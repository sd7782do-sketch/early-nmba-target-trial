## =============================================================================
## 03_initiation_counts.R  (v2)
## Conteggio initiation-only. Nessun filtro di durata nella definizione di esposizione.
##   A (primaria)    t0 = primo P/F <150 con PEEP >=5
##   B (sensibilità) t0 = misurazione di conferma (>=4h, <=24h dopo A), new-user
## Correzioni v2: percorsi risolti dinamicamente; derivazione del sangue arterioso
## da specimen_id quando la colonna specimen_status non è nel parquet.
## =============================================================================

suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source("_paths.R")

## ---- compatibilità con 01_config_corrected.R --------------------------------
## Il config corretto usa nomi diversi per alcuni PARAMS. Qui si riallineano,
## senza sovrascrivere nulla che sia già definito.
if (is.null(PARAMS$elig_window_h) && !is.null(PARAMS$intubation_to_t0_h))
  PARAMS$elig_window_h <- PARAMS$intubation_to_t0_h
if (is.null(PARAMS$elig_window_h)) PARAMS$elig_window_h <- 48

req_par <- c("age_min","pf_threshold","peep_min","pair_window_h",
             "elig_window_h","vent_gap_h","inf_gap_h","fu_days")
bad <- req_par[vapply(req_par, function(k)
  is.null(PARAMS[[k]]) || length(PARAMS[[k]]) != 1 || is.na(PARAMS[[k]]), TRUE)]
if (length(bad)) stop("PARAMS non validi (NULL, mancanti o di lunghezza != 1): ",
                      paste(bad, collapse = ", "))

req_items <- c("peep","fio2_ce","pao2_le","vent_proc","nmba_input","sedation_input")
bad_it <- req_items[vapply(req_items, function(k)
  is.null(ITEMS[[k]]) || !length(ITEMS[[k]]), TRUE)]
if (length(bad_it)) stop("ITEMS mancanti: ", paste(bad_it, collapse = ", "))

cat("Finestra eleggibilità:", PARAMS$elig_window_h, "ore dall'intubazione\n")

P         <- PARAMS
GRACE     <- c(6, 12, 24)
CONF_LO   <- 4
CONF_HI   <- 24
SPEC_ITEM <- 52033          # specimen type nel pannello emogas

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con, "SET memory_limit='12GB'"); dbExecute(con, "SET threads=6")
dbExecute(con, sprintf("SET temp_directory='%s'",
                       gsub("\\\\","/", file.path(normalizePath(tempdir()), "duckdb_spill"))))
dbExecute(con, "SET preserve_insertion_order=false")

## ---- sorgenti (riusa i parquet se ci sono, altrimenti li crea) ---------------
CE_PAR <- ensure_subset(con, "chartevents_oxygenation\\.parquet$", "chartevents_oxygenation.parquet",
                        "chartevents", c("stay_id","itemid","charttime","valuenum"),
                        unique(c(ITEMS$peep, ITEMS$fio2_ce)))
LE_PAR <- ensure_subset(con, "labevents_bloodgas\\.parquet$", "labevents_bloodgas.parquet",
                        "labevents",
                        c("subject_id","hadm_id","specimen_id","itemid","charttime","value","valuenum"),
                        unique(c(ITEMS$pao2_le, SPEC_ITEM)))

dbExecute(con, sprintf("CREATE VIEW ce AS SELECT * FROM read_parquet('%s')", CE_PAR))
dbExecute(con, sprintf("CREATE VIEW le AS SELECT * FROM read_parquet('%s')", LE_PAR))

for (t in c("icustays","patients","procedureevents","inputevents"))
  dbExecute(con, sprintf("CREATE TABLE %s AS SELECT * FROM %s",
                         ifelse(t == "inputevents", "inputev", t), src(t)))

## ---- PaO2 arterioso ----------------------------------------------------------
le_cols <- cols_of(con, sprintf("read_parquet('%s')", LE_PAR))

arterial_sql <- if ("specimen_status" %in% le_cols) {
  cat("Uso specimen_status già presente nel parquet.\n")
  sprintf("SELECT subject_id, hadm_id, charttime, valuenum
           FROM le WHERE itemid IN (%s) AND specimen_status = 'documented_arterial'",
          sqlids(ITEMS$pao2_le))
} else if (all(c("specimen_id","value") %in% le_cols)) {
  cat("specimen_status assente: derivo l'arterioso da specimen_id + itemid", SPEC_ITEM, "\n")
  sprintf("WITH spec AS (SELECT DISTINCT specimen_id, trim(value) AS spec_type
                         FROM le WHERE itemid = %d)
           SELECT p.subject_id, p.hadm_id, p.charttime, p.valuenum
           FROM le p JOIN spec s USING (specimen_id)
           WHERE p.itemid IN (%s) AND s.spec_type = 'ART.'",
          SPEC_ITEM, sqlids(ITEMS$pao2_le))
} else {
  stop("Il parquet di labevents non ha né specimen_status né specimen_id+value.\n",
       "Cancellalo e rilancia: ensure_subset lo ricrea con le colonne giuste.")
}

dbExecute(con, sprintf("
CREATE TABLE pao2 AS
SELECT i.stay_id, a.charttime, a.valuenum AS pao2
FROM (%s) a
JOIN icustays i ON i.hadm_id = a.hadm_id
WHERE a.valuenum BETWEEN 20 AND 700
  AND a.charttime BETWEEN i.intime - INTERVAL 6 HOUR AND i.outtime", arterial_sql))

cat("Righe PaO2 arteriose in ICU:", dbGetQuery(con, "SELECT count(*) n FROM pao2")$n, "\n")

## ---- FiO2 / PEEP / P/F -------------------------------------------------------
dbExecute(con, sprintf("
CREATE TABLE fio2 AS
SELECT * FROM (SELECT stay_id, charttime,
                      CASE WHEN valuenum <= 1 THEN valuenum*100 ELSE valuenum END AS fio2
               FROM ce WHERE itemid IN (%s))
WHERE fio2 BETWEEN 21 AND 100", sqlids(ITEMS$fio2_ce)))

dbExecute(con, sprintf("
CREATE TABLE peep AS SELECT stay_id, charttime, valuenum AS peep
FROM ce WHERE itemid IN (%s) AND valuenum BETWEEN 0 AND 40", sqlids(ITEMS$peep)))

dbExecute(con, sprintf("
CREATE TABLE pf AS
WITH a AS (SELECT p.stay_id, p.charttime, p.pao2, f.fio2, f.charttime f_time
           FROM pao2 p ASOF LEFT JOIN fio2 f
             ON p.stay_id = f.stay_id AND p.charttime >= f.charttime),
b AS (SELECT a.*, e.peep, e.charttime e_time
      FROM a ASOF LEFT JOIN peep e
        ON a.stay_id = e.stay_id AND a.charttime >= e.charttime)
SELECT stay_id, charttime, pao2, fio2, peep, 100.0*pao2/fio2 AS pf
FROM b
WHERE fio2 IS NOT NULL AND peep IS NOT NULL
  AND date_diff('minute', f_time, charttime) <= %d
  AND date_diff('minute', e_time, charttime) <= %d",
  P$pair_window_h*60, P$pair_window_h*60))

## ---- ventilazione invasiva ---------------------------------------------------
dbExecute(con, sprintf("
CREATE TABLE vent_ep AS
WITH raw AS (SELECT stay_id, starttime, endtime FROM procedureevents
             WHERE itemid IN (%s) AND endtime > starttime),
ord AS (SELECT *, LAG(endtime) OVER (PARTITION BY stay_id ORDER BY starttime) prev_end FROM raw),
flg AS (SELECT *, CASE WHEN prev_end IS NULL OR starttime > prev_end + INTERVAL %d HOUR
                       THEN 1 ELSE 0 END new_ep FROM ord),
grp AS (SELECT *, SUM(new_ep) OVER (PARTITION BY stay_id ORDER BY starttime
                                    ROWS UNBOUNDED PRECEDING) ep FROM flg)
SELECT stay_id, ep, MIN(starttime) vent_start, MAX(endtime) vent_end
FROM grp GROUP BY stay_id, ep", sqlids(ITEMS$vent_proc), P$vent_gap_h))

dbExecute(con, "CREATE TABLE vent_first AS
  SELECT * FROM (SELECT *, ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY vent_start) rn
                 FROM vent_ep) WHERE rn = 1")

dbExecute(con, sprintf("
CREATE TABLE qual AS
SELECT v.stay_id, pf.charttime, pf.pao2, pf.fio2, pf.peep, pf.pf, v.vent_start, v.vent_end
FROM vent_first v JOIN pf ON pf.stay_id = v.stay_id
WHERE pf.pf < %f AND pf.peep >= %f
  AND pf.charttime >= v.vent_start
  AND pf.charttime <= v.vent_start + INTERVAL %d HOUR
  AND pf.charttime <= v.vent_end", P$pf_threshold, P$peep_min, P$elig_window_h))

dbExecute(con, "
CREATE TABLE t0_A AS
SELECT stay_id, charttime AS time_zero, pao2, fio2, peep, pf, vent_start, vent_end
FROM (SELECT *, ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY charttime) rn FROM qual)
WHERE rn = 1")

dbExecute(con, sprintf("
CREATE TABLE t0_B AS
SELECT stay_id, time_zero, pao2, fio2, peep, pf, vent_start, vent_end FROM (
  SELECT q.stay_id, q.charttime AS time_zero, q.pao2, q.fio2, q.peep, q.pf,
         q.vent_start, q.vent_end,
         ROW_NUMBER() OVER (PARTITION BY q.stay_id ORDER BY q.charttime) rn
  FROM qual q JOIN t0_A a USING (stay_id)
  WHERE q.charttime >= a.time_zero + INTERVAL %d HOUR
    AND q.charttime <= a.time_zero + INTERVAL %d HOUR
) WHERE rn = 1", CONF_LO, CONF_HI))

## ---- infusioni NMBA: solo l'inizio -------------------------------------------
dbExecute(con, sprintf("
CREATE TABLE nmba_ep AS
WITH seg AS (SELECT stay_id, itemid, starttime, endtime FROM inputev
             WHERE itemid IN (%s) AND (%s) AND endtime > starttime),
ord AS (SELECT *, LAG(endtime) OVER (PARTITION BY stay_id, itemid ORDER BY starttime) prev_end FROM seg),
flg AS (SELECT *, CASE WHEN prev_end IS NULL OR starttime > prev_end + INTERVAL %d HOUR
                       THEN 1 ELSE 0 END new_ep FROM ord),
grp AS (SELECT *, SUM(new_ep) OVER (PARTITION BY stay_id, itemid ORDER BY starttime
                                    ROWS UNBOUNDED PRECEDING) ep FROM flg)
SELECT stay_id, itemid, MIN(starttime) inf_start, MAX(endtime) inf_end,
       date_diff('minute', MIN(starttime), MAX(endtime))/60.0 inf_hours
FROM grp GROUP BY stay_id, itemid, ep", sqlids(ITEMS$nmba_input), INFUSION_RULE, P$inf_gap_h))

## ---- coorti ------------------------------------------------------------------
build_cohort <- function(t0_tbl, suffix) {
  dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE coh_%s AS
  SELECT t.*, i.subject_id, i.hadm_id, i.intime, i.outtime,
         p.gender, p.dod, p.anchor_year_group,
         p.anchor_age + (date_part('year', i.intime) - p.anchor_year) AS age,
         (SELECT count(*) FROM nmba_ep n
           WHERE n.stay_id = t.stay_id
             AND n.inf_start <  t.time_zero AND n.inf_end > t.time_zero) AS prevalent_user,
         (SELECT count(*) FROM inputev s
           WHERE s.stay_id = t.stay_id AND s.itemid IN (%s) AND s.rate IS NOT NULL
             AND s.starttime <= t.time_zero AND s.endtime >= t.time_zero) AS sed_active_t0,
         (SELECT min(n.inf_start) FROM nmba_ep n
           WHERE n.stay_id = t.stay_id AND n.inf_start >= t.time_zero) AS first_init
  FROM %s t
  JOIN icustays i USING (stay_id)
  JOIN patients p USING (subject_id)
  WHERE p.anchor_age + (date_part('year', i.intime) - p.anchor_year) >= %d",
  suffix, sqlids(ITEMS$sedation_input), t0_tbl, P$age_min))

  dbExecute(con, sprintf("
  CREATE OR REPLACE TABLE coh_%s AS
  SELECT *, date_diff('minute', time_zero, first_init)/60.0 AS h_to_init,
         CASE WHEN dod IS NOT NULL AND dod <= CAST(time_zero AS DATE) + %d
              THEN 1 ELSE 0 END AS death_28
  FROM coh_%s", suffix, P$fu_days, suffix))
}
build_cohort("t0_A", "A"); build_cohort("t0_B", "B")

## ---- conteggi ----------------------------------------------------------------
count_variant <- function(suffix, label) do.call(rbind, lapply(GRACE, function(g)
  dbGetQuery(con, sprintf("
    SELECT '%s' AS coorte, %d AS grace_h,
           count(*) AS eleggibili,
           sum(CASE WHEN prevalent_user > 0 THEN 1 ELSE 0 END) AS prevalent_esclusi,
           sum(CASE WHEN prevalent_user = 0 THEN 1 ELSE 0 END) AS new_user,
           sum(CASE WHEN prevalent_user = 0 AND sed_active_t0 > 0 THEN 1 ELSE 0 END) AS nu_sedati,
           sum(CASE WHEN prevalent_user = 0 AND first_init IS NOT NULL THEN 1 ELSE 0 END) AS init_mai_filtro,
           sum(CASE WHEN prevalent_user = 0 AND h_to_init <= %d THEN 1 ELSE 0 END) AS init_entro_grace,
           sum(CASE WHEN prevalent_user = 0 AND sed_active_t0 > 0 AND h_to_init <= %d
                    THEN 1 ELSE 0 END) AS esposti_finali,
           sum(CASE WHEN prevalent_user = 0 AND sed_active_t0 > 0 AND h_to_init <= %d
                    AND death_28 = 1 THEN 1 ELSE 0 END) AS eventi_esposti
    FROM coh_%s", label, g, g, g, g, suffix))))

counts <- rbind(count_variant("A", "A_primo_punto"), count_variant("B", "B_conferma"))
print(counts)

cost <- dbGetQuery(con, "
SELECT (SELECT count(*) FROM coh_A) eleggibili_A,
       (SELECT count(*) FROM coh_B) eleggibili_B,
       (SELECT count(*) FROM coh_A) - (SELECT count(*) FROM coh_B) persi_per_conferma")
print(cost)

init_dist <- dbGetQuery(con, "
SELECT 'A' coorte, quantile_cont(h_to_init,0.25) q25, median(h_to_init) mediana,
       quantile_cont(h_to_init,0.75) q75, count(h_to_init) n
FROM coh_A WHERE prevalent_user = 0
UNION ALL
SELECT 'B', quantile_cont(h_to_init,0.25), median(h_to_init),
       quantile_cont(h_to_init,0.75), count(h_to_init)
FROM coh_B WHERE prevalent_user = 0")
print(init_dist)

## ---- ESS con anchor_year_group ----------------------------------------------
ess_kish <- function(w) sum(w)^2 / sum(w^2)
ess_tab <- do.call(rbind, lapply(c("A","B"), function(s) do.call(rbind, lapply(GRACE, function(g) {
  d <- dbGetQuery(con, sprintf("
    SELECT age, gender, anchor_year_group, pf, peep, fio2, death_28,
           CASE WHEN h_to_init <= %d THEN 1 ELSE 0 END AS treat
    FROM coh_%s WHERE prevalent_user = 0 AND sed_active_t0 > 0", g, s))
  if (sum(d$treat) < 20) return(NULL)
  d$gender <- factor(d$gender); d$anchor_year_group <- factor(d$anchor_year_group)
  ps <- fitted(glm(treat ~ age + gender + anchor_year_group + pf + peep + fio2,
                   family = binomial, data = d))
  w <- ifelse(d$treat == 1, 1 - ps, ps)
  data.frame(coorte = s, grace_h = g, n = nrow(d), n_treated = sum(d$treat),
             eventi_treated = sum(d$death_28[d$treat == 1]),
             ess_ow_totale = ess_kish(w), ess_ow_treated = ess_kish(w[d$treat == 1]),
             ps_tr_min = min(ps[d$treat == 1]), ps_tr_max = max(ps[d$treat == 1]))
}))))
print(ess_tab)

## ---- scrittura ---------------------------------------------------------------
od <- file.path(OUT_DIR, "initiation"); dir.create(od, showWarnings = FALSE, recursive = TRUE)
write.csv(counts,    file.path(od, "counts.csv"), row.names = FALSE)
write.csv(cost,      file.path(od, "confirmation_cost.csv"), row.names = FALSE)
write.csv(init_dist, file.path(od, "init_timing.csv"), row.names = FALSE)
if (!is.null(ess_tab)) write.csv(ess_tab, file.path(od, "ess.csv"), row.names = FALSE)
for (s in c("A","B"))
  dbExecute(con, sprintf("COPY (SELECT * FROM coh_%s) TO '%s' (FORMAT PARQUET)", s,
                         gsub("\\\\","/", file.path(od, sprintf("cohort_%s.parquet", s)))))

dbDisconnect(con, shutdown = TRUE)
cat("\nScritto in", od, "\n")
