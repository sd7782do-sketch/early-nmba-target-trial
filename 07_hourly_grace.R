## =============================================================================
## 07_hourly_grace.R
## Griglia oraria delle covariate nella finestra di grace, per i pesi di censura
## tempo-varianti dello script 06.
##
## Riusa le definizioni di 03 e 05 invece di riscriverle: stessi ITEMS, stesso
## PaO2 arterioso da specimen_id, stessi ASOF join per il P/F, stessa conversione
## in equivalenti noradrenalinici. Nessun itemid cablato qui dentro.
##
## Output: OUT_DIR/initiation/hourly_grace.parquet
##   stay_id, hour (1..24), pf, peep, plateau, ne_equiv
## =============================================================================

suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source("_paths.R")

H_MAX     <- 24L
SPEC_ITEM <- 52033          # specimen type nel pannello emogas, come in 03
COHORT    <- "A"

## stesso allineamento dei PARAMS fatto in 03
if (is.null(PARAMS$elig_window_h) && !is.null(PARAMS$intubation_to_t0_h))
  PARAMS$elig_window_h <- PARAMS$intubation_to_t0_h
if (is.null(PARAMS$elig_window_h)) PARAMS$elig_window_h <- 48
P <- PARAMS

## ITEMS di gravita', definiti in 05 e non nel config
if (is.null(ITEMS$vaso_input) || !length(ITEMS$vaso_input))
  ITEMS$vaso_input <- c(221906, 221289, 221749, 222315, 221662)

req <- c("peep", "fio2_ce", "pao2_le", "plateau", "vaso_input")
bad <- req[vapply(req, function(k) is.null(ITEMS[[k]]) || !length(ITEMS[[k]]), TRUE)]
if (length(bad)) stop("ITEMS mancanti: ", paste(bad, collapse = ", "))

COH <- file.path(OUT_DIR, "initiation", sprintf("cohort_%s_sev.parquet", COHORT))
if (!file.exists(COH)) stop("Manca ", COH)
OUT <- file.path(OUT_DIR, "initiation", "hourly_grace.parquet")
as_path <- function(p) gsub("\\\\", "/", normalizePath(p, mustWork = FALSE))
lbl <- function(ids) if (is.null(ids) || !length(ids)) "NULL" else paste(ids, collapse = ",")

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
if (!inherits(con, "duckdb_connection")) stop("Connessione DuckDB non creata.")
ddb_exec <- function(sql) dbExecute(con, sql)
ddb_get  <- function(sql) dbGetQuery(con, sql)

ddb_exec("SET memory_limit='12GB'")
ddb_exec("SET threads=6")
ddb_exec(sprintf("SET temp_directory='%s'", gsub("\\\\", "/", file.path(normalizePath(tempdir()), "duckdb_spill"))))
ddb_exec("SET preserve_insertion_order=false")

## ---- sorgenti: gli stessi sottoinsiemi usati da 03 e 05 ---------------------
CE_PAR <- ensure_subset(con, "chartevents_respiratory\\.parquet$", "chartevents_respiratory.parquet",
                        "chartevents", c("stay_id", "itemid", "charttime", "valuenum"),
                        unique(c(ITEMS$peep, ITEMS$fio2_ce, ITEMS$plateau, ITEMS$vt)))
LE_PAR <- ensure_subset(con, "labevents_bloodgas\\.parquet$", "labevents_bloodgas.parquet",
                        "labevents",
                        c("subject_id", "hadm_id", "specimen_id", "itemid", "charttime", "value", "valuenum"),
                        unique(c(ITEMS$pao2_le, SPEC_ITEM)))
cat("chartevents:", CE_PAR, "\nlabevents  :", LE_PAR, "\n")

ddb_exec(sprintf("CREATE VIEW ce AS SELECT * FROM read_parquet('%s')", CE_PAR))
ddb_exec(sprintf("CREATE VIEW le AS SELECT * FROM read_parquet('%s')", LE_PAR))
ddb_exec(sprintf("CREATE TABLE icustays AS SELECT * FROM %s", src("icustays")))
ddb_exec(sprintf("CREATE TABLE inputev AS SELECT * FROM %s", src("inputevents")))
ddb_exec(sprintf("CREATE TABLE coh AS SELECT * FROM read_parquet('%s')", as_path(COH)))
cat("Coorte:", ddb_get("SELECT COUNT(*) AS n FROM coh")$n, "stay\n")

## ---- scheletro stay per ora --------------------------------------------------
ddb_exec(sprintf("CREATE TABLE skel AS SELECT c.stay_id, h.hour, c.time_zero + to_hours(h.hour - 1) AS h_start, c.time_zero + to_hours(h.hour) AS h_end, c.pf AS pf0, c.peep AS peep0, c.plateau AS plateau0, COALESCE(c.ne_equiv, 0) AS ne0 FROM coh c CROSS JOIN (SELECT UNNEST(generate_series(1, %d)) AS hour) h", H_MAX))

## ---- 1. P/F: stessa costruzione di 03 ---------------------------------------
le_cols <- cols_of(con, sprintf("read_parquet('%s')", LE_PAR))
arterial_sql <- if ("specimen_status" %in% le_cols) {
  sprintf("SELECT subject_id, hadm_id, charttime, valuenum FROM le WHERE itemid IN (%s) AND specimen_status = 'documented_arterial'", lbl(ITEMS$pao2_le))
} else {
  sprintf("WITH spec AS (SELECT DISTINCT specimen_id, trim(value) AS spec_type FROM le WHERE itemid = %d) SELECT p.subject_id, p.hadm_id, p.charttime, p.valuenum FROM le p JOIN spec s USING (specimen_id) WHERE p.itemid IN (%s) AND s.spec_type = 'ART.'", SPEC_ITEM, lbl(ITEMS$pao2_le))
}
ddb_exec(sprintf("CREATE TABLE pao2 AS SELECT i.stay_id, a.charttime, a.valuenum AS pao2 FROM (%s) a JOIN icustays i ON i.hadm_id = a.hadm_id WHERE a.valuenum BETWEEN 20 AND 700 AND a.charttime BETWEEN i.intime - INTERVAL 6 HOUR AND i.outtime", arterial_sql))
ddb_exec(sprintf("CREATE TABLE fio2 AS SELECT * FROM (SELECT stay_id, charttime, CASE WHEN valuenum <= 1 THEN valuenum*100 ELSE valuenum END AS fio2 FROM ce WHERE itemid IN (%s)) WHERE fio2 BETWEEN 21 AND 100", lbl(ITEMS$fio2_ce)))
ddb_exec(sprintf("CREATE TABLE peep AS SELECT stay_id, charttime, valuenum AS peep FROM ce WHERE itemid IN (%s) AND valuenum BETWEEN 0 AND 40", lbl(ITEMS$peep)))
ddb_exec(sprintf("CREATE TABLE pf AS WITH a AS (SELECT p.stay_id, p.charttime, p.pao2, f.fio2, f.charttime f_time FROM pao2 p ASOF LEFT JOIN fio2 f ON p.stay_id = f.stay_id AND p.charttime >= f.charttime), b AS (SELECT a.*, e.peep, e.charttime e_time FROM a ASOF LEFT JOIN peep e ON a.stay_id = e.stay_id AND a.charttime >= e.charttime) SELECT stay_id, charttime, 100.0*pao2/fio2 AS pf FROM b WHERE fio2 IS NOT NULL AND peep IS NOT NULL AND date_diff('minute', f_time, charttime) <= %d AND date_diff('minute', e_time, charttime) <= %d", P$pair_window_h*60, P$pair_window_h*60))
cat("Emogas con P/F calcolabile:", ddb_get("SELECT COUNT(*) AS n FROM pf")$n, "\n")

## ---- 2. valori orari: ultimo dentro ciascuna ora ----------------------------
ddb_exec("CREATE TABLE hv_pf AS SELECT s.stay_id, s.hour, (SELECT x.pf FROM pf x WHERE x.stay_id = s.stay_id AND x.charttime > s.h_start AND x.charttime <= s.h_end ORDER BY x.charttime DESC LIMIT 1) AS v FROM skel s")
ddb_exec("CREATE TABLE hv_peep AS SELECT s.stay_id, s.hour, (SELECT x.peep FROM peep x WHERE x.stay_id = s.stay_id AND x.charttime > s.h_start AND x.charttime <= s.h_end ORDER BY x.charttime DESC LIMIT 1) AS v FROM skel s")
ddb_exec(sprintf("CREATE TABLE plat AS SELECT stay_id, charttime, valuenum AS plateau FROM ce WHERE itemid IN (%s) AND valuenum IS NOT NULL", lbl(ITEMS$plateau)))
ddb_exec("CREATE TABLE hv_plateau AS SELECT s.stay_id, s.hour, (SELECT x.plateau FROM plat x WHERE x.stay_id = s.stay_id AND x.charttime > s.h_start AND x.charttime <= s.h_end ORDER BY x.charttime DESC LIMIT 1) AS v FROM skel s")

## ---- 3. equivalenti noradrenalinici: stessa conversione di 05 ---------------
ddb_exec(sprintf("CREATE TABLE vaso_run AS SELECT i.stay_id, i.itemid, i.starttime, i.endtime, CASE WHEN lower(i.rateuom) LIKE 'mcg/kg/min%%' THEN i.rate WHEN lower(i.rateuom) LIKE 'mcg/kg/hour%%' THEN i.rate/60.0 WHEN lower(i.rateuom) LIKE 'mcg/min%%' THEN i.rate/COALESCE(i.patientweight, 80) WHEN lower(i.rateuom) LIKE 'mg/kg/min%%' THEN i.rate*1000 WHEN lower(i.rateuom) LIKE 'mg/min%%' THEN i.rate*1000/COALESCE(i.patientweight, 80) WHEN lower(i.rateuom) LIKE 'units/min%%' THEN i.rate*2.5 WHEN lower(i.rateuom) LIKE 'units/hour%%' THEN i.rate/60.0*2.5 ELSE NULL END AS dose_kgmin FROM inputev i SEMI JOIN coh c ON c.stay_id = i.stay_id WHERE i.itemid IN (%s) AND i.rate IS NOT NULL AND i.endtime > i.starttime", lbl(ITEMS$vaso_input)))
ddb_exec("CREATE TABLE hv_ne AS SELECT s.stay_id, s.hour, COALESCE((SELECT SUM(v.dose_kgmin) FROM vaso_run v WHERE v.stay_id = s.stay_id AND v.dose_kgmin IS NOT NULL AND v.starttime < s.h_end AND v.endtime > s.h_start), 0) AS v FROM skel s")

## ---- 4. assemblaggio, ancoraggio al basale, riporto in avanti ---------------
ddb_exec("CREATE TABLE raw_grid AS SELECT s.stay_id, s.hour, CASE WHEN s.hour = 1 THEN COALESCE(pf.v, s.pf0) ELSE pf.v END AS pf, CASE WHEN s.hour = 1 THEN COALESCE(pe.v, s.peep0) ELSE pe.v END AS peep, CASE WHEN s.hour = 1 THEN COALESCE(pl.v, s.plateau0) ELSE pl.v END AS plateau, CASE WHEN s.hour = 1 THEN COALESCE(ne.v, s.ne0) ELSE ne.v END AS ne_equiv FROM skel s LEFT JOIN hv_pf pf ON pf.stay_id = s.stay_id AND pf.hour = s.hour LEFT JOIN hv_peep pe ON pe.stay_id = s.stay_id AND pe.hour = s.hour LEFT JOIN hv_plateau pl ON pl.stay_id = s.stay_id AND pl.hour = s.hour LEFT JOIN hv_ne ne ON ne.stay_id = s.stay_id AND ne.hour = s.hour")

ddb_exec(sprintf("COPY (SELECT stay_id, hour, last_value(pf IGNORE NULLS) OVER w AS pf, last_value(peep IGNORE NULLS) OVER w AS peep, last_value(plateau IGNORE NULLS) OVER w AS plateau, last_value(ne_equiv IGNORE NULLS) OVER w AS ne_equiv FROM raw_grid WINDOW w AS (PARTITION BY stay_id ORDER BY hour ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW) ORDER BY stay_id, hour) TO '%s' (FORMAT PARQUET)", as_path(OUT)))

## ---- 5. controlli ------------------------------------------------------------
cat("\n=== 1. Misure NUOVE per ora, prima del riporto in avanti ===\n")
print(ddb_get("SELECT hour, ROUND(AVG(CASE WHEN pf IS NOT NULL THEN 1 ELSE 0 END),3) pf, ROUND(AVG(CASE WHEN peep IS NOT NULL THEN 1 ELSE 0 END),3) peep, ROUND(AVG(CASE WHEN plateau IS NOT NULL THEN 1 ELSE 0 END),3) plateau FROM raw_grid WHERE hour <= 12 GROUP BY hour ORDER BY hour"))

cat("\n=== 2. Ora 1 contro basale: deve essere vicino a 1 ===\n")
print(ddb_get(sprintf("SELECT ROUND(AVG(CASE WHEN ABS(g.pf - c.pf) < 1e-6 THEN 1 ELSE 0 END),3) pf, ROUND(AVG(CASE WHEN ABS(g.peep - c.peep) < 1e-6 THEN 1 ELSE 0 END),3) peep FROM read_parquet('%s') g JOIN coh c USING (stay_id) WHERE g.hour = 1", as_path(OUT))))

cat("\n=== 3. Quota di stay in cui il valore cambia dentro la finestra ===\n")
print(ddb_get(sprintf("SELECT ROUND(AVG(CASE WHEN n_pf > 1 THEN 1 ELSE 0 END),3) pf, ROUND(AVG(CASE WHEN n_peep > 1 THEN 1 ELSE 0 END),3) peep, ROUND(AVG(CASE WHEN n_plat > 1 THEN 1 ELSE 0 END),3) plateau, ROUND(AVG(CASE WHEN n_ne > 1 THEN 1 ELSE 0 END),3) ne FROM (SELECT stay_id, COUNT(DISTINCT pf) n_pf, COUNT(DISTINCT peep) n_peep, COUNT(DISTINCT plateau) n_plat, COUNT(DISTINCT ne_equiv) n_ne FROM read_parquet('%s') WHERE hour <= 12 GROUP BY stay_id)", as_path(OUT))))

cat("\n=== 4. Variazione ora 12 meno ora 1, per finestra di inizio ===\n")
print(ddb_get(sprintf("WITH a AS (SELECT * FROM read_parquet('%s') WHERE hour IN (1,12)), d AS (SELECT x.stay_id, x.pf - y.pf dpf, x.peep - y.peep dpeep, x.plateau - y.plateau dplat, x.ne_equiv - y.ne_equiv dne FROM a x JOIN a y USING (stay_id) WHERE x.hour = 12 AND y.hour = 1) SELECT CASE WHEN c.h_to_init <= 12 THEN 'early' ELSE 'no early' END grp, COUNT(*) n, ROUND(MEDIAN(d.dpf),1) d_pf, ROUND(MEDIAN(d.dpeep),1) d_peep, ROUND(MEDIAN(d.dplat),1) d_plateau, ROUND(MEDIAN(d.dne),3) d_ne FROM d JOIN coh c USING (stay_id) WHERE c.prevalent_user = 0 AND c.sed_active_t0 > 0 GROUP BY 1", as_path(OUT))))

dbDisconnect(con, shutdown = TRUE)
cat("\nScritto:", OUT, "\n")
