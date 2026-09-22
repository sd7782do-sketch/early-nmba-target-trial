## =============================================================================
## 05_baseline_severity.R
## Aggancia a time zero le covariate di gravità mancanti e riscrive la coorte:
##   - vasopressori attivi + dose in equivalenti di noradrenalina (mcg/kg/min)
##   - inotropi attivi
##   - lattato, creatinina, bilirubina, piastrine (ultimo valore prima di t0)
##   - plateau pressure e driving pressure quando disponibili
## Tutto valutato SOLO fino a t0: nessuna informazione futura.
## Output: cohort_<X>_sev.parquet, usato dagli script 07, 08 e 06.
## =============================================================================

suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source("_paths.R")

## =============================================================================
## ITEMS di gravità, definiti QUI e non nel config: 01_config_corrected.R
## riassegna ITEMS al suo interno e cancellerebbe qualunque riga aggiunta in fondo.
## I valori sotto sono gli itemid abituali di MIMIC-IV e vengono VERIFICATI
## contro d_items / d_labitems più sotto: se un'etichetta non combacia, lo script
## si ferma invece di procedere con un item sbagliato.
## =============================================================================
ITEMS$vaso_input     <- c(221906,  # norepinephrine
                          221289,  # epinephrine
                          221749,  # phenylephrine
                          222315,  # vasopressin
                          221662)  # dopamine
ITEMS$inotrope_input <- c(221653,  # dobutamine
                          221986)  # milrinone
ITEMS$lactate_le     <- 50813      # lactate
ITEMS$creat_le       <- 50912      # creatinine
ITEMS$bili_le        <- 50885      # bilirubin, total
ITEMS$plt_le         <- 51265      # platelet count

VERIFY_ITEMS <- TRUE

COHORT <- "A"
LAC_WINDOW_H  <- 6     # lattato: quanto indietro accettare
CHEM_WINDOW_H <- 24    # creatinina/bilirubina/piastrine
VENT_WINDOW_H <- 6     # plateau / tidal volume

need <- c("vaso_input","lactate_le")
miss <- need[vapply(need, function(k) is.null(ITEMS[[k]]) || !length(ITEMS[[k]]), TRUE)]
if (length(miss))
  stop("Mancano ITEMS: ", paste(miss, collapse = ", "),
       "\nEseguire 00b_inventory_severity.R e completare il config.")

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con, "SET memory_limit='12GB'"); dbExecute(con, "SET threads=6")
dbExecute(con, sprintf("SET temp_directory='%s'",
                       gsub("\\\\","/", file.path(normalizePath(tempdir()), "duckdb_spill"))))
dbExecute(con, "SET preserve_insertion_order=false")

## ---- verifica delle etichette prima di usare gli itemid ---------------------
if (VERIFY_ITEMS) {
  want_ce <- c("221906"="norepinephrine", "221289"="epinephrine",
               "221749"="phenylephrine",  "222315"="vasopressin",
               "221662"="dopamine",       "221653"="dobutamine",
               "221986"="milrinone")
  want_le <- c("50813"="lactate", "50912"="creatinine",
               "50885"="bilirubin", "51265"="platelet")

  di <- dbGetQuery(con, sprintf("SELECT itemid, label FROM %s WHERE itemid IN (%s)",
                                src("d_items"), paste(names(want_ce), collapse = ",")))
  dl <- dbGetQuery(con, sprintf("SELECT itemid, label FROM %s WHERE itemid IN (%s)",
                                src("d_labitems"), paste(names(want_le), collapse = ",")))
  names(di) <- tolower(names(di)); names(dl) <- tolower(names(dl))

  chk <- function(tab, want, what) {
    bad <- character(0)
    for (id in names(want)) {
      lab <- tab$label[tab$itemid == as.integer(id)]
      if (!length(lab) || !grepl(want[[id]], lab, ignore.case = TRUE))
        bad <- c(bad, sprintf("  %s: atteso '%s', trovato '%s'", id, want[[id]],
                              if (length(lab)) lab[1] else "<assente>"))
    }
    if (length(bad))
      stop("itemid ", what, " non verificati:\n", paste(bad, collapse = "\n"),
           "\nCorreggere in cima a 05_baseline_severity.R usando\n",
           "  feasibility_out/inventory_severity/items_severity_candidates.csv")
    cat("itemid", what, "verificati.\n")
  }
  chk(di, want_ce, "ICU (inputevents)")
  chk(dl, want_le, "lab")
}

COH <- file.path(OUT_DIR, "initiation", sprintf("cohort_%s.parquet", COHORT))
if (!file.exists(COH)) stop("Manca ", COH)
dbExecute(con, sprintf("CREATE TABLE coh AS SELECT * FROM read_parquet('%s')",
                       gsub("\\\\","/", normalizePath(COH))))
dbExecute(con, sprintf("CREATE TABLE inputev AS SELECT * FROM %s", src("inputevents")))

## ---- sottoinsieme labevents per gli analiti di gravità ----------------------
lab_ids <- unique(c(ITEMS$lactate_le, ITEMS$creat_le, ITEMS$bili_le, ITEMS$plt_le))
LAB_PAR <- ensure_subset(con, "labevents_severity.*\\.parquet$",
                         "labevents_severity.parquet", "labevents",
                         c("subject_id","hadm_id","itemid","charttime","valuenum"),
                         lab_ids)
dbExecute(con, sprintf("CREATE VIEW lab AS SELECT * FROM read_parquet('%s')", LAB_PAR))

## =============================================================================
## 1. Vasopressori attivi a t0, in equivalenti di noradrenalina
##    ATTENZIONE: le conversioni sotto vanno confrontate con vaso_rateuom_audit.csv.
##    Se le unità osservate non sono quelle previste, il CASE va riscritto.
## =============================================================================
lbl <- function(ids) if (is.null(ids) || !length(ids)) "NULL" else paste(ids, collapse = ",")

dbExecute(con, sprintf("
CREATE TABLE vaso_t0 AS
WITH run AS (
  SELECT c.stay_id, i.itemid, i.rate, i.rateuom,
         COALESCE(i.patientweight, 80) AS wt
  FROM coh c JOIN inputev i ON i.stay_id = c.stay_id
  WHERE i.itemid IN (%s) AND i.rate IS NOT NULL
    AND i.starttime <= c.time_zero AND i.endtime >= c.time_zero
), norm AS (
  SELECT stay_id, itemid,
         CASE
           WHEN lower(rateuom) LIKE 'mcg/kg/min%%'  THEN rate
           WHEN lower(rateuom) LIKE 'mcg/kg/hour%%' THEN rate/60.0
           WHEN lower(rateuom) LIKE 'mcg/min%%'     THEN rate/wt
           WHEN lower(rateuom) LIKE 'mg/kg/min%%'   THEN rate*1000
           WHEN lower(rateuom) LIKE 'mg/min%%'      THEN rate*1000/wt
           WHEN lower(rateuom) LIKE 'units/min%%'   THEN rate*2.5     -- vasopressina
           WHEN lower(rateuom) LIKE 'units/hour%%'  THEN rate/60.0*2.5
           ELSE NULL END AS dose_kgmin
  FROM run
)
SELECT stay_id,
       1 AS vaso_active,
       SUM(dose_kgmin) AS ne_equiv_raw,
       count(DISTINCT itemid) AS n_vaso
FROM norm WHERE dose_kgmin IS NOT NULL
GROUP BY stay_id", lbl(ITEMS$vaso_input)))

if (!is.null(ITEMS$inotrope_input) && length(ITEMS$inotrope_input)) {
  dbExecute(con, sprintf("
  CREATE TABLE ino_t0 AS
  SELECT DISTINCT c.stay_id, 1 AS inotrope_active
  FROM coh c JOIN inputev i ON i.stay_id = c.stay_id
  WHERE i.itemid IN (%s) AND i.rate IS NOT NULL
    AND i.starttime <= c.time_zero AND i.endtime >= c.time_zero",
  lbl(ITEMS$inotrope_input)))
} else {
  dbExecute(con, "CREATE TABLE ino_t0 AS SELECT NULL::BIGINT stay_id, NULL::INT inotrope_active WHERE FALSE")
}

## =============================================================================
## 2. Analiti di laboratorio: ultimo valore prima di t0 (ASOF)
## =============================================================================
lab_last <- function(alias, ids, win_h) {
  if (is.null(ids) || !length(ids)) {
    dbExecute(con, sprintf("CREATE TABLE lab_%s AS
      SELECT NULL::BIGINT stay_id, NULL::DOUBLE %s WHERE FALSE", alias, alias))
    return(invisible(NULL))
  }
  dbExecute(con, sprintf("
  CREATE TABLE lab_%s AS
  WITH v AS (
    SELECT c.stay_id, l.charttime, l.valuenum
    FROM coh c JOIN lab l ON l.hadm_id = c.hadm_id
    WHERE l.itemid IN (%s) AND l.valuenum IS NOT NULL
      AND l.charttime <= c.time_zero
      AND l.charttime >= c.time_zero - INTERVAL %d HOUR
  )
  SELECT stay_id, valuenum AS %s FROM (
    SELECT *, ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY charttime DESC) rn FROM v
  ) WHERE rn = 1", alias, lbl(ids), win_h, alias))
}

lab_last("lactate", ITEMS$lactate_le, LAC_WINDOW_H)
lab_last("creat",   ITEMS$creat_le,   CHEM_WINDOW_H)
lab_last("bili",    ITEMS$bili_le,    CHEM_WINDOW_H)
lab_last("plt",     ITEMS$plt_le,     CHEM_WINDOW_H)

## =============================================================================
## 3. Meccanica ventilatoria: plateau e driving pressure
## =============================================================================
CE_PAR <- ensure_subset(con, "chartevents_respiratory\\.parquet$",
                        "chartevents_respiratory.parquet", "chartevents",
                        c("stay_id","itemid","charttime","valuenum"),
                        unique(c(ITEMS$peep, ITEMS$fio2_ce, ITEMS$plateau, ITEMS$vt)))
dbExecute(con, sprintf("CREATE VIEW ce AS SELECT * FROM read_parquet('%s')", CE_PAR))

ce_last <- function(alias, ids) {
  if (is.null(ids) || !length(ids)) {
    dbExecute(con, sprintf("CREATE TABLE ce_%s AS
      SELECT NULL::BIGINT stay_id, NULL::DOUBLE %s WHERE FALSE", alias, alias))
    return(invisible(NULL))
  }
  dbExecute(con, sprintf("
  CREATE TABLE ce_%s AS
  SELECT stay_id, valuenum AS %s FROM (
    SELECT c.stay_id, x.valuenum,
           ROW_NUMBER() OVER (PARTITION BY c.stay_id ORDER BY x.charttime DESC) rn
    FROM coh c JOIN ce x ON x.stay_id = c.stay_id
    WHERE x.itemid IN (%s) AND x.valuenum IS NOT NULL
      AND x.charttime <= c.time_zero
      AND x.charttime >= c.time_zero - INTERVAL %d HOUR
  ) WHERE rn = 1", alias, alias, lbl(ids), VENT_WINDOW_H))
}
ce_last("plateau", ITEMS$plateau)
ce_last("vt",      ITEMS$vt)

## =============================================================================
## 4. Coorte arricchita
## =============================================================================
dbExecute(con, "
CREATE TABLE coh_sev AS
SELECT c.*,
       COALESCE(v.vaso_active, 0)    AS vaso_active,
       v.ne_equiv_raw                AS ne_equiv,
       COALESCE(v.n_vaso, 0)         AS n_vaso,
       COALESCE(i.inotrope_active,0) AS inotrope_active,
       l1.lactate, l2.creat, l3.bili, l4.plt,
       p.plateau, t.vt,
       CASE WHEN p.plateau IS NOT NULL THEN p.plateau - c.peep END AS driving_pressure
FROM coh c
LEFT JOIN vaso_t0  v  USING (stay_id)
LEFT JOIN ino_t0   i  USING (stay_id)
LEFT JOIN lab_lactate l1 USING (stay_id)
LEFT JOIN lab_creat   l2 USING (stay_id)
LEFT JOIN lab_bili    l3 USING (stay_id)
LEFT JOIN lab_plt     l4 USING (stay_id)
LEFT JOIN ce_plateau  p  USING (stay_id)
LEFT JOIN ce_vt       t  USING (stay_id)")

## ---- completezza e squilibrio grezzo ----------------------------------------
compl <- dbGetQuery(con, sprintf("
SELECT count(*) n,
       avg(vaso_active) p_vaso,
       avg(CASE WHEN ne_equiv IS NOT NULL THEN 1.0 ELSE 0 END) p_dose_calcolabile,
       avg(CASE WHEN lactate IS NOT NULL THEN 1.0 ELSE 0 END)  p_lattato,
       avg(CASE WHEN creat   IS NOT NULL THEN 1.0 ELSE 0 END)  p_creatinina,
       avg(CASE WHEN bili    IS NOT NULL THEN 1.0 ELSE 0 END)  p_bilirubina,
       avg(CASE WHEN plt     IS NOT NULL THEN 1.0 ELSE 0 END)  p_piastrine,
       avg(CASE WHEN plateau IS NOT NULL THEN 1.0 ELSE 0 END)  p_plateau
FROM coh_sev WHERE prevalent_user = 0 AND sed_active_t0 > 0"))
print(t(compl))

bal <- dbGetQuery(con, sprintf("
SELECT CASE WHEN h_to_init <= 12 THEN 1 ELSE 0 END AS treat,
       count(*) n,
       avg(vaso_active) p_vaso, median(ne_equiv) ne_mediana,
       median(lactate) lattato_mediana, median(plateau) plateau_mediano,
       median(driving_pressure) dp_mediana, avg(death_28) mort28
FROM coh_sev WHERE prevalent_user = 0 AND sed_active_t0 > 0
GROUP BY 1"))
cat("\n--- squilibrio grezzo sulle nuove covariate ---\n"); print(bal)

out <- file.path(OUT_DIR, "initiation", sprintf("cohort_%s_sev.parquet", COHORT))
dbExecute(con, sprintf("COPY (SELECT * FROM coh_sev) TO '%s' (FORMAT PARQUET)",
                       gsub("\\\\","/", out)))
write.csv(compl, file.path(OUT_DIR, "initiation", "severity_completeness.csv"), row.names = FALSE)
write.csv(bal,   file.path(OUT_DIR, "initiation", "severity_balance.csv"), row.names = FALSE)

dbDisconnect(con, shutdown = TRUE)
cat("\nScritto:", out, "\nGli script 07, 08 e 06 lo useranno come input.\n")
