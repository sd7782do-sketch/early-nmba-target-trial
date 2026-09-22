## =============================================================================
## 08_deathtime_and_audit.R
## Due compiti separati, entrambi rapidi:
##   A. aggiunge deathtime al parquet di coorte -> cohort_A_sev_dt.parquet
##   B. cerca l'origine della discrepanza 6.370 contro 6.360 del vecchio audit
##
## Da lanciare dalla radice del progetto dopo aver configurato `_paths.R`:
##   source("R/08_deathtime_and_audit.R")
## =============================================================================

suppressPackageStartupMessages({ library(DBI); library(duckdb) })
source("_paths.R")

COH_IN  <- file.path(OUT_DIR, "initiation", "cohort_A_sev.parquet")
COH_OUT <- file.path(OUT_DIR, "initiation", "cohort_A_sev_dt.parquet")
AUDIT   <- file.path(OUT_DIR, "audit"); dir.create(AUDIT, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(COH_IN)) stop("Manca ", COH_IN)
as_path <- function(p) gsub("\\\\", "/", normalizePath(p, mustWork = FALSE))

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
if (!inherits(con, "duckdb_connection")) stop("Connessione DuckDB non creata.")
ddb_exec <- function(sql) dbExecute(con, sql)
ddb_get  <- function(sql) dbGetQuery(con, sql)
ddb_exec("SET threads=6")
ddb_exec("SET preserve_insertion_order=false")

## =============================================================================
## A. deathtime
## =============================================================================
cat("\n########## A. join di deathtime ##########\n")
ddb_exec(sprintf("CREATE TABLE coh AS SELECT * FROM read_parquet('%s')", as_path(COH_IN)))
ddb_exec(sprintf("CREATE TABLE adm AS SELECT * FROM %s", src("admissions")))

## deathtime e' non nullo solo nel ricovero in cui il decesso e' avvenuto,
## che puo' non essere quello indice: si aggancia per subject_id.
ddb_exec("CREATE TABLE dt AS SELECT subject_id, MIN(deathtime) AS deathtime FROM adm WHERE deathtime IS NOT NULL GROUP BY subject_id")

ddb_exec(sprintf("COPY (SELECT c.*, d.deathtime FROM coh c LEFT JOIN dt d ON d.subject_id = c.subject_id) TO '%s' (FORMAT PARQUET)", as_path(COH_OUT)))

cat("\n--- copertura e coerenza con dod\n")
print(ddb_get(sprintf("SELECT COUNT(*) AS n_stay, SUM(CASE WHEN dod IS NOT NULL THEN 1 ELSE 0 END) AS con_dod, SUM(CASE WHEN deathtime IS NOT NULL THEN 1 ELSE 0 END) AS con_deathtime, SUM(CASE WHEN dod IS NOT NULL AND deathtime IS NULL THEN 1 ELSE 0 END) AS solo_data, SUM(CASE WHEN deathtime IS NOT NULL AND CAST(deathtime AS DATE) <> CAST(dod AS DATE) THEN 1 ELSE 0 END) AS date_discordanti FROM read_parquet('%s')", as_path(COH_OUT))))

cat("\n--- decessi entro 28 giorni: quanti guadagnano l'ora esatta\n")
print(ddb_get(sprintf("SELECT SUM(CASE WHEN deathtime IS NOT NULL THEN 1 ELSE 0 END) AS con_ora, SUM(CASE WHEN deathtime IS NULL AND dod IS NOT NULL THEN 1 ELSE 0 END) AS solo_data FROM read_parquet('%s') WHERE dod IS NOT NULL AND date_diff('day', CAST(time_zero AS DATE), CAST(dod AS DATE)) <= 28", as_path(COH_OUT))))

cat("\n--- controllo di segno: deathtime prima di time zero non deve esistere\n")
print(ddb_get(sprintf("SELECT SUM(CASE WHEN deathtime <= time_zero THEN 1 ELSE 0 END) AS deathtime_prima_di_t0, ROUND(MIN(date_diff('hour', time_zero, deathtime)), 1) AS ore_min FROM read_parquet('%s') WHERE deathtime IS NOT NULL", as_path(COH_OUT))))

cat("\nScritto:", COH_OUT, "\n")
cat("Eventuali deathtime <= time zero devono essere auditati; lo script 06 li tratta come mancanti e ricade su dod.\n")
cat("Poi in 06 basta che esista: viene preferito automaticamente a cohort_A_sev.parquet.\n")

## =============================================================================
## B. riconciliazione 6.370 contro 6.360
## =============================================================================
cat("\n########## B. audit della coorte eleggibile ##########\n")

## B1. identificativi della build corrente, per confronti futuri
ids <- ddb_get(sprintf("SELECT stay_id, subject_id, hadm_id, time_zero FROM read_parquet('%s') ORDER BY stay_id", as_path(COH_IN)))
saveRDS(ids, file.path(AUDIT, "eligible_ids_current.rds"))
write.csv(ids["stay_id"], file.path(AUDIT, "eligible_stay_ids_current.csv"), row.names = FALSE)
cat("Identificativi salvati:", nrow(ids), "stay in", file.path(AUDIT, "eligible_ids_current.rds"), "\n")

## B2. dove compare il numero 6360 nei file di output e negli script
cat("\n--- ricerca di 6360 nei file di testo di OUT_DIR e nella cartella degli script\n")
targets <- c(list.files(OUT_DIR, pattern = "\\.(csv|txt|log|md)$", recursive = TRUE, full.names = TRUE),
             list.files(getwd(), pattern = "\\.(csv|txt|log|md|R)$", full.names = TRUE))
hits <- Filter(function(f) {
    x <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    any(grepl("6360|6,360", x))
}, targets)
if (length(hits)) {
    for (h in hits) {
        cat("\n==", h, "\n")
        x <- readLines(h, warn = FALSE)
        cat(grep("6360|6,360", x, value = TRUE), sep = "\n")
    }
} else cat("Nessun file contiene 6360: il conteggio del vecchio audit non e' stato salvato.\n")

## B3. altri parquet di coorte presenti, con conteggi e date
cat("\n--- altri file di coorte presenti\n")
cfiles <- list.files(file.path(OUT_DIR, "initiation"), pattern = "\\.parquet$", full.names = TRUE)
if (length(cfiles)) {
    info <- do.call(rbind, lapply(cfiles, function(f) {
        n <- tryCatch(ddb_get(sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s')", as_path(f)))$n,
                      error = function(e) NA_integer_)
        data.frame(file = basename(f), righe = n, modificato = file.info(f)$mtime)
    }))
    print(info, row.names = FALSE)
}

## B4. ipotesi da testare: la vecchia definizione di PaO2 arterioso
##     (solo itemid, senza spec_type = 'ART.') dava un conteggio diverso
cat("\n--- conteggio di PaO2 arteriosi con le due definizioni\n")
LE_PAR <- ensure_subset(con, "labevents_bloodgas\\.parquet$", "labevents_bloodgas.parquet",
                        "labevents",
                        c("subject_id", "hadm_id", "specimen_id", "itemid", "charttime", "value", "valuenum"),
                        unique(c(ITEMS$pao2_le, 52033)))
ddb_exec(sprintf("CREATE VIEW le AS SELECT * FROM read_parquet('%s')", LE_PAR))
n_old <- ddb_get(sprintf("SELECT COUNT(*) AS n FROM le WHERE itemid IN (%s) AND valuenum BETWEEN 20 AND 700", sqlids(ITEMS$pao2_le)))$n
n_new <- ddb_get(sprintf("WITH spec AS (SELECT DISTINCT specimen_id, trim(value) AS t FROM le WHERE itemid = 52033) SELECT COUNT(*) AS n FROM le p JOIN spec s USING (specimen_id) WHERE p.itemid IN (%s) AND s.t = 'ART.' AND p.valuenum BETWEEN 20 AND 700", sqlids(ITEMS$pao2_le)))$n
cat("PaO2 senza filtro arterioso (vecchia definizione):", n_old, "\n")
cat("PaO2 con spec_type = 'ART.' (definizione attuale):", n_new, "\n")
cat("Se le due differiscono, la discrepanza 6.370/6.360 e' spiegata dal cambio di definizione.\n")

dbDisconnect(con, shutdown = TRUE)
cat("\nFatto.\n")
