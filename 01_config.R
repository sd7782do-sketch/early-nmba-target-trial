## =============================================================================
## 01_config.R
## Portable public configuration for MIMIC-IV v3.1.
##
## Before running, set MIMIC_DIR to the local MIMIC-IV root. The standard
## hosp/ and icu/ layout is detected for Parquet, CSV.GZ, or CSV files.
## Optionally set NMBA_OUT_DIR; otherwise outputs go to ./feasibility_out.
## =============================================================================

CONFIG_REVISION <- "public_v1.0.0"
MIMIC_VERSION   <- "3.1"

## ---- directories ------------------------------------------------------------
if (!exists("MIMIC_DIR")) MIMIC_DIR <- Sys.getenv("MIMIC_DIR", unset = "")
if (!nzchar(MIMIC_DIR) || !dir.exists(MIMIC_DIR)) {
  stop(
    "MIMIC_DIR is not set or does not exist. Set it before running, for example:\n",
    "  Sys.setenv(MIMIC_DIR = 'D:/data/mimiciv/3.1')"
  )
}
MIMIC_DIR <- normalizePath(MIMIC_DIR, winslash = "/", mustWork = TRUE)

if (!exists("OUT_DIR")) {
  out_env <- Sys.getenv("NMBA_OUT_DIR", unset = "")
  OUT_DIR <- if (nzchar(out_env)) out_env else "feasibility_out"
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_DIR <- normalizePath(OUT_DIR, winslash = "/", mustWork = TRUE)

## ---- source-table resolution ------------------------------------------------
TABLE_MODULE <- c(
  admissions      = "hosp",
  patients        = "hosp",
  labevents       = "hosp",
  d_labitems      = "hosp",
  icustays        = "icu",
  procedureevents = "icu",
  inputevents     = "icu",
  chartevents     = "icu",
  d_items         = "icu"
)

sql_quote <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")

resolve_table_path <- function(name) {
  name <- tolower(name)
  module <- unname(TABLE_MODULE[[name]])
  if (is.null(module)) stop("Unknown MIMIC-IV table: ", name)

  filenames <- paste0(name, c(".parquet", ".csv.gz", ".csv"))
  roots <- unique(c(file.path(MIMIC_DIR, module), MIMIC_DIR))
  direct <- unlist(lapply(roots, function(root) file.path(root, filenames)),
                   use.names = FALSE)
  hits <- direct[file.exists(direct)]

  ## Also support a partitioned Parquet directory: <root>/<module>/<table>/*.parquet
  part_dir <- file.path(MIMIC_DIR, module, name)
  if (!length(hits) && dir.exists(part_dir) &&
      length(list.files(part_dir, pattern = "\\.parquet$", ignore.case = TRUE))) {
    return(paste0(normalizePath(part_dir, winslash = "/", mustWork = TRUE),
                  "/*.parquet"))
  }

  if (!length(hits)) {
    pattern <- paste0("^", name, "\\.(parquet|csv|csv\\.gz)$")
    hits <- list.files(MIMIC_DIR, pattern = pattern, recursive = TRUE,
                       full.names = TRUE, ignore.case = TRUE)
    if (length(hits)) {
      norm <- gsub("\\\\", "/", hits)
      in_module <- grepl(paste0("/", module, "/"), norm, fixed = TRUE)
      if (any(in_module)) hits <- hits[in_module]
    }
  }

  hits <- unique(hits)
  if (!length(hits)) {
    stop("Table not found under MIMIC_DIR: ", module, "/", name,
         ". Expected .parquet, .csv.gz, or .csv.")
  }
  if (length(hits) > 1L) {
    stop("Ambiguous table ", name, "; candidates: ",
         paste(normalizePath(hits, winslash = "/", mustWork = TRUE),
               collapse = "; "))
  }
  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}

src <- function(name, module = NULL) {
  ## `module` is accepted for compatibility; the public mapping above is fixed.
  p <- resolve_table_path(name)
  if (grepl("\\.parquet$", p, ignore.case = TRUE)) {
    return(sprintf("read_parquet(%s)", sql_quote(p)))
  }
  sprintf("read_csv_auto(%s, header=true, sample_size=50000)", sql_quote(p))
}

sqlids <- function(x) {
  x <- suppressWarnings(as.integer(as.character(x)))
  x <- x[!is.na(x)]
  if (!length(x)) "NULL" else paste(unique(x), collapse = ",")
}

## ---- confirmed item mapping -------------------------------------------------
NMBA_ITEM_MAP <- data.frame(
  itemid = c(221555L, 229233L, 222062L, 227213L, 227214L, 229788L),
  molecule = c("cisatracurium", "rocuronium", "vecuronium",
               "vecuronium", "cisatracurium", "rocuronium"),
  source_table = c(rep("inputevents", 3), rep("chartevents", 3)),
  exposure_role = c(rep("primary", 3), rep("intubation_secondary", 3)),
  label = c("Cisatracurium", "Rocuronium", "Vecuronium",
            "Vecuronium (Intubation)", "Cis-atracurium (Intubation)",
            "Rocuronium (Intubation)"),
  stringsAsFactors = FALSE
)

NMBA_INPUT_ITEMIDS <- NMBA_ITEM_MAP$itemid[
  NMBA_ITEM_MAP$source_table == "inputevents"
]
NMBA_INTUBATION_ITEMIDS <- NMBA_ITEM_MAP$itemid[
  NMBA_ITEM_MAP$source_table == "chartevents"
]

## These 11 inputevents itemids are frozen from the manually reviewed inventory
## used for the reported run. The cohort query additionally requires a non-null
## rate and an infusion overlapping time zero.
SEDATION_INPUT_ITEMIDS <- c(
  222168L, # Propofol
  221744L, # Fentanyl
  225154L, # Morphine Sulfate
  221833L, # Hydromorphone (Dilaudid)
  225942L, # Fentanyl (Concentrate)
  221668L, # Midazolam (Versed)
  221385L, # Lorazepam (Ativan)
  229420L, # Dexmedetomidine (Precedex)
  225150L, # Dexmedetomidine (Precedex)
  221712L, # Ketamine
  225972L  # Fentanyl (Push); rate/overlap criteria still apply
)

ITEMS <- list(
  nmba_input      = NMBA_INPUT_ITEMIDS,
  nmba_intubation = NMBA_INTUBATION_ITEMIDS,
  sedation_input  = SEDATION_INPUT_ITEMIDS,
  vent_proc       = c(225792L), # Invasive Ventilation interval
  peep            = c(220339L), # PEEP set
  fio2_ce         = c(223835L), # Inspired O2 Fraction
  fio2_le         = NULL,
  pao2_le         = c(50821L),
  specimen_type_le = c(52033L),
  plateau         = c(224696L),
  vt              = c(224685L),
  height          = NULL,
  weight          = NULL
)

## ---- design parameters ------------------------------------------------------
PARAMS <- list(
  age_min = 18,
  pf_threshold = 150,
  peep_min = 5,
  pair_window_h = 2,
  intubation_to_t0_h = 48,
  elig_window_h = 48,
  vent_gap_h = 14,
  inf_gap_h = 1,
  grace_h = c(6, 12, 24),
  fu_days = 28,
  pao2_provenance_mode = "specimen_arterial_only",
  arterial_specimen_values = c("ART."),
  primary_grace_h = 12,
  use_intubation_chartevents_as_primary = FALSE,
  sedation_window = "pre_time_zero_only",
  rass_window = "pre_time_zero_only",
  rass_during_nmba = "not_interpretable"
)

TIME_ZERO_DEFINITION <- paste(
  "first documented arterial (52033 = ART.) pO2/FiO2 <150 with PEEP >=5",
  "within 48 hours after invasive ventilation starts"
)

## Continuous infusions only; bolus rules are retained for transparency but are
## not used in the primary exposure definition.
INFUSION_RULE <- paste(
  "rate IS NOT NULL",
  "AND (lower(trim(ordercategorydescription)) = 'continuous med'",
  "     OR lower(trim(ordercategoryname)) = '01-drips')",
  "AND (statusdescription IS NULL",
  "     OR lower(trim(statusdescription)) <> 'rewritten')",
  sep = "\n"
)

BOLUS_RULE <- paste(
  "rate IS NULL",
  "AND (lower(trim(ordercategorydescription)) = 'drug push'",
  "     OR lower(trim(ordercategoryname)) = '05-med bolus')",
  sep = "\n"
)

write.csv(NMBA_ITEM_MAP, file.path(OUT_DIR, "nmba_item_map_confirmed.csv"),
          row.names = FALSE)

cat("Config loaded.\n")
cat("MIMIC-IV version:", MIMIC_VERSION, "\n")
cat("MIMIC_DIR:", MIMIC_DIR, "\n")
cat("OUT_DIR:", OUT_DIR, "\n")
cat("Primary NMBA itemids:", paste(NMBA_INPUT_ITEMIDS, collapse = ", "), "\n")
cat("Sedation itemids:", paste(SEDATION_INPUT_ITEMIDS, collapse = ", "), "\n")
cat("Time zero:", TIME_ZERO_DEFINITION, "\n")
