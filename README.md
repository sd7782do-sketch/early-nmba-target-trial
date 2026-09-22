# R scripts

Run the scripts from the repository root in this order:

| Script | Input | Output |
|---|---|---|
| `03_initiation_counts.R` | Authorized MIMIC-IV tables | `initiation/cohort_A.parquet` and feasibility summaries |
| `05_baseline_severity.R` | `cohort_A.parquet` plus MIMIC-IV severity variables | `initiation/cohort_A_sev.parquet` |
| `08_deathtime_and_audit.R` | `cohort_A_sev.parquet` plus admissions | `initiation/cohort_A_sev_dt.parquet` and restricted local audit files |
| `07_hourly_grace.R` | Enriched cohort plus MIMIC-IV respiratory data | `initiation/hourly_grace.parquet` |
| `06_primary_analysis.R` | Enriched cohort, exact death times, and hourly grid | `primary/` estimates, diagnostics, tables, figures, and manifest |

```r
source("_paths.R")
source("R/03_initiation_counts.R")
source("R/05_baseline_severity.R")
source("R/08_deathtime_and_audit.R")
source("R/07_hourly_grace.R")
source("R/06_primary_analysis.R")
```

`01_config.R` freezes the item mapping and design constants; `_paths.R`
provides output-path and intermediate-file helpers. Script 06 automatically
detects `cohort_A_sev_dt.parquet` and `hourly_grace.parquet` and records their
use in `run_manifest.txt`.

Set `QUICK <- TRUE` in script 06 for a structural check before the full run.
See the root `README.md` for setup, data-governance requirements, outputs, and
the reproducibility boundary.
