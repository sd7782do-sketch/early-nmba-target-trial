# Release notes

## v1.0.0 — 2026-09-22

Initial public code release for the target trial emulation of early continuous
neuromuscular blockade in severe hypoxemic respiratory failure using MIMIC-IV
v3.1.

Included:

- portable configuration with the frozen, manually reviewed item mapping;
- cohort construction and initiation counts (`R/03_initiation_counts.R`);
- baseline severity enrichment (`R/05_baseline_severity.R`);
- exact death-time integration and cohort audit (`R/08_deathtime_and_audit.R`);
- hourly grace-period covariate extraction (`R/07_hourly_grace.R`);
- final primary clone-censor-weight analysis (`R/06_primary_analysis.R`);
- target trial protocol and statistical analysis plan;
- complete analysis decision and correction log; and
- citation, licensing, deposit, and data-governance documentation.

The public configuration removes personal paths and local inventory-index
dependencies while preserving the selected itemids and analytic definitions
used for the reported run.

This release contains no credentialed MIMIC-IV data, derived patient-level
files, credentials, or generated study outputs.
