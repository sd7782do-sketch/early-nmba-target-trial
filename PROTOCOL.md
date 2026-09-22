# Target trial protocol and statistical analysis plan

Early continuous neuromuscular blockade in severe hypoxemic respiratory failure:
a target trial emulation in MIMIC-IV.

This document states the trial the analysis emulates, the emulation, and the analysis
plan, in the form the code implements. Where a decision was taken after the data were
first examined, `ANALYSIS_LOG.md` records when and why; nothing here is presented as
having been prespecified when it was not.

---

## 1. The target trial

| Component | Specification |
|---|---|
| **Eligibility** | Adults receiving invasive mechanical ventilation with a documented arterial PaO₂/FiO₂ below 150 mm Hg and PEEP of at least 5 cm H₂O, occurring within 48 hours of the start of invasive ventilation; an active sedative infusion at eligibility; no continuous neuromuscular blocking agent (NMBA) infusion already running at eligibility. |
| **Treatment strategies** | (1) Start a continuous NMBA infusion within 12 hours of eligibility. (2) Do not start a continuous NMBA infusion within that window. No constraint on agent, dose, or duration; treatment after the window is unrestricted in both arms. |
| **Assignment** | Random, with the strategies revealed at time zero. |
| **Time zero** | Randomization, coinciding with the moment eligibility is met. |
| **Outcome** | Death from any cause within 28 days of time zero. |
| **Causal contrast** | Adherence to the assigned initiation strategy during the grace period. Subsequent treatment is not constrained, so this is a per-protocol contrast limited to the initiation decision, not an intention-to-treat effect of assignment. |
| **Analysis** | Comparison of 28-day mortality risks under the two strategies. |

## 2. The emulation

**Data source.** MIMIC-IV, single centre, deidentified, credentialed access. Patient anchor
year groups 2008–2010 through 2020–2022. Dates are shifted per patient, so within-patient
intervals are preserved but calendar dates are not interpretable.

**Unit of analysis.** The patient. When a patient has more than one intensive care stay
meeting eligibility, only the first is retained, and the restriction is applied to the raw
eligible set **before** any exclusion, so that later stays are not selected conditionally on
exclusion of the first.

**Time zero.** The first time point within 48 hours of the start of invasive ventilation at
which an arterial PaO₂/FiO₂ below 150 mm Hg coincides with PEEP of at least 5 cm H₂O.
Arterial samples are identified through the specimen record (`specimen_id` linked to the
blood-gas specimen type), not by analyte alone. PaO₂ is paired with the most recent
preceding FiO₂ and PEEP within a two-hour window.

**Exposure.** Start of a continuous NMBA infusion (cisatracurium, rocuronium, or vecuronium
recorded as infusions) after time zero and within the grace period. Boluses given around
intubation are recorded separately in the source data and are not treated as exposure. No
minimum duration is required, because conditioning on treatment persistence would use
information unavailable at the moment of initiation.

**Comparator.** No continuous NMBA infusion started within the same window. This includes
patients who start later and patients who never start.

**Outcome.** Death from any cause within 28 days of time zero, from the recorded time of
death where available and otherwise from the recorded date of death, which also captures
deaths after hospital discharge. Hospital discharge is not a censoring event. Deaths known
only as a date are placed at the end of the recorded day, a convention that keeps event
times positive; alternative placements are examined in sensitivity analyses.

**Assignment.** Not observed. Each patient is cloned into both strategy arms at time zero,
clones are censored when their observed treatment becomes incompatible with their assigned
strategy, and the resulting selection is addressed by inverse probability of censoring
weighting (section 4).

## 3. Covariates

Entered into the propensity score model, measured at or before time zero: age, sex, anchor
year group, PaO₂/FiO₂ ratio, PEEP, inspired oxygen fraction, driving pressure, tidal volume
(in millilitres; height is unavailable, so indexing to predicted body weight is not
possible), active vasopressor infusion and its dose in norepinephrine equivalents, active
inotrope infusion, lactate, creatinine, platelet count.

Assessed for balance but not entered: plateau pressure (the arithmetic sum of driving
pressure and PEEP, both included) and PaO₂ (a component of the oxygenation ratio).

**Range filtering.** Values outside the analysis bounds below are set to missing before
imputation and handled by multiple imputation. These are data-processing rules, not claims
about physiological impossibility; an analysis retaining the original values is reported as
a sensitivity.

| Variable | Retained range |
|---|---|
| PaO₂/FiO₂ | 20 to 600 mm Hg |
| PEEP | 0 to 30 cm H₂O |
| FiO₂ | 21 to 100% |
| PaO₂ | 20 to 600 mm Hg |
| Plateau pressure | 5 to 60 cm H₂O |
| Driving pressure | 1 to 45 cm H₂O |
| Tidal volume | 100 to 1,200 mL |
| Lactate | 0.2 to 30 mmol/L |
| Creatinine | 0.1 to 20 mg/dL |
| Platelets | 1 to 1,500 ×10³/µL |
| Norepinephrine equivalents | 0 to 5 µg/kg/min |
| Age | 18 to 100 years |

A vasopressor dose of zero in a patient without a vasopressor infusion is a structural
zero, not a missing value.

**Missing data.** Multiple imputation by chained equations, ten datasets, predictive mean
matching. The observed initiation indicator and the outcome are predictors and are not
imputed. Estimates are combined by Rubin's rules, ratios on the log scale.

## 4. Analysis

**Cloning and censoring.** Each patient enters both arms at time zero. The no-early clone is
censored at the moment a continuous NMBA infusion starts, if that occurs within the grace
period. The early clone is censored at the end of the grace period if no infusion has
started. A death contributes to both clones only if neither has deviated before it.

**Two propensity score models, with distinct roles.**

- *Tilt*: P(initiation within the grace period | X), fitted in the full analytic cohort.
  This defines the target population and, fitted by logistic regression on the sample to
  which the weights apply, preserves the exact balance property of overlap weighting.
- *Adherence*: the same specification, fitted among patients still at risk at the grace
  period deadline. This is the denominator of the early clone's censoring weight; patients
  who died earlier had no opportunity to initiate and do not belong to it.

**Weights.** For the early clone, the censoring weight is 1 before the deadline — no
censoring is possible there — and the reciprocal of the estimated adherence probability
afterwards. For the no-early clone, it is the reciprocal of the cumulative probability of
remaining uninitiated across the hourly intervals of the grace period, held constant
afterwards; the hourly initiation hazard is modelled on baseline covariates plus the
time-varying PaO₂/FiO₂, PEEP and plateau pressure of the current hour, because the artificial
censoring is triggered by deterioration occurring inside the window. The baseline tilt and
adherence probabilities are bounded to 0.02 and 0.98; the hourly initiation hazard is capped
at 0.98, and cumulative adherence probabilities are floored at 0.02 before inversion.

Censoring weights are multiplied by a baseline overlap tilt, e(X)[1 − e(X)], applied
identically to both clones. Because tilt and adherence are estimated on different samples,
the early clone's final weight is e(X)[1 − e(X)] divided by the estimated adherence
probability rather than 1 − e(X) exactly.

The censoring weights are **not** truncated in the primary analysis. Truncation at the 99th
percentile within arm is reported as a sensitivity; it degrades balance under the final
weights and is not preferred (see the v8 section of `ANALYSIS_LOG.md`).

**Estimand.** The difference in 28-day risk under the two initiation strategies in the
overlap-targeted population. Identification requires consistency, positivity, and
conditional exchangeability for the artificial censoring.

**Outcome model.** Follow-up is divided into hourly intervals within the grace period and
daily intervals thereafter, ending at exactly 28 days. Discrete-time hazards are modelled by
weighted pooled binomial regression with a complementary log-log link and a log
interval-width offset, including strategy arm, a natural cubic spline of time with three
degrees of freedom, and their interaction. Cumulative risks follow from the fitted hazards.

**Uncertainty.** Two hundred nonparametric bootstrap resamples at the patient level within
each imputed dataset, drawn before cloning, with propensity scores, censoring weights and
outcome model re-estimated in each resample. Within-imputation bootstrap variance and
between-imputation variance are combined by Rubin's rules. Curve intervals are pointwise,
not simultaneous.

**Prespecified checks.**

1. *Calibration.* The parametric curves are compared with a saturated weighted estimator
   fitting a separate hazard in each interval and arm. Both use the same event times, so
   agreement says nothing about uncertainty in death timing.
2. *End of grace period.* The risk difference at the deadline is reported descriptively. It
   is **not** a falsification test: the no-early clone is censored at the moment of
   initiation, so the strategies diverge inside the window and the contrast is not
   structurally null.
3. *Balance under the final weights*, in the post-grace risk set, alongside baseline balance.
4. *Severity by initiation window.* Baseline characteristics and observed mortality across
   strata of time to first infusion, to separate severity present at time zero from severity
   developing afterwards.

**Sensitivity analyses.** Grace periods of 6 and 24 hours; no range filtering of covariates;
truncation of censoring weights at the 99th percentile; alternative within-day placement of
date-only deaths (12 and 18 hours); omission of the overlap tilt, which targets the full
eligible cohort rather than the overlap population. E-values quantify the strength of
unmeasured confounding that would explain away the point estimate or move the interval to
the null.

## 5. What this design does not do

It does not estimate an intention-to-treat effect, because assignment is unobserved. It does
not identify a causal effect without the assumptions above, and the censoring model relies on
covariates measured at time zero and hourly within the grace period only. It cannot separate
harm from blockade from harm from the deterioration that prompts it; the analysis of severity
by initiation window locates the problem but does not resolve it.
