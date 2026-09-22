# Analysis log

A record of the design decisions taken during this analysis, the errors found and corrected,
and what each change did to the estimate. It exists so that a reader can tell which choices
were prespecified, which were forced by the data, and which were corrections of mistakes.

Every run writes a `run_manifest.txt` with the full parameter set and `sessionInfo()`. The
estimates below are the 28-day risk ratio for early initiation against no early initiation,
under the overlap tilt.

---

## Summary across versions

| Version | Cohort | Principal change | RR (95% CI) |
|---|---|---|---|
| v1 | 5,128 | First complete run | 1.41 (1.26–1.59) |
| v2 | 5,128 | Censoring weight of the early clone corrected | 1.53 (1.37–1.70) |
| v3 | 5,128 | Range filtering of covariates; tidal volume added to the model | 1.47 (1.32–1.63) |
| v4 | 4,918 | One qualifying stay per patient | 1.45 (1.30–1.62) |
| v5 | 4,918 | Two propensity scores; deaths placed at end of day | 1.46 (1.31–1.64) |
| v7 | 4,918 | Time-varying censoring covariates | 1.38 (1.22–1.56) |
| v8 | 4,918 | Weight truncation at the 99th percentile (rejected) | 1.52 (1.36–1.69) |
| **v10** | **4,918** | **Exact death times; truncation moved to sensitivity** | **1.37 (1.22–1.55)** |

Sensitivity analyses in the final run span 1.37 to 1.53, except the analysis without the
overlap tilt (1.90), which addresses a different target population. No discretionary choice
reversed the direction or changed the order of magnitude. The reduction from 1.46 to 1.37
when time-varying censoring covariates were introduced is itself a finding: roughly a
quarter of the excess previously attributed to exposure was confounding by deterioration
occurring inside the grace period.

---

## v1 — first complete run

Design as in the feasibility protocol: exposure defined solely by initiation within a
12-hour grace period with no minimum duration; cohort restricted to patients with an active
sedative infusion at time zero; 28-day mortality; clone-censor-weight with overlap
weighting; multiple imputation with Rubin's rules; bootstrap re-estimating the weights.

## v2 — censoring weight of the early clone

**Error found.** The early clone received a weight of 1/e(X) across its whole follow-up,
including the grace period. It cannot be censored before the deadline, so its censoring
weight there must be 1. The same deaths in the first 12 hours therefore carried roughly
1/e(X) times more weight in the early arm than in the comparator, and the grace-period rows
dominated the early arm's hazard.

**Correction.** Weight 1 during the grace period, 1/e(X) afterwards. **Effect:** 1.41 → 1.53,
entirely through the exposed arm.

Two further problems were corrected in the same pass. Landmark estimates at 7 and 14 days had
been read from the interval whose midpoint was nearest the landmark, returning 6.5 and 13.5
days; grid boundaries were moved onto whole days. A regular expression intended to silence
the benign "non-integer successes" warning did not match, so genuine warnings could have been
hidden; warnings are now filtered correctly and any non-benign warning is written to a log.

## v3 — range filtering, model additions, checks

**Data problem found.** Physiologically impossible values in covariates entering the
propensity score: negative driving pressures, creatinine above 20 mg/dL, and 173
norepinephrine-equivalent values above 5 µg/kg/min, the largest 35.8, which is a unit error
rather than a dose. These are set to missing and imputed; the bounds are data-processing
rules and an analysis retaining the original values is reported as a sensitivity.

**Model addition.** Tidal volume added to the propensity score: it was the only covariate
above |SMD| 0.10 after weighting. It enters in millilitres because height is unavailable.

**Checks added.** A saturated weighted estimator alongside the parametric model, introduced
because the fitted curves overstated mortality in the first 12 hours about fivefold relative
to the crude count; the comparison showed the discrepancy confined to the first day.
Balance under the final clone-censor weights reported alongside baseline balance.

**Effect:** 1.53 → 1.47.

**Interpretation error corrected.** The risk difference at the end of the grace period had
been described as a falsification test, on the grounds that no deviation is observable before
the deadline. This is false and contradicted the code: the no-early clone is censored at the
moment of initiation, usually well before the deadline, so the strategies diverge inside the
window and the contrast is not structurally null. The quantity is retained as descriptive and
labelled as such in the script, the manifest, and the manuscript.

## v4 — unit of analysis

**Problem found.** The cohort file contained 6,370 qualifying stays for 6,112 patients, with
73 patients contributing two stays within a single hospitalization. The bootstrap therefore
resampled stays while the Methods claimed patients; the same death could enter twice; and a
patient treated during an earlier qualifying stay re-entered as a new user, because prevalent
use was assessed only at each stay's own time zero.

**Correction.** Restriction to the first qualifying stay per patient, applied to the raw
cohort **before** any exclusion. Applying it afterwards would select the second episode
conditionally on exclusion of the first, which is informative selection.

**Effect:** cohort 5,128 → 4,918; RR 1.47 → 1.45. The conclusion was unaffected, which is the
useful finding: the correction was necessary for correctness, not because it changed the
answer.

**Descriptive analysis added.** Baseline severity and observed mortality across strata of
time to first infusion, introduced to distinguish severity present at time zero from severity
developing afterwards. It produced the most informative single result of the study: measured
baseline severity *decreases* as the delay to initiation increases, while observed mortality
does not.

## v5 — two propensity scores, death-time convention

**Problem found (introduced in v4).** The propensity score had been fitted only among
patients at risk at the deadline, correctly reasoning that patients who died earlier had no
opportunity to initiate. But the same score was also used for the overlap tilt, and the exact
balance property of overlap weighting holds only when the score is fitted on the sample to
which the weights apply. Baseline balance for lactate degraded to 0.109 with the correction
overshooting in sign.

**Correction.** Two scores with distinct roles: the tilt fitted in the full analytic cohort,
the adherence probability in the risk set at the deadline. A consequence to state in any
write-up: the early clone's final weight is e(X)[1 − e(X)] divided by the estimated adherence
probability, not 1 − e(X) exactly. **Effect:** maximum baseline |SMD| 0.109 → 0.018.

**Second problem found (also from v4).** Deaths known only as a date had been placed at
midday. For a patient whose time zero falls in the afternoon this places death before time
zero: 41 events had to be anchored to the first interval, and deaths inside the grace period
rose from 44 to 179. Placement at the end of the day keeps all event times positive.

**Effect:** 1.45 → 1.46, with placement sensitivities spanning 1.46 to 1.48.

## v6–v7 — time-varying censoring weights

**Rationale.** The artificial censoring of the no-early clone is triggered by the start of an
infusion, which typically follows deterioration occurring inside the grace period. Censoring
weights built only on baseline covariates cannot account for it, and the comparison across
initiation windows suggested such a determinant exists.

An hourly grid of PaO₂/FiO₂, PEEP, plateau pressure and vasopressor dose was extracted for
the first 24 hours (script 07), reusing the cohort-construction definitions rather than
restating them. Norepinephrine equivalents were excluded from the censoring model: 5.2% of
hourly values exceed 5 µg/kg/min, with a 99th percentile of 23 and a maximum of 749 — the
same unit problem as at baseline, amplified by summing overlapping infusion segments.

Two coding faults were found and fixed. The censoring model copied only `PS_VARS` into its
design frame, but plateau pressure is deliberately not in `PS_VARS`, so the column was NULL.
And with plateau missing in 23% of patients, `glm` would have dropped rows and returned a
shorter `fitted()` vector than the person-time frame, misaligning the weights **without an
error**; the fallback chain (hourly value, then baseline, then driving pressure + PEEP) and an
explicit length check now prevent this.

**Effect:** 1.46 → 1.38. The exposed arm's risk was unchanged (0.5007) while the comparator's
rose from 0.343 to 0.363 — the direction expected if the censoring had been informative.

**Cost.** The comparator arm's effective sample size fell from 1,600 to 816, with a few
weights above 100, and balance under the final weights rose from 0.038 to 0.080.

## v8 — weight truncation, tested and rejected

**Attempt.** Truncate censoring weights at the 99th percentile within arm, to tame the tail.

**Result.** The cap in the comparator arm fell to 1.614. The distribution of the censoring
weight is concentrated near 1, so the 1% above the cap carried the entire reweighting for
informative censoring. Truncating it made the weights nearly uniform: the comparator risk fell
to 0.327, *below* the value obtained without time-varying covariates at all, the risk ratio
rose to 1.52, and balance under the final weights degraded from 0.080 to 0.171 — worse than
any other version.

**Decision.** Rejected as the primary specification. Balance is the arbiter, not effective
sample size: a high effective sample size with nearly uniform weights is a symptom, not a
virtue. Truncation is retained as a sensitivity analysis, where it is reported for
transparency.

## v10 — exact death times

**Change.** `admissions.deathtime` joined into the cohort by `subject_id` (the field is
populated only for the admission in which death occurred, which need not be the index one).
Of deaths within 28 days of time zero, 94.6% gained an exact recorded time. Across the
cohort, 628 deaths remain date-only. The two sources agree on the calendar date in every
case. Thirteen records in the enriched pre-restriction cohort—and six in the final analytic
cohort—had a recorded time at or before time zero and were treated as having no exact time.

**Effect:** 1.38 → 1.37, and three concrete improvements. The sensitivity to within-day
placement of deaths collapsed from a 1.46–1.48 span to 1.37 in every variant, so that
limitation leaves the manuscript. The maximum calibration discrepancy beyond 24 hours fell
from 2.2 to 1.3 percentage points. No event required anchoring to the first interval.

**Final specification:** no weight truncation; time-varying covariates PaO₂/FiO₂, PEEP and
plateau pressure; exact death times where available. RR 1.374 (1.215–1.554), risk difference
13.6 percentage points (8.3–18.9).

---

## Decisions taken and not revisited

- **Physiological rather than Berlin definition of the population.** The radiological
  criterion cannot be harmonised across the databases intended for a transportability
  analysis, so eligibility rests on oxygenation and PEEP. Stated as a limitation.
- **No minimum infusion duration in the exposure definition.** Requiring persistence would
  condition on information unavailable at initiation.
- **Restriction to patients with an active sedative infusion at time zero.** Neuromuscular
  blockade without concurrent sedation is not a clinically admissible strategy. This narrows
  the population and does not establish comparable sedation depth.
- **Plateau pressure not entered into the propensity score.** It is the arithmetic sum of
  driving pressure and PEEP, both included. It is reported as a balance variable, and for
  that reason is consistently among the worst-balanced.
- **No censoring at hospital discharge.** The mortality field covers deaths after discharge.
- **Overlap weighting as the primary estimand.** With 8.8% exposed and markedly asymmetric
  propensity score distributions, the overlap population is the one the data support. The
  analysis without the tilt addresses a different target population, not a stronger version
  of the same result.
- **No instrumental variable, no negative control outcome.** Considered at design and dropped
  for lack of a credible instrument and a suitable control outcome.
- **eICU-CRD transportability deferred.** The exposure and outcome definitions do not
  transfer without substantial rework (offset times, free-text drug records, no
  post-discharge mortality), so it belongs to separate work rather than to this analysis.

## What this analysis does not establish

The estimate is discordant with the randomized evidence in both direction and magnitude. The
internal evidence — measured baseline severity decreasing with delay to initiation while
mortality does not — points to a determinant arising after time zero rather than to a
pharmacological effect, but cannot distinguish deterioration that prompts treatment from harm
caused by it. Nothing in this repository should be read as evidence that early neuromuscular
blockade increases mortality.
