# GitHub and Figshare deposit guide

Use the same version label (`v1.0.0`) and the same archive on both platforms.
Do not upload MIMIC-IV tables, derived patient-level files, local configuration,
credentials, inventory paths, or generated audit identifiers.

## GitHub

Create a public repository named:

```text
early-nmba-target-trial
```

Suggested description:

> R code for a MIMIC-IV target trial emulation of continuous neuromuscular
> blockade initiated within 12 hours of physiological eligibility.

Create the repository without auto-generating a README, `.gitignore`, or license,
because all three are included. Upload the **contents** of the
`early-nmba-target-trial` directory so that `README.md` is at the repository root.

Suggested first commit message:

```text
Initial public analysis code release
```

After the files are visible, create a release:

- tag: `v1.0.0`
- release title: `v1.0.0 — Initial public release`
- release notes: copy the `v1.0.0` section of `RELEASE_NOTES.md`
- asset: attach the unchanged distribution ZIP

## Figshare

Create a **new item**, not a Project. Use the following metadata.

**Item type**

```text
Software
```

**Title**

```text
Early continuous neuromuscular blockade in severe hypoxemic respiratory failure: target trial emulation code
```

**Authors**

```text
Daniele Orso
```

Add an ORCID only if it is verified in the author's Figshare/ORCID profile.

**Description**

> R code and design documentation for a target trial emulation in MIMIC-IV
> v3.1 comparing initiation of a continuous neuromuscular blocking agent
> infusion within 12 hours of physiological eligibility with no initiation
> during that window. The analysis uses cloning, artificial censoring, inverse
> probability of censoring weighting with time-varying respiratory covariates,
> a baseline overlap tilt, multiple imputation, and a patient-level bootstrap.
> The archive includes cohort construction, baseline severity enrichment, exact
> death-time integration, hourly grace-period covariate extraction, the final
> primary analysis, the target trial protocol, and a complete analysis log. No
> MIMIC-IV source or patient-level derived data are distributed.

**Keywords**

```text
target trial emulation
neuromuscular blockade
acute hypoxemic respiratory failure
clone-censor-weight
overlap weighting
MIMIC-IV
critical care
```

**License**

```text
MIT
```

**Related material**

```text
https://github.com/sd7782do-sketch/early-nmba-target-trial
```

Choose the relationship identifying the GitHub page as the source-code
repository or related software, according to the labels shown by the interface.

**File**

Upload the unchanged `early-nmba-target-trial_v1.0.0.zip` archive.

Before publishing, reserve the DOI and add it to the manuscript's code
availability statement. Publish only after the author list, title, license,
file, description, and GitHub URL have been checked: later changes to files,
title, or authors may create a new Figshare version.

## After the Figshare DOI is reserved

The initial archive is valid without a self-referential DOI. After reservation,
add the DOI to the GitHub README and `CITATION.cff` in a metadata-only follow-up
commit. Do not replace the deposited v1.0.0 ZIP solely to insert its own DOI
unless a self-contained archive is required, because changing the deposited
file may create another Figshare version.

Suggested manuscript statement:

> Analysis code and design documentation are available on GitHub
> (https://github.com/sd7782do-sketch/early-nmba-target-trial) and are archived
> on Figshare (DOI: [insert reserved DOI]). MIMIC-IV is available to credentialed
> users through PhysioNet and cannot be redistributed.
