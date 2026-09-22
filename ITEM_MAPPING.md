# Frozen MIMIC-IV item mapping

This mapping reproduces the manually reviewed selections used for the reported
run. Item labels should still be checked against `d_items` and `d_labitems` when
porting the code to a different MIMIC-IV version.

## Primary NMBA exposure (`inputevents`)

| itemid | Label |
|---:|---|
| 221555 | Cisatracurium |
| 229233 | Rocuronium |
| 222062 | Vecuronium |

The primary exposure additionally requires a documented continuous infusion
according to `INFUSION_RULE`. Peri-intubation `chartevents` items 227213,
227214, and 229788 are retained for documentation only and do not define the
primary exposure.

## Active sedation at time zero (`inputevents`)

| itemid | Label |
|---:|---|
| 222168 | Propofol |
| 221744 | Fentanyl |
| 225154 | Morphine Sulfate |
| 221833 | Hydromorphone (Dilaudid) |
| 225942 | Fentanyl (Concentrate) |
| 221668 | Midazolam (Versed) |
| 221385 | Lorazepam (Ativan) |
| 229420 | Dexmedetomidine (Precedex) |
| 225150 | Dexmedetomidine (Precedex) |
| 221712 | Ketamine |
| 225972 | Fentanyl (Push) |

Inclusion in this list does not by itself establish active sedation. The cohort
query requires a non-null infusion rate and an administration interval
overlapping time zero. Thus an item labelled as a push does not qualify unless
the source record satisfies those infusion/overlap criteria.

## Respiratory and laboratory definitions

| Role | Table | itemid | Label/definition |
|---|---|---:|---|
| Invasive ventilation interval | `procedureevents` | 225792 | Invasive Ventilation |
| PEEP | `chartevents` | 220339 | PEEP set |
| FiO2 | `chartevents` | 223835 | Inspired O2 Fraction |
| PaO2 | `labevents` | 50821 | pO2, blood gas |
| Blood-gas specimen type | `labevents` | 52033 | Specimen Type |
| Plateau pressure | `chartevents` | 224696 | Plateau Pressure |
| Tidal volume | `chartevents` | 224685 | Tidal Volume (observed) |

Arterial PaO2 requires item 52033 to equal `ART.` on the same `specimen_id`.

## Baseline severity items added by script 05

| Role | itemids |
|---|---|
| Vasopressors | 221906, 221289, 221749, 222315, 221662 |
| Inotropes | 221653, 221986 |
| Lactate | 50813 |
| Creatinine | 50912 |
| Bilirubin, total | 50885 |
| Platelet count | 51265 |
