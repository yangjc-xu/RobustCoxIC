# Data Dictionary

This dictionary documents the two UK Biobank–derived analysis data files used in
the data application (Section 7 of the manuscript). These files contain
individual-level UK Biobank data and are **not** included in the reproducibility
materials; access requires a formal application to the UK Biobank
(http://www.ukbiobank.ac.uk). The dictionary describes their structure so that an
approved researcher can reconstruct and interpret them. The analysis scripts
expect both files in the `data/` subdirectory of the reproducibility materials.

## Files

| File | Endpoint | Rows (participants) |
|------|----------|---------------------|
| `diabetes.tab`     | Age at onset of diabetes     | 465,565 |
| `hypertension.tab` | Age at onset of hypertension | 483,525 |

Both are tab-delimited plain-text files with a header row, read with
`data.table::fread()`. Each has **236 columns**: a participant identifier, the two
interval-censoring endpoints, 13 baseline covariates, and 220 genetic variants.

## Variables

| Variable | Type | Description | Coding / units |
|----------|------|-------------|----------------|
| `f.eid` | integer | Participant identifier. The analysis scripts replace it with a sequential de-identified index (`1:n`). | — |
| `left` | numeric | Left endpoint of the interval-censoring window: the largest age at which the participant is known to be event-free (`0` if none). | age in days (`diabetes.tab`) or years (`hypertension.tab`) |
| `right` | numeric | Right endpoint of the window: the smallest age at which the event is known to have occurred; `Inf` if the event was not observed (right-censored). | age in days (`diabetes.tab`) or years (`hypertension.tab`); `Inf` = right-censored |
| `sex` | binary | Biological sex. | `1` = male, `0` = female |
| `white` | binary | Race. | `1` = white, `0` = non-white |
| `England` | binary | Birthplace. | `1` = England, `0` = other |
| `pc1`–`pc10` | numeric | Top ten genetic principal components (population substructure). | continuous |
| 220 variant columns | numeric | Pre-screened genetic variants; one column per variant, named by its variant identifier (available via `colnames`). | additive genotype dosage (0/1/2 minor alleles) |

## Column order

The two files place the covariate and variant blocks in a different order; the
analysis code selects the covariates by name and the variants by these fixed
column positions (`variant_cols` in `application/functions.R`).

- **`diabetes.tab`:** column 1 = `f.eid`; columns 2–3 = `left`, `right`;
  columns 4–16 = `sex`, `white`, `England`, `pc1`–`pc10`;
  columns 17–236 = the 220 genetic variants.
- **`hypertension.tab`:** column 1 = `f.eid`; columns 2–3 = `left`, `right`;
  columns 4–223 = the 220 genetic variants;
  columns 224–236 = `sex`, `white`, `England`, `pc1`–`pc10`.

## Derivation

The endpoints were derived by converting ICD-10 diagnosis codes to PheWAS codes
and taking the difference between the date of diagnosis and the date of birth
(age at onset), yielding interval-censored ages; see Section 7 of the manuscript
and `application/fit_null.R`. With the `rounded` option, `application/fit_null.R`
rounds `left` and `right` down to whole years (after dividing by 365 for
`diabetes.tab`).
