# Reproducibility Materials

**Robust Inference for the Cox Proportional Hazards Model with Interval-Censored Data**

This repository contains the R and C++ code needed to reproduce all numerical
results in the manuscript and the supplementary materials: the simulation
studies (Tables 1, 2, S1, S2, S3, and Figure S1) and the UK Biobank
data application (Figure 1, Figures S2 and S3).

---

## 1. Software requirements

* R version 4.5.3 with a working C++ toolchain (the C++ routines are compiled on the
  fly through `Rcpp` / `RcppArmadillo`).
* R packages (versions used in parentheses):
  `Rcpp` (1.1.1-1.1), `RcppArmadillo` (15.2.4-1), `tidyverse` (2.0.0),
  `dplyr` (1.2.0), `MASS` (7.3-65), `data.table` (1.18.2.1), `ggpubr` (0.6.3),
  `ggplot2` (4.0.2), `latex2exp` (0.9.6).

```r
install.packages(c("Rcpp","RcppArmadillo","tidyverse","dplyr","MASS",
                   "data.table","ggpubr","ggplot2","latex2exp"))
```

The C++ files (`src/*.cpp`) are compiled automatically by `Rcpp::sourceCpp()`
inside the scripts; they do not need to be compiled manually.

---

## 2. Directory structure

```
reproducibility_materials/
├── README.md
├── DATA_DICTIONARY.md
├── src/                             # Code shared by all scripts
│   ├── unireg_indep_improved.cpp    # NPMLE (EM) + robust/model-based variance
│   ├── unireg_indep_UKB.cpp         # Memory-efficient NPMLE for n = 400,000
│   ├── score_test_improved.cpp      # Score test with a time-independent covariate
│   ├── score_test_dep_improved.cpp  # Score test with a time-dependent covariate
│   └── utils.R                      # Helper function used by all R scripts
├── data/                            # Place the UK Biobank data files here
├── figures/                         # Reference copies of the published figures
├── simulation/                      # Simulation studies
│   ├── functions.R                  # Data generation, model fitting, summaries
│   ├── run_table1.R                 # Table 1  (type I error, n = 200, 400)
│   ├── run_tableS1.R                # Table S1 (type I error, n = 400,000)
│   ├── run_table2.R                 # Tables 2, S2, S3 (parameter estimation)
│   ├── plot_m.R                     # Figure S1: m-curve (supplement, Section S.11)
│   ├── summarize_simulation.R       # Example: how to read the saved outputs
│   └── output/                      # Simulation outputs are written here
└── application/                     # UK Biobank data application
    ├── functions.R                  # Endpoint settings and plotting function
    ├── fit_null.R                   # Fit the null Cox models
    ├── score_test.R                 # Score tests for the 220 genetic variants
    ├── gof_test.R                   # Goodness-of-fit residual processes
    ├── Figure1.R                    # Figure 1
    ├── FigureS2.R / FigureS3.R      # Figures S2 and S3
    └── output/                      # Application outputs are written here
```

Each implementation appears in a single place: the scripts `run_*.R`,
`fit_null.R`, `score_test.R`, and `gof_test.R` take the inputs that vary (true
model, sample size, endpoint, covariate) as command-line arguments and call the
shared functions.

---

## 3. Working directory and paths

All paths in the scripts are relative to the root directory of the unzipped
materials (the folder containing this README), and no path needs to be edited.
Run every script with `Rscript` from the root directory, e.g.,

```
cd reproducibility_materials
Rscript simulation/run_table1.R 1 200
```

The inputs that vary (true model, sample size, endpoint, covariate) are passed
as command-line arguments, as described below.

---

## 4. Data availability

The application uses individual-level UK Biobank data, which **cannot be
redistributed** and are therefore **not included** here. Access requires a
formal application to the UK Biobank (https://www.ukbiobank.ac.uk). The
application scripts expect two derived data files in the `data/` subdirectory:

* `data/diabetes.tab`   — columns `left`, `right` (interval endpoints, in days) plus
  the covariates, for the diabetes endpoint;
* `data/hypertension.tab` — the same for the hypertension endpoint, with the
  interval endpoints in years.

These were constructed as described in Section 7 of the manuscript (converting
ICD-10 codes to PheWAS codes and taking the age at diagnosis); see
`DATA_DICTIONARY.md`. No simulation step requires any external data.

---

## 5. Reproducing the simulation studies

All simulation data are generated internally, so no external files are needed.
Each script begins with a random seed drawn at run time (saved with the
output), so the reported numbers are Monte Carlo averages over 10,000
replicates; small deviations (in the last reported digit) from the published
values are expected. The last, optional argument of each script is the number
of replicates, which can be reduced for a quick check.

### 5.1 Table 1 — type I error rates, n = 200 and 400

```
Rscript simulation/run_table1.R [model] [n] [replicates]
```

`model` is 1–12 and `n` is the sample size; without arguments, all twelve
models are run at n = 200 and 400. Each run saves
`simulation/output/table1_model{model}_n{n}.RData`.

The true models (matching the rows of Table 1) are:

| Model | True model |
|-------|------------|
| 1  | λ(t)=(t/8)·exp(0.5X₂+0.5X₃) |
| 2  | λ(t)=(t/8)·exp(0.8X₁²) |
| 3  | λ(t)=(t/8)·exp(0.2X₁⁴+0.5X₂) |
| 4  | λ(t)=(t/8)·exp(0.8X₁²+0.5X₂+0.5X₃) |
| 5  | λ(t)=(t/8)·(1+0.5X₂²) |
| 6  | λ(t)=(t/8)·exp(X₁²)(1+0.5X₂²) |
| 7  | λ(t)=(t/8)·(2+X₂)² |
| 8  | λ(t)=(t/8)·{1+0.5·exp(X₁²)+X₂}² |
| 9  | T = exp(0.5X₂) + φ |
| 10 | T = exp(0.8X₁² + 0.5X₂) + φ |
| 11 | log T = 0.5X₂ + ε |
| 12 | log T = 0.8X₁² + 0.5X₂ + ε |

In each saved object, `summary_stats_wald` and `summary_stats_score` report the
empirical type I error: `alpha` is the **robust** test (sandwich variance) and
`alpha1`/`alpha2` are the **model-based** tests.

### 5.2 Table S1 — type I error rates, n = 400,000

```
Rscript simulation/run_tableS1.R [model] [replicates]
```

`model` is 1–8 (the first eight true models above), evaluated at n = 400,000
using the memory-efficient `src/unireg_indep_UKB.cpp`. Each run saves
`simulation/output/tableS1_model{model}.RData`; the type I error rates are
`colMeans(rej_ind)`.

### 5.3 Tables 2, S2, and S3 — parameter estimation and the ratio γ₀

```
Rscript simulation/run_table2.R [model] [n] [replicates]
```

`model` is 1–6 and `n` is the sample size; without arguments, all six models
are run at n = 200, 400, 2000, and 4000. Each run saves
`simulation/output/table2_model{model}_n{n}.RData` with an object
`summary_stats` (`Est`, `SE`, `SEE`, `Perc`).

| Model | True model |
|-------|------------|
| 1 | λ(t)=(t/2)·expit(0.5X₁+X₂) |
| 2 | λ(t)=(t/2)·Φ(0.5X₁+X₂) |
| 3 | λ(t)=(t/2)·exp(0.5X₁+X₂+W), W~N(0,1) |
| 4 | T = exp{−(0.5X₁+X₂)} + φ |
| 5 | log T = −(0.5X₁+X₂) + ε₁, ε₁~N(0,2) |
| 6 | log T = −(0.5X₁+X₂) + ε₂, ε₂~N(0,1/2) |

* **Table 2** uses n = 200 and 400; **Table S2** uses n = 2000 and 4000.
* **Table S3** reports the proportionality constant α\* as the empirical mean of
  `beta1_hat/0.5` and `beta2_hat/1`; see `simulation/summarize_simulation.R`.

### 5.4 Figure S1: m-curve (supplement, Section S.11)

```
Rscript simulation/plot_m.R
```

generates the diagnostic plot of `m(y | U, K)` for the two
accelerated-failure-time models and writes `simulation/output/m_curve.pdf` and
`simulation/output/m_curve.eps`.

### 5.5 Reading the outputs

`simulation/summarize_simulation.R` shows, with worked examples, how to load
the saved `.RData` files and recover each reported quantity.

---

## 6. Reproducing the data application (requires UK Biobank data)

Place `diabetes.tab` and `hypertension.tab` in the `data/` subdirectory.
Intermediate `.RData` objects and the figures are written to, and read from,
`application/output/`. Run the scripts in the following order, where
`{endpoint}` is `diabetes` or `hypertension`.

**Figure 1 (score statistics for 220 genetic variants):**
1. `Rscript application/fit_null.R {endpoint}`
   → `nullmodel_{endpoint}.RData`.
2. `Rscript application/score_test.R {endpoint}`
   → `score_test_{endpoint}.RData`.
3. `Rscript application/Figure1.R` → `Figure1.pdf` / `Figure1.eps`. This script
   also prints the numbers of significant genetic variants reported in
   Section 7.

**Figures S2 and S3 (goodness-of-fit residual processes):**
1. Null models: `Rscript application/fit_null.R {endpoint}` (used by the
   functional-form checks) and `Rscript application/fit_null.R {endpoint} rounded`
   (used by the proportional-hazards checks, with the interval endpoints
   rounded down to whole years)
   → `nullmodel_{endpoint}.RData`, `nullmodel_{endpoint}_rounded.RData`.
2. Functional-form checks: `Rscript application/gof_test.R {endpoint} ff [covariate]`
   → `ff_{covariate}_{endpoint}.RData`.
3. Proportional-hazards checks: `Rscript application/gof_test.R {endpoint} prop [covariate]`
   → `prop_{covariate}_{endpoint}.RData`.
4. Assemble: `Rscript application/FigureS2.R` (diabetes) and
   `Rscript application/FigureS3.R` (hypertension)
   → `FigureS2.*`, `FigureS3.*`.

Without the `covariate` argument, `gof_test.R` runs all covariates examined for
that endpoint in turn; the covariates are listed in `gof_config` in
`application/functions.R`. The p-values reported in Section 7 are printed by
`gof_test.R`; Figures S2 and S3 display the significant residual processes.
The `.eps`/`.pdf` files in `figures/` are reference copies of the published
figures.

---

## 7. C++ routines

| File | Purpose |
|------|---------|
| `src/unireg_indep_improved.cpp` | NPMLE for the (misspecified) Cox model via the EM algorithm; profile likelihood; robust (sandwich) and model-based variance estimators. |
| `src/unireg_indep_UKB.cpp` | Memory-efficient version of the EM step used for the n = 400,000 simulations. |
| `src/score_test_improved.cpp` | Robust and model-based score statistics when the covariate of interest is time-independent. |
| `src/score_test_dep_improved.cpp` | Score statistics when the covariate of interest is time-dependent (used for the application's goodness-of-fit tests). |

---

## 8. Notes

* Random seeds are drawn at run time, so results are Monte Carlo estimates over
  10,000 replicates; expect small variation in the last reported digit.
* The n = 400,000 simulations and the full 220-variant application are
  computationally intensive and were run on a shared 16-core research Linux
  server without a job scheduler. The code contains no explicit
  parallelization: each script runs as a single R process, and independent
  runs were launched concurrently with `Rscript` in separate `screen`
  sessions on the same server, e.g.,
  `screen -dmS tableS1_model1 Rscript simulation/run_tableS1.R 1`. The eight
  models of Table S1 were run concurrently in eight sessions and
  took approximately two days in total; the UK Biobank analysis was run in the
  same way and also took approximately two days.
