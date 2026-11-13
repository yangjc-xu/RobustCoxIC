## ---------------------------------------------------------------------------
## Example helper for inspecting the saved simulation outputs.
##
## run_table1.R, run_tableS1.R, and run_table2.R save their results to
## simulation/output/. This script shows how to read those files and recover
## the quantities reported in the manuscript and the supplementary materials.
## Run it from the root directory of the reproducibility materials.
## ---------------------------------------------------------------------------

library(tidyverse)

## --- Table 1 (n = 200, 400): type I error rates ---
## run_table1.R saves summary_stats_wald and summary_stats_score, where
##   $alpha           = robust test (sandwich variance)
##   $alpha1, $alpha2 = model-based tests (two perturbation-based variance
##                      estimators of the profile likelihood)
## Example: true model 5 (Table 1) at n = 400.
load("simulation/output/table1_model5_n400.RData")
summary_stats_wald
summary_stats_score

## --- Table S1 (n = 400,000): type I error rates ---
## run_tableS1.R saves rej_ind (rejection indicators) and p_value; the type I
## error rates are the column means of rej_ind.
## Example:
## load("simulation/output/tableS1_model5.RData"); colMeans(rej_ind)

## --- Table 2 (n = 200, 400) and Table S2 (n = 2000, 4000): estimation ---
## run_table2.R saves summary_stats with the estimate (Est), empirical standard
## error (SE), mean standard error estimate (SEE), and sign/coverage percentage
## (Perc). Example: true model 6 (Table 2) at n = 4000.
load("simulation/output/table2_model6_n4000.RData")
summary_stats

## --- Table S3: proportionality constant alpha* ---
## alpha* is the empirical mean of the estimate divided by the true coefficient
## (beta_01 = 0.5 and beta_02 = 1).
summary_alpha <- simulation_results %>%
  summarize(
    alpha_beta1 = mean(beta1_hat / 0.5),
    alpha_beta2 = mean(beta2_hat / 1)
  )
summary_alpha
