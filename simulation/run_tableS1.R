## ---------------------------------------------------------------------------
## Table S1: empirical Type I error rates of the model-based and robust score
## tests for the null effect of X1 with n = 400,000 (first eight models of
## Table 1), using the memory-efficient EM routine in src/unireg_indep_UKB.cpp.
##
## Run from the root directory of the reproducibility materials:
##   Rscript simulation/run_tableS1.R                # all 8 models, one after another
##   Rscript simulation/run_tableS1.R 3              # model 3
##   Rscript simulation/run_tableS1.R 3 100          # model 3, 100 replicates
##
## Arguments: [model (1-8)] [number of replicates (default 10000)].
## Each model is computationally intensive; the models can be run concurrently
## as separate processes (see README.md).
## Output: simulation/output/tableS1_model{model}.RData
## ---------------------------------------------------------------------------

if (!file.exists("src/unireg_indep_UKB.cpp")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
seed <- as.integer(abs(rnorm(1) * 100000))
set.seed(seed)
library(tidyverse)
library(Rcpp)
sourceCpp("src/unireg_indep_UKB.cpp")
sourceCpp("src/score_test_improved.cpp")
source("src/utils.R")
source("simulation/functions.R")

args = commandArgs(trailingOnly = TRUE)
models = if (length(args) >= 1) as.integer(args[1]) else 1:8
# Number of simulations
N <- if (length(args) >= 2) as.integer(args[2]) else 10000

for (model_id in models) {
  if (!(model_id %in% 1:8)) stop("model must be an integer between 1 and 8")
  
  p_value = matrix(NA, nrow = N, ncol = 3)
  rej_ind = matrix(NA, nrow = N, ncol = 3)
  for (i in 1:N) {
    cat("model =", model_id, ", i =", i,"\n")
    cat("Model fitting!\n")
    startT = Sys.time()
    simulated_data = sim_data_type1(n_sample = 400000, model = model_id, exam_scheme = "grid")
    res = fit_type1(simulated_data, eps = 0.0001, wald = FALSE)
    cat("Time difference = ",Sys.time()-startT,"\n\n")
    p_value[i,1] = 1-pchisq(res$score_stat^2/res$var1, df = 1)
    p_value[i,2] = 1-pchisq(res$score_stat^2/res$var2, df = 1)
    p_value[i,3] = 1-pchisq(res$score_stat^2/res$var, df = 1)
    rej_ind[i,1] = (p_value[i,1] < 0.05)
    rej_ind[i,2] = (p_value[i,2] < 0.05)
    rej_ind[i,3] = (p_value[i,3] < 0.05)
  }
  
  # Columns 1 and 2: model-based score tests; column 3: robust score test
  print(mean(rej_ind[,1]))
  print(mean(rej_ind[,2]))
  print(mean(rej_ind[,3]))
  
  save(rej_ind, p_value, seed, file = sprintf("simulation/output/tableS1_model%d.RData", model_id))
}
