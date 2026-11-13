## ---------------------------------------------------------------------------
## Table 1: empirical Type I error rates of the model-based and robust Wald and
## score tests for the null effect of X1 (n = 200 and 400).
##
## Run from the root directory of the reproducibility materials:
##   Rscript simulation/run_table1.R                 # all 12 models, n = 200 and 400
##   Rscript simulation/run_table1.R 3               # model 3, n = 200 and 400
##   Rscript simulation/run_table1.R 3 200           # model 3, n = 200
##   Rscript simulation/run_table1.R 3 200 100       # model 3, n = 200, 100 replicates
##
## Arguments: [model (1-12)] [sample size] [number of replicates (default 10000)].
## The models are numbered in the order of the rows of Table 1.
## Output: simulation/output/table1_model{model}_n{n}.RData
## ---------------------------------------------------------------------------

if (!file.exists("src/unireg_indep_improved.cpp")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
seed <- as.integer(abs(rnorm(1) * 100000))
set.seed(seed)
library(tidyverse)
library(Rcpp)
sourceCpp("src/unireg_indep_improved.cpp")
sourceCpp("src/score_test_improved.cpp")
source("src/utils.R")
source("simulation/functions.R")

args = commandArgs(trailingOnly = TRUE)
models = if (length(args) >= 1) as.integer(args[1]) else 1:12
sizes = if (length(args) >= 2) as.integer(args[2]) else c(200, 400)
# Number of simulations
N <- if (length(args) >= 3) as.integer(args[3]) else 10000

for (model_id in models) {
  for (n_sample in sizes) {
    cat("Table 1: model", model_id, ", n =", n_sample, "\n")
    
    # Generate datasets
    simulated_data <- tibble(sim_id = 1:N)
    simulated_data$data <- map(simulated_data$sim_id,
                               function(id) sim_data_type1(n_sample = n_sample, model = model_id))
    
    simulation_results <- simulated_data %>%
      mutate(results = map(data, fit_type1)) %>%
      unnest(results)
    simulation_results$data <- NULL
    
    summary_stats_wald <- summarize_wald(simulation_results)
    summary_stats_score <- summarize_score(simulation_results)
    print(summary_stats_wald)
    print(summary_stats_score)
    
    save(simulation_results, summary_stats_wald, summary_stats_score, seed,
         file = sprintf("simulation/output/table1_model%d_n%d.RData", model_id, n_sample))
  }
}
