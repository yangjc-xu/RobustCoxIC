## ---------------------------------------------------------------------------
## Tables 2, S2, and S3: parameter estimation under the misspecified Cox model
## (Table 2: n = 200 and 400; Table S2: n = 2000 and 4000; Table S3: estimates
## of alpha* at all four sample sizes).
##
## Run from the root directory of the reproducibility materials:
##   Rscript simulation/run_table2.R                 # all 6 models, n = 200, 400, 2000, 4000
##   Rscript simulation/run_table2.R 3               # model 3, all four sample sizes
##   Rscript simulation/run_table2.R 3 200           # model 3, n = 200
##   Rscript simulation/run_table2.R 3 200 100       # model 3, n = 200, 100 replicates
##
## Arguments: [model (1-6)] [sample size] [number of replicates (default 10000)].
## The models are numbered in the order of the rows of Table 2.
## Output: simulation/output/table2_model{model}_n{n}.RData
## ---------------------------------------------------------------------------

if (!file.exists("src/unireg_indep_improved.cpp")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
seed <- as.integer(abs(rnorm(1) * 100000))
set.seed(seed)
library(tidyverse)
library(Rcpp)
library(MASS)
sourceCpp("src/unireg_indep_improved.cpp")
source("src/utils.R")
source("simulation/functions.R")

args = commandArgs(trailingOnly = TRUE)
models = if (length(args) >= 1) as.integer(args[1]) else 1:6
sizes = if (length(args) >= 2) as.integer(args[2]) else c(200, 400, 2000, 4000)
# Number of simulations
N <- if (length(args) >= 3) as.integer(args[3]) else 10000

for (model_id in models) {
  for (n_sample in sizes) {
    cat("Table 2: model", model_id, ", n =", n_sample, "\n")
    
    # Generate datasets
    simulated_data <- tibble(sim_id = 1:N)
    simulated_data$data <- map(simulated_data$sim_id,
                               function(id) sim_data_ratio(n_sample = n_sample, model = model_id))
    
    simulation_results <- simulated_data %>%
      mutate(results = map(data, fit_ratio)) %>%
      unnest(results)
    simulation_results$data <- NULL
    
    summary_stats <- summarize_ratio(simulation_results)
    print(summary_stats)
    
    save(simulation_results, summary_stats, seed,
         file = sprintf("simulation/output/table2_model%d_n%d.RData", model_id, n_sample))
  }
}
