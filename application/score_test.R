## ---------------------------------------------------------------------------
## Model-based and robust score tests for the 220 genetic variants.
##
## Run from the root directory of the reproducibility materials, after
## application/fit_null.R has been run for the same endpoint:
##   Rscript application/score_test.R diabetes
##   Rscript application/score_test.R hypertension
##
## Argument: <endpoint (diabetes or hypertension)>.
## Input:  data/{endpoint}.tab, application/output/nullmodel_{endpoint}.RData
## Output: application/output/score_test_{endpoint}.RData
## ---------------------------------------------------------------------------

if (!file.exists("src/score_test_improved.cpp")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
library(data.table)
library(Rcpp)
library(tidyverse)
sourceCpp("src/score_test_improved.cpp")
source("application/functions.R")

args = commandArgs(trailingOnly = TRUE)
endpoint = check_endpoint(args[1])
config = endpoint_config[[endpoint]]

load(sprintf("application/output/nullmodel_%s.RData", endpoint))
df = fread(config$data_file)
df$f.eid = 1:nrow(df)
variant = as.matrix(df[, config$variant_cols, with = FALSE])
variant_name = colnames(variant)
statistic = rep(NA,length(variant_name))
variance = matrix(NA, nrow = length(variant_name), ncol = 3)
p_value = matrix(NA, nrow = length(variant_name), ncol = 3)

# Score
# Columns 1 and 2 of variance and p_value: model-based; column 3: robust
p = length(result$beta)
beta_initial = c(result$beta, 0)
lambda_initial = result$lambda
for (i in 1:length(variant_name)) {
  cat("i =", i, "\n")
  startT = Sys.time()
  X_new = as.matrix(cbind(X, variant[,i]))
  score = score_test_improved(beta_initial, X_new, lambda_initial, simpdat, 0.001, Time, 2.5)
  statistic[i] = score$Stat
  variance[i,1] = as.numeric(score$Info1[p+1,p+1]-score$Info1[p+1,1:p]%*%solve(score$Info1[1:p,1:p])%*%score$Info1[1:p,p+1])
  variance[i,2] = as.numeric(score$Info2[p+1,p+1]-score$Info2[p+1,1:p]%*%solve(score$Info2[1:p,1:p])%*%score$Info2[1:p,p+1])
  variance[i,3] = mean((score$Dpl[,p+1]-score$Info2[p+1,1:p]%*%solve(score$Info2[1:p,1:p])%*%t(score$Dpl[,1:p]))^2)
  p_value[i,1] = 1-pchisq(statistic[i]^2/variance[i,1], df = 1)
  p_value[i,2] = 1-pchisq(statistic[i]^2/variance[i,2], df = 1)
  p_value[i,3] = 1-pchisq(statistic[i]^2/variance[i,3], df = 1)
  cat("Time difference = ",Sys.time()-startT,"\n\n")
}

save(statistic, variance, p_value, variant_name,
     file = sprintf("application/output/score_test_%s.RData", endpoint))
