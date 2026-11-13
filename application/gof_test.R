## ---------------------------------------------------------------------------
## Goodness-of-fit tests of the null Cox models (supremum tests based on
## cumulative sums of residuals): functional forms of covariates ("ff") and the
## proportional hazards assumption ("prop").
##
## Run from the root directory of the reproducibility materials, after
## application/fit_null.R has been run for the same endpoint (without "rounded"
## for the functional-form checks and with "rounded" for the proportional
## hazards checks):
##   Rscript application/gof_test.R hypertension ff pc4      # one covariate
##   Rscript application/gof_test.R hypertension prop sex
##   Rscript application/gof_test.R hypertension prop        # all covariates in turn
##
## Arguments: <endpoint (diabetes or hypertension)> <ff or prop> [covariate].
## The covariates examined for each endpoint are listed in gof_config in
## application/functions.R.
## Input:  application/output/nullmodel_{endpoint}.RData (ff) or
##         application/output/nullmodel_{endpoint}_rounded.RData (prop)
## Output: application/output/{ff or prop}_{covariate}_{endpoint}.RData
## ---------------------------------------------------------------------------

if (!file.exists("src/score_test_improved.cpp")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/score_test_improved.cpp")
sourceCpp("src/score_test_dep_improved.cpp")
source("application/functions.R")

gof_test = function(endpoint, type, covariate){
  cat("Endpoint:", endpoint, "; test:", type, "; covariate:", covariate, "\n")
  
  # loads result, X, simpdat, and Time
  if (type == "ff") {
    load(sprintf("application/output/nullmodel_%s.RData", endpoint))
  } else {
    load(sprintf("application/output/nullmodel_%s_rounded.RData", endpoint))
  }
  
  n = nrow(X)
  constant = 5
  n_sim = 1000
  p = length(result$beta)
  beta_new = c(result$beta, 0)
  col = match(covariate, covariate_names)
  
  if (type == "ff") {
    grid_range = gof_config[[endpoint]]$ff[[covariate]]
    if (is.null(grid_range)) stop("no functional-form check is defined for this covariate")
    x = seq(grid_range[1], grid_range[2], length.out = 101)
  } else {
    if (!(covariate %in% gof_config[[endpoint]]$prop)) stop("no proportional hazards check is defined for this covariate")
    x = sort(Time, decreasing = F)
  }
  W = numeric(length(x))
  W_simu1 = matrix(data = 0, nrow = n_sim, ncol = length(x))
  G = matrix(data = rnorm(n = n_sim*n), nrow = n_sim, ncol = n)
  
  for (i in 1:length(x)) {
    cat("i =", i, "\n")
    cat("x =", x[i], "\n")
    
    if (type == "ff") {
      X_extra = ifelse(X[,col] <= x[i], 1, 0)
      X_new = cbind(X, X_extra)
      res = score_test_improved(beta_new, X_new, result$lambda, simpdat, 0.001, Time, constant)
    } else {
      # the genetic principal components are standardized; the binary covariates are not
      if (covariate %in% c("sex", "white", "england")) {
        X_dep = (X[,col]) %x% t(Time <= x[i])
      } else {
        X_dep = (scale(X[,col])) %x% t(Time <= x[i])
      }
      res = score_test_dep_improved(beta_new, X, X_dep, result$lambda, simpdat, 0.001, Time, constant)
    }
    W[i] = res$Stat
    
    sum_unit1 = as.vector(res$Dpl[,p+1] - 
                            res$Info1[p+1,1:p]%*%
                            solve(res$Info1[1:p,1:p])%*%
                            t(res$Dpl[,1:p]))
    for (j in 1:n_sim) { W_simu1[j,i] = sum(G[j,] * sum_unit1)/sqrt(n) }
  }
  
  pval = mean(max(abs(W)) <= apply(abs(W_simu1), MARGIN = 1, FUN = max))
  cat("P-value =", pval, "\n\n")
  
  save(W, W_simu1, x, pval,
       file = sprintf("application/output/%s_%s_%s.RData", type, covariate, endpoint))
}

args = commandArgs(trailingOnly = TRUE)
endpoint = check_endpoint(args[1])
type = args[2]
if (is.na(type) || !(type %in% c("ff", "prop"))) stop("the second argument must be 'ff' or 'prop'")
if (length(args) >= 3) {
  covariates = args[3]
} else if (type == "ff") {
  covariates = names(gof_config[[endpoint]]$ff)
} else {
  covariates = gof_config[[endpoint]]$prop
}

for (covariate in covariates) {
  gof_test(endpoint, type, covariate)
}
