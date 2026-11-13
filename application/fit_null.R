## ---------------------------------------------------------------------------
## Fit the null Cox model (baseline covariates only) for one endpoint.
##
## Run from the root directory of the reproducibility materials:
##   Rscript application/fit_null.R diabetes
##   Rscript application/fit_null.R diabetes rounded
##   Rscript application/fit_null.R hypertension
##   Rscript application/fit_null.R hypertension rounded
##
## Arguments: <endpoint (diabetes or hypertension)> [rounded].
## With "rounded", the interval endpoints are rounded down to whole years; this
## version of the null model is used for the proportional hazards checks.
## Input:  data/{endpoint}.tab
## Output: application/output/nullmodel_{endpoint}.RData or
##         application/output/nullmodel_{endpoint}_rounded.RData
## ---------------------------------------------------------------------------

if (!file.exists("src/unireg_indep_improved.cpp")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
library(data.table)
library(Rcpp)
library(dplyr)
sourceCpp("src/unireg_indep_improved.cpp")
source("src/utils.R")
source("application/functions.R")

args = commandArgs(trailingOnly = TRUE)
endpoint = check_endpoint(args[1])
rounded = length(args) >= 2 && args[2] == "rounded"
config = endpoint_config[[endpoint]]

df = fread(config$data_file)
print(sum(df$left>0 & df$right<Inf))
print(sum(df$left==0))
print(sum(df$right==Inf))
df$f.eid = 1:nrow(df)
if (rounded) {
  round_divisor = config$round_divisor
  df = df %>%
    mutate(left = floor(left/round_divisor),
           right = floor(right/round_divisor))
}

## Data processing
n = dim(df)[1]
simpdat = data.frame(Subject_ID = df$f.eid, 
                     df[,c(2,3)],
                     R_star = ifelse(df$right==Inf, df$left, df$right))
simpdat$new_ID = df$f.eid
simpdat = simpdat[order(simpdat$R_star, decreasing = TRUE),]
Time = find_time(simpdat)
simpdat = as.matrix(simpdat)

X = model.matrix(~ sex+white+England+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, df)[,-1]
## Model fitting
cat("Model fitting!\n")
startT = Sys.time()
result = unireg_indep_EM_improved(X, simpdat, 0.001, Time)
cat("Time difference = ",Sys.time()-startT,"\n\n")

print(result$beta)
print(sqrt(diag(result$Cov1)))
print(sqrt(diag(result$Cov2)))
print(sqrt(diag(result$Cov_robust)))

save(result, X, simpdat, Time,
     file = sprintf("application/output/nullmodel_%s%s.RData", endpoint, if (rounded) "_rounded" else ""))
