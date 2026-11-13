## ---------------------------------------------------------------------------
## Shared functions for the simulation studies (Tables 1, 2, S1, S2, and S3
## and Figure S1). This file only defines functions; it is sourced by
## run_table1.R, run_tableS1.R, run_table2.R, and plot_m.R.
## ---------------------------------------------------------------------------

##################### Data generation #####################################

## Event times with hazard function (2 * t / scale^2) * rate, i.e.,
## (t/8) * rate for scale = 4 and (t/2) * rate for scale = 2.
hazard_event_time = function(rate, scale){
  N = length(rate)
  EventTime = rep(NA, N)
  for (i in 1:N) {
    Time = rexp(n = 1, rate = rate[i])
    EventTime[i] = scale*sqrt(Time)
  }
  return(EventTime)
}

## Six examination times per subject and the resulting interval (Left_Time,
## Right_Time] that brackets the event time.
##   exam_scheme = "continuous": six time points from Un(0,5), with the kth
##     ordered time point shifted by Un(0.1(k-1), 0.1k) (Tables 1, 2, S2, S3).
##   exam_scheme = "grid": six time points drawn with replacement from the grid
##     0.02, 0.06, ..., 0.98 (Table S1, n = 400,000).
interval_censor = function(EventTime, exam_scheme = "continuous"){
  N = length(EventTime)
  ExTime = matrix(NA, nrow = N, ncol = 6)
  Left_Time = rep(NA, N)
  Right_Time = rep(NA, N)
  for (i in 1:N) {
    if (exam_scheme == "continuous") {
      tmp_ExTime = sort(runif(n = 6, min = 0, max = 5), decreasing = F)
      ExTime[i,1] = tmp_ExTime[1] + runif(1, min = 0, max = 0.1)
      ExTime[i,2] = tmp_ExTime[2] + runif(1, min = 0.1, max = 0.2)
      ExTime[i,3] = tmp_ExTime[3] + runif(1, min = 0.2, max = 0.3)
      ExTime[i,4] = tmp_ExTime[4] + runif(1, min = 0.3, max = 0.4)
      ExTime[i,5] = tmp_ExTime[5] + runif(1, min = 0.4, max = 0.5)
      ExTime[i,6] = tmp_ExTime[6] + runif(1, min = 0.5, max = 0.6)
    } else {
      ExTime[i,] = sort(sample(x = seq(0.02,1,0.04), size = 6, replace = TRUE))
    }
    
    if(EventTime[i] > ExTime[i,6]){
      Left_Time[i] = ExTime[i,6]
      Right_Time[i] = Inf
    }else if(EventTime[i] < ExTime[i,1]){
      Left_Time[i] = 0
      Right_Time[i] = ExTime[i,1]
    }else{
      index = which(ExTime[i,] >= EventTime[i])[1]
      Left_Time[i] = ExTime[i,index-1]
      Right_Time[i] = ExTime[i,index]
    }
  }
  return(list(Left_Time = Left_Time, Right_Time = Right_Time))
}

## True models of Table 1 (model = 1, ..., 12, in the order of the rows of
## Table 1). The first eight models are also used in Table S1.
type1_event_time = function(model, X1, X2, X3){
  N = length(X1)
  switch(as.character(model),
         "1"  = hazard_event_time(exp(0.5 * X2 + 0.5 * X3), scale = 4),
         "2"  = hazard_event_time(exp(0.8 * X1^2), scale = 4),
         "3"  = hazard_event_time(exp(0.2 * X1^4 + 0.5 * X2), scale = 4),
         "4"  = hazard_event_time(exp(0.8 * X1^2 + 0.5 * X2 + 0.5 * X3), scale = 4),
         "5"  = hazard_event_time(1 + 0.5 * X2^2, scale = 4),
         "6"  = hazard_event_time(exp(X1^2) * (1 + 0.5 * X2^2), scale = 4),
         "7"  = hazard_event_time((2 + X2)^2, scale = 4),
         "8"  = hazard_event_time((1 + 0.5*exp(X1^2) + X2)^2, scale = 4),
         "9"  = exp(0.5*X2) + rexp(n = N, rate = 1),
         "10" = exp(0.8*X1^2 + 0.5*X2) + rexp(n = N, rate = 1),
         "11" = exp(0.5*X2 + rnorm(n = N, mean = 0, sd = 1)),
         "12" = exp(0.8*X1^2 + 0.5*X2 + rnorm(n = N, mean = 0, sd = 1)),
         stop("model must be an integer between 1 and 12"))
}

## Simulate one dataset for the first series of studies (Tables 1 and S1).
sim_data_type1 = function(n_sample, model, exam_scheme = "continuous"){
  N = n_sample
  Subject_ID = 1:N
  X1 = rnorm(n = N, mean = 0, sd = 1)
  X2 = rnorm(n = N, mean = 0, sd = 1)
  X3 = runif(n = N, min = -3, max = 3)
  EventTime = type1_event_time(model, X1, X2, X3)
  
  censored = interval_censor(EventTime, exam_scheme)
  SimData = data.frame(Subject_ID = Subject_ID,
                       Left_Time = censored$Left_Time,
                       Right_Time = censored$Right_Time,
                       X1 = X1, X2 = X2, Event_Time = EventTime)
  return(SimData)
}

## True models of Table 2 (model = 1, ..., 6, in the order of the rows of
## Table 2), where X is the n x 2 matrix of covariates.
ratio_event_time = function(model, X){
  N = nrow(X)
  if (model == 1) {
    linear_term = 0.5 * X[,1] + X[,2]
    EventTime = hazard_event_time(exp(linear_term)/(1+exp(linear_term)), scale = 2)
  } else if (model == 2) {
    linear_term = 0.5 * X[,1] + X[,2]
    EventTime = hazard_event_time(pnorm(linear_term), scale = 2)
  } else if (model == 3) {
    b = rnorm(n = N, mean = 0, sd = 1)
    linear_term = 0.5 * X[,1] + X[,2] + b
    EventTime = hazard_event_time(exp(linear_term), scale = 2)
  } else if (model == 4) {
    b = rexp(n = N, rate = 1)
    linear_term = 0.5 * X[,1] + X[,2]
    EventTime = exp(-linear_term) + b
  } else if (model == 5) {
    b = rnorm(n = N, mean = 0, sd = sqrt(2))
    linear_term = 0.5 * X[,1] + X[,2]
    EventTime = exp(-linear_term + b)
  } else if (model == 6) {
    b = rnorm(n = N, mean = 0, sd = sqrt(0.5))
    linear_term = 0.5 * X[,1] + X[,2]
    EventTime = exp(-linear_term + b)
  } else {
    stop("model must be an integer between 1 and 6")
  }
  return(EventTime)
}

## Simulate one dataset for the second series of studies (Tables 2, S2, S3).
sim_data_ratio = function(n_sample, model){
  N = n_sample
  Subject_ID = 1:N
  
  mu = c(0,0)
  Sigma <- matrix(c(1, 0.3,
                    0.3, 1), nrow = 2)
  
  X = mvrnorm(n = N, mu = mu, Sigma = Sigma)
  EventTime = ratio_event_time(model, X)
  
  censored = interval_censor(EventTime, "continuous")
  SimData = data.frame(Subject_ID = Subject_ID,
                       Left_Time = censored$Left_Time,
                       Right_Time = censored$Right_Time,
                       X1 = X[,1], X2 = X[,2], Event_Time = EventTime)
  return(SimData)
}

##################### Model fitting #######################################

## Convert a simulated dataset into the inputs of the C++ routines.
prepare_data = function(SimData){
  simpdat_simulated = data.frame(Subject_ID = 1:nrow(SimData), 
                                 SimData[,c(2,3)],
                                 R_star = ifelse(SimData$Right_Time==Inf, SimData$Left_Time, SimData$Right_Time))
  simpdat_simulated$new_ID = 1:nrow(simpdat_simulated)
  simpdat_simulated = simpdat_simulated[order(simpdat_simulated$R_star, decreasing = TRUE),]
  Time_simulated = find_time(simpdat_simulated)
  simpdat_simulated = as.matrix(simpdat_simulated)
  return(list(simpdat = simpdat_simulated, Time = Time_simulated))
}

## Fit the working Cox model with X1 and X2 and test the null effect of X1.
##   eps:  convergence threshold of the EM algorithm.
##   wald: if TRUE, the Wald statistics are also computed (Table 1); if FALSE,
##         only the score statistics are computed (Table S1).
fit_type1 = function(SimData, eps = 0.001, wald = TRUE) {
  prepared = prepare_data(SimData)
  simpdat_simulated = prepared$simpdat
  Time_simulated = prepared$Time
  
  # Wald
  if (wald) {
    X_simulated = as.matrix(model.matrix(~ X2 + X1, SimData)[,-1])
    model = unireg_indep_EM_improved(X_simulated, simpdat_simulated, eps, Time_simulated)
  }
  
  # Score
  X_rest = as.matrix(model.matrix(~ X2, SimData)[,-1])
  model_rest = unireg_indep_EM_improved(X_rest, simpdat_simulated, eps, Time_simulated)
  beta_initial = c(model_rest$beta, 0)
  lambda_initial = model_rest$lambda
  X_new_score = cbind(X_rest, SimData$X1)
  score = score_test_improved(beta_initial, X_new_score, lambda_initial, simpdat_simulated, eps, Time_simulated, 2.5)
  var1 = score$Info1[2,2]-score$Info1[1,2]^2/score$Info1[1,1]
  var2 = score$Info2[2,2]-score$Info2[1,2]^2/score$Info2[1,1]
  var = mean((score$Dpl[,2]-score$Info2[1,2]/score$Info2[1,1]*score$Dpl[,1])^2)
  
  results = tibble(
    score_stat = score$Stat,
    var1 = var1,
    var2 = var2,
    var = var
  )
  if (wald) {
    results = bind_cols(
      tibble(
        beta_hat = model$beta[2,1],
        beta_se = sqrt(model$Cov_robust[2,2]),
        beta_se1 = sqrt(model$Cov1[2,2]),
        beta_se2 = sqrt(model$Cov2[2,2])
      ),
      results
    )
  }
  return(results)
}

## Fit the working Cox model with X1 and X2 (Tables 2, S2, S3).
fit_ratio <- function(SimData) {
  prepared = prepare_data(SimData)
  simpdat_simulated = prepared$simpdat
  Time_simulated = prepared$Time
  X_simulated = as.matrix(model.matrix(~ X1 + X2, SimData)[,-1])
  model = unireg_indep_EM_improved(X_simulated, simpdat_simulated, 0.001, Time_simulated)
  tibble(
    beta1_hat = model$beta[1,1],
    beta2_hat = model$beta[2,1],
    beta1_se = sqrt(model$Cov_robust[1,1]),
    beta2_se = sqrt(model$Cov_robust[2,2]),
    beta_cov = model$Cov_robust[1,2],
    beta1_se1 = sqrt(model$Cov1[1,1]),
    beta1_se2 = sqrt(model$Cov2[1,1]),
    beta2_se1 = sqrt(model$Cov1[2,2]),
    beta2_se2 = sqrt(model$Cov2[2,2])
  )
}

##################### Summaries ###########################################

## Wald tests in Table 1: alpha is the robust test (sandwich variance), and
## alpha1 and alpha2 are the model-based tests.
summarize_wald = function(simulation_results){
  simulation_results %>%
    summarize(
      beta_mean = mean(beta_hat),
      beta_bias = mean(beta_hat),
      beta_sd = sd(beta_hat),
      beta_see = mean(beta_se),
      alpha = mean(abs(beta_hat)/beta_se > qnorm(0.975)),
      beta_see1 = mean(beta_se1),
      alpha1 = mean(abs(beta_hat)/beta_se1 > qnorm(0.975)),
      beta_see2 = mean(beta_se2),
      alpha2 = mean(abs(beta_hat)/beta_se2 > qnorm(0.975))
    )
}

## Score tests in Table 1: alpha is the robust test, and alpha1 and alpha2 are
## the model-based tests.
summarize_score = function(simulation_results){
  simulation_results %>%
    summarize(
      alpha = mean(score_stat^2/var > qchisq(0.95,1)),
      alpha1 = mean(score_stat^2/var1 > qchisq(0.95,1)),
      alpha2 = mean(score_stat^2/var2 > qchisq(0.95,1))
    )
}

## Tables 2 and S2: the true ratio is gamma_0 = beta_01 / beta_02 = 0.5.
summarize_ratio = function(simulation_results){
  simulation_results %>%
    summarize(
      beta1_mean = mean(beta1_hat),
      beta2_mean = mean(beta2_hat),
      ratio_mean = mean(beta1_hat/beta2_hat),
      beta1_sd = sd(beta1_hat),
      beta2_sd = sd(beta2_hat),
      ratio_sd = sd(beta1_hat/beta2_hat),
      beta1_see = mean(beta1_se),
      beta2_see = mean(beta2_se),
      ratio_see = mean(sqrt((beta1_hat/beta2_hat)^2*beta2_se^2 + beta1_se^2-2*(beta1_hat/beta2_hat)*beta_cov) / beta2_hat),
      P1 = mean(beta1_hat > 0),
      P2 = mean(beta2_hat > 0),
      CP = mean(beta1_hat/beta2_hat - qnorm(0.975)*sqrt((beta1_hat/beta2_hat)^2*beta2_se^2 + beta1_se^2-2*(beta1_hat/beta2_hat)*beta_cov) / beta2_hat < 0.5 &
                  beta1_hat/beta2_hat + qnorm(0.975)*sqrt((beta1_hat/beta2_hat)^2*beta2_se^2 + beta1_se^2-2*(beta1_hat/beta2_hat)*beta_cov) / beta2_hat > 0.5)
    )
}
