# Functions to fit and summarize LowFR models as described in "Low-Rank
# Longitudinal Factor Regression" paper

library(tidyverse)
library(rstan)
source("R/helper.R")
source("R/sim_data.R")

###############################################################################
############# Fit a Stan model and return specified output ####################
###############################################################################
fit_LowFR <- function(y_obs, X_obs, p=10, k=NULL, TT=3,
                      output="all",
                      burnin=1000, samples=1000, chains=4,
                      random_seed=1234) {
  
  # Note: output can be either "stan_fit", "all", or a list of parameters 
  # for which to return posterior samples. This list can include items from
  #       "alpha_0"
  #       "alpha"
  #       "Gamma"
  #       "sigma2"
  #       "Sigma"
  #       "phi"
  #       "tau"
  #       "theta"
  #       "Omega"
  # If "all" is specified, samples for the full list above will be returned as
  # a list of arrays.
  
  # detect n_obs
  n_obs <- length(y_obs)
  
  # choose k using SVD if needed
  if (is.null(k)) {
    k <- k_svd_LowFR(X_obs=X_obs, p=p, TT=TT)
  }
  
  # compile stan model
  m <- stan_model("Stan/LowFR.stan")
  
  # fit model
  options(mc.cores = parallel::detectCores())
  fit <- sampling(m,
                  list(N=n_obs,
                       p=p,
                       k=k,
                       TT=TT,
                       H=min(p,TT),
                       X=X_obs,
                       y=y_obs),
                  chains=chains,
                  iter=burnin+samples,
                  warmup=burnin,
                  seed=random_seed,
                  init=0)
  
  # return the stan fit if specified as output
  if (output == "stan_fit") {
    return(fit)
  }
  
  # otherwise, extract posterior samples and return list of specified parameters
  post_samples <- rstan::extract(fit)
  results <- list()
  result_names <- c()
  if ("alpha_0" %in% output | output == "all") {
    #results <- append(results, post_samples$alpha_0)
    results <- c(results, list(post_samples$alpha_0))
    result_names <- c(result_names, "alpha_0")
  }
  if ("alpha" %in% output | output == "all") {
    #results <- append(results, post_samples$alpha)
    results <- c(results, list(post_samples$alpha))
    result_names <- c(result_names, "alpha")
  }
  if ("Gamma" %in% output | output == "all") {
    #results <- append(results, post_samples$Gamma)
    results <- c(results, list(post_samples$Gamma))
    result_names <- c(result_names, "Gamma")
  }
  if ("sigma2" %in% output | output == "all") {
    #results <- append(results, post_samples$sigma2)
    results <- c(results, list(post_samples$sigma2))
    result_names <- c(result_names, "sigma2")
  }
  if ("Sigma" %in% output | output == "all") {
    #results <- append(results, post_samples$Sigma)
    results <- c(results, list(post_samples$Sigma))
    result_names <- c(result_names, "Sigma")
  }
  if ("phi" %in% output | output == "all") {
    #results <- append(results, post_samples$phi)
    results <- c(results, list(post_samples$phi))
    result_names <- c(result_names, "phi")
  }
  if ("tau" %in% output | output == "all") {
    #results <- append(results, post_samples$tau)
    results <- c(results, list(post_samples$tau))
    result_names <- c(result_names, "tau")
  }
  if ("theta" %in% output | output == "all") {
    #results <- append(results, post_samples$theta)
    results <- c(results, list(post_samples$theta))
    result_names <- c(result_names, "theta")
  }
  if ("Omega" %in% output | output == "all") {
    #results <- append(results, post_samples$Omega)
    results <- c(results, list(post_samples$Omega))
    result_names <- c(result_names, "Omega")
  }
  names(results) <- result_names
  return(results)
}

###############################################################################
######################### Choose k based on SVD ###############################
###############################################################################
k_svd_LowFR <- function(X_obs, p, TT) {
  # create data matrix of x_it's to do SVD and determine number of factors
  Xit_mat <- matrix(nrow=(TT*nrow(X_obs)), ncol=p)
  
  # loop over data to populate matrix
  for (i in 1:nrow(X_obs)) {
    for (t in 1:TT) {
      Xit_mat[(TT*(i-1) + t),] <- X_obs[i,seq(from=t, to=p*TT, by=TT)]
    }
  }
  
  # do SVD to get singular values
  sing_vals <- svd(Xit_mat)$d
  frac_sing_vals <- rep(NA, p)
  for (i in 1:p) {
    frac_sing_vals[i] <- sum(sing_vals[1:i]) / sum(sing_vals)
  }
  
  # find index where frac first exceeds 0.9
  for (i in 1:p) {
    if (frac_sing_vals[i] > 0.9) {
      print(paste0("Setting k = ", i, " based on SVD."))
      return(i)
    }
  }
  print("k selection failed. Setting k = p.")
  return(p)
}

###############################################################################
############## Calculate main and interaction MSE of a fitted model ###########
###############################################################################
get_mse_LowFR <- function(data, fit) {
  # extract main and interaction samples and true values
  alpha_samps <- fit$alpha
  Gamma_samps <- fit$Gamma
  alpha_true <- data$alpha
  Gamma_true <- data$Gamma
  
  # summarize main effect mse
  pT <- length(alpha_samps[1,])
  main_mse_vec <- rep(NA, pT)
  for (i in 1:pT) {
    main_mse_vec[i] <- squared_error(true_val=alpha_true[i],
                                     est_vec=alpha_samps[,i])
  }
  main_mse <- mean(main_mse_vec)
  
  # summarize interaction mse
  int_mse_vec <- rep(NA, pT*(pT+1)/2)
  curr_index <- 1
  for (i in 1:pT) {
    for (j in i:pT) {
      if (i == j) {
        int_mse_vec[curr_index] <- squared_error(true_val=Gamma_true[i,j],
                                                 est_vec=Gamma_samps[,i,j])
        curr_index = curr_index + 1
      }
      else {
        int_mse_vec[curr_index] <- squared_error(true_val=Gamma_true[i,j] +
                                                   Gamma_true[j,i],
                                                 est_vec=Gamma_samps[,i,j] +
                                                   Gamma_samps[,j,i])
        curr_index <- curr_index + 1
      }
    }
  }
  int_mse <- mean(int_mse_vec)
  
  # return calculated values
  output <- list(main_mse, int_mse)
  names(output) <- c("main_mse", "int_mse")
  return(output)
}

###############################################################################
############ Calculate TP and TN rates for main and interaction effects #######
###############################################################################
get_TP_TN_rates_LowFR <- function(data, fit) {
  # extract main and interaction samples and true values
  alpha_samps <- fit$alpha
  Gamma_samps <- fit$Gamma
  alpha_true <- data$alpha
  Gamma_true <- data$Gamma
  
  # summarize main effect TP, TN, FP, and FNs
  pT <- length(alpha_samps[1,])
  main_ss_vec <- rep(NA, pT)
  for (i in 1:pT) {
    main_ss_vec[i] <- sens_spec(true_val=alpha_true[i],
                                     est_vec=alpha_samps[,i])
  }
  if (sum(main_ss_vec == "TP") > 0) {
    main_TP_rate <- sum(main_ss_vec=="TP") / (sum(main_ss_vec=="TP") +
                                                sum(main_ss_vec=="FN"))
  }
  else {
    main_TP_rate <- 0
  }
  if (sum(main_ss_vec == "TN") > 0) {
    main_TN_rate <- sum(main_ss_vec=="TN") / (sum(main_ss_vec=="TN") +
                                                sum(main_ss_vec=="FP"))
  }
  else {
    main_TN_rate <- 0
  }
  
  # summarize interaction TP, TN, FP, and FNs
  int_ss_vec <- rep(NA, pT*(pT+1)/2)
  curr_index <- 1
  for (i in 1:pT) {
    for (j in i:pT) {
      if (i == j) {
        int_ss_vec[curr_index] <- sens_spec(true_val=Gamma_true[i,j],
                                                 est_vec=Gamma_samps[,i,j])
        curr_index = curr_index + 1
      }
      else {
        int_ss_vec[curr_index] <- sens_spec(true_val=Gamma_true[i,j] +
                                                   Gamma_true[j,i],
                                                 est_vec=Gamma_samps[,i,j] +
                                                   Gamma_samps[,j,i])
        curr_index <- curr_index + 1
      }
    }
  }
  if (sum(int_ss_vec=="TP") > 0) {
    int_TP_rate <- sum(int_ss_vec=="TP") / (sum(int_ss_vec=="TP") +
                                              sum(int_ss_vec=="FN"))
  }
  else {
    int_TP_rate <- 0
  }
  if (sum(int_ss_vec=="TN") > 0) {
    int_TN_rate <- sum(int_ss_vec=="TN") / (sum(int_ss_vec=="TN") +
                                              sum(int_ss_vec=="FP"))
  }
  else {
    int_TN_rate <- 0
  }
  
  # return calculated values
  output <- list(main_TP_rate, main_TN_rate, int_TP_rate, int_TN_rate)
  names(output) <- c("main_TP_rate", "main_TN_rate", "int_TP_rate", "int_TN_rate")
  return(output)
}

###############################################################################
############ Calculate 95% CI coverage for main and interaction effects #######
###############################################################################
get_coverage_LowFR <- function(data, fit) {
  # extract main and interaction samples and true values
  alpha_samps <- fit$alpha
  Gamma_samps <- fit$Gamma
  alpha_true <- data$alpha
  Gamma_true <- data$Gamma
  
  # summarize main effect mse
  pT <- length(alpha_samps[1,])
  main_coverage_vec <- rep(NA, pT)
  for (i in 1:pT) {
    main_coverage_vec[i] <- cover_bool(true_val=alpha_true[i],
                                     est_vec=alpha_samps[,i])
  }
  main_coverage <- mean(main_coverage_vec)
  
  # summarize interaction mse
  int_coverage_vec <- rep(NA, pT*(pT+1)/2)
  curr_index <- 1
  for (i in 1:pT) {
    for (j in i:pT) {
      if (i == j) {
        int_coverage_vec[curr_index] <- cover_bool(true_val=Gamma_true[i,j],
                                                 est_vec=Gamma_samps[,i,j])
        curr_index = curr_index + 1
      }
      else {
        int_coverage_vec[curr_index] <- cover_bool(true_val=Gamma_true[i,j] +
                                                   Gamma_true[j,i],
                                                 est_vec=Gamma_samps[,i,j] +
                                                   Gamma_samps[,j,i])
        curr_index <- curr_index + 1
      }
    }
  }
  int_coverage <- mean(int_coverage_vec)
  
  # return calculated values
  output <- list(main_coverage, int_coverage)
  names(output) <- c("main_coverage", "int_coverage")
  return(output)
}

###############################################################################
######################### Summarize all accuracy results ######################
###############################################################################
summarize_accuracy_LowFR <- function(data, fit) {
  mse <- get_mse_LowFR(data, fit)
  coverage <- get_coverage_LowFR(data, fit)
  TP_TN_rate <- get_TP_TN_rates_LowFR(data, fit)
  output <- c(mse, coverage, TP_TN_rate)
  return(output)
}

###############################################################################
############ Summarize main and interaction effective sample sizes ############
###############################################################################
get_ESS_LowFR <- function(fit) {
  # extract main and interaction samples
  alpha_samps <- fit$alpha
  Gamma_samps <- fit$Gamma
  
  # create vector and matrix to store ESS values
  pT <- length(alpha_samps[1,])
  main_ESS <- rep(NA, pT)
  int_ESS <- matrix(data=rep(NA, pT*pT), nrow=pT)
  
  # populate with ESS values
  for (i in 1:pT) {
    main_ESS[i] <- get_ESS(alpha_samps[,i])
    for (j in 1:pT) {
      int_ESS[i,j] <- get_ESS(Gamma_samps[,i,j])
    }
  }
  
  # calculate min values
  main_ESS_min <- min(main_ESS)
  int_ESS_min <- min(int_ESS)
  
  # return output
  output <- list(main_ESS, int_ESS, main_ESS_min, int_ESS_min)
  names(output) <- c("main_ESS", "int_ESS", "main_ESS_min", "int_ESS_min")
  return(output)
}

###############################################################################
######## Summarize main and interaction Gelman-Rubin diagnostic (Rhat) ########
###############################################################################
get_Rhat_LowFR <- function(fit, num_samples=1000, num_chains=4) {
  # extract main and interaction samples
  alpha_samps <- fit$alpha
  Gamma_samps <- fit$Gamma
  
  # create vector and matrix to store ESS values
  pT <- length(alpha_samps[1,])
  main_Rhat <- rep(NA, pT)
  int_Rhat <- matrix(data=rep(NA, pT*pT), nrow=pT)
  
  # populate with ESS values
  for (i in 1:pT) {
    main_Rhat[i] <- get_Rhat(alpha_samps[,i], num_samples=num_samples, num_chains=num_chains)
    for (j in 1:pT) {
      int_Rhat[i,j] <- get_Rhat(Gamma_samps[,i,j], num_samples=num_samples, num_chains=num_chains)
    }
  }
  
  # calculate min values
  main_Rhat_max <- max(main_Rhat)
  int_Rhat_max <- max(int_Rhat)
  
  # return output
  output <- list(main_Rhat, int_Rhat, main_Rhat_max, int_Rhat_max)
  names(output) <- c("main_Rhat", "int_Rhat", "main_Rhat_max", "int_Rhat_max")
  return(output)
}

###############################################################################
########################### Summarize mixing diagnostics ######################
###############################################################################
summarize_mixing_LowFR <- function(fit, num_samples=1000, num_chains=4) {
  ESS <- get_ESS_LowFR(fit)
  Rhat <- get_Rhat_LowFR(fit, num_samples, num_chains)
  output <- c(ESS, Rhat)
  return(output)
}


###############################################################################
########################### Run a full LowFR simulation #######################
###############################################################################
run_sim_LowFR <- function(scenario,
                          save_output=FALSE,
                          path="",
                          datagen_seed=1234,
                          stan_seed=1234,
                          n_obs=200,
                          p=10,
                          k=5,
                          k_model=NULL,
                          TT=3,
                          parameter_output="all",
                          burnin=1000,
                          samples=1000,
                          chains=4,
                          simulation_output=c("post_samples",
                                              "accuracy_summary",
                                              "diagnostic_summary")) {
  
  # simulate data using given simulation
  if (scenario == 1) {
    data <- simulate_scenario1(random_seed=datagen_seed,
                               n_obs=n_obs,
                               p=p,
                               k=k,
                               TT=TT)
  }
  else if (scenario == 2) {
    data <- simulate_scenario2(random_seed=datagen_seed,
                               n_obs=n_obs,
                               p=p,
                               k=k,
                               TT=TT)
  }
  else if (scenario == 3) {
    data <- simulate_scenario3(random_seed=datagen_seed,
                               n_obs=n_obs,
                               p=p,
                               TT=TT)
  }
  else {
    print("'scenario' must be set to one of 1, 2, or 3 in run_sim_LowFR function")
    return(NULL)
  }
  
  # fit LowFR model to data
  post_samples <- fit_LowFR(y_obs=data$y_obs,
                   X_obs=data$X_obs,
                   p=p,
                   k=k_model,
                   TT=TT,
                   output=parameter_output,
                   burnin=burnin,
                   samples=samples,
                   chains=chains,
                   random_seed=stan_seed)
  
  # summarize accuracy
  accuracy_summary <- summarize_accuracy_LowFR(data=data,
                                               fit=post_samples)
  
  # summarize cumulative effects
  cumulative_effect_summary <- get_cumulative_accuracy(data=data,
                                                       post_samples=post_samples,
                                                       model="LowFR")
  
  # summarize sampling diagnostics
  diagnostic_summary <- summarize_mixing_LowFR(fit=post_samples,
                                        num_samples=samples,
                                        num_chains=chains)
  
  # combine output
  result <- list()
  result_names <- c()
  if ("post_samples" %in% simulation_output) {
    result <- c(result, list(post_samples))
    result_names <- c(result_names, "post_samples")
  }
  if ("accuracy_summary" %in% simulation_output) {
    result <- c(result, list(accuracy_summary, cumulative_effect_summary))
    result_names <- c(result_names, "accuracy_summary", "cumulative_effect_summary")
  }
  if ("diagnostic_summary" %in% simulation_output) {
    result <- c(result, list(diagnostic_summary))
    result_names <- c(result_names, "diagnostic_summary")
  }
  names(result) <- result_names
  
  # save to given filepath if specified
  if (save_output == TRUE) {
    file_name <- paste0(path, "LowFR_scenario", scenario, "_seed", datagen_seed, ".rds")
    saveRDS(result, file=file_name)
  }
  
  # return output
  return(result)
}



