# Functions to fit and summarize FIN models

library(tidyverse)
source("R/helper.R")
source("R/sim_data.R")
source("R/infinitefactor/R/interactionDL.R")

###############################################################################
############# Fit  model and return specified output ####################
###############################################################################
fit_FIN <- function(y_obs, X_obs, k=NULL,
                      burnin=5000, samples=10000,
                      verbose=TRUE,
                      random_seed=1234,
                      output=c("coefSamples")) {
  
  # choose k using SVD if needed
  if (is.null(k)) {
    k <- k_svd_FIN(X_obs=X_obs)
  }
  
  # set random seed for reproducibility
  set.seed(random_seed)
  
  # fit model
  fit <- interactionDL(y=y_obs,
                       X=X_obs,
                       k=k,
                       nrun=samples,
                       burn=burnin,
                       verbose=verbose,
                       output=output)
  
  
  # return samples
  names(fit) <- c("theta", "Omega", "alpha_0", "alpha", "Gamma")
  fit$theta <- aperm(fit$theta, c(2,1))
  fit$Omega <- aperm(fit$Omega, c(3,1,2))
  fit$alpha <- aperm(fit$alpha, c(2,1))
  fit$Gamma <- aperm(fit$Gamma, c(3,1,2))
  return(fit)
}

###############################################################################
######################### Choose k based on SVD ###############################
###############################################################################
k_svd_FIN <- function(X_obs) {
  
  # get p from dimension of X_obs
  p <- ncol(X_obs)
  
  # do SVD to get singular values
  sing_vals <- svd(X_obs)$d
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
get_mse_FIN <- function(data, fit) {
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
get_TP_TN_rates_FIN <- function(data, fit) {
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
get_coverage_FIN <- function(data, fit) {
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
summarize_accuracy_FIN <- function(data, fit) {
  mse <- get_mse_FIN(data, fit)
  coverage <- get_coverage_FIN(data, fit)
  TP_TN_rate <- get_TP_TN_rates_FIN(data, fit)
  output <- c(mse, coverage, TP_TN_rate)
  return(output)
}

###############################################################################
############ Summarize main and interaction effective sample sizes ############
###############################################################################
get_ESS_FIN <- function(fit) {
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
########################### Summarize mixing diagnostics ######################
###############################################################################
summarize_mixing_FIN <- function(fit) {
  ESS <- get_ESS_FIN(fit)
  output <- c(ESS)
  return(output)
}


###############################################################################
########################### Run a full FIN simulation #######################
###############################################################################
run_sim_FIN <- function(scenario,
                        n_obs=200,
                        p=10,
                        k=5,
                        TT=3,
                          save_output=FALSE,
                          path="",
                          datagen_seed=1234,
                          model_seed=1234,
                          k_model=NULL,
                          parameter_output=c("coefSamples"),
                          burnin=5000,
                          samples=10000,
                          verbose=TRUE,
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
    print("'scenario' must be set to one of 1, 2, or 3 in run_sim_FIN function")
    return(NULL)
  }
  
  # fit FIN model to data
  post_samples <- fit_FIN(y_obs=data$y_obs,
                          X_obs=data$X_obs,
                          k=k_model,
                          burnin=burnin,
                          samples=samples,
                          verbose=verbose,
                          random_seed=model_seed,
                          output=parameter_output)
  
  # summarize accuracy
  accuracy_summary <- summarize_accuracy_FIN(data=data,
                                               fit=post_samples)
  
  # summarize cumulative effects
  cumulative_effect_summary <- get_cumulative_accuracy(data=data,
                                                       post_samples=post_samples,
                                                       model="FIN")
  
  # summarize sampling diagnostics
  diagnostic_summary <- summarize_mixing_FIN(fit=post_samples)
  
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
    file_name <- paste0(path, "FIN_scenario", scenario, "_seed", datagen_seed, ".rds")
    saveRDS(result, file=file_name)
  }
  
  # return output
  return(result)
}
