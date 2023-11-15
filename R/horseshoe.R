# Functions to fit and summarize horseshoe regression models

library(tidyverse)
library(horseshoe)
source("R/helper.R")
source("R/sim_data.R")

###############################################################################
############# Fit  model and return specified output ####################
###############################################################################
fit_horseshoe <- function(y_obs, X_obs,
                    burnin=1000, samples=5000,
                    random_seed=1234) {
  
  # set random seed for reproducibility
  set.seed(random_seed)
  
  # detect pT and n_obs
  pT <- ncol(X_obs)
  n_obs <- nrow(X_obs)
  
  # expand X_obs to include an intercept column and all pairwise interactions
  X_wide <- X_obs
  for (i in 1:(pT)) {
    for (j in 1:i) {
      X_wide <- cbind(X_wide, X_wide[,i] * X_wide[,j])
    }
  }
  X_wide <- cbind(rep(1,n_obs), X_wide)
  
  # fit model
  fit <- horseshoe(y_obs,
                   X_wide,
                   method.tau="halfCauchy",
                   method.sigma="Jeffreys",
                   burn=burnin,
                   nmc=samples)
  
  # reformat data into alpha_0, alpha, and Gamma
  alpha_0 <- fit$BetaSamples[1,]
  mains <- fit$BetaSamples[(2:(pT+1)),]
  ints <- fit$BetaSamples[((pT+2):(1+pT + pT*(pT+1)/2)),]
  
  # rearrange main effects
  mains <- aperm(mains, c(2,1))
  
  # rearrange interactions
  int_mat <- array(dim=c(samples, pT, pT))
  curr_index <- 1
  for (i in 1:pT) {
    for (j in 1:i) {
      int_mat[,i,j] <- ints[curr_index,] / 2
      int_mat[,j,i] <- ints[curr_index,] / 2
      curr_index <- curr_index + 1
    }
  }
  
  # return output
  output <- list(alpha_0, mains, int_mat)
  names(output) <- c("alpha_0", "alpha", "Gamma")
  return(output)
}


###############################################################################
############## Calculate main and interaction MSE of a fitted model ###########
###############################################################################
get_mse_horseshoe <- function(data, fit) {
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
get_TP_TN_rates_horseshoe <- function(data, fit) {
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
get_coverage_horseshoe <- function(data, fit) {
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
summarize_accuracy_horseshoe <- function(data, fit) {
  mse <- get_mse_horseshoe(data, fit)
  coverage <- get_coverage_horseshoe(data, fit)
  TP_TN_rate <- get_TP_TN_rates_horseshoe(data, fit)
  output <- c(mse, coverage, TP_TN_rate)
  return(output)
}

###############################################################################
############ Summarize main and interaction effective sample sizes ############
###############################################################################
get_ESS_horseshoe <- function(fit) {
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
summarize_mixing_horseshoe <- function(fit) {
  ESS <- get_ESS_horseshoe(fit)
  output <- c(ESS)
  return(output)
}


###############################################################################
########################### Run a full horseshoe simulation #######################
###############################################################################
run_sim_horseshoe <- function(scenario,
                        n_obs=200,
                        p=10,
                        k=5,
                        TT=3,
                        save_output=FALSE,
                        path="",
                        datagen_seed=1234,
                        model_seed=1234,
                        burnin=1000,
                        samples=5000,
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
    print("'scenario' must be set to one of 1, 2, or 3 in run_sim_horseshoe function")
    return(NULL)
  }
  
  # fit horseshoe model to data
  post_samples <- fit_horseshoe(y_obs=data$y_obs,
                          X_obs=data$X_obs,
                          burnin=burnin,
                          samples=samples,
                          random_seed=model_seed)
  
  # summarize accuracy
  accuracy_summary <- summarize_accuracy_horseshoe(data=data,
                                             fit=post_samples)
  
  # summarize cumulative effects
  cumulative_effect_summary <- get_cumulative_accuracy(data=data,
                                                       post_samples=post_samples,
                                                       model="horseshoe")
  
  # summarize sampling diagnostics
  diagnostic_summary <- summarize_mixing_horseshoe(fit=post_samples)
  
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
    file_name <- paste0(path, "horseshoe_scenario", scenario, "_seed", datagen_seed, ".rds")
    saveRDS(result, file=file_name)
  }
  
  # return output
  return(result)
}


