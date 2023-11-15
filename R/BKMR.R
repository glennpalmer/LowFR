# Functions to fit and summarize BKMR models

library(tidyverse)
library(bkmr)
source("R/helper.R")
source("R/sim_data.R")

###############################################################################
########################## Fit and return model ###############################
###############################################################################
fit_BKMR <- function(y_obs, X_obs,
                          samples=8000,
                          random_seed=1234) {
  
  # set random seed for reproducibility
  set.seed(random_seed)
  
  # fit model
  fit <- kmbayes(y=y_obs, Z=X_obs, iter=samples, varsel=TRUE)
  
  # return fit
  return(fit)
}

###############################################################################
###### Fit and return model -- grouped variable selection by exposure #########
###############################################################################
fit_BKMR_group <- function(y_obs, X_obs, p=10, TT=3,
                     samples=8000,
                     random_seed=1234) {
  
  # set random seed for reproducibility
  set.seed(random_seed)
  
  # specify groups
  groups <- rep(1:p, each=TT)
  
  # fit model
  fit <- kmbayes(y=y_obs, Z=X_obs, iter=samples, varsel=TRUE, groups=groups)
  
  # return fit
  return(fit)
}

###############################################################################
################ Summarize "main" and "interaction" effects ###################
###############################################################################
get_coefficients_BKMR_old <- function(fit) {
  
  # calculate data dimensions
  n_obs <- nrow(fit$Z)
  pT <- ncol(fit$Z)
  n_samples <- nrow(fit$beta)
  
  ########################## Main Effects ###########################
  
  # one more try at doing the whole thing at once
  Z_C <- matrix(nrow=(pT*(pT+1)/2), ncol=pT)
  Z_D <- matrix(nrow=(pT*(pT+1)/2), ncol=pT)
  curr_index <- 1
  same_indices <- c()
  for (i in 1:pT) {
    for (j in 1:i) {
      Z_C[curr_index,] <- rep(0,pT)
      Z_D[curr_index,] <- rep(0,pT)
      Z_C[curr_index,i] <- 0.5
      Z_D[curr_index,i] <- -0.5
      if (i != j) {
        Z_C[curr_index,j] <- 1
        Z_D[curr_index,j] <- 1
        same_indices <- c(same_indices, curr_index)
      }
      curr_index <- curr_index + 1
    }
  }
  
  # get predictions for these points -- this step is super slow...
  start_time <- Sys.time()
  pred_vals <- SamplePred(fit=fit, Znew=rbind(Z_C, Z_D), Xnew=cbind(0), sel=(floor(n_samples/2)+1):n_samples)
  end_time <- Sys.time()
  print(paste0("Took ", end_time - start_time, " hours to SamplePred for all 5000 posterior samples"))
  
  # initialize storage
  alpha <- matrix(nrow=floor(n_samples/2), ncol=pT)
  Gamma <- array(dim=c(floor(n_samples/2), pT, pT))
  
  # get main and interaction effects
  curr_indexC <- 1
  curr_indexD <- pT*(pT+1)/2 + 1
  for (i in 1:pT) {
    for (j in 1:i) {
      if (i == j) {
        alpha[,i] <- pred_vals[,curr_indexC] - pred_vals[,curr_indexD]
        Gamma[,i,j] <- (pred_vals[,curr_indexC] + pred_vals[,curr_indexD]) * 2
        curr_indexC <- curr_indexC + 1
        curr_indexD <- curr_indexD + 1
      }
      else {
        Gamma[,i,j] <- ((pred_vals[,curr_indexC] - pred_vals[,curr_indexD]) -
                          (pred_vals[,same_indices[i]] - pred_vals[,(same_indices[i] + pT*(pT+1)/2)])) / 2
        Gamma[,j,i] <- Gamma[,i,j]
        curr_indexC <- curr_indexC + 1
        curr_indexD <- curr_indexD + 1
      }
    }
  }
  
  # return output 
  output <- list(alpha, Gamma)
  names(output) <- c("alpha", "Gamma")
  return(output)
}


###############################################################################
################ Summarize "main" and "interaction" effects ###################
###############################################################################
get_coefficients_BKMR <- function(fit) {
  
  # calculate data dimensions
  n_obs <- nrow(fit$Z)
  pT <- ncol(fit$Z)
  n_samples <- nrow(fit$beta)
  
  # initialize storage
  alpha <- matrix(nrow=floor(n_samples/2), ncol=pT)
  Gamma <- array(dim=c(floor(n_samples/2), pT, pT))
  
  ########################## Main and Interaction Effects #####################
  
  start_time <- Sys.time()
  # loop over main effects and interactions and compute them one at a time
  for (i in 1:pT) {
    print(paste0("Starting for i=", i, " at:"))
    print(Sys.time())
    for (j in 1:i) {
      Z_curr <- matrix(rep(0, 4*pT), nrow=4, ncol=pT)
      Z_curr[1,i] <- 0.5
      Z_curr[2,i] <- -0.5
      Z_curr[3,i] <- 0.5
      Z_curr[4,i] <- -0.5
      if (i != j) {
        Z_curr[3,j] <- 1
        Z_curr[4,j] <- 1
      }
      # sample at these new points
      pred_vals <- SamplePred(fit=fit, Znew=Z_curr, Xnew=cbind(0))#, sel=(floor(n_samples/2)+1):n_samples)
      
      # record relevant main and interaction effects
      alpha[,i] <- pred_vals[,1] - pred_vals[,2]
      if (i != j) {
        Gamma[,i,j] <- ((pred_vals[,3] - pred_vals[,4]) - alpha[i]) / 2
        Gamma[,j,i] <- Gamma[,i,j]
      }
      else {
        Gamma[,i,j] <- (pred_vals[,1] + pred_vals[,2]) * 2
      }
    }
  }
  end_time <- Sys.time()
  print(paste0("Took ", end_time - start_time, " minutes to SamplePred for all 5000 posterior samples"))
  
  # return output 
  output <- list(alpha, Gamma)
  names(output) <- c("alpha", "Gamma")
  return(output)
}


###############################################################################
############## Calculate main and interaction MSE of a fitted model ###########
###############################################################################
get_mse_BKMR <- function(data, fit) {
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
get_TP_TN_rates_BKMR <- function(data, fit) {
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
get_coverage_BKMR <- function(data, fit) {
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
summarize_accuracy_BKMR <- function(data, fit) {
  mse <- get_mse_horseshoe(data, fit)
  coverage <- get_coverage_horseshoe(data, fit)
  TP_TN_rate <- get_TP_TN_rates_horseshoe(data, fit)
  output <- c(mse, coverage, TP_TN_rate)
  return(output)
}

###############################################################################
############ Summarize main and interaction effective sample sizes ############
###############################################################################
get_ESS_BKMR <- function(fit) {
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
summarize_mixing_BKMR <- function(fit) {
  ESS <- get_ESS_horseshoe(fit)
  output <- c(ESS)
  return(output)
}


###############################################################################
########################### Run a full BKMR simulation #######################
###############################################################################
run_sim_BKMR <- function(scenario,
                              n_obs=200,
                              p=10,
                              k=5,
                              TT=3,
                              save_output=FALSE,
                              path="",
                              datagen_seed=1234,
                              model_seed=1234,
                              samples=8000,
                              verbose=TRUE,
                              simulation_output=c("model_fit",
                                                  "post_samples",
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
    print("'scenario' must be set to one of 1, 2, or 3 in run_sim_BKMR function")
    return(NULL)
  }
  
  # fit BKMR model
  fit <- fit_BKMR(y_obs=data$y_obs,
                  X_obs=data$X_obs,
                       samples=samples,
                       random_seed=model_seed)
  
  # extract posterior samples
  post_samples <- get_coefficients_BKMR(fit=fit)
  
  # summarize accuracy
  accuracy_summary <- summarize_accuracy_horseshoe(data=data,
                                                   fit=post_samples)
  
  # summarize cumulative effects
  cumulative_effect_summary <- get_cumulative_accuracy(data=data,
                                                       post_samples=post_samples,
                                                       original_fit=fit,
                                                       model="BKMR")
  
  # summarize sampling diagnostics
  diagnostic_summary <- summarize_mixing_horseshoe(fit=post_samples)
  
  # combine output
  result <- list()
  result_names <- c()
  if ("model_fit" %in% simulation_output) {
    result <- c(result, list(fit))
    result_names <- c(result_names, "model_fit")
  }
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
    file_name <- paste0(path, "BKMR_scenario", scenario, "_seed", datagen_seed, ".rds")
    saveRDS(result, file=file_name)
  }
  
  # return output
  return(result)
}

###############################################################################
######################## Run a full BKMR_group simulation #####################
###############################################################################
run_sim_BKMR_group <- function(scenario,
                         n_obs=200,
                         p=10,
                         k=5,
                         TT=3,
                         save_output=FALSE,
                         path="",
                         datagen_seed=1234,
                         model_seed=1234,
                         samples=8000,
                         verbose=TRUE,
                         simulation_output=c("model_fit",
                                             "post_samples",
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
    print("'scenario' must be set to one of 1, 2, or 3 in run_sim_BKMR function")
    return(NULL)
  }
  
  # fit BKMR model
  fit <- fit_BKMR_group(y_obs=data$y_obs,
                  X_obs=data$X_obs,
                  p=p,
                  TT=TT,
                  samples=samples,
                  random_seed=model_seed)
  
  # extract posterior samples
  post_samples <- get_coefficients_BKMR(fit=fit)
  
  # summarize accuracy
  accuracy_summary <- summarize_accuracy_horseshoe(data=data,
                                                   fit=post_samples)
  
  # summarize cumulative effects
  cumulative_effect_summary <- get_cumulative_accuracy(data=data,
                                                       post_samples=post_samples,
                                                       original_fit=fit,
                                                       model="BKMR_group")
  
  # summarize sampling diagnostics
  diagnostic_summary <- summarize_mixing_horseshoe(fit=post_samples)
  
  # combine output
  result <- list()
  result_names <- c()
  if ("model_fit" %in% simulation_output) {
    result <- c(result, list(fit))
    result_names <- c(result_names, "model_fit")
  }
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
    file_name <- paste0(path, "BKMR_group_scenario", scenario, "_seed", datagen_seed, ".rds")
    saveRDS(result, file=file_name)
  }
  
  # return output
  return(result)
}
