## Summarize simulation or other results

###############################################################################
########## Summarize single simulation for all scenarios/models ###############
###############################################################################
simulation_1234_summary <- function() {
  path <- "simulations_1234/"
  model_list <- c("LowFR", "horseshoe", "FIN", "BKMR", "BKMR_group", "CorrQuadReg")
  for (m in model_list) {
    for (scenario in 1:3) {
      # load data for given model and scenario
      curr_path <- paste0(path, m, "/scenario", scenario, "/")
      curr_files <- list.files(path=curr_path)
      curr_path <- paste0(curr_path, curr_files[1])
      curr_data <- readRDS(curr_path)
      
      # summarize the results
      
    }
  }
}

###############################################################################
######## Summarize single simulation for a given scenario and model ###########
###############################################################################
single_sim_summary <- function(model, scenario, print_mixing=FALSE) {
  # load specified data
  path <- "simulations_1234/"
  curr_path <- paste0(path, model, "/scenario", scenario, "/")
  curr_file <- list.files(path=curr_path)
  curr_path <- paste0(curr_path, curr_file[1])
  curr_data <- readRDS(curr_path)
  
  # print summary
  print(paste0("Scenario ", scenario, " ", model, " summary:"))
  print("Accuracy:")
  print(curr_data$accuracy_summary)
  if (print_mixing) {
    print("Mixing:")
    print(curr_data$diagnostic_summary)
  }
  
  return(0)
}

###############################################################################
### Summarize list of simulations for a given scenario, seed_list, and model ##
###############################################################################
multi_sim_summary <- function(model, scenario, seed_list=c(1,5,13,17,21,24,25,33,
                                                           34,37,41,49,61,69,73,
                                                           77,85,89,97,1234),
                              path="simulations_11-6/", print_mixing=FALSE) {
  
  # initialize storage for all quantities of interest
  main_mse <- c()
  int_mse <- c()
  main_coverage <- c()
  int_coverage <- c()
  main_TP_rate <- c()
  main_TN_rate <- c()
  int_TP_rate <- c()
  int_TN_rate <- c()
  if (print_mixing==TRUE) {
    main_ESS_mins <- c()
    int_ESS_mins <- c()
    main_Rhat_maxes <- c()
    int_Rhat_maxes <- c()
  }
  
  # loop over model fits with given seeds and compile the overall results
  for (seed in seed_list) {
    # get path to given model, scenario, and seed
    full_path=paste0(path, model, "/scenario", scenario, "/",
                model, "_scenario", scenario, "_seed", seed, ".rds")
    
    # load the model
    curr_model <- readRDS(full_path)
    
    # add accuracy results to vectors
    main_mse <- c(main_mse, curr_model$accuracy_summary$main_mse)
    int_mse <- c(int_mse, curr_model$accuracy_summary$int_mse)
    main_coverage <- c(main_coverage, curr_model$accuracy_summary$main_coverage)
    int_coverage <- c(int_coverage, curr_model$accuracy_summary$int_coverage)
    main_TP_rate <- c(main_TP_rate, curr_model$accuracy_summary$main_TP_rate)
    int_TP_rate <- c(int_TP_rate, curr_model$accuracy_summary$int_TP_rate)
    main_TN_rate <- c(main_TN_rate, curr_model$accuracy_summary$main_TN_rate)
    int_TN_rate <- c(int_TN_rate, curr_model$accuracy_summary$int_TN_rate)
    
    # add mixing results if requested
    if (print_mixing==TRUE) {
      main_ESS_mins <- c(main_ESS_mins, curr_model$diagnostic_summary$main_ESS_min)
      int_ESS_mins <- c(int_ESS_mins, curr_model$diagnostic_summary$int_ESS_min)
      if (model %in% c("LowFR", "CorrQuadReg")) {
        main_Rhat_maxes <- c(main_Rhat_maxes, curr_model$diagnostic_summary$main_Rhat_max)
        int_Rhat_maxes <- c(int_Rhat_maxes, curr_model$diagnostic_summary$int_Rhat_max)
      }
    }
  }
  
  # print summaries
  print(paste0("Averaged accuracy summary for ", model, " Scenario ", scenario, ":"))
  print("Main MSE:")
  print(mean(main_mse))
  print("Interaction MSE:")
  print(mean(int_mse))
  print("Main coverage:")
  print(mean(main_coverage))
  print("Int coverage:")
  print(mean(int_coverage))
  print("Main TP rate:")
  print(mean(main_TP_rate))
  print("Main TN rate:")
  print(mean(main_TN_rate))
  print("Int TP rate:")
  print(mean(int_TP_rate))
  print("Int TN rate:")
  print(mean(int_TN_rate))
  
  # print mixing summaries if requested
  if (print_mixing) {
    print(paste0("Overall mixing summary for ", model, " Scenario ", scenario, ":"))
    print("Min main ESS's across models:")
    print(main_ESS_mins)
    print("Min int ESS's across models:")
    print(int_ESS_mins)
    print("Overall min main ESS:")
    print(min(main_ESS_mins))
    print("Overall min int ESS:")
    print(min(int_ESS_mins))
    if (model %in% c("LowFR", "CorrQuadReg")) {
      print("Overall max main Rhat:")
      print(max(main_Rhat_maxes))
      print("Overall max int Rhat:")
      print(max(int_Rhat_maxes))
    }
  }
}


###############################################################################
### Get accuracy of effects for cumulative exposure at all times -- i.e., the##
### expected increase in the outcome if an exposure increases from -1 to 1   ##
### at all three times.                                                      ##      
###############################################################################
get_cumulative_accuracy_postrun <- function(model, scenario, seed,
                                    path="simulations_11-6/") {
  
  # get the path for the given model results
  full_path=paste0(path, model, "/scenario", scenario, "/",
                   model, "_scenario", scenario, "_seed", seed, ".rds")
  
  # load the model
  curr_model <- readRDS(full_path)
  
  # get the dataset and true coefficients
  if (scenario == 1) {
    curr_data <- simulate_scenario1(random_seed=seed)
  }
  else if (scenario == 2) {
    curr_data <- simulate_scenario2(random_seed = seed)
  }
  else if (scenario == 3) {
    curr_data <- simulate_scenario3(random_seed = seed)
  }
  else {
    print("Need to specify scenario as one of 1, 2, or 3.")
    return(NULL)
  }
  
  # initialize storage
  cum_se <- c()
  cum_coverage_bool <- c()
  cum_sens_spec <- c()
  
  # calculate accuracy -- LowFR, FIN, CorrQuadReg, Horseshoe
  if (model %in% c("LowFR", "FIN", "CorrQuadReg", "Horseshoe")) {
    for (j in 1:10) {
      curr_samples <- 2*(curr_model$post_samples$alpha[,(j-1)*3 + 1] +
                               curr_model$post_samples$alpha[,(j-1)*3 + 2] +
                               curr_model$post_samples$alpha[,(j-1)*3 + 3])
      curr_true_val <- 2*(curr_data$alpha[(j-1)*3 + 1] +
                            curr_data$alpha[(j-1)*3 + 2] +
                            curr_data$alpha[(j-1)*3 + 3])
      
      # calculate accuracy summaries
      curr_mse <- squared_error(true_val=curr_true_val,
                                est_vec=curr_samples)
      
      curr_coverage <- cover_bool(true_val = curr_true_val,
                                  est_vec=curr_samples)
      
      curr_sens_spec <- sens_spec(true_val=curr_true_val,
                                  est_vec=curr_samples)
      
      # append to storage
      cum_se <- c(cum_se, curr_mse)
      cum_coverage_bool <- c(cum_coverage_bool, curr_coverage)
      cum_sens_spec <- c(cum_sens_spec, curr_sens_spec)
    }
  }
  
  # calculate accuracy -- BKMR and BKMR_group
  if (model %in% c("BKMR", "BKMR_group")) {
    curr_fit <- curr_model$model_fit
    
    cum_samples <- get_cum_samples_BKMR(curr_fit)
    
    for (j in 1:10) {
      curr_samples <- cum_samples[,j]
      
      curr_true_val <- 2*(curr_data$alpha[(j-1)*3 + 1] +
                            curr_data$alpha[(j-1)*3 + 2] +
                            curr_data$alpha[(j-1)*3 + 3])
      
      # calculate accuracy summaries
      curr_mse <- squared_error(true_val=curr_true_val,
                                est_vec=curr_samples)
      
      curr_coverage <- cover_bool(true_val = curr_true_val,
                                  est_vec=curr_samples)
      
      curr_sens_spec <- sens_spec(true_val=curr_true_val,
                                  est_vec=curr_samples)
      
      # append to storage
      cum_se <- c(cum_se, curr_mse)
      cum_coverage_bool <- c(cum_coverage_bool, curr_coverage)
      cum_sens_spec <- c(cum_sens_spec, curr_sens_spec)
    }
  }
  
  
  # calculate summaries
  cum_mse <- mean(cum_se)
  cum_coverage <- mean(cum_coverage_bool)
  cum_TP_rate <- sum(cum_sens_spec=="TP") / (sum(cum_sens_spec %in% c("TP", "FN")))
  cum_TN_rate <- sum(cum_sens_spec=="TN") / (sum(cum_sens_spec %in% c("TN", "FP")))
  
  # return summaries as a list
  output <- list(cum_mse, cum_coverage, cum_TP_rate, cum_TN_rate)
  names(output) <- c("cum_mse", "cum_coverage", "cum_TP_rate", "cum_TN_rate")
  return(output)
}


###############################################################################
############ Get cumulative effect samples for BKMR or BKMR_group #############
###############################################################################
get_cum_samples_BKMR_postrun <- function(fit) {
  # generate new points to sample
  Znew <- matrix(rep(0, 300), nrow=10, ncol=30)
  for (i in 1:10) {
    for (j in ((i-1)*3 + 1):(i*3)) {
      Znew[i,j] <- 1
    }
  }
  Znew <- rbind(Znew, -1*Znew)

  # sample new points
  n_samples <- fit$iter
  pred_vals <- SamplePred(fit=fit, Znew=Znew, Xnew=cbind(0), sel=(floor(n_samples/2)+1):n_samples)
  
  # summarize differences -- i.e., cumulative exposure effects for increase
  # of each exposure from -1 to 1
  alpha_mat <- pred_vals[,1:10] - pred_vals[,11:20]
  
  return(alpha_mat)
}


###############################################################################
#################### Get averaged cumulative effects ##########################     
###############################################################################
multi_sim_cum_effects <- function(model, scenario,
                                  seed_list=c(1,5,13,17,21,24,25,33,
                                              34,37,41,49,61,69,73,
                                              77,85,89,97,1234),
                                  path="simulations_11-6/") {
  
  # initialize storage
  mse <- c()
  coverage <- c()
  TP_rate <- c()
  TN_rate <- c()
  
  # add results for each random seed
  for (seed in seed_list) {
    curr_results <- get_cumulative_accuracy(model=model, scenario=scenario,
                                            seed=seed, path=path)
    
    mse <- c(mse, curr_results$cum_mse)
    coverage <- c(coverage, curr_results$cum_coverage)
    TP_rate <- c(TP_rate, curr_results$cum_TP_rate)
    TN_rate <- c(TN_rate, curr_results$cum_TN_rate)
  }
  
  mean_mse <- mean(mse)
  mean_coverage <- mean(coverage)
  mean_TP_rate <- mean(TP_rate)
  mean_TN_rate <- mean(TN_rate)
  
  # return results
  output <- list(mean_mse, mean_coverage, mean_TP_rate, mean_TN_rate)
  names(output) <- c("cumulative_exp_mse", "cumulative_exp_coverage",
                     "cumulative_exp_TP_rate", "cumulative_exp_TN_rate")
  return(output)
}
  



