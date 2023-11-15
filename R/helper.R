# helper functions used across different models
library(coda)

###############################################################################
############# Calculate squared estimation error of posterior mean ############
###############################################################################
squared_error <- function(true_val, est_vec, estimator="mean") {
  # returns the squared estimation error of mean(est_vec) for true_val
  if (estimator == "mean") {
    return((mean(est_vec) - true_val)^2)
  }
  else {
    print("No valid estimator for squared error specified. Returning NA.")
    return(NA)
  }
}

###############################################################################
################### Calculate credible interval of given width ################
###############################################################################
get_CI <- function(est_vec, interval_width) {
  lower <- (1 - interval_width) / 2
  upper <- 1 - lower
  quants <- quantile(x=est_vec, probs=c(lower, upper))
  return(c(quants[[1]], quants[[2]]))
}

###############################################################################
############ Return whether or not the true value is in a given CI ############
###############################################################################
cover_bool <- function(true_val, est_vec, interval_width=0.95) {
  # returns whether the specified interval of est_vec contains true_val
  quants <- get_CI(est_vec=est_vec, interval_width=interval_width)
  if (quants[1] < true_val & quants[2] > true_val) {
    return(TRUE)
  }
  return(FALSE)
}

###############################################################################
############################ Detect TP, TN, FP, or FN #########################
###############################################################################
sens_spec <- function(true_val, est_vec, interval_width=0.95) {
  # if true_val=0 returns:
    # "TN" for 95% CI of est_vec including 0
    # "FP" for 95% CI of est_vec excluding 0
  # if true_val!=0 returns:
    # "TP" for 95% CI of est_vec excluding 0 and correct sign
    # "FN" for 95% CI of est_vec including 0 or incorrect sign
  
  # compute given quantiles
  quants <- get_CI(est_vec=est_vec, interval_width=interval_width)
  
  # check for TN or FP:
  if (true_val == 0) {
    if (quants[1] <= 0 & quants[2] >= 0) {
      return("TN")
    }
    else {
      return("FP")
    }
  }
  
  # check for TP or FN
  else {
    if (quants[1]*true_val > 0 & quants[2]*true_val > 0) {
      return("TP")
    }
    else {
      return("FN")
    }
  }
}

###############################################################################
############################ Calculate ESS of a vector ########################
###############################################################################
get_ESS <- function(est_vec) {
  samps <- mcmc(est_vec)
  return(effectiveSize(samps))
}

###############################################################################
############################ Calculate Rhat of a vector #######################
###############################################################################
get_Rhat <- function(est_vec, num_samples, num_chains) {
  # check for correct dimensionality
  if (length(est_vec) != num_samples*num_chains) {
    print("in get_Rhat function, length(est_vec) must equal num_samples*num_chains")
    return(NA)
  }
  
  # create mcmc.list object
  chain_list <- list()
  for (i in 1:num_chains) {
    chain_list <- c(chain_list, list(mcmc(est_vec[((i-1)*num_samples + 1):(i*num_samples)])))
  }
  samps <- mcmc.list(chain_list)
  return(gelman.diag(samps)$psrf[1,1])
}


###############################################################################
### Get accuracy of effects for cumulative exposure at all times -- i.e., the##
### expected increase in the outcome if an exposure increases from -1 to 1   ##
### at all three times.                                                      ##      
###############################################################################
get_cumulative_accuracy <- function(data, post_samples, model, original_fit=NULL) {
  
  # load the data
  curr_data <- data
  
  # initialize storage
  cum_se <- c()
  cum_coverage_bool <- c()
  cum_sens_spec <- c()
  
  # calculate accuracy -- LowFR, FIN, CorrQuadReg, Horseshoe
  if (model %in% c("LowFR", "FIN", "CorrQuadReg", "Horseshoe")) {
    for (j in 1:10) {
      curr_samples <- 2*(post_samples$alpha[,(j-1)*3 + 1] +
                           post_samples$alpha[,(j-1)*3 + 2] +
                           post_samples$alpha[,(j-1)*3 + 3])
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
    curr_fit <- original_fit
    
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
get_cum_samples_BKMR <- function(fit) {
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




