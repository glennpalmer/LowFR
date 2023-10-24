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





