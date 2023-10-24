# Functions to fit and summarize LowFR models as described in "Low-Rank
# Longitudinal Factor Regression" paper

library(tidyverse)
library(rstan)

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
                       H=min(p,TT),
                       X=X_obs,
                       y=y_obs),
                  chains=chains,
                  iter=burnin+samples,
                  warmup=burnin,
                  seed=random_seed)
  
  # return the stan fit if specified as output
  if (output == "stan_fit") {
    return(fit)
  }
  
  # otherwise, extract posterior samples and return list of specified parameters
  post_samples <- extract(fit)
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







