# Functions to fit and summarize LowFR models as described in "Low-Rank
# Longitudinal Factor Regression" paper

library(tidyverse)
library(rstan)

###############################################################################
############# Fit a Stan model and return specified output ####################
###############################################################################
fit_LowFR <- function(y_obs, X_obs, p, k=NULL, TT,
                      output="stan_fit",
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
  
  # choose k using SVD if needed
  if (is.null(k)) {
    # add function call here to compute k using SVD method
  }
  
  # compile stan model
  m <- stan_model("../Stan/LowFR.stan")
  
  # fit model
  options(mc.cores = parallel::detectCores())
  fit <- sampling(m,
                  list(N=n_obs,
                       p=p,
                       k=k,
                       H=min(p,TT),
                       X=X_obs,
                       y=y_obs,
                       chains=chains,
                       iter=burnin+samples,
                       warmup=burnin,
                       seed=random_seed))
  
  # return the stan fit if specified as output
  if (output == "stan_fit") {
    return(fit)
  }
  
  # otherwise, extract posterior samples and return list of specified parameters
  post_samples <- extract(fit)
  results <- list()
  result_names <- c()
  if ("alpha_0" %in% output | output == "all") {
    results <- list.append(results, post_samples$alpha_0)
    result_names <- c(result_names, "alpha_0")
  }
  if ("alpha" %in% output | output == "all") {
    results <- list.append(results, post_samples$alpha)
    result_names <- c(result_names, "alpha")
  }
  if ("Gamma" %in% output | output == "all") {
    results <- list.append(results, post_samples$Gamma)
    result_names <- c(result_names, "Gamma")
  }
  if ("sigma2" %in% output | output == "all") {
    results <- list.append(results, post_samples$sigma2)
    result_names <- c(result_names, "sigma2")
  }
  if ("Sigma" %in% output | output == "all") {
    results <- list.append(results, post_samples$Sigma)
    result_names <- c(result_names, "Sigma")
  }
  if ("phi" %in% output | output == "all") {
    results <- list.append(results, post_samples$phi)
    result_names <- c(result_names, "phi")
  }
  if ("tau" %in% output | output == "all") {
    results <- list.append(results, post_samples$tau)
    result_names <- c(result_names, "tau")
  }
  if ("theta" %in% output | output == "all") {
    results <- list.append(results, post_samples$theta)
    result_names <- c(result_names, "theta")
  }
  if ("Omega" %in% output | output == "all") {
    results <- list.append(results, post_samples$Omega)
    result_names <- c(result_names, "Omega")
  }
  names(results) <- result_names
  return(results)
}



