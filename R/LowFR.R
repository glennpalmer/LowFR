# Functions to fit and summarize LowFR models as described in "Low-Rank
# Longitudinal Factor Regression" paper

library(tidyverse)
library(rstan)

###############################################################################
############# Fit a Stan model and return specified output ####################
###############################################################################
fit_LowFR <- function(y_obs, X_obs, p, k, TT,
                      output="stan_model",
                      burnin=1000, samples=1000, chains=4) {
  
  # Note: output can be either "stan_model", "all", or a list of parameters 
  # for which to return posterior samples. The list can include items from
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
  
  
}