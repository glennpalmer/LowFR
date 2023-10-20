# Functions to generate data needed for simulations in "Low-Rank
# Longitudinal Factor Regression"

# import necessary libraries
library(dirmult)
library(mvtnorm)
library(psych)

###############################################################################
########################### Scenario 1 (Section 4) ############################
########## (Data generated from LowFR model with rank-1 theta) ################
###############################################################################
simulate_scenario1 <- function(set_seed = TRUE,
                               random_seed=1234,
                               n_obs=200,
                               p=10,
                               k=5,
                               TT=3, 
                               phi=0.5,
                               sigma2=1,
                               Sigma_diag=rep(0.25,p),
                               mu=0,
                               nonzero_main=2,
                               nonzero_int=3,
                               Lambda_sparsity=1) {
  
  # set random seed
  if (set_seed) {
    set.seed(random_seed)
  }
  
  # create Sigma matrix from vector input
  Sigma = diag(Sigma_diag)
  
  # specify covariance matrix for eta_i's over time
  Phi_mat <- matrix(data = rep(1, TT*TT), nrow=TT)
  for (i in 1:TT) {
    for (j in 1:TT) {
      if (j != i) {
        Phi_mat[i,j] <- phi
      }
    }
  }
  
  # specify beta and omega (main effects)
  beta1 <- rep(0, k)
  for (j in 1:nonzero_main) {
    beta1[j] <- runif(1,1,2)
    if (runif(1) > 0.5) {
      beta1[j] <- beta1[j] * -1
    }
  }
  omega1 <- rdirichlet(n=1, alpha=rep(2,TT))[1,]
  
  # specify B and W (interactions)
  B <- matrix(rep(0, k^2), nrow=k)
  nonzero_B <- sample(1:(k^2), nonzero_int, replace=FALSE)
  for (j in 1:nonzero_int) {
    B[ceiling(nonzero_B[j] / k), nonzero_B[j] %% k] <- runif(1,1,2)
    if (runif(1) > 0.5) {
      B[ceiling(nonzero_B[j] / k), nonzero_B[j] %% k] <- B[ceiling(nonzero_B[j] / k), nonzero_B[j] %% k] * -1
    }
  }
  W <- matrix(rdirichlet(n=1, alpha=rep(2,TT*TT))[1,], nrow=TT)
  
  # specify vector of main-effect coefficients based on the above to generate y's
  theta <- as.vector(omega1 %*% t(beta1))
  
  # specify matrix for quadratic regression coefficients
  Omega <- kronecker(B, W)
  
  # make Omega symmetric for easier interpretability
  Omega <- (t(Omega) + Omega) / 2
  
  # specify factor loadings matrix Lambda
  Lambda = matrix(rep(0,p*k), nrow=p)
  for (i in 1:p) {
    for (j in 1:k) {
      if (runif(1) < Lambda_sparsity) {
        Lambda[i,j] <- rnorm(1)
      }
    }
  }
  Lambda_kron_I = kronecker(Lambda, diag(TT))
  
  # specify Sigma kronecker matrix (covariance matrix for x_i)
  Sigma_kron_phi = kronecker(Sigma, Phi_mat)
  
  ################################## Simulate Data ############################
  
  # generate etas
  Eta_obs <- matrix(data=rep(NA, n_obs * k * TT), nrow=n_obs)
  for (i in 1:n_obs) {
    for (j in 1:k) {
      Eta_obs[i, (1 + (j-1) * TT):(j * TT)] <- rmvnorm(n=1, sigma = Phi_mat)
    }
  }
  # generate ys using etas
  y_obs <- rep(NA, n_obs)
  for (i in 1:n_obs) {
    y_obs[i] <- mu + (t(theta) %*% Eta_obs[i,])[1] + 
      (t(Eta_obs[i,]) %*% Omega %*% Eta_obs[i,])[1] +
      rnorm(n=1, mean=0, sd=sqrt(sigma2))
  }
  
  # generate xs using etas
  X_obs = matrix(nrow=n_obs, ncol=p*TT)
  for (i in 1:n_obs) {
    X_obs[i,] = rmvnorm(n=1, mean=Lambda_kron_I %*% Eta_obs[i,], sigma=Sigma_kron_phi)
  }
  
  # calculate y_i | x_i using formula from paper
  # define needed matrices
  V_regression <- kronecker(solve(t(Lambda) %*% solve(Sigma) %*% Lambda + diag(k)), Phi_mat)
  A_regression <- V_regression %*% kronecker(t(Lambda) %*% solve(Sigma), solve(Phi_mat))
  
  # calculate induced coefficients
  alpha_0 <- tr(Omega %*% V_regression)
  alpha <- as.vector(t(theta) %*% A_regression)
  Gamma <- t(A_regression) %*% Omega %*% A_regression
  
  # return results
  output = list(y_obs, X_obs, alpha_0, alpha, Gamma, theta, Omega)
  names(output) = c("y_obs", "X_obs", "alpha_0", "alpha", "Gamma", "theta", "Omega")
  return(output)
}

###############################################################################
########################### Scenario 2 (Section 4) ############################
########## (Data generated from LowFR model with rank-2 theta) ################
###############################################################################
simulate_scenario2 <- function(set_seed = TRUE,
                               random_seed=1234,
                               n_obs=200,
                               p=10,
                               k=5,
                               TT=3, 
                               phi=0.5,
                               sigma2=1,
                               Sigma_diag=rep(0.25,p),
                               mu=0,
                               nonzero_main=2,
                               nonzero_int=3,
                               Lambda_sparsity=1) {
  
  # set random seed
  if (set_seed) {
    set.seed(random_seed)
  }
  
  # create Sigma matrix from vector input
  Sigma = diag(Sigma_diag)
  
  # specify covariance matrix for eta_i's over time
  Phi_mat <- matrix(data = rep(1, TT*TT), nrow=TT)
  for (i in 1:TT) {
    for (j in 1:TT) {
      if (j != i) {
        Phi_mat[i,j] <- phi
      }
    }
  }
  
  # specify beta and omega (main effects)
  beta1 <- rep(0, k)
  beta2 <- rep(0, k)
  nonzero_beta1 <- sample(1:k, nonzero_main, replace=FALSE)
  nonzero_beta2 <- sample(1:k, nonzero_main, replace=FALSE)
  for (j in 1:nonzero_main) {
    beta1[nonzero_beta1[j]] <- runif(1,1,2)
    if (runif(1) > 0.5) {
      beta1[nonzero_beta1[j]] <- beta1[nonzero_beta1[j]] * -1
    }
    beta2[nonzero_beta2[j]] <- runif(1,1,2)
    if (runif(1) > 0.5) {
      beta2[nonzero_beta2[j]] <- beta2[nonzero_beta2[j]] * -1
    }
  }
  omega1 <- rdirichlet(n=1, alpha=rep(2,TT))[1,]
  omega2 <- rdirichlet(n=1, alpha=rep(2,TT))[1,]
  
  # specify B and W (interactions)
  B <- matrix(rep(0, k^2), nrow=k)
  nonzero_B <- sample(1:(k^2), nonzero_int, replace=FALSE)
  for (j in 1:nonzero_int) {
    B[ceiling(nonzero_B[j] / k), nonzero_B[j] %% k] <- runif(1,1,2)
    if (runif(1) > 0.5) {
      B[ceiling(nonzero_B[j] / k), nonzero_B[j] %% k] <- B[ceiling(nonzero_B[j] / k), nonzero_B[j] %% k] * -1
    }
  }
  W <- matrix(rdirichlet(n=1, alpha=rep(2,TT*TT))[1,], nrow=TT)
  
  # specify vector of main-effect coefficients based on the above to generate y's
  theta <- as.vector(omega1 %*% t(beta1) + omega2 %*% t(beta2))
  
  # specify matrix for quadratic regression coefficients
  Omega <- kronecker(B, W)
  
  # make Omega symmetric for easier interpretability
  Omega <- (t(Omega) + Omega) / 2
  
  # specify factor loadings matrix Lambda
  Lambda = matrix(rep(0,p*k), nrow=p)
  for (i in 1:p) {
    for (j in 1:k) {
      if (runif(1) < Lambda_sparsity) {
        Lambda[i,j] <- rnorm(1)
      }
    }
  }
  Lambda_kron_I = kronecker(Lambda, diag(TT))
  
  # specify Sigma kronecker matrix (covariance matrix for x_i)
  Sigma_kron_phi = kronecker(Sigma, Phi_mat)
  
  ################################## Simulate Data ############################
  
  # generate etas
  Eta_obs <- matrix(data=rep(NA, n_obs * k * TT), nrow=n_obs)
  for (i in 1:n_obs) {
    for (j in 1:k) {
      Eta_obs[i, (1 + (j-1) * TT):(j * TT)] <- rmvnorm(n=1, sigma = Phi_mat)
    }
  }
  # generate ys using etas
  y_obs <- rep(NA, n_obs)
  for (i in 1:n_obs) {
    y_obs[i] <- mu + (t(theta) %*% Eta_obs[i,])[1] + 
      (t(Eta_obs[i,]) %*% Omega %*% Eta_obs[i,])[1] +
      rnorm(n=1, mean=0, sd=sqrt(sigma2))
  }
  
  # generate xs using etas
  X_obs = matrix(nrow=n_obs, ncol=p*TT)
  for (i in 1:n_obs) {
    X_obs[i,] = rmvnorm(n=1, mean=Lambda_kron_I %*% Eta_obs[i,], sigma=Sigma_kron_phi)
  }
  
  # calculate y_i | x_i using formula from paper
  # define needed matrices
  V_regression <- kronecker(solve(t(Lambda) %*% solve(Sigma) %*% Lambda + diag(k)), Phi_mat)
  A_regression <- V_regression %*% kronecker(t(Lambda) %*% solve(Sigma), solve(Phi_mat))
  
  # calculate actual coefficients
  alpha_0 <- tr(Omega %*% V_regression)
  alpha <- as.vector(t(theta) %*% A_regression)
  Gamma <- t(A_regression) %*% Omega %*% A_regression
  
  # return results
  output = list(y_obs, X_obs, alpha_0, alpha, Gamma, theta, Omega)
  names(output) = c("y_obs", "X_obs", "alpha_0", "alpha", "Gamma", "theta", "Omega")
  return(output)
}

###############################################################################
########################### Scenario 3 (Section 4) ############################
########## (Data generated from misspecified model with kronecker #############
########## covariance for x_i and coefficients within times are   #############
########## generated as constant plus noise (so full rank).)      #############
###############################################################################
simulate_scenario3 <- function(set_seed = TRUE,
                               random_seed=1234,
                               n_obs=200,
                               p=10,
                               k=5,
                               TT=3, 
                               phi1=0.7,
                               phi2=0.7,
                               sigma2=1,
                               Sigma_diag=rep(0.25,p),
                               mu=0,
                               nonzero_main=4,
                               nonzero_int=10) {
  
  # set random seed
  if (set_seed) {
    set.seed(random_seed)
  }
  
  # create Sigma matrix from vector input
  Sigma = diag(Sigma_diag)
  
  # specify covariance for exposures over time
  Phi_mat1 <- matrix(data = rep(1, TT*TT), nrow=TT)
  for (i in 1:TT) {
    for (j in 1:TT) {
      if (j != i) {
        Phi_mat1[i,j] <- phi1
      }
    }
  }
  
  # specify covariance across exposures within time
  Phi_mat2 <- matrix(data = rep(1, p*p), nrow=p)
  for (i in 1:p) {
    for (j in 1:p) {
      if (j != i) {
        Phi_mat2[i,j] <- phi2
      }
    }
  }
  
  # specify covariance of x
  X_cov <- kronecker(Phi_mat2, Phi_mat1)
  
  # specify beta (mean for each exposure main effects)
  # and alpha (actual main effects)
  alpha <- rep(0, p*TT)
  beta1 <- rep(0, p)
  nonzero_beta1 <- sample(1:p, nonzero_main, replace=FALSE)
  for (j in 1:nonzero_main) {
    beta1[nonzero_beta1[j]] <- runif(1,0.2,0.4)
    if (runif(1) > 0.5) {
      beta1[nonzero_beta1[j]] <- beta1[nonzero_beta1[j]] * -1
    }
    for (t in 1:TT) {
      alpha[TT * (nonzero_beta1[j] - 1) + t] <- rnorm(n=1, mean=beta1[nonzero_beta1[j]], sd=0.05)
    }
  }
  
  # specify B (mean for each interaction)
  B <- matrix(rep(0, p^2), nrow=p)
  nonzero_B <- sample(1:(p^2), nonzero_int, replace=FALSE)
  for (j in 1:nonzero_int) {
    B[ceiling(nonzero_B[j] / p), nonzero_B[j] %% p] <- runif(1,0.05,0.15)
    if (runif(1) > 0.5) {
      B[ceiling(nonzero_B[j] / p), nonzero_B[j] %% p] <- B[ceiling(nonzero_B[j] / p), nonzero_B[j] %% p] * -1
    }
  }
  
  # specify Gamma (interactions) using B
  Gamma <- matrix(rep(0, p*p*TT*TT), nrow=p*TT)
  for (j1 in 1:p) {
    for (j2 in 1:p) {
      if (B[j1,j2] != 0) {
        for (t1 in 1:TT) {
          for (t2 in 1:TT) {
            Gamma[TT*(j1-1) + t1, TT*(j2-1) + t2] <- rnorm(n=1, mean=B[j1,j2], sd=0.01)
          }
        }
      }
    }
  }
  
  # make Gamma symmetric for easier interpretability
  Gamma <- (t(Gamma) + Gamma) / 2
  
  ################################## Simulate Data ############################
  
  # generate Xs
  X_obs = matrix(nrow=n_obs, ncol=p*TT)
  for (i in 1:n_obs) {
    X_obs[i,] <- rmvnorm(n=1, sigma=X_cov)
  }
  
  # generate ys using Xs
  y_obs <- rep(NA, n_obs)
  for (i in 1:n_obs) {
    y_obs[i] <- mu + (t(alpha) %*% X_obs[i,])[1] + 
      (t(X_obs[i,]) %*% Gamma %*% X_obs[i,])[1] +
      rnorm(n=1, mean=0, sd=sqrt(sigma2))
  }
  
  # here no need to compute induced coefficients since no factor model
  alpha_0 <- mu
  
  # return results
  output = list(y_obs, X_obs, alpha_0, alpha, Gamma)
  names(output) = c("y_obs", "X_obs", "alpha_0", "alpha", "Gamma")
  return(output)
}










