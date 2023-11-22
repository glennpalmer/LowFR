// Final model for fitting cross-validation splits to ELEMENT data.

functions{
  matrix kronecker(matrix A, matrix B) {
    matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
    int m;
    int n;
    int p;
    int q;
    m = rows(A);
    n = cols(A);
    p = rows(B);
    q = cols(B);
    for (i in 1:m) {
      for (j in 1:n) {
        int row_start;
        int row_end;
        int col_start;
        int col_end;
        row_start = (i - 1) * p + 1;
        row_end = (i - 1) * p + p;
        col_start = (j - 1) * q + 1;
        col_end = (j - 1) * q + q;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    return C;
  }
}

data {
  int<lower=0> N_train;
  int<lower=0> N_test;
  int<lower=0> p;
  int<lower=0> k;
  int<lower=0> TT;
  vector[N_train] y_train;
  matrix[N_train, p*TT] X_train;
  vector[N_test] y_test;
  matrix[N_test, p*TT] X_test;
  vector[N_train] sex_train;
  vector[N_train] age_train;
  vector[N_train] bmi_train;
  vector[N_train] onset_train;
  vector[N_test] sex_test;
  vector[N_test] age_test;
  vector[N_test] bmi_test;
  vector[N_test] onset_test;
}

parameters {
  // parameters for regression of y on latent factors
  // intercept
  real mu;
  // main effects
  vector[k] beta1;
  vector[k] beta2;
  vector[k] beta3;
  vector[TT] omega1;
  vector[TT] omega2;
  vector[TT] omega3;
  // interactions
  matrix[k,k] B;
  matrix[TT,TT] W;
  // variance
  real<lower=0> sigma2;
  
  // parameters for main effects and interaction for sex and age
  real alpha_sex;
  real alpha_age;
  real alpha_bmi;
  real alpha_onset;
  
  // parameters for factor model
  matrix[p,k] Lambda;
  vector<lower=0>[p] Sigma;
  real<lower=0,upper=1> phi;
  matrix[N_train, k*TT] Eta;
  
  // parameters for multiplicative gamma process prior
  real<lower=0> delta[3];
  real<lower=0> xi1[k+TT];
  real<lower=0> xi2[k+TT];
  real<lower=0> xi3[k+TT];
  real<lower=0> a1;
  real<lower=0> a2;
  
  // parameters for multiplicative gamma process prior -- interactions
  real<lower=0> tau_int;
  real<lower=0> xi_int[k*k+TT*TT];
  real<lower=0> a1_int;
}

transformed parameters {
  // declare theta
  vector[k*TT] theta;
  
  // declare Omega
  matrix[k*TT, k*TT] Omega;
  
  // declare Phi_mat
  matrix[TT, TT] Phi_mat;
  
  // declare tau (for MGP)
  real<lower=0> tau[3];
  
  // define theta in terms of beta and omega
  for (j in 1:k) {
    for (t in 1:TT) {
      theta[t + (j-1) * TT] = beta1[j] * omega1[t] +
      beta2[j] * omega2[t] +
      beta3[j] * omega3[t];
    }
  }
  
  // define Phi_mat
  Phi_mat = rep_matrix(0,TT,TT);
  for (i in 1:TT) {
    for (j in 1:TT) {
      if (i == j) {
        Phi_mat[i, j] = 1;
      }
      else {
        Phi_mat[i, j] = phi;
      }
    }
  }
  
  // define Omega in terms of the beta_int's and omega_int's
  Omega = kronecker(B, W);
  Omega = (Omega + Omega') / 2;
  
  // define tau in terms of delta
  tau[1] = delta[1];
  tau[2] = delta[1] * delta[2];
  tau[3] = tau[2] * delta[3];
}

model {
  // define I_kron_phi for use in the model for eta_i
  matrix[k*TT, k*TT] I_kron_phi = kronecker(diag_matrix(rep_vector(1,k)), Phi_mat);
  
  // define Sigma_kron_phi for use in the model for x_i
  matrix[p*TT, p*TT] Sigma_kron_phi = kronecker(diag_matrix(Sigma), Phi_mat);
  
  // declare Lambda_kron_I in terms of Lambda
  matrix[TT*p, TT*k] Lambda_kron_I = kronecker(Lambda, diag_matrix(rep_vector(1,TT)));
  
  // priors for variance terms
  Sigma ~ inv_gamma(1,1);
  sigma2 ~ inv_gamma(1,1);
  
  // intercept
  mu ~ normal(0, sqrt(10));
  
  // beta and omega multiplicative gamma priors
  for (j in 1:k) {
    beta1[j] ~ normal(0, 1 / sqrt(xi1[j] * tau[1]));
    beta2[j] ~ normal(0, 1 / sqrt(xi2[j] * tau[2]));
    beta3[j] ~ normal(0, 1 / sqrt(xi3[j] * tau[3]));
  }
  
  for (t in 1:TT) {
    omega1[t] ~ normal(0, 1 / sqrt(xi1[k + t] * tau[1]));
    omega2[t] ~ normal(0, 1 / sqrt(xi2[k + t] * tau[2]));
    omega3[t] ~ normal(0, 1 / sqrt(xi3[k + t] * tau[3]));
  }
  
  // priors for sex and age covariates
  alpha_age ~ normal(0, sqrt(10));
  alpha_sex ~ normal(0, sqrt(10));
  alpha_bmi ~ normal(0, sqrt(10));
  alpha_onset ~ normal(0, sqrt(10));
  
  // priors for multiplicative gamma hyperparameters
  delta[1] ~ gamma(a1, 1);
  delta[2] ~ gamma(a2, 1);
  delta[3] ~ gamma(a2, 1);
  xi1 ~ gamma(1.5, 1.5);
  xi2 ~ gamma(1.5, 1.5);
  xi3 ~ gamma(1.5, 1.5);
  a1 ~ gamma(2,1);
  a2 ~ gamma(2,1);
  
  // priors for multiplicative gamma hyperparameters -- interactions
  tau_int ~ gamma(a1_int, 1);
  xi_int ~ gamma(1.5, 1.5);
  a1_int ~ gamma(2,1);
  
  // quadratic regression terms
  for (i in 1:k) {
    for (j in 1:k) {
      B[i,j] ~ normal(0, 1 / sqrt(xi_int[(i-1)*k + j] * tau_int));
    }
  }
  for (i in 1:TT) {
    for (j in 1:TT) {
      W[i,j] ~ normal(0, 1 / sqrt(xi_int[k * k + (i-1)*TT + j] * tau_int));
    }
  }
  
  // prior distributions for factor model parameters
  for (i in 1:N_train) {
    Eta[i,] ~ multi_normal(rep_vector(0,k*TT), I_kron_phi);
  }
  for (i in 1:p) {
    for (j in 1:k) {
      Lambda[i,j] ~ normal(0, sqrt(10));
    }
  }
  phi ~ uniform(0,1);
  
  // model for X and y
  for (i in 1:N_train) {
    X_train[i,] ~ multi_normal(Lambda_kron_I * to_vector(Eta[i,]), Sigma_kron_phi);
    
    y_train[i] ~ normal(mu + Eta[i] * theta + quad_form(Omega, Eta[i,]') +
                  alpha_sex * sex_train[i] + alpha_age * age_train[i] +
                  alpha_bmi * bmi_train[i] + alpha_onset * onset_train[i], sqrt(sigma2));
  }
}


generated quantities {
  // generate the induced intercept, main effects, and interactions
  matrix[k*TT, p*TT] A;
  matrix[k*TT, k*TT] V;
  matrix[p,p] Sigma_inv;
  real alpha_0; // induced intercept
  vector[p*TT] alpha; // induced linear coefficients of y on x
  matrix[p*TT, p*TT] Gamma; // induced symmetric interaction matrix for regression of y on x
  vector[N_test] y_test_pred;
  vector[N_test] pred_errors;
  
  // Define Sigma_inv
  Sigma_inv = rep_matrix(0,p,p);
  for (j in 1:p) {
    Sigma_inv[j,j] = 1 / Sigma[j];
  }
  
  // Define A and V
  V = kronecker(inverse_spd(Lambda' * Sigma_inv * Lambda + diag_matrix(rep_vector(1,k))), Phi_mat);
  A = V * kronecker(Lambda' * Sigma_inv, inverse_spd(Phi_mat));
  
  // Define induced regression terms
  alpha_0 = mu + trace(Omega * V);
  alpha = A' * theta;
  Gamma = A' * Omega * A;
  
  // compute predicted y_test values
  for (i in 1:N_test) {
    y_test_pred[i] = alpha_0 + X_test[i]*alpha + quad_form(Gamma, X_test[i,]') +
                  alpha_sex * sex_test[i] + alpha_age * age_test[i] +
                  alpha_bmi * bmi_test[i] + alpha_onset * onset_test[i];
    
  }
  pred_errors = y_test_pred - y_test;
}

