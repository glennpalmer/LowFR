// Finalized model from "Low-rank longitudinal factor regression" paper
// (This implementation does not include additional covariates and cannot
// handle missing data.)

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
  int<lower=0> N;
  int<lower=0> p;
  int<lower=0> k;
  int<lower=0> TT;
  int<lower=0> H;
  vector[N] y;
  matrix[N, p*TT] X;
}

parameters {
  // parameters for regressoion of y on latent factors
  // intercept
  real mu;
  // main effects
  matrix[k,H] beta;
  matrix[TT,H] omega;
  // interactions
  matrix[k,k] B;
  matrix[TT,TT] W;
  // variance
  real<lower=0> sigma2;
  
  // parameters for factor model
  matrix[p,k] Lambda;
  vector<lower=0>[p] Sigma;
  real<lower=0, upper=1> phi;
  matrix[N, k*TT] Eta;
  
  // parameters for multiplicative gamma process prior
  real<lower=0> delta[H];
  real<lower=0> xi[k+TT,H];
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
  real<lower=0> tau[H];
  
  // define theta in terms of beta and omega
  theta = to_vector(omega * beta');
  
  // define Omega in terms of B and W
  Omega = kronecker(B, W);
  Omega = (Omega + Omega') / 2;
  
  // define Phi_mat in terms of phi
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

  // construct tau in terms of delta
  tau[1] = delta[1];
  for (l in 2:H) {
    tau[l] = tau[l-1] * delta[l];
  }
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
    for (l in 1:H) {
      beta[j,l] ~ normal(0, 1 / sqrt(xi[j,l] * tau[l]));
    }
  }
  
  for (t in 1:TT) {
    for (l in 1:H) {
      omega[t,l] ~ normal(0, 1 / sqrt(xi[k+t,l] * tau[l]));
    }
  }
  
  // priors for multiplicative gamma hyperparameters
  delta[1] ~ gamma(a1, 1);
  for (l in 2:H) {
    delta[l] ~ gamma(a2,1);
  }
  for (j in 1:(k+TT)) {
    for (l in 1:H) {
      xi[j,l] ~ gamma(1.5,1.5);
    }
  }
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
  for (i in 1:N) {
    Eta[i,] ~ multi_normal(rep_vector(0,k*TT), I_kron_phi);
  }
  for (i in 1:p) {
    for (j in 1:k) {
      Lambda[i,j] ~ normal(0, sqrt(10));
    }
  }
  phi ~ uniform(0,1);
  
  // model for X and y
  for (i in 1:N) {
    X[i,] ~ multi_normal(Lambda_kron_I * to_vector(Eta[i,]), Sigma_kron_phi);
    
    y[i] ~ normal(mu + Eta[i] * theta + quad_form(Omega, Eta[i,]'), sqrt(sigma2));
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
}


