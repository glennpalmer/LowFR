// *** This version of CorrQuadReg is the same as the one in CorrQuadReg.stan
// except that here we add an additional n x r matrix of covariates which are
// included in the model via linear regression with N(0,10) priors on the
// coefficients.
//
// A naive longitudinal quadratic regression model that learns a correlation
// parameter for coefficients within exposures and pairs of exposures across
// time. Referred to in the paper as CorrQuadReg.

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
  
  matrix blocks_to_symmetric_mat(matrix off_diag, matrix diag, int p, int TT) {
    matrix[p*TT,p*TT] Gamma;
    int curr_time_index = 1;
    int curr_exp_index = 1;
    for (i in 1:p) {
      for (j in 1:p) {
        curr_time_index = 1;
        if (i == j) {
          for (t1 in 1:TT) {
            for (t2 in 1:TT) {
              if (t2 > t1) {
                Gamma[(i-1)*TT + t1, (j-1)*TT + t2] = 0;
              }
              else {
                Gamma[(i-1)*TT + t1, (j-1)*TT + t2] = diag[i,curr_time_index];
                curr_time_index = curr_time_index + 1;
              }
            }
          }
        }
        else if (i >= j) {
          for (t1 in 1:TT) {
            for (t2 in 1:TT) {
              Gamma[(i-1)*TT + t1, (j-1)*TT + t2] = off_diag[curr_exp_index, curr_time_index];
              curr_time_index = curr_time_index + 1;
            }
          }
          curr_exp_index = curr_exp_index + 1;
        }
        else {
          for (t1 in 1:TT) {
            for (t2 in 1:TT) {
              Gamma[(i-1)*TT + t1, (j-1)*TT + t2] = 0;
            }
          }
        }
      }
    }
    Gamma = (Gamma + Gamma') / 2;
    return Gamma;
  }
}

data {
  int<lower=0> N;
  int<lower=0> p;
  int<lower=0> TT;
  int<lower=0> r;
  vector[N] y;
  matrix[N, p*TT] X;
  matrix[N, r] Z;
}

parameters {
  // regression coefficients
  real alpha_0;
  vector[p*TT] alpha;
  matrix[p, (TT*(TT+1))/2] Gamma_diag_blocks;
  matrix[(p*(p-1))/2, TT*TT] Gamma_off_diag_blocks;
  
  // coefficients for additional covariates
  vector[r] beta;
  
  // hyperparameters
  real<lower=0> phi2_alpha;
  real<lower=0> phi2_gamma;
  real<lower=0,upper=1> psi_alpha;
  real<lower=0,upper=1> psi_gamma;
  
  // variance of outcome
  real<lower=0> sigma2;
}

transformed parameters {
  // declare Psi_mat_alpha
  matrix[TT,TT] Psi_mat_alpha;
  
  // declare diag and off-diag Psi_mat_gammas
  matrix[TT*TT,TT*TT] Psi_mat_gamma_off_diag;
  matrix[(TT*(TT+1))/2, (TT*(TT+1))/2] Psi_mat_gamma_diag;
  
  // declare Gamma matrix
  matrix[p*TT,p*TT] Gamma;
  
  // define Psi_mat_alpha
  for (i in 1:TT) {
    for (j in 1:i) {
      if (i == j) {
        Psi_mat_alpha[i,j] = 1;
      }
      else {
        Psi_mat_alpha[i,j] = psi_alpha;
        Psi_mat_alpha[j,i] = psi_alpha;
      }
    }
  }
  
  // define Psi_mat_gamma matrices
  for (i in 1:(TT*TT)) {
    for (j in 1:i) {
      if (i == j) {
        Psi_mat_gamma_off_diag[i,j] = 1;
      }
      else {
        Psi_mat_gamma_off_diag[i,j] = psi_gamma;
        Psi_mat_gamma_off_diag[j,i] = psi_gamma;
      }
    }
  }
  for (i in 1:((TT*(TT+1))/2)) {
    for (j in 1:i) {
      if (i == j) {
        Psi_mat_gamma_diag[i,j] = 1;
      }
      else {
        Psi_mat_gamma_diag[i,j] = psi_gamma;
        Psi_mat_gamma_diag[j,i] = psi_gamma;
      }
    }
  }
  
  // define symmetric Gamma matrix using off-diagonal and diagonal blocks defined above
  Gamma = blocks_to_symmetric_mat(Gamma_off_diag_blocks, Gamma_diag_blocks, p, TT);
  
}


model {
  // define covariance for alpha
  matrix[p*TT,p*TT] alpha_cov = phi2_alpha * kronecker(diag_matrix(rep_vector(1,p)), Psi_mat_alpha);
  
  // priors
  alpha_0 ~ normal(0,sqrt(10));
  alpha ~ multi_normal(rep_vector(0,p*TT), alpha_cov);
  sigma2 ~ inv_gamma(1,1);
  phi2_alpha ~ inv_gamma(1,1);
  phi2_gamma ~ inv_gamma(1,1);
  psi_alpha ~ uniform(0,1);
  psi_gamma ~ uniform(0,1);
  
  // additional covariates
  beta ~ normal(0, sqrt(10));
  
  // loop over diag and off-diag blocks to assign priors to Gamma elements
  for (i in 1:p) {
    Gamma_diag_blocks[i,] ~ multi_normal(rep_vector(0, (TT*(TT+1))/2), phi2_gamma*Psi_mat_gamma_diag);
  }
  for (i in 1:((p*(p-1))/2)) {
    Gamma_off_diag_blocks[i,] ~ multi_normal(rep_vector(0,TT*TT), phi2_gamma*Psi_mat_gamma_off_diag);
  }
  
  // likelihood
  for (i in 1:N) {
    y[i] ~ normal(alpha_0 + X[i]*alpha + quad_form(Gamma, X[i,]') + Z[i]*beta, sigma2);
  }
}

generated quantities {
  vector[N] y_rep;
  for (i in 1:N) {
    y_rep[i] = normal_rng(alpha_0 + X[i]*alpha + quad_form(Gamma, X[i,]') + Z[i]*beta, sigma2);
  }
}


