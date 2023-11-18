// rank-1 model for motivating intro figure

data {
  int<lower=0> N;
  int<lower=0> p;
  int<lower=0> TT;
  vector[N] y;
  matrix[N, p*TT] X;
}

parameters {
  vector[p] beta1;
  vector[TT] omega1;
  real<lower=0> sigma2;
}

transformed parameters {
  // declare theta
  vector[p*TT] theta;
  
  // define theta in terms of beta and omega
  for (j in 1:p) {
    for (t in 1:TT) {
      theta[t + (j-1) * TT] = beta1[j] * omega1[t];
    }
  }
}

model {
  beta1 ~ normal(0, sqrt(10));
  omega1 ~ normal(0, sqrt(10));
  sigma2 ~ inv_gamma(1,1);
  y ~ normal(X * theta, sqrt(sigma2));
}

