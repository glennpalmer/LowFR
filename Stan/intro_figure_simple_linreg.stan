// simple linear regression model for motivating intro figure

data {
  int<lower=0> N;
  int<lower=0> p;
  int<lower=0> TT;
  vector[N] y;
  matrix[N, p*TT] X;
}

parameters {
  vector[p*TT] beta;
  real<lower=0> sigma2;
}

model {
  beta ~ normal(0, sqrt(10));
  sigma2 ~ inv_gamma(1,1);
  y ~ normal(X * beta, sqrt(sigma2));
}


