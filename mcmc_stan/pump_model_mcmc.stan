data {
  int<lower=0> N;  // total number of observations
  int<lower=0> s[N];       // Failure count
  real<lower=0> t[N];       // Time to failure

  // Prior parameters
  real<lower=0> beta_prior_shape;
  real<lower=0> beta_prior_scale;
  real<lower=0> omega_0;
}
parameters {
  real<lower=0> lambda[N];
  real<lower=0> beta;
}
model {
  //data variance priors
  beta ~ inv_gamma(beta_prior_shape, 1 / beta_prior_scale);

  for (n in 1:N) {
    lambda[n] ~ gamma(omega_0, 1 / beta);
    s[n] ~ poisson(t[n] * lambda[n]);
  }
}
generated quantities {
    real log_prior;
    real beta_log_prior;
    real lambda_log_prior[N];

    log_prior = 0;
    beta_log_prior = inv_gamma_lpdf(beta | beta_prior_shape, 1 / beta_prior_scale);
    log_prior = log_prior + beta_log_prior; 
    for (n in 1:N) {
      lambda_log_prior[n] = gamma_lpdf(lambda[n] | omega_0, 1 / beta);
      log_prior = log_prior + lambda_log_prior[n]; 
    }
}
