// The log prior for the Hessian.
data {
  int<lower=0> N;  // total number of observations
}
parameters {
  // Variational parameters.
  real<lower=0> lambda_shape[N];
  real<lower=0> lambda_rate[N];
  real beta_mean;
  real<lower=0> beta_var;

  // Prior parameters unconstrained so that grads are in the original space.
  real beta_prior_shape;
  real beta_prior_scale;
  real omega_0;
}
transformed parameters {
  // Expected variational parameters.
  real e_lambda[N];
  real e_log_lambda[N];
  real e_beta_inv;
  real e_log_beta;

  for (n in 1:N) {
    e_lambda[n] = lambda_shape[n] / lambda_rate[n];
    e_log_lambda[n] = digamma(lambda_shape[n]) - log(lambda_rate[n]);
  }

  e_beta_inv = exp(0.5 * beta_var - beta_mean);
  e_log_beta = beta_mean;
}
model {
  // priors
  // beta ~ inv_gamma(beta_prior_shape, 1 / beta_prior_scale);
  target += -1 * (beta_prior_shape + 1) * e_log_beta - e_beta_inv / beta_prior_scale;

  // Data
  for (n in 1:N) {
    // lambda[n] ~ gamma(omega_0, 1 / beta);
    target += (omega_0 - 1) * e_log_lambda[n] - e_lambda[n] * e_beta_inv - omega_0 * e_log_beta;
  }
}
