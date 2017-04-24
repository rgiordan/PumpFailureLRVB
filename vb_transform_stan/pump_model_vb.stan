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
  // Variational parameters.
  real<lower=0> lambda_star_shape[N];
  real<lower=0> lambda_star_rate[N];
  real beta_mean;
  real<lower=0> beta_var;
}
transformed parameters {
  // Expected variational parameters.
  real e_lambda_star[N];
  real e_log_lambda_star[N];
  real e_beta;
  real e_beta_inv;
  real e_log_beta;
  
  for (n in 1:N) {
    e_lambda_star[n] = lambda_star_shape[n] / lambda_star_rate[n];
    e_log_lambda_star[n] = digamma(lambda_star_shape[n]) - log(lambda_star_rate[n]);
  }

  e_beta_inv = exp(0.5 * beta_var - beta_mean);
  e_beta = exp(0.5 * beta_var + beta_mean);
  e_log_beta = beta_mean;
}
model {
  // priors
  // beta ~ inv_gamma(beta_prior_shape, 1 / beta_prior_scale);
  target += -1 * (beta_prior_shape + 1) * e_log_beta - e_beta_inv / beta_prior_scale;

  // Data
  for (n in 1:N) {
    // s[n] ~ poisson(t[n] * (lambda_star[n] * beta));
    target += -1 * e_beta * e_lambda_star[n] * t[n] +
              s[n] * (e_log_lambda_star[n] + e_log_beta);
    
    // lambda[n] ~ gamma(omega_0, 1);
    target += (omega_0 - 1) * e_log_lambda_star[n] - e_lambda_star[n];
  }
  
  // Entropy
  // Beta
  target += 0.5 * log(beta_var) + exp(0.5 + beta_mean);

  // Lambda
  for (n in 1:N) {
    target += lambda_star_shape[n] - log(lambda_star_rate[n]) + lgamma(lambda_star_shape[n]) +
              (1 - lambda_star_shape[n]) * digamma(lambda_star_shape[n]);
  }
}
