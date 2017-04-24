data {
  int<lower=0> N;  // total number of observations
  real<lower=0> t[N];       // Time to failure
}
parameters {
// Prior parameters unconstrained so that grads are in the original space.
  real beta_prior_shape;
  real beta_prior_scale;
  real omega_0;

  real<lower=0> lambda[N];
  real<lower=0> beta;
}
model {
  // The log prior only.
  beta ~ inv_gamma(beta_prior_shape, 1 / beta_prior_scale);

  for (n in 1:N) {
    lambda[n] ~ gamma(omega_0, 1 / beta);
  }
}
