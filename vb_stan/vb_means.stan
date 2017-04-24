data {
    int N;
    int n_mean;  // We will calculate this mean.  If n_mean == -1, calculate the beta mean.
}
parameters {
  // Variational parameters.  These must be identical to pump_model_vb.stan.
  real<lower=0> lambda_shape[N];
  real<lower=0> lambda_rate[N];
  real beta_mean;
  real<lower=0> beta_var;
}
transformed parameters {
  // Expected variational parameters.
  real e_lambda[N];
  real e_beta;
  
  for (n in 1:N) {
    e_lambda[n] = lambda_shape[n] / lambda_rate[n];
  }

  e_beta = exp(beta_mean + 0.5 * beta_var);
}
model {
    if (n_mean > 0) {
        target += e_lambda[n_mean];
    } else if (n_mean == -1) {
        target += e_beta;
    } else {
        reject("Invalid value for n_mean.");
    }
}
