data {
    int N;

    // We will calculate this lambda mean,
    // unless n_mean == -1, in which case calculate the beta mean.
    int n_mean;  
}
parameters {
  // Variational parameters.  These must be identical to pump_model_vb.stan.
  real<lower=0> lambda_star_shape[N];
  real<lower=0> lambda_star_rate[N];
  real beta_mean;
  real<lower=0> beta_var;
}
transformed parameters {
  // Expected variational parameters.
  real e_lambda[N];
  real e_beta;
  
  e_beta = exp(beta_mean + 0.5 * beta_var);
  for (n in 1:N) {
    e_lambda[n] = e_beta * lambda_star_shape[n] / lambda_star_rate[n];
  }
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
