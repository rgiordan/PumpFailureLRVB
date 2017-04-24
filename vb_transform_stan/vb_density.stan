data {
  int<lower=0> N;  // total number of observations

  // If n == -1, then get the log density of beta.  Otherwise get the
  // density for the corresponding element of lambda.
  int n;
}
parameters {
  // Variational parameters.
  real<lower=0> lambda_star_shape[N];
  real<lower=0> lambda_star_rate[N];
  real beta_mean;
  real<lower=0> beta_var;
  
  // The values at which we evaluate the variational log density
  real lambda_star[N];
  real beta;
  }
model {
  if (n > 0) {
      target += gamma_lpdf(lambda_star[n] | lambda_star_shape[n], lambda_star_rate[n]);
  } else if (n == -1) {
      target += lognormal_lpdf(beta | beta_mean, sqrt(beta_var));
  } else {
      reject("Invalid value for n.");
  }
}
