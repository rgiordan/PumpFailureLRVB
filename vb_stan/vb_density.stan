data {
  int<lower=0> N;  // total number of observations

  // If n == -1, then get the log density of beta.  Otherwise get the
  // density for the corresponding element of lambda.
  int n;
}
parameters {
  // Variational parameters.
  real<lower=0> lambda_shape[N];
  real<lower=0> lambda_rate[N];
  real beta_mean;
  real<lower=0> beta_var;
  
  // The values at which we evaluate the variational log density
  real lambda[N];
  real beta;
  }
model {
  if (n > 0) {
      target += gamma_lpdf(lambda[n] | lambda_shape, lambda_rate);
  } else if (n == -1) {
      target += lognormal_lpdf(beta | beta_mean, sqrt(beta_var));
  } else {
      reject("Invalid value for n.");
  }
}
