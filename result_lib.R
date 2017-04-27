LoadInEnvironment <- function(rdata_file) {
  new_env <- new.env()
  load(rdata_file, envir=new_env)
  return(new_env)
}


GetParamRow <- function(value, par, metric, method, group=-1) {
  data.frame(value=value, par=par, metric=metric, method=method, group=group)
}

GetMCMCLambdaGroupResults <- function(ind, mcmc_draws) {
  return(rbind(
    GetParamRow(value=mean(mcmc_draws$lambda[, ind]),
                par="lambda", group=ind, method="mcmc", metric="mean"),
    GetParamRow(value=sd(mcmc_draws$lambda[, ind]),
                par="lambda", group=ind, method="mcmc", metric="sd")
  ))
}

GetMCMCBetaResults <- function(mcmc_draws) {
  return(rbind(
    GetParamRow(value=mean(mcmc_draws$beta), par="beta", method="mcmc", metric="mean"),
    GetParamRow(value=sd(mcmc_draws$beta), par="beta", method="mcmc", metric="sd")
  ))
}

GetMCMCLambdaSensitivityResults <- function(mcmc_draws, log_prior_grad_mat, prior_param_names) {
  lambda_draws <- mcmc_draws$lambda - colMeans(mcmc_draws$lambda)
  sens <- cov(lambda_draws, log_prior_grad_mat)

  # normalize the sensitivity by the posterior standard deviation.
  lambda_sd_scale <- sqrt(diag(cov(lambda_draws)))
  result_list <- list()
  for (ind in 1:ncol(lambda_draws)) {
    for (prior_ind in 1:length(prior_param_names)) {
      result_list[[length(result_list) + 1]] <-
        GetParamRow(value=sens[ind, prior_ind] / lambda_sd_scale[ind],
                    par="lambda", group=ind, method="mcmc",
                    metric=prior_param_names[prior_ind])
    }
  }
  return(do.call(rbind, result_list))
}

GetMCMCBetaSensitivityResults <- function(mcmc_draws, log_prior_grad_mat, prior_param_names) {
  beta_draws <- mcmc_draws$beta - mean(mcmc_draws$beta)
  sens <- cov(beta_draws, log_prior_grad_mat)

  # normalize the sensitivity by the posterior standard deviation.
  beta_sd_scale <- sd(beta_draws)
  result_list <- list()
  for (prior_ind in 1:length(prior_param_names)) {
    result_list[[length(result_list) + 1]] <-
      GetParamRow(value=sens[1, prior_ind] / beta_sd_scale,
                  par="beta", method="mcmc", metric=prior_param_names[prior_ind])
  }
  return(do.call(rbind, result_list))
}

GetMCMCBetaResults <- function(mcmc_draws) {
  return(rbind(
    GetParamRow(value=mean(mcmc_draws$beta), par="beta", method="mcmc", metric="mean"),
    GetParamRow(value=sd(mcmc_draws$beta), par="beta", method="mcmc", metric="sd")
  ))
}

GetVBLambdaShape <- function(ind, vb_fit) {
    vb_fit$par[[paste("lambda_shape[", ind, "]", sep="")]]
}

GetVBLambdaRate <- function(ind, vb_fit) {
    vb_fit$par[[paste("lambda_rate[", ind, "]", sep="")]]
}

GetVBLambdaGroupResults <- function(ind, vb_fit) {
  lambda_shape <- GetVBLambdaShape(ind, vb_fit)
  lambda_rate <- GetVBLambdaRate(ind, vb_fit)
  return(rbind(
    GetParamRow(value=lambda_shape / lambda_rate,
                par="lambda", group=ind, method="mfvb", metric="mean"),
    GetParamRow(value=sqrt(lambda_shape) / lambda_rate,
                par="lambda", group=ind, method="mfvb", metric="sd")
  ))
}

GetVBBetaResults <- function(vb_fit) {
  # They are actually parameters of a lognormal, not the mean and variance.
  beta_mu <- vb_fit$par[["beta_mean"]]
  beta_sigma2 <- vb_fit$par[["beta_var"]]
  beta_mean <- exp(beta_mu + 0.5 * beta_sigma2)
  beta_sd <- beta_mean * sqrt(exp(beta_sigma2) - 1)
  return(rbind(
    GetParamRow(value=beta_mean, par="beta", method="mfvb", metric="mean"),
    GetParamRow(value=beta_sd, par="beta", method="mfvb", metric="sd")
  ))
}


GetVBLambdaStarShape <- function(ind, vb_fit) {
    vb_fit$par[[paste("lambda_star_shape[", ind, "]", sep="")]]
}

GetVBLambdaStarRate <- function(ind, vb_fit) {
    vb_fit$par[[paste("lambda_star_rate[", ind, "]", sep="")]]
}

GetVBLambdaStarGroupResults <- function(ind, vb_fit) {
  lambda_star_shape <- GetVBLambdaStarShape(ind, vb_fit)
  lambda_star_rate <- GetVBLambdaStarRate(ind, vb_fit)
  beta_mu <- vb_fit$par[["beta_mean"]]
  beta_sigma2 <- vb_fit$par[["beta_var"]]

  e_beta <- exp(beta_mu + 0.5 * beta_sigma2)
  var_beta <- e_beta^2 * (exp(beta_sigma2) - 1)

  e_lambda_star <- lambda_star_shape / lambda_star_rate
  var_lambda_star <- lambda_star_shape / (lambda_star_rate ^ 2)

  e_lambda <- e_beta * e_lambda_star
  var_lambda <- var_beta * var_lambda_star + var_beta * (e_lambda_star^2) + (e_beta^2) * var_lambda_star

  return(rbind(
    GetParamRow(value=e_lambda, par="lambda", group=ind, method="mfvb", metric="mean"),
    GetParamRow(value=sqrt(var_lambda), par="lambda", group=ind, method="mfvb", metric="sd")
  ))
}
