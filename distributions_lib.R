# Requires data_dat in the global environment
BetaLogPrior <- function(beta) {
  # beta ~ inv_gamma(beta_prior_shape, 1 / beta_prior_scale);
  dinvgamma(beta, shape=data_dat$beta_prior_shape, scale=data_dat$beta_prior_scale, log=TRUE)
}

# Requires data_dat in the global environment
LambdaStarLogPrior <- function(lambda_star) {
  # lambda_star[n] | beta ~ gamma(omega_0, 1);
  return((data_dat$omega_0 - 1) * log(lambda_star) - lambda_star - lgamma(data_dat$omega_0))
}



########################################
# MCMC functions


# Use the LRVBUtils version instead
# GetMCMCInfluenceFunctions <- function(param_draws, GetLogPrior) {
#   GetEYGivenX <- function(x_draws, y_draws) {
#     e_y_given_x <- loess.smooth(x_draws, y_draws)
#     return(approx(e_y_given_x$x, e_y_given_x$y, xout=x_draws)$y)
#   }
#   
#   mcmc_dens <- density(param_draws)
#   dens_at_draws <- approx(mcmc_dens$x, mcmc_dens$y, xout=param_draws)
#   mcmc_influence_ratio <- exp(log(dens_at_draws$y) - GetLogPrior(param_draws))
#   mcmc_importance_ratio <- exp(GetLogPrior(param_draws) - log(dens_at_draws$y))
#   
#   GetMCMCInfluence <- function(g_draws) {
#     conditional_mean_diff <- GetEYGivenX(x_draws=param_draws, y_draws=g_draws) - mean(g_draws)
#     return(conditional_mean_diff * mcmc_influence_ratio)
#   }
#   
#   GetMCMCWorstCase <- function(g_draws, correct=TRUE) {
#     mcmc_influence <- GetMCMCInfluence(g_draws)
#     if (correct) {
#       mcmc_influence_pos <- mcmc_influence * (mcmc_influence > 0)
#       mcmc_influence_neg <- -1 * mcmc_influence * (mcmc_influence < 0)
#       mcmc_wc <- max(
#         sqrt(mean((mcmc_influence_pos^2) * mcmc_importance_ratio)),
#         sqrt(mean((mcmc_influence_neg^2) * mcmc_importance_ratio)))
#     } else {
#       # This incorrect method appears to be what the prior marginals paper used.
#       mcmc_wc <- sqrt(mean((mcmc_influence^2) * mcmc_importance_ratio))
#     }
#     return(mcmc_wc)  
#   }
#   
#   return(list(GetMCMCInfluence=GetMCMCInfluence, GetMCMCWorstCase=GetMCMCWorstCase, dens_at_draws=dens_at_draws))
# }


GetMCMCWorstCaseResults <- function(MCMCWorstCaseFunction, prior_metric) {
  mcmc_beta_worst_case_list <- list()
  cat("Getting sensitivities to ", prior_metric, " \n")
  cat("   ... of beta\n")
  mcmc_beta_worst_case_list[[length(mcmc_beta_worst_case_list) + 1]] <-
    GetParamRow(value=MCMCWorstCaseFunction(mcmc_draws$beta),
                par="beta", metric=prior_metric, method="mcmc")

  for (ind in 1:data_dat$N) {
    cat("   ... of lambda ", ind, "\n")
    mcmc_beta_worst_case_list[[length(mcmc_beta_worst_case_list) + 1]] <-
      GetParamRow(value=MCMCWorstCaseFunction(mcmc_draws$lambda[, ind]), par="lambda", group=ind,
                  metric=prior_metric, method="mcmc")
  }
  return(do.call(rbind, mcmc_beta_worst_case_list))  
}


################################################
# VB functions

# Use the LRVBUtils one instead
# GetVariationalInfluenceResults <- function(
#   num_draws,
#   metric,
#   DrawImportanceSamples,
#   GetImportanceLogProb,
#   GetLogQGradTerm,
#   GetLogQ,
#   GetLogPrior,
#   moment_df=moment_df) {
#   
#   u_draws <- DrawImportanceSamples(num_draws)
#   
#   GetImportanceLogRatio <- function(u) { GetLogPrior(u) - GetImportanceLogProb(u) }
#   GetInfluenceLogRatio <- function(u) { GetLogQ(u) - GetLogPrior(u) }
#   importance_lp_ratio <- sapply(u_draws, GetImportanceLogRatio)
#   influence_lp_ratio <- sapply(u_draws, GetInfluenceLogRatio)
#   influence_fun  <-
#     do.call(rbind, lapply(u_draws, function(u) { GetLogQGradTerm(u) })) * exp(influence_lp_ratio)
#   u_influence_mat <- (influence_fun ^ 2) * exp(importance_lp_ratio)
#   u_influence_mat_pos <- ((influence_fun > 0) * influence_fun ^ 2) * exp(importance_lp_ratio)
#   u_influence_mat_neg <- ((influence_fun < 0) * influence_fun ^ 2) * exp(importance_lp_ratio)
#   
#   worst_case <-
#     sapply(1:ncol(influence_fun),
#            function(ind) { sqrt(max(mean(u_influence_mat_pos[, ind]),
#                                     mean(u_influence_mat_neg[, ind]))) })
#   
#   return(list(
#     df=GetParamRow(value=worst_case, par=moment_df$par, metric=metric, method="lrvb", group=moment_df$group),
#     u_draws=u_draws,
#     influence_fun=influence_fun,
#     importance_lp_ratio=importance_lp_ratio,
#     influence_lp_ratio=influence_lp_ratio))
# }
# 


BetaLogQ <- function(beta) {
  dlnorm(beta, meanlog=vb_fit$par["beta_mean"], sdlog=sqrt(vb_fit$par["beta_var"]), log=TRUE)
}

GetGetBetaLogQGradTermFunction <- function() {
  # Don't actually fit -- just get a fit object.
  vb_log_q_beta_fitobj <- sampling(vb_log_density, list(N=data_dat$N, n=-1), iter=1, chains=1)
  beta_draw_dat <- vb_fit_dat
  beta_draw_dat$lambda_star <- rep(0, data_dat$N) # not used
  vb_pars_free <- unconstrain_pars(vb_means_fitobj, vb_fit_dat)
  param_ids <- 1:length(vb_pars_free)
  
  GetSingleBetaLogQGradTerm <- function(beta) {
    beta_draw_dat$beta <- beta
    vb_pars_free <- unconstrain_pars(vb_log_q_beta_fitobj, beta_draw_dat)
    grad_log_q <- grad_log_prob(vb_log_q_beta_fitobj, vb_pars_free, adjust_transform=FALSE)[param_ids]
    return(grad_log_q)
  }
  
  GetBetaLogQGradTerm <- function(beta_vec) {
    grad_log_q_mat <- sapply(beta_vec, GetSingleBetaLogQGradTerm)
    return(as.matrix(-1 * moment_jac %*% solve(vb_fit$hessian, grad_log_q_mat)))
  }
  
  return(GetBetaLogQGradTerm)
}

beta_importance_scale <- 1
GetBetaImportanceLogProb <- function(u) {
  dlnorm(u, meanlog=vb_fit$par["beta_mean"], sdlog=beta_importance_scale * sqrt(vb_fit$par["beta_var"]), log=TRUE)
}

GetBetaImportanceDraws <- function(num_draws) {
  rlnorm(num_draws, meanlog=vb_fit$par["beta_mean"], sdlog=beta_importance_scale * sqrt(vb_fit$par["beta_var"]))
}



# Lambda star

LambdaStarLogQ <- function(lambda, ind) {
  lambda_star_shape <- GetVBLambdaStarShape(ind, vb_fit)
  lambda_star_rate <- GetVBLambdaStarRate(ind, vb_fit)
  return(dgamma(lambda, shape=lambda_star_shape, rate=lambda_star_rate, log=TRUE))
}

GetLambdaStarLogQGradTermFunction <- function(ind) {
  vb_log_q_lambda_fitobj <- sampling(vb_log_density, list(N=data_dat$N, n=ind), iter=1, chains=1)
  lambda_draw_dat <- vb_fit_dat
  lambda_draw_dat$beta <- 0 # not used
  lambda_draw_dat$lambda_star <- rep(0, data_dat$N) # Only the [ind] element is used.
  vb_pars_free <- unconstrain_pars(vb_means_fitobj, vb_fit_dat)
  param_ids <- 1:length(vb_pars_free)
  
  # TODO: subtract the importance-weighted means?
  GetSingleLambdaStarLogQGradTerm <- function(lambda_star) {
    lambda_draw_dat$lambda_star[ind] <- lambda_star
    vb_pars_free <- unconstrain_pars(vb_log_q_lambda_fitobj, lambda_draw_dat)
    grad_log_q <- grad_log_prob(vb_log_q_lambda_fitobj, vb_pars_free, adjust_transform=FALSE)[param_ids]
    return(grad_log_q)
  }
  
  GetLambdaStarLogQGradTerm <- function(lambda_star_vec) {
    grad_log_q_mat <- sapply(lambda_star_vec, GetSingleLambdaStarLogQGradTerm)
    return(as.matrix(-1 * moment_jac %*% solve(vb_fit$hessian, grad_log_q_mat)))
  }
  
  return(GetLambdaStarLogQGradTerm)
}

lambda_importance_scale <- 2
GetLambdaStarImportanceLogProb <- function(u) {
  dgamma(u, shape=GetVBLambdaStarShape(ind, vb_fit) / (lambda_importance_scale^2),
         rate=GetVBLambdaStarRate(ind, vb_fit) / (lambda_importance_scale^2), log=TRUE)
}
GetLambdaStarImportanceDraws <- function(num_draws) {
  rgamma(num_draws, shape=GetVBLambdaStarShape(ind, vb_fit) / (lambda_importance_scale^2),
         rate=GetVBLambdaStarRate(ind, vb_fit) / (lambda_importance_scale^2))
}
GetLambdaStarImportanceLogProbRatio <- function(u) {
  return(LambdaStarLogPrior(u) - GetLambdaStarImportanceDensity(u))
}





