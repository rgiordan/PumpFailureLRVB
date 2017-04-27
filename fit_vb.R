library(rstan)
library(dplyr)
library(reshape2)
library(ggplot2)
library(numDeriv)
library(invgamma)
library(gridExtra)
library(LRVBUtils)

project_directory <- file.path(Sys.getenv("GIT_REPO_LOC"), "PumpFailureLRVB")
data_directory <- file.path(project_directory, "data")
source(file.path(project_directory, "result_lib.R"))
source(file.path(project_directory, "distributions_lib.R"))

model_name <- "pump_model_transform"

save_filename <- paste(model_name, "vb_analysis.Rdata", sep="_")

rdata_file <- file.path(data_directory, paste(model_name, "stan.Rdata", sep="_"))
if (file.exists(rdata_file)) {
  print("Loading from file.")
  load(rdata_file)
} else {
  print("Compiling model.")
  vb_means <- stan_model(file.path(project_directory, "vb_transform_stan/vb_means.stan"))
  vb_model <- stan_model(file.path(project_directory, "vb_transform_stan/pump_model_vb.stan"))
  vb_log_prior <- stan_model(file.path(project_directory, "vb_transform_stan/pump_model_log_prior_vb.stan"))
  vb_log_density <- stan_model(file.path(project_directory, "vb_transform_stan/vb_density.stan"))

  mcmc_model <- stan_model(file.path(project_directory, "mcmc_stan/pump_model_mcmc.stan"))
  mcmc_log_prior <- stan_model(file.path(project_directory, "mcmc_stan/pump_model_mcmc_log_prior.stan"))

  save(vb_model, vb_means, vb_log_prior, vb_log_density,
       mcmc_model, mcmc_log_prior, file=rdata_file)
}

mcmc_results_file <- paste(model_name, "mcmc.Rdata", sep="_")
mcmc_env <- LoadInEnvironment(file.path(project_directory, "data", mcmc_results_file))

pump_data <- read.csv(file.path(project_directory, "data/pump_failure_data.csv"))
data_dat <- list(N=nrow(pump_data), s=pump_data$failures, t=pump_data$thousand_hours,
                 beta_prior_shape=0.01, beta_prior_scale=1.0, omega_0=1.802)


#####################
# VB fit

lambda_star_mean <- colMeans(mcmc_env$mcmc_draws$lambda)
vb_dat <- list(lambda_star_shape=lambda_star_mean , lambda_star_rate=rep(1.0, data_dat$N),
               beta_mean=mean(mcmc_env$mcmc_draws$beta), beta_var=var(mcmc_env$mcmc_draws$beta))
vb_fit <- optimizing(vb_model, data=data_dat, init=vb_dat, algorithm="BFGS",
                     hessian=TRUE, verbose=TRUE, iter=1000, tol_rel_grad=1e-16)
vb_dat <- list(lambda_star_shape=vb_fit$par[sprintf("lambda_star_shape[%d]", 1:10)],
               lambda_star_rate=vb_fit$par[sprintf("lambda_star_rate[%d]", 1:10)],
               beta_mean=vb_fit$par["beta_mean"],
               beta_var=vb_fit$par["beta_var"])
vb_fit <- optimizing(vb_model, data=data_dat, init=vb_dat, algorithm="Newton", hessian=TRUE, verbose=TRUE, iter=100)

vb_results <- rbind(
  do.call(rbind, lapply(1:data_dat$N, function(ind) { GetVBLambdaStarGroupResults(ind, vb_fit) } )),
  GetVBBetaResults(vb_fit))


################################
# LRVB

# Get the Jacobian of the mean transform.  Check that the means match.

lambda_star_shapes <- sapply(1:data_dat$N, function(ind) { vb_fit$par[[paste("lambda_star_shape[", ind, "]", sep="")]] })
vb_fit_dat <- list(lambda_star_shape=sapply(1:data_dat$N, function(ind) GetVBLambdaStarShape(ind, vb_fit)),
                   lambda_star_rate=sapply(1:data_dat$N, function(ind) GetVBLambdaStarRate(ind, vb_fit)),
                   beta_mean=vb_fit$par[["beta_mean"]], beta_var=vb_fit$par[["beta_var"]])

# Don't actually fit -- just get a fit object.  You can ignore sampling errors.
vb_means_fitobj <- sampling(vb_means, list(N=data_dat$N, n_mean=1), iter=1, chains=1)
vb_pars_free <- unconstrain_pars(vb_means_fitobj, vb_fit_dat)

# Get beta, and initialize the moment objects.
vb_means_fitobj <- sampling(vb_means, list(N=data_dat$N, n_mean=-1), iter=1, chains=1)
beta_mean <- log_prob(vb_means_fitobj, vb_pars_free, adjust_transform=FALSE, gradient=TRUE)
moment_jac <- t(attr(beta_mean, "gradient"))
moment_df <- GetParamRow(value=as.numeric(beta_mean), par="beta", method="stan_mfvb", metric="mean")

# Get the lambda means.
for (n in 1:data_dat$N) {
  print(n)
  vb_means_fitobj <- sampling(vb_means, list(N=data_dat$N, n_mean=n), iter=1, chains=1)
  lambda_mean <- log_prob(vb_means_fitobj, vb_pars_free, adjust_transform=FALSE, gradient=TRUE)
  moment_df <- rbind(moment_df,
                     GetParamRow(value=as.numeric(lambda_mean), par="lambda", group=n, method="stan_mfvb", metric="mean"))
  moment_jac <- rbind(moment_jac, t(attr(lambda_mean, "gradient")))
}

# These should match.
check_df <- rbind(filter(vb_results, metric=="mean"), moment_df) %>%
  dcast(par + group ~ method, value.var="value") %>%
  mutate(diff=mfvb - stan_mfvb)
stopifnot(all(abs(check_df$diff < 1e-8)))

# Sanity check.
vb_means_fitobj@.MISC$stan_fit_instance$param_names()
vb_means_fitobj@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)

# Calculate LRVB covariance.
# Dig that the hessian is actually now calculated numerically, not with autodiff.
lrvb_cov <- -1 * moment_jac %*% solve(vb_fit$hessian, t(moment_jac))
moment_lrvb_df <- moment_df %>% mutate(method="lrvb", metric="sd", value=sqrt(diag(lrvb_cov)))


#####################
# Look at the estimation results.
results <-
  rbind(vb_results, moment_lrvb_df, mcmc_env$mcmc_results) %>%
  dcast(par + metric + group ~ method, value.var="value")

if (FALSE) {
grid.arrange(
  ggplot(filter(results, metric=="mean")) +
    geom_point(aes(x=mcmc, y=mfvb, color=par)) +
    geom_abline(aes(intercept=0, slope=1))
,
  ggplot(filter(results, metric=="sd")) +
    geom_point(aes(x=mcmc, y=mfvb, color="mfvb", shape=par), size=3) +
    geom_point(aes(x=mcmc, y=lrvb, color="lrvb", shape=par), size=3) +
    geom_abline(aes(intercept=0, slope=1))
, ncol=2  
)
}




##############################################
# Sensitivity
##############################################

########################
# Get the VB parametric sensitivity.

# Normalize all the sensitivity measures by the LRVB standard deviation.
lrvb_sd_scale <- sqrt(diag(lrvb_cov))

vb_log_prior_dat <- vb_fit_dat
vb_log_prior_dat$beta_prior_shape <- data_dat$beta_prior_shape
vb_log_prior_dat$beta_prior_scale <- data_dat$beta_prior_scale
vb_log_prior_dat$omega_0 <- data_dat$omega_0

vb_log_prior_fitobj <- sampling(vb_log_prior, list(N=data_dat$N), iter=1, chains=1)
vb_log_prior_pars_free <- unconstrain_pars(vb_log_prior_fitobj, vb_log_prior_dat)

log_prior_grad <- function(vb_log_prior_pars_free) {
  grad_log_prob(vb_log_prior_fitobj, vb_log_prior_pars_free, adjust_transform=FALSE)
}

log_prior_hess <- jacobian(log_prior_grad, vb_log_prior_pars_free)
stopifnot(max(abs(log_prior_hess - t(log_prior_hess))) < 1e-8)
vb_pars_free <- unconstrain_pars(vb_means_fitobj, vb_fit_dat)
param_ids <- 1:length(vb_pars_free)
prior_ids <- setdiff(1:length(vb_log_prior_pars_free), param_ids)
log_prior_jac <- (0.5 * (log_prior_hess + t(log_prior_hess)))[param_ids, prior_ids]

prior_param_names <- vb_log_prior_fitobj@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)[prior_ids]
log_prior_sens <- -1 * moment_jac %*% solve(vb_fit$hessian, log_prior_jac) / lrvb_sd_scale

vb_sensitivity_list <- list()
for (prior_ind in 1:length(prior_param_names)) {
  vb_sensitivity_list[[length(vb_sensitivity_list) + 1]] <-
    moment_df %>% mutate(method="lrvb", metric=prior_param_names[prior_ind], value=log_prior_sens[, prior_ind])
}
vb_sensitivity_df <- do.call(rbind, vb_sensitivity_list)


###########################################
# Combine parametric sensitivity results

sens_results <-
  rbind(mcmc_env$mcmc_sensitivity_df, vb_sensitivity_df) %>%
  dcast(par + metric + group ~ method, value.var="value")


#####################################
# Non-parametric sensitivity

num_draws <- 5000

# Get the beta influence function results
# Note: you can ignore the Stan warning.
GetBetaLogQGradTerm <- GetGetBetaLogQGradTermFunction()

beta_influence_results <- GetVariationalInfluenceResults(
  num_draws,
  DrawImportanceSamples=GetBetaImportanceDraws,
  GetImportanceLogProb=GetBetaImportanceLogProb,
  GetLogQGradTerms=GetBetaLogQGradTerm,
  GetLogQ=BetaLogQ,
  GetLogPrior=BetaLogPrior)

vb_beta_worst_case_df <-
  with(beta_influence_results,
       GetParamRow(value=worst_case / lrvb_sd_scale, par=moment_df$par,
                   metric="beta", method="lrvb", group=moment_df$group))


# Get the lambda influence functions
lambda_worst_case_list <- list()
for (ind in 1:data_dat$N) {
  cat("Ind ", ind, "\n")
  lambda_influence_results <- GetVariationalInfluenceResults(
    num_draws,
    DrawImportanceSamples=GetLambdaStarImportanceDraws,
    GetImportanceLogProb=GetLambdaStarImportanceLogProb,
    GetLogQGradTerm=GetLambdaStarLogQGradTermFunction(ind),
    GetLogQ=function(u) { LambdaStarLogQ(u, ind) },
    GetLogPrior=LambdaStarLogPrior)

  lambda_influence_results$df <-
    with(lambda_influence_results,
         GetParamRow(value=worst_case / lrvb_sd_scale, par=moment_df$par,
                     metric=paste("lambda_star[", ind, "]", sep=""),
                     method="lrvb", group=moment_df$group))

  lambda_worst_case_list[[length(lambda_worst_case_list) + 1]] <- lambda_influence_results
}

vb_lambda_worst_case_df <- do.call(rbind, lapply(lambda_worst_case_list, function(x) { x$df }))

# Also normalize the gustafson result
gustafson_result_df <- mcmc_env$gustafson_result_df %>%
  mutate(value=value / sd) %>% select(-sd)

mcmc_worst_case_df <- mcmc_env$mcmc_worst_case_df %>%
  mutate(value=value / sd) %>% select(-sd)

worst_case_df <-
  rbind(gustafson_result_df,
        mcmc_worst_case_df,
        vb_beta_worst_case_df, vb_lambda_worst_case_df) %>%
  dcast(par + group + metric ~ method, value.var="value") %>%
  filter(metric != "lambda_star")


##########################
# Graph individual influence function

# Make a matrix of MCMC draws that matches the indices of moment_df
draws_mat <- matrix(NaN, length(mcmc_env$mcmc_draws$beta), nrow(moment_df))
draws_mat[, which(moment_df$par == "beta")] <- mcmc_env$mcmc_draws$beta
for (g in 1:max(moment_df$group)) {
  ind <- which(moment_df$par == "lambda" & moment_df$group==g)
  draws_mat[, ind] <- mcmc_env$mcmc_draws$lambda[, g] / mcmc_env$mcmc_draws$beta
}


GetMCMCInfluenceFunctionsDataFrame <- function(mcmc_inf_funs, ind, GetLogPrior, metric, num_mcmc_draws=5000) {
  param_draws <- mcmc_inf_funs$param_draws
  g_draws <- draws_mat[, ind]
  worst_case_results <- mcmc_inf_funs$GetMCMCWorstCaseResults(g_draws)
  mcmc_rows <- sample(length(param_draws), num_mcmc_draws)
  mcmc_inf_df <-
    data.frame(u=param_draws[mcmc_rows],
               log_prior=GetLogPrior(param_draws[mcmc_rows]),
               log_dens=log(mcmc_inf_funs$dens_at_draws[mcmc_rows]),
               worst_case_u=worst_case_results$worst_u[mcmc_rows],
               influence=worst_case_results$mcmc_influence[mcmc_rows],
               gbar=mcmc_inf_funs$GetConditionalMeanDiff(g_draws)[mcmc_rows],
               method="mcmc",
               metric=metric)
  return(mcmc_inf_df)
}


GetVBInfluenceFunctionsDataFrame <- function(vb_influence_results, ind, metric) {
  vb_inf_df <- with(vb_influence_results,
                    data.frame(u=u_draws, log_dens=log_q,
                               log_prior=vb_influence_results$log_prior,
                               worst_case_u=worst_case_u[, ind],
                               influence=influence_fun[, ind],
                               gbar=log_q_grad_terms[ind, ], method="lrvb",
                               metric=metric))
  return(vb_inf_df)
}




ind <- which(moment_df$par == "beta")
vb_inf_df <- GetVBInfluenceFunctionsDataFrame(beta_influence_results, ind, metric=sprintf("beta_on_%d", ind))
mcmc_inf_df <- GetMCMCInfluenceFunctionsDataFrame(
  mcmc_env$beta_mcmc_influence_funs, ind, BetaLogPrior, metric=sprintf("beta_on_%d", ind))


###########################
# save for paper

save(results, sens_results, worst_case_df, vb_inf_df, mcmc_inf_df,
     file=file.path(project_directory, "data", save_filename))

stop("Graphs follow, not executing")

# Basics




##################################
# Look at influence functions


grid.arrange(
  ggplot() +
    geom_point(data=vb_inf_df, aes(x=u, y=influence, color="lrvb")) +
    geom_point(data=mcmc_inf_df, aes(x=u, y=influence, color="mcmc"))
  ,
  ggplot() +
    geom_point(data=vb_inf_df, aes(x=u, y=gbar, color="lrvb")) +
    geom_point(data=mcmc_inf_df, aes(x=u, y=gbar, color="mcmc"))
  ,
  ggplot() +
    geom_point(data=vb_inf_df, aes(x=u, y=exp(log_dens), color="lrvb")) +
    geom_point(data=mcmc_inf_df, aes(x=u, y=exp(log_dens), color="mcmc"))
  ,
  # ggplot() +
  #   geom_point(data=vb_inf_df, aes(x=u, y=exp(log_prior), color="lrvb")) +
  #   geom_point(data=mcmc_inf_df, aes(x=u, y=exp(log_prior), color="mcmc"))
  # ,
  ggplot() +
    geom_point(data=vb_inf_df, aes(x=u, y=worst_case_u, color="lrvb")) +
    geom_point(data=mcmc_inf_df, aes(x=u, y=worst_case_u, color="mcmc")) +
    geom_point(data=mcmc_inf_df, aes(x=u, y=exp(log_prior), color="prior"))
  , ncol=4
)

# The prior is so diffuse that this looks crazy.
num_prior_points <- 500
prior_draws <- qinvgamma(seq(1 / (num_prior_points + 1), 1 - 1 / (num_prior_points + 1), length.out=num_prior_points),
                         shape=data_dat$beta_prior_shape, scale=data_dat$beta_prior_scale)
prior_draws <- seq(0, 2, length.out=num_prior_points)
log_prior_vals <- BetaLogPrior(prior_draws)
ggplot() +
  geom_point(data=vb_inf_df, aes(x=u, y=worst_case_u, color="lrvb")) +
  geom_point(data=mcmc_inf_df, aes(x=u, y=worst_case_u, color="mcmc")) +
  geom_point(aes(x=prior_draws, y=exp(log_prior_vals)))


# Compare g terms
beta_draws <- mcmc_env$mcmc_draws$beta
g_draws <- beta_draws
e_g_cond <- beta_mcmc_funs$GetConditionalMeanDiff(g_draws)
mcmc_beta_mean <- filter(mcmc_env$mcmc_results, par=="beta", metric=="mean")$value
mfvb_beta_mean <- filter(moment_df, par=="beta", metric=="mean")$value

ggplot() +
  geom_point(aes(x=beta_influence_results$u_draws, y=beta_influence_results$log_q_grad_terms[ind, ], color="lrvb")) +
  geom_point(aes(x=beta_draws[mcmc_rows], y=e_g_cond[mcmc_rows], color="mcmc")) +
  geom_vline(aes(xintercept=mfvb_beta_mean, color="lrvb")) +
  geom_vline(aes(xintercept=mcmc_beta_mean, color="mcmc"))


# Sanity checking that the densities are actually densitites
vb_dens <- exp(BetaLogQ(beta_influence_results$u_draws))
sorting_order <- order(beta_influence_results$u_draws)
sum(diff(beta_influence_results$u_draws[sorting_order]) * vb_dens[sorting_order[-1]])
sorting_order <- order(mcmc_env$mcmc_draws$beta)
sum(diff(mcmc_env$mcmc_draws$beta[sorting_order]) * beta_mcmc_funs$dens_at_draws[sorting_order[-1]])

# Lambda influence
lambda_ind <- 3
vb_influence_results <- lambda_worst_case_list[[lambda_ind]]
lambda_star_draws <- mcmc_env$mcmc_draws$lambda[, lambda_ind] / mcmc_env$mcmc_draws$beta
mcmc_funs <- GetMCMCInfluenceFunctions(lambda_star_draws, LambdaStarLogPrior)

# ... on beta
vb_ind <- which(moment_df$par == "beta")
mcmc_influence <- mcmc_funs$GetMCMCInfluence(mcmc_env$mcmc_draws$beta)
mcmc_rows <- sample(1:length(mcmc_influence), 5000)

ggplot() +
  geom_point(aes(x=vb_influence_results$u_draws, y=vb_influence_results$influence_fun[, vb_ind], color="lrvb")) +
  geom_point(aes(x=lambda_star_draws[mcmc_rows], y=mcmc_influence[mcmc_rows], color="mcmc"))

# ... and on itself.
vb_ind <- with(moment_df, which(par == "lambda" & group == lambda_ind))
mcmc_influence <- mcmc_funs$GetMCMCInfluence(mcmc_env$mcmc_draws$lambda[, lambda_ind])
mcmc_rows <- sample(1:length(mcmc_influence), 5000)

ggplot() +
  geom_point(aes(x=vb_influence_results$u_draws, y=vb_influence_results$influence_fun[, vb_ind], color="lrvb")) +
  geom_point(aes(x=lambda_star_draws[mcmc_rows], y=mcmc_influence[mcmc_rows], color="mcmc"))

# Compare densities
vb_dens <- exp(LambdaStarLogQ(vb_influence_results$u_draws, lambda_ind))
ggplot() +
  geom_point(aes(x=vb_influence_results$u_draws, y=vb_dens, color="lrvb")) +
  geom_point(aes(x=lambda_star_draws[mcmc_rows], y=mcmc_funs$dens_at_draws[mcmc_rows], color="mcmc"))


####################
# Graphs

ggplot(filter(worst_case_df, metric == "beta")) +
  geom_point(aes(x=gustafson, y=lrvb, color=paste(par, group)), size=3) +
  geom_abline(aes(intercept=0, slope=1))

ggplot(filter(worst_case_df, metric != "beta")) +
  geom_point(aes(x=gustafson, y=lrvb, color=metric), size=3) +
  geom_abline(aes(intercept=0, slope=1))

# Check Gustafson
ggplot(filter(worst_case_df, metric == "beta")) +
  geom_point(aes(x=gustafson, y=mcmc, color=paste(par, group)), size=3) +
  geom_abline(aes(intercept=0, slope=1))

ggplot(filter(worst_case_df, metric != "beta")) +
  geom_point(aes(x=gustafson, y=mcmc, color=paste(par, group)), size=3) +
  geom_abline(aes(intercept=0, slope=1))


ggplot(sens_results) +
  geom_point(aes(x=mcmc, y=lrvb, shape=par), size=3) +
  geom_abline(aes(intercept=0, slope=1)) +
  facet_grid(~ metric)
