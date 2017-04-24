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

mcmc_results_file <- paste(model_name, "mcmc.Rdata", sep="_")

#################
# Load or compile the stan model.

rdata_file <- file.path(data_directory, paste(model_name, "stan.Rdata", sep="_")) 
if (file.exists(rdata_file)) {
  print("Loading from file.")
  load(rdata_file)
} else {
  print("Compiling model.")
  mcmc_model <- stan_model(file.path(project_directory, "mcmc_stan/pump_model_mcmc.stan"))
  mcmc_log_prior <- stan_model(file.path(project_directory, "mcmc_stan/pump_model_mcmc_log_prior.stan"))

  vb_means <- stan_model(file.path(project_directory, "vb_transform_stan/vb_means.stan"))
  vb_model <- stan_model(file.path(project_directory, "vb_transform_stan/pump_model_vb.stan"))
  vb_log_prior <- stan_model(file.path(project_directory, "vb_transform_stan/pump_model_log_prior_vb.stan"))
  vb_log_density <- stan_model(file.path(project_directory, "vb_transform_stan/vb_density.stan"))
  
  save(vb_model, vb_means, vb_log_prior, vb_log_density,
       mcmc_model, mcmc_log_prior, file=rdata_file)
}


#################
# Load the data.

pump_data <- read.csv(file.path(data_directory, "pump_failure_data.csv"))
data_dat <- list(N=nrow(pump_data), s=pump_data$failures, t=pump_data$thousand_hours,
                 beta_prior_shape=0.01, beta_prior_scale=1.0, omega_0=1.802)
gustafson_result_original_df <-
  read.csv(file.path(data_directory, "gustafson_worst_case.csv")) %>%
  rename(prior=to.of) %>%
  melt(id.vars="prior") %>%
  rename(moment=variable) %>%
  mutate(par=gsub("[0-9]*", "", moment)) %>%
  mutate(par=sub("l", "lambda", par)) %>%
  mutate(group=ifelse(par=="beta", -1, as.numeric(sub("l", "", moment)))) %>%
  mutate(metric=sub("l", "lambda_star", prior)) %>%
  mutate(metric=sub("([0-9]+)", "[\\1]", metric, fixed=FALSE))

# Get Gustafson's results in our format.
gustafson_result_df <-
  GetParamRow(value=gustafson_result_original_df$value,
              par=gustafson_result_original_df$par,
              group=gustafson_result_original_df$group,
              metric=gustafson_result_original_df$metric,
              method="gustafson")


################
# Run MCMC

fit <- sampling(mcmc_model, data_dat, iter=500000, chains=1)
mcmc_draws <- extract(fit)
mcmc_means <- get_posterior_mean(fit)
lambda_means <- sapply(1:data_dat$N, function(ind) { mean(mcmc_draws$lambda[, ind]) } )
lambda_sds <- sapply(1:data_dat$N, function(ind) { sd(mcmc_draws$lambda[, ind]) } )
pump_data$mcmc_mean <- lambda_means
pump_data$mcmc_sd <- lambda_sds

if (FALSE) {
  ggplot(pump_data) +
    geom_point(aes(x=gustafson_mean, y=mcmc_mean)) +
    geom_abline(aes(intercept=0, slope=1))
  ggplot(pump_data) +
    geom_point(aes(x=gustafson_sd, y=mcmc_sd)) +
    geom_abline(aes(intercept=0, slope=1))
}

mcmc_results <- rbind(
  do.call(rbind, lapply(1:data_dat$N, function(ind) { GetMCMCLambdaGroupResults(ind, mcmc_draws) } )),
  GetMCMCBetaResults(mcmc_draws))


######################
# MCMC sensitivity

num_draws <- length(mcmc_draws$beta)
log_prior_grad_mat <- matrix(NA, nrow=num_draws, ncol=3)
mcmc_log_prior_fitobj <- sampling(mcmc_log_prior, list(N=data_dat$N, t=data_dat$t), iter=1, chains=1)
# These should be the prior param names
prior_ind <- 1:3
prior_param_names <- mcmc_log_prior_fitobj@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)[prior_ind]

prog_bar <- txtProgressBar(min=1, max=num_draws, style=3)
for (draw in 1:length(mcmc_draws$beta)) {
  setTxtProgressBar(prog_bar, value=draw)
  mcmc_draw_dat <- list(lambda=mcmc_draws$lambda[draw, ], beta=mcmc_draws$beta[draw],
                        beta_prior_shape=data_dat$beta_prior_shape, beta_prior_scale=data_dat$beta_prior_scale,
                        omega_0=data_dat$omega_0)
  prior_pars_free <- unconstrain_pars(mcmc_log_prior_fitobj, mcmc_draw_dat)
  log_prior_grad_mat[draw, ] <- grad_log_prob(mcmc_log_prior_fitobj, prior_pars_free, adjust_transform=FALSE)[prior_ind]
}
close(prog_bar)

mcmc_sensitivity_df <- rbind(
  GetMCMCLambdaSensitivityResults(mcmc_draws=mcmc_draws, log_prior_grad_mat, prior_param_names),
  GetMCMCBetaSensitivityResults(mcmc_draws=mcmc_draws, log_prior_grad_mat, prior_param_names))



##############################
# Our own MCMC versions of the worst case analysis

# Beta
beta_mcmc_influence_funs <- GetMCMCInfluenceFunctions(mcmc_draws$beta, BetaLogPrior)
GetMCMCWorstCase <- beta_mcmc_influence_funs$GetMCMCWorstCase
mcmc_beta_worst_case_df <-GetMCMCWorstCaseResults(GetMCMCWorstCase, "beta")

# Lambda
lambda_mcmc_influence_funs_list <- list()
lambda_worst_case_list <- list()
for (prior_ind in 1:data_dat$N) {
  param_draws <- mcmc_draws$lambda[, prior_ind] / mcmc_draws$beta
  lambda_influence_funs <- GetMCMCInfluenceFunctions(param_draws, LambdaStarLogPrior)
  lambda_mcmc_influence_funs_list[[prior_ind]] <- lambda_influence_funs
  GetMCMCWorstCase <- lambda_influence_funs$GetMCMCWorstCase
  prior_metric <- paste("lambda_star[", prior_ind, "]", sep="")
  lambda_worst_case_list[[length(lambda_worst_case_list) + 1]] <-
    GetMCMCWorstCaseResults(GetMCMCWorstCase, prior_metric)
}

mcmc_worst_case_df <- rbind(do.call(rbind, lambda_worst_case_list), mcmc_beta_worst_case_df)

#######################################
# Save

save(
  mcmc_results,
  mcmc_sensitivity_df,
  mcmc_worst_case_df,
  gustafson_result_df,
  mcmc_draws,
  beta_mcmc_influence_funs,
  lambda_worst_case_list,
  file=file.path(data_directory, "mcmc_results_file"))
