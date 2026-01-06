library(here)
library(brms)
library(bayesplot)
library(ggplot2)

setwd(paste(here::here(), "/Model_effect_size/data", sep = "", collapse = ""))

fit_post  <- readRDS("multivariate_bayes_model_dec_23.rds")
fit_prior <- readRDS("multivariate_bayes_model_prior_only.rds")

df <- fit_post$data  # use ONE object name

hist(df$rat)

pp_rat <- posterior_predict(fit_prior, resp = "rat", ndraws = 500)
pp_sin <- posterior_predict(fit_prior, resp = "sinuosity", ndraws = 500)

# Confirm matrix dims align with data rows
print(dim(pp_rat))
print(nrow(df))


ppc_dens_overlay(y = df$rat, yrep = pp_rat) +
  ggtitle("Prior predictive check: rat")

ppc_dens_overlay(y = df$sinuosity, yrep = pp_sin) +
  ggtitle("Prior predictive check: sinuosity")


