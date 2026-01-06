library(here)
library(brms)
library(cmdstanr)
library(ggplot2)
library(ggridges)
setwd(paste(here::here(),"/Model_effect_size/data", sep="",collapse=""))

df<-read.csv("2025_12_10_113234_targ0.1-0.2-0.3-0.4-0.5-0.6-0.7-0.8-0.9_cl_size5_r_size5.csv")
df$rat<-df$total_cost/df$euclidean

hist(log(df$rat))
hist(df$sinuosity)
##
## 1. Define formulas with appropriate distributions
##

f1 <- bf(
  sinuosity ~ cost1 +(1+cost1|c1_fraction),
  family = Gamma(link="log"))


f2 <- bf(
  rat ~ cost1 +(1+cost1|c1_fraction),
  family = Gamma(link="log"))
    
    

##
## 2. Priors
##

priors <- c(
  
  # Population-level effects (log scale)
  set_prior("normal(0, 1)", class = "Intercept", resp = "sinuosity"),
  set_prior("normal(0, 1)", class = "b", resp = "sinuosity"),
  
  set_prior("normal(4, 2)", class = "Intercept", resp = "rat"),
  set_prior("normal(0, 1)", class = "b", resp = "rat"),
  
  # Group-level standard deviations
  set_prior("exponential(1)", class = "sd", resp = "sinuosity"),
  set_prior("exponential(1)", class = "sd", resp = "rat"),
  
  # Correlation structure (LKJ prior)
  set_prior("lkj(2)", class = "cor"),
  
  # Gamma shape parameters
  set_prior("exponential(1)", class = "shape", resp = "sinuosity"),
  set_prior("exponential(1)", class = "shape", resp = "rat")
)


##
## 3. Fit the model
##

fit <- brm(
  f1 + f2,
  data = df,
  prior = priors,
  chains = 4, cores = 4, iter = 4000,
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)

summary(fit)
saveRDS(fit, "multivariate_bayes_model_dec_23.rds")
getwd()


# PRIOR-ONLY MODEL
fit_prior <- brm(
  f1 + f2,
  data = df,
  prior = priors,
  sample_prior = "only",
  chains = 4,
  cores = 4,
  iter = 2000,
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)

saveRDS(fit_prior, "multivariate_bayes_model_prior_only.rds")
