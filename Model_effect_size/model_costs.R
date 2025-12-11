library(here)
library(brms)
library(cmdstanr)
library(ggplot2)
library(ggridges)
setwd(paste(here::here(),"/Model_effect_size/data", sep="",collapse=""))

df<-read.csv("2025_12_10_113234_targ0.1-0.2-0.3-0.4-0.5-0.6-0.7-0.8-0.9_cl_size5_r_size5.csv")
head(df)
df$rat<-df$total_cost/df$euclidean

ggplot(df, aes(x=rat, group=cost1,fill=cost1))+geom_density(alpha=.2)+
  facet_wrap(~c1_fraction) +ggtitle("accumulated cost as a function of distance to \n target by open cost")

ggplot(df, aes(x=sinuosity, group=cost1,fill=cost1))+geom_density(alpha=.2) +
  ggtitle("Sinuousity in different open substrate cost conditions")


ggplot(df, aes(x=rat,y=as.factor(c1_fraction),fill=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1)

ggplot(df, aes(x=total_cost,y=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1) +
  ggtitle("Total cost/euclidian distance to target \n divided by open area cost\n and percent open area")  
  

ggplot(df, aes(x=sinuosity,y=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1) +ggtitle("Sinuosity divided by open area cost and percent open area")

hist(df$total_cost)
d<-rgamma(10000, 10,.0001)
hist(d)


hist(df$sinuosity)
##
## 1. Define formulas with appropriate distributions
##

hist(df$sinuousity)
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
saveRDS(fit, "multivariate_bayes_model_dec.rds")
getwd()


