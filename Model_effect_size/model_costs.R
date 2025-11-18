library(here)
library(brms)
library(cmdstanr)
library(ggplot2)
setwd(paste(here::here(),"/Model_effect_size/data", sep="",collapse=""))

df<-read.csv("2025_11_18_11targ0.6_cl_size5_r_size5targ0.7_cl_size5_r_size5targ0.8_cl_size5_r_size5targ0.9_cl_size5_r_size5.csv")
head(df)
df$rat<-df$total_cost/df$total_length

ggplot(df, aes(x=total_cost/total_length, group=cost2,fill=cost2))+geom_density(alpha=.2)
ggplot(df, aes(x=sinuosity, group=cost2,fill=cost2))+geom_density(alpha=.2)

d<-rgamma(10000, 3,.01)
hist(d)
#hist(df$total_cost/df$total_length)
#hist(df$sinuosity)
names(df)
#A gamma distribution is most appropriate. 
m <- brm(
  formula = rat ~ cost2 * target_closed_fraction + (1 + cost2 | replicate),
  data = df,
  family = Gamma(link=log),
  cores = 4,
  chains = 4,
  seed = 123,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  backend = "cmdstanr"
)

summary(m)

setwd(paste(here::here(),"/Model_effect_size/model", sep="",collapse=""))

saveRDS(m, "gamma_m_nov.RDS")
m <- readRDS("gamma_m_nov.RDS")


