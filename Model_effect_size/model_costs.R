library(here)
library(brms)
library(cmdstanr)
library(ggplot2)
library(ggridges)
setwd(paste(here::here(),"/Model_effect_size/data", sep="",collapse=""))

df<-read.csv("2025_11_19_153724_targ0.1-0.2-0.3-0.4_cl_size5_r_size5.csv")
head(df)
df$rat<-df$total_cost/df$total_length

ggplot(df, aes(x=rat, group=cost1,fill=cost1))+geom_density(alpha=.2)+facet_wrap(~c1_fraction)
ggplot(df, aes(x=sinuosity, group=cost1,fill=cost1))+geom_density(alpha=.2)


ggplot(df, aes(x=rat,y=as.factor(c1_fraction),fill=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1)

ggplot(df, aes(x=total_cost,y=as.factor(c1_fraction),fill=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1)

ggplot(df, aes(x=sinuosity,y=as.factor(c1_fraction),fill=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1)






d<-rgamma(10000, 3,.01)
hist(d)
#hist(df$total_cost/df$total_length)
#hist(df$sinuosity)
names(df)
#A gamma distribution is most appropriate. 
m <- brm(
  formula = rat ~ cost1 * c1_fraction + (1 + cost1 | replicate),
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


