library(ggplot2)
library(dplyr)
library(here)


setwd(file.path(here::here(),"Model_effect_size/data"))
df <- read.csv("2025_11_18_11targ0.6_cl_size5_r_size5targ0.7_cl_size5_r_size5targ0.8_cl_size5_r_size5targ0.9_cl_size5_r_size5.csv")

df <- df[df$total_length > 0, ]
df <- df[which(df$total_cost>0),]
df$rat <- df$total_cost / df$total_length



#Effects seem counterintuitive


df$rat<-df$cost2/df$total_length
ggplot(df, aes(x=rat,y=as.factor(cost2),fill=as.factor(cost2)))+geom_density_ridges()+
  facet_wrap(~target_closed_fraction)

df$rat<-df$cost2/df$total_length
ggplot(df, aes(x=rat,y=as.factor(target_closed_fraction),fill=as.factor(target_closed_fraction)))+geom_density_ridges()+
  facet_wrap(~cost2)


ggplot(df, aes(x=sinuosity,y=as.factor(target_closed_fraction),fill=as.factor(target_closed_fraction)))+geom_density_ridges()+
  facet_wrap(~cost2)