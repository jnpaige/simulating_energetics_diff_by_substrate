library(ggplot2)
library(dplyr)
library(here)
library(ggridges)

setwd(file.path(here::here(),"Model_effect_size/data"))
df <- read.csv("2025_12_02_121620_targ0.1-0.2-0.3-0.4-0.5-0.6-0.7-0.8-0.9_cl_size5_r_size5.csv")

df <- df[df$total_length > 0, ]
df <- df[which(df$total_cost>0),]
df$rat <- df$total_cost / df$euclidean


library(dplyr)

median_table <- df %>%
  group_by(c1_fraction, cost1) %>%
  summarize(
    median_rat = median(rat, na.rm = TRUE),
    median_total_cost = median(total_cost, na.rm = TRUE),
    .groups = "drop"
  )

median_table

ggplot(median_table, aes(x=median_rat,y=median_total_cost,col=cost1))+
  geom_point() +
  #facet_wrap(~cost1)+
  geom_smooth(method="lm")

median_table
library(ggplot2)
library(ggridges)

ggplot(df, aes(x = rat, y = as.factor(cost1), fill = as.factor(cost1))) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(data = median_table,
             aes(xintercept = median_rat, color = as.factor(cost1)),
             linetype = "dashed", size = 0.6, show.legend = FALSE) +
  facet_wrap(~c1_fraction) +
  labs(x = "Cost / Length (rat)", y = "Cost1",
       fill = "Cost1", color = "Cost1",
       title = "Density of rat values by cost1 with median reference lines") +
  theme_minimal()







#Effects seem counterintuitive


ggplot(df, aes(x=rat,y=as.factor(cost1),fill=as.factor(cost1)))+geom_density_ridges()+
  facet_wrap(~c1_fraction)


ggplot(df, aes(x=total_cost,y=as.factor(cost1),fill=as.factor(cost1)))+geom_density_ridges()+
  facet_wrap(~c1_fraction)

ggplot(df, aes(x=sinuosity,y=as.factor(cost1),fill=as.factor(cost1)))+geom_density_ridges()+
  facet_wrap(~c1_fraction)





