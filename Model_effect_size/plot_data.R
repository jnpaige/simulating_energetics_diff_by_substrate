library(ggplot2)
library(dplyr)
library(here)
library(ggridges)

setwd(file.path(here::here(),"Model_effect_size/data"))
df <- read.csv("2025_12_10_113234_targ0.1-0.2-0.3-0.4-0.5-0.6-0.7-0.8-0.9_cl_size5_r_size5.csv")

df <- df[df$total_length > 0, ]
df <- df[which(df$total_cost>0),]
df$rat <- df$total_cost / df$euclidean


### RAT is cost per straight line distance.

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


library(ggplot2)
library(ggridges)

p<-ggplot(df, aes(x = rat, y = as.factor(cost1), fill = as.factor(cost1))) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(data = median_table,
             aes(xintercept = median_rat, color = as.factor(cost1)),
             linetype = "dashed", size = 0.6, show.legend = FALSE) +
  facet_wrap(~c1_fraction) +
  labs(x = "Cost / sinuosity", y = "Cost1",
       fill = "Cost1", color = "Cost1",
       title = "Density of rat values by cost1 with median reference lines") +
  theme_minimal()
p
ggsave("density of ratio by cost.png", p,  dpi=600, width=10, height=6)

#Ridge plots. Summarizing model outputs


p<-ggplot(df, aes(y=sinuosity, x=cost1))+facet_wrap(~c1_fraction) +geom_point(alpha=.04) +theme_bw() +
  labs(x = "Open substrate cost", y = "Sinuosity",
     title = "Sinuosity as a function of open substrate \n cost across open substrate proportion conditions") 

p
ggsave("open_cost_on_sinuosity_sim_output.png", dpi=600)

p<-ggplot(df, aes(y=rat, x=cost1))+facet_wrap(~c1_fraction) +geom_point(alpha=.04) +theme_bw() +
  labs(x = "Open substrate cost", y = "Cost/Euclidean Distance",
       title = "Cost/euclidian distance as a function of open substrate \n cost across open substrate proportion conditions") 
p
ggsave("open_cost_on_ratio_sim_output.png", dpi=600)


p<-ggplot(df, aes(x=cells_open, y=rat))+geom_point(alpha=.3)+
  facet_wrap(~c1_fraction) +theme_bw()

p

p<-ggplot(df, aes(x=cells_open, y=rat))+geom_point(alpha=.3)+
  facet_wrap(~c1_fraction) +theme_bw()
p

