library(ggplot2)
library(brms)
library(bayesplot)
library(tidybayes)
library(dplyr)
library(here)

setwd(file.path(here::here(),"Model_effect_size/model"))
m <- readRDS("multivariate_bayes_model_dec.RDS")

setwd(file.path(here::here(),"Model_effect_size/data"))
df<-read.csv("2025_12_10_113234_targ0.1-0.2-0.3-0.4-0.5-0.6-0.7-0.8-0.9_cl_size5_r_size5.csv")




library(brms)
library(ggplot2)
library(dplyr)

# Population-level effect of cost1 on rat
ce_rat <- conditional_effects(
  m,
  resp       = "rat",
  effects    = "cost1",
  re_formula = NA    # drop group-level terms
)

# Population-level effect of cost1 on sinuosity
ce_sin <- conditional_effects(
  m,
  resp       = "sinuosity",
  effects    = "cost1",
  re_formula = NA
)

# These are already ggplot objects in a list of length 1
p<-plot(ce_rat, points = FALSE)[[1]] +
  labs(
    title = "Population-level effect of open substrate cost on cost/euclidean distance",
    x = "relative cost of open substrate",
    y = "Predicted cost/euclidean distance (Gamma mean)"
  )


p<-p+theme_bw()
p
ggsave("pop_effect_substrate_rat.png",p,dpi=600)


p<-plot(ce_sin, points = FALSE)[[1]] +
  labs(
    title = "Population-level effect of open substrate cost on sinuosity",
    x = "relative cost of open substrate",
    y = "Predicted sinuosity (Gamma mean)"
  )

p<-p+theme_bw()
p
ggsave("pop_effect_substrate_sinu.png",p,dpi=600)



### Focus on subset of open area costs at three levels of open area 

cost_seq    <- seq(min(df$cost1), max(df$cost1), length.out = 100)
open_levels <- c(0.1, 0.4, 0.7)   # 10%, 40%, 70% open area

newdat <- expand.grid(
  cost1       = cost_seq,
  c1_fraction = open_levels
)



f_rat <- fitted(
  m,
  newdata    = newdat,
  resp       = "rat",
  re_formula = NULL,   # include group-level terms
  summary    = TRUE,   # gives mean and quantiles
  probs      = c(0.055, 0.945)  # ~89% CI (or change to c(0.025,0.975) for 95%)
)

rat_pred <- cbind(newdat, as.data.frame(f_rat)) %>%
  rename(
    rat_est   = Estimate,
    rat_lower = Q5.5,   # or Q2.5 if you used probs = c(0.025, 0.975)
    rat_upper = Q94.5
  )


f_sin <- fitted(
  m,
  newdata    = newdat,
  resp       = "sinuosity",
  re_formula = NULL,
  summary    = TRUE,
  probs      = c(0.055, 0.945)
)

sin_pred <- cbind(newdat, as.data.frame(f_sin)) %>%
  rename(
    sin_est   = Estimate,
    sin_lower = Q5.5,
    sin_upper = Q94.5
  )


library(tidyr)

rat_long <- rat_pred %>%
  mutate(
    response = "rat",
    est      = rat_est,
    lower    = rat_lower,
    upper    = rat_upper
  ) %>%
  select(cost1, c1_fraction, response, est, lower, upper)

sin_long <- sin_pred %>%
  mutate(
    response = "sinuosity",
    est      = sin_est,
    lower    = sin_lower,
    upper    = sin_upper
  ) %>%
  select(cost1, c1_fraction, response, est, lower, upper)

plot_dat <- bind_rows(rat_long, sin_long)

p<-ggplot(
  plot_dat,
  aes(
    x    = cost1,
    y    = est,
    ymin = lower,
    ymax = upper,
    colour = factor(c1_fraction),
    fill   = factor(c1_fraction)
  )
) +
  geom_ribbon(alpha = 0.18, colour = NA) +
  geom_line() +
  facet_wrap(~ response, scales = "free_y") +
  labs(
    x        = "Relative cost of open substrate",
    y        = "Posterior mean (Gamma scale)",
    colour   = "Open area fraction",
    fill     = "Open area fraction",
    title    = "Effect of relative cost of open substrate \n on transport cost/euclidean distance and sinuosity\nat 10%, 40%, 70% open area"
  ) +
  theme_bw()

p

ggsave("cond_effect_open_area.png",p,dpi=600)









