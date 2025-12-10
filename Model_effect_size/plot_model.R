library(ggplot2)
library(brms)
library(bayesplot)
library(tidybayes)
library(dplyr)
library(here)

setwd(file.path(here::here(),"Model_effect_size/model"))
m <- readRDS("gamma_m_nov.RDS")

setwd(file.path(here::here(),"Model_effect_size/data"))
df<-read.csv("2025_11_19_153724_targ0.1-0.2-0.3-0.4_cl_size5_r_size5.csv")

df <- df[df$total_length > 0, ]
df <- df[which(df$total_cost>0),]
df$rat <- df$total_cost / df$total_length
names(df)

pred_data <- expand.grid(
  cost1 = seq(min(df$cost1), max(df$cost1), length.out = 20),
  c1_fraction = seq(min(df$c1_fraction),
                               max(df$c1_fraction), length.out = 20),
  replicate = unique(df$replicate)
)
df
library(tidyr)
spaghetti_pred <- fitted(
  m, 
  newdata = pred_data,
  re_formula = NULL,
  summary = FALSE
)

# Convert to long format: each row = one posterior draw for one row of pred_data
spaghetti_df <- cbind(pred_data, as.data.frame(t(spaghetti_pred)))


spaghetti_df <- spaghetti_df %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "draw",
    values_to = "estimate"
  )


i<-sample(seq(1,length(spaghetti_df$cost1),by=1),size=10000, replace=FALSE)
d<-sample(spaghetti_df[i,], )


p <- ggplot(d, aes(x = cost1, y = estimate,
                             group = interaction(c1_fraction, draw))) +
  geom_line(alpha = 0.4, linewidth = 0.3, color = "steelblue") +
  labs(x = "cost2",
       y = "total_cost / total_length",
       title = "Posterior Least-Cost Predictions Across Closed Substrate Levels") +
  theme_minimal() 
p


### Pop info

pred_data <- data.frame(
  cost1 = seq(min(df$cost1), max(df$cost1), length.out = 100),
  c1_fraction = mean(df$c1_fraction)
)

pop_pred <- fitted(
  m, newdata = pred_data, re_formula = NA,
  probs = c(0.10, 0.90)
)

pred_interval <- predict(
  m, newdata = pred_data, re_formula = NA,
  probs = c(0.10, 0.90)
)

plot_df <- data.frame(
  cost1 = pred_data$cost1,
  estimate = pop_pred[, "Estimate"],
  lower_mean = pop_pred[, "Q10"],
  upper_mean = pop_pred[, "Q90"],
  lower_pred = pred_interval[, "Q10"],
  upper_pred = pred_interval[, "Q90"]
)

p <- ggplot() +
  geom_point(data = df, aes(x = cost1, y = rat),
             alpha = 0.8, size = 1, color = "gray50") +
  geom_ribbon(data = plot_df, aes(x = cost1, ymin = lower_pred, ymax = upper_pred),
              alpha = 0.2, fill = "red") +
  geom_ribbon(data = plot_df, aes(x = cost1, ymin = lower_mean, ymax = upper_mean),
              alpha = 0.4, fill = "steelblue") +
  geom_line(data = plot_df, aes(x = cost1, y = estimate),
            linewidth = 1.5, color = "darkblue") +
  theme_bw() +
  labs(x = "cost1", y = "rat",
       title = "80% HDI & Prediction Intervals\n(at mean closed coverage)")
p


## Posterior distribution of the effect of cost
p <- m %>%
  as_draws_df() %>%
  ggplot(aes(x = b_cost1)) +
  geom_density(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(x = "cost2 coefficient",
       y = "Density",
       title = "Posterior Distribution of cost2 Effect")
p



library(tidybayes)
library(dplyr)
library(ggplot2)


### Plot the
# Choose closed fractions to visualize
c1_fraction_tape <- c(0.6, 0.7, 0.8, 0.9)

# Extract posterior draws
draws <- m %>% as_draws_df()

# Compute marginal effect of cost2 at each closed fraction
marginal_df <- lapply(c1_fraction_tape, function(cf) {
  tibble(
    c1_fraction = cf,
    marginal_effect = draws$b_cost1 + draws$`b_cost1:c1_fraction` * cf
  )
}) %>% bind_rows()


# Faceted density plot
p <- ggplot(marginal_df, aes(x = marginal_effect, group=c1_fraction, fill=c1_fraction)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    x = "Marginal effect of cost1",
    y = "Density",
    title = "Posterior Distribution of the Effect of cost2\nAcross Closed Substrate Fractions"
  )

p






library(tidyverse)
library(tidybayes)
library(ggridges)

c1_fraction_tape <- c(0.6, 0.7, 0.8, 0.9)
costs <- c(1, 1.5, 3.25, 5)   # include 1 explicitly as baseline

# Create prediction dataset
pred_grid <- expand_grid(
  c1_fraction = c1_fraction_tape,
  cost1 = costs
)

# Draw from posterior predictive expected values
draws_pred <- m %>%
  epred_draws(
    newdata = pred_grid,
    re_formula = NA
  )

# Compute effect sizes relative to cost = 1 *within each draw and closed level*
effects <- draws_pred %>%
  group_by(.draw, c1_fraction) %>%
  mutate(
    baseline = .epred[cost1 == 1],
    effect_size = .epred / baseline
  ) %>%
  ungroup()

# Make cost a nicer factor for ordering in ridge plot

effects <- effects %>%
  mutate(cost_label = paste0("cost = ", cost1),
         cost_label = factor(cost_label, levels = unique(cost_label)))

effects
# Ridge plot: densities of effect sizes, stacked by closed
ggplot(effects, aes(x = effect_size, y = cost_label, fill = cost_label)) +
  geom_density_ridges(alpha = 0.8, color = "white") +
  facet_wrap(~ c1_fraction, ncol = 1, scales = "free_y") +
  labs(
    x = "Effect size (relative to cost = 1)",
    y = "Cost",
    title = "Posterior Effect Sizes by Closed Level",
    subtitle = "Effect = RAT(closed, cost) / RAT(closed, 1)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
v

