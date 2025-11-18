library(ggplot2)
library(brms)
library(bayesplot)
library(tidybayes)
library(dplyr)
m <- readRDS("gamma_m_nov.RDS")
df<-read.csv("2025_11_11_20targ0.8_cl_size5_r_size5.csv")

df$rat<-df$total_cost/df$total_length
hist(df$rat)
# Create prediction data for plotting
pred_data <- expand.grid(
  cost2 = seq(min(df$cost2), max(df$cost2), length.out = 100),
  replicate = unique(df$replicate)
)

# Get predictions for spaghetti plot
spaghetti_pred <- fitted(m, newdata = pred_data, re_formula = NULL, 
                         summary = FALSE)  # Keep all posterior samples

# Convert to dataframe for plotting
spaghetti_df <- as.data.frame(t(spaghetti_pred))
spaghetti_df$cost2 <- pred_data$cost2
spaghetti_df$replicate <- pred_data$replicate

# Create spaghetti plot
spaghetti_plot <- ggplot(spaghetti_df, aes(x = cost2, y = V1, group = replicate)) +
  geom_line(alpha = 0.3, linewidth = 0.5, color = "steelblue") +
  labs(x = "cost2", y = "total_cost/total_length")+
  theme_minimal()

print(spaghetti_plot)
ggsave("spaghetti_plot_rat_cost.jpeg",spaghetti_plot, dpi=300)



library(ggplot2)
library(tidybayes)

# Create population-level predictions
pop_data <- data.frame(
  cost2 = seq(min(df$cost2), max(df$cost2), length.out = 100)
)

pop_pred <- fitted(m, newdata = pop_data, re_formula = NA, 
                   probs = c(0.1, 0.90))

pop_df <- data.frame(
  cost2 = pop_data$cost2,
  estimate = pop_pred[, "Estimate"],
  lower = pop_pred[, "Q10"],
  upper = pop_pred[, "Q90"]
)


# Create plot with raw data and population-level effect
p <- ggplot() +
  # Raw data points with some transparency to handle overplotting
  geom_point(data = df, aes(x = cost2, y = rat), 
             alpha = 0.3, size = 1, color = "gray50") +
  # HDI envelope
  geom_ribbon(data = pop_df, aes(x = cost2, ymin = lower, ymax = upper),
              alpha = 0.3, fill = "steelblue") +
  # Population mean line
  geom_line(data = pop_df, aes(x = cost2, y = estimate),
            linewidth = 1.5, color = "darkblue") +
  labs(x = "cost2", y = "rat (total_cost/total_length)")+
  theme_minimal()

p


library(ggplot2)
library(tidybayes)

# Create prediction data
pred_data <- data.frame(
  cost2 = seq(min(df$cost2), max(df$cost2), length.out = 100)
)

# Get population-level predictions (conditional mean) with 80% HDI
pop_pred <- fitted(m, newdata = pred_data, re_formula = NA, 
                   probs = c(0.10, 0.90))  # 80% HDI for population mean

pop_df <- data.frame(
  cost2 = pred_data$cost2,
  estimate = pop_pred[, "Estimate"],
  lower_mean = pop_pred[, "Q10"],
  upper_mean = pop_pred[, "Q90"]
)

# Get PREDICTION intervals (includes residual uncertainty)
pred_interval <- predict(m, newdata = pred_data, re_formula = NA,
                         probs = c(0.10, 0.90))  # 80% prediction interval

pred_df <- data.frame(
  cost2 = pred_data$cost2,
  lower_pred = pred_interval[, "Q10"],
  upper_pred = pred_interval[, "Q90"]
)

# Combine both
plot_df <- cbind(pop_df, pred_df[, c("lower_pred", "upper_pred")])

# Create plot with both 80% HDI for mean and 80% prediction interval
prediction_plot <- ggplot() +
  # Raw data points
  geom_point(data = df, aes(x = cost2, y = rat), 
             alpha = 0.8, size = 1, color = "gray50") +
  # 80% Prediction interval (includes residual uncertainty)
  geom_ribbon(data = plot_df, aes(x = cost2, ymin = lower_pred, ymax = upper_pred),
              alpha = 0.2, fill = "red", color = NA) +
  # 80% HDI for population mean
  geom_ribbon(data = plot_df, aes(x = cost2, ymin = lower_mean, ymax = upper_mean),
              alpha = 0.4, fill = "steelblue") +
  # Population mean line
  geom_line(data = plot_df, aes(x = cost2, y = estimate),
            linewidth = 1.5, color = "darkblue") +
  labs(x = "cost2", y = "rat (total_cost/total_length)")+
  theme_minimal()+theme_bw() +ggtitle("80% HDI population level effects of substrate cost\n on total cost/total length")

prediction_plot

ggsave("80_hdi_cost_on_rat.jpeg",prediction_plot,dpi=400)


# If you want to see the actual posterior distribution of the cost2 coefficient
posterior_plot <- m %>%
  as_draws_df() %>%
  ggplot(aes(x = b_cost2)) +
  geom_density(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "cost2 coefficient", y = "Density",
       title = "Posterior distribution of cost2 effect") +
  theme_minimal()

ggsave("cost2_posterior.jpeg",posterior_plot,dpi=400)

print(posterior_plot)
