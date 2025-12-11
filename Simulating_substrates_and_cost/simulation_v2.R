# sim_fast.R
library(terra)
library(raster)
library(gdistance)
library(here)
library(doParallel)    # optional for parallel runs
library(foreach)


# Source the fast functions file
source('/Users/jpaige/Desktop/Research repositories/simulating_energetics_diff_by_substrate/Simulating_substrates_and_cost/functions/functions_V4.R') 

# ----------- parameters -------------
n <- 100
mod <- 5
k <- 100
cost_tape0 <- c(1)
cost_tape1 <- c(1, 1.5, 2, 2.5,3,3.5,4,4.5,5,5.5,6)
c1_fraction <- c(.1, .2, .3, .4,.5,.6,.7,.8,.9)  # earlier you used 'c1_fraction' naming
c1_cluster_size <- 5
resource_cluster_size <- 5

# metadata
metadata <- paste0("targ", paste(c1_fraction, collapse="-"),
                   "_cl_size", c1_cluster_size, "_r_size", resource_cluster_size)

# parameter grid
param_grid <- expand.grid(cost0 = cost_tape0,
                          cost1 = cost_tape1,
                          c1_fraction = c1_fraction,
                          replicate = seq_len(k),
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)
param_grid$metadata <- metadata
param_grid$total_length <- NA_real_
param_grid$euclidean <- NA_real_
param_grid$sinuosity <- NA_real_
param_grid$total_cost <- NA_real_
param_grid$cells_open <- NA_real_
param_grid$cells_closed <- NA_real_


# ----------- precompute the raster_template for the CROPPED domain -------------
# cropped domain dimensions
n_crop <- n - 2 * mod
if (n_crop <= 0) stop("mod too large relative to grid size")

# build raster::RasterLayer template that matches the cropped extent
r_template <- raster::raster(nrows = n_crop, ncols = n_crop,
                             xmn = mod, xmx = n - mod, ymn = mod, ymx = n - mod)
# set an initial value (not necessary) but ensures geometry is correct
r_template <- raster::setValues(r_template, rep(NA_real_, raster::ncell(r_template)))

# ---------- optional: parallel setup ----------
use_parallel <- FALSE
if (use_parallel) {
  cores <- parallel::detectCores()
  cl <- makeCluster(max(1, cores - 1))
  doParallel::registerDoParallel(cl)
}

# ---------- main loop (can be run with foreach %dopar% if use_parallel = TRUE) ----------
iter_seq <- seq_len(nrow(param_grid))

if (!use_parallel) {
  for (i in iter_seq) {
    message(sprintf("iter %d/%d : cost0=%g cost1=%g c1_fraction=%g replicate=%d",
                    i, nrow(param_grid), param_grid$cost0[i], param_grid$cost1[i],
                    param_grid$c1_fraction[i], param_grid$replicate[i]))
    
    out <- model_landscape_and_movement_fast(grid_size = n,
                                             mod = mod,
                                             cost0 = param_grid$cost0[i],
                                             cost1 = param_grid$cost1[i],
                                             c1_fraction = param_grid$c1_fraction[i],
                                             c1_cluster_size = c1_cluster_size,
                                             resource_cluster_size = resource_cluster_size,
                                             raster_template = r_template)
    
    if (is.null(out)) next
    
    # calculate metrics
    metrics <- calculate_path_metrics(out[[1]], out[[2]],
                                      substrate_r = out[[5]],
                                      travel_cost_r = out[[6]])
    
    param_grid$total_length[i] <- metrics$total_length
    param_grid$euclidean[i] <- metrics$euclidean
    param_grid$sinuosity[i] <- metrics$sinuosity
    param_grid$total_cost[i] <- metrics$total_cost
    param_grid$cells_open[i] <- metrics$cells_open
    param_grid$cells_closed[i] <- metrics$cells_closed
    
  }
} else {
  # parallel version (returns list rows -> combine back)
  res_list <- foreach(i = iter_seq, .packages = c("terra","raster","gdistance","sp","sf")) %dopar% {
    out <- model_landscape_and_movement_fast(grid_size = n,
                                             mod = mod,
                                             cost0 = param_grid$cost0[i],
                                             cost1 = param_grid$cost1[i],
                                             c1_fraction = param_grid$c1_fraction[i],
                                             c1_cluster_size = c1_cluster_size,
                                             resource_cluster_size = resource_cluster_size,
                                             raster_template = r_template)
    if (is.null(out)) return(NULL)
    metrics <- calculate_path_metrics(out[[1]], out[[2]], substrate_r = out[[3]])
    data.frame(total_length = metrics$total_length,
               euclidean = metrics$euclidean,
               sinuosity = metrics$sinuosity,
               total_cost = metrics$total_cost,
               i = i)
  }
  # merge back results
  for (r in res_list) {
    if (is.null(r)) next
    row_i <- r$i
    param_grid$total_length[row_i] <- r$total_length
    param_grid$euclidean[row_i] <- r$euclidean
    param_grid$sinuosity[row_i] <- r$sinuosity
    param_grid$total_cost[row_i] <- r$total_cost
  }
  stopCluster(cl)
}

#save results
file_name <- paste0(format(Sys.time(), "%Y_%m_%d_%H%M%S_"), metadata, ".csv")
write.csv(param_grid, file_name, row.names = FALSE)
message("Saved results to: ", file_name)


library(ggplot2)
library(tidybayes)
library(ggridges)
#Basic plots. 
df<-param_grid
length(df$cost0)


df <- df[df$total_length > 0, ]
df <- df[which(df$total_cost>0),]
df$rat <- df$total_cost / df$total_length


ggplot(df, aes(x=sinuosity, y=cells_open))+geom_point()+facet_wrap(~c1_fraction)
ggplot(df, aes(x=total_cost, y=total_length))+geom_point()+facet_wrap(~c1_fraction)
ggplot(df, aes(x=total_cost, y=total_length))+geom_point()+facet_wrap(~cost1)


ggplot(df, aes(x=rat,y=as.factor(c1_fraction),fill=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1)

ggplot(df, aes(x=rat,y=as.factor(cost1),fill=as.factor(cost1)))+geom_density_ridges()+
  facet_wrap(~c1_fraction)


ggplot(df, aes(x=sinuosity,y=as.factor(c1_fraction),fill=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1)

library(dplyr)
med1 <- df %>%
  group_by(c1_fraction, cost1) %>%
  summarize(median = median(rat, na.rm = TRUE), .groups = "drop")

ggplot(df, aes(x = rat, y = as.factor(cost1), fill = as.factor(cost1))) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(
    data = med1,
    aes(xintercept = median),
    color = "red", linewidth = 0.7
  ) +
  facet_wrap(~ c1_fraction) +
  theme_ridges() +
  theme(legend.position = "none") +ggtitle()




# Compute per-facet medians
med2 <- df %>%
  group_by(cost1, c1_fraction) %>%
  summarize(median = median(sinuosity, na.rm = TRUE), .groups = "drop")

ggplot(df, aes(x = sinuosity, y = as.factor(c1_fraction), fill = as.factor(c1_fraction))) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(
    data = med2,
    aes(xintercept = median),
    color = "red", linewidth = 0.7
  ) +
  facet_wrap(~ cost1) +
  theme_ridges() +
  theme(legend.position = "none")



