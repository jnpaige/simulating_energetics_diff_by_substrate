# plot_grids_openCoverage_cost1low_high.R
library(terra)
library(raster)
library(gdistance)
library(here)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)

source('/Users/jpaige/Desktop/Research repositories/simulating_energetics_diff_by_substrate/Simulating_substrates_and_cost/functions/functions_V4.R')

set.seed(133)

# =========================
# CONFIG (mirror sim_fast.R semantics)
# =========================
n <- 100
mod <- 5
c1_cluster_size <- 5
resource_cluster_size <- 5

cost0 <- 1          # CLOSED cost
cost1_low  <- 2     # OPEN cost (low)
cost1_high <- 6     # OPEN cost (high)

# You asked: vary OPEN coverage between .2 and .6 in the same conditions
open_coverages <- c(0.1, 0.7)

# IMPORTANT: your sim uses c1_fraction; we keep the name.
# Here we set c1_fraction to mean "open coverage" (to match your request wording),
# i.e., c1_fraction == open_coverage.
# If in your function c1_fraction actually means CLOSED fraction, change line below to:
#   c1_fraction = 1 - open_coverage
scenarios <- tibble::tribble(
  ~scenario,               ~open_coverage, ~cost1,
  "open0p2_costLOW",       0.1,            cost1_low,
  "open0p6_costLOW",       0.7,            cost1_low,
  "open0p2_costHIGH",      0.1,            cost1_high,
  "open0p6_costHIGH",      0.7,            cost1_high
)

out_dir <- file.path(getwd(), "example_plots_grids")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# raster template (cropped domain)
# =========================
n_crop <- n - 2 * mod
r_template <- raster::raster(
  nrows = n_crop, ncols = n_crop,
  xmn = mod, xmx = n - mod,
  ymn = mod, ymx = n - mod
)
r_template <- raster::setValues(r_template, rep(NA_real_, raster::ncell(r_template)))

max_tries <- 50

# store plots
sub_plots  <- vector("list", nrow(scenarios))
cost_plots <- vector("list", nrow(scenarios))

# =========================
# helper: run one sim with retries (inline logic, no wrapper functions)
# =========================
for (i in seq_len(nrow(scenarios))) {
  
  scenario_name <- scenarios$scenario[i]
  open_cov      <- scenarios$open_coverage[i]
  cost1         <- scenarios$cost1[i]
  
  # NOTE: set c1_fraction according to your requested meaning (open coverage).
  c1_fraction <- open_cov
  
  message(sprintf("[%d/%d] %s | open_coverage=%.2f | cost1(open)=%.2f",
                  i, nrow(scenarios), scenario_name, open_cov, cost1))
  
  out <- NULL
  for (attempt in seq_len(max_tries)) {
    out <- model_landscape_and_movement_fast(
      grid_size = n,
      mod = mod,
      cost0 = cost0,
      cost1 = cost1,
      c1_fraction = c1_fraction,
      c1_cluster_size = c1_cluster_size,
      resource_cluster_size = resource_cluster_size,
      raster_template = r_template
    )
    if (!is.null(out)) break
  }
  
  # unpack
  cost_surface <- out[[1]]
  lcp          <- out[[2]]
  xy_start     <- out[[3]]
  xy_target    <- out[[4]]
  substrate_r  <- out[[5]]
  travel_cost  <- out[[6]]
  
  if (inherits(cost_surface, "RasterLayer")) cost_surface <- rast(cost_surface)
  if (inherits(substrate_r,  "RasterLayer")) substrate_r  <- rast(substrate_r)
  if (inherits(travel_cost,  "RasterLayer")) travel_cost  <- rast(travel_cost)
  if (inherits(lcp, "SpatialLines") || inherits(lcp, "SpatialLinesDataFrame")) {
    lcp <- terra::vect(lcp)
  }
  
  # data frames
  df_sub <- as.data.frame(substrate_r, xy = TRUE, na.rm = FALSE)
  names(df_sub)[3] <- "substrate"
  
  df_tc <- as.data.frame(travel_cost, xy = TRUE, na.rm = FALSE)
  names(df_tc)[3] <- "travel_cost"
  
  df_cost <- as.data.frame(cost_surface, xy = TRUE, na.rm = FALSE)
  names(df_cost)[3] <- "acc_cost"
  
  lcp_df <- as.data.frame(terra::geom(lcp)) %>%
    filter(!is.na(x) & !is.na(y)) %>%
    mutate(part = as.factor(part)) %>%
    select(x, y, part)
  
  start_df  <- data.frame(x = xy_start[1],  y = xy_start[2])
  target_df <- data.frame(x = xy_target[1], y = xy_target[2])
  
  # classify substrate by mean travel cost vs cost0/cost1 (your clean fix)
  df_join <- df_sub %>%
    left_join(df_tc, by = c("x", "y")) %>%
    filter(!is.na(substrate), !is.na(travel_cost))
  
  tab_cost <- df_join %>%
    group_by(substrate) %>%
    summarise(mean_travel_cost = mean(travel_cost, na.rm = TRUE), .groups = "drop")
  
  open_value   <- tab_cost$substrate[which.min(abs(tab_cost$mean_travel_cost - cost1))]
  closed_value <- tab_cost$substrate[which.min(abs(tab_cost$mean_travel_cost - cost0))]
  
  df_sub <- df_sub %>%
    mutate(substrate_class = case_when(
      substrate == open_value   ~ "open substrate",
      substrate == closed_value ~ "closed substrate"
    ))
  
  title_txt <- sprintf("open_coverage=%.2f | cost1(open)=%.2f", open_cov, cost1)
  
  # ---- substrate plot ----
  sub_plots[[i]] <- ggplot() +
    geom_raster(data = df_sub, aes(x = x, y = y, fill = substrate_class), alpha = .5) +
    geom_path(data = lcp_df, aes(x = x, y = y, group = part), linewidth = 1.0, col = "red") +
    geom_point(data = start_df,  aes(x = x, y = y), shape = 21, size = 2.5) +
    geom_point(data = target_df, aes(x = x, y = y), shape = 21, size = 2.5) +
    scale_fill_manual(
      values = c("open substrate" = "gray85", "closed substrate" = "gray20"),
      name = "",
      breaks = c("open substrate", "closed substrate")
    ) +
    coord_equal() +
    theme_minimal() +
    theme_bw() +
    labs(title = title_txt, x = NULL, y = NULL)
  
  # ---- accumulated cost plot ----
  cost_plots[[i]] <- ggplot() +
    geom_raster(data = df_cost, aes(x = x, y = y, fill = acc_cost), alpha = .85) +
    geom_path(data = lcp_df, aes(x = x, y = y, group = part), linewidth = 1.0, col = "red") +
    geom_point(data = start_df,  aes(x = x, y = y), shape = 21, size = 2.5) +
    geom_point(data = target_df, aes(x = x, y = y), shape = 21, size = 2.5) +
    scale_fill_viridis_c(name = "acc. cost") +
    coord_equal() +
    theme_minimal() +
    theme_bw() +
    labs(title = title_txt, x = NULL, y = NULL)
}

# =========================
# order into 2x2 grids:
# top row = costLOW (open0.2, open0.6)
# bottom  = costHIGH (open0.2, open0.6)
# =========================
# =========================
# make 2x2 grids with tight spacing and panel labels
# =========================
grid_sub  <- (sub_plots[[1]] | sub_plots[[2]]) / (sub_plots[[3]] | sub_plots[[4]])
grid_cost <- (cost_plots[[1]] | cost_plots[[2]]) / (cost_plots[[3]] | cost_plots[[4]])

# add A–D tags and reduce padding
grid_sub  <- grid_sub  + plot_annotation(tag_levels = "A") &
  theme(plot.margin = margin(2, 2, 2, 1))

grid_cost <- grid_cost + plot_annotation(tag_levels = "A") &
  theme(plot.margin = margin(2, 2, 2, 1))

# save wide, high-res PNGs
sub_file  <- file.path(out_dir, "GRID_substrate_4scenarios.png")
cost_file <- file.path(out_dir, "GRID_accCost_4scenarios.png")

ggsave(sub_file,  grid_sub,  width = 8, height = 8, dpi = 600)
ggsave(cost_file, grid_cost, width = 8, height = 8, dpi = 600)

message("Saved grids:\n  ", sub_file, "\n  ", cost_file)

