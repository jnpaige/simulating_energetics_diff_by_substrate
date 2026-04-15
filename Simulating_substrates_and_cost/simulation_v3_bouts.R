# simulation_v3_bouts.R
# =============================================================================
# Bout-based simulation of sifaka movement on heterogeneous landscapes
# =============================================================================
#
# KEY CHANGES FROM V2:
# - Landscape is 3x3 km with 5m resolution (600x600 cells)
# - Instead of single start-to-resource LCPs, we simulate iterated BOUTS
# - Bout distances drawn from a Gamma distribution matching real sifaka GPS data
# - Each bout: random direction -> target point -> LCP -> record metrics
# - We record per-bout: LCP distance, euclidean distance, sinuosity,
#   cumulative cost, distance in open vs closed substrate
#
# PERFORMANCE ARCHITECTURE:
# - The transition matrix is built ONCE per landscape (~15-40 sec for 600x600)
# - Each shortestPath call is cheap (~0.01-0.1 sec) since bouts are short
# - Metric extraction uses pre-extracted matrices (no raster::extract per bout)
# - Total per landscape: ~1-3 min for 1000 bouts
# - ~150 landscapes total → estimated 2.5-7.5 hours
#
# For faster iteration during development, reduce:
#   - grid_size to 300 (1.5km @ 5m, ~4x faster transition build)
#   - n_bouts to 200
#   - k_replicates to 2
#
# =============================================================================

library(terra)
library(raster)
library(gdistance)
library(sf)
library(sp)
library(dplyr)
library(doParallel)
library(foreach)

# =============================================================================
# 1. BOUT DISTANCE DISTRIBUTION
# =============================================================================
# From real sifaka GPS data (Travel_bout_distance_density_plot.pdf):
#   - Most bouts very short (0-100m), right-skewed, tail to ~1500m
#   - Gamma(shape=1.2, rate=0.015) gives:
#     mean~80m, median~56m, 95th%~260m
#
# TODO: Replace with MLE fit from actual bout data:
#   library(fitdistrplus)
#   fit <- fitdist(bout_distances, "gamma")
#   bout_distance_shape <- fit$estimate["shape"]
#   bout_distance_rate  <- fit$estimate["rate"]

BOUT_DIST_SHAPE <- 1.2
BOUT_DIST_RATE  <- 0.015

sample_bout_distance <- function(n = 1) {
  d <- rgamma(n, shape = BOUT_DIST_SHAPE, rate = BOUT_DIST_RATE)
  pmax(d, 5)  # floor at 5m (one cell)
}


# =============================================================================
# 2. LANDSCAPE GENERATION
# =============================================================================

generate_landscape <- function(grid_size = 600,
                               cell_size = 5,
                               buffer_cells = 50,
                               cost0 = 1,
                               cost1 = 1,
                               c1_fraction = 0.3,
                               cluster_size_m = 50) {
  
  cluster_cells <- max(3, round(cluster_size_m / cell_size))
  if (cluster_cells %% 2 == 0) cluster_cells <- cluster_cells + 1
  
  extent_m <- grid_size * cell_size
  
  # Full grid
  land <- terra::rast(nrows = grid_size, ncols = grid_size,
                      xmin = 0, xmax = extent_m,
                      ymin = 0, ymax = extent_m)
  
  # Spatially autocorrelated binary substrate
  substrate <- terra::rast(land)
  terra::values(substrate) <- runif(terra::ncell(substrate))
  substrate <- terra::focal(substrate, w = cluster_cells, fun = mean, na.rm = TRUE)
  
  smn <- terra::global(substrate, "min", na.rm = TRUE)[[1]]
  smx <- terra::global(substrate, "max", na.rm = TRUE)[[1]]
  if (smx > smn) {
    substrate <- (substrate - smn) / (smx - smn)
  } else {
    terra::values(substrate) <- 0
  }
  
  v <- terra::values(substrate)
  v[is.na(v)] <- 0
  threshold <- stats::quantile(v, probs = c1_fraction, na.rm = TRUE)
  v_bin <- ifelse(v <= threshold, 1L, 0L)  # 1=open, 0=closed
  terra::values(substrate) <- v_bin
  
  # Travel cost raster
  vals_sub <- terra::values(substrate)
  vals_sub[is.na(vals_sub)] <- 0
  vals_cost <- ifelse(vals_sub == 1, cost1, cost0)
  
  travel_cost_terra <- terra::rast(land)
  terra::values(travel_cost_terra) <- vals_cost
  
  # Crop to interior
  buf_m <- buffer_cells * cell_size
  crop_ext <- terra::ext(buf_m, extent_m - buf_m,
                         buf_m, extent_m - buf_m)
  
  travel_cost_crop <- terra::crop(travel_cost_terra, crop_ext)
  substrate_crop   <- terra::crop(substrate, crop_ext)
  
  # Convert to raster::RasterLayer for gdistance
  interior_size <- grid_size - 2 * buffer_cells
  r_template <- raster::raster(nrows = interior_size, ncols = interior_size,
                               xmn = buf_m, xmx = extent_m - buf_m,
                               ymn = buf_m, ymx = extent_m - buf_m)
  
  travel_cost_r <- raster::setValues(r_template, terra::values(travel_cost_crop))
  substrate_r   <- raster::setValues(r_template, terra::values(substrate_crop))
  
  # ---------------------------------------------------------------
  # BUILD TRANSITION LAYER — the one expensive step per landscape
  # ---------------------------------------------------------------
  # Creates sparse ~NxN adjacency matrix (8-connected).
  # All subsequent shortestPath() calls reuse this.
  # ---------------------------------------------------------------
  tr <- gdistance::transition(1 / travel_cost_r, transitionFunction = mean, directions = 8)
  trC <- gdistance::geoCorrection(tr, type = "c")
  
  # Pre-extract substrate and cost as matrices for fast lookup
  sub_matrix  <- raster::as.matrix(substrate_r)
  cost_matrix <- raster::as.matrix(travel_cost_r)
  
  return(list(
    trC           = trC,
    travel_cost_r = travel_cost_r,
    substrate_r   = substrate_r,
    sub_matrix    = sub_matrix,
    cost_matrix   = cost_matrix,
    cell_size     = cell_size,
    interior_ext  = raster::extent(r_template),
    nrow          = interior_size,
    ncol          = interior_size
  ))
}


# =============================================================================
# 3. FAST METRIC EXTRACTION
# =============================================================================
# Uses pre-extracted matrices instead of raster::extract per bout.
# Converts path coordinates to row/col indices for direct matrix lookup.

extract_bout_metrics <- function(path_coords, landscape) {
  
  ext <- landscape$interior_ext
  nr  <- landscape$nrow
  nc  <- landscape$ncol
  cs  <- landscape$cell_size
  
  n_pts <- nrow(path_coords)
  if (n_pts < 2) return(NULL)
  
  # Segment lengths (vectorized)
  dx <- diff(path_coords[, 1])
  dy <- diff(path_coords[, 2])
  seg_lengths <- sqrt(dx^2 + dy^2)
  lcp_length  <- sum(seg_lengths)
  
  # Euclidean: first to last point
  euclid <- sqrt((path_coords[n_pts, 1] - path_coords[1, 1])^2 +
                   (path_coords[n_pts, 2] - path_coords[1, 2])^2)
  
  sinuosity <- if (euclid > 0) lcp_length / euclid else NA_real_
  
  # Convert coordinates to row/col for matrix lookup
  # raster row 1 = top (ymax), increases downward
  col_idx <- pmin(nc, pmax(1, floor((path_coords[, 1] - ext@xmin) / cs) + 1))
  row_idx <- pmin(nr, pmax(1, nr - floor((path_coords[, 2] - ext@ymin) / cs)))
  
  # Look up substrate and cost at each vertex
  sub_vals  <- landscape$sub_matrix[cbind(row_idx, col_idx)]
  cost_vals <- landscape$cost_matrix[cbind(row_idx, col_idx)]
  
  # Per-segment attribution (start vertex determines substrate/cost)
  n_segs  <- length(seg_lengths)
  is_open <- (sub_vals[1:n_segs] == 1)
  is_open[is.na(is_open)] <- FALSE
  
  dist_open   <- sum(seg_lengths[is_open])
  dist_closed <- sum(seg_lengths[!is_open])
  
  seg_costs <- cost_vals[1:n_segs]
  seg_costs[is.na(seg_costs)] <- 1
  total_cost <- sum(seg_costs * seg_lengths)
  
  n_cells_open   <- sum(sub_vals == 1, na.rm = TRUE)
  n_cells_closed <- sum(sub_vals == 0, na.rm = TRUE)
  
  list(
    lcp_length     = lcp_length,
    euclidean      = euclid,
    sinuosity      = sinuosity,
    total_cost     = total_cost,
    dist_open      = dist_open,
    dist_closed    = dist_closed,
    n_cells_open   = n_cells_open,
    n_cells_closed = n_cells_closed
  )
}


# =============================================================================
# 4. SINGLE BOUT
# =============================================================================

simulate_one_bout <- function(xy_start, bout_dist_m, landscape) {
  
  ext    <- landscape$interior_ext
  cs     <- landscape$cell_size
  margin <- cs * 3
  
  # Random direction
  theta <- runif(1, 0, 2 * pi)
  
  xy_target <- c(
    xy_start[1] + bout_dist_m * cos(theta),
    xy_start[2] + bout_dist_m * sin(theta)
  )
  
  # Clamp to interior
  xy_target[1] <- max(ext@xmin + margin, min(ext@xmax - margin, xy_target[1]))
  xy_target[2] <- max(ext@ymin + margin, min(ext@ymax - margin, xy_target[2]))
  
  euclid <- sqrt(sum((xy_target - xy_start)^2))
  if (euclid < cs) return(list(metrics = NULL, xy_end = xy_target))
  
  # LCP (fast on precomputed transition layer)
  lcp_sp <- tryCatch(
    gdistance::shortestPath(landscape$trC,
                            matrix(xy_start, nrow = 1),
                            matrix(xy_target, nrow = 1),
                            output = "SpatialLines"),
    error = function(e) NULL
  )
  if (is.null(lcp_sp)) return(list(metrics = NULL, xy_end = xy_target))
  
  # Extract path coords from sp object (no terra conversion needed)
  path_coords <- do.call(rbind, lapply(lcp_sp@lines, function(L) {
    do.call(rbind, lapply(L@Lines, function(ll) ll@coords))
  }))
  
  if (is.null(path_coords) || nrow(path_coords) < 2) {
    return(list(metrics = NULL, xy_end = xy_target))
  }
  
  # Metrics via matrix lookup
  m <- extract_bout_metrics(path_coords, landscape)
  if (is.null(m)) return(list(metrics = NULL, xy_end = xy_target))
  
  xy_end <- path_coords[nrow(path_coords), ]
  
  metrics <- data.frame(
    euclidean           = m$euclidean,
    lcp_length          = m$lcp_length,
    sinuosity           = m$sinuosity,
    total_cost          = m$total_cost,
    dist_open           = m$dist_open,
    dist_closed         = m$dist_closed,
    n_cells_open        = m$n_cells_open,
    n_cells_closed      = m$n_cells_closed,
    frac_dist_open      = m$dist_open / (m$dist_open + m$dist_closed + 1e-12),
    bout_dist_requested = bout_dist_m,
    start_x             = xy_start[1],
    start_y             = xy_start[2],
    end_x               = xy_end[1],
    end_y               = xy_end[2],
    stringsAsFactors    = FALSE
  )
  
  list(metrics = metrics, xy_end = xy_end)
}


# =============================================================================
# 5. ITERATED BOUTS ON ONE LANDSCAPE
# =============================================================================

simulate_bouts_on_landscape <- function(landscape, n_bouts = 1000,
                                        verbose = TRUE) {
  
  ext    <- landscape$interior_ext
  cs     <- landscape$cell_size
  margin <- cs * 10
  
  # Random start
  xy_current <- c(
    runif(1, ext@xmin + margin, ext@xmax - margin),
    runif(1, ext@ymin + margin, ext@ymax - margin)
  )
  
  results   <- vector("list", n_bouts)
  n_success <- 0L
  n_fail    <- 0L
  
  for (b in seq_len(n_bouts)) {
    bout_dist <- sample_bout_distance(1)
    bout <- simulate_one_bout(xy_current, bout_dist, landscape)
    
    if (!is.null(bout$metrics)) {
      bout$metrics$bout_id <- b
      results[[b]] <- bout$metrics
      n_success <- n_success + 1L
    } else {
      n_fail <- n_fail + 1L
    }
    
    xy_current <- bout$xy_end
    
    # Re-center if near edge
    if (xy_current[1] < ext@xmin + margin ||
        xy_current[1] > ext@xmax - margin ||
        xy_current[2] < ext@ymin + margin ||
        xy_current[2] > ext@ymax - margin) {
      xy_current <- c(
        runif(1, ext@xmin + margin, ext@xmax - margin),
        runif(1, ext@ymin + margin, ext@ymax - margin)
      )
    }
    
    if (verbose && b %% 200 == 0) {
      message(sprintf("  bout %d/%d (ok: %d, fail: %d)",
                      b, n_bouts, n_success, n_fail))
    }
  }
  
  bout_df <- do.call(rbind, results[!sapply(results, is.null)])
  
  if (verbose) {
    message(sprintf("  Complete: %d ok, %d failed", n_success, n_fail))
  }
  
  bout_df
}


# =============================================================================
# 6. MAIN SIMULATION LOOP (PARALLEL)
# =============================================================================
#
# Each landscape is fully independent, so we parallelize across landscapes.
# The transition matrix (~100MB per landscape) lives in each worker's memory.
# With 8 cores, expect ~8x speedup → ~20-40 min instead of 2.5-5 hours.
#
# Memory estimate: each worker holds one landscape at a time.
#   - 600x600 transition matrix (sparse): ~100-200 MB
#   - Substrate/cost matrices: ~3 MB
#   - Bout results: ~5 MB
#   Total per worker: ~200-300 MB
#   With 8 workers: ~2 GB (fine for most machines with 16+ GB RAM)
#
# If memory is tight, reduce n_cores below.
# =============================================================================

library(doParallel)
library(foreach)

# --- Parameters ---
grid_size      <- 600       # 600 cells × 5m = 3km
cell_size      <- 5         # meters per cell
buffer_cells   <- 50        # 250m buffer each side
cluster_size_m <- 50        # substrate patch size (meters)
n_bouts        <- 1000      # bouts per landscape
k_replicates   <- 5         # replicates per condition

cost0          <- 1         # closed (arboreal) cost — always 1
cost1_values   <- c(1, 1.25, 1.5, 2, 3, 4)
c1_fraction_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# --- FOR QUICK TESTING, uncomment these: ---
# grid_size    <- 300
# n_bouts      <- 200
# k_replicates <- 2
# cost1_values <- c(1, 2, 4)
# c1_fraction_values <- c(0.3, 0.7)

# --- Parallel setup ---
# Set to FALSE to run sequentially (useful for debugging)
use_parallel <- TRUE
n_cores      <- max(1, parallel::detectCores() - 1)  # leave 1 core free

param_grid <- expand.grid(
  cost1       = cost1_values,
  c1_fraction = c1_fraction_values,
  replicate   = seq_len(k_replicates),
  stringsAsFactors = FALSE
)

message(sprintf("=== Bout simulation v3 ==="))
message(sprintf("Grid: %dx%d @ %dm = %dx%dm landscape",
                grid_size, grid_size, cell_size,
                grid_size * cell_size, grid_size * cell_size))
message(sprintf("Interior: %dx%d cells after %dm buffer",
                grid_size - 2*buffer_cells, grid_size - 2*buffer_cells,
                buffer_cells * cell_size))
message(sprintf("Conditions: %d cost x %d fraction x %d reps = %d landscapes",
                length(cost1_values), length(c1_fraction_values),
                k_replicates, nrow(param_grid)))
message(sprintf("Bouts per landscape: %d", n_bouts))
message(sprintf("Total bouts: %d", nrow(param_grid) * n_bouts))
message(sprintf("Parallel: %s (%d cores)\n",
                ifelse(use_parallel, "YES", "NO"), n_cores))

# --- Worker function: one complete landscape + all its bouts ---
# This is what each core runs independently.
run_one_landscape <- function(i, param_grid, grid_size, cell_size,
                              buffer_cells, cost0, cluster_size_m, n_bouts) {
  
  c1_cost <- param_grid$cost1[i]
  c1_frac <- param_grid$c1_fraction[i]
  rep_id  <- param_grid$replicate[i]
  
  # Generate landscape + transition layer
  landscape <- tryCatch(
    generate_landscape(
      grid_size      = grid_size,
      cell_size      = cell_size,
      buffer_cells   = buffer_cells,
      cost0          = cost0,
      cost1          = c1_cost,
      c1_fraction    = c1_frac,
      cluster_size_m = cluster_size_m
    ),
    error = function(e) NULL
  )
  if (is.null(landscape)) return(NULL)
  
  # Simulate bouts (verbose = FALSE in parallel to avoid garbled output)
  bout_df <- simulate_bouts_on_landscape(landscape, n_bouts = n_bouts,
                                         verbose = FALSE)
  
  if (!is.null(bout_df) && nrow(bout_df) > 0) {
    bout_df$cost0       <- cost0
    bout_df$cost1       <- c1_cost
    bout_df$c1_fraction <- c1_frac
    bout_df$replicate   <- rep_id
    bout_df$landscape_id <- i
  }
  
  return(bout_df)
}

t_start_all <- proc.time()

if (use_parallel) {
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  # Export all functions and globals to workers
  parallel::clusterExport(cl, c(
    "generate_landscape",
    "simulate_bouts_on_landscape",
    "simulate_one_bout",
    "extract_bout_metrics",
    "sample_bout_distance",
    "run_one_landscape",
    "BOUT_DIST_SHAPE",
    "BOUT_DIST_RATE"
  ))
  
  # Run in parallel with progress tracking via .verbose
  all_results <- foreach(
    i = seq_len(nrow(param_grid)),
    .packages = c("terra", "raster", "gdistance", "sp"),
    .errorhandling = "pass",     # don't kill everything if one landscape fails
    .verbose = FALSE
  ) %dopar% {
    run_one_landscape(i, param_grid, grid_size, cell_size,
                      buffer_cells, cost0, cluster_size_m, n_bouts)
  }
  
  parallel::stopCluster(cl)
  
  # Check for errors
  errors <- sapply(all_results, function(x) inherits(x, "error") || inherits(x, "simpleError"))
  if (any(errors)) {
    message(sprintf("WARNING: %d/%d landscapes had errors:", sum(errors), length(errors)))
    for (j in which(errors)) {
      message(sprintf("  [%d] cost1=%.2f, frac=%.1f, rep=%d: %s",
                      j, param_grid$cost1[j], param_grid$c1_fraction[j],
                      param_grid$replicate[j],
                      conditionMessage(all_results[[j]])))
    }
    all_results[errors] <- list(NULL)
  }
  
} else {
  # Sequential fallback (for debugging)
  all_results <- vector("list", nrow(param_grid))
  
  for (i in seq_len(nrow(param_grid))) {
    message(sprintf("[%d/%d] cost1=%.2f | frac=%.1f | rep=%d",
                    i, nrow(param_grid),
                    param_grid$cost1[i], param_grid$c1_fraction[i],
                    param_grid$replicate[i]))
    
    t0 <- proc.time()
    all_results[[i]] <- run_one_landscape(
      i, param_grid, grid_size, cell_size,
      buffer_cells, cost0, cluster_size_m, n_bouts
    )
    t1 <- proc.time()
    
    nr <- if (!is.null(all_results[[i]])) nrow(all_results[[i]]) else 0
    message(sprintf("  %d bouts in %.1f sec", nr, (t1 - t0)[3]))
  }
}

t_end_all <- proc.time()
message(sprintf("\n=== Total runtime: %.1f min ===",
                (t_end_all - t_start_all)[3] / 60))

# =============================================================================
# 7. COMBINE AND SAVE
# =============================================================================

results_df <- do.call(rbind, all_results[!sapply(all_results, is.null)])
results_df$cost_per_euclid <- results_df$total_cost / results_df$euclidean

metadata <- sprintf("bout_sim_v3_%s_n%d_k%d",
                    format(Sys.time(), "%Y_%m_%d_%H%M%S"),
                    n_bouts, k_replicates)

file_name <- paste0(metadata, ".csv")
write.csv(results_df, file_name, row.names = FALSE)
message("Saved: ", file_name)
message(sprintf("Total rows: %d", nrow(results_df)))

# =============================================================================
# 8. DIAGNOSTIC PLOTS
# =============================================================================

library(ggplot2)

# Bout distance distribution
p1 <- ggplot(results_df, aes(x = euclidean)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(title = "Simulated bout Euclidean distances",
       x = "Euclidean distance (m)", y = "Count") +
  theme_bw()
ggsave("bout_v3_dist_histogram.png", p1, dpi = 300, width = 8, height = 5)

# Sinuosity by cost × open fraction
p2 <- ggplot(results_df, aes(x = sinuosity, fill = as.factor(cost1))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ c1_fraction, labeller = label_both) +
  labs(title = "Bout sinuosity by open substrate cost",
       x = "Sinuosity", fill = "Open cost") +
  theme_bw() +
  coord_cartesian(xlim = c(1, 2.5))
ggsave("bout_v3_sinuosity_density.png", p2, dpi = 300, width = 12, height = 8)

# Fraction distance in open
p3 <- ggplot(results_df, aes(x = frac_dist_open, fill = as.factor(cost1))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ c1_fraction, labeller = label_both) +
  labs(title = "Fraction of bout distance in open substrate",
       x = "Fraction open", fill = "Open cost") +
  theme_bw()
ggsave("bout_v3_frac_open_density.png", p3, dpi = 300, width = 12, height = 8)

# Cost per euclidean
p4 <- ggplot(results_df, aes(x = cost_per_euclid, fill = as.factor(cost1))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ c1_fraction, labeller = label_both) +
  labs(title = "Cost per Euclidean distance",
       x = "Cost / Euclidean", fill = "Open cost") +
  theme_bw()
ggsave("bout_v3_cost_per_euclid.png", p4, dpi = 300, width = 12, height = 8)

# Sinuosity vs distance — the key scale-dependent relationship
p5 <- ggplot(results_df %>% filter(c1_fraction %in% c(0.3, 0.5, 0.7)),
             aes(x = euclidean, y = sinuosity, color = as.factor(cost1))) +
  geom_point(alpha = 0.05, size = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ c1_fraction, labeller = label_both) +
  labs(title = "Sinuosity vs bout distance by condition",
       x = "Euclidean distance (m)", y = "Sinuosity", color = "Open cost") +
  theme_bw() +
  coord_cartesian(ylim = c(1, 2))
ggsave("bout_v3_sinuosity_vs_dist.png", p5, dpi = 300, width = 14, height = 5)

# Summary table
summary_table <- results_df %>%
  group_by(cost1, c1_fraction) %>%
  summarise(
    n = n(),
    med_sinuosity    = median(sinuosity, na.rm = TRUE),
    med_euclid       = median(euclidean, na.rm = TRUE),
    med_lcp          = median(lcp_length, na.rm = TRUE),
    med_cost_euclid  = median(cost_per_euclid, na.rm = TRUE),
    med_frac_open    = median(frac_dist_open, na.rm = TRUE),
    mean_frac_open   = mean(frac_dist_open, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_table, paste0(metadata, "_summary.csv"), row.names = FALSE)
print(summary_table, n = 50)

message("\nDone!")