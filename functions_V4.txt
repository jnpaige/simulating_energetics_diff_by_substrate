# functions_fast.R
library(terra)
library(raster)
library(gdistance)
library(sf)

# This function returns a list like v3.
# list(cost_surface_terra, lcp_terra, total_cum_cost, xy_start, xy_nearest_resource)
# but it tries to avoid expensive conversions and repeated accCost calls.

model_landscape_and_movement_fast <- function(grid_size, mod, cost0, cost1,
                                              c1_fraction,
                                              c1_cluster_size,
                                              resource_cluster_size,
                                              raster_template = NULL) {
  #defensive checks
  if (!is.numeric(c1_fraction) || length(c1_fraction) != 1 || is.na(c1_fraction)) {
    stop("c1_fraction must be a single numeric in (0,1).")
  }
  if (missing(raster_template) || is.null(raster_template)) {
    stop("raster_template (a raster::RasterLayer for the CROPPED domain) must be provided.")
  }
  
  #1. Create base terra raster (full grid)
  land <- terra::rast(nrows = grid_size, ncols = grid_size,
                      xmin = 0, xmax = grid_size,
                      ymin = 0, ymax = grid_size)
  
  #2. Substrate generation (vectorized, minimal terra ops)
  substrate <- terra::rast(land)
  terra::values(substrate) <- runif(terra::ncell(substrate))
  substrate <- terra::focal(substrate, w = c1_cluster_size, fun = mean, na.rm = TRUE)
  
  # normalize 0..1
  smn <- terra::global(substrate, "min", na.rm = TRUE)[[1]]
  smx <- terra::global(substrate, "max", na.rm = TRUE)[[1]]
  if (smx == smn) {
    terra::values(substrate) <- 0
  } else {
    substrate <- (substrate - smn) / (smx - smn)
  }
  
  v <- terra::values(substrate)
  v[is.na(v)] <- 0
  
  # threshold -> produce binary substrate (1 = open, 0 = closed)
  threshold <- stats::quantile(v, probs = c1_fraction, na.rm = TRUE)
  
  # FIXED: 1=open should be the *lowest* c1_fraction of terrain,
  # not the highest, unless you explicitly want high==open.
  v_bin <- ifelse(v <= threshold, 1L, 0L)
  
  substrate_binary <- substrate
  terra::values(substrate_binary) <- v_bin
  substrate <- substrate_binary
  
  #3. Resource raster (single seeded cluster) ----------
  resource <- terra::rast(land)
  random_cell <- sample(terra::ncell(resource), 1)
  terra::values(resource) <- 0
  terra::values(resource)[random_cell] <- 1
  resource <- terra::focal(resource, w = resource_cluster_size, fun = mean, na.rm = TRUE)
  
  rvals <- terra::values(resource)
  rvals[is.na(rvals)] <- 0
  if (any(rvals > 0, na.rm = TRUE)) {
    rth <- stats::quantile(rvals[rvals > 0], 0.7, na.rm = TRUE)
    rvals <- ifelse(rvals < rth, 0L, 1L)
  } else {
    rvals[] <- 0L
  }
  terra::values(resource) <- rvals
  
  # ---------- 4. Build travel cost raster for full grid, then crop to template ----------
  # substrate is 0/1 (closed/open). Map 1 -> cost1, 0 -> cost0
  vals_sub <- terra::values(substrate)
  vals_sub[is.na(vals_sub)] <- 0
  vals_cost_full <- ifelse(vals_sub >= 1, cost1, cost0)
  
  travel_cost_full <- substrate
  terra::values(travel_cost_full) <- vals_cost_full
  
  # Crop to central area (matching raster_template extent)
  crop_ext <- terra::ext(mod, grid_size - mod, mod, grid_size - mod)
  travel_cost_crop <- terra::crop(travel_cost_full, crop_ext)
  resource_crop <- terra::crop(resource, crop_ext)
  substrate_crop <- terra::crop(substrate, crop_ext)
  
  # ---------- 5. Convert to raster::RasterLayer using provided template ----------
  # raster_template must match the cropped domain geometry (nrows/ncols/xmin/xmax/...)
  # We copy values into that template using raster::setValues (fast)
  vals_tc <- terra::values(travel_cost_crop)
  vals_res <- terra::values(resource_crop)
  vals_sub_crop <- terra::values(substrate_crop)
  
  # sanity checks
  if (length(vals_tc) != raster::ncell(raster_template)) {
    stop("raster_template cells do not match cropped domain -- check mod / grid_size.")
  }
  
  travel_cost_r <- raster::setValues(raster_template, vals_tc)
  resource_r <- raster::setValues(raster_template, vals_res)
  substrate_r <- raster::setValues(raster_template, vals_sub_crop)
  
  # Guard against non-positive costs
  if (any(!is.finite(raster::values(travel_cost_r)))) stop("Non-finite travel costs.")
  if (any(raster::values(travel_cost_r) <= 0)) stop("travel_cost_r must be > 0 everywhere.")
  
  # ---------- 6. Build gdistance transition & geo-correct ----------
  tr <- gdistance::transition(1 / travel_cost_r, transitionFunction = mean, directions = 8)
  trC <- gdistance::geoCorrection(tr, type = "c")
  
  # ---------- 7. cost surface to resources (single accCost call) ----------
  resource_cells_idx <- which(raster::values(resource_r) == 1)
  if (length(resource_cells_idx) == 0) {
    warning("No resource cells found in cropped domain; skipping.")
    return(NULL)
  }
  resource_coords <- raster::xyFromCell(resource_r, resource_cells_idx)
  
  cost_to_resources <- tryCatch({
    gdistance::accCost(trC, resource_coords)
  }, error = function(e) {
    warning("accCost failed: ", e$message)
    return(NULL)
  })
  if (is.null(cost_to_resources)) return(NULL)
  
  # ---------- 8. Choose start: farthest non-resource cell (use cost_to_resources values) ----------
  dist_vals <- raster::values(cost_to_resources)
  dist_vals[resource_cells_idx] <- NA  # ignore resource cells for choosing start
  if (all(is.na(dist_vals))) {
    warning("All distances NA (maybe disconnected) — skipping.")
    return(NULL)
  }
  start_idx <- which.max(dist_vals)
  xy_start <- raster::xyFromCell(cost_to_resources, start_idx)
  
  # ---------- 9. Nearest reachable resource to start (use cost surface values at resource coords) ----------
  # resource_cost_values: cost_to_resources at resource_coords
  resource_cost_values <- raster::extract(cost_to_resources, resource_coords)
  if (all(is.na(resource_cost_values))) {
    warning("No reachable resources from start.")
    return(NULL)
  }
  nearest_resource_idx <- which.min(resource_cost_values)
  xy_nearest_resource <- resource_coords[nearest_resource_idx, , drop = FALSE]
  
  # ---------- 10. Shortest path from start to nearest resource ----------
  lcp <- tryCatch({
    gdistance::shortestPath(trC, xy_start, xy_nearest_resource, output = "SpatialLines")
  }, error = function(e) {
    warning("shortestPath failed: ", e$message)
    return(NULL)
  })
  if (is.null(lcp)) return(NULL)
  
  # Convert SpatialLines to terra vect for compatibility with your downstream code
  # shortestPath returns SpatialLines; convert via sf then to terra vect
  lcp_sf <- sf::st_as_sfc(sf::st_as_text(sf::st_as_sfc(lcp)))
  # simpler conversion: use sp -> sf -> terra
  lcp_sp <- lcp
  lcp_sf <- sf::st_as_sf(lcp_sp)
  lcp_terra <- terra::vect(lcp_sp)
  
  # ---------- 11. Total cumulative cost along path (vectorized) ----------
  # extract coordinates of lcp vertices
  geom_mat <- do.call(rbind, lapply(lcp_sp@lines, function(L) {
    do.call(rbind, lapply(L@Lines, function(ll) ll@coords))
  }))
  # convert to SpatialPoints (raster::extract will take coords)
  pts <- sp::SpatialPoints(coords = geom_mat, proj4string = raster::crs(travel_cost_r))
  values_on_path <- raster::extract(cost_to_resources, pts)
  total_cum_cost <- sum(values_on_path, na.rm = TRUE)
  
  # Return objects in familiar structure
  cost_surface_terra <- terra::rast(cost_to_resources)  # small conversion once
  return(list(cost_surface_terra, lcp_terra, total_cum_cost, xy_start, xy_nearest_resource))
}






prepare_for_plotting<-function(substrate, cost_surface,lcp,resources,xy_start,xy_target){
  l<-list()
  # Convert rasters for ggplot
  df_substrate <- as.data.frame(rast(substrate), xy = TRUE)
  colnames(df_substrate)[3] <- "substrate"
  
  df_cost <- as.data.frame(cost_surface, xy = TRUE)
  colnames(df_cost)[3] <- "cost"
  
  # Convert SpatVector (least-cost path) to sf
  lcp_sf <- st_as_sf(lcp)
  
  # ✅ Convert resources to SpatRaster and then extract points
  resources_spat <- rast(resources)
  res_points <- as.points(resources_spat, values = TRUE)
  res_points <- res_points[res_points$lyr.1 == 1, ]
  res_sf <- st_as_sf(res_points)
  
  # Start/target points as sf
  start_sf <- st_as_sf(data.frame(x = xy_start[1], y = xy_start[2]), coords = c("x", "y"), crs = crs(resources_spat))
  target_sf <- st_as_sf(data.frame(x = xy_target[1], y = xy_target[2]), coords = c("x", "y"), crs = crs(resources_spat))
  
  l[[1]]<-df_substrate
  l[[2]]<-lcp_sf
  l[[3]]<-res_sf
  l[[4]]<-start_sf
  l[[5]]<-target_sf
  l[[6]]<-df_cost
  return(l)}






### Work with results:

library(terra)


calculate_path_metrics <- function(cost_surface,lcp) {
  g <- geom(lcp)  # matrix: geom, part, x, y, hole
  g <- as.data.frame(g)
  g$x <- as.numeric(g$x)
  g$y <- as.numeric(g$y)
  g$part <- as.numeric(g$part)
  
  total_length <- 0
  straight_start_end <- c()
  parts <- unique(g$part)
  
  for (p in parts) {
    rows <- which(g$part == p)
    if (length(rows) < 2) next
    coords <- cbind(g$x[rows], g$y[rows])
    dx <- diff(coords[,1])
    dy <- diff(coords[,2])
    segs <- sqrt(dx^2 + dy^2)
    total_length <- total_length + sum(segs, na.rm = TRUE)
    
    straight_start_end <- rbind(straight_start_end,
                                c(x1 = coords[1,1], y1 = coords[1,2],
                                  x2 = coords[nrow(coords),1], y2 = coords[nrow(coords),2]))
  }
  
  # Straight-line distance (start–end)
  if (nrow(straight_start_end) >= 1) {
    start <- straight_start_end[1, c("x1", "y1")]
    last  <- straight_start_end[nrow(straight_start_end), c("x2", "y2")]
    euclid <- sqrt((last[1] - start[1])^2 + (last[2] - start[2])^2)
  } else {
    euclid <- NA_real_
  }
  
  sinuosity <- ifelse(!is.na(euclid) && euclid > 0, total_length / euclid, NA_real_)
  
  total_cum_cost<-sum(terra::extract(cost_surface,lcp, touches=TRUE,cells=TRUE,xy=TRUE,method="bilinear")$layer)
  
  data.frame(
    total_length = total_length,
    euclidean = euclid,
    sinuosity = sinuosity,
    total_cost = total_cum_cost
  )
}

# Usage:
#metrics <- calculate_path_metrics(t2[[2]],t2[[3]][2])  # or lcp



