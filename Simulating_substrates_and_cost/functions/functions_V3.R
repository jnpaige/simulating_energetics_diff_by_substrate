library(terra)
library(gdistance)
library(raster)
## Functions
#Comprehensive function that generats rasters and also builds cost surf and least path.



#In this version, only one resource is planted on the landscape. This avoids a lot of stochasticity
##we are not interested in testing. Of course, on a super resource dense landscape you don't need to move much.
##So locomotor type efficiency not that big a deal. 

model_landscape_and_movement <- function(grid_size, mod, cost0, cost1,
                                         c1_fraction,
                                         c1_cluster_size,
                                         resource_cluster_size) {
  
  
  # 1. Base raster
  land <- terra::rast(nrows = grid_size, ncols = grid_size,
                      xmin = 0, xmax = grid_size,
                      ymin = 0, ymax = grid_size)
  

  ### 2. Improved substrate raster generation with controlled fraction
  substrate <- terra::rast(land)
  terra::values(substrate) <- runif(terra::ncell(substrate))
  
  # Create clusters
  substrate <- terra::focal(substrate, w = c1_cluster_size, fun = mean, na.rm = TRUE)

  # Normalize
  substrate_min <- terra::global(substrate, "min", na.rm = TRUE)[[1]]
  substrate_max <- terra::global(substrate, "max", na.rm = TRUE)[[1]]
  #normalize the clusters
  
  substrate <- (substrate - substrate_min) / (substrate_max - substrate_min)

  # Pull values to memory safely
  v <- terra::values(substrate)
  v[is.na(v)] <- 0
  
  
  # Quantile-based threshold
  # If you less than or qual the threshold you get a 1 (this is open area/cost2), otherwise 0 (its closed/cost1)
  threshold <- stats::quantile(v, probs = c1_fraction, na.rm = TRUE)
  v_bin <- ifelse(v >= threshold, 0, 1)
  
  
  # Assign back
  substrate_binary <- substrate
  terra::values(substrate_binary) <- v_bin
  
  substrate <- substrate_binary

  ### 3. Improved resource raster generation
  resource <- terra::rast(land)
  # Ensure at least one resource
  random_cell <- sample(terra::ncell(resource), 1)
  resource[random_cell] <- 1
  
  resource <- terra::focal(resource, w = resource_cluster_size, fun = mean, na.rm = TRUE)
  
  rvals <- terra::values(resource)
  rvals[is.na(rvals)] <- 0
  
  if (any(rvals > 0, na.rm = TRUE)) {
    resource_threshold <- stats::quantile(rvals[rvals > 0], 0.7, na.rm = TRUE)
    rvals <- ifelse(rvals < resource_threshold, 0, 1)
  } else {
    rvals[] <- 0
  }
  
  terra::values(resource) <- rvals
  
  
  ### 4. Travel cost raster (robust version)
  
  substrate <- terra::ifel(substrate >= 1, 1, 0)
  substrate[is.na(substrate)] <- 0
  
  
  
  # Force raster to materialize in memory
  #substrate <- terra::deepcopy(substrate)
  #substrate <- terra::setValues(substrate, terra::values(substrate))
  

  # If you less than or qual the threshold you get a 1 (this is open area/cost2), otherwise 0 (its closed/cost1)
  # Now safely build travel cost layer
  travel_cost <- substrate
  vals <- terra::values(substrate)
  vals[vals == 1] <- cost1
  vals[vals == 0] <- cost0
  travel_cost <- terra::setValues(travel_cost, vals)

  
  ### 5. Crop extent
  crop_ext <- terra::ext(mod, grid_size - mod, mod, grid_size - mod)
  substrate_crop <- terra::crop(substrate, crop_ext)
  resource_crop <- terra::crop(resource, crop_ext)
  travel_cost_crop <- terra::crop(travel_cost, crop_ext)
  
  ### 6. Convert to raster for gdistance
  travel_cost_raster <- raster::raster(travel_cost_crop)
  resource_raster <- raster::raster(resource_crop)
  
  
  ### 7. Build gdistance transition object
  tr <- gdistance::transition(1 / travel_cost_raster, transitionFunction = mean, directions = 8)
  trC <- gdistance::geoCorrection(tr, type = "c")
  
  ### 8. Cost surface to nearest resource
  resource_coords <- raster::xyFromCell(resource_raster,
                                        which(raster::values(resource_raster) == 1))
  if (nrow(resource_coords) == 0) {
    warning("No resource cells found, skipping iteration.")
    return(NULL)
  }
  
  cost_to_resources <- tryCatch({
    gdistance::accCost(trC, resource_coords)
  }, error = function(e) return(NULL))
  if (is.null(cost_to_resources)) {
    warning("accCost failed for resource_coords, skipping iteration.")
    return(NULL)
  }
  
  cost_surface_terra <- terra::rast(cost_to_resources)
  ### 9. Pick start point far from resources
  non_resource_cells <- which(raster::values(resource_raster) == 0)
  if (length(non_resource_cells) == 0) {
    warning("No non-resource cells found, skipping iteration.")
    return(NULL)
  }
  
  sample_cells <- sample(non_resource_cells, min(100, length(non_resource_cells)))
  sample_coords <- raster::xyFromCell(resource_raster, sample_cells)
  
  distances <- sapply(seq_len(nrow(sample_coords)), function(i) {
    xy <- sample_coords[i, , drop = FALSE]
    cost_from_xy <- tryCatch(gdistance::accCost(trC, xy), error = function(e) return(NULL))
    if (is.null(cost_from_xy)) return(NA)
    min(raster::values(cost_from_xy)[raster::values(resource_raster) == 1], na.rm = TRUE)
  })
  
  if (all(is.na(distances))) {
    warning("All distances NA — disconnected landscape. Skipping iteration.")
    return(NULL)
  }
  
  distances[!is.finite(distances)] <- NA
  start_cell <- sample_cells[which.max(distances)]
  xy_start <- raster::xyFromCell(resource_raster, start_cell)
  
  ### 10. Least-cost path
  cost_from_start <- tryCatch({
    gdistance::accCost(trC, xy_start)
  }, error = function(e) return(NULL))
  if (is.null(cost_from_start)) {
    warning("Failed to compute accCost from start, skipping iteration.")
    return(NULL)
  }
  
  resource_cost_values <- terra::extract(terra::rast(cost_from_start), resource_coords)
  if (all(is.na(resource_cost_values[,1]))) {
    warning("No reachable resource from start, skipping iteration.")
    return(NULL)
  }
  
  nearest_resource_idx <- which.min(resource_cost_values[,1])
  xy_nearest_resource <- resource_coords[nearest_resource_idx, , drop = FALSE]
  
  lcp <- tryCatch({
    gdistance::shortestPath(trC, xy_start, xy_nearest_resource, output = "SpatialLines")
  }, error = function(e) return(NULL))
  if (is.null(lcp)) {
    warning("shortestPath failed, skipping iteration.")
    return(NULL)
  }
  
  lcp_terra <- terra::vect(lcp)
  total_cum_cost <- terra::extract(terra::rast(cost_from_start),
                                   terra::vect(matrix(xy_nearest_resource, ncol = 2))) |>
    as.numeric()
  

  return(list(cost_surface_terra,lcp_terra, total_cum_cost, xy_start, xy_nearest_resource))
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


