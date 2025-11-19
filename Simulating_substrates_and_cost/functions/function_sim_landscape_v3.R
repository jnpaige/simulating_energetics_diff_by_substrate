library(terra)
library(gdistance)
library(raster)


model_landscape_and_movement <- function(grid_size, mod, cost0, cost1,
                                         c1_fraction,
                                         c1_cluster_size,
                                         resource_cluster_size) {
  
  
  
  # 1. Base raster
  land <- terra::rast(nrows = grid_size, ncols = grid_size,
                      xmin = 0, xmax = grid_size,
                      ymin = 0, ymax = grid_size)
  
  plot(land)
  l1 <- list()
  
  ### 2. Improved substrate raster generation with controlled fraction
  substrate <- terra::rast(land)
  terra::values(substrate) <- runif(terra::ncell(substrate))
  
  # Create clusters
  substrate <- terra::focal(substrate, w = c1_cluster_size, fun = mean, na.rm = TRUE)
  plot(substrate)
  
  # Normalize
  substrate_min <- terra::global(substrate, "min", na.rm = TRUE)[[1]]
  substrate_max <- terra::global(substrate, "max", na.rm = TRUE)[[1]]
  #normalize the clusters
  
  substrate <- (substrate - substrate_min) / (substrate_max - substrate_min)
  plot(substrate)
  
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
  plot(substrate)
  
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
  
  plot(substrate)
  
  # Force raster to materialize in memory
  #substrate <- terra::deepcopy(substrate)
  #substrate <- terra::setValues(substrate, terra::values(substrate))
  
  plot(substrate)
  
  # If you less than or qual the threshold you get a 1 (this is open area/cost2), otherwise 0 (its closed/cost1)
  # Now safely build travel cost layer
  travel_cost <- substrate
  plot(travel_cost)
  vals <- terra::values(substrate)
  vals[vals == 1] <- cost1
  vals[vals == 0] <- cost0
  travel_cost <- terra::setValues(travel_cost, vals)
  plot(travel_cost)
  
  
  ### 5. Crop extent
  crop_ext <- terra::ext(mod, grid_size - mod, mod, grid_size - mod)
  substrate_crop <- terra::crop(substrate, crop_ext)
  resource_crop <- terra::crop(resource, crop_ext)
  travel_cost_crop <- terra::crop(travel_cost, crop_ext)
  
  ### 6. Convert to raster for gdistance
  travel_cost_raster <- raster::raster(travel_cost_crop)
  resource_raster <- raster::raster(resource_crop)
  
  plot(travel_cost_raster)
  l1[[1]] <- travel_cost_raster
  l1[[2]] <- resource_raster
  
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
  plot(cost_to_resources)
  plot(cost_surface_terra)
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
  
  
  l2 <- list(cost_surface_terra, lcp_terra, total_cum_cost, xy_start, xy_nearest_resource)
  final_list <- list(l1, l2)
  return(final_list)
}





##testing


### 1. Create base landscape
grid_size <- 100  # grid size
mod <- 5 # number of cells to remove from each side
cost0<-1
cost1<-2

c1_fraction=.05
c1_cluster_size=5
resource_prob=0.01
resource_cluster_size=5

  
































t<-model_landscape_and_movement(grid_size= n,
                                mod=mod,
                                cost1=cost1,
                                cost2=cost2,  
                                c1_fraction=c1_fraction,
                                c1_cluster_size=c1_cluster_size, 
                                resource_cluster_size=resource_cluster_size )



# Extract components
substrate <- t[[1]][[1]]        # raster (may be RasterLayer)
resources <- t[[1]][[2]]        # raster (may be RasterLayer)
cost_surface <- t[[2]][[1]]     # SpatRaster
lcp <- t[[2]][[2]]              # SpatVector
xy_start <- t[[2]][[4]]
xy_target <- t[[2]][[5]]

t<-prepare_for_plotting(substrate, cost_surface,lcp,resources,xy_start,xy_target)

### Now plot data. 
df_substrate<-t[[1]]
lcp_sf<-t[[2]]
res_sf<-t[[3]]
start_sf<-t[[4]]
target_sf<-t[[5]]
df_cost<-t[[6]]



fig1 <- ggplot() +
  geom_raster(data = df_substrate, aes(x, y, fill = factor(substrate))) +
  scale_fill_manual(values = c("gray20", "gray80"),
                    name = "Substrate",
                    labels = c("Closed (0)", "Open (1)")) +
  geom_sf(data = lcp_sf, color = "red", size = 1.1) +
  geom_sf(data = res_sf, color = "white", size = 1.2, shape = 21, fill = "white") +
  geom_sf(data = start_sf, color = "black", fill = "blue", shape = 21, size = 3) +
  geom_sf(data = target_sf, color = "black", fill = "yellow", shape = 21, size = 3) +
  labs(title = "(A) Substrate Map with Resources and Least-Cost Path") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"))

fig1

fig2 <- ggplot() +
  geom_raster(data = df_cost, aes(x, y, fill = cost)) +
  scale_fill_viridis_c(option = "C", name = "Cost to nearest resource") +
  geom_sf(data = lcp_sf, color = "red", size = 1.1) +
  geom_sf(data = res_sf, color = "white", size = 1.2, shape = 21, fill = "white") +
  geom_sf(data = start_sf, color = "black", fill = "blue", shape = 21, size = 3) +
  geom_sf(data = target_sf, color = "black", fill = "yellow", shape = 21, size = 3) +
  labs(title = "(B) Cost Surface with Resources and Least-Cost Path") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"))

fig2
fig1







