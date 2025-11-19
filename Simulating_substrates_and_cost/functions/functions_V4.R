library(terra)
library(gdistance)

model_landscape_and_movement <- function(grid_size, cost_open, cost_closed, open_fraction) {
  # Create landscape
  land <- rast(nrows = grid_size, ncols = grid_size, xmin = 0, xmax = grid_size, ymin = 0, ymax = grid_size)
  
  # Simple substrate with random clusters
  substrate <- rast(land)
  values(substrate) <- runif(ncell(substrate))
  substrate <- focal(substrate, w = 5, fun = mean)
  substrate <- substrate > quantile(values(substrate), open_fraction, na.rm = TRUE)
  
  # Simple resource placement
  resource <- rast(land)
  resource[sample(ncell(resource), 10)] <- 1
  resource <- focal(resource, w = 3, fun = mean) > 0.3
  
  # Cost surface
  cost_surface <- ifel(substrate, cost_open, cost_closed)
  
  # Find path
  cost_raster <- raster(cost_surface)
  tr <- transition(1 / cost_raster, mean, 8)
  tr <- geoCorrection(tr, "c")
  
  resource_coords <- xyFromCell(raster(resource), which(values(raster(resource)) == 1))
  if (nrow(resource_coords) == 0) return(NULL)
  
  # Start from edge
  start_point <- c(1, 1)
  costs <- accCost(tr, start_point)
  target_idx <- which.min(extract(rast(costs), resource_coords)[,1])
  target_point <- resource_coords[target_idx, ]
  
  path <- shortestPath(tr, start_point, target_point, "SpatialLines")
  
  return(list(
    cost_surface = rast(costs),
    path = vect(path),
    start = start_point,
    target = target_point
  ))
}

calculate_path_metrics <- function(result) {
  if (is.null(result)) return(NULL)
  
  coords <- geom(result$path)[, c("x", "y")]
  if (nrow(coords) < 2) return(NULL)
  
  segments <- sqrt(diff(coords[,1])^2 + diff(coords[,2])^2)
  path_length <- sum(segments)
  euclidean <- sqrt(sum((coords[1,] - coords[nrow(coords),])^2))
  
  path_cost <- sum(extract(result$cost_surface, result$path)$layer, na.rm = TRUE)
  
  data.frame(
    path_length = path_length,
    euclidean = euclidean,
    sinuosity = path_length / euclidean,
    total_cost = path_cost
  )
}
