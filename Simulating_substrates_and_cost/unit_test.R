library(terra)
library(ggplot2)
library(sf)
library(dplyr)

setwd(paste(here::here(),"/Simulating_substrates_and_cost/functions",sep="",collapse=""))
source("functions_V3.R")
#set.seed(42)

### 1. Create base landscape
n <- 100  # grid size
mod <- 5 # number of cells to remove from each side
cost1<-1
cost2<-2

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





