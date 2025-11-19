## Iterated Djikstra sim.

##

# install.packages(c("terra", "gdistance"))
library(terra)
library(gdistance)
library(here)


setwd(paste(here::here(),"/Simulating_substrates_and_cost/functions",sep="",collapse=""))


source("functions_V3.R")
#set.seed(42)

### 1. Create base landscape
n <- 100  # grid size
mod <- 5 # number of cells to remove from each side
#cost1<-1
#cost2<-3
cost0=1 ## Keep this as a baseline this is the closed substrate
cost1=1.5## Vary the cost of second substrate (open)
c1_fraction=.2 

c1_cluster_size=5
resource_cluster_size=5



output <- model_landscape_and_movement(grid_size=n, 
                                       mod=mod, 
                                       cost0=cost0,
                                       cost1=cost1,
                                       c1_fraction =c1_fraction[i], 
                                       c1_cluster_size = c1_cluster_size,
                                       resource_cluster_size=resource_cluster_size)

calculate_path_metrics(output[[1]],output[[2]])
prepare_for_plotting()
