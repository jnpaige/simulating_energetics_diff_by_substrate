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
k<-5   # Simulations per variable combo
cost_tape0<-c(1) ## Keep this as a baseline this is the closed substrate
cost_tape1<-c(1.5,2.5,3.5,4.5) ## Vary the cost of second substrate (open)
target_open_fraction=c(.1,.2,.3,.4) ## 

closed_cluster_size=5
resource_cluster_size=5
metadata<-paste("targ", target_open_fraction, 
                "_cl_size", closed_cluster_size,
#                "_r_p", resource_prob,
                "_r_size", resource_cluster_size,collapse="",sep="")

set.seed(100)
results <- data.frame(
  total_length = numeric(),
  euclidean = numeric(),
  sinuosity = numeric(),
  total_cost=numeric(),
  cost1 = numeric(),
  cost2 = numeric(),
  replicate = integer()
)


#preallocate grid
param_grid <- expand.grid(cost0 = cost_tape0,
                          cost1 = cost_tape1,
                          c1_fraction=target_open_fraction,
                          replicate = seq_len(k))
param_grid
n_total <- nrow(param_grid)
param_grid$metadata<-metadata
# Example data frame to store results
param_grid$total_length <- NA
param_grid$euclidean <- NA
param_grid$sinuosity <- NA
param_grid$total_cost <- NA


for(i in 1:nrow(param_grid)){
    #output <- model_landscape_and_movement_v8(n, mod, param_grid$cost1[i],param_grid$cost2[i])
  output <- model_landscape_and_movement(grid_size=n, 
                                         mod=mod, 
                                         cost0=param_grid$cost0[i],
                                         cost1=param_grid$cost1[i],
                                         c1_fraction =param_grid$c1_fraction[i], 
                                         c1_cluster_size = closed_cluster_size,
                                         resource_cluster_size=resource_cluster_size)
    
  if (is.null(output)) {
    next  # skip this iteration
  }
  
    # calculate metrics
    metrics <- calculate_path_metrics(output[[1]], output[[2]])
    
    # store results
    param_grid$total_length[i] <- metrics$total_length
    param_grid$euclidean[i] <- metrics$euclidean
    param_grid$sinuosity[i] <- metrics$sinuosity
    param_grid$total_cost[i] <- metrics$total_cost
    
}
here::here()
param_grid 


plot(output[[2]][[1]])
plot(output[[1]][[1]])
plot(output[[2]][[2]])
file_name<-paste0(format(Sys.time(), "%Y_%m_%d_%H"),metadata,".csv",sep="",collapse="_")


write.csv(param_grid,file_name)

library(ggplot2)
library(tidybayes)
library(ggridges)
#Basic plots. 
df<-param_grid
length(df$cost0)

df <- df[df$total_length > 0, ]
df <- df[which(df$total_cost>0),]
df$rat <- df$total_cost / df$total_length



ggplot(df, aes(x=rat,y=as.factor(c1_fraction),fill=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1)


ggplot(df, aes(x=sinuosity,y=as.factor(c1_fraction),fill=as.factor(c1_fraction)))+geom_density_ridges()+
  facet_wrap(~cost1)


