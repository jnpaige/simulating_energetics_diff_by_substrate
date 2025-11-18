### Iterated Djikstra sim v5.
## This version has a trycatch in the simulation loop so that we skip over
## the common edge cases where a lcp fails. Haven't diagnosed why it fails yet
## But its in 1 in 20 or so simulations, so highly frequent. 
## I've been developing different versions of the modelling function
## To try to remove edge effects (doesn't seem to be the cause v5 function)
## It might be because the xy coordinate is too close to the resource. 
## But, I tried but still regular failures (v6)
## For now I keep skipping them. 
## Ok figures it out. Ocassionally there are substrates with no variability, its all either closed or its open.
## This causes problems in calculating the lcp. So in the sim we ensure that there are 
## some proportions of open/closed areas. This is hacky, and future versions should
## Produce a more elegant way of simulating landscape. 
## The v8 of the model function works well. cleaning it up though to simplify. 

# install.packages(c("terra", "gdistance"))
library(terra)
library(gdistance)
library(here)

setwd(paste(here::here(),"/Functions",sep="",collapse=""))
source("Iterated_djikstra_approach_functions_V3.R")
#set.seed(42)

### 1. Create base landscape
n <- 100  # grid size
mod <- 5 # number of cells to remove from each side
#cost1<-1
#cost2<-3
k<-40   # Simulations per variable combo
cost_tape1<-c(1)
cost_tape2<-c(1,2,4,8)
target_closed_fraction=.8
closed_cluster_size=5
resource_cluster_size=5
metadata<-paste("targ", target_closed_fraction, 
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
param_grid <- expand.grid(cost1 = cost_tape1,
                          cost2 = cost_tape2,
                          replicate = seq_len(k))
n_total <- nrow(param_grid)
param_grid$metadata<-metadata

# Example data frame to store results
param_grid$total_length <- NA
param_grid$euclidean <- NA
param_grid$sinuosity <- NA
param_grid$total_cost <- NA

for(i in 1:nrow(param_grid)){
    #output <- model_landscape_and_movement_v8(n, mod, param_grid$cost1[i],param_grid$cost2[i])
  output <- model_landscape_and_movement(n, mod, param_grid$cost1[i],param_grid$cost2[i],
                                            target_closed_fraction=target_closed_fraction, 
                                            closed_cluster_size = closed_cluster_size,
                                            resource_cluster_size=resource_cluster_size)
    
  if (is.null(output)) {
    next  # skip this iteration
  }
  
    # calculate metrics
    metrics <- calculate_path_metrics(output[[2]][[1]], output[[2]][[2]])
    
    # store results
    param_grid$total_length[i] <- metrics$total_length
    param_grid$euclidean[i] <- metrics$euclidean
    param_grid$sinuosity[i] <- metrics$sinuosity
    param_grid$total_cost[i] <- metrics$total_cost
    
  }
param_grid
i
plot(output[[2]][[1]])
plot(output[[1]][[1]])
plot(output[[2]][[2]])
file_name<-paste0(format(Sys.time(), "%Y_%m_%d_%H"),metadata,".csv",sep="",collapse="_")


write.csv(param_grid,file_name)


#Basic plots. 
df<-param_grid
ggplot(df,aes(x=total_length, y=sinuosity))+geom_point() +facet_wrap(~cost2) 
ggplot(df,aes(x=total_length, y=total_cost))+geom_point() +facet_wrap(~cost2) 

  #xlab("difference in cost between substrates")


ggplot(df,aes(x=total_length, group=cost2, fill=cost2))+geom_density(alpha=.2) 
ggplot(df,aes(x=euclidean, group=cost2, fill=cost2))+geom_density(alpha=.2) 
ggplot(df,aes(x=total_length, group=cost2, fill=cost2))+geom_density(alpha=.2) 
ggplot(df,aes(x=sinuosity, group=cost2, fill=cost2))+geom_density(alpha=.2) 
ggplot(df,aes(x=total_cost, group=cost2, fill=cost2))+geom_density(alpha=.2) 
ggplot(df,aes(x=euclidean,y=total_length, group=cost2, color=cost2))+geom_point(alpha=.9) 


ggplot(df,aes(x=total_cost/total_length, group=cost2, fill=cost2))+geom_density(alpha=.2) 








