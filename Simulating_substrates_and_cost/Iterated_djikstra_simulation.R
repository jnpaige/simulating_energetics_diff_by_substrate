## Iterated Djikstra sim.

##

# install.packages(c("terra", "gdistance"))
library(terra)
library(gdistance)
library(here)


setwd(paste(here::here(),"/Simulating_substrates_and_cost/functions",sep="",collapse=""))


source("Iterated_djikstra_approach_functions_V3.R")
#set.seed(42)

### 1. Create base landscape
n <- 100  # grid size
mod <- 5 # number of cells to remove from each side
#cost1<-1
#cost2<-3
k<-20   # Simulations per variable combo
cost_tape1<-c(1) ## Keep this as a baseline
cost_tape2<-c(1.5,2,4) ## Vary the cost of second substrate (open)
target_closed_fraction=c(.6,.7,.8,.9) ## Vary the proportion of first substrate (closed).

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
                          target_closed_fraction=target_closed_fraction,
                          replicate = seq_len(k))
n_total <- nrow(param_grid)
param_grid$metadata<-metadata
# Example data frame to store results
param_grid$total_length <- NA
param_grid$euclidean <- NA
param_grid$sinuosity <- NA
param_grid$total_cost <- NA
param_grid

for(i in 1:nrow(param_grid)){
    #output <- model_landscape_and_movement_v8(n, mod, param_grid$cost1[i],param_grid$cost2[i])
  output <- model_landscape_and_movement(n, mod, param_grid$cost1[i],param_grid$cost2[i],
                                            target_closed_fraction=param_grid$target_closed_fraction[i], 
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
here::here()
param_grid
i
plot(output[[2]][[1]])
plot(output[[1]][[1]])
plot(output[[2]][[2]])
file_name<-paste0(format(Sys.time(), "%Y_%m_%d_%H"),metadata,".csv",sep="",collapse="_")


write.csv(param_grid,file_name)

library(ggplot2)
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








