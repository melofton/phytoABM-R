# Title: Workflow script for HABs ABM with one PFT
# Author: Mary Lofton
# Date: 16MAR24

# Purpose: Workflow to run R-based phytoplankton ABM

# =============================================================================
# Load packages and custom functions
# =============================================================================
library(tidyverse)
library(rLakeAnalyzer)
source("./code/HABs_ABM_functions_1group.R")

# =============================================================================
# Simulate phytoplankton dynamics
# =============================================================================

# ----- Initialise depths
lake_depths <- seq(from = 0.1, to = 9.3, by = 0.1)

# ----- Initialise individuals (phytos)
inds <- initialize_phytos(depths = lake_depths)

# ----- Initialise environment (water temperature)
lake_env <- initialize_env(depths = lake_depths)


# ---- Assign traits
## define traits for PFTs
traits_lst <- list(
  #create cell diameter trait
  diam = 20, #um
  
  #create cell density trait
  dens = 1050, # In range for green algae
  
  #create cell shape trait
  shape = 1, # make them all spherical for now
  
  #create temperature sensitivity traits
  Tmin = c(4),
  Topt = c(23),
  Tmax = c(30),
  I_S = c(300)
)


# ---- Start the simulation 
env <- lake_env$env_init
ts         <- 0;
time_steps <- 60*72;
inds_hist  <- NULL;
start_time <- Sys.time()
while(ts < time_steps){
  inds            <- movement(inds, env); 
  inds            <- growth(inds, repr_col = 7, traits = traits_lst, growth_env = env);
  inds            <- death(inds);
  ts              <- ts + 1; 
  env             <- update_env(env, lake_env$wtemp, lake_env$swr, time_steps = time_steps, tstep = ts, depths = lake_depths);
  inds_hist[[ts]] <- inds;
  print(ts)
  print(length(inds[,1]))
}
end_time <- Sys.time()
run_time = end_time - start_time
print(run_time)




# =============================================================================
# Print the results
# =============================================================================

ind_yloc <- array(data = NA, dim = c(time_steps,length(lake_depths)+1))
colnames(ind_yloc) <- c("timestep",lake_depths)

for(i in 1:time_steps){
  for(j in 1:length(lake_depths)){
    
    #generic PFT
    ind_yloc[i, 1] <- i;                      # Save the time step
    ind_yloc[i, j+1] <- length(inds_hist[[i]][which(round(inds_hist[[i]][,2],1) == round(lake_depths[j],1) & inds_hist[[i]][,1] == 1),1]); # Save the number of individuals at each depth
  }
}
#print(ind_yloc);




# =============================================================================
# Wrangle and save the results
# =============================================================================

##plot IBM output
plot_yloc <- data.frame(ind_yloc) %>%
  gather(X0.1:X9.3, key = "Depth_m",value = "num_agents") %>%
  mutate(Depth_m = as.double(substring(Depth_m, 2)))

plot_ts <- plot_yloc %>%
  filter(Depth_m >= 0.1 & Depth_m <= 2) %>%
  group_by(timestep) %>%
  summarize(surface_agents = mean(num_agents))

# #Write output to file for plotting for proposal
write.csv(plot_yloc, file = "./model_output/ABM_depthByTimestep_47hr.csv", row.names = FALSE)
write.csv(plot_ts, file = "./model_output/ABM_surfaceTimeseries_47hr.csv", row.names = FALSE)

saveRDS(inds_hist, file = "./model_output/ABM_output_47hr.rds")
