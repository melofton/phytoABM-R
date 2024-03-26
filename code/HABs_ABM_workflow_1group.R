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
lake_env <- initialize_env(depths = lake_depths, n_days = 14)


# ---- Assign traits
## define traits for PFTs
traits_lst <- list(
  #create cell diameter trait
  diam = 20, #um; not sure how to do this for filaments
  
  #create cell density trait
  dens = 950, # In range for cyanobacteria
  
  #create cell shape trait
  shape = 1, # this represents spherical and we need filament eventually
  
  #create temperature and light sensitivity traits
  # Tmin = c(4),
  # Topt = c(23),
  # Tmax = c(30),
  T_0 = c(22), # optimal temperature for Anabaena Prokopkin et al. https://doi.org/10.1016/j.ecolmodel.2005.05.011
  q = c(5), # thermal dispersion for Anabaena Prokopkin et al. https://doi.org/10.1016/j.ecolmodel.2005.05.011
  I_S = c(300),
  
  #create maximum growth rate trait
  umax = 0.78/1440, #maximum specific growth rate for Anabaena reported in Reynolds 2006
  
  #respiration traits
  R_resp = 0.08/1440, #this is from Cayelan/Kamilla's GLM-AED calibration converted to minute scale
  theta_resp = 1.08,
  
  #nutrient uptake traits
  K_N = 2,
  N_0 = 0,
  K_P = 0.05,
  P_0 = 0
)


# ---- Start the simulation 
env <- lake_env$env_init
ts         <- 0;
time_steps <- 60*330;
inds_hist  <- NULL;
start_time <- Sys.time()
while(ts < time_steps){
  inds            <- movement(inds, env, yloc = 2, ymax = 9.3, wnd = unname(env[1,6]), delta_z = 9); 
  inds            <- growth(inds, repr_col = 7, traits = traits_lst, growth_env = env);
  inds            <- death(inds, traits = traits_lst);
  ts              <- ts + 1; 
  env             <- update_env(env, lake_env$wtemp, lake_env$met, lake_env$din, lake_env$frp, time_steps = time_steps, tstep = ts, depths = lake_depths);
  inds_hist[[ts]] <- inds;
  print(ts)
  print(length(inds[,1]))
}
end_time <- Sys.time()
run_time = end_time - start_time
print(run_time)




# =============================================================================
# Wrangle and save the results
# =============================================================================

# load in previous model run if desired
inds_hist <- readRDS("./model_output/ABM_output_330hr.rds")

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

##dataframes for plotting
plot_yloc <- data.frame(ind_yloc) %>%
  gather(X0.1:X9.3, key = "Depth_m",value = "num_agents") %>%
  mutate(Depth_m = as.double(substring(Depth_m, 2)))

plot_ts <- plot_yloc %>%
  group_by(timestep) %>%
  summarize(surface_agents = sum(num_agents))

# #Write output to file for plotting 
write.csv(plot_yloc, file = "./model_output/ABM_depthByTimestep_330hr.csv", row.names = FALSE)
write.csv(plot_ts, file = "./model_output/ABM_agentTimeseries_330hr.csv", row.names = FALSE)

saveRDS(inds_hist, file = "./model_output/ABM_output_330hr.rds")
