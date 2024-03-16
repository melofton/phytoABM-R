# Title: Workflow script for HABs ABM
# Author: Mary Lofton
# Date: 16MAR24

# Purpose: Workflow to run R-based phytoplankton ABM

# =============================================================================
# Load packages and custom functions
# =============================================================================
library(tidyverse)
library(rLakeAnalyzer)
source("./code/HABs_ABM_functions.R")

# =============================================================================
# Simulate phytoplankton dynamics
# =============================================================================

# ----- Initialise individuals (phytos)

fp <- read_csv("./data/FluoroProbe_2021-08-02_FCR_50.csv")

# Create inds by looping through each depth increment and then creating that number of individuals at that depth based on ug/L
inds <- array(data = 0, dim = c(0, 8))

depths <- seq(from = 0.1, to = 9.3, by = 0.1)

for(i in 1:length(depths)){
  
  #isolate a particular depth
  temp <- fp[i,]
  
  #create cell diameter trait
  min_diam.d = 10 # 10 um 
  max_diam.d = 25 # 25 um in meters
  
  min_diam.c = 10 # 10 um 
  max_diam.c = 25 # 25 um in meters
  
  min_diam.g = 10 # 10 um 
  max_diam.g = 25 # 25 um in meters
  
  #create cell density trait
  min_dens.d = 1070 # In range for diatoms
  max_dens.d = 1130 # To make sure they all sink
  
  min_dens.c = 920 # In range for cyanos
  max_dens.c = 980 # To make sure they all float
  
  min_dens.g = 1020 # In range for chlorophytes
  max_dens.g = 1080 # To make sure they all sink
  
  #create cell shape trait
  shape.d = 1 # make them all spherical for now
  
  shape.c = 1 # make them all spherical for now
  
  shape.g = 1 # make them all spherical for now
  
  #populate depths with traits
  temp.df.d <- array(data = 0, dim = c(round(temp$BrownAlgae_ugL,0),8))
  temp.df.d[,1] <- 1 #placeholder for taxon ID or some other trait
  temp.df.d[,2] <- depths[i]
  temp.df.d[,3] <- runif(dim(temp.df.d)[1], min = min_diam.d, max = max_diam.d)
  temp.df.d[,4] <- runif(dim(temp.df.d)[1], min = min_dens.d, max = max_dens.d)
  temp.df.d[,5] <- shape.d
  
  temp.df.c <- array(data = 0, dim = c(round(temp$Bluegreens_ugL,0),8))
  temp.df.c[,1] <- 2 #placeholder for taxon ID or some other trait
  temp.df.c[,2] <- depths[i]
  temp.df.c[,3] <- runif(dim(temp.df.c)[1], min = min_diam.c, max = max_diam.c)
  temp.df.c[,4] <- runif(dim(temp.df.c)[1], min = min_dens.c, max = max_dens.c)
  temp.df.c[,5] <- shape.c
  
  temp.df.g <- array(data = 0, dim = c(round(temp$GreenAlgae_ugL,0),8))
  temp.df.g[,1] <- 3 #placeholder for taxon ID or some other trait
  temp.df.g[,2] <- depths[i]
  temp.df.g[,3] <- runif(dim(temp.df.g)[1], min = min_diam.g, max = max_diam.g)
  temp.df.g[,4] <- runif(dim(temp.df.g)[1], min = min_dens.g, max = max_dens.g)
  temp.df.g[,5] <- shape.g
  
  temp.df <- rbind(temp.df.d, temp.df.c, temp.df.g)
  
  inds <- rbind(inds, temp.df)
}

colnames(inds) <- c("PFG","yloc","cell_diam","cell_dens","cell_shape","dcol","repr","light_exp")





# ----- Initialise environment (water temperature)

env <- array(data = 0, dim = c(length(depths),5))

######### MANUALLY SIMULATED WATER TEMPERATURE ###########
# env[,1] <- 5 #water temperature;  no stratification
# env[,1] <- c(rev(seq(from = 15, to = 20, by = ((20-15)/(40-1)))),
#              rev(seq(from = 12, to = 15, by = ((15-12)/(10-1)))),
#              rev(seq(from = 10, to = 12, by = ((12-10)/(43-1))))) #water temperature; stronger stratification
# env[,1] <- c(rev(seq(from = 13, to = 15, by = ((15-13)/(40-1)))),
#              rev(seq(from = 11, to = 13, by = ((13-11)/(10-1)))),
#              rev(seq(from = 10, to = 11, by = ((11-10)/(43-1))))) #water temperature; weaker stratification
######### MANUALLY SIMULATED WATER TEMPERATURE ###########

######### GLM WATER TEMPERATURE ##########################
wtemp <- read_csv("./model_inputs/cal_wtemp_GLM.csv") 
env[,1] <- unlist(wtemp[,7])
######### GLM WATER TEMPERATURE ##########################

# second column is depth
env[,2] <- round(depths,1)

# fifth column is light
swr <- read_csv("./data/cal_met_GLM.csv")
env[,5] <- exp(as.double(log(swr[6,2])) - 0.5*depths)
colnames(env) <- c("wt","yloc","dens","visc","light")





# ---- Start the simulation 
ts         <- 0;
time_steps <- 60*18;
inds_hist  <- NULL;
start_time <- Sys.time()
while(ts < time_steps){
  inds            <- movement(inds, env); 
  inds            <- birth(inds);
  inds            <- death(inds);
  ts              <- ts + 1; 
  env             <- update_env(env, wtemp, swr, time_steps = time_steps, tstep = ts);
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

ind_yloc <- NULL
ind_yloc[[1]] <- array(data = NA, dim = c(time_steps,length(depths)+1))
colnames(ind_yloc[[1]]) <- c("timestep",depths)

ind_yloc[[2]] <- array(data = NA, dim = c(time_steps,length(depths)+1))
colnames(ind_yloc[[2]]) <- c("timestep",depths)

ind_yloc[[3]] <- array(data = NA, dim = c(time_steps,length(depths)+1))
colnames(ind_yloc[[3]]) <- c("timestep",depths)

for(i in 1:time_steps){
  for(j in 1:length(depths)){
    
    #diatom
    ind_yloc[[1]][i, 1] <- i;                      # Save the time step
    ind_yloc[[1]][i, j+1] <- length(inds_hist[[i]][which(round(inds_hist[[i]][,2],1) == round(depths[j],1) & inds_hist[[i]][,1] == 1),1]); # Save the number of individuals at each depth
    
    #cyano
    ind_yloc[[2]][i, 1] <- i;                      # Save the time step
    ind_yloc[[2]][i, j+1] <- length(inds_hist[[i]][which(round(inds_hist[[i]][,2],1) == round(depths[j],1) & inds_hist[[i]][,1] == 2),1]); # Save the number of individuals at each depth
    
    #green
    ind_yloc[[3]][i, 1] <- i;                      # Save the time step
    ind_yloc[[3]][i, j+1] <- length(inds_hist[[i]][which(round(inds_hist[[i]][,2],1) == round(depths[j],1) & inds_hist[[i]][,1] == 3),1]); # Save the number of individuals at each depth
    
  }
}
#print(ind_yloc);




# =============================================================================
# Wrangle and save the results
# =============================================================================

##plot IBM output
plot_yloc.d <- data.frame(ind_yloc[[1]]) %>%
  gather(X0.1:X9.3, key = "Depth_m",value = "num_agents") %>%
  mutate(Depth_m = as.double(substring(Depth_m, 2)))
#head(plot_yloc.d)

plot_ts.d <- plot_yloc.d %>%
  filter(Depth_m >= 0.1 & Depth_m <= 2) %>%
  group_by(timestep) %>%
  summarize(surface_agents = mean(num_agents))

plot_yloc.c <- data.frame(ind_yloc[[2]]) %>%
  gather(X0.1:X9.3, key = "Depth_m",value = "num_agents") %>%
  mutate(Depth_m = as.double(substring(Depth_m, 2)))
#head(plot_yloc.c)

plot_ts.c <- plot_yloc.c %>%
  filter(Depth_m >= 0.1 & Depth_m <= 2) %>%
  group_by(timestep) %>%
  summarize(surface_agents = mean(num_agents))

plot_yloc.g <- data.frame(ind_yloc[[3]]) %>%
  gather(X0.1:X9.3, key = "Depth_m",value = "num_agents") %>%
  mutate(Depth_m = as.double(substring(Depth_m, 2)))
#head(plot_yloc.g)

plot_ts.g <- plot_yloc.g %>%
  filter(Depth_m >= 0.1 & Depth_m <= 2) %>%
  group_by(timestep) %>%
  summarize(surface_agents = mean(num_agents))

plot_ts.all <- left_join(plot_ts.d, plot_ts.c, by = c("timestep")) %>%
  left_join(plot_ts.g, by = c("timestep")) %>%
  rename(diatom = surface_agents.x,
         cyano = surface_agents.y,
         green = surface_agents) %>%
  mutate(total = diatom + cyano + green) %>%
  gather(diatom:total, key = pft, value = num_agents)

plot_ts.prop <- left_join(plot_ts.d, plot_ts.c, by = c("timestep")) %>%
  left_join(plot_ts.g, by = c("timestep")) %>%
  rename(diatom = surface_agents.x,
         cyano = surface_agents.y,
         green = surface_agents) %>%
  mutate(total = diatom + cyano + green) %>%
  mutate(total_prop = total/total*100,
         diatom_prop = diatom/total*100,
         cyano_prop = cyano/total*100,
         green_prop = green/total*100) %>%
  select(timestep, diatom_prop, cyano_prop, green_prop, total_prop) %>%
  gather(diatom_prop:total_prop, key = pft, value = prop_agents)

# #Write output to file for plotting for proposal
write.csv(plot_yloc.d, file = "./ABM_output/ABM_diatom_depthByTimestep_18hr.csv", row.names = FALSE)
write.csv(plot_ts.d, file = "./ABM_output/ABM_diatom_surfaceTimeseries_18hr.csv", row.names = FALSE)

write.csv(plot_yloc.c, file = "./ABM_output/ABM_cyano_depthByTimestep_18hr.csv", row.names = FALSE)
write.csv(plot_ts.c, file = "./ABM_output/ABM_cyano_surfaceTimeseries_18hr.csv", row.names = FALSE)

write.csv(plot_yloc.g, file = "./ABM_output/ABM_green_depthByTimestep_18hr.csv", row.names = FALSE)
write.csv(plot_ts.g, file = "./ABM_output/ABM_green_surfaceTimeseries_18hr.csv", row.names = FALSE)

write.csv(plot_ts.all, file = "./ABM_output/ABM_all_surfaceTimeseries_18hr.csv",row.names = FALSE)
write.csv(plot_ts.prop, file = "./ABM_output/ABM_allProp_surfaceTimeseries_18hr.csv",row.names = FALSE)

saveRDS(inds_hist, file = "ABM_output_18hr.rds")
