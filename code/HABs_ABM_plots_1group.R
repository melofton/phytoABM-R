#Title: Plot for Societal RoL proposal
#Date: 02FEB23
#Author: MEL

#load packages
library(tidyverse)
library(lubridate)
library(cowplot)


# =============================================================================
# Plot the results
# =============================================================================

# Plot ABM output

ts_keep <- seq(from = 0, to = 120, by = 1)

plot_yloc <- read_csv("./model_output/ABM_depthByTimestep_2hr.csv") %>%
  filter(timestep %in% ts_keep)

DepthTime <- ggplot(data = plot_yloc, aes(x = num_agents, y = Depth_m, group = timestep, color = timestep))+
  geom_path()+
  scale_y_reverse()+
  theme_classic()+
  ggtitle("")+
  xlab("Individuals")+
  ylab("Depth (m)")+
  labs(color = "Time (minutes)")
  #theme(legend.position = "none")
DepthTime

ggsave(DepthTime, filename = "./plot_output/depth_by_time_1pft.tif",height = 4, width = 3.5,
       units = "in", dpi = 300, dev = "tiff")

dates <- seq(from=as.POSIXct("2021-08-09 05:01"),to=as.POSIXct("2021-08-09 05:01")+120*60,by="min", tz = "UTC") 
plot_dates <- tibble(dates[ts_keep])

plot_ts <- read_csv("./model_output/ABM_agentTimeseries_2hr.csv") %>%
  filter(timestep %in% ts_keep) %>%
  bind_cols(plot_dates) %>%
  rename(datetime = `dates[ts_keep]`)

agents_ts <- ggplot(data = plot_ts, aes(x = datetime, y = surface_agents)) +
  geom_line()+
  theme_bw()+
  ylab("Number of agents")+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "1 hours")+
  xlab("Time of day")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
agents_ts

ggsave(agents_ts, filename = "./plot_output/agent_timeseries_1pft.tif",height = 3, width = 10,
       units = "in", dpi = 300, dev = "tiff")

# get distribution of vertical movement over time
inds_hist <- readRDS("./model_output/ABM_output_2hr.rds")

lake_depths <- seq(from = 0.1, to = 9.3, by = 0.1)
time_steps <- 60*2;

ind_zmove <- array(data = NA, dim = c(time_steps,length(lake_depths)+1))
colnames(ind_zmove) <- c("timestep",lake_depths)

for(i in 1:time_steps){
  for(j in 1:length(lake_depths)){
    
    #generic PFT
    ind_zmove[i, 1] <- i;                      # Save the time step
    ind_zmove[i, j+1] <- mean(inds_hist[[i]][which(round(inds_hist[[i]][,2],1) == round(lake_depths[j],1) & inds_hist[[i]][,1] == 1),9], na.rm = TRUE); # Save the number of individuals at each depth
  }
}

#print(ind_yloc);

##dataframes for plotting
plot_zmove <- data.frame(ind_zmove) %>%
  gather(X0.1:X9.3, key = "Depth_m",value = "mean_zmove") %>%
  mutate(Depth_m = as.double(substring(Depth_m, 2))) %>%
  group_by(Depth_m) %>%
  summarize(grand_mean_zmove = mean(mean_zmove, na.rm = TRUE))

zmove_plot <- ggplot(data = plot_zmove, aes(x = Depth_m, y = grand_mean_zmove))+
  geom_bar(stat = "identity", color = "darkblue", fill = "lightblue")+
  xlab("Depth (m)")+
  ylab("Average vertical distance moved")+
  theme_classic()
zmove_plot

ggsave(zmove_plot, filename = "./plot_output/mean_zmove.tif",height = 3, width = 4,
       units = "in", dpi = 300, dev = "tiff")
