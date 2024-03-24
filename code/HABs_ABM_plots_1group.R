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

ts_keep <- seq(from = 0, to = 19800, by = 60)

plot_yloc <- read_csv("./model_output/ABM_depthByTimestep_330hr.csv") %>%
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

dates <- seq(from=as.POSIXct("2021-08-09 05:01"),to=as.POSIXct("2021-08-09 05:01")+19800*60,by="min", tz = "UTC") 
plot_dates <- tibble(dates[ts_keep])

plot_ts <- read_csv("./model_output/ABM_agentTimeseries_330hr.csv") %>%
  filter(timestep %in% ts_keep) %>%
  bind_cols(plot_dates) %>%
  rename(datetime = `dates[ts_keep]`)

agents_ts <- ggplot(data = plot_ts, aes(x = datetime, y = surface_agents)) +
  geom_line()+
  theme_bw()+
  ylab("Number of agents")+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "8 hours")+
  xlab("Time of day")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
agents_ts

ggsave(agents_ts, filename = "./plot_output/agent_timeseries_1pft.tif",height = 3, width = 10,
       units = "in", dpi = 300, dev = "tiff")
