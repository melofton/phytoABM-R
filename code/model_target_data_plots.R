# Title: Model target data plots
# Author: Mary Lofton
# Date: 25MAR24

# Purpose: develop visualizations of target data

# load packages
library(tidyverse)
library(lubridate)

fp <- read_csv("./data/FluoroProbe_2014_2023.csv") %>%
  filter(Reservoir == "FCR" & Site == 50 & month(DateTime) == 8 & year(DateTime) == 2021) %>%
  select(DateTime, Depth_m, TotalConc_ugL) %>%
  mutate(Date = as.factor(date(DateTime)))

august_dcm <- ggplot(data = fp, aes(x = TotalConc_ugL, y = Depth_m, group = Date, color = Date))+
  geom_path()+
  scale_y_reverse()+
  xlab("Phytoplankton biomass (ug/L)")+
  ylab("Depth (m)")+
  theme_bw()+
  theme(legend.position='none')+
  ggtitle("FCR August depth profiles in 2021")
august_dcm

ggsave(august_dcm, filename = "./plot_output/August_DCM_in_FCR.tif",height = 4, width = 3.5,
       units = "in", dpi = 300, dev = "tiff")
