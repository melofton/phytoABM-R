# Title: Plots of ABM model driver data
# Author: Mary Lofton
# Date: 04MAR23

# Purpose: visualize model driver data (which is GLM-AED driver data + model output) 
# for ABM model runs

# load packages
library(tidyverse)
library(lubridate)
library(akima)
library(colorRamps)

# read in and wrangle data

#GLM temperature and shortwave
temp <- read_csv("./data/cal_wtemp_GLM.csv")

swr <- read_csv("./data/cal_met_GLM.csv") %>%
  mutate(Hour = hour(time))

# plotting functions
temp_heatmap <- function(temp_data){
  
  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  temp_new <- temp_data %>%
    pivot_longer(`0`:`23`, names_to = "Hour", values_to = "Temp_C") %>%
    rename(Depth_m = depth) %>%
    mutate(Hour = as.numeric(Hour))
  
  interp <- akima::interp(x=temp_new$Hour, y = temp_new$Depth_m, z = temp_new$Temp_C,
                   xo = seq(min(temp_new$Hour), max(temp_new$Hour), by = 0.1), 
                   yo = seq(min(temp_new$Depth_m), max(temp_new$Depth_m), by = 0.01),
                   extrap = T, linear = T, duplicate = "strip")
  interp <- akima::interp2xyz(interp, data.frame=T)
  
  fig_title <- "GLM-AED modeled temperature output for 2021-08-09"
  
  p1 <- ggplot(interp, aes(x=x, y=y))+
    geom_raster(aes(fill=z))+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = colorRamps::blue2green2red(60), na.value="gray")+
    labs(x = "Hour of day", y = "Depth (m)", title = fig_title,fill='Water temperature (°C)')+
    theme_bw()
  
  print(p1)
  
}

temp_plot <- temp_heatmap(temp_data = temp)
ggsave(temp_plot, filename = "./plot_output/temp_environment.png", device = "png",
       width = 6, height = 4, units = "in")

swr_plot <- ggplot(data = swr, aes(x = Hour, y = ShortWave))+
  geom_line(linewidth = 1)+
  theme_classic()+
  xlab("Hour of day")+
  ylab("Shortwave radiation (W/m2)")
ggsave(swr_plot, filename = "./plot_output/swr_environment.png", device = "png",
       width = 6, height = 4, units = "in")
