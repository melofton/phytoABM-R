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
library(rLakeAnalyzer)

# read in and wrangle data

#GLM temperature, met, and nutrients
temp <- read_csv("./data/cal_wtemp_GLM.csv")
din <- read_csv("./data/cal_din_GLM.csv")
frp <- read_csv("./data/cal_frp_GLM.csv")

met <- read_csv("./data/cal_met_GLM.csv") %>%
  mutate(Hour = hour(time))

# plotting functions
temp_heatmap <- function(temp_data){
  
  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  temp_new <- temp_data %>%
    pivot_longer(`0`:`23`, names_to = "Hour", values_to = "Temp_C") %>%
    rename(Depth_m = depth) %>%
    mutate(Hour = as.numeric(Hour))
  
  td <- temp_data %>%
    mutate(across(`0`:`23`, ~ thermo.depth(.x,depths = depth))) %>%
    first() %>%
    select(-depth) %>%
    pivot_longer(`0`:`23`,names_to = "hour", values_to = "thermo.depth.m") %>%
    mutate(hour = as.numeric(hour))
  
  interp <- akima::interp(x=temp_new$Hour, y = temp_new$Depth_m, z = temp_new$Temp_C,
                   xo = seq(min(temp_new$Hour), max(temp_new$Hour), by = 0.1), 
                   yo = seq(min(temp_new$Depth_m), max(temp_new$Depth_m), by = 0.01),
                   extrap = T, linear = T, duplicate = "strip")
  interp <- akima::interp2xyz(interp, data.frame=T)
  
  fig_title <- "GLM-AED modeled temperature output for 2021-08-09"
  
  p1 <- ggplot()+
    geom_raster(data = interp, aes(x=x, y=y, fill = z))+
    geom_line(data = td, aes(x=hour, y = thermo.depth.m, color = "Thermocline"), linewidth = 1)+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = colorRamps::blue2green2red(60), na.value="gray")+
    scale_color_manual(values = c("Thermocline" = "black"), name = "")+
    labs(x = "Hour of day", y = "Depth (m)", title = fig_title,fill='Water temperature (°C)')+
    theme_bw()
  
  print(p1)
  
}

temp_plot <- temp_heatmap(temp_data = temp)
ggsave(temp_plot, filename = "./plot_output/temp_environment.png", device = "png",
       width = 6, height = 4, units = "in")

swr_plot <- ggplot(data = met, aes(x = Hour, y = ShortWave))+
  geom_line(linewidth = 1)+
  theme_classic()+
  xlab("Hour of day")+
  ylab("Shortwave radiation (W/m2)")
swr_plot
ggsave(swr_plot, filename = "./plot_output/swr_environment.png", device = "png",
       width = 6, height = 4, units = "in")

wnd_plot <- ggplot(data = met, aes(x = Hour, y = WindSpeed))+
  geom_line(linewidth = 1)+
  theme_classic()+
  xlab("Hour of day")+
  ylab("WindSpeed (m/s)")
wnd_plot
ggsave(wnd_plot, filename = "./plot_output/wnd_environment.png", device = "png",
       width = 6, height = 4, units = "in")

#density plot
density_heatmap <- function(temp_data = temp){
  
  density <- temp_data %>%
    mutate(across(!depth, water.density))
  
  td <- temp_data %>%
    mutate(across(`0`:`23`, ~ thermo.depth(.x,depths = depth))) %>%
    first() %>%
    select(-depth) %>%
    pivot_longer(`0`:`23`,names_to = "hour", values_to = "thermo.depth.m") %>%
    mutate(hour = as.numeric(hour))
  
  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  density_new <- density %>%
    pivot_longer(`0`:`23`, names_to = "Hour", values_to = "Density_kgm3") %>%
    rename(Depth_m = depth) %>%
    mutate(Hour = as.numeric(Hour))
  
  interp <- akima::interp(x=density_new$Hour, y = density_new$Depth_m, z = density_new$Density_kgm3,
                          xo = seq(min(density_new$Hour), max(density_new$Hour), by = 0.1), 
                          yo = seq(min(density_new$Depth_m), max(density_new$Depth_m), by = 0.01),
                          extrap = T, linear = T, duplicate = "strip")
  interp <- akima::interp2xyz(interp, data.frame=T)
  
  fig_title <- "Calculated water density for 2021-08-09"
  
  p1 <- ggplot()+
    geom_raster(data = interp, aes(x=x, y=y, fill = z))+
    geom_line(data = td, aes(x=hour, y = thermo.depth.m, color = "Thermocline"), linewidth = 1)+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(colorRamps::blue2yellow(60)), na.value="gray")+
    scale_color_manual(values = c("Thermocline" = "black"), name = "")+
    labs(x = "Hour of day", y = "Depth (m)", title = fig_title,fill='Water density (kg/m3)')+
    theme_bw()
  
  print(p1)
  
}

density_plot <- density_heatmap(temp_data = temp)
ggsave(density_plot, filename = "./plot_output/density_environment.png", device = "png",
       width = 6, height = 4, units = "in")

#diffusivity plot
wnd <- met %>%
  select(-ShortWave)
ymax = 9.3

K_heatmap <- function(temp_data = temp, wnd, ymax){
  
  density <- temp_data %>%
    mutate(across(!depth, water.density))
  
  td <- temp_data %>%
    mutate(across(`0`:`23`, ~ thermo.depth(.x,depths = depth))) %>%
    first() %>%
    select(-depth) %>%
    pivot_longer(`0`:`23`,names_to = "hour", values_to = "thermo.depth.m") %>%
    mutate(hour = as.numeric(hour))
  
  avgEpiDense = density %>%
    pivot_longer(`0`:`23`,names_to = "hour", values_to = "density") %>%
    mutate(hour = as.numeric(hour)) %>%
    left_join(td) %>%
    group_by(hour) %>%
    filter(any(depth <= thermo.depth.m)) %>%
    summarize(averageEpiDense = mean(density, na.rm = TRUE))
  
  u_star = uStar(wndSpeed = unlist(wnd$WindSpeed), wndHeight = 3, averageEpiDense = unlist(avgEpiDense$averageEpiDense))
  
  K_data <- array(data = 0, dim = c(length(temp_data$depth), 25))
  K_data[,1] <- temp_data$depth
  
  for(i in 1:length(u_star)){
    K_data[,i+1] <- 0.4*u_star[i]*ymax*((ymax - K_data[,1])/ymax)*(1-((ymax - K_data[,1])/ymax)) + 0.00001
  }
  
  K_data <- data.frame(K_data)
  colnames(K_data) <- c("depth",c(0:23))

  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  K_data_new <- K_data %>%
    pivot_longer(`0`:`23`, names_to = "Hour", values_to = "K_perm2pers") %>%
    rename(Depth_m = depth) %>%
    mutate(Hour = as.numeric(Hour))
  
  interp <- akima::interp(x=K_data_new$Hour, y = K_data_new$Depth_m, z = K_data_new$K_perm2pers,
                          xo = seq(min(K_data_new$Hour), max(K_data_new$Hour), by = 0.1), 
                          yo = seq(min(K_data_new$Depth_m), max(K_data_new$Depth_m), by = 0.01),
                          extrap = T, linear = T, duplicate = "strip")
  interp <- akima::interp2xyz(interp, data.frame=T)
  
  fig_title <- "Calculated K for 2021-08-09"
  
  p1 <- ggplot()+
    geom_raster(data = interp, aes(x=x, y=y, fill = z))+
    geom_line(data = td, aes(x=hour, y = thermo.depth.m, color = "Thermocline"), linewidth = 1)+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(colorRamps::magenta2green(60)), na.value="gray")+
    scale_color_manual(values = c("Thermocline" = "black"), name = "")+
    labs(x = "Hour of day", y = "Depth (m)", title = fig_title,fill='K (m-2 s-1)')+
    theme_bw()
  
  print(p1)
  
}

K_plot <- K_heatmap(temp_data = temp, wnd, ymax)
ggsave(K_plot, filename = "./plot_output/K_environment.png", device = "png",
       width = 6, height = 4, units = "in")


# ustar plot

ustar_plot <- function(temp_data = temp, wnd, ymax){
  
  density <- temp_data %>%
    mutate(across(!Depth_m, water.density))
  
  td <- temp_data %>%
    mutate(across(`0`:`23`, ~ thermo.depth(.x,depths = Depth_m))) %>%
    first() %>%
    select(-Depth_m) %>%
    pivot_longer(`0`:`23`,names_to = "hour", values_to = "thermo.depth.m") %>%
    mutate(hour = as.numeric(hour))
  
  avgEpiDense = density %>%
    pivot_longer(`0`:`23`,names_to = "hour", values_to = "density") %>%
    mutate(hour = as.numeric(hour)) %>%
    left_join(td) %>%
    group_by(hour) %>%
    filter(any(Depth_m <= thermo.depth.m)) %>%
    summarize(averageEpiDense = mean(density, na.rm = TRUE))
  
  u_star = uStar(wndSpeed = unlist(wnd$WindSpeed), wndHeight = 3, averageEpiDense = unlist(avgEpiDense$averageEpiDense))
  
  K_data <- array(data = 0, dim = c(length(temp_data$depth), 25))
  K_data[,1] <- temp_data$depth
  
  for(i in 1:length(u_star)){
    K_data[,i+1] <- 0.4*u_star[i]*ymax*((ymax - K_data[,1])/ymax)*(1-((ymax - K_data[,1])/ymax)) + 0.00001
  }
  
  K_data <- data.frame(K_data)
  colnames(K_data) <- c("depth",c(0:23))
  
  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  K_data_new <- K_data %>%
    pivot_longer(`0`:`23`, names_to = "Hour", values_to = "K_perm2pers") %>%
    rename(Depth_m = depth) %>%
    mutate(Hour = as.numeric(Hour))
  
  interp <- akima::interp(x=K_data_new$Hour, y = K_data_new$Depth_m, z = K_data_new$K_perm2pers,
                          xo = seq(min(K_data_new$Hour), max(K_data_new$Hour), by = 0.1), 
                          yo = seq(min(K_data_new$Depth_m), max(K_data_new$Depth_m), by = 0.01),
                          extrap = T, linear = T, duplicate = "strip")
  interp <- akima::interp2xyz(interp, data.frame=T)
  
  fig_title <- "Calculated K for 2021-08-09"
  
  p1 <- ggplot()+
    geom_raster(data = interp, aes(x=x, y=y, fill = z))+
    geom_line(data = td, aes(x=hour, y = thermo.depth.m, color = "Thermocline"), linewidth = 1)+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(colorRamps::magenta2green(60)), na.value="gray")+
    scale_color_manual(values = c("Thermocline" = "black"), name = "")+
    labs(x = "Hour of day", y = "Depth (m)", title = fig_title,fill='K (m-2 s-1)')+
    theme_bw()
  
  print(p1)
  
}

K_plot <- K_heatmap(temp_data = temp, wnd, ymax)
ggsave(K_plot, filename = "./plot_output/K_environment.png", device = "png",
       width = 6, height = 4, units = "in")

# plotting functions
din_heatmap <- function(din_data, temp_data){
  
  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  din_new <- din_data %>%
    pivot_longer(`0`:`23`, names_to = "Hour", values_to = "DIN") %>%
    mutate(Hour = as.numeric(Hour))
  
  td <- temp_data %>%
    mutate(across(`0`:`23`, ~ thermo.depth(.x,depths = Depth_m))) %>%
    first() %>%
    select(-Depth_m) %>%
    pivot_longer(`0`:`23`,names_to = "hour", values_to = "thermo.depth.m") %>%
    mutate(hour = as.numeric(hour))
  
  interp <- akima::interp(x=din_new$Hour, y = din_new$Depth_m, z = din_new$DIN,
                          xo = seq(min(din_new$Hour), max(din_new$Hour), by = 0.1), 
                          yo = seq(min(din_new$Depth_m), max(din_new$Depth_m), by = 0.01),
                          extrap = T, linear = T, duplicate = "strip")
  interp <- akima::interp2xyz(interp, data.frame=T)
  
  fig_title <- "GLM-AED modeled DIN output for 2021-08-09"
  
  p1 <- ggplot()+
    geom_raster(data = interp, aes(x=x, y=y, fill = z))+
    geom_line(data = td, aes(x=hour, y = thermo.depth.m, color = "Thermocline"), linewidth = 1)+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = colorRamps::blue2green2red(60), na.value="gray")+
    scale_color_manual(values = c("Thermocline" = "black"), name = "")+
    labs(x = "Hour of day", y = "Depth (m)", title = fig_title,fill='DIN (mmol N m-3)')+
    theme_bw()
  
  print(p1)
  
}

din_plot <- din_heatmap(din_data = din, temp_data = temp)
ggsave(din_plot, filename = "./plot_output/din_environment.png", device = "png",
       width = 6, height = 4, units = "in")

# plotting functions
frp_heatmap <- function(frp_data, temp_data){
  
  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  frp_new <- frp_data %>%
    pivot_longer(`0`:`23`, names_to = "Hour", values_to = "FRP") %>%
    mutate(Hour = as.numeric(Hour))
  
  td <- temp_data %>%
    mutate(across(`0`:`23`, ~ thermo.depth(.x,depths = Depth_m))) %>%
    first() %>%
    select(-Depth_m) %>%
    pivot_longer(`0`:`23`,names_to = "hour", values_to = "thermo.depth.m") %>%
    mutate(hour = as.numeric(hour))
  
  interp <- akima::interp(x=frp_new$Hour, y = frp_new$Depth_m, z = frp_new$FRP,
                          xo = seq(min(frp_new$Hour), max(frp_new$Hour), by = 0.1), 
                          yo = seq(min(frp_new$Depth_m), max(frp_new$Depth_m), by = 0.01),
                          extrap = T, linear = T, duplicate = "strip")
  interp <- akima::interp2xyz(interp, data.frame=T)
  
  fig_title <- "GLM-AED modeled FRP output for 2021-08-09"
  
  p1 <- ggplot()+
    geom_raster(data = interp, aes(x=x, y=y, fill = z))+
    geom_line(data = td, aes(x=hour, y = thermo.depth.m, color = "Thermocline"), linewidth = 1)+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = colorRamps::blue2green2red(60), na.value="gray")+
    scale_color_manual(values = c("Thermocline" = "black"), name = "")+
    labs(x = "Hour of day", y = "Depth (m)", title = fig_title,fill='FRP (mmol P m-3)')+
    theme_bw()
  
  print(p1)
  
}

frp_plot <- frp_heatmap(frp_data = frp, temp_data = temp)
ggsave(frp_plot, filename = "./plot_output/frp_environment.png", device = "png",
       width = 6, height = 4, units = "in")



