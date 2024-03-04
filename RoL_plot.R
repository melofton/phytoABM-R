#Title: Plot for Societal RoL proposal
#Date: 02FEB23
#Author: MEL

#load packages
library(tidyverse)
library(lubridate)

#read in and wrangle data

#GLM
glm <- read_csv("GLM_output_summer_2021.csv")

GLM_day <- glm %>%
  filter(date(DateTime) == "2021-08-09") %>%
  select(-PHY_tchla_1.6, -temp_1.6) %>%
  gather(PHY_tphy_1.6:PHY_cyano_1.6, key = pft, value = mmolCm3)

glm_prop <- glm %>%
  filter(date(DateTime) == "2021-08-09") %>%
  mutate(total_prop = PHY_tphy_1.6/PHY_tphy_1.6*100,
         diatom_prop = PHY_diatom_1.6/PHY_tphy_1.6*100,
         cyano_prop = PHY_cyano_1.6/PHY_tphy_1.6*100,
         green_prop = PHY_green_1.6/PHY_tphy_1.6*100) %>%
  select(DateTime, total_prop:green_prop) %>%
  gather(total_prop:green_prop, key = "pft", value = "prop_mmolCm3")

glm_join <- glm %>%
  filter(date(DateTime) == "2021-08-09") %>%
  select(DateTime, PHY_tchla_1.6)

#FP
fp <- read_csv("./HABs-ABM/data/FluoroProbe_2014_2021.csv") %>%
  filter(Reservoir == "FCR" & Site == 50 & year(DateTime) == 2021 & month(DateTime) %in% c(6:9))

fp2 <- fp %>%
  filter(Depth_m <= 9.3) %>%
  mutate(TotalConcNoMixed_ugL = TotalConc_ugL - MixedAlgae_ugL) %>%
  select(DateTime, Depth_m, GreenAlgae_ugL, Bluegreens_ugL, BrownAlgae_ugL, TotalConcNoMixed_ugL, Temp_degC)

prof <- fp2 %>%
  filter(date(DateTime) == "2021-08-09") %>%
  gather(GreenAlgae_ugL:TotalConcNoMixed_ugL, key = "spectral_group", value = "ugL")

fp_surface <- prof %>%
  filter(Depth_m <= 2) %>%
  group_by(spectral_group) %>%
  summarize(surface_ugL = mean(ugL, na.rm = TRUE))
fp_surface[,"surface_prop"] <- fp_surface[,"surface_ugL"]/unlist(fp_surface[4,2])
fp_surface[,"pft"] <- c("cyano_prop","diatom_prop","green_prop","total_prop")

#EXO
exo <- read_csv("./HABs-ABM/data/FCR_catwalk_EDI_2018_2022.csv") %>%
  filter(year(DateTime) == 2021 & month(DateTime) %in% c(6:9)) %>%
  select(DateTime, EXOChla_ugL_1)

exo_day <- exo %>%
  filter(date(DateTime) == "2021-08-09")
exo_no_outliers <- exo %>%
  filter(date(DateTime) == "2021-08-09") %>%
  filter(EXOChla_ugL_1 < 15)

#ABM

#function to get model datetimes
dates <- rep(seq(from=as.POSIXct("2021-08-09 05:01"),to=as.POSIXct("2021-08-09 23:00"),by="min", tz = "UTC"),4)

abm <- read_csv("./HABs-ABM/ABM_output/ABM_all_surfaceTimeseries_18hr.csv") %>%
  mutate(DateTime = dates)

abm_prop <- read_csv("./HABs-ABM/ABM_output/ABM_allProp_surfaceTimeseries.csv") %>%
  mutate(DateTime = dates) %>%
  mutate(prop_agents = prop_agents/100) %>%
  mutate(facet = "pretty")

abm_join <- abm %>%
  filter(pft == "total") %>%
  mutate(DateTime = force_tz(DateTime, tzone = "UTC")) %>%
  select(DateTime, num_agents)

#plotting

ggplot(data = GLM_day, aes(x = DateTime, y = mmolCm3, group = pft, color = pft))+
  geom_line()+
  theme_bw()

ggplot(data = exo_day, aes(x = DateTime, y = EXOChla_ugL_1))+
  geom_line()+
  #ylim(5,15)+
  theme_bw()

ggplot(data = abm, aes(x = timestep, y = num_agents, group = pft, color = pft)) +
  geom_line()+
  theme_bw()

ggplot() +
  geom_line(data = abm_prop, aes(x = DateTime, y = prop_agents, group = pft, color = pft))+
  geom_point(aes(x = rep(as.POSIXct("2021-08-09 13:45"),4), y = fp_surface$surface_prop, group = fp_surface$pft, color = fp_surface$pft))+
  theme_bw()

ggplot(data = glm_prop, aes(x = DateTime, y = prop_mmolCm3, group = pft, color = pft)) +
  geom_line()+
  theme_bw()

money_df <- tibble(DateTime = seq(from=as.POSIXct("2021-08-09 06:00"),to=as.POSIXct("2021-08-09 20:00"),by="min", tz = "UTC")) %>%
  mutate(DateTime = force_tz(DateTime, tzone = "UTC")) %>%
  left_join(abm_join, by = "DateTime") %>%
  left_join(exo_no_outliers, by = "DateTime") %>%
  left_join(glm_join, by = "DateTime") %>%
  rename(ABM = num_agents) %>%
  gather(ABM:PHY_tchla_1.6, key = data_source, value = phyto) %>%
  mutate(data_source = factor(data_source, levels = c("EXOChla_ugL_1","ABM","PHY_tchla_1.6")))

date.breaks = as.POSIXct(c("2021-08-09 06:00","2021-08-09 09:00",
                           "2021-08-09 12:00","2021-08-09 15:00",
                           "2021-08-09 18:00"))

base_money_df <- tibble(DateTime = seq(from=as.POSIXct("2021-08-09 06:00"),to=as.POSIXct("2021-08-09 19:00"),by="min", tz = "UTC")) %>%
  mutate(DateTime = force_tz(DateTime, tzone = "UTC")) %>%
  left_join(abm_join, by = "DateTime") %>%
  left_join(exo_no_outliers, by = "DateTime") %>%
  left_join(glm_join, by = "DateTime") %>%
  rename(ABM = num_agents) 

jpeg("./HABs-ABM/ABM_vs_obs.jpeg", res = 300, width = 5, height = 3, units = "in")
par(mar = c(4, 4, 1, 4) + 0.1, mgp = c(2.4,1,0),cex.lab = 1.3)              # Additional space for second y-axis
plot(base_money_df$DateTime, base_money_df$EXOChla_ugL_1, 
     type = "p", pch = 16, col = "green", xlab = "Time of day", 
     ylab = expression(paste("Observed chl-a"," (",mu,"g ",L^-1,")")), 
     yaxt = "n")
axis(2, las = 2)
par(new = TRUE)      
par(mar = c(4, 4, 1, 4) + 0.1, mgp = c(2.4,1,0),cex.lab = 1.3)              # Additional space for second y-axis# Add new plot
plot(base_money_df$DateTime, base_money_df$ABM, axes = FALSE, xlab = "", 
     ylab = "", type = "l", lwd = 2)
axis(side = 4, at = pretty(range(base_money_df$ABM)), las = 2)      # Add second axis
mtext(expression("IBM (individuals)"), side = 4, line = 3, cex = 1.3)             # Add second axis label
legend("topleft", lty = c(NA,1), pch = c(16,NA), lwd = c(NA,2), col = c("green","black"),legend = c("Observed","Modeled"),bty = "n")
dev.off()

par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(x, y1, pch = 16, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(x, y2, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(y2)))      # Add second axis
mtext("y2", side = 4, line = 3)             # Add second axis label


fac.labs <- c("IBM","Observed","Process-based")
names(fac.labs) <- c("ABM","EXOChla_ugL_1","PHY_tchla_1.6")

ann_text <- data.frame(DateTime = rep(as.POSIXct("2021-08-09 02:00"),times = 3),
                       phyto = c(12.5, 162.5, 4.65),
                       lab = c("A", "B", "C"),
                       data_source = factor(c("EXOChla_ugL_1","ABM","PHY_tchla_1.6"),levels = c("EXOChla_ugL_1","ABM","PHY_tchla_1.6")))

compare <- ggplot(data = money_df, aes(x = DateTime, y = phyto))+
  geom_point(size = 0.5)+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours")+
  facet_grid(rows = vars(data_source), scales = "free_y", switch = "y",
             labeller = labeller(data_source = fac.labs))+
  theme_bw()+
  theme(strip.placement = "outside", legend.position = "none",
        strip.text.y = element_text(
          size = 12
        ),
        strip.background = element_blank(
        ),
        axis.text.x = element_text(size = 9),
        plot.margin = margin(b = 0, r = 5, t = 5))+
        #axis.ticks.x = element_blank()
  ylab("")+
  xlab("")+
  geom_text(data = ann_text,label = ann_text$lab, size = 6)
compare

ggsave(compare, filename = "./HABs-ABM/IBM_plot_v2.tif",height = 4.5, width = 3.5,
       units = "in", dpi = 300, dev = "tiff")

fac.labs2 <- c("Proportion")
names(fac.labs2) <- c("pretty")

groups <- ggplot() +
  geom_line(data = abm_prop, aes(x = DateTime, y = prop_agents, group = pft, color = pft, linetype = "IBM"))+
  geom_point(aes(x = rep(as.POSIXct("2021-08-09 13:45"),4), y = fp_surface$surface_prop, group = fp_surface$pft, color = fp_surface$pft, shape = "Observed"))+
  theme_bw()+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours")+
  ylab("Proportion of phytoplankton")+
  xlab("Time of Day")+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        strip.placement = "outside",
        strip.text.x = element_text(
          size = 12, face = "bold"
        ),
        strip.background = element_blank(
        ),
        legend.box = "horizontal",
        legend.margin=margin(),
        plot.margin = margin(t = 0,r = 5,b = 5),
        legend.spacing.y = unit(0, 'cm'),
        plot.title = element_blank())+
  ylab("")+
  scale_color_manual(labels = c("Cyanobacteria","Brown algae","Green algae","Total"),
                    values = c("darkcyan","brown","olivedrab3","darkblue"))+
  scale_linetype_manual(name = "", values = c(IBM = "solid"))+
  scale_shape_manual(name = "", values = c(Observed = 16))+
  facet_grid(rows = vars(facet), scales = "free_y", switch = "y",
             labeller = labeller(facet = fac.labs2))+
  guides(color=guide_legend(nrow=4,byrow=TRUE))
  
groups

library(cowplot)

final_plot<-plot_grid(compare, groups, scale = 1,
                nrow = 2, ncol = 1, rel_heights = c(2.5,1.7),
                align = "v")
ggsave(final_plot, filename = "./HABs-ABM/IBM_plot.tif",height = 7, width = 4,
       units = "in", dpi = 300, dev = "tiff")

