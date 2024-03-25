#Title: Data munging
#Author: Mary Lofton
#Date: last updated 20DEC22

#Purpose: Download all published data from Environmental Data Initiative (EDI) 
#and tidy it so it can be pulled in for analysis

library(tidyverse)
library(lubridate)
library(data.table)

##Data download----
options(timeout = 1000)
#download FP data from EDI 
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")
dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "Reservoir",     
                 "Site",     
                 "DateTime",     
                 "CastID",     
                 "Depth_m",     
                 "GreenAlgae_ugL",     
                 "Bluegreens_ugL",     
                 "BrownAlgae_ugL",     
                 "MixedAlgae_ugL",     
                 "TotalConc_ugL",     
                 "YellowSubstances_ugL",     
                 "Temp_C",     
                 "Transmission_perc",     
                 "RFU_370nm",     
                 "RFU_470nm",     
                 "RFU_525nm",     
                 "RFU_570nm",     
                 "RFU_590nm",     
                 "RFU_610nm",     
                 "Flag_GreenAlgae_ugL",     
                 "Flag_Bluegreens_ugL",     
                 "Flag_BrownAlgae_ugL",     
                 "Flag_MixedAlgae_ugL",     
                 "Flag_YellowSubstances_ugL",     
                 "Flag_TotalConc_ugL",     
                 "Flag_Temp_C",     
                 "Flag_Transmission_perc",     
                 "Flag_RFU_525nm",     
                 "Flag_RFU_570nm",     
                 "Flag_RFU_610nm",     
                 "Flag_RFU_370nm",     
                 "Flag_RFU_590nm",     
                 "Flag_RFU_470nm"    ), check.names=TRUE)

unlink(infile1)
write.csv(dt1, "./data/FluoroProbe_2014_2023.csv", row.names = FALSE)

#download EXO data from EDI 
data  <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.271.7&entityid=71e6b946b751aa1b966ab5653b01077f"

download.file(data,destfile = "./HABs-ABM/data/FCR_catwalk_EDI_2018_2022.csv", method='libcurl')

##Data munging----

#FP----
fp <- read_csv("./HABs-ABM/data/FluoroProbe_2014_2021.csv") %>%
  filter(Reservoir == "FCR" & Site == 50 & year(DateTime) == 2021 & month(DateTime) %in% c(8:12))

fp2 <- fp %>%
  filter(date(DateTime) == "2021-08-02" & Depth_m <= 9.3) %>%
  mutate(TotalConcNoMixed_ugL = TotalConc_ugL - MixedAlgae_ugL) %>%
  select(DateTime, Depth_m, GreenAlgae_ugL, Bluegreens_ugL, BrownAlgae_ugL, TotalConcNoMixed_ugL, Temp_degC)

# Assign all depths in between 0.1 m increments to the relevant increment (e.g., anything 
# between 0.06 - 0.15 is 0.1; anything between 0.16 and 0.25 is 0.2; etc.)
fp3 <- fp2 %>%
  mutate(Depth_inc = round(Depth_m, 1)) %>%
  group_by(Depth_inc) %>%
  summarize(GreenAlgae_ugL = as.double(mean(GreenAlgae_ugL, na.rm = TRUE)),
            Bluegreens_ugL = as.double(mean(Bluegreens_ugL, na.rm = TRUE)),
            BrownAlgae_ugL = as.double(mean(BrownAlgae_ugL, na.rm = TRUE)),
            TotalConcNoMixed_ugL = as.double(mean(TotalConcNoMixed_ugL, na.rm = TRUE)))

#actual missing depths: 0.1, 1.6, 2.9, 6.5, 8.5, 8.6, 9.3
depths <- data.frame(Depth_inc = c(0.1, 1.6, 2.9, 6.5, 8.5, 8.6, 9.3),
                     GreenAlgae_ugL = rep(NA,7),
                     Bluegreens_ugL = rep(NA,7),
                     BrownAlgae_ugL = rep(NA,7),
                     TotalConcNoMixed_ugL = rep(NA,7))
depths <- tibble(depths)

fp4 <- bind_rows(depths,fp3) %>%
  arrange(Depth_inc)
fp4[1,c(2:5)] <- fp4[2,c(2:5)]
fp4[16,c(2:5)] <- fp4[17,c(2:5)]
fp4[29,c(2:5)] <- fp4[30,c(2:5)]
fp4[65,c(2:5)] <- fp4[66,c(2:5)]
fp4[c(85:86),c(2:5)] <- fp4[87,c(2:5)]
fp4[93,c(2:5)] <- fp4[92,c(2:5)]

ggplot(data = fp4, aes(x = TotalConcNoMixed_ugL, y = Depth_inc))+
  geom_path()+
  scale_y_reverse()+
  theme_classic()

write.csv(fp4,"./HABs-ABM/data/FluoroProbe_2021-08-02_FCR_50.csv",row.names = FALSE)

#plot 2021 data
fp <- read_csv("./HABs-ABM/data/FluoroProbe_2014_2021.csv") %>%
  filter(Reservoir == "FCR" & Site == 50 & year(DateTime) == 2021 & month(DateTime) %in% c(6:9))

fp2 <- fp %>%
  filter(Depth_m <= 9.3) %>%
  mutate(TotalConcNoMixed_ugL = TotalConc_ugL - MixedAlgae_ugL) %>%
  select(DateTime, Depth_m, GreenAlgae_ugL, Bluegreens_ugL, BrownAlgae_ugL, TotalConcNoMixed_ugL, Temp_degC)

fp_ts <- fp2 %>%
  filter(Depth_m >= 1.3 & Depth_m <= 1.8) %>%
  mutate(Date = date(DateTime)) %>%
  group_by(Date) %>%
  summarize(TotalConcNoMixed_ugL = mean(TotalConcNoMixed_ugL, na.rm = TRUE),
            GreenAlgae_ugL = mean(GreenAlgae_ugL, na.rm = TRUE),
            Bluegreens_ugL = mean(Bluegreens_ugL, na.rm = TRUE),
            BrownAlgae_ugL = mean(BrownAlgae_ugL, na.rm = TRUE)) %>%
  gather(TotalConcNoMixed_ugL:BrownAlgae_ugL, key = "spectral_group", value = "ugL")

ggplot(data = fp_ts, aes(x = Date, y = ugL, group = spectral_group, color = spectral_group))+
  geom_line()+
  geom_point()+
  theme_bw()

prof <- fp2 %>%
  filter(date(DateTime) == "2021-08-09") %>%
  gather(GreenAlgae_ugL:TotalConcNoMixed_ugL, key = "spectral_group", value = "ugL")

ggplot(data = prof, aes(x = ugL, y = Depth_m, group = spectral_group, color = spectral_group))+
  geom_path()+
  scale_y_reverse()+
  theme_bw()

fp_surface <- prof %>%
  filter(Depth_m <= 2) %>%
  group_by(spectral_group) %>%
  summarize(surface_ugL = mean(ugL, na.rm = TRUE))

#EXO----
exo <- read_csv("./HABs-ABM/data/FCR_catwalk_EDI_2018_2022.csv") %>%
  filter(year(DateTime) == 2021 & month(DateTime) %in% c(6:9)) %>%
  select(DateTime, EXOChla_ugL_1)

ggplot(data = exo, aes(x = DateTime, y = EXOChla_ugL_1))+
  geom_line()+
  theme_bw()

exo_day <- exo %>%
  filter(date(DateTime) == "2021-08-09")

ggplot(data = exo_day, aes(x = DateTime, y = EXOChla_ugL_1))+
  geom_line()+
  theme_bw()
##Get model output ----

# Load packages
library(sp) 
library(devtools)
library(glmtools) 
library(GLMr) 
library(tidyverse)
library(lubridate)
# If this worked, GLMr should load without error messages. Hooray!

# See what version of GLM you are running- should be v.2.x.x
glm_version() 

# Set nml file name
sim_folder <- paste0("/Users/MaryLofton/RProjects/FCR-phyto-modeling/Eco-KGML-transfer-learning/data/data_raw/ModelOutputFCR") 
nc_file <- file.path(sim_folder, 'output.nc') 

# Get vars
#sim_vars(nc_file)
depths = 1.6
vars = c("temp","PHY_tchla","PHY_tphy","PHY_green","PHY_diatom","PHY_cyano")

var <- get_var(nc_file, var_name = vars[1], reference="surface", z_out=depths) %>%
  filter(year(DateTime) == 2021 & month(DateTime) %in% c(6:9)) 

temp <- var

for(k in 2:length(vars)){
  
var <- get_var(nc_file, var_name = vars[k], reference="surface", z_out=depths) %>%
  filter(year(DateTime) == 2021 & month(DateTime) %in% c(6:9)) 

temp <- left_join(temp,var,by = "DateTime")

}

write.csv(temp, file = "GLM_output_summer_2021.csv",row.names = FALSE)

#Plot vars
vars2 <- paste0(vars, "_1.6")

for (j in 1:length(vars2)){
  
  p <- ggplot(data = temp, aes(x = DateTime, y = .data[[vars2[j]]]))+
    geom_line()+
    theme_bw()
  print(p)
  
}

GLM_day <- temp %>%
  filter(date(DateTime) == "2021-08-09") %>%
  gather(PHY_tphy_1.6:PHY_cyano_1.6, key = pft, value = mmolCm3)

ggplot(data = GLM_day, aes(x = DateTime, y = mmolCm3, group = pft, color = pft))+
  geom_line()+
  theme_bw()

#Pull water temp profiles for calibration date (2021-08-09)
depths <- seq(from = 0.1, to = 9.3, by = 0.1)

wtemp <- get_var(nc_file, var_name = "temp", reference="surface", z_out=depths) %>%
  filter(date(DateTime) == "2021-08-09")

cal_wtemp <- wtemp %>%
  gather(temp_0.1:temp_9.3, key = "Depth_m", value = "Temp_C") %>%
  separate(Depth_m, sep = "_", c("var","depth")) %>%
  select(-var)

cal_wtemp2 <- cal_wtemp %>%
  mutate(hour = hour(DateTime)) %>%
  select(-DateTime) %>%
  spread(key = hour, value = Temp_C)
cal_wtemp2[93,c(2:25)] <- cal_wtemp2[92,c(2:25)]

write.csv(cal_wtemp2, file = "./HABs-ABM/data/cal_wtemp_GLM.csv",row.names = FALSE)

##Get met data ----
met <- read_csv("./data/met_avg_filtered.csv") %>%
  select(time, ShortWave, WindSpeed) %>%
  filter(date(time) == "2021-08-09")
met[6,2] <- (met[5,2]+met[7,2])/2
met[20,2] <- (met[19,2]+met[21,2])/2

write.csv(met, file = "./data/cal_met_GLM.csv",row.names = FALSE)

##Statistics ----

#RMSE: GLM-AED vs. EXO chl-a
exo_daily <- exo %>%
  mutate(Date = date(DateTime)) %>%
  group_by(Date) %>%
  summarize(exo_chla = mean(EXOChla_ugL_1, na.rm = TRUE))

glm_chla <- temp %>%
  mutate(Date = date(DateTime)) %>%
  group_by(Date) %>%
  summarize(glm_chla = mean(PHY_tchla_1.6, na.rm = TRUE))

#define function
rmse = function(m, o){
  sqrt(mean((m - o)^2, na.rm = TRUE))
}

#plot
plot(x = exo_daily$Date, y = exo_daily$exo_chla, type = "p",pch =16)
lines(x = exo_daily$Date, y = glm_chla$glm_chla, col = "red")

RMSE <- rmse(glm_chla$glm_chla, exo_daily$exo_chla)
