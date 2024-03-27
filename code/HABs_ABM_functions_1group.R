# Title: HABs AMB functions for model with 1 PFT
# Author: Mary Lofton
# Date: 16MAR24

# Purpose: Library of custom functions to run R-based HABs ABM

# =============================================================================
# Movement function
# =============================================================================
movement <- function(inds, env, yloc = 2, ymax = 9.3, wnd = wnd, delta_z = 9, depths){
  
  total_inds <- dim(inds)[1]
  
  cell_diam <- inds[,3] # Get the cell diameters of all individuals
  cell_dens <- inds[,4] # Get the cell densities of all individuals
  cell_shape <- inds[,5] # Get the cell shapes of all individuals
  
  #calculate water density (Kell, 1975, eq. 16); need to replace this with water.density function from rLakeAnalyzer https://github.com/GLEON/rLakeAnalyzer/blob/master/R/water.density.R
  env[,3] <- water.density(wtr = env[,1], sal = env[,1]*0)

  #calculate water viscosity (Kestin et al. 1978, eq. 15)
  env[,4] <- (10^(((20 - env[,1])/(env[,1] + 96))*(1.2364 - 0.00137*(20 - env[,1]) + 0.0000057*(20 - env[,1])^2)))*1000 #this actually calculates the ratio of viscosity at target temp/viscosity at 20 degrees C, but viscosity at 20 degrees C is ~ 1
  
  #calculate buoyancy frequency (rLakeAnalyzer)
  buoyancy_freq <- buoyancy.freq(wtr = env[,1], depths = lake_depths)
  env[,10] <- c(first(buoyancy_freq),buoyancy_freq)
  
  #calculate thermocline depth (rLakeAnalyzer)
  td = thermo.depth(wtr = env[,1], depths = env[,2])
  
  #calculate N
  N <- sqrt(env[,10]) 
  
  #specify K (vertical eddy diffusivity) follow Yeates & Imberger 1994 eq 5
  K <- NULL
  
  for(i in 1:length(N)){
    if(N[i] == 0 & depths[i] < td){
      K[i] = 1e-2
    } else if(N[i] == 0 & depths[i] > td) {
      K[i] = 1e-6
      
    } else {
      K[i] = 0.0000000000049 * ((9.80665*281.0694) / N[i])  #Area in km2 is 0.079 and we convert to square meters take the square root of that
    }
  }
  
  K_prime = c(first(diff(K)),diff(K))
  
  K_df <- data.frame(depth = as.numeric(depths),
                     K = K,
                     K_prime = K_prime)
  K_df <- tibble(K_df)
  
  #specify model timestep
  del_t = 60 #in seconds
  
  # #calculate u_star (shear velocity)
  # td = thermo.depth(wtr = env[,1], depths = env[,2])
  # avgEpiDense = mean(env[which(env[,2] <= td),3])
  # u_star = uStar(wndSpeed = wnd, wndHeight = 3, averageEpiDense = avgEpiDense)
  # 
  #loop through and calculate individual movement
  for(j in 1:total_inds){
    
    # curr_dens <- env[which(env[,2] == round(inds[j,2],1)),3]
    # curr_visc <- env[which(env[,2] == round(inds[j,2],1)),4]
    
    w_s <- -0.13/1440 #this is what's in Cayelan/Kamilla's GLM-AED calibration for cyanobacteria at FCR converted to meters per minute instead of per day
    # w_s <- ((9.8081*cell_diam[j]^2*(cell_dens[j] - curr_dens))/(18*cell_shape[j]*curr_visc))/100000 # Define the cell velocity given cell diameter, density, and shape
    
    
    ####### Attempt at random walk ###############################################
    
    #pull individual's depth
    z <- round(unname(inds[j, yloc]),1)
    K_z <- K[which.min(depths - z)]
    K_prime_z <- K_prime[which.min(depths - z)]
    
    #specify random walk following Visser 1997 https://www.int-res.com/articles/meps/158/m158p275.pdf
    # between Ross & Sharples paper and Visser 1997 paper, a lot of confusion re: what z is!
    # is it height above sediments, or depth with a negative sign, or depth with no negative sign
    # elev <- ymax - z
    
    z_t1 <- z + K_prime_z*z*del_t + runif(1, min = -0.1,max = 0.1)*sqrt((2*K_z*(z + 0.5*K_prime_z*z*del_t)*del_t)/(1/3)) + w_s # w_s*del_t # this is needed if you use Stokes
    
    ####### end test code #######################################################
    inds[j, yloc] <- z_t1
    inds[j, delta_z] <- z_t1 - z
    
  }
  
  
  # =========   The reflecting boundary is added below
  for(i in 1:total_inds){ 
    if(inds[i, yloc] > ymax){         # If it moved past the maximum depth
      inds[i, yloc] <- ymax;        # Then move it back to the maximum depth
    }
    if(inds[i, yloc] < 0.1){            # If it moved below 0.1 (above surface)
      inds[i, yloc] <- 0.1;           # Then move it back to 0.1 (surface)
    }
  } 
  # =========  Now all individuals should stay on the landscape
  return(inds);
}

# =============================================================================
# Growth function
# =============================================================================

growth <- function(inds, repr_col = 7, traits = traits_lst, growth_env = env){
  total_inds       <- dim(inds)[1]; # Get the number of individuals in inds
  ind_cols         <- dim(inds)[2]; # Total inds columns
  
  #unpack traits
  diam = traits$diam
  dens = traits$dens
  shape = traits$shape
  # Tmin = traits$Tmin
  # Topt = traits$Topt
  # Tmax = traits$Tmax
  T_0 = traits$T_0
  q = traits$q
  I_S = traits$I_S
  umax = traits$umax
  N_0 = traits$N_0
  K_N = traits$K_N
  P_0 = traits$P_0
  K_P = traits$K_P
  
  ###Temp-dependent growth
  for(i in 1:total_inds){
    
    #get water temp and SWR where individual is located
    TEMP <- env[which(env[,2] == round(inds[i,2],1)),1]
    SWR <- env[which(env[,2] == round(inds[i,2],1)),5]
    DIN <- env[which(env[,2] == round(inds[i,2],1)),7]
    FRP <- env[which(env[,2] == round(inds[i,2],1)),8]
    

      #temp-dependent growth following Hellweger et al. 2008 https://aslopubs.onlinelibrary.wiley.com/doi/epdf/10.4319/lo.2008.53.4.1227?src=getftr
      # fT = ((TEMP - Tmin) / (Topt - Tmin)) *((Tmax - TEMP) / (Tmax - Topt)) ^((Tmax - Topt) / (Topt - Tmin))
      # if(fT < 0 | is.na(fT)){fT <- 0}
      fT = exp(-( (TEMP - T_0) / q )^2)
      
      #light-dependent growth with photoinhibition (Steele 1962) https://doi.org/10.4319/lo.1962.7.2.0137
      fI = (SWR/I_S) * exp(1 - (SWR/I_S))
      if(SWR < 5e-5 | fI < 5e-5){fI = 0.0}
      
      #nutrient dependent growth (Monod model) for DIN
      fN = (DIN - N_0) / (DIN - N_0 + K_N)
      
      #nutrient dependent growth (Monod model) for FRP
      fP = (FRP - P_0) / (FRP - P_0 + K_P)
      
    
    inds[i, repr_col] <- rbinom(n = 1, size = 1, prob = umax*min(c(fT,fI,fN,fP)))
    
  }
  
  #make array for offspring
  new_inds <- array(data = 0, dim = c(0, ind_cols))
  
  #figure out how many offspring of each PFT per depth
  for(j in 1:length(env[,2])){
    
    total_off <- sum(inds[which(inds[,1] == 1 & round(inds[,2],1) == env[j,2]), repr_col])

    # ---- We now have the total number of new offspring for each PFT at this depth; now add traits
    
    #generic PFT
    temp.df <- array(data = 0, dim = c(total_off, ind_cols))
    temp.df[,1] <- 1 #placeholder for taxon ID or some other trait
    temp.df[,2] <- env[j,2]
    temp.df[,3] <- diam
    temp.df[,4] <- dens
    temp.df[,5] <- shape
    
    new_inds <- rbind(new_inds, temp.df)
  }
  
  # ---- Our new offspring can now be attached in the inds array
  inds <- rbind(inds, new_inds);
  return(inds);
}

# =============================================================================
# Death function
# =============================================================================
death <- function(inds, dcol = 6, yloc = 2, ymax = 9.3, pft = 1, traits = traits_lst, light_exp = 8, death_env = env){
  total_inds <- dim(inds)[1] # Get the number of individuals in inds
  
  # unpack paramters
  R_resp = traits$R_resp
  theta_resp = traits$theta_resp
  
  for(i in 1:total_inds){ 
    
    # Conduct bernoulli draws for death based on grazing rates for individual pfts
    TEMP <- env[which(env[,2] == round(inds[i,2],1)),1]
    resp = R_resp*theta_resp^(TEMP-20)
    inds[i, dcol] <- rbinom(1, 1, resp)
    
    # #update light exposure; individuals with no light for 24 hrs or more die
    # SWR <- env[which(env[,2] == round(inds[i,2],1)),5]
    # if(SWR < 5e-5){inds[i, light_exp] <- inds[i, light_exp] + 1}
    # if(inds[i, light_exp] >= 24*60){inds[i, dcol] <- 1}
    
    #anything at sediments dies
    if(inds[i, yloc] >= ymax){# If it moved to or past the maximum depth
      inds[i, dcol] <- 1;        # Then it dies
    }
  }
  
  inds            <- inds[inds[, dcol] == 0,]; # Retain living
  
  return(inds);
}

# ========================================================
# Update environment function
# ========================================================
update_env <- function(env, wtemp, met, din, frp, time_steps = 60*2, tstep = ts, depths = lake_depths){
  
  ts = tstep
  
  update_times <- seq(from = 60, to = time_steps, by = 60)
  
  if(ts %in% update_times){
    
    #calculate which column of wtemp we want
    col <- 7+ts/60
    
    #update water temp column
    env[,1] <- unlist(wtemp[,col])
    
    #calculate which row of swr and wnd we want
    row <- 6+ts/60
    
    #update light column
    env[,5] <- exp(as.double(log(met[row,2])) - 0.5*depths)
    
    #update wind column
    env[,6] <- unlist(met[row,3])
    
    #update din column
    env[,7] <- unlist(din[,col])
    
    #update frp column
    env[,8] <- unlist(frp[,col])
    
  }
  
  return(env)
  
}

# ========================================================
# Initialize phytoplankton function
# ========================================================

initialize_phytos <- function(depths){
  # ----- Initialise individuals (phytos)
  
  # fp <- read_csv("./data/FluoroProbe_2021-08-02_FCR_50.csv")
  fp = tibble(Depth_inc = depths,
                  GreenAlgae_ugL = 10,
                  Bluegreens_ugL = 10,
                  BrownAlgae_ugL = 10,
                  TotalConcNoMixed_ugL = 10)
  
  # Create inds by looping through each depth increment and then creating that number of individuals at that depth based on ug/L
  inds <- array(data = 0, dim = c(0, 9))
  
  for(i in 1:length(depths)){
    
    #isolate a particular depth
    temp <- fp[i,]
    
    #create cell diameter trait
    diam = 20 # 20 um 
    
    #create cell density trait
    dens = 1050 # In range for chlorophytes
    
    #create cell shape trait
    shape = 1 # make them all spherical for now
    
    #populate depths with traits
    temp.df <- array(data = 0, dim = c(round(temp$TotalConcNoMixed_ugL,0),9)) # THIS IS WHERE YOU 
    # NEED TO EVENTUALLY ACCOUNT FOR LAYER VOLUME!! (MULTIPLY BY LITERS IN THAT LAYER)
    temp.df[,1] <- 1 #placeholder for taxon ID or some other trait
    temp.df[,2] <- depths[i]
    temp.df[,3] <- diam
    temp.df[,4] <- dens
    temp.df[,5] <- shape

    inds <- rbind(inds, temp.df)
  }
  
  colnames(inds) <- c("PFG","yloc","cell_diam","cell_dens","cell_shape","dcol","repr","light_exp","delta_z")
  
  return(inds)
}

# ========================================================
# Initialize environment function
# ========================================================

initialize_env <- function(depths, n_days){
  
  env <- array(data = 0, dim = c(length(depths),10))
  
  # first column is temperature
  wtemp <- read_csv("./data/cal_wtemp_GLM.csv") 
  wtemp_depths <- wtemp[,1]
  wtemp <- wtemp[,-1]
  wtemp <- do.call("cbind", replicate(n_days, wtemp, simplify = FALSE))
  wtemp <- cbind(depths, wtemp)
  colnames(wtemp) <- c("depth",seq(1:(ncol(wtemp)-1)))
  env[,1] <- unlist(wtemp[,7])
  
  # second column is depth
  env[,2] <- round(depths,1)
  
  met <- read_csv("./data/cal_met_GLM.csv")
  met <- do.call("rbind", replicate(n_days, met, simplify = FALSE)) 
  met <- met %>%
    mutate(time = seq(1:nrow(met)))
  
  # fifth column is light
  Kd = 0.4 # this is from Kamilla's PEST calibrated GLM
  env[,5] <- exp(as.double(log(met[6,2])) - Kd*depths) # get light attenuation with depth

  # sixth column is wind
  env[,6] <- unlist(met[6,3])
  
  # seventh column is din
  din <- read_csv("./data/cal_din_GLM.csv") 
  din_depths <- din[,1]
  din <- din[,-1]
  din <- do.call("cbind", replicate(n_days, din, simplify = FALSE))
  din <- cbind(depths, din)
  colnames(din) <- c("depth",seq(1:(ncol(din)-1)))
  env[,7] <- unlist(din[,7])
  
  # eighth column is frp
  frp <- read_csv("./data/cal_frp_GLM.csv") 
  frp_depths <- frp[,1]
  frp <- frp[,-1]
  frp <- do.call("cbind", replicate(n_days, frp, simplify = FALSE))
  frp <- cbind(depths, frp)
  colnames(frp) <- c("depth",seq(1:(ncol(frp)-1)))
  env[,8] <- unlist(frp[,7])
  
  colnames(env) <- c("wt","yloc","dens","visc","light","wnd","din","frp","delta_z","buoyancy_freq")
  
  return(list(env_init = env, wtemp = wtemp, met = met, din = din, frp = frp))
  
}
