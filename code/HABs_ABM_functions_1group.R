# Title: HABs AMB functions for model with 1 PFT
# Author: Mary Lofton
# Date: 16MAR24

# Purpose: Library of custom functions to run R-based HABs ABM

# =============================================================================
# Movement function
# =============================================================================
movement <- function(inds, env, yloc = 2, ymax = 9.3, wnd = 3){
  
  total_inds <- dim(inds)[1]
  
  cell_diam <- inds[,3] # Get the cell diameters of all individuals
  cell_dens <- inds[,4] # Get the cell densities of all individuals
  cell_shape <- inds[,5] # Get the cell shapes of all individuals
  
  #calculate water density
  env[,3] <- (999.83952 + 16.945176*env[,1] - 0.0079870401*env[,1]^2 - 0.000046170461*env[,1]^3 + 0.00000010556302*env[,1]^4 - 0.00000000028054253*env[,1]^5) / (1 + 0.016879850*env[,1]) #Kell equation for water density (Kell, 1975)
  
  #calculate water viscosity
  env[,4] <- (10^(((20 - env[,1])/(env[,1] + 96))*(1.2364 - 0.00137*(20 - env[,1]) + 0.0000057*(20 - env[,1])^2)))*1000 #this actually calculates the ratio of viscosity at target temp/viscosity at 20 degrees C, but viscosity at 20 degrees C is ~ 1
  
  #specify model timestep
  del_t = 60 #in seconds
  
  #calculate u_star (shear velocity)
  td = thermo.depth(wtr = env[,1], depths = env[,2])
  avgEpiDense = mean(env[which(env[,2] <= td),3])
  u_star = uStar(wndSpeed = wnd, wndHeight = 3, averageEpiDense = avgEpiDense)
  
  #loop through and calculate individual movement
  for(j in 1:total_inds){
    
    curr_dens <- env[which(env[,2] == round(inds[j,2],1)),3]
    curr_visc <- env[which(env[,2] == round(inds[j,2],1)),4]
    
    w_s <- 0.00009028 #this is what's in GLM-AED calibration converted to meters per minute instead of per day
    # w_s <- ((9.8081*cell_diam[j]^2*(cell_dens[j] - curr_dens))/(18*cell_shape[j]*curr_visc))/100000 # Define the cell velocity given cell diameter, density, and shape
    
    
    ####### Attempt at random walk ###############################################
    
    #pull individual's depth
    z <- inds[j, yloc]
    
    #specify K (vertical eddy diffusivity)
    K = 0.4*u_star*9.3*((9.3 - z)/9.3)*(1-((9.3 - z)/9.3)) + 0.00001
    
    #specify K_prime
    K_prime = 0.4 * u_star * 9.3 * ((9.3 - z)/9.3) * (1/9.3) - 0.4 * u_star * 
      9.3 * (1/9.3) * (1 - ((9.3 - z)/9.3))
    
    # ## calculation of function for K_prime ################
    # f = expression(0.4*u_star*9.3*((9.3 - z)/9.3)*(1-((9.3 - z)/9.3)) + 0.00001)
    # K_prime_eq = D(f, "z")
    # K_prime_eq
    # K_prime2_eq = D(K_prime_eq, "z")
    # K_prime2 = -(0.4 * u_star * 9.3 * (1/9.3) * (1/9.3) + 0.4 * u_star * 9.3 * 
    #   (1/9.3) * (1/9.3))
    # #######################################################
    
    #specify random walk
    elev <- ymax - z
    
    z_t1 <- elev + K_prime*elev*del_t + runif(1, min = -1,max = 1)*sqrt((2*K*(elev + 0.5*K_prime*elev*del_t)*del_t)/(1/3)) + w_s # w_s*del_t # this is needed if you use Stokes
    
    ####### end test code #######################################################
    inds[j, yloc] <- 9.3 - z_t1
    
  }
  
  
  # =========   The reflecting boundary is added below
  for(i in 1:total_inds){ 
    if(inds[i, yloc] > ymax){         # If it moved past the maximum depth
      inds[i, yloc] <- ymax;        # Then move it back to the maximum depth
    }
    if(inds[i, yloc] < 0.5){            # If it is close to top boundary
      inds[i, yloc] <- runif(1, min = 0.1, max = 0.5);           # Create random mixed layer
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
  Tmin = traits$Tmin
  Topt = traits$Topt
  Tmax = traits$Tmax
  I_S = traits$I_S
  
  ###Temp-dependent growth
  for(i in 1:total_inds){
    
    #get water temp and SWR where individual is located
    TEMP <- env[which(env[,2] == round(inds[i,2],1)),1]
    SWR <- env[which(env[,2] == round(inds[i,2],1)),5]
    

      #temp-dependent growth
      fT = ((TEMP - Tmin) / (Topt - Tmin)) *((Tmax - TEMP) / (Tmax - Topt)) ^((Tmax - Topt) / (Topt - Tmin))
      if(fT < 0 | is.na(fT)){fT <- 0}
      
      #light-dependent growth with photoinhibition
      fI = (SWR/I_S) * exp(1 - (SWR/I_S))
      if(SWR < 5e-5 | fI < 5e-5){fI = 0.0}
    
    inds[i, repr_col] <- rbinom(n = 1, size = 1, prob = fT*fI/60)
    
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
death <- function(inds, dcol = 6, yloc = 2, ymax = 9.3, pft = 1, graze = 0.35/60, light_exp = 8, death_env = env){
  total_inds <- dim(inds)[1] # Get the number of individuals in inds
  
  for(i in 1:total_inds){ 
    
    # Conduct bernoulli draws for death based on grazing rates for individual pfts
    inds[i, dcol] <- rbinom(1, 1, graze)
    
    #update light exposure; individuals with no light for 24 hrs or more die
    SWR <- env[which(env[,2] == round(inds[i,2],1)),5]
    if(SWR < 1){inds[i, light_exp] <- inds[i, light_exp] + 1}
    if(inds[i, light_exp] >= 24*60){inds[i, dcol] <- 1}
    
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
update_env <- function(env, wtemp, swr, time_steps = 60*2, tstep = ts, depths = lake_depths){
  
  ts = tstep
  
  update_times <- seq(from = 60, to = time_steps, by = 60)
  
  if(ts %in% update_times){
    
    #calculate which column of wtemp we want
    col <- 7+ts/60
    
    #update water temp column
    env[,1] <- unlist(wtemp[,col])
    
    #calculate which row of swr we want
    row <- 6+ts/60
    
    #update light column
    env[,5] <- exp(as.double(log(swr[row,2])) - 0.5*depths)
    
  }
  
  return(env)
  
}

# ========================================================
# Initialize phytoplankton function
# ========================================================

initialize_phytos <- function(depths){
  # ----- Initialise individuals (phytos)
  
  fp <- read_csv("./data/FluoroProbe_2021-08-02_FCR_50.csv")
  
  # Create inds by looping through each depth increment and then creating that number of individuals at that depth based on ug/L
  inds <- array(data = 0, dim = c(0, 8))
  
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
    temp.df <- array(data = 0, dim = c(round(temp$TotalConcNoMixed_ugL,0),8)) # THIS IS WHERE YOU 
    # NEED TO EVENTUALLY ACCOUNT FOR LAYER VOLUME!! (MULTIPLY BY LITERS IN THAT LAYER)
    temp.df[,1] <- 1 #placeholder for taxon ID or some other trait
    temp.df[,2] <- depths[i]
    temp.df[,3] <- diam
    temp.df[,4] <- dens
    temp.df[,5] <- shape
    
    inds <- rbind(inds, temp.df)
  }
  
  colnames(inds) <- c("PFG","yloc","cell_diam","cell_dens","cell_shape","dcol","repr","light_exp")
  
  return(inds)
}

# ========================================================
# Initialize environment function
# ========================================================

initialize_env <- function(depths){
  
  env <- array(data = 0, dim = c(length(depths),5))
  
  wtemp <- read_csv("./data/cal_wtemp_GLM.csv") 
  wtemp_depths <- wtemp[,1]
  wtemp <- wtemp[,-1]
  n=7
  wtemp <- do.call("cbind", replicate(n, wtemp, simplify = FALSE))
  wtemp <- cbind(depths, wtemp)
  colnames(wtemp) <- c("depth",seq(1:(ncol(wtemp)-1)))
  env[,1] <- unlist(wtemp[,7])
  
  # second column is depth
  env[,2] <- round(depths,1)
  
  # fifth column is light
  swr <- read_csv("./data/cal_met_GLM.csv")
  n=7
  swr <- do.call("rbind", replicate(n, swr, simplify = FALSE)) 
  swr <- swr %>%
    mutate(time = seq(1:nrow(swr)))
  env[,5] <- exp(as.double(log(swr[6,2])) - 0.5*depths)
  colnames(env) <- c("wt","yloc","dens","visc","light")
  
  return(list(env_init = env, wtemp = wtemp, swr = swr))
  
  }