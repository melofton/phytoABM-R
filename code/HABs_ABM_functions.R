# Title: HABs AMB functions
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
    
    w_s <- ((9.8081*cell_diam[j]^2*(cell_dens[j] - curr_dens))/(18*cell_shape[j]*curr_visc))/100000 # Define the cell velocity given cell diameter, density, and shape
    
    
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
    
    z_t1 <- elev + K_prime*elev*del_t + runif(1, min = -1,max = 1)*sqrt((2*K*(elev + 0.5*K_prime*elev*del_t)*del_t)/(1/3)) + w_s*del_t
    
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
## define traits for PFTs
traits_lst <- list(
  #create cell diameter trait
  min_diam.d = 10, # 10 um 
  max_diam.d = 25, # 25 um in meters
  
  min_diam.c = 10, # 10 um 
  max_diam.c = 25, # 25 um in meters
  
  min_diam.g = 10, # 10 um 
  max_diam.g = 25, # 25 um in meters
  
  #create cell density trait
  min_dens.d = 1070, # In range for diatoms
  max_dens.d = 1130, # To make sure they all sink
  
  min_dens.c = 920, # In range for cyanos
  max_dens.c = 980, # To make sure they all float
  
  min_dens.g = 1020, # In range for chlorophytes
  max_dens.g = 1080, # To make sure they all sink
  
  #create cell shape trait
  shape.d = 1, # make them all spherical for now
  
  shape.c = 1, # make them all spherical for now
  
  shape.g = 1, # make them all spherical for now
  
  #create temperature sensitivity traits
  Tmin = c(4,4,4),
  Topt = c(23,23,23),
  Tmax = c(30,30,30),
  I_S = c(300,300,300)
)

birth <- function(inds, repr_col = 7, traits = traits_lst){
  total_inds       <- dim(inds)[1]; # Get the number of individuals in inds
  ind_cols         <- dim(inds)[2]; # Total inds columns
  
  #unpack traits
  min_diam.d = traits$min_diam.d
  max_diam.d = traits$max_diam.d
  min_diam.c = traits$min_diam.c
  max_diam.c = traits$max_diam.c
  min_diam.g = traits$min_diam.g
  max_diam.g = traits$max_diam.g
  min_dens.d = traits$min_dens.d
  max_dens.d = traits$max_dens.d
  min_dens.c = traits$min_dens.c
  max_dens.c = traits$max_dens.c
  min_dens.g = traits$min_dens.g
  max_dens.g = traits$max_dens.g
  shape.d = traits$shape.d
  shape.c = traits$shape.c
  shape.g = traits$shape.g
  Tmin = traits$Tmin
  Topt = traits$Topt
  Tmax = traits$Tmax
  I_S = traits$I_S
  
  ###Temp-dependent growth
  for(i in 1:total_inds){
    
    #get water temp and SWR where individual is located
    TEMP <- env[which(env[,2] == round(inds[i,2],1)),1]
    SWR <- env[which(env[,2] == round(inds[i,2],1)),5]
    
    if(inds[i,1] == 1){ #diatoms
      
      #temp-dependent growth
      fT = ((TEMP - Tmin[1]) / (Topt[1] - Tmin[1])) *((Tmax[1] - TEMP) / (Tmax[1] - Topt[1])) ^((Tmax[1] - Topt[1]) / (Topt[1] - Tmin[1]))
      if(fT < 0 | is.na(fT)){fT <- 0}
      
      #light-dependent growth with photoinhibition
      fI = (SWR/I_S[1]) * exp(1 - (SWR/I_S[1]))
      if(SWR < 5e-5 | fI < 5e-5){fI = 0.0}
      
    }
    
    if(inds[i,1] == 2){ #cyanos
      
      #temp-dependent growth
      fT = ((TEMP - Tmin[2]) / (Topt[2] - Tmin[2])) *((Tmax[2] - TEMP) / (Tmax[2] - Topt[2])) ^((Tmax[2] - Topt[2]) / (Topt[2] - Tmin[2]))
      if(fT < 0 | is.na(fT)){fT <- 0}
      
      #light-dependent growth
      fI = (SWR/I_S[2]) * exp(1 - (SWR/I_S[2]))
      if(SWR < 5e-5 | fI < 5e-5){fI = 0.0}
      
    }
    
    if(inds[i,1] == 3){ #green algae
      
      #temp-dependent growth
      fT = ((TEMP - Tmin[3]) / (Topt[3] - Tmin[3])) *((Tmax[3] - TEMP) / (Tmax[3] - Topt[3])) ^((Tmax[3] - Topt[3]) / (Topt[3] - Tmin[3]))
      if(fT < 0 | is.na(fT)){fT <- 0}
      
      #light-dependent growth
      fI = (SWR/I_S[3]) * exp(1 - (SWR/I_S[3]))
      if(SWR < 5e-5 | fI < 5e-5){fI = 0.0}
      
    }
    
    inds[i, repr_col] <- rbinom(n = 1, size = 1, prob = fT*fI/60)
    
  }
  
  #make array for offspring
  new_inds <- array(data = 0, dim = c(0, ind_cols))
  
  #figure out how many offspring of each PFT per depth
  for(j in 1:length(env[,2])){
    
    total_off_d <- sum(inds[which(inds[,1] == 1 & round(inds[,2],1) == env[j,2]), repr_col])
    total_off_c <- sum(inds[which(inds[,1] == 2 & round(inds[,2],1) == env[j,2]), repr_col])
    total_off_g <- sum(inds[which(inds[,1] == 3 & round(inds[,2],1) == env[j,2]), repr_col])
    
    # ---- We now have the total number of new offspring for each PFT at this depth; now add traits
    
    #diatoms
    temp.df.d <- array(data = 0, dim = c(total_off_d, ind_cols))
    temp.df.d[,1] <- 1 #placeholder for taxon ID or some other trait
    temp.df.d[,2] <- env[j,2]
    temp.df.d[,3] <- runif(dim(temp.df.d)[1], min = min_diam.d, max = max_diam.d)
    temp.df.d[,4] <- runif(dim(temp.df.d)[1], min = min_dens.d, max = max_dens.d)
    temp.df.d[,5] <- shape.d
    
    #cyanos
    temp.df.c <- array(data = 0, dim = c(total_off_c, ind_cols))
    temp.df.c[,1] <- 2 #placeholder for taxon ID or some other trait
    temp.df.c[,2] <- env[j,2]
    temp.df.c[,3] <- runif(dim(temp.df.c)[1], min = min_diam.c, max = max_diam.c)
    temp.df.c[,4] <- runif(dim(temp.df.c)[1], min = min_dens.c, max = max_dens.c)
    temp.df.c[,5] <- shape.c
    
    #green algae
    temp.df.g <- array(data = 0, dim = c(total_off_g, ind_cols))
    temp.df.g[,1] <- 3 #placeholder for taxon ID or some other trait
    temp.df.g[,2] <- env[j,2]
    temp.df.g[,3] <- runif(dim(temp.df.g)[1], min = min_diam.g, max = max_diam.g)
    temp.df.g[,4] <- runif(dim(temp.df.g)[1], min = min_dens.g, max = max_dens.g)
    temp.df.g[,5] <- shape.g
    
    temp.df <- rbind(temp.df.d, temp.df.c, temp.df.g)
    new_inds <- rbind(new_inds, temp.df)
  }
  
  # ---- Our new offspring can now be attached in the inds array
  inds <- rbind(inds, new_inds);
  return(inds);
}

# =============================================================================
# Death function
# =============================================================================
death <- function(inds, dcol = 6, yloc = 2, ymax = 9.3, pft = 1, graze_d = 0.35/60, graze_c = 0.35/60, graze_g = 0.35/60, light_exp = 8){
  total_inds <- dim(inds)[1] # Get the number of individuals in inds
  
  for(i in 1:total_inds){ 
    
    # Conduct bernoulli draws for death based on grazing rates for individual pfts
    if(inds[i, pft] == 1){ 
      inds[i, dcol] <- rbinom(1, 1, graze_d)
    }
    if(inds[i, pft] == 2){
      inds[i, dcol] <- rbinom(1, 1, graze_c)
    }
    if(inds[i, pft] == 3){
      inds[i, dcol] <- rbinom(1, 1, graze_g)
    }
    
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
update_env <- function(env, wtemp, swr, time_steps = 60*2, tstep = ts){
  
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