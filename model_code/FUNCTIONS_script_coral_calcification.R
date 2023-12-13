
# ---
# title: "functions_coral_calcification"
# author: "Anne Kruijt"
# date: '2023-04-04'
# ---


#--------------FUNCTIONS USED IN THE CORAL CALCIFICATION MODEL-------------------


#-------------------- Irradiance calculations ---------------------------------

Irradiance_func <- function(lats, days, G_daily_lats, omega_s, declination, daylength) {
  
  for( i in 1:length(lats)){
    
    lat <- lats[i]
    
    for (j in days){   # j is the Julian day number
      
      # ----- declination angle in radians
      d <- 23.45 * sin( 2 * pi * (284 + j)/ 365) * deg_to_rad # simple formula by Cooper (1969), gives result in degrees
      
      declination[j,i] <- d * (1/deg_to_rad) # declination in radians
      
      # ----- sun hour angle in radians ----   following Berger (1978) and Duffie and Beckman (2013)
      if ( abs(lat -d/deg_to_rad ) > 90){
        omega_s[j,i] <- 0 * deg_to_rad
      }
      else if (abs(lat + d/deg_to_rad ) > 90){
        omega_s[j,i] <- 180 * deg_to_rad
      }
      else {
        omega_s[j,i] <- acos(-tan(deg_to_rad*lat)*tan( d )) 
      }  
      
      ##--- day length calculation----  following Duffie and Beckman (2013)
      
      daylength[j,i] <- 2 * omega_s[j,i] * (1/ deg_to_rad) / 15    
      
      ##------- Daily irradiance, method following Berger (1978), Duffie and Beckman (2013) 
      
      Ecorr <- 1 + 0.033 * cos( (2 * pi * j) / 365 ) # correction factor: Duffie and Beckman (2013)
      
      G_daily_lats[j,i] <- Isc * Ecorr * (omega_s[j,i]*sin(d)*sin(deg_to_rad * lat) + cos(d)*cos(deg_to_rad*lat)*sin(omega_s[j,i]))
      
    }
    
  }
  return(list(omega_s = omega_s, daylength = daylength, G_daily_lats = G_daily_lats, declination = declination))
  
}

#-------------------- Light refraction function --------------------------------

Reflection <- function(latitude, declination) {
  na <- 1 #
  nw <- 1.335 # incidence of refraction water (Shifrin)
  Phi_a <- latitude - abs(declination) #angle of light beam with the normal at 12h
  Phi_w <- asin( na * sin(Phi_a * deg_to_rad) / nw ) * (1/deg_to_rad) # angle of refracted beam in the water with normal, calculated using Snell's law
  
  # Light reflection, with Fresnel equations
  Rperp <- (sin((Phi_a - Phi_w)*deg_to_rad))^2 / (sin((Phi_a + Phi_w) * deg_to_rad))^2
  Rparalell <- (tan((Phi_a - Phi_w)* deg_to_rad))^2 / (tan((Phi_a + Phi_w)* deg_to_rad))^2
  Reff <- 0.5*(Rperp + Rparalell) # effective fraction of reflection
  
  
  return(Reff)
}


#--------------------- Functions to compute the DEATH of a coral due to darkness--------------------------
DARKDAYS <- function(daylength, irradiance, threshold_setting){
  threshold =  threshold_setting
  daylength[irradiance < threshold] <- 0   # minimum average irradiance on a day that is needed for calcification. Below this, coral is only respiring (Gattuso, ...)
  darkdays <- colCounts(daylength, value = 0)  # number of days in a year a coral is below the threshold value, for each latitude 
  return(darkdays)
}

DEATH <- function(darkdays, latitude_index, maxdays){   
  
  if (darkdays[latitude_index] > maxdays){  # if the amount of dark days at a certain latitude is more than the max amount of dark days a coral can handle, it will die (value =0)
    val <- 0 }
  else {
    val <- 1}
  return(val)
}

#-------- Light penetration into water----------------

LambertBeer <- function(E0, z, kpar){
  E <- E0 * exp(-kpar *z)
  return(E)}

f_light <- function(E, Emax, Ek){
  f <- tanh(E/Ek)   # tangent curve
  return(f)}

#--------------- temperature dependence curve------------------------


f_temp_block <- function(Temp_cold, Temp_warm,Tboundary_low, Tboundary_high){
  
  # The combined curve
  
  if (Temp_cold < Tboundary_low){
    f_fraction <- 0
  }
  else if (Temp_warm > (Tboundary_high)){
    f_fraction <- 0
  }
  else{
    f_fraction <- 1
  }
  
  #from rate to fraction:
  return(f_fraction)
}


#---------------- calcification  --------------------#
#Note: in the final comparison of scenarios we scale all calcification rates to the maximum occurring in the scenario, so decreasing
# rates become a percentage of the maximum (100%). Since we assume omega to be constant and unlimiting, in that sense the actual 
# calcification rate coming out of the function CalcLinearReef does not matter. We however keep this function here in case we later 
# would want to study changes under different saturationc states (not the purpose of this study).

# Linear relationship, based on experimental and observational studies: at omega =1 calcification is still positive
CalcLinearReef <- function(Omega){
  G <- (15.5*Omega + 41.69) /(24*60*60) #millimol/m2/s at constant light and temp, following from Eyre 2018
  return(G)}


# Calcification function based on the linear relation with Omega
calcfunc <- function(Omega, lightdep, tempdep){
  # calcification rate
  calc_reef <-  CalcLinearReef(Omega) * lightdep * tempdep # mmol/m2/s (lightdep and tempdep are determined per day)
  calc_reef <- pmax(calc_reef, 0)  # making sure there is no negative calcification
  calc_reef[is.nan(calc_reef) == TRUE] <- 0
  calc_reef_year <- colSums(calc_reef * 60 * 60 *24, dim = 1) 
  # converting from /s to /day and then summing value for each day to obtain yearly total,
  # we realize that the instantaneous rate is not constant the whole day long, but since the light and tempdep are the only
  # limiting factors we consider in this study, and they are calculated per day, this simplification does not affect our conclusions
  return(calc_reef_year)
}


# ------ purely mathematical functions ----- #

interpolating_temperatures <- function(start_lat, final_lat, length_original, temperatures_original, lat_step_new){
  
  #Given latitude values
  original_latitudes <- seq(start_lat, final_lat, length.out = length_original)
  
  # New latitude values
  new_latitudes <- seq(start_lat, final_lat, by = lat_step_new)
  
  # Linear interpolation to get temperature values at new latitudes
  new_temperatures <- interp1(original_latitudes, temperatures_original, new_latitudes)
  
  return(new_temperatures)
}

interpolating_temperatures_month_to_year <- function(start_time, final_time, length_original, temperatures_original, time_step_new){
  
  #Given time
  original_time <- seq(start_time, final_time, length.out = length_original)
  
  # New time values
  new_times <- seq(start_time, final_time, by = time_step_new)
  
  # Linear interpolation to get temperature values at new times
  new_temperatures <- interp1(original_time, temperatures_original, new_times)
  
  return(new_temperatures)
}


# function for scaling data to 0-100% 
scaling <- function(data){
  datanew <- data/max(data) *100
  return(datanew)
}