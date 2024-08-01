# ---
# title: "Master function calcification"
# author: "Anne Kruijt"
# date: '2023-04-04'
# ---


## SCRIPT CONTAINING FUNCTION FOR COMPUTING CORAL CALCIFICATION ##
#-----------------------------------------------------------------------#

Model_simulation <- function(darkday_tolerance_ms, 
                             kpar_ms,
                             Ek_ms,
                             threshold_setting_ms,
                             Tcurve_ms,
                             Tscenario_ms1,
                             Tscenario_ms2,
                             day = days, lat = lats,
                             boxdepth_ms = boxdepth, 
                             Emax_ms = Emax, 
                             Tset_low_ms = Tset_low, Tset_high_ms = Tset_high,
                             width_low_ms = width_low, width_high_ms = width_high) {
  
  timelen <- length(days)
  latlen <-  length(lats)
  
  ##-------------START OF CALCULATIONS------------------------
  #-----------------------------------------------------------#
  
  ##------------ IRRADIANCE CALCULATIONS-------------------
  
  # Empty matrices for storing irradiance data
  
  G_daily_lats <- matrix(, nrow=timelen, ncol=latlen) # empty matrix for storing computed daily irradiance 
  omega_s <- matrix(, nrow=timelen, ncol=latlen) # empty matrix for storing computed sun hour angles
  declination <- matrix(, nrow=timelen, ncol=latlen)  # empty matrix for storing declination angles
  daylength <- matrix(, nrow=timelen, ncol=latlen)  # empty matrix

  
  # Calculating the irradiance data

  Irradiance_data <- Irradiance_func(lats, days, G_daily_lats, omega_s, declination, daylength)
  omega_s <- Irradiance_data$omega_s
  daylength <- Irradiance_data$daylength
  G_daily_lats <- Irradiance_data$G_daily_lats
  declination <- Irradiance_data$declination

  #------------ Refraction calculation--------------------
  
  R_fraction <- matrix(, nrow=timelen, ncol=latlen)
  
  for (i in 1:length(lats)){
    lat <- lats[i]
    
    for(j in days){
      dec <- declination[j,i] # in degrees
      R_fraction[j,i] <- Reflection(lat, dec)
    }
  }
  
 
  #Conversion to umol/m2/s, by scaling with prescribed Emax 
  Light1 <- as.data.frame(array(dim= c(timelen,latlen ), data = G_daily_lats ))
  
  factor <- 1/ max(Light1)
  Light1_scaled <- Light1 * factor  # scaling the irradiance curve, setting the maximum irradiance to 1
  
  Light2_scaled <- Light1_scaled * Emax 
  
  #------------ Correcting irradiance for reflection from surface 
  Light_corrected <- Light2_scaled - Light2_scaled* R_fraction
  

  ##----------- LIGHT LIMITATION CALCULATION----------------
  
  #----------- Creating dataframe to store the calculated value for light limitation
  
  column.names <- lats
  row.names <- days
  light_limitation <- array(,dim = c(timelen, latlen), dimnames = list(row.names,column.names))
  light_not_limiting <- array(1 ,dim = c(timelen, latlen), dimnames = list(row.names,column.names))
  
  #----------- Computing at which latitudes corals die due to excess dark days -----
  
  darkdays <- DARKDAYS(daylength, Light_corrected, threshold_setting_ms)
  
  
  #---- Computing light limitation for each latitude and day of the year
  
  for (i in 1:length(lats)) {
    
    death <- DEATH(darkdays, i, darkday_tolerance_ms)
    
    for (j in days){
      
      E0 <- Light_corrected[j,i] # surface light
      E <- LambertBeer(E0, boxdepth_ms, kpar_ms) # bottom light, mol/m2/day
      light_limitation[j,i] <- f_light(E, Emax_ms, Ek_ms) * death  # limitation by light for each latitude and day of year
    }
  }
  
  ##------- TEMPERATURE LIMITATION CALCULATIONS ----------##
  
  # creating empty matrix
  Temp_limitation <- matrix(, nrow = timelen, ncol =latlen)
  
  f_temp_function <- Tcurve_ms
  
  
  # in the case of daily varying temperatures, we use the following:
  # for (i in 1:length(lats)) {
  #   for (j in 1:timelen){
  #     Temp_limitation[j,i] <- f_temp_block(Tscenario_ms[i,j], Tset_low_ms, Tset_high_ms, width_low_ms, width_high_ms, Tboundary_low, Tboundary_high)
  #   }
  # 
  # }
  
  # calculating the temperature dependence for each latitude
  for (i in 1:length(lats)) {
    Temp_limitation[,i] <- f_temp_function(Tscenario_ms1[i], Tscenario_ms2[i],Tboundary_low, Tboundary_high)
    }
  
  #-----------  CALCIFICATION CALCULATION -----------------
  
  #### In case of inspecting the influence of actual meaningful calcification rates (based on known relationships with Omega)
  # # average calcification rate for each day in 
  # calc_reef <-  CalcLinearReef(OmegaArag) * Temp_limitation * light_limitation
  # calc_reef <- pmax(calc_reef, 0)  # making sure there is no negative calcification
  # calc_reef[is.nan(calc_reef) == TRUE] <- 0
  # 
  # # total calcification per year
  # calc_reef_year <- calcfunc(OmegaArag, light_limitation, Temp_limitation)
  # calc_reef_year_nolightlim <- calcfunc(OmegaArag, light_not_limiting, Temp_limitation)
  
  ##### In case of inspecting calcification potential solely depending on temperature and light limitation
  calc_reef <-  1 * Temp_limitation * light_not_limiting#light_limitation # calcification potential per day
  calc_reef <- pmax(calc_reef, 0)  # making sure there is no negative calcification
  calc_reef[is.nan(calc_reef) == TRUE] <- 0
  
  # total calcification potential per year
  calc_reef_year <- colSums(calc_reef, dim =1)
  calc_reef_year_nolightlim <- colSums(1 * Temp_limitation, dim =1 )
  
  #---------- END OF CALCULATIONS--------------------------------
  --------------------------------------------------------------

   return(list(Irradiance_data = Irradiance_data, light_limitation = light_limitation, calc_reef = calc_reef,
               calc_reef_year =calc_reef_year, calc_reef_year_nolightlim = calc_reef_year_nolightlim, darkdays = darkdays, Esurf = Light_corrected, templim =Temp_limitation))
  
}
