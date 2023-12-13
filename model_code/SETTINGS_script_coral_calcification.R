# ---
# title: "Settings_coral_calcification"
# author: "Anne Kruijt"
# date: '2023-04-04'
# ---

#---------- DEFINING THE MODEL SETTINGS----------------------#


##------------ Domain of the model simulation ---------------------
boxdepth  <- 15 # unit= meters
lats <- seq(0, 90, 1)     # 90 degrees on Northern hemisphere
days <- seq(1, 365, 1)  # One year
# indices for lengths of rows and columns
latlen <- length(lats)
timelen <- length(days)


##------------ Constants for irradiance calculations -----------------
Isc <- 1.367 # Solar constant, kWatt/m2 
T <- 24 # hours in a day
deg_to_rad <- pi/180 # conversion from degrees to radians

##------------ Constants for the light regime ------------------------
Emax <-  625 # umol/m2/s # maximum light at earth surface 
kpar_standard <- 0.05 # Very oligotrophic systems have even lower kpars (close to 0)

## ----------- Prescribed settings for calcification and temperature tolerance -------------
OmegaArag <- 3 # Value for present study doesnt really matter since we assume omega to be unlimiting everywhere
Topt_high <- 26
Tboundary_high <- 50 # set to 36 for normal upper boundary and 50 or higher in case of studying the effect of lower temperature boundary only
Tboundary_low <- 16



