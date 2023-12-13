# ---
# title: "Temperature scenario creation"
# author: "Anne Kruijt"
# date: '2023-11-6'
# ---


## SCRIPT CONTAINING CODE TO CREATE THE DIFFERENT TEMPERATURE SCENARIOS ##
#-----------------------------------------------------------------------#
#  To be used in the master coral calcification script

file_location <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(file_location)

library(here)

source("FUNCTIONS_script_coral_calcification.R")


# --------- reading in CESM files
# reference:	Zhu, J., Poulsen, C. J., & Tierney, J. E. (2019). 
# Simulation of Eocene extreme warmth and high climate sensitivity through cloud feedbacks. 
# Science Advances, 5(9). https://doi.org/10.1126/sciadv.aax1874

CESM_temperature_data3x_monthly <- read.csv("../data_files/SST_monthly_3x.csv") 
CESM_temperature_data6x_monthly <- read.csv("../data_files/SST_monthly_6x.csv")
CESM_temperature_data9x_monthly <- read.csv("../data_files/SST_monthly_9x.csv")

CESM_latitudes <- seq(90,-90, length.out = length(CESM_temperature_data3x_monthly[,2])) # resolution of CESM is 2.5 degrees


#---------- reading in MODERN TEMPERATURE files:----------------------#
# reference: 22.	Locarnini, R. A. , Mishonov, A. V., Baranova, O. K., Boyer, T. P., Zweng, M. M., García, H., Reagan, J., Seidov, D., Weathers, K., Paver, C., & Smolyar, I. (2019).
# World Ocean Atlas 2018, Volume 1: Temperature. (Vol. 1). NOAA Atlas NESDIS.

ModernSST_monthly <- read.csv('../data_files/modern_SST_data_monthly.csv')
ModernSST_monthly[,1] <- seq(-90,90, length.out = length(ModernSST_monthly[,1]))


# ---- MONTHLY TO DAILY TEMPERATURE VALUES ------#

df_temperaturesCESM_3xM <- data.frame(CESM_latitudes, CESM_temperature_data3x_monthly[,3:14]) # creating df of relevant SST data together with latitude
df_temperaturesCESM_6xM <- data.frame(CESM_latitudes, CESM_temperature_data6x_monthly[,3:14])
df_temperaturesCESM_9xM <- data.frame(CESM_latitudes, CESM_temperature_data9x_monthly[,3:14])


temp3x_interpolated_M <- array(data = NaN, dim = c(181,13)) # empty array to be filled with values in for-loop
temp6x_interpolated_M <- array(data = NaN, dim = c(181,13))
temp9x_interpolated_M <- array(data = NaN, dim = c(181,13))
ModernSST_M <- array(data = NaN, dim = c(181,13))
temp3x_interpolated_D <- array(data = NaN, dim = c(181,365))
temp6x_interpolated_D <- array(data = NaN, dim = c(181,365))
temp9x_interpolated_D <- array(data = NaN, dim = c(181,365))
ModernSST_D <- array(data = NaN, dim = c(181,365))

#interpolating between latitudes
for (i in 1:13){
  Temp3xM <- rev(df_temperaturesCESM_3xM[,i]) # reversing order to make data start with the value at 90degS
  temp3x_interpolated_M[,i]<- interpolating_temperatures(-90, 90, length(Temp3xM), Temp3xM, 1)
  
  Temp6xM <- rev(df_temperaturesCESM_6xM[,i])
  temp6x_interpolated_M[,i]<- interpolating_temperatures(-90, 90, length(Temp6xM), Temp6xM, 1)
  
  Temp9xM <- rev(df_temperaturesCESM_9xM[,i])
  temp9x_interpolated_M[,i]<- interpolating_temperatures(-90, 90, length(Temp9xM), Temp9xM, 1)
  
  ModernSST_M[,i] <-  interpolating_temperatures(-90, 90, length(as.numeric(ModernSST_monthly[,i])), as.numeric(ModernSST_monthly[,i]), 1)
}

#interpolating between months 
for (i in 11:181){ # skipping first 10 rows, these contain NaNs
  temp3x_interpolated_D[i,]<- interpolating_temperatures_month_to_year(1, 365, length(temp3x_interpolated_M[i,2:13]), temp3x_interpolated_M[i,2:13], 1) # from 2:13, first column contains the latitudes
  temp6x_interpolated_D[i,]<- interpolating_temperatures_month_to_year(1, 365, length(temp6x_interpolated_M[i,2:13]), temp6x_interpolated_M[i,2:13], 1)
  temp9x_interpolated_D[i,]<- interpolating_temperatures_month_to_year(1, 365, length(temp9x_interpolated_M[i,2:13]), temp9x_interpolated_M[i,2:13], 1)
}

for (i in 14:181){ # skipping first 13 rows, these contain NaNs
  ModernSST_D[i,]<- interpolating_temperatures_month_to_year(1, 365, length(ModernSST_M[i,2:13]), ModernSST_M[i,2:13], 1) 
  
}

# 0:90 degrees N (extracting data for half sphere (focus on one hemishpere in the paper, better coverage in the geological record)
temp3x_interpolated_D_NH <- temp3x_interpolated_D[90:181,]
temp6x_interpolated_D_NH <- temp6x_interpolated_D[90:181,]
temp9x_interpolated_D_NH <- temp9x_interpolated_D[90:181,]
ModernSST_D_NH <- ModernSST_D[90:181,]
ModernSST_M_NH <- data.matrix(ModernSST_M[90:181,2:13])

# ------------ YEARLY coldest DAY  ------------ #

# empty arrays, to be filled with minima. 
Tempx3Min <- array(, dim = c(181))
Tempx6Min <- array(, dim = c(181))
Tempx9Min <- array(, dim = c(181))
TempModernMin <- array(, dim = c(181))

Tempx3Min[11:181] <- rowMins(temp3x_interpolated_D[11:181,2:365])
Tempx6Min[11:181] <- rowMins(temp6x_interpolated_D[11:181,2:365])
Tempx9Min[11:181] <- rowMins(temp9x_interpolated_D[11:181,2:365])
TempModernMin[13:181] <- rowMins(ModernSST_D[13:181,2:365]) 

# at 0-90 degrees N
Tempx3Min_NH <- Tempx3Min[90:181]
Tempx6Min_NH<- Tempx6Min[90:181]
Tempx9Min_NH <- Tempx9Min[90:181]
TempModernMin_NH <- TempModernMin[90:181]

#--------  YEARLY coldest MONTH  -------------------------- #
Tempx3MinM <- array(, dim = c(181))
Tempx6MinM <- array(, dim = c(181))
Tempx9MinM <- array(, dim = c(181))
TempModernMinM <- array(, dim = c(181))

Tempx3MinM[11:181] <- rowMins(temp3x_interpolated_M[11:181,2:13])
Tempx6MinM[11:181] <- rowMins(temp6x_interpolated_M[11:181,2:13])
Tempx9MinM[11:181] <- rowMins(temp9x_interpolated_M[11:181,2:13])
TempModernMinM[13:181] <- rowMins(ModernSST_M[13:181,2:13])

# at 0-90 degrees NH
Tempx3MinM_NH <- Tempx3MinM[90:181]
Tempx6MinM_NH<- Tempx6MinM[90:181]
Tempx9MinM_NH <- Tempx9MinM[90:181]
TempModernMinM_NH <- TempModernMinM[90:181]

#---- YEARLY warmest MONTH -------------#
Tempx3MaxM <- array(, dim = c(181))
Tempx6MaxM <- array(, dim = c(181))
Tempx9MaxM <- array(, dim = c(181))
TempModernMaxM <- array(, dim = c(181))

Tempx3MaxM[11:181] <- rowMaxs(temp3x_interpolated_M[11:181,2:13])
Tempx6MaxM[11:181] <- rowMaxs(temp6x_interpolated_M[11:181,2:13])
Tempx9MaxM[11:181] <- rowMaxs(temp9x_interpolated_M[11:181,2:13])
TempModernMaxM[13:181] <- rowMaxs(ModernSST_M[13:181,2:13])

# at 0-90 degrees NH
Tempx3MaxM_NH <- Tempx3MaxM[90:181]
Tempx6MaxM_NH<- Tempx6MaxM[90:181]
Tempx9MaxM_NH <- Tempx9MaxM[90:181]
TempModernMaxM_NH <- TempModernMaxM[90:181]



# ------ YEARLY AVERAGE WINTER TEMPERARTURE (3 MONTH AVERAGE) ---------#
Tempx3winter <- array(, dim = c(181))
Tempx6winter <- array(, dim = c(181))
Tempx9winter <- array(, dim = c(181))
TempModern_winter  <- array(, dim = c(181))


Tempx3winter[11:89] <- rowMeans(subset((temp3x_interpolated_M[11:89,]), select = c(9,10,11))) # coldest SST months SH (August, September, Oktober)
Tempx3winter[90:181] <- rowMeans(subset((temp3x_interpolated_M[90:181,]), select = c(3,4,5))) # coldest SST months NH (February, March, April)
Tempx6winter[11:89] <- rowMeans(subset((temp6x_interpolated_M[11:89,]), select = c(9,10,11)))
Tempx6winter[90:181] <- rowMeans(subset((temp6x_interpolated_M[90:181,]), select = c(3,4,5)))
Tempx9winter[11:89] <- rowMeans(subset((temp9x_interpolated_M[11:89,]), select = c(9,10,11)))
Tempx9winter[90:181] <- rowMeans(subset((temp9x_interpolated_M[90:181,]), select = c(3,4,5)))

TempModern_winter[13:89] <- rowMeans(subset((ModernSST_M[13:89,]), select = c(9,10,11)))
TempModern_winter[90:181] <- rowMeans(subset((ModernSST_M[90:181,]), select = c(3,4,5)))

# 0-90 degrees
Tempx3winter_NH <- Tempx3winter[90:181]
Tempx6winter_NH<- Tempx6winter[90:181]
Tempx9winter_NH <- Tempx9winter[90:181]
TempModern_winter_NH <- TempModern_winter[90:181]

#------ writing to file (only those files used for plotting in figure 1 of the paper) --------#

write.csv(Tempx3MinM, paste("../model_output/", "Tempx3MinM", ".csv", sep = ""))
write.csv(Tempx6MinM, paste("../model_output/", "Tempx6MinM", ".csv", sep = ""))
write.csv(Tempx9MinM, paste("../model_output/", "Tempx9MinM", ".csv", sep = ""))
write.csv(TempModernMinM, paste("../model_output/", "TempModernMinM", ".csv", sep = ""))


write.csv(Tempx3MaxM, paste("../model_output/", "Tempx3MaxM", ".csv", sep = ""))
write.csv(Tempx6MaxM, paste("../model_output/", "Tempx6MaxM", ".csv", sep = ""))
write.csv(Tempx9MaxM, paste("../model_output/", "Tempx9MaxM", ".csv", sep = ""))
write.csv(TempModernMaxM, paste("../model_output/", "TempModernMaxM", ".csv", sep = ""))