##---- SCRIPT TO EXTRACT OCEAN TEMPERATURE DATA FROM NOAA ATLAS ---- ##

# Following the step described on github page of Martin Schobben: https://github.com/MartinSchobben/oceanexplorer


# ------------------------------------------------------------------- #

# Install oceanexplorer from CRAN: 
install.packages("oceanexplorer")


library(oceanexplorer)

# obtain the NOAA world ocean atlas for temperature
temp_global <- get_NOAA('temperature', 1, "annual" )

# filter for a certain depth:

depth_layer <- 10
temp_depthlayer <- filter_NOAA(temp_global, depth = depth_layer)

## extract from NOAA framework
temp_extracted <- temp_depthlayer$t_an

## calculate zonal average
zonal_av <- colMeans(temp_extracted, na.rm = TRUE)  # removing Nans

## exporting temperature data to csv file

filelocation <- "./model_output/"

# write.csv(zonal_av,file= file.path(file_location, 'modern_SST_data.csv'))


## plotting temperature data

plot(zonal_av)

plot_NOAA(temp_depthlayer, depth = NULL)

#--------- Monthly values --------------#

months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
zonal_av_months <- array(data = NaN, dim = c(180,12))
for (i in 1: length(months)){
  temp_global_month <- get_NOAA('temperature', 1, months[i] )
  # filter for a certain depth:
  depth_layer <- 10
  temp_depthlayer <- filter_NOAA(temp_global_month, depth = depth_layer)
  # extract from NOAA framework
  temp_extracted <- temp_depthlayer$t_an
  # calculate zonal average
  zonal_av_months[,i] <- colMeans(temp_extracted, na.rm = TRUE)  # removing Nans
}

# exporting temperature data to csv file
write.csv(zonal_av_months,file= file.path(filelocation, 'modern_SST_data_monthly.csv'))


