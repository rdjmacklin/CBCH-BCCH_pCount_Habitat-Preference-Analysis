################################################################################

# Script Title: Applying Amundson et al.'s (2014) model to point count data
# collected on the Lower Mainland and Vancouver Island, British Columbia, Canada.

# Script Author: Rory Macklin (macklin@zoology.ubc.ca)

# Date: March 23, 2023

################################################################################

## Check working directory. This should be the head directory of the repository.

getwd()

# setwd("../")

## Install required packages (if not already installed).

packages <- c("tidyverse", "REdaS", "unmarked", "suncalc", "lubridate", "measurements")

lib.path <- if(Sys.info()["login"] == "root") {
  .libPaths()
} else {
  "/home/macklin/R_pkg_lib"
}

installed <- as.data.frame(installed.packages(lib = lib.path))

for(i in packages){
  if(!(i %in% installed$Package)) {
    
    print(paste0("Package ", i, " isn't installed. Attempting installation."))
    
    install.packages(i, lib = lib.path, dependencies = TRUE)
    
    installed <- as.data.frame(installed.packages(lib = lib.path))
    
    if(i %in% installed$Package) {
      
      suppressPackageStartupMessages(library(i, lib.loc = lib.path, character.only = TRUE))
      
      print(paste0("Package ", i, " installed successfully and loaded!"))
      
    } else {
      
      print(paste0("Package ", i, " installation unsuccessful."))
      
    }
    
  } else {
    
    suppressPackageStartupMessages(library(i, lib.loc = lib.path, character.only = TRUE))
    print(paste0("Package ", i, " was already installed and has been loaded!"))
    
  }
  
  if(i == packages[length(packages)]) {
    
    rm(lib.path)

  }
}

## Read in our Bird Data

birds <- read_csv("./Data/Raw/Field_Data/BirdData.csv")

## Create some useful columns in the birds dataframe.

birds$SitePoint <- paste0(birds$Site, "-", birds$Point)

## Convert lat and long as read from the GPS (degree decimal minutes) to decimal degrees.

birds$Latitude <- as.numeric(measurements::conv_unit(paste0(birds$LatitudeHr,
                                                 " ",
                                                 ifelse(nchar(birds$LatitudeMin) == 2,
                                                       birds$LatitudeMin,
                                                       paste0(0, birds$LatitudeMin)),
                                                 ".",
                                                 ifelse(nchar(birds$LatitudeSec) == 3,
                                                        birds$LatitudeSec,
                                                        paste0(0, birds$LatitudeSec))
                                                 ), from = "deg_dec_min", to = "dec_deg"))


birds$Longitude <- as.numeric(measurements::conv_unit(paste0(birds$LongitudeHr,
                                                 " ",
                                                 ifelse(nchar(birds$LongitudeMin) == 2,
                                                        birds$LongitudeMin,
                                                        paste0(0, birds$LongitudeMin)),
                                                 ".",
                                                 ifelse(nchar(birds$LongitudeSec) == 3,
                                                        birds$LongitudeSec,
                                                        paste0(0, birds$LongitudeSec))
                                                 ), from = "deg_dec_min", to = "dec_deg")) * -1                

## Use the package "suncalc" to approaximate sunrise times at our points.

birds$SunriseHr <- as.numeric(substr(getSunlightTimes(data = birds %>%
                                             mutate(date = as.Date(paste0(Year, "-", 
                                                                          ifelse(nchar(Month) == 1,
                                                                                 paste0("0", Month),
                                                                                 Month),
                                                                          "-",
                                                                          ifelse(nchar(Day) == 1,
                                                                                 paste0("0", Day),
                                                                                 Day))),
                                                    lat = Latitude,
                                                    lon = Longitude) %>%
                                             select(date, lat, lon), keep = "sunrise", tz = "MST")$sunrise,
                          start = 12, stop= 13)) -1

birds$SunriseHr <- case_when(birds$Month == 2 ~ birds$SunriseHr,
                             birds$Month == 3 & birds$Day < 12 ~ birds$SunriseHr,
                             birds$Month == 3 & birds$Day > 11 ~ birds$SunriseHr + 1)

birds$SunriseMin <- as.numeric(substr(getSunlightTimes(data = birds %>%
                                                        mutate(date = as.Date(paste0(Year, "-", 
                                                                                     ifelse(nchar(Month) == 1,
                                                                                            paste0("0", Month),
                                                                                            Month),
                                                                                     "-",
                                                                                     ifelse(nchar(Day) == 1,
                                                                                            paste0("0", Day),
                                                                                            Day))),
                                                               lat = Latitude,
                                                               lon = Longitude) %>%
                                                        select(date, lat, lon), keep = "sunrise", tz = "MST")$sunrise,
                                     start = 15, stop= 16))
                          
## Subtract sunrise time (in minutes since midnight) from survey start time 
## (in minutes since midnight) to get time since sunrise (in minutes).
birds$TimeSinceSunrise <- ((birds$StartTimeHr * 60) + birds$StartTimeMin) - ((birds$SunriseHr*60) + birds$SunriseMin)

## Start conforming bird data to unmarkedFrameGDR

yDistance <- matrix()
