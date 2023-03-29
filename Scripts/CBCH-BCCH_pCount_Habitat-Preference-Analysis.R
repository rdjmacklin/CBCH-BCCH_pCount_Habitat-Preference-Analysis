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

packages <- c("tidyverse", "REdaS", "unmarked", "suncalc", "lubridate", 
              "measurements", "sf", "mapview", "ggplot2")

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
    rm(installed)
    
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

## Make an sf object to visualize points.

birds_sf <- st_as_sf(birds, coords = c("Longitude", "Latitude"),
                     crs = 4326)
## Load in tree data

trees <- read_csv("./Data/Raw/Field_Data/TreeData.csv")

trees$AngleToBase <- deg2rad(trees$AngleToBase)
trees$AngleToTop <- deg2rad(trees$AngleToTop)

trees$RawAdjustment[is.na(trees$RawAdjustment) == TRUE] <- 0

trees$CanopyHeight <- case_when(trees$Method == "LevelDistanceLRF" | trees$Method == "LevelDistanceEstimate" ~ round(sqrt(trees$HypotenuseDistance^2-trees$HorizontalDistance^2) + trees$RawAdjustment,2),
                                trees$Method == "DownslopeDistanceLRF" ~ round(sqrt(trees$HypotenuseDistance^2-trees$HorizontalDistance^2),2) + round(sqrt(trees$DownslopeHypotenuseDistance^2-trees$HorizontalDistance^2),2),
                                trees$Method == "DirectLRF" | trees$Method == "DirectEstimate" ~ trees$DirectHeight + trees$RawAdjustment,
                                trees$Method == "LevelAngleSmartphoneDistanceEstimate" | trees$Method == "LevelAngleSmartphoneDistanceLRF" ~ round((trees$HorizontalDistance*tan(trees$AngleToTop)) + trees$RawAdjustment, 2),
                                trees$Method == "UpslopeAngleSmartphoneDistanceEstimate" ~ round((trees$HorizontalDistance*tan(trees$AngleToTop)) - (trees$HorizontalDistance*trees$AngleToBase), 2),
                                trees$Method == "DownslopeAngleSmartphoneDistanceEstimate" ~ round((trees$HorizontalDistance*tan(trees$AngleToTop)) + (trees$HorizontalDistance*tan(trees$AngleToBase)),2)
)

treesBySite <- as.data.frame(matrix(ncol = 6, nrow = length(unique(trees$SheetNo))))
names(treesBySite) <- c("Site", "Region", "Elev", "CanopyHeight", "PropCon", "Density")

treesBySite$Site <- unique(paste0(trees$Site, "-", trees$Point))

for(i in treesBySite$Site) {
  treesBySite$Region[treesBySite$Site == i] <- unique(birds$Region[birds$SitePoint == i])
  treesBySite$Elev[treesBySite$Site == i] <- unique(trees$Elevation[trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])
  treesBySite$PropCon[treesBySite$Site == i] <- round(length(na.omit(trees$TreeType[trees$TreeType == "Con" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))]))/length(na.omit(trees$TreeType[trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])),2)
  treesBySite$CanopyHeight[treesBySite$Site == i] <- round(mean(na.omit(trees$CanopyHeight[trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])),2)
  treesBySite$Density[treesBySite$Site == i] <- round(mean(c(1000*length(trees$Direction[trees$Direction == "N" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])/(na.omit(trees$DistanceToFinalTree[trees$Direction == "N" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])*2),
                                                             1000*length(trees$Direction[trees$Direction == "E" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])/(na.omit(trees$DistanceToFinalTree[trees$Direction == "E" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])*2),
                                                             1000*length(trees$Direction[trees$Direction == "S" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])/(na.omit(trees$DistanceToFinalTree[trees$Direction == "S" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])*2),
                                                             1000*length(trees$Direction[trees$Direction == "W" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])/(na.omit(trees$DistanceToFinalTree[trees$Direction == "W" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])*2),
                                                             1000*length(trees$Direction[trees$Direction == "NW" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])/(na.omit(trees$DistanceToFinalTree[trees$Direction == "NW" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])*2),
                                                             1000*length(trees$Direction[trees$Direction == "NE" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])/(na.omit(trees$DistanceToFinalTree[trees$Direction == "NE" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])*2),
                                                             1000*length(trees$Direction[trees$Direction == "SW" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])/(na.omit(trees$DistanceToFinalTree[trees$Direction == "SW" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])*2),
                                                             1000*length(trees$Direction[trees$Direction == "SE" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])/(na.omit(trees$DistanceToFinalTree[trees$Direction == "SE" & trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])*2))
  ), 2)
  treesBySite$CanopyCover[treesBySite$Site == i] <- unique(trees$CanopyCover[trees$Site == substr(i, 1, 4) & as.numeric(trees$Point == substr(i, 6, 6))])
}

###### Lower Mainland

## Start conforming bird data to unmarkedFrameGDR

yDistanceLM <- matrix(nrow = length(unique(birds$SitePoint[birds$Region == "LowerMainland"])), ncol = 2)
rownames(yDistanceLM) <- unique(birds$SitePoint[birds$Region == "LowerMainland"])
colnames(yDistanceLM) <- c("A", "B")

for(i in rownames(yDistanceLM)) {
  
  for(j in colnames(yDistanceLM)) {
    
    yDistanceLM[i,j] <- length(na.omit(birds$BirdNo[birds$Species == "CBCH" & 
                                              birds$SitePoint == i & 
                                              birds$DistanceBin == j]))
    
  }
}

yRemovalLM <- matrix(nrow = length(unique(birds$SitePoint[birds$Region == "LowerMainland"])), ncol = 4)
rownames(yRemovalLM) <- unique(birds$SitePoint[birds$Region == "LowerMainland"])
colnames(yRemovalLM) <- c("A", "B", "C", "D")

for(i in rownames(yRemovalLM)) {
  
  for(j in colnames(yRemovalLM)) {
    
    yRemovalLM[i,j] <- length(na.omit(birds$BirdNo[birds$Species == "CBCH" & 
                                             birds$SitePoint == i & 
                                             birds$TimeBin == j &
                                             birds$DistanceBin %in% c("A", "B")]))
    
  }
}

## Gather site covariates

siteCovsLM <- as.data.frame(matrix(nrow = length(unique(birds$SitePoint[birds$Region == "LowerMainland"])), ncol = 4))
rownames(siteCovsLM) <- unique(birds$SitePoint[birds$Region == "LowerMainland"])
colnames(siteCovsLM) <- c("PropCon", "CanopyCover", "Density", "TimeSinceSunrise")

for(i in rownames(siteCovsLM)) {
  
  for(j in colnames(siteCovsLM)[1:3]) {
    
    siteCovsLM[i,j] <- treesBySite[treesBySite$Site == i, j]
    
  }
  
  siteCovsLM[i, "TimeSinceSunrise"] <- unique(birds$TimeSinceSunrise[birds$SitePoint == i])

}

## Fit model!

CBCH_umFrameLM <- unmarkedFrameGDR(yDistanceLM[c(1:(which(rownames(yDistanceLM) == "HALU-4")-1), (which(rownames(yDistanceLM) == "HALU-4") + 1):nrow(yDistanceLM)),],
                                 yRemovalLM[c(1:(which(rownames(yRemovalLM) == "HALU-4")-1), (which(rownames(yRemovalLM) == "HALU-4") + 1 ):nrow(yRemovalLM)),],
                                 numPrimary = 1, dist.breaks = c(0, 25, 50), unitsIn = "m", 
                                 siteCovs = siteCovsLM[c(1:(which(rownames(siteCovsLM) == "HALU-4")-1), (which(rownames(siteCovsLM) == "HALU-4") + 1):nrow(siteCovsLM)),],
                                 period.lengths = c(2.5, 2.5, 2.5, 2.5))

CBCH_umLM <- gdistremoval(lambdaformula = ~ scale(PropCon) + scale(CanopyCover) + 1 , phiformula = ~ 1, removalformula = ~scale(TimeSinceSunrise) + 1, distanceformula = ~ 1,
                          CBCH_umFrameLM, output = "density", unitsOut = "ha", mixture = "ZIP")


###### Vancouver Island

## Start conforming bird data to unmarkedFrameGDR

yDistanceVI <- matrix(nrow = length(unique(birds$SitePoint[birds$Region == "VancouverIsland"])), ncol = 2)
rownames(yDistanceVI) <- unique(birds$SitePoint[birds$Region == "VancouverIsland"])
colnames(yDistanceVI) <- c("A", "B")

for(i in rownames(yDistanceVI)) {
  
  for(j in colnames(yDistanceVI)) {
    
    yDistanceVI[i,j] <- length(na.omit(birds$BirdNo[birds$Species == "CBCH" & 
                                              birds$SitePoint == i & 
                                              birds$DistanceBin == j]))
    
  }
}

yRemovalVI <- matrix(nrow = length(unique(birds$SitePoint[birds$Region == "VancouverIsland"])), ncol = 4)
rownames(yRemovalVI) <- unique(birds$SitePoint[birds$Region == "VancouverIsland"])
colnames(yRemovalVI) <- c("A", "B", "C", "D")

for(i in rownames(yRemovalVI)) {
  
  for(j in colnames(yRemovalVI)) {
    
    yRemovalVI[i,j] <- length(na.omit(birds$BirdNo[birds$Species == "CBCH" & 
                                             birds$SitePoint == i & 
                                             birds$TimeBin == j &
                                             birds$DistanceBin %in% c("A", "B")]))
    
  }
}

## Gather site covariates

siteCovsVI <- as.data.frame(matrix(nrow = length(unique(birds$SitePoint[birds$Region == "VancouverIsland"])), ncol = 3))
rownames(siteCovsVI) <- unique(birds$SitePoint[birds$Region == "VancouverIsland"])
colnames(siteCovsVI) <- c("PropCon", "CanopyCover", "TimeSinceSunrise")

for(i in rownames(siteCovsVI)) {
  
  for(j in colnames(siteCovsVI)[1:2]) {
    
    siteCovsVI[i,j] <- treesBySite[treesBySite$Site == i, j]
    
  }
  
  siteCovsVI[i, "TimeSinceSunrise"] <- unique(birds$TimeSinceSunrise[birds$SitePoint == i])
  
}

## Fit model!

CBCH_umFrameVI <- unmarkedFrameGDR(yDistanceVI,
                                   yRemovalVI,
                                   numPrimary = 1, dist.breaks = c(0, 25, 50), unitsIn = "m", 
                                   siteCovs = siteCovsVI,
                                   period.lengths = c(2.5, 2.5, 2.5, 2.5))

CBCH_umVI <- gdistremoval(lambdaformula = ~ scale(PropCon) + scale(CanopyCover) + 1 , phiformula = ~  1, removalformula = ~ scale(TimeSinceSunrise) + 1, distanceformula = ~1,
                          CBCH_umFrameVI, output = "density", unitsOut = "ha", mixture = "ZIP")

## Prediction surface

pred_surface_pCon <- as.data.frame(cbind(seq(from = 0, to = 1, by = 0.02), 100))
colnames(pred_surface_pCon) <- c("PropCon", "CanopyCover")

pred_surface_phi <- as.data.frame(seq(from = 0, to = 400, by = 5))
colnames(pred_surface_phi) <- "TimeSinceSunrise"

LM_Pred <- cbind(predict(CBCH_umLM, newdata = pred_surface_pCon, type = "lambda", appendData = TRUE), "LowerMainland")
names(LM_Pred)[length(names(LM_Pred))] <- "Region"

VI_Pred <- cbind(predict(CBCH_umVI, newdata = pred_surface_pCon, type = "lambda", appendData = TRUE), "VancouverIsland")
names(VI_Pred)[length(names(VI_Pred))] <- "Region"

predictions_PropCon <- rbind(LM_Pred, VI_Pred)

ggplot(data = predictions_PropCon, aes(x= PropCon, y = Predicted, group = Region)) +
  geom_jitter(aes(color = Region)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Region, colour = Region)) +
  theme_classic() +
  ylim(0,100)

pred_surface_canCov <- as.data.frame(cbind(0.5, seq(from = 0, to = 100, by = 5)))
names(pred_surface_canCov) <- names(pred_surface_pCon)

LM_Pred_canCov <- cbind(predict(CBCH_umLM, newdata = pred_surface_canCov, type = "lambda", appendData = TRUE), "LowerMainland")
names(LM_Pred_canCov)[length(names(LM_Pred_canCov))] <- "Region"

VI_Pred_canCov <- cbind(predict(CBCH_umVI, newdata = pred_surface_canCov, type = "lambda", appendData = TRUE), "VancouverIsland")
names(VI_Pred_canCov)[length(names(VI_Pred_canCov))] <- "Region"

predictions_canCov <- rbind(LM_Pred_canCov, VI_Pred_canCov)

ggplot(data = predictions_canCov, aes(x= CanopyCover, y = Predicted, group = Region)) +
  geom_jitter(aes(color = Region)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Region, colour = Region)) +
  theme_classic() +
  ylim(0,100)
