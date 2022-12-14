####################################
## Using GBIF occurrences, determine ranges as area of occupancy (AOO) and extent of occurrences (alpha hulls) 
# This is the r-script to test in PC. The formal run was done in HPC by submitting array jobs


rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(dplyr)
library(rgbif)
library(sp)
library(rgeos)
library(rgdal)
library(cleangeo)
library(raster)
library(maptools)
library(alphahull)


# path directory
path_gbif <- "data/gbif"

# load checklist and land/ocean boundary
load("data/gbif/data_to_get_occurrences_20220205.RDATA")

# load functions to generate alpha hulls or buffed  points
source("rscripts/99_ah2sp.R")
source("rscripts/99_ah_range.R")

# GBIF download keys
key <- keys[1] 

# load GBIF occurrences
load(paste(path_gbif, paste0("distribution_records/occ_", key, ".RDATA"),sep="/"))


# transform land and ocean shapefiles in behrmann equal area projection
behr <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
buffland_behr <- spTransform(buffland, behr)
buffocean_behr <- spTransform(buffocean, behr)


##########
## to calculate area of occupany (AOO), defined as number of grid cells containing records

## extent of raster
ras <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90)
ext <- extent(projectExtent(ras, behr))

## rasters at different resolutions
ras10 <- raster(xmn=ext[1],xmx=ext[2],ymn=ext[3],ymx=ext[4],crs=behr)
res(ras10) <- 10
ras50 <- raster(xmn=ext[1],xmx=ext[2],ymn=ext[3],ymx=ext[4],crs=behr)
res(ras50) <- 50
ras100 <- raster(xmn=ext[1],xmx=ext[2],ymn=ext[3],ymx=ext[4],crs=behr)
res(ras100) <- 100


## Calculate AOO
# use all occurrences
xy_behr <- project(as.matrix(occ[, c("decimallongitude","decimallatitude")]), proj = as.character(behr)) 
occ_rasid <- occ %>% 
  mutate(ras10id = xy_behr %>% cellFromXY(ras10,.)) %>% 
  mutate(ras50id = xy_behr %>% cellFromXY(ras50,.)) %>% 
  mutate(ras100id = xy_behr %>% cellFromXY(ras100,.)) %>% 
  dplyr::select(c(species, specieskey, ras10id:ras100id)) %>% 
  distinct()

aoo <- distinct(occ_rasid, specieskey, ras10id) %>% count(specieskey) %>% rename(aoo10 = n) %>%
  full_join(distinct(occ_rasid, specieskey, ras50id) %>% count(specieskey) %>% rename(aoo50 = n)) %>%
  full_join(distinct(occ_rasid, specieskey, ras100id) %>% count(specieskey) %>% rename(aoo100 = n))
spsuma <- spsuma %>% left_join(aoo) 


#############
## to calculate extent of occurrences (EOO), defined as alpha hulls

## Alpha hulls with alpha =6
ahull6 <- list()
spkeys <-  pull(occ, specieskey) %>% unique()

for(i in 1:length(spkeys)){
  keep <-  pull(spsuma[spsuma$specieskey == spkeys[i],], "keep_new")
  exclude_map <- NULL
  if(keep == "land" & !is.na(keep)) exclude_map <- buffland
  if(keep == "ocean" & !is.na(keep)) exclude_map <- buffocean
  xy <- subset(occ, specieskey == spkeys[i]) %>% 
    dplyr::select("decimallongitude","decimallatitude","specieskey")
  print(i)
  
  # for few species with numerous occurrences, round coordinates to 1 decimal and randomly select 10000 occurrences to construct alpha hulls
  if(nrow(xy) > 10000) {
    xy <- xy %>%
      mutate(dplyr::across(c(decimallatitude, decimallongitude), round, 1)) %>% 
      distinct()
    xy <- xy[sample(1:nrow(xy), min(10000, nrow(xy))), ]
  }
  ahull6[[i]] <- ah_range(xy=xy, alpha=6, buff=10, exclude_map = exclude_map, is.write=FALSE)
}



# summary of alpha hulls
ahull6.suma <- bind_cols(specieskey = spkeys, do.call(rbind, lapply(ahull6, function(x) x[[2]]))) %>% as_tibble()

# ranges as alpha hulls
ahull6 <- lapply(ahull6, function(x) x[[1]])
ahull6 <- do.call("rbind", ahull6)


# EOO: areas of alpha hulls. Here I apply a zero-width buffer to avoid orphaned hole problem
eoo <- data.frame(specieskey = spkeys, ahull6=NA)
eoo[,2] <- spTransform(ahull6, behr) %>% gBuffer(byid=TRUE, width=0) %>% gArea(byid=TRUE)
spsuma <- spsuma %>% left_join(eoo) 


# save results
save(spsuma, file = paste(path_gbif,paste0("ranges/rangesize_", key, ".RDATA"),sep="/"))
save(ras10, ras50, ras100, occ_rasid, ahull6.suma, ahull6, 
     file = paste(path_gbif,paste0("ranges/range_details_", key, ".RDATA"),sep="/"))

