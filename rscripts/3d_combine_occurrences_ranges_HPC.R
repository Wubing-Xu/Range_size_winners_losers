####################################
## Combine GBIF occurrences (raw and rasterized), alpha hulls and estimates of range size from individual jobs run HPC

#install and load other packages
packages <- c("tidyverse","dplyr","rgbif","sp","rgeos","rgdal","cleangeo","raster","maptools","alphahull")

for(x in packages){
    if(!require(x,character.only=TRUE,lib.loc = "/gpfs0/home/wubing/R/library")){
		install.packages(x,repos = "http://cran.us.r-project.org", lib="/gpfs0/home/wubing/R/library", dependencies = TRUE)
		require(x,character.only=TRUE,lib.loc = "/gpfs0/home/wubing/R/library")
  }
}

# Set user dependent working directories
path2wd <- "/work/wubing/homogenization_occupancy"
setwd(path2wd)

# path directory
path_gbif <- "data/gbif"

# load checklist and land/ocean boundary
load("data/gbif/data_to_get_occurrences.RDATA")

## combine species metadata and range size 
spsuma <- lapply(keys, function(x) get(load(paste(path_gbif,paste0("ranges/rangesize_", x, ".RDATA"), sep="/"))) )
spsuma <- do.call(bind_rows, spsuma)
spsuma <- distinct(spsuma, species, .keep_all = TRUE)

# combine species distributions at four resolutions and species ranges
ranges <- lapply(keys, function(x) mget(load(paste(path_gbif,paste0("ranges/range_details_", x, ".RDATA"), sep="/"))))
ras10 <- ranges[[1]][[1]]
ras50 <- ranges[[1]][[2]]
ras100 <- ranges[[1]][[3]]
occ_rasid <- do.call(bind_rows, lapply(ranges, "[[", 4))
ahull6_suma <- do.call(bind_rows, lapply(ranges, "[[", 5))
ahull6 <- do.call(rbind, lapply(ranges, "[[", 6))


## combine GBIF occurrences
occ <- lapply(keys, function(x) get(load(paste(path_gbif, paste0("distribution_records/occ_", x, ".RDATA"),sep="/"))) )
occ <- do.call(bind_rows, occ)

# save the combined results
save(spsuma, file = paste(path_gbif, paste0("Species_rangesize.RDATA"), sep="/"))
save(ras10, ras50, ras100, occ_rasid, file = paste(path_gbif, paste0("Species_distributions_gridcells.RDATA"), sep="/"))
save(ahull6, ahull6_suma, file = paste(path_gbif, paste0("Species_distributions_rangemaps.RDATA"), sep="/"))
save(occ, file = paste(path_gbif, paste0("GBIF_occurrences.RDATA"), sep="/"))

