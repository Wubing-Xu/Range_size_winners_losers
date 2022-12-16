## submit GBIF occurrence download requests in PC
# debug r-scripts to download occurrences, and then clean occurrence, which will be run in HPC

rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(rgbif)
library(CoordinateCleaner)
library(rnaturalearth)
library(rgeos)

# gbif account 
user <- "wubing" # gbif.org username 
pwd <- "xxxxxxxxx" # gbif.org password
email <- "wbingxu@gmail.com" # email 

# path directory for store gbif occurrences
path_gbif <- "data/gbif/distribution_records"


# download world land and ocean from natural earth
land  <- rnaturalearth::ne_download(scale=10, type="land", category='physical', load=TRUE)
ocean  <- rnaturalearth::ne_download(scale=10, type="ocean", category='physical', load=TRUE)
bound  <- rnaturalearth::ne_download(scale=10, type="wgs84_bounding_box", category='physical', load=TRUE)

# add a 0.1 degree buffer, but keep the extent within -180 - 180 degrees
# buffered land
buffland <- rgeos::gBuffer(land, width=0.1, byid=TRUE)
buffland <- gIntersection(buffland, bound, byid=TRUE)
buffland <- SpatialPolygonsDataFrame(buffland, land@data, match.ID=FALSE)
# buffered ocean
buffocean <- rgeos::gBuffer(ocean, width=0.1, byid=TRUE)
buffocean <- gIntersection(buffocean, bound, byid=TRUE)
buffocean <- SpatialPolygonsDataFrame(buffocean, ocean@data, match.ID=FALSE)

# load checklist
load("data/combined_checklists/spgbif_habitat.RDATA")

# the specieskey from gbif which are used to download occurrences 
gbif_specieskey <- spgbif %>% pull(specieskey) %>% unique()


## submit download request for GBIF occurrence data
# divide all species into 40 groups, which are used for data requests separately
nsp <- length(gbif_specieskey) # 19111 species
nsp.bin <- ceiling(nsp/40) # number of species in each group
bins <- rep(1:40, times =nsp.bin) 
bins <- bins[1:nsp]

# submit data requests
keys <- vector()
for (i in 41:42){
  x <- occ_download(
    pred_in("taxonKey", gbif_specieskey[bins==i]),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user=user, pwd=pwd, email=email
  )
  keys[i] <- x[1]
  print(i)
  Sys.sleep(400)
}

save(spgbif, gbif_specieskey, buffland, buffocean, keys, file="data/gbif/data_to_get_occurrences.RDATA")

## submit all species together. get the doi for citation. DOI: 10.15468/dl.6vdkbn
x2 <- occ_download(
  pred_in("taxonKey", gbif_specieskey),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user, pwd=pwd, email=email
)


####################
## the below lines were just used for debugging in PC
# to formal computation, run the 3b_get_gbif_occurrences_HPC_array.R in HPC

####################
# get GBIF occurrences for those species that have not been included in previous extraction

# the species that have been included in previous extraction
load("data/gbif/data_to_get_occurrences_20220205.RDATA")
load("data/Species_meta_rangesize.RDATA")

# load species list
load("data/combined_checklists/spgbif_habitat.RDATA")

# all specieskey from gbif
gbif_specieskey <- spgbif %>% pull(specieskey) %>% unique()

# the species key with no estimates of range size in the previous extraction
gbif_specieskey <- gbif_specieskey[!gbif_specieskey %in% pull(spsuma, specieskey)]

## submit download request for GBIF occurrence data
# divide all species into several groups, which are used for data requests separately
nsp <- length(gbif_specieskey)
nsp.bin <- ceiling(nsp/10) # number of species in each group
bins <- rep(111:120, times =nsp.bin) 
bins <- bins[1:nsp]

# submit data requests
# keys <- vector()
for (i in 111:120){
  x <- occ_download(
    pred_in("taxonKey", gbif_specieskey[bins==i]),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user=user, pwd=pwd, email=email
  )
  keys[i] <- x[1]
  print(i)
  Sys.sleep(400)
}

save(spgbif, gbif_specieskey, buffland, buffocean, keys, file="data/gbif/data_to_get_occurrences_xxxxxxxx.RDATA")


table(pull(spsuma, specieskey) %in% gbif_specieskey)
table(gbif_specieskey %in% pull(spsuma, specieskey))
length(gbif_specieskey) - nrow(spsuma)



#################################
## download and clean GBIF occurrences

# load species list and land/ocean boundary
load("data/gbif/data_to_get_occurrences_20220205.RDATA")

# GBIF download keys
keys <- occ_download_list(user=user, pwd=pwd, limit=310, start=0)
# key <- keys$results$key[1]
key <- keys$results$key[keys$results$totalRecords <100000 & keys$results$totalRecords >10000][1]

# get occurrences from GBIF
occ <- try(occ_download_get(key, path = path_gbif), silent = TRUE)
occ <- occ_download_import(key=key, path = path_gbif)
occ <- occ %>% 
  setNames(tolower(names(.))) %>%
  filter(specieskey %in% gbif_specieskey)


# summary of species information
spsuma_temp <- occ %>% count(specieskey) %>% rename(n_raw = n)
spsuma <- spgbif %>% inner_join(spsuma_temp)

# clean gbif occurrences
occ <- occ %>% 
  filter(occurrencestatus  == "PRESENT") %>%
  filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")) %>%
  filter(!basisofrecord == "MATERIAL_SAMPLE") %>%  #remove *most* metagenomics records 
  filter(!publishingorgkey == "ab733144-7043-4e88-bd4f-fca7bf858880") %>% 
  filter(!taxonrank == "UNRANKED") %>% 
  filter(coordinateuncertaintyinmeters <= 100000 | is.na(coordinateuncertaintyinmeters)) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>%
  cc_gbif(buffer = 1000) %>% 
  clean_coordinates(tests = c("capitals","centroids","equal","institutions","zeros"), 
                    capitals_rad = 1000, centroids_detail = "country", value = "clean") %>%
  mutate(dplyr::across(c(decimallatitude, decimallongitude), round,2)) %>%
  dplyr::select(specieskey, species, decimallongitude, decimallatitude)

spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_clean = n))  

#remove duplicates
occ <- distinct(occ)
spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_nodupl = n))

# flag records in land using buffered land boundary
occ$inland <- cc_sea(occ, ref = buffland, value = "flagged")
# flag records in ocean using buffered ocean boundary  
occ$inocean <- cc_sea(occ, ref = buffocean, value = "flagged")

# number of records in buffered land and ocean
spsuma_temp <- occ %>% count(specieskey, inland) %>% filter(inland) %>% dplyr::select(specieskey, n_land = n)
spsuma <- spsuma %>% left_join(spsuma_temp)
spsuma_temp <- occ %>% count(specieskey, inocean) %>% filter(inocean) %>% dplyr::select(specieskey, n_ocean = n)
spsuma <- spsuma %>% left_join(spsuma_temp)

# some species have no information of land or sea habitats.
#If 95%> records in one habitat and more than in the other habitat, we assume it as the main habitat
spsuma <- spsuma %>% 
  mutate(keep_new = ifelse(n_land/n_nodupl > 0.95  & n_land > n_ocean & is.na(keep) ,"land", keep)) %>%
  mutate(keep_new = ifelse(n_ocean/n_nodupl > 0.95 & n_ocean > n_land & is.na(keep_new) ,"ocean", keep_new)) %>%
  relocate(keep_new, .after = keep)

# remove records in sea for land species and remove records in land for sea species
occ <- occ %>% 
  filter(!(!inland & specieskey %in% pull(spsuma[spsuma$keep_new == "land",], specieskey))) %>%
  filter(!(!inocean & specieskey %in% pull(spsuma[spsuma$keep_new == "ocean",], specieskey)))

spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_final = n))


save(occ, spsuma, file = paste(path_gbif, paste0("occ_", key, ".RDATA"),sep="/"))
