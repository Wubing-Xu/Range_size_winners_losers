# download occurrences, and then clean occurrence
# I submit array jobs in HPC and save clean occurrences for each individual job


#install and load packages
packages <- c("tidyverse","dplyr","CoordinateCleaner","rnaturalearth","rgeos","rgbif")

for(x in packages){
    if(!require(x,character.only=TRUE, lib.loc = "/gpfs0/home/wubing/R/library")){
		install.packages(x, repos = "http://cran.us.r-project.org", lib="/gpfs0/home/wubing/R/library", dependencies = TRUE)
		require(x,character.only=TRUE, lib.loc = "/gpfs0/home/wubing/R/library")
  }
}

# Set user dependent working directories
path2wd <- "/work/wubing/homogenization_occupancy"
setwd(path2wd)

# gbif account 
user <- "wubing" # gbif.org username 
pwd <- "xxxxxxxxxx" # gbif.org password
email <- "wbingxu@gmail.com" # email

# path directory for store gbif occurrences
path_gbif <- "data/gbif/distribution_records"

# species list and land/ocean boundary
load("data/gbif/data_to_get_occurrences_20220205.RDATA")

# download keys
## get task id from array jobs submitted in HPC cluster 
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
key <- keys[task]


# A function to unzip larg files (>4G): from https://stackoverflow.com/questions/42740206/r-possible-truncation-of-4gb-file
decompress_file <- function(directory, file, .file_cache = FALSE) {
  
  if (.file_cache == TRUE) {
    print("decompression skipped")
  } else {
    
    # Set working directory for decompression
    wd <- getwd()
    setwd(directory)
    
    # Run decompression
    decompression <-
      system2("unzip",
              args = c("-o", # include override flag
                       file),
              stdout = TRUE)
    
    # Reset working directory
    setwd(wd); rm(wd)
    
    # Test for success criteria
    # change the search depending on 
    # your implementation
    if (grepl("Warning message", tail(decompression, 1))) {
      print(decompression)
    }
  }
}    

# A self_defined function of occ_download_import for large files
occ_download_import_self <- function (key = NULL, path = ".", fill = FALSE, 
          encoding = "UTF-8", ...) 
{
  decompress_file(directory = path, file = paste0(key, ".zip"))
  
  targetpath <- sprintf("%s/%s.csv", path, key)
  df <- data.table::fread(targetpath, data.table = FALSE, fill = fill, 
                          encoding = encoding)
  df$countryCode[is.na(df$countryCode)] <- "NA"
  df <- structure(tibble::as_tibble(df), type = "single")
  
  file.remove(targetpath)
  return(df)
}

# get occurrences from GBIF
occ <- try(occ_download_get(key, path = path_gbif), silent = TRUE)
occ   <- occ_download_import_self(key=key, path = path_gbif)
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
  dplyr::select(specieskey, species, decimallongitude, decimallatitude, year)

spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_clean = n))  


# add period for each occurrences (a species in a location) as 1 if data were sampled before 2000, and as 2 if after 2000, 
# and as 3 (1+2) if both before and after 2000
# and remove duplicates
occ <- occ %>% 
  mutate(period = ifelse(year>2000 | is.na(year), 2, 1)) %>%
  dplyr::select(-year) %>% 
  distinct() %>%
  group_by(specieskey, decimallongitude, decimallatitude) %>% 
  mutate(period = sum(period)) %>% 
  ungroup() %>%
  distinct()

spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_nodupl = n))  

# flag records in land sing buffered land boundary
occ$inland <- cc_sea(occ, ref = buffland, value = "flagged")
# flag records in ocean using buffered ocean boundary  
occ$inocean <- cc_sea(occ, ref = buffocean, value = "flagged")

# number of records in buffered land and ocean
spsuma_temp <- occ %>% count(specieskey, inland) %>% filter(inland) %>% select(specieskey, n_land = n)
spsuma <- spsuma %>% left_join(spsuma_temp)
spsuma_temp <- occ %>% count(specieskey, inocean) %>% filter(inocean) %>% select(specieskey, n_ocean = n)
spsuma <- spsuma %>% left_join(spsuma_temp)

# remove records in sea for land species and remove records in land for sea species
occ <- occ %>% filter(!(!inland & specieskey %in% pull(spsuma[spsuma$keep == "land",], specieskey))) %>%
  filter(!(!inocean & specieskey %in% pull(spsuma[spsuma$keep == "ocean",], specieskey)))
spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_final = n)) 

# some species have no information of land or sea habitats.
#If 95%> records in one habitat and more than in the other habitat, we assume it as the main habitat
spsuma <- spsuma %>% 
  mutate(keep_new = ifelse(n_land/n_nodupl > 0.95 & n_land > n_ocean & is.na(keep) ,"land", keep)) %>%
  mutate(keep_new = ifelse(n_ocean/n_nodupl > 0.95 & n_ocean > n_land & is.na(keep_new) ,"ocean", keep_new))

save(occ, spsuma, file = paste(path_gbif, paste0("occ_", key, ".RDATA"),sep="/"))
