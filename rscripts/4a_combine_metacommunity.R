## Combine the filtered datasets from four database (BioTIME, RivFishTime, InsectChange, Metacommunity Resurvey) and unify metadata  
## and add information such as how many species in each database have range size estimates, the number of occurrence of each species in assemblage data and GBIF


rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","vegan", "reshape2", "sf")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


load("data/Metacommunity_Resurvey/MetacommunityResurvey_filtered.RDATA")
load("data/BioTIME/BioTIME_4locations_10years_filtered_final.RDATA")
load("data/RivFishTIME/RivFishTIME_4sites_10years_filtered.RDATA")
load("data/Role_insects/inset_metacommunities_filtered.RDATA")
load("data/combined_checklists/splist_gbif.RDATA")
load("data/Species_rangesize.RDATA")



##############
## BioTIME data

# update species with GBIF names and indicate whether estimates of range size were calculated
bt <- bt_filtered %>% 
  # filter(!STUDY_ID %in% c(70, 458:465)) %>% # remove problematic studies 
  left_join(splist.gbif %>% dplyr::select(input.species, species, specieskey), by = c("GENUS_SPECIES" = "input.species")) %>% 
  mutate(species = ifelse(is.na(species), GENUS_SPECIES, species)) %>%
  dplyr::select(-GENUS_SPECIES) %>%
  filter((sum.allrawdata.ABUNDANCE > 0 & !is.na(sum.allrawdata.ABUNDANCE)) | (sum.allrawdata.BIOMASS > 0 & !is.na(sum.allrawdata.BIOMASS))) %>% 
  left_join(spsuma %>% dplyr::select(specieskey, aoo10), by ="specieskey") %>%
  mutate(has_range = ifelse(is.na(aoo10), "no", "yes")) %>%
  dplyr::select(-aoo10)

## prepare the dataset used to calculate occupancy
# select the needed columns and calculate number of samples and duration
# choose years with at least 4 samples, studies with >=2 years, and duration > 10 years,
bt_input <- bt %>% dplyr::select(STUDY_ID, YEAR, location, species, specieskey, has_range, LATITUDE, LONGITUDE) %>% 
  rename(studyID = STUDY_ID, year = YEAR, sample = location, latitude = LATITUDE, longitude = LONGITUDE) %>% 
  distinct() %>%
  group_by(studyID, year) %>%
  mutate(n_samp = n_distinct(sample)) %>%
  filter(n_samp > 3) %>%
  group_by(studyID) %>%
  mutate(all_samp = n_distinct(sample),
         min_samp = min(n_samp),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)


# check how many studies and their attributes
bt_studies <- bt_input %>% distinct(studyID, all_samp, min_samp, n_years, duration) #108 studies
table(bt_studies$all_samp)
table(bt_studies$min_samp) # 90 studies >= 10
table(bt_studies$n_years)  # 71 studies >= 10
table(bt_studies$duration)


## meta data for selected studies from BioTIME
bt_input_meta <- bt %>% 
  rename(studyID = STUDY_ID, year = YEAR, sample = location, study_name = TITLE) %>% 
  unite(col = event, studyID, year, sample, remove=FALSE) %>%
  filter(event %in% (bt_input %>%  unite(col = event, studyID, year, sample) %>% pull(event))) %>%
  distinct(studyID, sample, LATITUDE, LONGITUDE, CENT_LONG, CENT_LAT, REALM, TAXA, CLIMATE, GRAIN_SQ_KM, AREA_SQ_KM, study_name) %>%
  distinct(studyID, sample, .keep_all =TRUE) %>%
  group_by(studyID) %>%
  mutate(min_lat = min(LATITUDE, na.rm = TRUE),
         max_lat = max(LATITUDE, na.rm = TRUE),
         min_long = min(LONGITUDE, na.rm = TRUE),
         max_long = max(LONGITUDE, na.rm = TRUE),
         cent_lat = mean(LATITUDE, na.rm = TRUE),
         cent_long = mean(LONGITUDE, na.rm = TRUE)) %>%
  ungroup()

# update central longitude for several studies spanning 180 degree
bt_input_meta <- bt_input_meta %>% 
  mutate(long.new = ifelse((max_long - min_long) > 180, LONGITUDE, NA)) %>% 
  mutate(long.new = ifelse(long.new < 0, long.new + 180, long.new - 180)) %>% 
  group_by(studyID) %>%
  mutate(cent_long_new = mean(long.new, na.rm = TRUE)) %>% 
  ungroup() %>%
  mutate(cent_long_new = ifelse(cent_long_new < 0, cent_long_new + 180, cent_long_new - 180)) %>%
  # distinct(studyID, CENT_LONG, cent_long, cent_long_new)
  mutate(cent_long = ifelse(!is.na(cent_long_new), cent_long_new, cent_long)) %>%
  dplyr::select(-c(long.new, cent_long_new))


## Calculate the area of samples within studies
# the spatial points
bt_input_extent <- bt_input_meta %>% 
  distinct(studyID, LATITUDE , LONGITUDE, min_long, max_long ) %>%
  # update longitude for several studies spanning 180 degree
  mutate(long.new = ifelse((max_long - min_long) > 180, LONGITUDE, NA)) %>% 
  mutate(long.new = ifelse(long.new < 0, long.new + 180, long.new - 180)) %>% 
  mutate(LONGITUDE = ifelse(!is.na(long.new), long.new, LONGITUDE)) %>% 
  # determine spatial objects
  distinct(studyID, LATITUDE , LONGITUDE) %>%
  st_as_sf(coords = c('LONGITUDE', 'LATITUDE'))  %>% 
  st_set_crs(4326) %>%
  group_by(studyID) %>% 
  summarise(geometry = st_combine(geometry))

# areas are defined with convex hull
bt_input_area <- bt_input_extent %>% 
  st_convex_hull() %>% 
  st_area() %>%
  as.numeric()
bt_input_area <-bind_cols(studyID = bt_input_extent$studyID, area = round(bt_input_area/10^6,1)) #the unit of area is km2


# Compare the area and coordinates between the determined based on sites and the raw values
bt_input_meta_check <-  bt_input_meta %>% 
  distinct(studyID, min_lat, max_lat, min_long, max_long, CENT_LONG, CENT_LAT, cent_lat, cent_long, GRAIN_SQ_KM, AREA_SQ_KM) %>% 
  left_join(bt_input_area)

ggplot(bt_input_meta_check) + geom_point(aes(area, AREA_SQ_KM)) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_x_log10() + scale_y_log10()

ggplot(bt_input_meta_check) + geom_point(aes(cent_lat, CENT_LAT)) + 
  geom_abline(intercept = 0, slope = 1) 

ggplot(bt_input_meta_check) + geom_point(aes(cent_long, CENT_LONG)) + 
  geom_abline(intercept = 0, slope = 1) 

ggplot(bt_input_meta_check) + 
  geom_text(aes(cent_long, cent_lat,  label = studyID, colour = "new")) +
  geom_text(aes(CENT_LONG, cent_lat, label = studyID, colour = "raw")) +
  xlim(-180, 180) + ylim(-75, 75)


# choose the determined area based on sites 
bt_input_area <- bt_input_meta_check %>% 
  mutate(mean_grain_m2 = GRAIN_SQ_KM*10^6,
         sd_grain_m2 = NA,
         extent_km2 = area) %>%
  mutate(mean_grain_m2 = ifelse(mean_grain_m2 == 0, NA, mean_grain_m2)) %>%
  dplyr::select(studyID, mean_grain_m2, sd_grain_m2, extent_km2)

# the start and end years of studies
bt_input_year <- bt_input %>% 
  group_by(studyID) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

# combined meta data with area and start and end years
bt_input_meta <- bt_input_meta %>% 
  dplyr::select(studyID, taxon = TAXA, realm = REALM, climate = CLIMATE, study_name, min_lat, max_lat, min_long, max_long, cent_lat, cent_long) %>%
  distinct() %>%
  inner_join(bt_input_area) %>%
  inner_join(bt_input_year)

# 7 studies with taxa recorded as "All"
bt_studies_taxon_all <- bt_input_meta %>% filter(taxon == "All") %>% distinct(studyID) %>% pull()
bt %>% filter(STUDY_ID %in% bt_studies_taxon_all) %>% distinct(STUDY_ID, TAXA, ORGANISMS)

# 11 studies with taxa recorded as "Benthos"
bt_studies_taxon_benthos <- bt_input_meta %>% filter(taxon == "Benthos") %>% distinct(studyID) %>% pull()
bt %>% filter(STUDY_ID %in% bt_studies_taxon_benthos) %>% distinct(STUDY_ID, TAXA, ORGANISMS)



##########
# metacommunity resurvey data

#  update species with GBIF names and indicate whether estimates of range size were calculated
mr <- mr_filtered %>% 
  rename(input_species = species) %>%
  left_join(splist.gbif %>% dplyr::select(input.species, species, specieskey), by = c("input_species" = "input.species")) %>% 
  mutate(species = ifelse(is.na(species), input_species, species)) %>%
  dplyr::select(-input_species) %>%
  left_join(spsuma %>% dplyr::select(specieskey, aoo10), by ="specieskey") %>%
  mutate(has_range = ifelse(is.na(aoo10), "no", "yes")) %>%
  dplyr::select(-aoo10) 


## prepare the dataset used to calculate occupancy
# select the needed columns and calculate number of samples and duration
# choose years with at least 4 samples, studies with >=2 years, and duration > 10 years,
mr_input <- mr %>% 
  dplyr::select(studyID, year, sample, species, specieskey, has_range, latitude, longitude) %>% 
  distinct() %>%
  group_by(studyID, year) %>%
  mutate(n_samp = n_distinct(sample)) %>%
  filter(n_samp > 3) %>%
  group_by(studyID) %>%
  mutate(all_samp = n_distinct(sample),
         min_samp = min(n_samp),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)


# check how many studies and their attributes
mr_studies <- mr_input %>% distinct(studyID, all_samp, min_samp, n_years, duration) #97 studies
table(mr_studies$all_samp)
table(mr_studies$min_samp) # 57 studies >= 10
table(mr_studies$n_years)  # 25 studies with n_years =2, 28 studies >= 10
table(mr_studies$duration)


#### meta data for selected studies from bh
mr_input_meta <- mr %>% 
  unite(col = event, studyID, year, sample, remove=FALSE) %>%
  filter(event %in% (mr_input %>%  unite(col = event, studyID, year, sample) %>% pull(event))) %>%
  distinct(studyID, study_name, dataset_id, sample, latitude, longitude, realm, taxon, 
           alpha_grain_m2, gamma_extent_km2 = gamma_bounding_box_km2, alpha_grain_type, gamma_extent_type = gamma_bounding_box_type) %>% 
  distinct(studyID, sample, .keep_all =TRUE) %>%
  group_by(studyID) %>%
  mutate(min_lat = min(latitude, na.rm = TRUE),
         max_lat = max(latitude, na.rm = TRUE),
         min_long = min(longitude, na.rm = TRUE),
         max_long = max(longitude, na.rm = TRUE),
         cent_lat = mean(latitude, na.rm = TRUE),
         cent_long = mean(longitude, na.rm = TRUE)) %>%
  ungroup()


# no studies with longitudinal spans >180
mr_input_meta %>% filter(max_long -min_long > 180) %>% distinct(studyID)


## Calculate the area of samples within studies
# the spatial points
mr_input_extent <- mr_input_meta %>% 
  distinct(studyID, latitude , longitude, min_long, max_long) %>%
  # update longitude for several studies spanning 180 degree
  mutate(long.new = ifelse((max_long - min_long) > 180 & 
                             (!studyID %in% c("rosenblad_2016a_global" ,"rosenblad_2016b_global")), longitude, NA)) %>% 
  mutate(long.new = ifelse(long.new < 0, long.new + 180, long.new - 180)) %>% 
  mutate(longitude = ifelse(!is.na(long.new), long.new, longitude)) %>% 
  # determine spatial objects
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  distinct(studyID, latitude, longitude) %>%
  st_as_sf(coords = c('longitude', 'latitude'))  %>% 
  st_set_crs(4326) %>%
  group_by(studyID) %>% 
  summarise(geometry = st_combine(geometry))

# areas are defined with convex hull
mr_input_area <- mr_input_extent %>% 
  st_convex_hull() %>% 
  st_area() %>%
  as.numeric()
mr_input_area <-bind_cols(studyID = mr_input_extent$studyID, area = round(mr_input_area/10^6,1)) #the unit of area is km2

# add gamma extent from raw meta data for studies missing values
mr_input_area_inmeta <- mr_input_meta %>% 
  group_by(studyID, study_name) %>% 
  summarise(mean_grain_m2 = mean(alpha_grain_m2, na.rm=TRUE),
            sd_grain_m2 = sd(alpha_grain_m2, na.rm=TRUE),
            extent_km2 = mean(gamma_extent_km2, na.rm=TRUE))
mr_input_area <- full_join(mr_input_area, mr_input_area_inmeta)

# Because many studies didn't provide coordinates and thus can't calculated area based on locations, use the extent provided by database
table(mr_input_area$area == 0) #31 FALSE, 66 TRUE
mr_input_area <- mr_input_area %>% 
  dplyr::select(-area)

# all studies have  extent.
mr_input_area %>% filter(is.na(extent_km2))

# the start and end years of studies
mr_input_year <- mr_input %>% 
  group_by(studyID) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

# combined meta data with area and start and end years
mr_input_meta <- mr_input_meta %>% 
  dplyr::select(studyID, study_name, dataset_id, taxon, realm, min_lat, max_lat, min_long, max_long, cent_lat, cent_long) %>%
  distinct() %>%
  inner_join(mr_input_area) %>%
  inner_join(mr_input_year)

# check whether some studies have no central coordinates
mr_input_meta %>% 
  filter(is.na(cent_lat) | is.na(cent_long)) %>% 
  distinct( studyID, study_name, dataset_id)


#########################
## RivFishTIME data

# update species, and determine whether range size data exists
ft <- ft_filtered %>% 
  left_join(splist.gbif %>% dplyr::select(input.species, species, specieskey), by = c("Species" = "input.species")) %>% 
  mutate(species = ifelse(is.na(species), Species, species)) %>%
  dplyr::select(-Species) %>%
  left_join(spsuma %>% dplyr::select(specieskey, aoo10), by ="specieskey") %>%
  mutate(has_range = ifelse(is.na(aoo10), "no", "yes")) %>%
  dplyr::select(-aoo10) %>% 
  relocate(species, has_range, .after = Year)
 

## prepare the dataset used to calculate occupancy
# select the needed columns and calculate number of samples and duration
# choose years with at least 4 samples, studies with >=2 years, and duration > 10 years,
ft_input <- ft %>% dplyr::select(SourceID, Year, TimeSeriesID, species, specieskey, has_range, Latitude, Longitude) %>%  #, Latitude, Longitude
  rename(studyID = SourceID, year = Year, sample = TimeSeriesID, latitude = Latitude, longitude = Longitude) %>%
  distinct() %>%
  group_by(studyID, year) %>%
  mutate(n_samp = n_distinct(sample)) %>%
  filter(n_samp > 3) %>%
  group_by(studyID) %>%
  mutate(all_samp = n_distinct(sample),
         min_samp = min(n_samp),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)

# check how many studies and their attributes
ft_studies <- ft_input %>% distinct(studyID, all_samp, min_samp, n_years, duration) #36 studies
table(ft_studies$all_samp)
table(ft_studies$min_samp) # # 23 studies >= 10
table(ft_studies$n_years)  # 21 studies >= 10
table(ft_studies$duration)


#### meta data for selected studies from ft
ft_input_meta <- ft %>% 
  rename(studyID = SourceID, year = Year, sample = TimeSeriesID) %>%
  unite(col = event, studyID, year, sample, remove=FALSE) %>%
  filter(event %in% (ft_input %>%  unite(col = event, studyID, year, sample) %>% pull(event))) %>%
  distinct(studyID, sample, Latitude , Longitude) %>%
  filter(!duplicated(.[,c("studyID", "sample")])) %>%
  group_by(studyID) %>%
  mutate(min_lat = min(Latitude, na.rm = TRUE),
         max_lat = max(Latitude, na.rm = TRUE),
         min_long = min(Longitude, na.rm = TRUE),
         max_long = max(Longitude, na.rm = TRUE),
         cent_lat = mean(Latitude, na.rm = TRUE),
         cent_long = mean(Longitude, na.rm =TRUE)) %>%
  ungroup() %>%
  mutate(taxon = "fish",
         realm = "freshwater")

# no study have longitudinal spans >180
ft_input_meta %>% filter(max_long -min_long > 180) %>% distinct(studyID)


## Calculate the area of samples within studies
# the spatial points
ft_input_extent <- ft_input_meta %>% 
  distinct(studyID, Latitude, Longitude) %>%
  st_as_sf(coords = c('Longitude', 'Latitude'))  %>% 
  st_set_crs(4326) %>%
  group_by(studyID) %>% 
  summarise(geometry = st_combine(geometry))

# areas are defined with convex hull
ft_input_area <- ft_input_extent %>% 
  st_convex_hull() %>% 
  st_area() %>%
  as.numeric()
ft_input_area <- bind_cols(studyID = ft_input_extent$studyID, area = round(ft_input_area/10^6,1)) #the unit of area is km2

# add columns in the meta data in other two datasets
ft_input_area <- ft_input_area %>% 
  mutate(mean_grain_m2 =NA,
         sd_grain_m2 =NA,
         extent_km2 = area) %>%
  dplyr::select(-area)
  
# the start and end years of studies
ft_input_year <- ft_input %>% 
  group_by(studyID) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

# combined area with other meta data
ft_input_meta <- ft_input_meta %>% 
  dplyr::select(studyID, taxon, realm, min_lat, max_lat, min_long, max_long, cent_lat, cent_long) %>%
  distinct() %>%
  inner_join(ft_input_area) %>%
  inner_join(ft_input_year)


#########################
## Roles' InsectChange dataset

# update species, and determine whether range size data exists
it <- it_filtered %>% 
  rename(raw_species = species) %>%
  left_join(splist.gbif %>% dplyr::select(input.species, species, specieskey), by = c("raw_species" = "input.species")) %>% 
  mutate(species = ifelse(is.na(species), raw_species, species)) %>%
  dplyr::select(-raw_species) %>%
  left_join(spsuma %>% dplyr::select(specieskey, aoo10), by ="specieskey") %>%
  mutate(has_range = ifelse(is.na(aoo10), "no", "yes")) %>%
  dplyr::select(-aoo10) %>% 
  relocate(species, has_range, .after = Year)

## prepare the dataset used to calculate occupancy
# select the needed columns and calculate number of samples and duration
# choose years with at least 4 samples, studies with >=2 years, and duration > 10 years,
it_input <- it %>% dplyr::select(Datasource_ID, Year, Plot_ID, species, specieskey, has_range, Latitude, Longitude) %>%  #, Latitude, Longitude
  rename(studyID = Datasource_ID, year = Year, sample = Plot_ID, latitude = Latitude, longitude = Longitude) %>%
  distinct() %>%
  group_by(studyID, year) %>%
  mutate(n_samp = n_distinct(sample)) %>%
  filter(n_samp > 3) %>%
  group_by(studyID) %>%
  mutate(all_samp = n_distinct(sample),
         min_samp = min(n_samp),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)

# check how many studies and their attributes
it_studies <- it_input %>% distinct(studyID, all_samp, min_samp, n_years, duration) # 23 studies
table(it_studies$all_samp)
table(it_studies$min_samp)
table(it_studies$n_years)
table(it_studies$duration)


#### meta data for selected studies from it
it_input_meta <- it %>% 
  rename(studyID = Datasource_ID, study_name = Datasource_name, year = Year, sample = Plot_ID) %>%
  unite(col = event, studyID, year, sample, remove=FALSE) %>%
  filter(event %in% (it_input %>%  unite(col = event, studyID, year, sample) %>% pull(event))) %>%
  distinct(studyID, study_name, sample, Latitude , Longitude, Realm, Climate_zone) %>%
  filter(!duplicated(.[,c("studyID", "sample")])) %>%
  group_by(studyID) %>%
  mutate(min_lat = min(Latitude, na.rm = TRUE),
         max_lat = max(Latitude, na.rm = TRUE),
         min_long = min(Longitude, na.rm = TRUE),
         max_long = max(Longitude, na.rm = TRUE),
         cent_lat = mean(Latitude, na.rm = TRUE),
         cent_long = mean(Longitude, na.rm =TRUE)) %>%
  ungroup() %>%
  mutate(taxon = "invertebrates")

# no study have longitudinal spans >180
it_input_meta %>% filter(max_long -min_long > 180) %>% distinct(studyID)


## Calculate the area of samples within studies
# the spatial points
it_input_extent <- it_input_meta %>% 
  distinct(studyID, Latitude, Longitude) %>%
  st_as_sf(coords = c('Longitude', 'Latitude'))  %>% 
  st_set_crs(4326) %>%
  group_by(studyID) %>% 
  summarise(geometry = st_combine(geometry))

# areas are defined with convex hull
it_input_area <- it_input_extent %>% 
  st_convex_hull() %>% 
  st_area() %>%
  as.numeric()
it_input_area <- bind_cols(studyID = it_input_extent$studyID, area = round(it_input_area/10^6,1)) #the unit of area is km2

# add columns in the meta data in other datasets
it_input_area <- it_input_area %>% 
  mutate(mean_grain_m2 =NA,
         sd_grain_m2 =NA,
         extent_km2 = area) %>%
  dplyr::select(-area)


# two studies have no extent. Find the extents manually
it_input_area %>% filter(extent_km2 == 0)

# the area of 1367: a coarse estimate based on map with locations;
# the area of 1408: a coarse estimate based on summed area of catchments
it_input_area_manul <- data.frame(studyID = c("1367", "1408"),
                                  extent_km2 = c(10, 350))
id <- match(it_input_area_manul[, "studyID"], pull(it_input_area, studyID))
it_input_area[id[!is.na(id)],c("extent_km2")] <- it_input_area_manul[!is.na(id), c("extent_km2")]


# the start and end years of studies
it_input_year <- it_input %>% 
  group_by(studyID) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

# combined area with other meta data
it_input_meta <- it_input_meta %>% 
  dplyr::select(studyID, study_name, taxon, realm = Realm, climate = Climate_zone, min_lat, max_lat, min_long, max_long, cent_lat, cent_long) %>%
  distinct() %>%
  inner_join(it_input_area) %>%
  inner_join(it_input_year)



#######
## combine datasets from four database: BioTIME, RivFishTIME, InsectChange, Metacommunity Resurvey

# community data
dat <- bind_rows(
  bt_input %>% mutate(database = "bt", studyID = as.character(studyID)),
  ft_input %>% mutate(database = "ft", studyID = as.character(studyID)),
  it_input %>% mutate(database = "ic", studyID = as.character(studyID), sample = as.character(sample)),
  mr_input %>% mutate(database = "mr", studyID = as.character(studyID)),
  ) %>%
  unite(col = "study", database, studyID, remove = FALSE) %>% 
  relocate(database, .before = studyID)


# meta data
dat_meta <- bind_rows(
  bt_input_meta %>% 
    mutate(database = "bt",
           studyID = as.character(studyID)),
  ft_input_meta %>% 
    mutate(database = "ft",
           studyID = as.character(studyID)),
  it_input_meta %>% 
    mutate(database = "ic",
           studyID = as.character(studyID)),
   mr_input_meta %>% 
    mutate(database = "mr",
           studyID = as.character(studyID)), 
  ) %>%
  relocate(database) %>%
  relocate(climate, .after = realm)  %>%
  mutate(realm = tolower(realm)) %>%
  unite(col = "study", database, studyID, remove = FALSE)


# determine the climate based on central coordinates
dat_meta <- dat_meta %>% 
  mutate(climate = ifelse(abs(cent_lat) <= 23.5, "Tropical", NA),
         climate = ifelse(abs(cent_lat) > 23.5 & abs(cent_lat)<= 60, "Temperate", climate),
         climate = ifelse(abs(cent_lat) > 60, "Polar", climate)) %>%
  relocate(climate, .after = realm)


# determine the proportion of species that have estimates of ranges for each study
data_meta_sprich <- dat %>% 
  group_by(study) %>%
  summarise(sprich = n_distinct(species),
            sprich_range = n_distinct(species[has_range == "yes"]),
            Psprich_range = round(sprich_range/sprich, 4))

# 25studies have < 10 species with estimates of range size
data_meta_sprich %>% filter(sprich_range < 10)


#  total number of species: 25,607 species 
dat %>% distinct(species) %>% nrow() 
# number of species that are matched to GBIF backbone: 19,110
dat %>% filter(species %in% splist.gbif$species) %>% distinct(species) %>% nrow()
# number of species that have estimates of range size: 18,914
dat %>% filter(has_range == "yes") %>% distinct(species) %>% nrow() 


# Update community: select only species that have estimates of range size 
dat <- dat %>% 
  filter(has_range == "yes") %>%
  group_by(study, year) %>%
  mutate(n_samp = n_distinct(sample)) %>% 
  filter(n_samp > 3) %>%
  group_by(study) %>%
  mutate(all_samp = n_distinct(sample),
         min_samp = min(n_samp),
         n_years = n_distinct(year),
         duration = max(year) - min(year) +1,
         sprich = n_distinct(species)) %>% 
  # choose studies with number of species >= 10
  filter(n_years >= 2 & duration >= 10 & sprich >= 10) %>%
  group_by(study) %>%
  mutate(period = ifelse(year == min(year), "first", 
                         ifelse(year == max(year), "last", "intermediate"))) %>%
  relocate(period, .after = year) %>%
  ungroup() %>% 
  # add study type
  inner_join(dat_meta %>% dplyr::select(study, study_name))

# add number of occurrences (unique in the resolution of 0.01 degree) from community data and GBIF, which will be used in sensitive analyses
# number of occurrences from community data
nocc_species <- dat %>% ungroup() %>%
  dplyr::select(species, specieskey, latitude, longitude) %>%
  mutate(dplyr::across(c(latitude, longitude), round, 2)) %>%
  distinct() %>% 
  count(specieskey)
# add number of occurrences from GBIF
nocc_species <- nocc_species %>% 
  rename(nocc_community = n) %>%
  inner_join(spsuma %>% dplyr::select(specieskey, nocc_gbif = n_final)) %>%
  mutate(ratio_nocc_gbif = round(nocc_gbif/nocc_community, 2))

table(nocc_species$ratio_nocc_gbif < 1) # FALSE: 17880, TRUE: 835

# combine community data with number of occurrences
dat <- dat %>% 
  left_join(nocc_species)


# update meta data using updated community data
data_meta_update <- dat %>% 
  group_by(study) %>%
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1,
            min_samp = min(n_samp),
            max_samp = max(n_samp),
            nsamp_same = (min_samp == max_samp))

dat_meta <- dat_meta %>% 
  dplyr::select(study:extent_km2) %>% 
  inner_join(data_meta_update) %>%
  inner_join(data_meta_sprich)


## generalize taxon groups
dat_meta <- dat_meta %>%
  mutate(taxon = tolower(taxon))
table(dat_meta$taxon)

# label taxonomic groups manually
taxon_manual <- data.frame("raw" = c("amphibians", "all", 'benthos', 'birds', 'fish', 
                                     'freshwater invertebrates', 'herpetofauna', 'invertebrates', 'mammals', 'marine invertebrates', 
                                     'marine plants', 'multiple taxa', 'plants', 'terrestrial invertebrates', 'terrestrial plants'),
                           "new" = c("Amphibians and reptiles", "Multiple taxa", "Benthos", "Birds", "Fish", 
                                     "Invertebrates", "Amphibians and reptiles", "Invertebrates", "Mammals", "Invertebrates", 
                                     "Plants", "Multiple taxa", "Plants", "Invertebrates", "Plants")) 
unique(dat_meta$taxon)[!unique(dat_meta$taxon) %in% taxon_manual[,1]]

# add updated taxonomic groups
id <- match(pull(dat_meta, taxon), taxon_manual[,1], nomatch = NA)
table(is.na(id)) #all matched

dat_meta <- dat_meta %>%
  mutate(taxon_new = taxon_manual[id, 2]) %>%
  relocate(taxon_new, .after = taxon)

# transfer the first letter of realms to be capital 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
dat_meta <- dat_meta %>% mutate(realm = firstup(realm))


# check number of studies across database, study_types,realms, taxonomic groups, climate
table(dat_meta[,c( "database")])
table(dat_meta[,c("realm", "database")])
table(dat_meta[,c("climate", "database")])
table(dat_meta[,c("taxon_new", "realm")])
table(dat_meta[,c("taxon_new", "realm", "database")])
table(dat_meta$n_years == 2) # 34 studies vs. 204
table(dat_meta$n_years < 10) # 112 studies vs. 126


# found some errors on the central latitudes; update them manually  
cent_coords_manul <- data.frame(study = c("bt_187", "bt_97"),
                                cent_long = c(-75, 176.667),
                                cent_lat = c(37, 72.81345))
id <- match(cent_coords_manul[, "study"], pull(dat_meta, study))
dat_meta[id, c("cent_lat")] <- cent_coords_manul[, c("cent_lat")]
dat_meta[id, c("cent_long")] <- cent_coords_manul[, c("cent_long")]


save(dat, dat_meta,
     file = "data/Combined_assemblages.RDATA")
