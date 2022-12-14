## filter metacommunity-resurvey dataset
## to choose subsets of datasets with number of locations in each year >=4, duration >=10, 
## and matching locations across years to keep locations in similar configuration across years

rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(raster)
library(vegan)
library(reshape2)

mr <- read_csv("data/Metacommunity_Resurvey/metacommunity-survey-communities.csv")
mr_meta <- read_csv("data/Metacommunity_Resurvey/metacommunity-survey-metadata.csv")
mr_species <- read_csv("data/Metacommunity_Resurvey/manual_community_species_filled_20221003.csv")

mr <- mr %>% dplyr::select(dataset_id, regional, local,  year, species, species_original, value, metric, unit)


#########
# correct some errors

## use the manually checked species names to replace the raw one
mr <- mr %>% left_join(mr_species %>% distinct(dataset_id, species_original, species.new))
# all species are matched
 mr %>% filter(is.na(species.new)) %>% distinct(dataset_id, species, species_original, species.new)

mr <- mr %>%
  mutate(species = ifelse(is.na(species.new), species, species.new)) %>%
  dplyr::select(- species.new)


# magalhaes_2020: the baseline survey of 5 sites were performed at 4 years. Use the average year to replace it.
mr %>% filter(dataset_id == "magalhaes_2020") %>% distinct(regional, local, year)
id <- mr$dataset_id == "magalhaes_2020" & mr$year %in% c(2003:2006)
mr$year[id] <- 2005

id <- mr_meta$dataset_id == "magalhaes_2020" & mr_meta$year %in% c(2003:2006)
mr_meta$year[id] <- 2005

# anderson_2019b: the extent is missing; add it manually
id <- mr_meta$dataset_id == "anderson_2019b"
mr_meta$gamma_bounding_box_km2[id] <- 1.73

# christensen_2021: the extent seems incorrect; use the estimates from data description
id <- mr_meta$dataset_id == "christensen_2021"
mr_meta$gamma_bounding_box_km2[id] <- 782


# szydlowski_2022_snails: the extent is missing
id <- mr_meta$dataset_id == "szydlowski_2022_snails"
mr_meta$gamma_bounding_box_km2[id] <- 2640


# self_define studyID for combining with other database 
mr <- mr %>% unite(col = study_name, dataset_id, regional, sep="_",remove=FALSE)

mr_meta <- mr_meta %>% 
  unite(col = study_name, dataset_id, regional, sep="_", remove=FALSE) %>%
  mutate(studyID = factor(study_name, labels = paste0("sfd_", 1:n_distinct(study_name))),
         studyID =as.character(studyID)) %>%
  relocate(studyID)

# check whether study_name can not be matched between main data and meta data
mr %>%
  filter(!study_name %in% mr_meta$study_name) %>%
  distinct(study_name)

mr_meta %>% filter(!study_name %in% mr$study_name) %>%
  distinct(study_name)

# check whether some datasets have multiple taxon types and reams in the same region 
mr_meta %>% distinct(dataset_id, regional, realm, taxon) %>% filter(duplicated(.[,c("dataset_id", "regional")]))
mr_meta %>% distinct(dataset_id, realm, taxon) %>% filter(duplicated(.[,c("dataset_id")])) 

# check where some samples have duplicated information in the metadata
mr_meta %>% distinct() %>% filter(duplicated(.[,c("dataset_id","year", "regional" ,"local")])) 
mr_meta %>% distinct() %>% filter(duplicated(.[,c("dataset_id","year", "regional" ,"local")]))  %>% distinct(dataset_id)


# combine community and meta data
mr <- mr %>% left_join(mr_meta %>% 
                         dplyr::select(-c(data_pooled_by_authors, data_pooled_by_authors_comment, sampling_years, alpha_grain_comment,
                                          gamma_bounding_box_comment,  gamma_sum_grains_comment, comment, comment_standardisation)) %>% 
                         distinct(), 
                       by =c("study_name", "dataset_id", "year", "regional", "local")) %>%
  relocate(studyID)


# remove studies that have been included in other database
mr_meta <- mr_meta %>%
  # remove three studies that have included in InsectChange
  # valtonen_2018 = Hungary moths;  schuch_2011 = Germany Marchand Schuch, magnuson_2020 = LTER NTL Macroinvertebrates
  # the study "willig_2010" has been included in BioTIME (StudyID = 54), but BioTIME doesn't provide coordinates or plotID for this study. keep it here 
  filter(! dataset_id %in% c("valtonen_2018", "schuch_2011", "magnuson_2020")) 

mr <- mr %>% 
  filter(dataset_id %in% mr_meta$dataset_id) %>%
  rename(sample = local)


# remove studies with few sites and short duration
mr_4loc <- mr %>%
  group_by(studyID, year) %>%
  mutate(n_samp = n_distinct(sample)) %>%
  filter(n_samp > 3) %>%
  group_by(studyID) %>%
  mutate(duration = max(year) - min(year) +1) %>%
  ungroup() %>%
  filter(duration >= 10) %>%
  dplyr::select(!c(n_samp, duration))


############
## the number and locations of sites are somewhat different across years
# match sites across years based on grid-cells to keep similar configurations of sites

# Calculate minimum, maximum latitude and longitude, and spans of longitude and longitude
meta_spat <- mr_4loc %>% 
  distinct(studyID, sample, year, latitude, longitude) %>%
  group_by(studyID) %>%
  mutate(min_latitude =  min(latitude),
         max_latitude =  max(latitude),
         min_longitude =  min(longitude),
         max_longitude =  max(longitude)) %>%
  # If locations cross 180 degree, update the minimum longitude to calculate spans of longitude
  mutate(longitude = ifelse((max_longitude - min_longitude) >180 & longitude < 0, longitude + 360, longitude)) %>%
  mutate(span_latitude = max(latitude) - min(latitude),
         span_longitude = max(longitude) - min(longitude)) %>%
  ungroup() %>%
  distinct(studyID, min_latitude, max_latitude, min_longitude, max_longitude, span_latitude, span_longitude)

# Check whether some  studies have no coordinates. 
# If no coordinates, set latitudinal and longitudinal span as zero
meta_spat$span_latitude[is.na(meta_spat$span_latitude)] <- 0
meta_spat$span_longitude[is.na(meta_spat$span_longitude)] <- 0


# Check spatial distributions of samples in grid-cells.
i = 5
study_spat <- meta_spat[i, c("studyID", "min_longitude", "max_longitude", "min_latitude", "max_latitude", "span_longitude", "span_latitude")] %>% as.data.frame()
# set the resolution as the mean 1/5 spans of longitude and longitude
res <- mean(c(study_spat[, 6]/5, study_spat[, 7]/5)) 
ras <- raster(xmn = study_spat[1, 2] - res, xmx = study_spat[1, 3] + res, 
              ymn = study_spat[1, 4] - res, ymx = study_spat[1, 5] + res, 
              resolution = res)
values(ras) <- 1:length(ras)
plot(ras)
points(mr_4loc %>% 
         distinct(studyID, sample, latitude, longitude) %>% 
         filter(studyID == study_spat[1, 1]) %>% dplyr::select(longitude, latitude))


# add cells for each study. cells are only comparable within studies
mr_4loc <- mr_4loc %>% mutate(cell = NA)
for(i in 1:nrow(meta_spat)){
  study_spat <- meta_spat[i, c("studyID", "min_longitude", "max_longitude", "min_latitude", "max_latitude", "span_longitude", "span_latitude")] %>% as.data.frame()
  
  if(study_spat$span_longitude > 0){
    # set the resolution as the mean 1/5 spans of longitude and longitude
    res <- mean(c(study_spat[, 6]/5, study_spat[, 7]/5))
    ras <- raster(xmn = study_spat[1, 2] - res, xmx = study_spat[1, 3] + res, 
                  ymn = study_spat[1, 4] - res, ymx = study_spat[1, 5] + res, 
                  resolution = res)
    values(ras) <- 1:length(ras)
    id <- mr_4loc$studyID == study_spat[1, 1]
    mr_4loc$cell[id] <- cellFromXY(ras, mr_4loc[id, c("longitude", "latitude")] %>% as.data.frame())
  }
  
  #  studies with no coordinates or only with regional central coordinates. Use the identities of plots as the cells 
  if(study_spat$span_longitude == 0){
    id <- mr_4loc$studyID == study_spat[1, 1]
    mr_4loc$cell[id] <- mr_4loc$sample[id]
  }
}


# filter the dataset to remove rare cells and years with small extent
mr_4loc_filtered <- NULL
mr_years_max_loc <- NULL
for(i in 1:nrow(meta_spat)){
  # perform loop for each study
  study <- mr_4loc %>% 
    filter(studyID == meta_spat$studyID[i])
  
  # calculate number of locations in each cell of each year: rows are years, columns are cells, elements are number of locations
  year_cell <- as.matrix(xtabs( ~ year + cell, data = unique(study[, c("sample","year","cell")]), sparse=TRUE)) 
  
  
  # keep cells with locations in all years or cells with density of locations grater than 50% of the mean value
  cell_loc <- data.frame("cell" = colnames(year_cell),
                         "p_years" = colMeans(year_cell>0),
                         "mean_loc" = colMeans(year_cell)) %>% 
    mutate(p_loc = round(mean_loc/sum(mean_loc),3)) #relative density of locations
  
  # keep cells with locations in all years or with density of locations > the half of the mean
  # if the remaining cells have >= 4 locations
  id_cell <- with(cell_loc, p_years==1 | p_loc > 0.5*1/nrow(cell_loc))
  if(max(rowSums(year_cell[, id_cell, drop=FALSE]))>=4){
    year_cell <-  year_cell[, id_cell, drop=FALSE]
  }
  
  # number of co-occurred (in the same cells) locations between years
  co_loc <- as.matrix(designdist(year_cell, method = "J", terms= "minimum"))
  
  # find which two years (year-pair) have the maximum number of co-occurred locations and duration >=10 years
  # these two years have priority to be kept, and other years will be compared to the two years
  max_co_loc <- melt(co_loc) %>% 
    as_tibble() %>%
    set_names("year1","year2","n_loc") %>%
    mutate(year1 = as.numeric(year1),
           year2 = as.numeric(year2),
           duration = year2 - year1 + 1) %>%
    filter(duration > 9 & n_loc > 3)
  if(nrow(max_co_loc) == 0) {next}
  max_co_loc <-  filter(max_co_loc, n_loc >= 0.9*max(n_loc))
  max_co_loc <- filter(max_co_loc, duration == max(duration))
  max_co_loc <-  filter(max_co_loc, n_loc == max(n_loc))
  
  # keep the cells that have locations in the both determined years
  cell_year1 <- year_cell[rownames(year_cell) %in% unlist(max_co_loc[1,1]),]
  cell_year2 <- year_cell[rownames(year_cell) %in% unlist(max_co_loc[1,2]),]
  cell_shared <- cell_year1 > 0 & cell_year2 > 0
  year_cell <- year_cell[,cell_shared, drop=FALSE]
  
  # other years except the two priority years will be compared to the two years
  # weight cells based on mean number of locations and calculate the sums of weights of cells with locations for each year
  weight <- colMeans(year_cell)/sum(colMeans(year_cell))
  year_cell_binomial <- (year_cell > 0)*1
  cum_weight_year <- year_cell_binomial %*% weight
  
  # remove years with the sum weight less than 80%: the years with small extent
  id_year <- cum_weight_year > 0.9
  year_cell <-  year_cell[id_year, , drop=FALSE]
  
  
  # remove years with less than 50% of mean number of locations, but keep the duration >= 10 years
  id_year1 <- rowSums(year_cell) >= 0.5*mean(rowSums(year_cell))
  year_keep <- as.numeric(rownames(year_cell))[id_year1]
  if((max(year_keep) - min(year_keep) + 1) > 9){
    year_cell <-  year_cell[id_year1, , drop=FALSE]
  }
  
  # keep the locations is the chosen cells and years
  rare_study <- study %>% 
    filter(cell %in% colnames(year_cell) & year %in% rownames(year_cell))
  
  mr_years_max_loc <- bind_rows(mr_years_max_loc , bind_cols(studyID = meta_spat$studyID[i], max_co_loc[1, ]))
  mr_4loc_filtered  <- bind_rows(mr_4loc_filtered, rare_study)
}

ggplot(data = rare_study, aes(longitude, latitude)) + 
  facet_wrap(~year) +
  geom_point(size = 1, alpha = 0.7) + 
  coord_fixed()


# keep only years with at least 4 samples,
# and keep studies with at least 2 time points and duration >10 years
mr_4loc_filtered <- mr_4loc_filtered %>%
  group_by(studyID, year) %>%
  mutate(n_samp = n_distinct(sample)) %>%
  filter(n_samp > 3) %>%
  group_by(studyID) %>%
  mutate(all_samp = n_distinct(sample),
         min_samp = min(n_samp),
         n_years = n_distinct(year),
         duration = max(year) - min(year) +1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)


# check how many studies and their attributes
mr_studies <- mr_4loc_filtered %>% distinct(studyID, all_samp, min_samp, n_years, duration) #97 studies
table(mr_studies$all_samp)
table(mr_studies$min_samp) # 57 studies >= 10
table(mr_studies$n_years)  # 25 studies with n_years =2, 28 studies >= 10
table(mr_studies$duration)

# save the filtered data
mr_filtered <- mr_4loc_filtered %>% dplyr::select(-(n_samp:duration)) 
save(mr_filtered, file = "data/Metacommunity_Resurvey/metacommunityResurvey_filtered.RDATA")



# locations of the filtered dataset
mr_loc_filtered <- mr_4loc_filtered %>% 
  distinct(studyID, sample, year, latitude, longitude)

# locations of the full dataset
mr_loc <- mr_4loc %>% 
  distinct(studyID, sample, year, latitude, longitude) %>% 
  left_join(mr_loc_filtered %>%
              mutate(keep = "yes")) %>%
  mutate(keep = ifelse(is.na(keep), "no", keep)) %>% 
  group_by(studyID) %>%
  mutate(n_coord = n_distinct(latitude, longitude)) %>% 
  ungroup() %>%
  filter(!is.na(latitude) & n_coord >2)


# plot distributions of samples and indicate which records will be removed
pdf('data/Metacommunity_Resurvey/metacommunityResurvey_filtered.pdf', width = 12, height = 10)
id_study <- unique(mr_loc$studyID)
for(i in 1:length(id_study)){
  study <- mr_loc %>% 
    filter(studyID %in% id_study[i])
  
  p <- ggplot(data = study, aes(longitude, latitude)) + 
    facet_wrap(~year) +
    geom_point( aes(colour = keep), size = 1, alpha = 0.7) + 
    labs(title = id_study[i]) +
    theme(legend.position = "top", legend.text = element_text(size=12)) + 
    coord_fixed() +
    scale_color_manual(values=c("yes" = "deepskyblue", "no" = "coral"))
  
  print(p)
}
dev.off()
