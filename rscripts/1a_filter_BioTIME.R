## filter BioTIME database to choose subsets of datasets with number of locations in each year >=4, duration >=10, 
## and matching locations across years to keep locations in similar configuration across years

rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(raster)
library(vegan)
library(reshape2)


###############
# load public and private BioTIME datasets and combine them
bt_public <- read_csv("/iDiv/Data/BioTIME/BioTIMEQuery02_04_2018.csv")
bt_meta_public <- read_csv('/iDiv/Data/BioTIME/BioTIMEMetadata_02_04_2018.csv')
load("/iDiv/Data/BioTIME/biotime_extra_private_20210712.Rdata")

bt <- bind_rows(bt_public, bt) 
bt_meta <- bt_meta_public %>% 
  dplyr::select(- c(GENERAL_TREAT:TREAT_DATE, CONTACT_1:COMMENTS)) %>% 
  bind_rows(bt_meta %>% dplyr::select(- c(MIN_SAMPLE, AREA)))

save(bt, bt_meta, file = "data/BioTIME/biotime_2021.RDATA")


###############
## Choose subsets with number of locations in each year >=4, duration >=10, and number of years >2

load("data/BioTIME/biotime_2021.RDATA")

bt_loc <- bt %>% dplyr::select(STUDY_ID, YEAR, LATITUDE, LONGITUDE) %>% distinct()

# remove years with number of locations < 4.
bt_loc <-  bt_loc %>% 
  group_by(STUDY_ID, YEAR) %>%
  mutate(n_loc = n_distinct(LATITUDE, LONGITUDE)) %>%
  ungroup() %>%
  filter(n_loc >= 4) %>%
  dplyr::select(-n_loc)

# calculate number of locations, number of years and duration
# and remove studies with number of years < 2 and duration <10 years
meta_year <- bt_loc %>%
  group_by(STUDY_ID, YEAR) %>%
  summarise(n_loc = n_distinct(LATITUDE, LONGITUDE, na.rm = TRUE)) %>%
  group_by(STUDY_ID) %>%
  mutate(t_loc = sum(n_loc), 
         mean_loc = round(mean(n_loc),1), 
         min_loc = min(n_loc), 
         max_loc = max(n_loc),
         n_years = n_distinct(YEAR, na.rm = TRUE),
         duration = max(YEAR) - min(YEAR) +1) %>% 
  ungroup() %>%
  filter(n_years >=2 & duration >=10)

meta <- meta_year %>% 
  dplyr::select(-c(YEAR, n_loc)) %>% 
  distinct()


# frequency distribution of number of time points, duration, locations in a study
ggplot(meta) + geom_histogram(aes(n_years)) +  theme(axis.text=element_text(size=14), axis.title=element_text(size=20))

ggplot(meta) + geom_histogram(aes(duration)) + theme(axis.text=element_text(size=14), axis.title=element_text(size=20))

ggplot(meta) + geom_histogram(aes(mean_loc)) + scale_x_log10() + theme(axis.text=element_text(size=14), axis.title=element_text(size=20))


# remove studies with number of years < 2 and duration <10 years
bt_4loc_10yr <- bt_loc %>% 
  unite(col = location, LATITUDE, LONGITUDE, remove = FALSE) %>%
  unite(col = study_year, STUDY_ID, YEAR, remove = FALSE) %>%
  filter(study_year %in% (meta_year %>% unite(col = study_year, STUDY_ID, YEAR) %>% pull(study_year))) %>%
  dplyr::select(-study_year)



############
## the number and locations of sites are somewhat different across years
# match sites across years based on grid-cells to keep similar configurations of sites

# Calculate minimum, maximum latitude and longitude, and spans of longitude and longitude
meta_spat <- bt_4loc_10yr %>% 
  group_by(STUDY_ID) %>%
  mutate(min_LATITUDE =  min(LATITUDE),
         max_LATITUDE =  max(LATITUDE),
         min_LONGITUDE =  min(LONGITUDE),
         max_LONGITUDE =  max(LONGITUDE)) %>%
  # If locations cross 180 degree, update the minimum longitude to calculate spans of longitude
  mutate(LONGITUDE = ifelse((max_LONGITUDE - min_LONGITUDE) >180 & LONGITUDE < 0, LONGITUDE + 360, LONGITUDE)) %>%
  mutate(span_LATITUDE = max(LATITUDE) - min(LATITUDE),
         span_LONGITUDE = max(LONGITUDE) - min(LONGITUDE)) %>%
  ungroup() %>%
  distinct(STUDY_ID, min_LATITUDE, max_LATITUDE, min_LONGITUDE, max_LONGITUDE, span_LATITUDE, span_LONGITUDE)


# Check spatial distributions of samples in grid-cells.
i = 2
study_spat <- meta_spat[i, c("STUDY_ID", "min_LONGITUDE", "max_LONGITUDE", "min_LATITUDE", "max_LATITUDE", "span_LONGITUDE", "span_LATITUDE")] %>% as.data.frame()
# set the resolution as the mean 1/5 spans of longitude and longitude
res <- mean(c(study_spat[, 6]/5, study_spat[, 7]/5)) 
ras <- raster(xmn = study_spat[1, 2] - res, xmx = study_spat[1, 3] + res, 
              ymn = study_spat[1, 4] - res, ymx = study_spat[1, 5] + res, 
              resolution = res)
values(ras) <- 1:length(ras)
plot(ras)
points(bt_4loc_10yr %>% filter(STUDY_ID == study_spat[1, 1]) %>% dplyr::select(LONGITUDE, LATITUDE))


# add cells for each study. cells are only comparable within studies
bt_4loc_10yr <- bt_4loc_10yr %>% mutate(cell = NA)
for(i in 1:nrow(meta_spat)){
  study_spat <- meta_spat[i, c("STUDY_ID", "min_LONGITUDE", "max_LONGITUDE", "min_LATITUDE", "max_LATITUDE", "span_LONGITUDE", "span_LATITUDE")] %>% as.data.frame()
  # set the resolution as the mean 1/5 spans of longitude and longitude
  res <- mean(c(study_spat[, 6]/5, study_spat[, 7]/5))
  ras <- raster(xmn = study_spat[1, 2] - res, xmx = study_spat[1, 3] + res, 
                ymn = study_spat[1, 4] - res, ymx = study_spat[1, 5] + res, 
                resolution = res)
  values(ras) <- 1:length(ras)
  id <- bt_4loc_10yr$STUDY_ID == study_spat[1, 1]
  bt_4loc_10yr$cell[id] <- cellFromXY(ras, bt_4loc_10yr[id, c("LONGITUDE", "LATITUDE")] %>% as.data.frame())
}

# filter the dataset to remove rare cells and years with small extent
bt_4loc_10yr_filtered <- NULL
bt_years_max_loc <- NULL
for(i in 1:nrow(meta_spat)){
  # perform loop for each study
  study <- bt_4loc_10yr %>% filter(STUDY_ID == meta_spat$STUDY_ID[i]) # %>% distinct(year, sample)
  
  # calculate number of locations in each cell of each year: rows are years, columns are cells, elements are number of locations
  year_cell <- as.matrix(xtabs( ~ YEAR + cell, data = study[, c("YEAR","cell")], sparse=TRUE)) 
  
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
    filter(cell %in% colnames(year_cell) & YEAR %in% rownames(year_cell))
  
  bt_years_max_loc <- bind_rows(bt_years_max_loc , bind_cols(STUDY_ID = meta_spat$STUDY_ID[i], max_co_loc[1, ]))
  bt_4loc_10yr_filtered  <- bind_rows(bt_4loc_10yr_filtered , rare_study)
}

ggplot(data = rare_study, aes(LONGITUDE, LATITUDE)) + 
  facet_wrap(~YEAR) +
  geom_point(size = 1, alpha = 0.7) + 
  coord_fixed()


# remove years with < 4 locations and studies with number of years < 2 and duration < 10 years
bt_4loc_10yr_filtered <- bt_4loc_10yr_filtered %>% 
  group_by(STUDY_ID, YEAR) %>%
  mutate(n_loc = n_distinct(LATITUDE, LONGITUDE)) %>% 
  filter(n_loc >= 4) %>%
  group_by(STUDY_ID) %>%
  mutate(n_years = n_distinct(YEAR, na.rm = TRUE),
         duration = max(YEAR) - min(YEAR) + 1) %>%
  filter(n_years >= 2 & duration >= 10) %>%
  ungroup() %>%
  dplyr::select(-c(n_loc, n_years, duration))

# add the column "keep" to distinguish the locations that should be kept or removed
bt_4loc_10yr <- bt_4loc_10yr %>% 
  left_join(bt_4loc_10yr_filtered %>% mutate(keep = "yes")) %>%
  mutate(keep = ifelse(is.na(keep), "no", keep))

bt_years_max_loc <- bt_years_max_loc %>% 
  left_join(bt_meta %>% dplyr::select(STUDY_ID, REALM, TAXA))


# plot distributions of locations and indicate which locations will be removed
pdf('data/BioTIME/BioTIME_4locations_10years_filtered.pdf', width = 12, height = 10)
for(i in 1:nrow(bt_years_max_loc)){
  study <- bt_4loc_10yr %>% filter(STUDY_ID %in% bt_years_max_loc$STUDY_ID[i])
  #study_title <- meta$STUDY_ID[i]
  study_title <- paste(bt_years_max_loc[i,1], 
      "_year1=", bt_years_max_loc[i,2], "_year2=", bt_years_max_loc[i,3], 
      "_max.shared.samples=", bt_years_max_loc[i,4], 
      "_duration=", bt_years_max_loc[i,5], 
      "_realm=", bt_years_max_loc[i,6], "_taxa=", bt_years_max_loc[i,7], sep="")
  
  p <- ggplot(data = study, aes(LONGITUDE, LATITUDE)) + 
    facet_wrap(~YEAR) +
    geom_point( aes(colour = keep), size = 1, alpha = 0.7) + 
    labs(title = study_title) +
    theme(legend.position = "top", legend.text = element_text(size=12)) + 
    coord_fixed() +
    scale_color_manual(values = c("yes" = "deepskyblue", "no" = "coral"))

  print(p)
}
dev.off()


# occurrence data for filtered locations
bt_filtered <- bt_4loc_10yr_filtered %>% 
  inner_join(bt %>% dplyr::select(-X1), by = c("STUDY_ID",  "YEAR", "LATITUDE", "LONGITUDE")) %>%
  left_join(bt_meta, by = "STUDY_ID") %>%
  unite(col = "location", LATITUDE, LONGITUDE, remove = FALSE)

# Save data
save(bt_filtered, file = "data/BioTIME/BioTIME_4locations_10years_filtered.RDATA")


# check how many samples for each location
location_nsamples <- bt_filtered %>% 
  dplyr::select(STUDY_ID, YEAR, location, SAMPLE_DESC) %>% 
  distinct() %>%
  group_by(STUDY_ID, YEAR, location) %>%
  summarise(n_samp = n_distinct(SAMPLE_DESC)) %>%
  ungroup()

# ~10% of locations have more than one sample
table(location_nsamples$n_samp > 1)

# the minimum number of samples across locations within studies
study_location_nsamples <- location_nsamples %>% 
  group_by(STUDY_ID) %>% 
  summarise(n_location = n_distinct(location),
            min_samp = min(n_samp),
            max_samp = max(n_samp),
            median_samp = median(n_samp),
            mean_samp = mean(n_samp)) 

# 9 studies have multiple samples at all locations
table(study_location_nsamples$min_samp > 1)

# keep minimum number of samples across locations within studies for each location
set.seed(10)
bt_filtered_location <- bt_filtered %>% 
  dplyr::select(STUDY_ID, YEAR, location, SAMPLE_DESC) %>% 
  distinct() %>%
  group_by(STUDY_ID, YEAR, location) %>% 
  mutate(n_samp = n_distinct(SAMPLE_DESC)) %>%
  group_by(STUDY_ID) %>% 
  mutate(min_samp = min(n_samp)) %>%
  dplyr::select(-n_samp) %>%
  group_by(STUDY_ID, YEAR, location) %>% 
  sample_n(size = unique(min_samp))

bt_filtered <- bt_filtered %>% inner_join(bt_filtered_location %>% dplyr::select(-min_samp))

save(bt_filtered, file = "data/BioTIME/BioTIME_4locations_10years_filtered_final.RDATA")


# the StudyID used
bt_studyID_filtered <- bt_filtered %>% distinct(STUDY_ID, HABITAT, TAXA, ORGANISMS, TITLE )
write.csv(bt_studyID_filtered, file = "data/BioTIME/bt_studyID_filtered.csv", row.names = FALSE)

