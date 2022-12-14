## filter InsectChange database to choose subsets of datasets with number of locations in each year >=4, duration >=10, 
## and matching locations across years to keep locations in similar configuration across years

rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(raster)
library(vegan)
library(reshape2)

it <- readRDS("data/Role_insects/Raw data insect metacommunities.Rdata")

# the additional insect studies after the publication of InsectChange
load("data/Role_insects/mosquito_enc_standardized.RDATA")


# to combine these new studies with the main insect studies, change the column names
mosquito_enc_std <- mosquito_enc_std %>% 
  rename(Datasource_name = study, Plot_ID = plot, Latitude = Latitudes, Longitude = Longitudes, 
         Year = year, species = Species, Reference = Citations, Datasource_nameREDUNDANT = Tag , Number = N, Realm = realm) %>%
  mutate(Abundance.Biomass = "A") %>%
  relocate(Datasource_ID)


# use only records that were identified as species, and filtered out studies with < 4 locations
it_4loc <- it %>% 
  as_tibble() %>% 
  filter(Level == "Species") %>% 
  mutate(species = paste(Genus, Species, sep = " ")) %>% 
  # add the new studies
  mutate(Datasource_ID = as.character(Datasource_ID),
         Plot_ID = as.character(Plot_ID)) %>%
  bind_rows(mosquito_enc_std %>% 
              mutate(Datasource_ID = as.character(Datasource_ID))) %>%
  group_by(Datasource_ID) %>%
  mutate(n_loc = n_distinct(Plot_ID)) %>% 
  ungroup() %>%
  filter(n_loc > 3) %>%
  dplyr::select(!n_loc)

# the study 1367 survey in 1980 and 1981, and resurvey in 2008 and 2009; change the years as two periods
id1 <- it_4loc$Datasource_ID == "1367" & it_4loc$Year == 1981
id2 <- it_4loc$Datasource_ID == "1367" & it_4loc$Year == 2009
it_4loc$Year[id1] <- 1980
it_4loc$Year[id2] <- 2008


# check whether different plots have same coordinates.
# only few plots have same coordinates with the others; use the Plot_ID to indicate locations (sites)
it_4loc %>% 
  distinct(Datasource_ID, Plot_ID, Longitude, Latitude) %>%
  filter(duplicated(.[c("Longitude", "Latitude")]))

# 2 studies (ID = 1367, 1408) provide only regional central coordinates
it_4loc %>% 
  group_by(Datasource_ID) %>%
  summarise(n_loc = n_distinct(Longitude, Latitude)) %>%
  filter(n_loc == 1)


############
## the number and locations of sites are somewhat different across years
# match sites across years based on grid-cells to keep similar configurations of sites

# Calculate minimum, maximum Latitude and Longitude, and spans of Longitude and Longitude
meta_spat <- it_4loc %>% 
  distinct(Datasource_ID, Plot_ID, Year, Latitude, Longitude) %>%
  group_by(Datasource_ID) %>%
  mutate(min_Latitude =  min(Latitude),
         max_Latitude =  max(Latitude),
         min_Longitude =  min(Longitude),
         max_Longitude =  max(Longitude)) %>%
  # If locations cross 180 degree, update the minimum Longitude to calculate spans of Longitude
  mutate(Longitude = ifelse((max_Longitude - min_Longitude) >180 & Longitude < 0, Longitude + 360, Longitude)) %>%
  mutate(span_Latitude = max(Latitude) - min(Latitude),
         span_Longitude = max(Longitude) - min(Longitude)) %>%
  ungroup() %>%
  distinct(Datasource_ID, min_Latitude, max_Latitude, min_Longitude, max_Longitude, span_Latitude, span_Longitude)


# Check spatial distributions of samples in grid-cells.
i = 2
study_spat <- meta_spat[i, c("Datasource_ID", "min_Longitude", "max_Longitude", "min_Latitude", "max_Latitude", "span_Longitude", "span_Latitude")] %>% as.data.frame()
# set the resolution as the mean 1/5 spans of Longitude and Longitude
res <- mean(c(study_spat[, 6]/5, study_spat[, 7]/5)) 
ras <- raster(xmn = study_spat[1, 2] - res, xmx = study_spat[1, 3] + res, 
              ymn = study_spat[1, 4] - res, ymx = study_spat[1, 5] + res, 
              resolution = res)
values(ras) <- 1:length(ras)
plot(ras)
points(it_4loc %>% 
         distinct(Datasource_ID, Plot_ID, Latitude, Longitude) %>% 
         filter(Datasource_ID == study_spat[1, 1]) %>% dplyr::select(Longitude, Latitude))


# add cells for each study. cells are only comparable within studies
it_4loc <- it_4loc %>% mutate(cell = NA)
for(i in 1:nrow(meta_spat)){
  study_spat <- meta_spat[i, c("Datasource_ID", "min_Longitude", "max_Longitude", "min_Latitude", "max_Latitude", "span_Longitude", "span_Latitude")] %>% as.data.frame()
  
  if(study_spat$span_Longitude > 0){
    # set the resolution as the mean 1/5 spans of Longitude and Longitude
    res <- mean(c(study_spat[, 6]/5, study_spat[, 7]/5))
    ras <- raster(xmn = study_spat[1, 2] - res, xmx = study_spat[1, 3] + res, 
                  ymn = study_spat[1, 4] - res, ymx = study_spat[1, 5] + res, 
                  resolution = res)
    values(ras) <- 1:length(ras)
    id <- it_4loc$Datasource_ID == study_spat[1, 1]
    it_4loc$cell[id] <- cellFromXY(ras, it_4loc[id, c("Longitude", "Latitude")] %>% as.data.frame())
  }
  
  # two studies only have regional central coordinates. Use the PLOT_ID as the cells 
  if(study_spat$span_Longitude == 0){
    id <- it_4loc$Datasource_ID == study_spat[1, 1]
    it_4loc$cell[id] <- it_4loc$Plot_ID[id]
  }
}


# filter the dataset to remove rare cells and years with small extent
it_4loc_filtered <- NULL
it_years_max_loc <- NULL
for(i in 1:nrow(meta_spat)){
  # perform loop for each study
  study <- it_4loc %>% 
    filter(Datasource_ID == meta_spat$Datasource_ID[i])
  
  # calculate number of locations in each cell of each year: rows are years, columns are cells, elements are number of locations
  year_cell <- as.matrix(xtabs( ~ Year + cell, data = unique(study[, c("Plot_ID","Year","cell")]), sparse=TRUE)) 
  
  
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
    filter(cell %in% colnames(year_cell) & Year %in% rownames(year_cell))
  
  it_years_max_loc <- bind_rows(it_years_max_loc , bind_cols(Datasource_ID = meta_spat$Datasource_ID[i], max_co_loc[1, ]))
  it_4loc_filtered  <- bind_rows(it_4loc_filtered, rare_study)
}

ggplot(data = rare_study, aes(Longitude, Latitude)) + 
  facet_wrap(~Year) +
  geom_point(size = 1, alpha = 0.7) + 
  coord_fixed()


# keep only years with at least 4 plots,
# and keep studies with at least 2 time points and duration >10 years
it_4loc_filtered <- it_4loc_filtered %>%
  group_by(Datasource_ID, Year) %>%
  mutate(n_samp = n_distinct(Plot_ID)) %>%
  filter(n_samp > 3) %>%
  group_by(Datasource_ID) %>%
  mutate(all_samp = n_distinct(Plot_ID),
         min_samp = min(n_samp),
         n_years = n_distinct(Year),
         duration = max(Year) - min(Year) +1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)

# check how many studies and their attributes
it_studies <- it_4loc_filtered %>% distinct(Datasource_ID, all_samp, min_samp, n_years, duration) # 23 studies
table(it_studies$all_samp)
table(it_studies$min_samp)
table(it_studies$n_years)  
table(it_studies$duration)

# save the filtered data
it_filtered <- it_4loc_filtered %>% dplyr::select(-(n_samp:duration)) 
save(it_filtered, file = "data/Role_insects/inset_metacommunities_filtered.RDATA")


# locations of the filtered dataset
it_loc_filtered <- it_4loc_filtered %>% 
  distinct(Datasource_ID, Plot_ID, Year, Latitude, Longitude )

# locations of the full dataset
it_loc <- it_4loc %>% 
  distinct(Datasource_ID, Plot_ID, Year, Latitude, Longitude ) %>% 
  left_join(it_loc_filtered %>%
              mutate(keep = "yes")) %>%
  mutate(keep = ifelse(is.na(keep), "no", keep))

# plot distributions of samples and indicate which records will be removed
pdf('data/Role_insects/inset_metacommunities_filtered.pdf', width = 12, height = 10)
id_study <- unique(it_loc$Datasource_ID)
for(i in 1:length(id_study)){
  study <- it_loc %>% 
    filter(Datasource_ID %in% id_study[i])
  
  p <- ggplot(data = study, aes(Longitude, Latitude)) + 
    facet_wrap(~Year) +
    geom_point( aes(colour = keep), size = 1, alpha = 0.7) + 
    labs(title = id_study[i]) +
    theme(legend.position = "top", legend.text = element_text(size=12)) + 
    coord_fixed() +
    scale_color_manual(values=c("yes" = "deepskyblue", "no" = "coral"))
  
  print(p)
}
dev.off()

