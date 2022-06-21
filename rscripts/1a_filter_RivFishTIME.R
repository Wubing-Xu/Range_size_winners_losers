## filter RivFishTIME database to choose subsets of datasets with number of locations in each year >=4, duration >=10, 
## and matching locations across years to keep locations in similar configuration across years


rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(raster)
library(vegan)
library(reshape2)


ft <- read_csv("data/RivFishTIME/1873_10_1873_2_RivFishTIME_SurveyTable.csv")
meta <- read_csv('data/RivFishTIME/1873_10_1873_2_RivFishTIME_TimeseriesTable.csv')[,-13]


# studies with at least 4 time series
study_4loc <- meta %>%
  group_by(SourceID) %>%
  summarise(n_loc = n_distinct(TimeSeriesID),
            n_protocol = n_distinct(Protocol)) %>% 
  ungroup() %>% 
  filter(n_loc > 3)

# 10 studies with multiple protocols
study_mp <- study_4loc  %>% 
  filter(n_protocol > 1)

# Choose the protocol with the maximum number of time series
study_mp_selected <- meta %>% 
  filter(SourceID %in% pull(study_mp, SourceID)) %>%
  group_by(SourceID, Protocol) %>%
  summarise(n_loc = n_distinct(TimeSeriesID)) %>%
  filter(n_loc > 3) %>%
  group_by(SourceID) %>%
  filter(n_loc == max(n_loc)) %>%
  ungroup() %>% 
  distinct(SourceID, .keep_all = TRUE) #one study have 2 protocol with same number of time series. Keep one of them.

# the studies with at least 4 time series in the same protocol 
meta_4loc <- bind_rows(meta %>% filter(SourceID %in% 
                                         (study_4loc %>% filter(n_protocol == 1) %>% pull(SourceID))),
                       meta %>% inner_join(study_mp_selected %>% dplyr::select(-n_loc)))

# check number of sites and time series for each study, and whether they are equal
meta_4loc_nsites <- meta_4loc %>% 
  group_by(SourceID) %>% 
  summarise(n_site = n_distinct(SiteID),
            n_ts = n_distinct(TimeSeriesID))
table(meta_4loc_nsites$n_site == meta_4loc_nsites$n_ts) # all 39 studies is TRUE, meaning each time series has unique site

# communities for the selected time series
ft_4loc <- inner_join(ft, meta_4loc, by = "TimeSeriesID") %>% 
  relocate(SourceID)


############
## the number and locations of sites are somewhat different across years
# match sites across years based on grid-cells to keep similar configurations of sites

# Calculate minimum, maximum Latitude and Longitude, and spans of Longitude and Longitude
meta_spat <- ft_4loc %>% 
  distinct(SourceID, TimeSeriesID, Year, Latitude, Longitude) %>%
  group_by(SourceID) %>%
  mutate(min_Latitude =  min(Latitude),
         max_Latitude =  max(Latitude),
         min_Longitude =  min(Longitude),
         max_Longitude =  max(Longitude)) %>%
  # If locations cross 180 degree, update the minimum Longitude to calculate spans of Longitude
  mutate(Longitude = ifelse((max_Longitude - min_Longitude) >180 & Longitude < 0, Longitude + 360, Longitude)) %>%
  mutate(span_Latitude = max(Latitude) - min(Latitude),
         span_Longitude = max(Longitude) - min(Longitude)) %>%
  ungroup() %>%
  distinct(SourceID, min_Latitude, max_Latitude, min_Longitude, max_Longitude, span_Latitude, span_Longitude)


# Check spatial distributions of samples in grid-cells.
i = 2
study_spat <- meta_spat[i, c("SourceID", "min_Longitude", "max_Longitude", "min_Latitude", "max_Latitude", "span_Longitude", "span_Latitude")] %>% as.data.frame()
# set the resolution as the mean 1/5 spans of Longitude and Longitude
res <- mean(c(study_spat[, 6]/5, study_spat[, 7]/5)) 
ras <- raster(xmn = study_spat[1, 2] - res, xmx = study_spat[1, 3] + res, 
              ymn = study_spat[1, 4] - res, ymx = study_spat[1, 5] + res, 
              resolution = res)
values(ras) <- 1:length(ras)
plot(ras)
points(ft_4loc %>% 
         distinct(SourceID, TimeSeriesID, Latitude, Longitude) %>% 
         filter(SourceID == study_spat[1, 1]) %>% dplyr::select(Longitude, Latitude))


# add cells for each study. cells are only comparable within studies
ft_4loc <- ft_4loc %>% mutate(cell = NA)
for(i in 1:nrow(meta_spat)){
  study_spat <- meta_spat[i, c("SourceID", "min_Longitude", "max_Longitude", "min_Latitude", "max_Latitude", "span_Longitude", "span_Latitude")] %>% as.data.frame()
  # set the resolution as the mean 1/5 spans of Longitude and Longitude
  res <- mean(c(study_spat[, 6]/5, study_spat[, 7]/5))
  ras <- raster(xmn = study_spat[1, 2] - res, xmx = study_spat[1, 3] + res, 
                ymn = study_spat[1, 4] - res, ymx = study_spat[1, 5] + res, 
                resolution = res)
  values(ras) <- 1:length(ras)
  id <- ft_4loc$SourceID == study_spat[1, 1]
  ft_4loc$cell[id] <- cellFromXY(ras, ft_4loc[id, c("Longitude", "Latitude")] %>% as.data.frame())
}


# filter the dataset to remove rare cells and years with small extent
ft_4loc_filtered <- NULL
ft_years_max_loc <- NULL
for(i in 1:nrow(meta_spat)){
  # perform loop for each study
  study <- ft_4loc %>% 
    filter(SourceID == meta_spat$SourceID[i])
  
  # calculate number of locations in each cell of each year: rows are years, columns are cells, elements are number of locations
  year_cell <- as.matrix(xtabs( ~ Year + cell, data = unique(study[, c("TimeSeriesID","Year","cell")]), sparse=TRUE)) 
  
  
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
  
  ft_years_max_loc <- bind_rows(ft_years_max_loc , bind_cols(SourceID = meta_spat$SourceID[i], max_co_loc[1, ]))
  ft_4loc_filtered  <- bind_rows(ft_4loc_filtered, rare_study)
}

ggplot(data = rare_study, aes(Longitude, Latitude)) + 
  facet_wrap(~Year) +
  geom_point(size = 1, alpha = 0.7) + 
  coord_fixed()


# keep only years with at least 4 sites,
# and keep studies with at least 2 time points and duration >10 years
ft_4loc_filtered <- ft_4loc_filtered %>%
  group_by(SourceID, Year) %>%
  mutate(n_site = n_distinct(TimeSeriesID)) %>%
  filter(n_site > 3) %>%
  group_by(SourceID) %>%
  mutate(all_site = n_distinct(TimeSeriesID ),
         min_site = min(n_site),
         n_years = n_distinct(Year),
         duration = max(Year) - min(Year) +1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)


# check how many studies and their attributes
ft_studies <- ft_4loc_filtered %>% distinct(SourceID, all_site, min_site, n_years, duration) #36 studies
table(ft_studies$all_site == ft_studies$min_site) 
table(ft_studies$all_site)
table(ft_studies$min_site) # 23 studies >= 10
table(ft_studies$n_years)  # 21 studies >= 10
table(ft_studies$duration)


# add the column "keep" to distinguish the records that should be kept or removed
ft_4loc <- ft_4loc %>% 
  left_join(ft_4loc_filtered %>% 
              dplyr::select(- (n_site:duration)) %>%
              mutate(keep = "yes")) %>%
  mutate(keep = ifelse(is.na(keep), "no", keep))


# plot distributions of sites and indicate which records will be removed
pdf('data/RivFishTIME/RivFishTIME_4locations_filtered.pdf', width = 12, height = 10)
id_study <- unique(ft_4loc$SourceID)
for(i in 1:length(id_study)){
  study <- ft_4loc %>% 
    filter(SourceID %in% id_study[i]) %>% 
    distinct(SourceID, TimeSeriesID, Year, Latitude, Longitude, keep)
  
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


# check how many time series have different number of surveys across years
n_survey <- ft_4loc_filtered %>% 
  dplyr::select(-(n_site:duration)) %>%
  group_by(TimeSeriesID, Year) %>%
  summarise(n_survey = n_distinct(SurveyID)) %>%
  group_by(TimeSeriesID) %>% 
  summarise(min_survey = min(n_survey),
         max_survey = max(n_survey))

# ~40% of time series have different number of surveys between years
table(n_survey$min_survey == n_survey$max_survey)

# only 44 time series have multiple surveys at all year
table(n_survey$min_survey > 1) 

# keep one survey for each year of each time series
one_survey_perYear <- ft_4loc_filtered %>% 
  dplyr::select(SourceID, TimeSeriesID, Year, SurveyID) %>%
  distinct(SourceID, TimeSeriesID, Year, .keep_all=TRUE)
ft_4loc_filtered <-  ft_4loc_filtered %>% inner_join(one_survey_perYear)

# save the filtered data
ft_filtered <- ft_4loc_filtered %>% dplyr::select(-(n_site:duration)) 
save(ft_filtered, file = "data/RivFishTIME/RivFishTIME_4sites_10years_filtered.RDATA")

