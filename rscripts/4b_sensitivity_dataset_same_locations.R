##################################
## to test of sensitivity of data filtering, generate a dataset keeping same number and locations of sites through years 

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy",
                  "IDIVTS01" = "H:/wubing/iDiv/homogenization_occupancy")
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


load("data/Combined_assemblages.RDATA")


# filter the dataset to maximum the number of resurveyed samples (same locations) for at least two years within a region
dat_sloc <- NULL
for(i in 1:length(unique(dat$study))){
  # perform loop for each study
  study <- dat %>% 
    filter(study == unique(dat$study)[i])
  
  # get a matrix with row are years and colunmns is samples
  year_sample <- as.matrix(xtabs(~ year + sample, data = unique(study[,c("year","sample")]), sparse=TRUE))
  
  # number of co-occurred samples between years
  co_samp <- as.matrix(designdist(year_sample, method = "J", terms= "binary"))
  
  #  find which two years (year-pair) have the maximum number of co-occurred samples and duration
  max_co_samp <- melt(co_samp) %>% 
    as_tibble() %>%
    set_names("year1","year2","n_samp") %>%
    mutate(year1 = as.numeric(year1),
           year2 = as.numeric(year2),
           duration = year2 - year1 + 1) %>%
    filter(duration > 9 & n_samp > 3)
  if(nrow(max_co_samp) == 0) {next}
  max_co_samp <-  filter(max_co_samp, n_samp >= 0.9*max(n_samp))
  max_co_samp <- filter(max_co_samp, duration == max(duration))
  max_co_samp <-  filter(max_co_samp, n_samp == max(n_samp))
  
  # keep samples that are shared in the two determined years
  samps_year1 <- year_sample[rownames(year_sample) %in% unlist(max_co_samp[1,1]),]
  samps_year2 <- year_sample[rownames(year_sample) %in% unlist(max_co_samp[1,2]),]
  samps_shared <- samps_year1 > 0 & samps_year2 > 0
  year_sample <- year_sample[,samps_shared, drop=FALSE]
  
  # Other years except the two priority years will be compared to the two years. Keep only the years have the same sites with the priority years
  id <- rowSums(year_sample) == ncol(year_sample)
  year_sample <-  year_sample[id, ]
  
  #  only keep the selected co-occurred samples and years 
  study_filtered <- study %>% 
    filter(sample %in% colnames(year_sample) & year %in% rownames(year_sample))
  
  dat_sloc <- bind_rows(dat_sloc, study_filtered)
}


# update metadata
# get the start and end years, and durations of studies
dat_sloc_year <- dat_sloc %>% 
  group_by(study) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

dat_sloc_nsamp <- dat_sloc %>% 
  group_by(study) %>%
  summarise(n_samp = n_distinct(sample))


dat_sloc_meta <- dat_meta %>%
  filter(study %in% dat_sloc$study) %>%
  dplyr::select(c(study:climate, cent_lat, cent_long)) %>%
  left_join(dat_sloc_year) %>%
  left_join(dat_sloc_nsamp)


save(dat_sloc, dat_sloc_meta,
     file = "data/Combined_assemblages_same_locations.RDATA")

