## calculate occupancy in the first and last year and
## the average ones for early and late periods that pool multiple years
## and evaluate the uncertainty by randomly selecting same number of samples for both periods 200 times 


rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

library(tidyverse)


load("data/Combined_assemblages_20220208.RDATA")
load("data/Species_rangesize.RDATA")


# A function to calculate occupancy
get_occupancy <- function(data, fillspecies = FALSE, trimsamples = TRUE, occ_full = TRUE, occ_rarefy = TRUE, resamples = 1){
  # data should have five columns: study_id, year, sample and species, has_range
  colnames(data) <- c("study","year","sample","species")
  
  nsamples <- data %>%
    group_by(study, year) %>%
    summarise(n_samp = n_distinct(sample)) %>%
    group_by(study) %>%
    mutate(mean_samp = mean(n_samp),
           min_samp = min(n_samp)) %>% 
    ungroup()
  
  # delete years with number of locations < half of the mean across years
  if(trimsamples){
    nsamples <- nsamples %>% 
      filter(n_samp >= 0.5*mean_samp) %>%
      group_by(study) %>%
      mutate(mean_samp = mean(n_samp),
             min_samp = min(n_samp)) %>%
      ungroup()
  }
  
   data <- inner_join(data, nsamples) %>%
     group_by(study) %>% 
     filter(n_distinct(year) >= 2) %>% 
     ungroup()
   
   study_meta <- data %>%
     group_by(study) %>%
     summarise(
       sprich = n_distinct(species),
       n_years = n_distinct(year),
       start_year = min(year),
       end_year = max(year),
       duration = max(year) - min(year) + 1
     )

   n_occ_full <- NULL
   n_occ_rare <- NULL
   
   if(occ_full){
     n_occ_full <- data %>%
       group_by(study, year, species) %>% 
       summarise(n_occ = n_distinct(sample)) %>% 
       ungroup()
     
     if(fillspecies){
       n_occ_full <- group_by(n_occ_full, study) %>%
         complete(year, species, fill = list(n_occ = 0)) %>%
         ungroup()
     }
     n_occ_full <- n_occ_full %>% mutate(resample = 0)
   }
   
   if(occ_rarefy){
     data_nest <- data %>% 
       group_by(study, year, sample) %>%
       nest(data = species) %>%
       ungroup()
     
     for(i in 1:resamples){
       rare_data <- tibble()
       for(j in 1:length(unique(data_nest$study))){
          study <- data_nest %>% 
           filter(study == unique(data_nest$study)[j]) %>%
           group_by(study, year) %>% 
           sample_n(size = unique(min_samp))
         rare_data <- bind_rows(rare_data, study)
       }
       n_occ_rarefied <- rare_data %>% unnest(data) %>%
         group_by(study, year, species) %>% 
         summarise(n_occ = n_distinct(sample)) %>% 
         ungroup()
       
       if(fillspecies){
         n_occ_rarefied <- group_by(n_occ_rarefied, study) %>%
           complete(year, species, fill = list(n_occ = 0)) %>%
           ungroup()
       }
       n_occ_rarefied <- n_occ_rarefied %>% mutate(resample = i)
       n_occ_rare <- bind_rows(n_occ_rare, n_occ_rarefied)
       print(i)
     }
   }
   n_occ <- bind_rows(n_occ_full, n_occ_rare)
   
   occupancy <- study_meta %>% 
     left_join(nsamples) %>% 
     left_join(n_occ) %>%
     mutate(occupancy = ifelse(resample ==0, n_occ/n_samp, n_occ/min_samp)) %>%
     relocate(resample) %>%
     arrange(resample, study, year)
   
   return(occupancy)
}


# A function to calculate occupancy change
get_occup_change <- function(x){
  # x should have 3 columns, year, n_occ, occupancy
  x <- x[order(x[,1]),] %>% data.frame()
  presence <- ifelse(x$occupancy >0, 1, 0)
  
  nocc_first <- x[1, 2]
  nocc_last <- x[2, 2]
  occup_first <- x[1, 3]
  occup_last <- x[2, 3]
  occup_change <- occup_last - occup_first
  
  if(presence[1] == 1 & presence[2] == 1) dynamic = "persistent"
  if(presence[1] == 0 & presence[2] == 1) dynamic = "colonization"
  if(presence[1] == 1 & presence[2] == 0) dynamic = "extinction"
  if(presence[1] == 0 & presence[2] == 0) dynamic = "absent"
  
  res <- data.frame(dynamic, occup_first, occup_last, occup_change, nocc_first, nocc_last)
  return(res)
}


############
# calculate occupancy

#######
# Use the first and last years to calculate occupancy
dat_occupancy_2yr <- dat %>% 
  filter(period %in% c("first", "last")) %>% 
  dplyr::select(study, year, sample, specieskey)

# calculate occupancy
occupancy_2yr <- get_occupancy(data = dat_occupancy_2yr,
                           fillspecies = TRUE, trimsamples = FALSE, occ_full = FALSE, occ_rarefy = TRUE, resamples = 1)

occupancy_2yr <- occupancy_2yr %>% 
  rename(specieskey = species) %>%
  left_join(dat %>% distinct(study, database, studyID, study_name, year, period), by =c("study", "year")) %>%
  left_join(dat %>% distinct(species, specieskey), by =c("specieskey")) %>%
  mutate(nsamp_used = min_samp) %>%
  relocate(database, studyID, study_name, .after = study) %>% 
  relocate(period, .after = year) %>% 
  relocate(species, .after = specieskey) %>% 
  relocate(nsamp_used, .after = min_samp)



#######
# Combine same number of years before and after the median of start and end of years to calculate occupancy

# to determine which years will be combined
years_period <- dat %>% 
  distinct(study, duration, year) %>%
  group_by(study) %>%
  mutate(median_year = (max(year) + min(year))/2,
         nyears_early_median =  sum(year < unique(median_year)),
         nyears_late_median =  sum(year > unique(median_year)),
         nyears_combined = min(nyears_early_median, nyears_late_median),
         period_early = year <= sort(year)[nyears_combined],
         period_late = year >= sort(year, decreasing = TRUE)[nyears_combined],
         period_combined = ifelse(period_early, "early", ifelse(period_late, "late", "intermediate"))) %>%
  ungroup()

years_period <- years_period %>% 
  group_by(study) %>%
  mutate(year_period = case_when(period_combined == "early" ~ min(year),
                                 period_combined == "late" ~ max(year),
                                 period_combined == "intermediate" ~ mean(year))) %>%
  group_by(study, period_combined) %>%
  mutate(year_mean_period = mean(year)) %>% 
  group_by(study) %>%
  mutate(duration_mean = max(year_mean_period) - min(year_mean_period) + 1) %>%
  ungroup()

years_period <- years_period %>%
  dplyr::select(study, year, median_year, nyears_combined, period_combined, year_period, year_mean_period, duration_mean)

# add information about year combined to meta data
dat_meta <- dat_meta %>%
  left_join(years_period %>% distinct(study, nyears_combined, duration_mean))


## to calculate occupancy
# only use the time points to be pooled to calculate occupancy
dat_occupancy_period <- dat %>% 
  left_join(years_period, by = c("study", "year")) %>% 
  filter(period_combined %in% c("early", "late")) %>% 
  dplyr::select(study, year, sample, specieskey) %>%
  distinct()

# calculate occupancy
occupancy <- get_occupancy(data = dat_occupancy_period,
                           fillspecies = TRUE, trimsamples = FALSE, occ_full = FALSE, occ_rarefy = TRUE, resamples = 1)

occupancy <- occupancy %>% 
  rename(specieskey = species) %>%
  left_join(dat %>% distinct(study, database, studyID, study_name), by =c("study")) %>% 
  left_join(years_period, by = c("study", "year")) %>%
  left_join(dat %>% distinct(species, specieskey), by =c("specieskey")) %>%
  mutate(nsamp_used = min_samp*nyears_combined) %>%
  dplyr::select(-c(year)) %>%
  relocate(database, studyID, .after = study) %>% 
  relocate(year_period, period_combined, .after = duration) %>% 
  relocate(species, .after = specieskey) %>% 
  relocate(nsamp_used, .after = min_samp)

## the mean occupancy for the early and late periods
mean_occupancy_period  <- occupancy  %>% 
  group_by(resample, study, specieskey, period_combined) %>%
  summarise(n_occ = sum(n_occ),
            occupancy = mean(occupancy)) %>%
  ungroup()

occupancy_period <- occupancy %>% 
  dplyr::select(-c(n_samp, n_occ, occupancy)) %>% 
  distinct() %>% 
  left_join(mean_occupancy_period, by = c("resample","study", "specieskey", "period_combined"))




############
# calculate occupancy change

## use the data in first and last years
# nest different years of the same species
occupancy_2yr_nest <- occupancy_2yr %>% 
  dplyr::select(-period, - n_samp) %>%
  group_by(resample, study, species) %>%
  nest(data = c(year, n_occ, occupancy)) %>%
  ungroup()

occupancy_change_2yr <- tibble()
for(i in 1:nrow(occupancy_2yr_nest)){
  res <- get_occup_change(x = occupancy_2yr_nest[i,]$data[[1]])
  occupancy_change_2yr <- bind_rows(occupancy_change_2yr, res)
}

occupancy_change_2yr <- occupancy_2yr_nest %>% 
  dplyr::select(-data) %>%
  bind_cols(occupancy_change_2yr) %>% 
  left_join(dat %>% distinct(specieskey, nocc_community, nocc_gbif, ratio_nocc_gbif), by = c("specieskey")) %>%
  left_join(spsuma[,c(2,14:19)], by = c("specieskey")) %>%
  filter(dynamic != "absent")

table(occupancy_change_2yr[, c("dynamic", "database")])


## use the data combined years as periods
# nest different years of the same species
occupancy_period_nest <- occupancy_period %>% 
  dplyr::select(-period_combined, - year_mean_period) %>%
  group_by(resample, study, species) %>%
  nest(data = c(year_period, n_occ, occupancy)) %>%
  ungroup()

occupancy_change_period <- tibble()
for(i in 1:nrow(occupancy_period_nest)){
  res <- get_occup_change(x = occupancy_period_nest[i,]$data[[1]])
  occupancy_change_period <- bind_rows(occupancy_change_period, res)
}

occupancy_change_period <- occupancy_period_nest %>% 
  dplyr::select(-data) %>%
  bind_cols(occupancy_change_period) %>% 
  left_join(dat %>% distinct(specieskey, nocc_community, nocc_gbif, ratio_nocc_gbif), by = c("specieskey")) %>%
  left_join(spsuma[,c(2,14:19)], by = c("specieskey")) %>%
  filter(dynamic != "absent")

table(occupancy_change_period[, c("dynamic", "database")])


# save occupancy 
save(occupancy_2yr, occupancy_period, occupancy_change_2yr, occupancy_change_period, years_period, dat_meta, 
     file = "intermediate_results/occuapncy.RDATA")



#########################
## For sensitive analyses: resample samples for 200 times

# calculate occupancy
occupancy_resample <- get_occupancy(data = dat_occupancy_period,
                                    fillspecies = TRUE, trimsamples = FALSE, occ_full = FALSE, occ_rarefy = TRUE, resamples = 200)

occupancy_resample <- occupancy_resample %>% 
  rename(specieskey = species) %>%
  left_join(dat %>% distinct(study, database, studyID), by =c("study")) %>%
  left_join(years_period, by = c("study", "year")) %>%
  left_join(dat %>% distinct(species, specieskey), by =c("specieskey")) %>%
  mutate(nsamp_used = min_samp*nyears_combined) %>%
  select(-c(year)) %>%
  relocate(database, studyID, .after = study) %>% 
  relocate(year_period, period_combined, .after = duration) %>% 
  relocate(species, .after = specieskey) %>% 
  relocate(nsamp_used, .after = min_samp)

## the mean occupancy for the early and late periods
mean_occupancy_resample_period  <- occupancy_resample  %>% 
  group_by(resample, study, specieskey, period_combined) %>%
  summarise(n_occ = sum(n_occ),
            occupancy = mean(occupancy)) %>%
  ungroup()

occupancy_resample <- occupancy_resample %>% 
  dplyr::select(-c(n_samp, n_occ, occupancy)) %>% 
  distinct() %>% 
  left_join(mean_occupancy_resample_period, by = c("resample","study", "specieskey", "period_combined"))


#### calculate occupancy change
occupancy_resample_occup <- occupancy_resample %>% 
  dplyr::select(resample, study, database, studyID, period_combined, nsamp_used, specieskey, occupancy) %>%
  group_by(resample, study, specieskey) %>% 
  pivot_wider(names_from = period_combined, names_prefix = "occup_", values_from = occupancy) %>%
  ungroup()

occupancy_resample_nocc <- occupancy_resample %>% 
  dplyr::select(resample, study, database, studyID, period_combined, nsamp_used, specieskey, n_occ) %>%
  group_by(resample, study, specieskey) %>% 
  pivot_wider(names_from = period_combined, names_prefix = "nocc_", values_from = n_occ) %>%
  ungroup()

occupancy_change_resample <- occupancy_resample_occup %>% 
  left_join(occupancy_resample_nocc) %>% 
  mutate(dynamic = ifelse(occup_early > 0 & occup_late > 0, "persistent", 
                          ifelse(occup_early > 0 & occup_late == 0, "extinction", "colonization")),
         occup_change = occup_late - occup_early) %>%
  rename(occup_first = occup_early, occup_last = occup_late, nocc_first = nocc_early, nocc_last = nocc_late) %>%
  relocate(dynamic, occup_first, occup_last, occup_change, nocc_first, nocc_last, .after = last_col())

# add range size
occupancy_change_resample <- occupancy_change_resample %>%
  left_join(spsuma[,c(2, 14)], by = "specieskey")


save(occupancy_change_resample, dat_meta, file = "intermediate_results/occuapncy_resample.RDATA")

