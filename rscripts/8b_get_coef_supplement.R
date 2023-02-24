## get the fixed and random effects of range size and predicted/fitted occupancy and occupancy from the models 
## with interaction between range size with realm-taxa, or realm-region, and the models used for sensitivity analyse.
## These estimates will be used to draw supplement figures


rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ggplot2", "brms", "tidybayes", "abind")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)

# input the observation data and models
load("models/data_input_to_models.RDATA")
load("models/brm_oc_aoo50.RDATA")
load("models/brm_oc_aoo100.RDATA")
load("models/brm_oc_ahull6.RDATA")
load("models/brm_oc_aoo10_rgbif.RDATA")
load("models/brm_oc_aoo10_sloc.RDATA")
load("models/brm_oc_aoo10_realm_region.RDATA")
load("models/brm_oc_aoo10_realm_taxa.RDATA")


#############
## get coefficients from models using range size defined as AOO in the resolution of 50 km

## global fixed effects of range size
brm_oc_aoo50_fixed <- fixef(brm_oc_aoo50, probs = c(0.025, 0.975)) %>%
  as_tibble() %>% 
  mutate(term = c("intercept", "sigma_intercept" ,"slope", "sigma_nsamp")) 

## random effects of rang size for each study
brm_oc_aoo50_coef <- coef(brm_oc_aoo50, probs = c(0.025, 0.975))[[1]]
brm_oc_aoo50_coef <- as_tibble(brm_oc_aoo50_coef) %>%
  dplyr::select(1:8) %>%
   mutate(study = rownames(brm_oc_aoo50_coef)) 
colnames(brm_oc_aoo50_coef)[1:8] <- c("estimate_intercept", "se_intercept", "Q2.5_intercept", "Q97.5_intercept", 
                                      "estimate_slope", "se_slope", "Q2.5_slope", "Q97.5_slope")

brm_oc_aoo50_coef <- brm_oc_aoo50_coef %>% 
  mutate(sig_slope = ifelse(Q2.5_slope < 0 & Q97.5_slope >0, "neutral", ifelse(Q97.5_slope <= 0, "negative", "positive")),
         sig_slope = factor(sig_slope, levels = c("negative", "positive", "neutral"))) %>%
  left_join(oc_period %>% group_by(study) %>%
              summarise(xmin = min(aoo50),
                        xmax = max(aoo50),
                        cl.xmin = min(cl.aoo50),
                        cl.xmax = max(cl.aoo50))) %>%
  left_join(dat_meta) %>%
  relocate(study, database, studyID, study_name)


#############
## get coefficients from models using range size defined as AOO in the resolution of 100 km

## global fixed effects of range size
brm_oc_aoo100_fixed <- fixef(brm_oc_aoo100,  probs = c(0.025, 0.975)) %>%
  as_tibble() %>%
  mutate(term = c("intercept", "sigma_intercept" ,"slope", "sigma_nsamp")) 

## random effects of rang size for each study
brm_oc_aoo100_coef <- coef(brm_oc_aoo100, probs = c(0.025, 0.975))[[1]]
brm_oc_aoo100_coef <- as_tibble(brm_oc_aoo100_coef) %>%
  dplyr::select(1:8) %>%
  mutate(study = rownames(brm_oc_aoo100_coef)) 
colnames(brm_oc_aoo100_coef)[1:8] <- c("estimate_intercept", "se_intercept", "Q2.5_intercept", "Q97.5_intercept", 
                                      "estimate_slope", "se_slope", "Q2.5_slope", "Q97.5_slope")

brm_oc_aoo100_coef <- brm_oc_aoo100_coef %>% 
  mutate(sig_slope = ifelse(Q2.5_slope < 0 & Q97.5_slope >0, "neutral", ifelse(Q97.5_slope <= 0, "negative", "positive")),
         sig_slope = factor(sig_slope, levels = c("negative", "positive", "neutral"))) %>%
  left_join(oc_period %>% group_by(study) %>%
              summarise(xmin = min(aoo100),
                        xmax = max(aoo100),
                        cl.xmin = min(cl.aoo100),
                        cl.xmax = max(cl.aoo100))) %>%
  left_join(dat_meta) %>%
  relocate(study, database, studyID, study_name)


#############
## get coefficients from models using range size defined as ahull hulls (alpha = 6)

## global fixed effects of range size
brm_oc_ahull6_fixed <- fixef(brm_oc_ahull6, probs = c(0.025, 0.975)) %>%
  as_tibble() %>%
  mutate(term = c("intercept", "sigma_intercept" ,"slope", "sigma_nsamp")) 

colnames(brm_oc_ahull6_fixed)[3:4] <- c("Q2.5", "Q97.5")


## random effects of rang size for each study
brm_oc_ahull6_coef <- coef(brm_oc_ahull6, robust = TRUE, probs = c(0.025, 0.975))[[1]]
brm_oc_ahull6_coef <- as_tibble(brm_oc_ahull6_coef) %>%
  dplyr::select(1:8) %>%
  mutate(study = rownames(brm_oc_ahull6_coef)) 
colnames(brm_oc_ahull6_coef)[1:8] <- c("estimate_intercept", "se_intercept", "Q2.5_intercept", "Q97.5_intercept", 
                                      "estimate_slope", "se_slope", "Q2.5_slope", "Q97.5_slope")

brm_oc_ahull6_coef <- brm_oc_ahull6_coef %>% 
  mutate(sig_slope = ifelse(Q2.5_slope < 0 & Q97.5_slope >0, "neutral", ifelse(Q97.5_slope <= 0, "negative", "positive")),
         sig_slope = factor(sig_slope, levels = c("negative", "positive", "neutral"))) %>%
  left_join(oc_period %>% group_by(study) %>%
              summarise(xmin = min(ahull6),
                        xmax = max(ahull6),
                        cl.xmin = min(cl.ahull6),
                        cl.xmax = max(cl.ahull6))) %>%
  left_join(dat_meta) %>%
  relocate(study, database, studyID, study_name)


#############
## get coefficients from models using assemblage data with same locations through years 

## global fixed effects of range size
brm_oc_aoo10_sloc_fixed <- fixef(brm_oc_aoo10_sloc, probs = c(0.025, 0.975)) %>%
  as_tibble() %>%
  mutate(term = c("intercept", "sigma_intercept" ,"slope", "sigma_nsamp")) 


## random effects of rang size for each study
brm_oc_aoo10_sloc_coef <- coef(brm_oc_aoo10_sloc, probs = c(0.025, 0.975))[[1]]
brm_oc_aoo10_sloc_coef <- as_tibble(brm_oc_aoo10_sloc_coef) %>%
  dplyr::select(1:8) %>%
  mutate(study = rownames(brm_oc_aoo10_sloc_coef)) 
colnames(brm_oc_aoo10_sloc_coef)[1:8] <- c("estimate_intercept", "se_intercept", "Q2.5_intercept", "Q97.5_intercept", 
                                      "estimate_slope", "se_slope", "Q2.5_slope", "Q97.5_slope")

brm_oc_aoo10_sloc_coef <- brm_oc_aoo10_sloc_coef %>% 
  mutate(sig_slope = ifelse(Q2.5_slope < 0 & Q97.5_slope >0, "neutral", ifelse(Q97.5_slope <= 0, "negative", "positive")),
         sig_slope = factor(sig_slope, levels = c("negative", "positive", "neutral"))) %>%
  left_join(oc_period %>% group_by(study) %>%
              summarise(xmin = min(aoo10),
                        xmax = max(aoo10),
                        cl.xmin = min(cl.aoo10),
                        cl.xmax = max(cl.aoo10))) %>%
  left_join(dat_meta) %>%
  relocate(study, database, studyID, study_name)


#############
## get coefficients from models using species that have 5-times of occurrences in GBIF than in community data

## global fixed effects of range size
brm_oc_aoo10_rgbif_fixed <- fixef(brm_oc_aoo10_rgbif, probs = c(0.025, 0.975)) %>%
  as_tibble() %>%
  mutate(term = c("intercept", "sigma_intercept" ,"slope", "sigma_nsamp")) 

## random effects of rang size for each study
brm_oc_aoo10_rgbif_coef <- coef(brm_oc_aoo10_rgbif, probs = c(0.025, 0.975))[[1]]
brm_oc_aoo10_rgbif_coef <- as_tibble(brm_oc_aoo10_rgbif_coef) %>%
  dplyr::select(1:8) %>%
  mutate(study = rownames(brm_oc_aoo10_rgbif_coef)) 
colnames(brm_oc_aoo10_rgbif_coef)[1:8] <- c("estimate_intercept", "se_intercept", "Q2.5_intercept", "Q97.5_intercept", 
                                      "estimate_slope", "se_slope", "Q2.5_slope", "Q97.5_slope")

brm_oc_aoo10_rgbif_coef <- brm_oc_aoo10_rgbif_coef %>% 
  mutate(sig_slope = ifelse(Q2.5_slope < 0 & Q97.5_slope >0, "neutral", ifelse(Q97.5_slope <= 0, "negative", "positive")),
         sig_slope = factor(sig_slope, levels = c("negative", "positive", "neutral"))) %>%
  left_join(dat_meta) %>%
  relocate(study, database, studyID, study_name)


#############
## get global fixed effects of range size for each combination of realm and region
# number of studies for each realm_region
nstudy_realm_taxa<- oc_period %>% 
  distinct(study, realm, taxon_new) %>%
  group_by(realm, taxon_new) %>%
  summarise(nstudy = n_distinct(study)) %>%
  ungroup() %>%
  mutate(realm_taxa = paste(realm, taxon_new, sep = "_"),
         realm_taxa = gsub(" ", "", realm_taxa))

brm_oc_aoo10_realm_taxa_fixed <- fixef(brm_oc_aoo10_realm_taxa, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm_taxa") %>% 
  filter(! grepl("sigma", realm_taxa)) %>%
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
  mutate(realm_taxa = gsub("realm_taxa", "", realm_taxa),
         realm_taxa = gsub(":cl.aoo10", "", realm_taxa)) %>%
  left_join(nstudy_realm_taxa) %>%
  mutate(realm_taxa = paste(realm, taxon_new, sep = "_"))

# add ranges of range sizes for plotting
brm_oc_aoo10_realm_taxa_line <- brm_oc_aoo10_realm_taxa_fixed %>% 
  dplyr::select(realm_taxa, realm, taxon_new, term, Estimate) %>%
  pivot_wider(names_from = "term", values_from = Estimate) %>% 
  left_join(oc_period %>% group_by(realm_taxa) %>%
              summarise(xmin = min(aoo10),
                        xmax = max(aoo10),
                        cl.xmin = min(cl.aoo10),
                        cl.xmax = max(cl.aoo10)))

## get fitted values of occupancy change for each observed range size from models including 
# interaction between range size and realm_region

brm_oc_aoo10_realm_taxa_fitted <- fitted(brm_oc_aoo10_realm_taxa, re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(oc_period %>% 
              dplyr::select(study, aoo10, occup_change_sqroot, nsamp_used))

brm_oc_aoo10_realm_taxa_fitted <- brm_oc_aoo10_realm_taxa_fitted %>% 
  rename(oc_sqroot_pred = Estimate, oc_sqroot_Q2.5 = Q2.5, oc_sqroot_Q97.5 = Q97.5) %>%
  left_join(dat_meta %>% distinct(study, database, studyID, taxon_new, taxon_final, realm, region)) %>% 
  distinct(aoo10, realm, taxon_new , .keep_all = TRUE)



#############
## get global fixed effects of range size for each combination of realm and region
# number of studies for each realm_region
nstudy_realm_region <- oc_period %>% 
  distinct(study, realm, region) %>%
  group_by(realm, region) %>%
  summarise(nstudy = n_distinct(study)) %>%
  ungroup() %>%
  mutate(realm_region = paste(realm, region, sep = "_"),
         realm_region = gsub(" ", "", realm_region))

brm_oc_aoo10_realm_region_fixed <- fixef(brm_oc_aoo10_realm_region, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm_region") %>% 
  filter(! grepl("sigma", realm_region)) %>%
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
  mutate(realm_region = gsub("realm_region", "", realm_region),
         realm_region = gsub(":cl.aoo10", "", realm_region)) %>%
  left_join(nstudy_realm_region) %>%
  mutate(realm_region = paste(realm, region, sep = "_"))

# add ranges of range sizes for plotting
brm_oc_aoo10_realm_region_line <- brm_oc_aoo10_realm_region_fixed %>% 
  dplyr::select(realm_region, realm, region, term, Estimate) %>%
  pivot_wider(names_from = "term", values_from = Estimate) %>% 
  left_join(oc_period %>% group_by(realm_region) %>%
              summarise(xmin = min(aoo10),
                        xmax = max(aoo10),
                        cl.xmin = min(cl.aoo10),
                        cl.xmax = max(cl.aoo10)))

## get fitted values of  occupancy change for each observed range size from models including 
# interaction between range size and realm_region

brm_oc_aoo10_realm_region_fitted <- fitted(brm_oc_aoo10_realm_region, re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(oc_period %>% 
              dplyr::select(study, aoo10, occup_change_sqroot, nsamp_used))

brm_oc_aoo10_realm_region_fitted <- brm_oc_aoo10_realm_region_fitted %>% 
  rename(oc_sqroot_pred = Estimate, oc_sqroot_Q2.5 = Q2.5, oc_sqroot_Q97.5 = Q97.5) %>%
  left_join(dat_meta %>% distinct(study, database, studyID, taxon_new, taxon_final, realm, region)) %>% 
  distinct(aoo10, realm, region , .keep_all = TRUE)
  

# save results
save(brm_oc_aoo50_fixed, brm_oc_aoo50_coef, 
     brm_oc_aoo100_fixed, brm_oc_aoo100_coef, 
     brm_oc_ahull6_fixed, brm_oc_ahull6_coef, 
     brm_oc_aoo10_sloc_fixed, brm_oc_aoo10_sloc_coef,
     brm_oc_aoo10_rgbif_fixed, brm_oc_aoo10_rgbif_coef, 
     brm_oc_aoo10_realm_taxa_fixed, brm_oc_aoo10_realm_taxa_line, brm_oc_aoo10_realm_taxa_fitted,
     brm_oc_aoo10_realm_region_fixed, brm_oc_aoo10_realm_region_line, brm_oc_aoo10_realm_region_fitted,
     file = "results/coefs_supplement_brms.RDATA")
