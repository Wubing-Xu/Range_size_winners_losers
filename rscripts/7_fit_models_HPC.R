## fit models to assess relationships between range size and occupancy change using hierarchical generalized linear models
# models were fitted based on Bayesian inference using the R package brms
# models were run in HPC, submitted in array jobs

rm(list = ls())

# Set user dependent working directories
path2wd <- "/work/wubing/homogenization_occupancy"
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ggplot2", "lme4", "brms", "tidybayes")

for(x in needed_libs){
    if(!require(x,character.only=TRUE, lib.loc = "/gpfs0/home/wubing/R/library")){
		install.packages(x, repos = "http://cran.us.r-project.org", lib="/gpfs0/home/wubing/R/library", dependencies = TRUE)
		require(x, character.only = TRUE, lib.loc = "/gpfs0/home/wubing/R/library")
  }
}

# task index
index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# load occupancy change and meta data
load("intermediate_results/occuapncy.RDATA")
load("intermediate_results/occuapncy_same_locations.RDATA")
load("data/Assemblages_protection.RDATA")
load("data/Assemblages_regions.RDATA")
load("data/Assemblages_taxa.RDATA")


# changes in species richness between early and late periods for each study
sprich_dymamic <- occupancy_change_period %>% 
  group_by(study, database, studyID, sprich) %>%
  summarise(sprich_extin = sum(dynamic == "extinction"),
            sprich_pers = sum(dynamic == "persistent"),
            sprich_colon = sum(dynamic == "colonization"),
            sprich_first = sprich_extin + sprich_pers,
            sprich_last = sprich_colon + sprich_pers,
            sprich_total = sprich_extin + sprich_pers + sprich_colon, 
            sprich_change = log(sprich_last/sprich_first)) %>%
  ungroup() %>%
  mutate(p_sprich_extin = sprich_extin/sprich_total,
         p_sprich_pers = sprich_pers/sprich_total,
         p_sprich_colon = sprich_colon/sprich_total)


# combine meta data with percentage of samples in protected areas, regions, updated taxa,  
# and richness and its change and percentage of species in different dynamics 
dat_meta <- dat_meta %>% 
  left_join(dat_psamp_inPA %>% dplyr::select(study, psamp_inPA_bfstart, psamp_inPA_bflate, psamp_inPA_wdpa, coordinate_local)) %>%
  left_join(region %>% dplyr::select(study, region)) %>%
  left_join(sprich_dymamic %>% dplyr::select(study, sprich_first:p_sprich_colon)) %>%
  dplyr::select(!taxon_new) %>%
  left_join(dat_meta_taxa %>% dplyr::select(study, taxon_new, taxon_final)) %>%
  relocate(taxon_new, taxon_final, .after = taxon)

# for datasets including only samples in same locations across years
dat_sloc_meta <- dat_sloc_meta %>% 
  left_join(dat_psamp_inPA %>% dplyr::select(study, psamp_inPA_bfstart, psamp_inPA_bflate, psamp_inPA_wdpa, coordinate_local)) %>%
  left_join(region %>% dplyr::select(study, region)) %>%
  dplyr::select(!taxon_new) %>%
  left_join(dat_meta_taxa %>% dplyr::select(study, taxon_new, taxon_final)) %>%
  relocate(taxon_new, taxon_final, .after = taxon)

# transform variables and add metadata of studies
oc_period <- occupancy_change_period %>% 
  mutate(lost = ifelse(dynamic == "extinction", 1, 0),
         occup_first_no01 = ifelse(occup_first == 0, 0.01, occup_first),
         occup_first_no01 = ifelse(occup_first  == 1, 0.99, occup_first_no01),
         occup_last_no01 = ifelse(occup_last == 0, 0.01, occup_last),
         occup_last_no01 = ifelse(occup_last  == 1, 0.99, occup_last_no01),         
         occup_first_logit = qlogis(occup_first_no01),
         occup_last_logit = qlogis(occup_last_no01),
         occup_change_logit = occup_last_logit- occup_first_logit,
         cl.aoo10 = log10(aoo10)- mean(log10(aoo10)),
         cl.aoo50 = log10(aoo50)- mean(log10(aoo50)),
         cl.aoo100 = log10(aoo100)- mean(log10(aoo100)),
         cl.aoo10_early = log10(aoo10_early)- mean(log10(aoo10_early), na.rm=TRUE),
         cl.ahull6 = log10(ahull6)- mean(log10(ahull6))) %>%
  left_join(dat_meta %>% dplyr::select(study, taxon_new, taxon_final, realm, region))

oc_2yr <- occupancy_change_2yr %>% 
  mutate(lost = ifelse(dynamic == "extinction", 1, 0),
         occup_first_no01 = ifelse(occup_first == 0, 0.01, occup_first),
         occup_first_no01 = ifelse(occup_first == 1, 0.99, occup_first_no01),
         occup_last_no01 = ifelse(occup_last == 0, 0.01, occup_last),
         occup_last_no01 = ifelse(occup_last == 1, 0.99, occup_last_no01),         
         occup_first_logit = qlogis(occup_first_no01),
         occup_last_logit = qlogis(occup_last_no01),
         occup_change_logit = occup_last_logit- occup_first_logit,
         cl.aoo10 = log10(aoo10)- mean(log10(aoo10)),
         cl.aoo50 = log10(aoo50)- mean(log10(aoo50)),
         cl.aoo100 = log10(aoo100)- mean(log10(aoo100)),
         cl.ahull6 = log10(ahull6)- mean(log10(ahull6))) %>%
  left_join(dat_meta %>% dplyr::select(study, taxon_new, taxon_final, realm, region))

oc_period_sloc <- occupancy_change_sloc_period %>% 
  mutate(lost = ifelse(dynamic == "extinction", 1, 0),
         occup_first_no01 = ifelse(occup_first == 0, 0.01, occup_first),
         occup_first_no01 = ifelse(occup_first  == 1, 0.99, occup_first_no01),
         occup_last_no01 = ifelse(occup_last == 0, 0.01, occup_last),
         occup_last_no01 = ifelse(occup_last  == 1, 0.99, occup_last_no01),         
         occup_first_logit = qlogis(occup_first_no01),
         occup_last_logit = qlogis(occup_last_no01),
         occup_change_logit = occup_last_logit- occup_first_logit,
         cl.aoo10 = log10(aoo10)- mean(log10(aoo10)),
         cl.aoo50 = log10(aoo50)- mean(log10(aoo50)),
         cl.aoo100 = log10(aoo100)- mean(log10(aoo100)),
         cl.ahull6 = log10(ahull6)- mean(log10(ahull6))) %>%
  left_join(dat_meta %>% dplyr::select(study, taxon_new, taxon_final, realm, region))


# category variables are set as factors, which will be included in interaction term with rang size 
oc_period <- oc_period %>%
  mutate(taxon_new = factor(taxon_new, levels = c("Amphibians and reptiles", "Birds", "Fish", 
                                                  "Invertebrates", "Mammals", "Plants")),
         taxon_final = factor(taxon_final, levels = c("Amphibians and reptiles", "Birds", "Mammals",
                                                      "Terrestrial invertebrates", "Terrestrial plants", 
                                                      "Freshwater fish", "Freshwater invertebrates", "Freshwater plants", 
                                                      "Marine fish", "Marine invertebrates", "Marine plants")),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         region = factor(region, levels = c("Africa", "Asia", "Australia", "Europe", "North America", "South America",
                                              "Atlantic", "Pacific", "Arctic Ocean", "Indian Ocean", "Southern Ocean"),
                          labels = c("Africa", "Asia", "Australia", "Europe", "North America", "South America",
                                     "Atlantic Ocean", "Pacific Ocean", "Other Oceans", "Other Oceans", "Other Oceans"))) %>% 
  mutate(realm_taxa = factor(paste(realm, taxon_new, sep = "_")),
         realm_region = factor(paste(realm, region, sep = "_")))

# save(oc_period, oc_2yr, oc_period_sloc, dat_meta, dat_sloc_meta, file = "models/data_input_to_models.RDATA")


###############################
# fit models

###########
# relationship between range size (AOO 10 km) and occupancy change using all species
if(index == 1){
  t1 <- Sys.time()
  brm_oc_aoo10 <- brm(nocc_last|trials(nsamp_used) ~ cl.aoo10 + offset(occup_first_logit) + 
                        (1 + cl.aoo10|study), 
                      family = binomial, data = oc_period ,
                      chains = 4, cores = 4, iter = 8000, thin = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 12), 
                      file = "models/brm_output/brm_oc_aoo10")
  print(brm_oc_aoo10)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10, file = "models/brm_oc_aoo10.RDATA")
}

# sensitivity analyses using AOO 50 km
if(index == 2){
  t1 <- Sys.time()
  brm_oc_aoo50 <- brm(nocc_last|trials(nsamp_used) ~ cl.aoo50 + offset(occup_first_logit) + 
                        (1 + cl.aoo50|study), 
                      family = binomial, data = oc_period ,
                      chains = 4, cores = 4, iter = 8000, thin = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 12), 
                      file = "models/brm_output/brm_oc_aoo50")
  print(brm_oc_aoo50)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo50, file = "models/brm_oc_aoo50.RDATA")
}

# sensitivity analyses using AOO 100 km
if(index == 3){
  t1 <- Sys.time()
  brm_oc_aoo100 <- brm(nocc_last|trials(nsamp_used) ~ cl.aoo100 + offset(occup_first_logit) + 
                         (1 + cl.aoo100|study), 
                       family = binomial, data = oc_period ,
                       chains = 4, cores = 4, iter = 8000, thin = 4,
                       control = list(adapt_delta = 0.9, max_treedepth = 12), 
                       file = "models/brm_output/brm_oc_aoo100")
  print(brm_oc_aoo100)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo100, file = "models/brm_oc_aoo100.RDATA")
}

# sensitivity analyses using EOO (alpha hulls)
if(index == 4){
  t1 <- Sys.time()
  brm_oc_ahull6 <- brm(nocc_last|trials(nsamp_used) ~ cl.ahull6 + offset(occup_first_logit) + 
                         (1 + cl.ahull6|study), 
                       family = binomial, data = oc_period ,
                       chains = 4, cores = 4, iter = 8000, thin = 4,
                       control = list(adapt_delta = 0.9, max_treedepth = 12), 
                       file = "models/brm_output/brm_oc_ahull6")
  print(brm_oc_ahull6)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_ahull6, file = "models/brm_oc_ahull6.RDATA")
}

# sensitivity analyses using occupancy based on assemblages in the first and last years
if(index == 5){
  t1 <- Sys.time()
  brm_oc_aoo10_2yr <- brm(nocc_last|trials(nsamp_used) ~ cl.aoo10 + offset(occup_first_logit) + 
                            (1 + cl.aoo10|study), 
                          family = binomial, data = oc_2yr ,
                          chains = 4, cores = 4, iter = 8000, thin = 4,
                          control = list(adapt_delta = 0.9, max_treedepth = 12), 
                          file = "models/brm_output/brm_oc_aoo10_2yr")
  print(brm_oc_aoo10_2yr)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_2yr, file = "models/brm_oc_aoo10_2yr.RDATA")
}

# sensitivity analyses using only the species that had at least five times more occurrences in GBIF than in the assemblage dataset
if(index == 6){
  t1 <- Sys.time()
  brm_oc_aoo10_rgbif <- brm(nocc_last|trials(nsamp_used) ~ cl.aoo10 + offset(occup_first_logit) + 
                              (1 + cl.aoo10|study), 
                            family = binomial, 
                            data = oc_period %>% filter(ratio_nocc_gbif > 5),
                            chains = 4, cores = 4, iter = 8000, thin = 4,
                            control = list(adapt_delta = 0.9, max_treedepth = 12), 
                            file = "models/brm_output/brm_oc_aoo10_rgbif")
  print(brm_oc_aoo10_rgbif)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_rgbif, file = "models/brm_oc_aoo10_rgbif.RDATA")
}

# sensitivity analyses using species with initial occupancy > 0 and < 1
if(index == 7){
  t1 <- Sys.time()
  brm_oc_aoo10_no01 <- brm(nocc_last|trials(nsamp_used) ~ cl.aoo10 + offset(occup_first_logit) + 
                             (1 + cl.aoo10|study), 
                           family = binomial, 
                           data = oc_period %>% filter(occup_first > 0 & occup_first < 1),
                           chains = 4, cores = 4, iter = 8000, thin = 4,
                           control = list(adapt_delta = 0.9, max_treedepth = 12), 
                           file = "models/brm_output/brm_oc_aoo10_no01")
  print(brm_oc_aoo10_no01)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_no01, file = "models/brm_oc_aoo10_no01.RDATA")
}

# sensitivity analyses using assemblage data with sites in the same locations across years
if(index == 8){
  t1 <- Sys.time()
  brm_oc_aoo10_sloc <- brm(nocc_last|trials(nsamp_used) ~ cl.aoo10 + offset(occup_first_logit) + 
                        (1 + cl.aoo10|study), 
                      family = binomial, data = oc_period_sloc,
                      chains = 4, cores = 4, iter = 8000, thin = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 12), 
                      file = "models/brm_output/brm_oc_aoo10_sloc")
  print(brm_oc_aoo10_sloc)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_sloc, file = "models/brm_oc_aoo10_sloc.RDATA")
}


# include the interaction term between range size and one of study type, realm, taxa and region
# include the interaction term between range size and realm
if(index == 9){
  t1 <- Sys.time()
  brm_oc_aoo10_realm <- brm(nocc_last|trials(nsamp_used) ~ 0 +  realm + cl.aoo10:realm + offset(occup_first_logit) + 
                              (1 + cl.aoo10|study), 
                            family = binomial, data = oc_period ,
                            chains = 4, cores = 4, iter = 8000, thin = 4,
                            control = list(adapt_delta = 0.9, max_treedepth = 12), 
                            file = "models/brm_output/brm_oc_aoo10_realm")
  print(brm_oc_aoo10_realm)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm, file = "models/brm_oc_aoo10_realm.RDATA")
}

# include the interaction between range size and taxa
if(index == 10){
  t1 <- Sys.time()
  brm_oc_aoo10_taxa <- brm(nocc_last|trials(nsamp_used) ~ 0 +  taxon_new + cl.aoo10:taxon_new + offset(occup_first_logit) + 
                             (1 + cl.aoo10|study), 
                           family = binomial, data = oc_period ,
                           chains = 4, cores = 4, iter = 8000, thin = 4,
                           control = list(adapt_delta = 0.9, max_treedepth = 12), 
                           file = "models/brm_output/brm_oc_aoo10_taxa")
  print(brm_oc_aoo10_taxa)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_taxa, file = "models/brm_oc_aoo10_taxa.RDATA")
}

# include the interaction between range size and taxa_final (distinguish plants, fish in different realms)
if(index == 11){
  t1 <- Sys.time()
  brm_oc_aoo10_taxafinal <- brm(nocc_last|trials(nsamp_used) ~ 0 +  taxon_final + cl.aoo10:taxon_final + offset(occup_first_logit) + 
                                  (1 + cl.aoo10|study), 
                                family = binomial, data = oc_period ,
                                chains = 4, cores = 4, iter = 8000, thin = 4,
                                control = list(adapt_delta = 0.9, max_treedepth = 12), 
                                file = "models/brm_output/brm_oc_aoo10_taxafinal")
  print(brm_oc_aoo10_taxafinal)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_taxafinal, file = "models/brm_oc_aoo10_taxafinal.RDATA")
}

# include the interaction between range size and regions (continents and oceans)
if(index == 12){
  t1 <- Sys.time()
  brm_oc_aoo10_region <- brm(nocc_last|trials(nsamp_used) ~ 0 +  region + cl.aoo10:region + offset(occup_first_logit) + 
                                (1 + cl.aoo10|study), 
                              family = binomial, data = oc_period ,
                             chains = 4, cores = 4, iter = 8000, thin = 4,
                             control = list(adapt_delta = 0.9, max_treedepth = 12), 
                              file = "models/brm_output/brm_oc_aoo10_region")
  print(brm_oc_aoo10_region)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_region, file = "models/brm_oc_aoo10_region.RDATA")
}

# include the interaction between range size and the combination of realm and taxa
if(index == 13){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_taxa <- brm(nocc_last|trials(nsamp_used) ~ 0 +  realm_taxa + cl.aoo10:realm_taxa + offset(occup_first_logit) + 
                                   (1 + cl.aoo10|study), 
                                 family = binomial, data = oc_period ,
                                 chains = 4, cores = 4, iter = 8000, thin = 4,
                                 control = list(adapt_delta = 0.9, max_treedepth = 12), 
                                 file = "models/brm_output/brm_oc_aoo10_realm_taxa")
  print(brm_oc_aoo10_realm_taxa)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_taxa, file = "models/brm_oc_aoo10_realm_taxa.RDATA")
}

# include the interaction between range size and the combination of realm and regions
if(index == 14){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_region <- brm(nocc_last|trials(nsamp_used) ~ 0 +  realm_region + cl.aoo10:realm_region + offset(occup_first_logit) + 
                                     (1 + cl.aoo10|study), 
                                   family = binomial, data = oc_period ,
                                   chains = 4, cores = 4, iter = 8000, thin = 4,
                                   control = list(adapt_delta = 0.9, max_treedepth = 12), 
                                   file = "models/brm_output/brm_oc_aoo10_realm_region")
  print(brm_oc_aoo10_realm_region)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_region, file = "models/brm_oc_aoo10_realm_region.RDATA")
}
