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
load("intermediate_results/occupancy.RDATA")
load("intermediate_results/occupancy_same_locations.RDATA")
load("data/Assemblages_protection.RDATA")
load("data/Assemblages_regions.RDATA")
load("data/Assemblages_taxa.RDATA")


# combine meta data with percentage of samples in protected areas, regions, updated taxa,  
# and richness and its change and percentage of species in different dynamics 
dat_meta <- dat_meta %>% 
  mutate(nsamp_used = min_samp* nyears_combined) %>%
  left_join(dat_coverage_PA %>% dplyr::select(study, coordinate_local, psamp_inPA_bfstart, psamp_inPA_bflate, psamp_inPA_wdpa, 
                                              parea_inPA_bflate, parea_inPA_wdpa)) %>%
  left_join(region %>% dplyr::select(study, region)) %>%
  dplyr::select(!taxon_new) %>%
  left_join(dat_meta_taxa %>% dplyr::select(study, taxon_new, taxon_final)) %>%
  relocate(taxon_new, taxon_final, .after = taxon) %>%
  mutate(lat_cent = abs(cent_lat) - mean(abs(cent_lat)),
         sprich_cent = log10(sprich) - mean(log10(sprich)),       
         duration_cent = log10(duration) - mean(log10(duration)),
         startyr_cent = start_year - mean(start_year),
         extent_cent = log10(extent_km2) - mean(log10(extent_km2)),
         nsamp_cent = log10(nsamp_used) - mean(log10(nsamp_used)))
  
# for datasets including only samples in same locations across years
dat_sloc_meta <- dat_sloc_meta %>% 
  left_join(region %>% dplyr::select(study, region)) %>%
  dplyr::select(!taxon_new) %>%
  left_join(dat_meta_taxa %>% dplyr::select(study, taxon_new, taxon_final)) %>%
  relocate(taxon_new, taxon_final, .after = taxon)


# transform variables and add metadata of studies
oc_period <- occupancy_change_period %>% 
  mutate(occup_first_asin = asin(sqrt(occup_first)),
         occup_last_asin = asin(sqrt(occup_last)),
         occup_change_asin = occup_last_asin - occup_first_asin,
         occup_change_sqroot = ifelse(occup_change >= 0, sqrt(abs(occup_change)), -sqrt(abs(occup_change))),
         cl.aoo10 = log10(aoo10)- mean(log10(aoo10)),
         cl.aoo50 = log10(aoo50)- mean(log10(aoo50)),
         cl.aoo100 = log10(aoo100)- mean(log10(aoo100)),
         cl.ahull6 = log10(ahull6)- mean(log10(ahull6))) %>%
  left_join(dat_meta %>% dplyr::select(study, taxon_new, taxon_final, realm, region, 
                                       coordinate_local, psamp_inPA_bfstart, psamp_inPA_bflate, psamp_inPA_wdpa, parea_inPA_bflate, parea_inPA_wdpa,
                                       lat_cent, sprich_cent, duration_cent, startyr_cent, extent_cent, nsamp_cent))


oc_period_sloc <- occupancy_change_sloc_period %>% 
  mutate(occup_first_asin = asin(sqrt(occup_first)),
         occup_last_asin = asin(sqrt(occup_last)),
         occup_change_asin = occup_last_asin - occup_first_asin,
         occup_change_sqroot = ifelse(occup_change >= 0, sqrt(abs(occup_change)), -sqrt(abs(occup_change))),
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

# save(oc_period, oc_period_sloc, dat_meta, dat_sloc_meta, file = "models/data_input_to_models.RDATA")


###############################
# fit models

###########
# relationship between range size (AOO 10 km) and occupancy change using all species
if(index == 1){
  t1 <- Sys.time()
  brm_oc_aoo10 <- brm(bf(occup_change_sqroot ~ cl.aoo10 + (1 + cl.aoo10|study),
                         sigma ~ log10(nsamp_used)), 
                      family = gaussian(), data = oc_period,
                      chains = 4, cores = 4, iter = 8000, thin = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 10), 
                      file = "models/brm_output/brm_oc_aoo10")
  print(brm_oc_aoo10)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10, file = "models/brm_oc_aoo10.RDATA")
}

# sensitivity analyses using AOO 50 km
if(index == 2){
  t1 <- Sys.time()
  brm_oc_aoo50 <- brm(bf(occup_change_sqroot ~ cl.aoo50 + (1 + cl.aoo50|study),
                         sigma ~ log10(nsamp_used)), 
                      family = gaussian, data = oc_period,
                      chains = 4, cores = 4, iter = 8000, thin = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 10), 
                      file = "models/brm_output/brm_oc_aoo50")
  print(brm_oc_aoo50)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo50, file = "models/brm_oc_aoo50.RDATA")
}

# sensitivity analyses using AOO 100 km
if(index == 3){
  t1 <- Sys.time()
  brm_oc_aoo100 <- brm(bf(occup_change_sqroot ~ cl.aoo100 + (1 + cl.aoo100|study),
                          sigma ~ log10(nsamp_used)), 
                       family = gaussian, data = oc_period,
                       chains = 4, cores = 4, iter = 8000, thin = 4,
                       control = list(adapt_delta = 0.9, max_treedepth = 10), 
                       file = "models/brm_output/brm_oc_aoo100")
  print(brm_oc_aoo100)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo100, file = "models/brm_oc_aoo100.RDATA")
}

# sensitivity analyses using EOO (alpha hulls)
if(index == 4){
  t1 <- Sys.time()
  brm_oc_ahull6 <- brm(bf(occup_change_sqroot ~ cl.ahull6 + (1 + cl.ahull6|study),
                          sigma ~ log10(nsamp_used)), 
                       family = gaussian, data = oc_period,
                       chains = 4, cores = 4, iter = 8000, thin = 4,
                       control = list(adapt_delta = 0.9, max_treedepth = 10), 
                       file = "models/brm_output/brm_oc_ahull6")
  print(brm_oc_ahull6)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_ahull6, file = "models/brm_oc_ahull6.RDATA")
}


# sensitivity analyses using only the species that had at least five times more occurrences in GBIF than in the assemblage dataset
if(index == 5){
  t1 <- Sys.time()
  brm_oc_aoo10_rgbif <- brm(bf(occup_change_sqroot ~ cl.aoo10 + (1 + cl.aoo10|study),
                               sigma ~ log10(nsamp_used)), 
                            family = gaussian(), data = oc_period %>% filter(ratio_nocc_gbif > 5), 
                            chains = 4, cores = 4, iter = 8000, thin = 4,
                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                            file = "models/brm_output/brm_oc_aoo10_rgbif")
  print(brm_oc_aoo10_rgbif)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_rgbif, file = "models/brm_oc_aoo10_rgbif.RDATA")
}


# sensitivity analyses using assemblage data with sites in the same locations across years
if(index == 6){
  t1 <- Sys.time()
  brm_oc_aoo10_sloc <- brm(bf(occup_change_sqroot ~ cl.aoo10 + (1 + cl.aoo10|study),
                              sigma ~ log10(nsamp_used)), 
                           family = gaussian(), data = oc_period_sloc, 
                      chains = 4, cores = 4, iter = 8000, thin = 4,
                      control = list(adapt_delta = 0.9, max_treedepth = 10), 
                      file = "models/brm_output/brm_oc_aoo10_sloc")
  print(brm_oc_aoo10_sloc)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_sloc, file = "models/brm_oc_aoo10_sloc.RDATA")
}


# include the interaction term between range size and one of study type, realm, taxa and region
# include the interaction term between range size and realm
if(index == 7){
  t1 <- Sys.time()
  brm_oc_aoo10_realm <- brm(bf(occup_change_sqroot ~ 0 +  realm + cl.aoo10:realm + (1 + cl.aoo10|study),
                               sigma ~ log10(nsamp_used)), 
                            family = gaussian(), data = oc_period,
                            chains = 4, cores = 4, iter = 8000, thin = 4,
                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                            file = "models/brm_output/brm_oc_aoo10_realm")
  print(brm_oc_aoo10_realm)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm, file = "models/brm_oc_aoo10_realm.RDATA")
}

# include the interaction between range size and taxa
if(index == 8){
  t1 <- Sys.time()
  brm_oc_aoo10_taxa <- brm(bf(occup_change_sqroot ~ 0 +  taxon_new + cl.aoo10:taxon_new + (1 + cl.aoo10|study),
                              sigma ~ log10(nsamp_used)), 
                           family = gaussian(), data = oc_period, 
                           chains = 4, cores = 4, iter = 8000, thin = 4,
                           control = list(adapt_delta = 0.9, max_treedepth = 10), 
                           file = "models/brm_output/brm_oc_aoo10_taxa")
  print(brm_oc_aoo10_taxa)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_taxa, file = "models/brm_oc_aoo10_taxa.RDATA")
}

# include the interaction between range size and taxa_final (distinguish plants, fish in different realms)
if(index == 9){
  t1 <- Sys.time()
  brm_oc_aoo10_taxafinal <- brm(bf(occup_change_sqroot ~ 0 +  taxon_final + cl.aoo10:taxon_final + (1 + cl.aoo10|study),
                                   sigma ~ log10(nsamp_used)), 
                                family = gaussian(), data = oc_period, 
                                chains = 4, cores = 4, iter = 8000, thin = 4,
                                control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                file = "models/brm_output/brm_oc_aoo10_taxafinal")
  print(brm_oc_aoo10_taxafinal)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_taxafinal, file = "models/brm_oc_aoo10_taxafinal.RDATA")
}

# include the interaction between range size and regions (continents and oceans)
if(index == 10){
  t1 <- Sys.time()
  brm_oc_aoo10_region <- brm(bf(occup_change_sqroot ~ 0 +  region + cl.aoo10:region + (1 + cl.aoo10|study),
                                sigma ~ log10(nsamp_used)), 
                             family = gaussian(), data = oc_period, 
                             chains = 4, cores = 4, iter = 8000, thin = 4,
                             control = list(adapt_delta = 0.9, max_treedepth = 10), 
                              file = "models/brm_output/brm_oc_aoo10_region")
  print(brm_oc_aoo10_region)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_region, file = "models/brm_oc_aoo10_region.RDATA")
}

# include the interaction between range size and the combination of realm and taxa
if(index == 11){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_taxa <- brm(bf(occup_change_sqroot ~ 0 +  realm_taxa + cl.aoo10:realm_taxa + (1 + cl.aoo10|study),
                                    sigma ~ log10(nsamp_used)), 
                                 family = gaussian(), data = oc_period, 
                                 chains = 4, cores = 4, iter = 8000, thin = 4,
                                 control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                 file = "models/brm_output/brm_oc_aoo10_realm_taxa")
  print(brm_oc_aoo10_realm_taxa)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_taxa, file = "models/brm_oc_aoo10_realm_taxa.RDATA")
}

# include the interaction between range size and the combination of realm and regions
if(index == 12){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_region <- brm(bf(occup_change_sqroot ~ 0 +  realm_region + cl.aoo10:realm_region + (1 + cl.aoo10|study),
                                      sigma ~ log10(nsamp_used)), 
                                   family = gaussian(), data = oc_period, 
                                   chains = 4, cores = 4, iter = 8000, thin = 4,
                                   control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                   file = "models/brm_output/brm_oc_aoo10_realm_region")
  print(brm_oc_aoo10_realm_region)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_region, file = "models/brm_oc_aoo10_realm_region.RDATA")
}


# include the interaction between range size and protection status for realms
if(index == 13){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_psampwdpa <- brm(bf(occup_change_sqroot ~ 0 +  realm + (cl.aoo10*psamp_inPA_wdpa):realm + (1 + cl.aoo10|study),
                               sigma ~ log10(nsamp_used)), 
                            family = gaussian(), data = oc_period,
                            chains = 4, cores = 4, iter = 8000, thin = 4,
                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                            file = "models/brm_output/brm_oc_aoo10_realm_psampwdpa")
  print(brm_oc_aoo10_realm_psampwdpa)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_psampwdpa, file = "models/brm_oc_aoo10_realm_psampwdpa.RDATA")
}

if(index == 14){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_psampbflate <- brm(bf(occup_change_sqroot ~ 0 +  realm + (cl.aoo10*psamp_inPA_bflate):realm + (1 + cl.aoo10|study),
                                         sigma ~ log10(nsamp_used)), 
                                      family = gaussian(), data = oc_period,
                                      chains = 4, cores = 4, iter = 8000, thin = 4,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                      file = "models/brm_output/brm_oc_aoo10_realm_psampbflate")
  print(brm_oc_aoo10_realm_psampbflate)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_psampbflate, file = "models/brm_oc_aoo10_realm_psampbflate.RDATA")
}

if(index == 15){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_pareawdpa <- brm(bf(occup_change_sqroot ~ 0 +  realm + (cl.aoo10*parea_inPA_wdpa):realm + (1 + cl.aoo10|study),
                                         sigma ~ log10(nsamp_used)), 
                                      family = gaussian(), data = oc_period,
                                      chains = 4, cores = 4, iter = 8000, thin = 4,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                      file = "models/brm_output/brm_oc_aoo10_realm_pareawdpa")
  print(brm_oc_aoo10_realm_pareawdpa)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_pareawdpa, file = "models/brm_oc_aoo10_realm_pareawdpa.RDATA")
}

if(index == 16){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_pareabflate <- brm(bf(occup_change_sqroot ~ 0 +  realm + (cl.aoo10*parea_inPA_bflate):realm + (1 + cl.aoo10|study),
                                         sigma ~ log10(nsamp_used)), 
                                      family = gaussian(), data = oc_period,
                                      chains = 4, cores = 4, iter = 8000, thin = 4,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                      file = "models/brm_output/brm_oc_aoo10_realm_pareabflate")
  print(brm_oc_aoo10_realm_pareabflate)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_pareabflate, file = "models/brm_oc_aoo10_realm_pareabflate.RDATA")
}



if(index == 17){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_covariate <- brm(bf(occup_change_sqroot ~ 0 +  realm + (cl.aoo10 + psamp_inPA_wdpa + cl.aoo10:psamp_inPA_wdpa +
                                                                             lat_cent + cl.aoo10:lat_cent + 
                                                                             sprich_cent + cl.aoo10:sprich_cent + 
                                                                             duration_cent + cl.aoo10:duration_cent + 
                                                                             startyr_cent + cl.aoo10:startyr_cent + 
                                                                             extent_cent + cl.aoo10:extent_cent + 
                                                                             nsamp_cent + cl.aoo10:nsamp_cent):realm + (1 + cl.aoo10|study),
                                         sigma ~ log10(nsamp_used)), 
                                      family = gaussian(), data = oc_period,
                                      chains = 4, cores = 4, iter = 8000, thin = 4,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                      file = "models/brm_output/brm_oc_aoo10_realm_covariate")
  print(brm_oc_aoo10_realm_covariate)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_covariate, file = "models/brm_oc_aoo10_realm_covariate.RDATA")
}


if(index == 18){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_covariate1 <- brm(bf(occup_change_sqroot ~ 0 +  realm + (cl.aoo10 + psamp_inPA_bflate + cl.aoo10:psamp_inPA_bflate +
                                                                             lat_cent + cl.aoo10:lat_cent + 
                                                                             sprich_cent + cl.aoo10:sprich_cent + 
                                                                             duration_cent + cl.aoo10:duration_cent + 
                                                                             startyr_cent + cl.aoo10:startyr_cent + 
                                                                             extent_cent + cl.aoo10:extent_cent + 
                                                                             nsamp_cent + cl.aoo10:nsamp_cent):realm + (1 + cl.aoo10|study),
                                         sigma ~ log10(nsamp_used)), 
                                      family = gaussian(), data = oc_period,
                                      chains = 4, cores = 4, iter = 8000, thin = 4,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                      file = "models/brm_output/brm_oc_aoo10_realm_covariate1")
  print(brm_oc_aoo10_realm_covariate1)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_covariate1, file = "models/brm_oc_aoo10_realm_covariate1.RDATA")
}


if(index == 19){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_covariate2 <- brm(bf(occup_change_sqroot ~ 0 +  realm + (cl.aoo10 + parea_inPA_wdpa+ cl.aoo10:parea_inPA_wdpa +
                                                                              lat_cent + cl.aoo10:lat_cent + 
                                                                              sprich_cent + cl.aoo10:sprich_cent + 
                                                                              duration_cent + cl.aoo10:duration_cent + 
                                                                              startyr_cent + cl.aoo10:startyr_cent + 
                                                                              extent_cent + cl.aoo10:extent_cent + 
                                                                              nsamp_cent + cl.aoo10:nsamp_cent):realm + (1 + cl.aoo10|study),
                                          sigma ~ log10(nsamp_used)), 
                                       family = gaussian(), data = oc_period,
                                       chains = 4, cores = 4, iter = 8000, thin = 4,
                                       control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                       file = "models/brm_output/brm_oc_aoo10_realm_covariate2")
  print(brm_oc_aoo10_realm_covariate2)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_covariate2, file = "models/brm_oc_aoo10_realm_covariate2.RDATA")
}


if(index == 20){
  t1 <- Sys.time()
  brm_oc_aoo10_realm_covariate3 <- brm(bf(occup_change_sqroot ~ 0 +  realm + (cl.aoo10 + parea_inPA_bflate + cl.aoo10:parea_inPA_bflate +
                                                                              lat_cent + cl.aoo10:lat_cent + 
                                                                              sprich_cent + cl.aoo10:sprich_cent + 
                                                                              duration_cent + cl.aoo10:duration_cent + 
                                                                              startyr_cent + cl.aoo10:startyr_cent + 
                                                                              extent_cent + cl.aoo10:extent_cent + 
                                                                              nsamp_cent + cl.aoo10:nsamp_cent):realm + (1 + cl.aoo10|study),
                                          sigma ~ log10(nsamp_used)), 
                                       family = gaussian(), data = oc_period,
                                       chains = 4, cores = 4, iter = 8000, thin = 4,
                                       control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                       file = "models/brm_output/brm_oc_aoo10_realm_covariate3")
  print(brm_oc_aoo10_realm_covariate3)
  t2 <- Sys.time()
  print(t2 - t1)
  save(brm_oc_aoo10_realm_covariate3, file = "models/brm_oc_aoo10_realm_covariate3.RDATA")
}
