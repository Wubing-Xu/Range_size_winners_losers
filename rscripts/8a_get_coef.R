## get the fixed and random effects of range size for the overall pattern and different realms, 
## and predicted/fitted occupancy change, which will be used to draw main figures


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
load("models/brm_oc_aoo10.RDATA")
load("models/brm_oc_aoo10_realm.RDATA")



#############
## global fixed effects of range size
brm_oc_aoo10_fixed <- fixef(brm_oc_aoo10,  probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as_tibble() %>% 
  mutate(term = c("intercept", "sigma_intercept" ,"slope", "sigma_nsamp")) 

# add ranges of range sizes for plotting
brm_oc_aoo10_line <- brm_oc_aoo10_fixed %>% 
  dplyr::select(term, Estimate) %>%
  pivot_wider(names_from = "term", values_from = Estimate) %>%
  mutate(oc_period %>%
           summarise(xmin = min(aoo10),
                     xmax = max(aoo10),
                     cl.xmin = min(cl.aoo10),
                     cl.xmax = max(cl.aoo10)))

# save model summary
summary(brm_oc_aoo10)[["ngrps"]]
summary(brm_oc_aoo10)[["nobs"]]
brmsummary_oc_aoo10 <- summary(brm_oc_aoo10)[["fixed"]] %>%
  as.data.frame() %>%
  mutate(across(Estimate:'Rhat', round, 3))

write.csv(brmsummary_oc_aoo10, "Results/model_summary_brm_oc_aoo10.csv")



#############
## study-level effects of rang size 
brm_oc_aoo10_coef <- coef(brm_oc_aoo10, probs = c(0.025, 0.975))[[1]]
brm_oc_aoo10_coef <- as_tibble(brm_oc_aoo10_coef) %>%
  dplyr::select(1:8) %>%
   mutate(study = rownames(brm_oc_aoo10_coef))
colnames(brm_oc_aoo10_coef)[1:8] <- c("estimate_intercept", "se_intercept", "Q2.5_intercept", "Q97.5_intercept", 
                                      "estimate_slope", "se_slope", "Q2.5_slope", "Q97.5_slope")

# add study-level meta data
brm_oc_aoo10_coef <- brm_oc_aoo10_coef %>% 
  # indicate significance of slopes
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
## get fitted values of occupancy change for each observed range size

brm_oc_aoo10_fitted <- fitted(brm_oc_aoo10, re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(oc_period  %>% 
              dplyr::select(study, aoo10, occup_change_sqroot, nsamp_used))

brm_oc_aoo10_fitted <- brm_oc_aoo10_fitted %>% 
  rename(oc_sqroot_pred = Estimate, oc_sqroot_Q2.5 = Q2.5, oc_sqroot_Q97.5 = Q97.5) %>%
  left_join(dat_meta %>% distinct(study, database, studyID, taxon_new, taxon_final, realm, region)) %>% 
  distinct(aoo10, .keep_all = TRUE)



#############
## get global fixed effects of range size for each realm
brm_oc_aoo10_realm_fixed <- fixef(brm_oc_aoo10_realm, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  filter(! grepl("sigma", realm)) %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":cl.aoo10", "", realm)) %>%
  mutate(term =  rep(c("intercept", "slope"), each = n()/2))
  

# add ranges of range sizes for plotting
brm_oc_aoo10_realm_line <- brm_oc_aoo10_realm_fixed %>% 
  dplyr::select(term, realm, Estimate) %>%
  pivot_wider(names_from = "term", values_from = Estimate) %>% 
  left_join(oc_period %>% group_by(realm) %>%
              summarise(xmin = min(aoo10),
                        xmax = max(aoo10),
                        cl.xmin = min(cl.aoo10),
                        cl.xmax = max(cl.aoo10)))

# save model summary
summary(brm_oc_aoo10_realm)[["ngrps"]]
summary(brm_oc_aoo10_realm)[["nobs"]]
oc_period %>% distinct(realm, study) %>% group_by(realm) %>% summarise(n())
oc_period %>% group_by(realm) %>% summarise(n())

brmsummary_oc_aoo10_realm <- summary(brm_oc_aoo10_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', round, 3))

write.csv(brmsummary_oc_aoo10_realm, "Results/model_summary_brm_oc_aoo10_realm.csv")




#############
## get study-level effects of rang size from the model including the interaction between range size and realm

#  posterior of slopes for three realms
brm_oc_aoo10_realm_global_post <- posterior_samples(brm_oc_aoo10_realm, fixed = TRUE, subset = seq(1, 4000, by = 4),
                       pars = c('b_realmTerrestrial:cl.aoo10', "b_realmFreshwater:cl.aoo10", "b_realmMarine:cl.aoo10")) %>%
  as_tibble() %>%
  rename(slope_terrestrial = "b_realmTerrestrial:cl.aoo10", 
         slope_freshwater = "b_realmFreshwater:cl.aoo10",
         slope_marine = "b_realmMarine:cl.aoo10")

#  posterior of study_level random effects
study_levels <- brm_oc_aoo10_realm$data %>% 
  as_tibble() %>% 
  distinct(study) %>%
  mutate(level = study) %>%
  nest(level) 

brm_oc_aoo10_realm_study_post <- study_levels %>%
  mutate(r_slope = map(data, ~posterior_samples(brm_oc_aoo10_realm, 
                                              pars = paste('r_study[', as.character(.x$level), ',cl.aoo10]', sep=''),
                                              fixed = TRUE,
                                              subset = seq(1, 4000, by = 4)) %>% 
                       unlist() %>% as.numeric())) %>%
  select(-data) %>% 
  unnest(r_slope) 

# add the global and random effects
brm_oc_aoo10_realm_study_post <- brm_oc_aoo10_realm_study_post %>%
  mutate(slope_terrestrial = rep(brm_oc_aoo10_realm_global_post$slope_terrestrial, times = n_distinct(study)),
         slope_freshwater = rep(brm_oc_aoo10_realm_global_post$slope_freshwater, times = n_distinct(study)),
         slope_marine = rep(brm_oc_aoo10_realm_global_post$slope_marine, times = n_distinct(study))) %>% 
  left_join(dat_meta %>% dplyr::select(study, realm)) %>%
  mutate(slope = ifelse(realm == "Terrestrial", r_slope + slope_terrestrial, 
                        ifelse(realm == "Freshwater", r_slope + slope_freshwater,
                               r_slope + slope_marine)))

# calculate study-level effects of rang size
brm_oc_aoo10_realm_coef <- brm_oc_aoo10_realm_study_post %>%
  group_by(study, realm) %>%
  summarise(estimate_slope = mean(slope),
         Q2.5_slope = quantile(slope, prob = 0.05),
         Q97.5_slope = quantile(slope, prob = 0.95)) %>%
  ungroup() %>%
  mutate(sig_slope = ifelse(Q2.5_slope < 0 & Q97.5_slope >0, "neutral", ifelse(Q97.5_slope <= 0, "negative", "positive")),
         sig_slope = factor(sig_slope, levels = c("negative", "positive", "neutral"))) 
  


#############
## get fitted values of  occupancy change for each observed range size from models including 
# interaction between range size and realm

brm_oc_aoo10_realm_fitted <- fitted(brm_oc_aoo10_realm, re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(oc_period  %>% 
              dplyr::select(study, aoo10, occup_change_sqroot, nsamp_used))

brm_oc_aoo10_realm_fitted <- brm_oc_aoo10_realm_fitted %>% 
  rename(oc_sqroot_pred = Estimate, oc_sqroot_Q2.5 = Q2.5, oc_sqroot_Q97.5 = Q97.5) %>%
  left_join(dat_meta %>% distinct(study, database, studyID, taxon_new, taxon_final, realm, region)) %>% 
  distinct(aoo10, realm, .keep_all = TRUE)


save(brm_oc_aoo10_fixed, brm_oc_aoo10_line, brm_oc_aoo10_coef, brm_oc_aoo10_fitted, 
     brm_oc_aoo10_realm_fixed, brm_oc_aoo10_realm_line, brm_oc_aoo10_realm_coef, brm_oc_aoo10_realm_fitted,
     file = "results/coefs_main_brms.RDATA")

