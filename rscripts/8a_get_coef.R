## get the fixed and random effects of range size for the overall pattern and different realms, 
## and predicted/fitted occupancy and occupancy, which will be used to draw main figures


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
brm_oc_aoo10_fixed <- fixef(brm_oc_aoo10, robust = TRUE, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as_tibble() %>%
  mutate(term = c("intercept", "slope")) 

# add ranges of range sizes for plotting
brm_oc_aoo10_line <- brm_oc_aoo10_fixed %>% 
  dplyr::select(term, Estimate) %>%
  pivot_wider(names_from = "term", values_from = Estimate) %>%
  mutate(oc_period %>%
           summarise(xmin = min(aoo10),
                     xmax = max(aoo10),
                     cl.xmin = min(cl.aoo10),
                     cl.xmax = max(cl.aoo10)))


#############
## random effects of rang size for each study
brm_oc_aoo10_coef <- coef(brm_oc_aoo10, robust = TRUE, probs = c(0.025, 0.975))[[1]]
brm_oc_aoo10_coef <- as_tibble(brm_oc_aoo10_coef) %>%
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
## get fitted values of occupancy and occupancy change for each observed range size
# get the fitted values: the number of sites occupied at the second period
brm_oc_aoo10_fitted <- fitted(brm_oc_aoo10, re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(oc_period  %>% 
              dplyr::select(study, aoo10, occup_change_logit, nsamp_used, occup_first_logit, occup_last_logit ))

# calculate the predicted occupancy change
brm_oc_aoo10_fitted <- brm_oc_aoo10_fitted %>% 
  mutate(occup_last_logit_pred =  qlogis(Estimate/nsamp_used),
         occup_last_logit_Q2.5 =  qlogis(Q2.5/nsamp_used),
         occup_last_logit_Q97.5 =  qlogis(Q97.5/nsamp_used),
         occup_change_logit_pred = occup_last_logit_pred - occup_first_logit,
         occup_change_logit_Q2.5 = occup_last_logit_Q2.5 - occup_first_logit,
         occup_change_logit_Q97.5 = occup_last_logit_Q97.5 - occup_first_logit) %>%
  left_join(dat_meta %>% distinct(study, database, studyID, taxon_new, taxon_final, realm, climate)) %>% 
  distinct(aoo10, .keep_all = TRUE)


#############
## get global fixed effects of range size for each realm
brm_oc_aoo10_realm_fixed <- fixef(brm_oc_aoo10_realm, robust = TRUE, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":cl.aoo10", "", realm))

# add ranges of range sizes for plotting
brm_oc_aoo10_realm_line <- brm_oc_aoo10_realm_fixed %>% 
  dplyr::select(term, realm, Estimate) %>%
  pivot_wider(names_from = "term", values_from = Estimate) %>% 
  left_join(oc_period %>% group_by(realm) %>%
              summarise(xmin = min(aoo10),
                        xmax = max(aoo10),
                        cl.xmin = min(cl.aoo10),
                        cl.xmax = max(cl.aoo10)))


## get fitted values of occupancy and occupancy change for each observed range size from models including 
# interaction between range size and realm
# get the fitted values: the number of sites occupied at the second period
brm_oc_aoo10_realm_fitted <- fitted(brm_oc_aoo10_realm, re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(oc_period  %>% 
              dplyr::select(study, aoo10, occup_change_logit, nsamp_used, occup_first_logit, occup_last_logit))

# calculate the predicted occupancy change
brm_oc_aoo10_realm_fitted <- brm_oc_aoo10_realm_fitted %>% 
  mutate(occup_last_logit_pred =  qlogis(Estimate/nsamp_used),
         occup_last_logit_Q2.5 =  qlogis(Q2.5/nsamp_used),
         occup_last_logit_Q97.5 =  qlogis(Q97.5/nsamp_used),
         occup_change_logit_pred = occup_last_logit_pred - occup_first_logit,
         occup_change_logit_Q2.5 = occup_last_logit_Q2.5 - occup_first_logit,
         occup_change_logit_Q97.5 = occup_last_logit_Q97.5 - occup_first_logit) %>%
  left_join(dat_meta %>% distinct(study, database, studyID, taxon_new, taxon_final, realm, region)) %>% 
  distinct(aoo10, realm, .keep_all = TRUE)


save(brm_oc_aoo10_fixed, brm_oc_aoo10_line, brm_oc_aoo10_coef, brm_oc_aoo10_fitted, 
     brm_oc_aoo10_realm_fixed, brm_oc_aoo10_realm_line, brm_oc_aoo10_realm_fitted,
     file = "results/coefs_main_brms.RDATA")

