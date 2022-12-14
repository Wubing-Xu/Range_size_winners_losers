## get the effect of protection levels on occupancy change - range size relationships
#  get the predicted occupancy change across range size for different protection levels 

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
load("models/brm_oc_aoo10_realm_psampbflate.RDATA")
load("models/brm_oc_aoo10_realm_psampwdpa.RDATA")
load("models/brm_oc_aoo10_realm_pareabflate.RDATA")
load("models/brm_oc_aoo10_realm_pareawdpa.RDATA")
load("models/brm_oc_aoo10_realm_covariate1.RDATA")

#############
## get effect sizes of protection levels on occupancy change and interaction with range size

# use the proportion of sites in early-established PA
fixef_aoo10_psampbflate <- fixef(brm_oc_aoo10_realm_psampbflate, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  filter(grepl("psamp_inPA_bflate", realm)) %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":psamp_inPA_bflate", "", realm),
         realm = gsub(":cl.aoo10", "", realm)) %>%
  mutate(term =  rep(c("intercept", "slope"), each = n()/2))

# use the proportion of sites in all established PA
fixef_aoo10_psampwdpa <- fixef(brm_oc_aoo10_realm_psampwdpa, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  filter(grepl("psamp_inPA_wdpa", realm)) %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":psamp_inPA_wdpa", "", realm),
         realm = gsub(":cl.aoo10", "", realm)) %>%
  mutate(term =  rep(c("intercept", "slope"), each = n()/2))

# use the proportion of regional area in early-established PA
fixef_aoo10_pareabflate <- fixef(brm_oc_aoo10_realm_pareabflate, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  filter(grepl("parea_inPA_bflate", realm)) %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":parea_inPA_bflate", "", realm),
         realm = gsub(":cl.aoo10", "", realm)) %>%
  mutate(term =  rep(c("intercept", "slope"), each = n()/2))

# use the proportion of regional area in all established PA
fixef_aoo10_pareawdpa <- fixef(brm_oc_aoo10_realm_pareawdpa, probs = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  filter(grepl("parea_inPA_wdpa", realm)) %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":parea_inPA_wdpa", "", realm),
         realm = gsub(":cl.aoo10", "", realm)) %>%
  mutate(term =  rep(c("intercept", "slope"), each = n()/2))



###################
## get the predicted occupancy change across range size at different protection levels for three realm

# prepare the new data
nsamp_used <- oc_period %>% 
  distinct(study,realm, nsamp_used) %>%
  group_by(realm) %>%
  summarise(nsamp_used = round(10^mean(log10(nsamp_used))))

nd <- oc_period %>% 
  filter(!is.na(psamp_inPA_bflate)) %>%
  distinct(realm, aoo10, cl.aoo10) %>%
  left_join(nsamp_used) %>%
  slice(rep(1:n(), each=2)) %>%
  mutate(psamp_inPA_bflate = rep(c(0, 1), time =n()/2))

# predicted occupancy change from the model using proportion of sites in early-established PA
oc_aoo10_psampbflate_fitted <- fitted(brm_oc_aoo10_realm_psampbflate, newdata = nd,
                                      re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(nd)

# predicted occupancy change from the model using proportion of sites in all PA
oc_aoo10_psampwdpa_fitted <- fitted(brm_oc_aoo10_realm_psampwdpa, 
                                      newdata = nd %>% rename(psamp_inPA_wdpa = psamp_inPA_bflate),
                                      re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(nd)

# predicted occupancy change from the model using proportion of regional area in early-established PA
oc_aoo10_pareabflate_fitted <- fitted(brm_oc_aoo10_realm_pareabflate, 
                                      newdata = nd %>% rename(parea_inPA_bflate = psamp_inPA_bflate),
                                      re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(nd)

# predicted occupancy change from the model using proportion of regional area in all PA
oc_aoo10_pareawdpa_fitted <- fitted(brm_oc_aoo10_realm_pareawdpa, 
                                      newdata = nd %>% rename(parea_inPA_wdpa = psamp_inPA_bflate),
                                      re_formula = NA, nsamples = 1000) %>%
  as_tibble() %>% 
  bind_cols(nd)


save(fixef_aoo10_psampbflate, fixef_aoo10_psampwdpa, fixef_aoo10_pareabflate, fixef_aoo10_pareawdpa,
     oc_aoo10_psampbflate_fitted, oc_aoo10_psampwdpa_fitted, oc_aoo10_pareabflate_fitted, oc_aoo10_pareawdpa_fitted,
     file = "results/fixef_prediction_protection.RDATA")


####################
# save summary of the model used for main figure (fig. 3)  
summary(brm_oc_aoo10_realm_psampbflate)[["ngrps"]]
summary(brm_oc_aoo10_realm_psampbflate)[["nobs"]]
oc_period %>% filter(!is.na(psamp_inPA_bflate)) %>% distinct(realm, study) %>% group_by(realm) %>% summarise(n())
oc_period %>% filter(!is.na(psamp_inPA_bflate)) %>%  group_by(realm) %>% summarise(n())

brmsummary_oc_aoo10_realm_psampbflate <- summary(brm_oc_aoo10_realm_psampbflate)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', round, 3))

write.csv(brmsummary_oc_aoo10_realm_psampbflate, "Results/model_summary_brm_oc_aoo10_realm_protection.csv")



####################
# get the fixed effects from the model considering range size, realm, protection status and other six covariates, and save as a table 

fixef_aoo10_allCovariates <- fixef(brm_oc_aoo10_realm_covariate1, probs = c(0.025, 0.975, 0.05,0.95)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:Q95, round, 4))

write.csv(fixef_aoo10_allCovariates, file = "results/Table.S_fixef_aoo10_allCovariates.csv")

brmsummary_oc_aoo10_realm_allCovariates <- summary(brm_oc_aoo10_realm_covariate1)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'u-95% CI', round, 4),
         Rhat = round(Rhat, 3))

write.csv(brmsummary_oc_aoo10_realm_allCovariates, file = "results/model_summary_oc_aoo10_realm_allCovariates.csv")

