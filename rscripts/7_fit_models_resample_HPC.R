## test the sensitivity of our result (overall range size-occupancy change relationship) to the rarefaction process
## repeat the same model structure that was used in main text with 200 resampled datasets
# models were fitted based on Bayesian inference using the R package brms
# models were run in HPC


rm(list = ls())

# Set user dependent working directories
path2wd <- "/work/wubing/homogenization_occupancy"
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ggplot2", "lme4", "brms", "tidybayes", "parallel")

for(x in needed_libs){
    if(!require(x,character.only=TRUE, lib.loc = "/gpfs0/home/wubing/R/library")){
		install.packages(x, repos = "http://cran.us.r-project.org", lib="/gpfs0/home/wubing/R/library", dependencies = TRUE)
		require(x, character.only = TRUE, lib.loc = "/gpfs0/home/wubing/R/library")
  }
}


# load occupancy change and meta data
load("intermediate_results/occupancy_resample.RDATA")


# add variables that will be used in models
occupancy_change_resample <- occupancy_change_resample %>% 
  mutate(occup_first_asin = asin(sqrt(occup_first)),
         occup_last_asin = asin(sqrt(occup_last)),
         occup_change_asin = occup_last_asin - occup_first_asin,
         occup_change_sqroot = ifelse(occup_change >= 0, sqrt(abs(occup_change)), -sqrt(abs(occup_change))),
         cl.aoo10 = log10(aoo10)- mean(log10(aoo10)))

# set parallel computing
t1 <- Sys.time()
print("n_cores")
print(detectCores()) #the number of cores
no_cores <- 40 
cl <- makeCluster(no_cores, type="PSOCK") # Initiate cluster
clusterExport(cl, varlist = c("occupancy_change_resample"))
clusterEvalQ(cl, c(library(brms), library(dplyr), library(tibble)))


# a self_defined function to perform models and get fixed effect size
get_fixef_aoo10 <- function(i) {
  brm_oc_aoo10_resample <- brm(bf(occup_change_sqroot ~ cl.aoo10 + (1 + cl.aoo10|study),
                                  sigma ~ log10(nsamp_used)), 
                               family = gaussian(), data = occupancy_change_resample %>% filter(resample == i),
                               chains = 4, cores = 4, iter = 4000, thin = 2,
                               control = list(adapt_delta = 0.9, max_treedepth = 10))
  brm_oc_aoo10_fixef <- fixef(brm_oc_aoo10_resample) %>%
    as.data.frame() %>%
    rownames_to_column(var = "term") %>%
    mutate(resample = i)
  return(brm_oc_aoo10_fixef)
  
}

brm_oc_aoo10_fixef_resamples <- parLapply(cl, 1:200, get_fixef_aoo10)

save(brm_oc_aoo10_fixef_resamples,
     file="models/brm_oc_aoo10_fixef_resamples.RDATA")
t2 <- Sys.time()
t2 - t1

stopCluster(cl)
