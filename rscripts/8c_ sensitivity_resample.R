## test the sensitivity of our result to the rarefaction process
## repeat the same model structure that was used in main text with 200 resampled datasets using the package 'lme4'
## draw a figure to compare the effects of range size shown in the main text with those from 200 rarefaction processes

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
				  "IDIVTS01" = "H:/wubing/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ggplot2", "lme4")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)

load("intermediate_results/occuapncy_resample.RDATA")
load("models/data_input_to_models.RDATA")
load("results/coefs_main_brms.RDATA")


##########
# fit the same data set using GLMM
glmm_oc_aoo10 <- glmer(occup_last ~ cl.aoo10 + offset(occup_first_logit) + (1 +  cl.aoo10|study), 
                       weights = nsamp_used , family = binomial, 
                       data = oc_period)
summary(glmm_oc_aoo10)

# get the global slopes and study-level slopes
glmm_oc_aoo10_fixed <- fixef(glmm_oc_aoo10)

glmm_oc_aoo10_coef <- coef(glmm_oc_aoo10)[[1]] %>% 
  rownames_to_column(var = "study") %>%
  as_tibble() %>% 
  rename(glmm_intercept = "(Intercept)", glmm_slope = cl.aoo10)


##########
# refit models using 200 resampled datasets

# add variables that will be used in models
occupancy_change_resample <- occupancy_change_resample %>% 
  mutate(occup_first_no01 = ifelse(occup_first == 0, 0.01, occup_first),
         occup_first_no01 = ifelse(occup_first == 1, 0.99, occup_first_no01),
         occup_first_logit = qlogis(occup_first_no01),
         cl.aoo10 = log10(aoo10)- mean(log10(aoo10)))

# refit models and get the fixed effects
oc_aoo10_slopes_resample <- NULL 
for(i in unique(occupancy_change_resample$resample)){
  glmm_oc_aoo10_resample <- try(glmer(occup_last ~ cl.aoo10 + offset(occup_first_logit) + (1 +  cl.aoo10|study), 
                         weights = nsamp_used , family = binomial, 
                         data = occupancy_change_resample %>% filter(resample == i)), silent = TRUE)
  if(!inherits(glmm_oc_aoo10_resample, "try-error")){
    slope_fixef <- fixef(glmm_oc_aoo10_resample) %>% t()
    oc_aoo10_slopes_resample <- rbind(oc_aoo10_slopes_resample, cbind(resample = i, slope_fixef))
  }
}

dim(oc_aoo10_slopes_resample)
oc_aoo10_slopes_resample <- oc_aoo10_slopes_resample %>% 
  as_tibble() %>%
  rename(glmm_intercept = "(Intercept)", glmm_slope = cl.aoo10)

save(oc_aoo10_slopes_resample, file = "results/coefs_resample.Rdata")


####
# generate figures

# combine global and study-level slopes from brms and glm models
oc_aoo10_coef <- brm_oc_aoo10_coef %>%
  dplyr::select(study, brm_slope = estimate_slope, brm_Q2.5 = Q2.5_slope, brm_Q97.5 = Q97.5_slope) %>%
  left_join(glmm_oc_aoo10_coef %>% dplyr::select(study, glmm_slope))

oc_aoo10_fixed <- brm_oc_aoo10_fixed %>%
  dplyr::select(term, brm_slope = Estimate, brm_Q2.5 = Q2.5, brm_Q97.5 = Q97.5) %>%
  mutate(glmm_slope = glmm_oc_aoo10_fixed) %>% 
  filter(term == "slope")

# compare slopes from brms and glm models
plot_slops_models <- ggplot(data = oc_aoo10_coef) +
  geom_point(aes(x = brm_slope, y = glmm_slope), shape = 16, size = 0.8, alpha = 0.3) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  geom_point(data = oc_aoo10_fixed, aes(x =brm_slope, y = glmm_slope), 
             shape = 16, size = 1.5, alpha = 1, color = "red") + 
  geom_abline(intercept = 0, slope = 1, linetype =2, size = 0.3, colour = "blue", alpha = 0.5) + 
  labs(x = "brms model slopes", y = "GLM model slopes", tag = "A") + 
  theme_classic() + 
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'))
 
# frequency distribution of overall estimates of slopes from glmm 
plot_slopes_resample <- ggplot(oc_aoo10_slopes_resample) +
  geom_histogram(aes(glmm_slope)) + 
  geom_vline(xintercept = oc_aoo10_fixed$brm_slope, linetype =2, colour = "black") + 
  labs(x = "Overall estimates of slope", y = "Number of resamples", tag = "B") +
  theme_classic() + 
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'))

cowplot::plot_grid(plot_slops_models,plot_slopes_resample, nrow = 1, align = "hv")


ggsave("results/Fig.S_slopes_resample.png", width = 120, height = 60, units = 'mm', dpi = 600)

