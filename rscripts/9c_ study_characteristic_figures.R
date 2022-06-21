#######################################
##  frequency distribution plot of six study characteristics: latitude, regional species richness, 
## extent, number of samples in each period, duration and start year of sampling;
## plot showing correlations among study characteristics
## fit models assessing relationships between study-level slopes (occupancy change against range size) and study characteristics and draw plots

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ggplot2", "ggridges", "brms", "tidybayes", "cowplot", "RColorBrewer")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {  
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


load("models/data_input_to_models.RDATA")
load("results/coefs_supplement_brms.RDATA")
load("results/coefs_main_brms.RDATA")


dat_meta <- dat_meta %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  left_join(oc_period %>% distinct(study, nsamp_used))

brm_oc_aoo10_coef <- brm_oc_aoo10_coef %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  left_join(oc_period %>% distinct(study, nsamp_used))

# an error in the raw data: the year 2015 was incorrectly input as 2915. Correct here
id <- dat_meta$study_name == "price_2017_Illinois" & !is.na(dat_meta$study_name)
dat_meta$end_year[id] <- 2015
dat_meta$duration[id] <- 16
dat_meta$duration_mean[id] <- 16
id <- brm_oc_aoo10_coef$study_name == "price_2017_Illinois" & !is.na(brm_oc_aoo10_coef$study_name)
brm_oc_aoo10_coef$end_year[id] <- 2015
brm_oc_aoo10_coef$duration[id] <- 16
brm_oc_aoo10_coef$duration_mean[id] <- 16

realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#1B9E77", "Marine" = "#7570B3")


############
## Figure showing frequency distribution of study characteristics
# cent_lat, sprich, extent_km2, nsample_used, duration, start_year

plot_meta_latitude <- ggplot(data = dat_meta) +
  geom_histogram(aes(cent_lat, fill = realm)) + 
  #scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) + 
  labs(x = "Latitude (degree)", y = NULL, fill = NULL, tag = "A") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_richness <- ggplot(data = dat_meta) +
  geom_histogram(aes(sprich_total, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Regional species richness", y = NULL, tag = "B") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_extent <- ggplot(data = dat_meta) +
  geom_histogram(aes(extent_km2, fill = realm)) + 
  scale_x_log10(breaks = c(10^0, 10^3, 10^6)) + 
  labs(x = bquote('Extent of study sites ('*km^2*')'), y = NULL, tag = "C") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_nsamp <- ggplot(data = dat_meta) +
  geom_histogram(aes(nsamp_used, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Number of samples", y = NULL, tag = "D") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_duration <- ggplot(data = dat_meta) +
  geom_histogram(aes(duration, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Duration of sampling (years)", y = NULL, tag = "E") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_start <- ggplot(data = dat_meta) +
  geom_histogram(aes(start_year, fill = realm)) + 
  labs(x = "First year sampled", y = NULL, tag = "F") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = c(0.28, 0.85),
        legend.text = element_text(size = 8),
        legend.key.size = unit(4, "mm"),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

# save the figure
left <- cowplot::plot_grid(NULL)
right <- cowplot::plot_grid(plot_meta_latitude, plot_meta_richness, plot_meta_extent,
                   plot_meta_nsamp, plot_meta_duration, plot_meta_start,
                   nrow = 2, align = "hv") 

cowplot::plot_grid(left, right, nrow = 1, rel_widths = c(0.025, 1)) + 
  cowplot::draw_label("Number of studies", x = 0.01, angle = 90, size = 10)

ggsave("results/Fig.S_histogram_study_characteristics.png", width = 180, height = 125, units = 'mm')



####################
## Figure showing correlations among study characteristics (continuous variables)
library(GGally)

# choose 6 variables used above + protection and log-transform richness, extent, nsamp and duration
dat_meta_continuous <- dat_meta %>%
  dplyr::select(c(realm, study, cent_lat, sprich_total, extent_km2, nsamp_used, 
                  start_year, duration, psamp_inPA_bflate)) %>%
  mutate(cent_lat = abs(cent_lat),
         sprich_total = log10(sprich_total),
         extent_km2 = log10(extent_km2),
         nsamp_used = log10(nsamp_used),
         duration = log10(duration))

# generate figure
corrplot_meta <- ggpairs(dat_meta_continuous, columns = c('cent_lat', 'sprich_total', 'extent_km2', 'nsamp_used', 
                                                          'duration', 'start_year', 'psamp_inPA_bflate'),
                         upper = list(continuous = wrap("cor"), size = 2.81),
                         lower = list(continuous = wrap("points", size = 0.8, alpha = 0.5, shape = 16)),
                         columnLabels = c("|latitude|", "log10_richness", "log10_extent", "log10_n_samples", 
                                          "log10_duration", "start_year", "P_sites_PA")) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  theme_bw() +
  theme(axis.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7.5))

ggsave(corrplot_meta, file = "results/Fig.S_correlation_study_characteristics.png", width = 180, height = 165, units = 'mm', dpi = 600)



##################
## figures showing relationships between study-level slopes and study meta data
# consider: latitude, regional species richness, extent of sites, number of samples, duration of sampling

### firstly fit models considering interaction between realms and study-characteristics
# latitude 
brm_slope_realm_latitude <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + abs(cent_lat):realm + (1|study)),
                    data = brm_oc_aoo10_coef,
                    control = list(adapt_delta = 0.9, max_treedepth = 10),
                    cores = 4, chains = 4, iter = 8000, thin = 4,
                    file = "models/brm_output/brm_slope_realm_latitude")
brm_slope_realm_latitude
pp_check(brm_slope_realm_latitude)
fixef(brm_slope_realm_latitude, probs = c(0.025, 0.1, 0.9, 0.975))


####
# richness 
# use total number of richness with estimates of range sizes
brm_slope_realm_richness <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + log10(sprich_total):realm + (1|study)),
                         data = brm_oc_aoo10_coef,
                         control = list(adapt_delta = 0.9, max_treedepth = 10),
                         cores = 4, chains = 4, iter = 8000, thin = 4,
                         file = "models/brm_output/brm_slope_realm_richness")
brm_slope_realm_richness
pp_check(brm_slope_realm_richness)
fixef(brm_slope_realm_richness, probs = c(0.025, 0.1, 0.9, 0.975))

# check model with all regional species, including those without estimates of range size
# results are consistent. use the above model
brm_slope_realm_richnessAll <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + log10(sprich):realm + (1|study)),
                         data = brm_oc_aoo10_coef,
                         control = list(adapt_delta = 0.9, max_treedepth = 10),
                         cores = 4, chains = 4, iter = 8000, thin = 4,
                         file = "models/brm_output/brm_slope_realm_richnessAll")
brm_slope_realm_richnessAll
pp_check(brm_slope_realm_richnessAll)
fixef(brm_slope_realm_richnessAll, probs = c(0.025, 0.1, 0.9, 0.975))


####
# extent of sites
brm_slope_realm_extent <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + log10(extent_km2):realm + (1|study)),
                         data = brm_oc_aoo10_coef,
                         control = list(adapt_delta = 0.9, max_treedepth = 10),
                         cores = 4, chains = 4, iter = 8000, thin = 4,
                         file = "models/brm_output/brm_slope_realm_extent")
brm_slope_realm_extent
pp_check(brm_slope_realm_extent)
fixef(brm_slope_realm_extent, probs = c(0.025, 0.1, 0.9, 0.975))


####
# number of sites within study
brm_slope_realm_nsamp <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + log10(nsamp_used):realm + (1|study)),
                         data = brm_oc_aoo10_coef,
                         control = list(adapt_delta = 0.9, max_treedepth = 10),
                         cores = 4, chains = 4, iter = 8000, thin = 4,
                         file = "models/brm_output/brm_slope_realm_nsamp")
brm_slope_realm_nsamp
pp_check(brm_slope_realm_nsamp)
fixef(brm_slope_realm_nsamp, probs = c(0.025, 0.1, 0.9, 0.975))


####
# duration of sampling
# use the duration between the first and the last year
brm_slope_realm_duration_log <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + log10(duration):realm + (1|study)),
                         data = brm_oc_aoo10_coef,
                         control = list(adapt_delta = 0.9, max_treedepth = 10),
                         cores = 4, chains = 4, iter = 8000, thin = 4,
                         file = "models/brm_output/brm_slope_realm_duration_log")
brm_slope_realm_duration_log
pp_check(brm_slope_realm_duration_log)
fixef(brm_slope_realm_duration_log, probs = c(0.025, 0.1, 0.9, 0.975))


# use the duration between the mean year in the first and second periods
# results are consistent. use the above model
brm_slope_realm_durationMean_log <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + log10(duration_mean):realm + (1|study)),
                         data = brm_oc_aoo10_coef,
                         control = list(adapt_delta = 0.9, max_treedepth = 10),
                         cores = 4, chains = 4, iter = 8000, thin = 4,
                         file = "models/brm_output/brm_slope_realm_durationMean_log")
brm_slope_realm_durationMean_log
pp_check(brm_slope_realm_durationMean_log)
fixef(brm_slope_realm_durationMean_log, probs = c(0.025, 0.1, 0.9, 0.975))


####
# start year
brm_oc_aoo10_coef <- brm_oc_aoo10_coef %>%
  mutate(start_year1 = start_year - mean(start_year))

brm_slope_realm_startyr <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + start_year1:realm  + (1|study)),
                    data = brm_oc_aoo10_coef,
                    control = list(adapt_delta = 0.9, max_treedepth = 10),
                    cores = 4, chains = 4, iter = 8000, thin = 4,
                    file = "models/brm_output/brm_slope_realm_startyr")
brm_slope_realm_startyr
pp_check(brm_slope_realm_startyr)
fixef(brm_slope_realm_startyr, probs = c(0.025, 0.1, 0.9, 0.975))



## get overall slope and confidence interval, and fitted values
# latitude
# get the fixed effects and 95% and 80% confidence interval, and indicate the significance level
brm_slope_latitude_fixef <- fixef(brm_slope_realm_latitude, prob = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
    filter(term != "intercept") %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":abscent_lat", "", realm),
         sig_level= ifelse(Q2.5>0 | Q97.5 <0, "strong", ifelse(Q10>0 | Q90 <0, "weak", "no")))

# get the fitted values 
slope_latitude_fitted <- cbind(brm_slope_realm_latitude$data,
                               fitted(brm_slope_realm_latitude, re_formula = NA)) %>% 
  as_tibble() %>%
  left_join(brm_slope_latitude_fixef %>% dplyr::select(realm, sig_level))


# richness
# get the fixed effects and 95% and 80% confidence interval, and indicate the significance level
brm_slope_richness_fixef <- fixef(brm_slope_realm_richness, prob = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
  filter(term != "intercept") %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":log10sprich_total", "", realm),
         sig_level= ifelse(Q2.5>0 | Q97.5 <0, "strong", ifelse(Q10>0 | Q90 <0, "weak", "no")))

# get the fitted values 
slope_richness_fitted <- cbind(brm_slope_realm_richness$data,
                               fitted(brm_slope_realm_richness, re_formula = NA)) %>% 
  as_tibble() %>%
  left_join(brm_slope_richness_fixef %>% dplyr::select(realm, sig_level))


# extent
# get the fixed effects and 95% and 80% confidence interval, and indicate the significance level
brm_slope_extent_fixef <- fixef(brm_slope_realm_extent, prob = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
  filter(term != "intercept") %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":log10extent_km2", "", realm),
         sig_level= ifelse(Q2.5>0 | Q97.5 <0, "strong", ifelse(Q10>0 | Q90 <0, "weak", "no")))

# get the fitted values 
slope_extent_fitted <- cbind(brm_slope_realm_extent$data,
                               fitted(brm_slope_realm_extent, re_formula = NA)) %>% 
  as_tibble() %>%
  left_join(brm_slope_extent_fixef %>% dplyr::select(realm, sig_level))


# n samples used
# get the fixed effects and 95% and 80% confidence interval, and indicate the significance level
brm_slope_nsamp_fixef <- fixef(brm_slope_realm_nsamp, prob = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
  filter(term != "intercept") %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":log10nsamp_used", "", realm),
         sig_level= ifelse(Q2.5>0 | Q97.5 <0, "strong", ifelse(Q10>0 | Q90 <0, "weak", "no")))

# get the fitted values 
slope_nsamp_fitted <- cbind(brm_slope_realm_nsamp$data,
                               fitted(brm_slope_realm_nsamp, re_formula = NA)) %>% 
  as_tibble() %>%
  left_join(brm_slope_nsamp_fixef %>% dplyr::select(realm, sig_level))


# duration
# get the fixed effects and 95% and 80% confidence interval, and indicate the significance level
brm_slope_duration_fixef <- fixef(brm_slope_realm_duration_log, prob = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
  filter(term != "intercept") %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":log10duration", "", realm),
         sig_level= ifelse(Q2.5>0 | Q97.5 <0, "strong", ifelse(Q10>0 | Q90 <0, "weak", "no")))

# get the fitted values 
slope_duration_fitted <- cbind(brm_slope_realm_duration_log$data,
                               fitted(brm_slope_realm_duration_log, re_formula = NA)) %>% 
  as_tibble() %>%
  left_join(brm_slope_duration_fixef %>% dplyr::select(realm, sig_level))


# start year
# get the fixed effects and 95% and 80% confidence interval, and indicate the significance level
brm_slope_startyr_fixef <- fixef(brm_slope_realm_startyr, prob = c(0.025, 0.975, 0.1, 0.9)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "realm") %>% 
  mutate(term = rep(c("intercept", "slope"), each = n()/2)) %>%
  filter(term != "intercept") %>%
  mutate(realm = gsub("realm", "", realm),
         realm = gsub(":start_year1", "", realm),
         sig_level= ifelse(Q2.5>0 | Q97.5 <0, "strong", ifelse(Q10>0 | Q90 <0, "weak", "no")))

# get the fitted values 
slope_startyr_fitted <- cbind(brm_slope_realm_startyr$data,
                              fitted(brm_slope_realm_startyr, re_formula = NA)) %>% 
  as_tibble() %>%
  left_join(brm_slope_startyr_fixef %>% dplyr::select(realm, sig_level)) %>%
  mutate(start_year = start_year1 + mean(brm_oc_aoo10_coef$start_year),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


#### generate figure
plot_slope_latitude  <- ggplot() +
  geom_point(data = brm_oc_aoo10_coef, 
             aes(x = abs(cent_lat), y = estimate_slope, colour = realm), 
             shape = 19, size = 1, alpha = 0.8) +  
  geom_ribbon(data = slope_latitude_fitted, 
              aes(x = abs(cent_lat), ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.4) + 
  geom_line(data = slope_latitude_fitted, 
            aes(x = abs(cent_lat), y = Estimate, colour = realm), size = 1) +  
  geom_hline(yintercept = 0, linetype = 2, size = 0.3, colour = "black") +
  labs(x = "|Latitude|", y = NULL) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = "no",
        text = element_text(size = 8),
        axis.title = element_text(size = 9))

plot_slope_richness  <- ggplot() +
  geom_point(data = brm_oc_aoo10_coef, 
             aes(x = sprich_total, y = estimate_slope, colour = realm), 
             shape = 19, size = 1, alpha = 0.8) +  
  geom_ribbon(data = slope_richness_fitted, 
              aes(x = sprich_total, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.4) + 
  geom_line(data = slope_richness_fitted, 
            aes(x =sprich_total, y = Estimate, colour = realm), size = 1) +  
  geom_hline(yintercept = 0, linetype = 2, size = 0.3, colour = "black") +
  scale_x_log10() + 
  labs(x = "Regional species richness", y = NULL) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = "no",
        text = element_text(size = 8),
        axis.title = element_text(size = 9))

plot_slope_extent  <- ggplot() +
  geom_point(data = brm_oc_aoo10_coef, 
             aes(x = extent_km2, y = estimate_slope, colour = realm), 
             shape = 19, size = 1, alpha = 0.8) +  
  geom_ribbon(data = slope_extent_fitted, 
              aes(x = extent_km2, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.4) + 
  geom_line(data = slope_extent_fitted, 
            aes(x =extent_km2, y = Estimate, colour = realm), size = 1) +  
  geom_hline(yintercept = 0, linetype = 2, size = 0.3, colour = "black") +
  scale_x_log10() + 
  labs(x = bquote('Extent of study sites ('*km^2*')'), y = NULL) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = "no",
        text = element_text(size = 8),
        axis.title = element_text(size = 9))

plot_slope_nsamp  <- ggplot() +
  geom_point(data = brm_oc_aoo10_coef, 
             aes(x = nsamp_used, y = estimate_slope, colour = realm), 
             shape = 19, size = 1, alpha = 0.8) +  
  geom_ribbon(data = slope_nsamp_fitted, 
              aes(x = nsamp_used, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.4) + 
  geom_line(data = slope_nsamp_fitted, 
            aes(x =nsamp_used, y = Estimate, colour = realm), size = 1) +  
  geom_hline(yintercept = 0, linetype = 2, size = 0.3, colour = "black") +
  scale_x_log10() + 
  labs(x = "Number of samples", y = NULL) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = "no",
        text = element_text(size = 8),
        axis.title = element_text(size = 9))

plot_slope_duration  <- ggplot() +
  geom_point(data = brm_oc_aoo10_coef, 
             aes(x = duration, y = estimate_slope, colour = realm), 
             shape = 19, size = 1, alpha = 0.8) +  
  geom_ribbon(data = slope_duration_fitted, 
              aes(x = duration, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.4) + 
  geom_line(data = slope_duration_fitted, 
            aes(x =duration, y = Estimate, colour = realm), size = 1) +  
  geom_hline(yintercept = 0, linetype = 2, size = 0.3, colour = "black") +
  scale_x_log10() + 
  labs(x = "Duration of sampling (years)", y = NULL) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = "no",
        text = element_text(size = 8),
        axis.title = element_text(size = 9))

plot_slope_startyr  <- ggplot() +
  geom_point(data = brm_oc_aoo10_coef, 
             aes(x = start_year, y = estimate_slope, colour = realm), 
             shape = 19, size = 1, alpha = 0.8) +  
  geom_ribbon(data = slope_startyr_fitted, 
              aes(x = start_year, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.4) + 
  geom_line(data = slope_startyr_fitted, 
            aes(x =start_year, y = Estimate, colour = realm), size = 1) +  
  geom_hline(yintercept = 0, linetype = 2, size = 0.3, colour = "black") +
  #scale_x_log10() + 
  labs(x = "Start year", y = NULL) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = c(0.25, 0.85),
        legend.text = element_text(size = 8),
        legend.key.size = unit(4, "mm"),
        text = element_text(size = 8),
        axis.title = element_text(size = 9))

# save the figure 
left <- cowplot::plot_grid(NULL)
right <- cowplot::plot_grid(plot_slope_latitude, plot_slope_richness, plot_slope_extent,
                            plot_slope_nsamp, plot_slope_duration, plot_slope_startyr,
                            nrow = 2, align = "hv",
                            labels = c("A", "B", "C", "D", "E", "F"), label_size = 9) 

cowplot::plot_grid(left, right, nrow = 1, rel_widths = c(0.025, 1)) + 
  cowplot::draw_label("Effect of range size on occupancy change", x = 0.01, angle = 90, size = 10)


ggsave("results/Fig.S_slope_covariates.png", width = 180, height = 125, units = 'mm', dpi =600)

