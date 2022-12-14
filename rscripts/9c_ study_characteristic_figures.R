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


realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")


############
## Figure showing frequency distribution of study characteristics
# cent_lat, sprich, extent_km2, nsample_used, duration, start_year

plot_meta_latitude <- ggplot(data = dat_meta) +
  geom_histogram(aes(cent_lat, fill = realm)) + 
  #scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) + 
  labs(x = "Latitude (degree)", y = NULL, fill = NULL, tag = "a") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_richness <- ggplot(data = dat_meta) +
  geom_histogram(aes(sprich, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Regional species richness", y = NULL, tag = "b") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_extent <- ggplot(data = dat_meta) +
  geom_histogram(aes(extent_km2, fill = realm)) + 
  scale_x_log10(breaks = c(10^0, 10^3, 10^6)) + 
  labs(x = bquote('Extent of study sites ('*km^2*')'), y = NULL, tag = "c") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_nsamp <- ggplot(data = dat_meta) +
  geom_histogram(aes(nsamp_used, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Number of samples", y = NULL, tag = "d") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_duration <- ggplot(data = dat_meta) +
  geom_histogram(aes(duration, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Duration of sampling (years)", y = NULL, tag = "e") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_start <- ggplot(data = dat_meta) +
  geom_histogram(aes(start_year, fill = realm)) + 
  labs(x = "First year sampled", y = NULL, tag = "f") + 
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
  dplyr::select(c(realm, study, cent_lat, sprich, extent_km2, nsamp_used, 
                  start_year, duration, psamp_inPA_bflate)) %>%
  mutate(cent_lat = abs(cent_lat),
         sprich = log10(sprich),
         extent_km2 = log10(extent_km2),
         nsamp_used = log10(nsamp_used),
         duration = log10(duration))

# generate figure
corrplot_meta <- ggpairs(dat_meta_continuous, columns = c('cent_lat', 'sprich', 'extent_km2', 'nsamp_used', 
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

