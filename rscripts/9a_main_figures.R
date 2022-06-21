## draw Fig. 2 and Fig. 3, 
# and prepare the supplement figures & tables related to protection status


rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ggplot2", "ggridges", "brms", "tidybayes", "cowplot","sf", "rnaturalearth" )

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


load("models/data_input_to_models.RDATA")
load("results/coefs_main_brms.RDATA")

# an error in the raw data: the year 2015 was incorrectly input as 2915. Correct here
id <- dat_meta$study_name == "price_2017_Illinois" & !is.na(dat_meta$study_name)
dat_meta$end_year[id] <- 2015
dat_meta$duration[id] <- 16
dat_meta$duration_mean[id] <- 16
id <- brm_oc_aoo10_coef$study_name == "price_2017_Illinois" & !is.na(brm_oc_aoo10_coef$study_name)
brm_oc_aoo10_coef$end_year[id] <- 2015
brm_oc_aoo10_coef$duration[id] <- 16
brm_oc_aoo10_coef$duration_mean[id] <- 16

dat_meta <- dat_meta %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

brm_oc_aoo10_coef <- brm_oc_aoo10_coef %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  left_join(oc_period %>% distinct(study, nsamp_used))


# summary of used studies for table S1 
study_used <- brm_oc_aoo10_coef %>% dplyr::select(study, database, studyID, study_name, start_year, end_year, n_years, realm, taxon, taxon_new, region, psamp_inPA_bflate, estimate_slope)


# save study-level coefficients
write.csv(study_used, file = "results/Table_S1.csv", row.names=FALSE)
write.csv(brm_oc_aoo10_coef, file = "results/Study_summary.csv", row.names=FALSE)



############
## Figure 2a: the map showing study-level slopes, with zoomed North America and Europe
#  https://www.r-bloggers.com/2019/04/zooming-in-on-maps-with-sf-and-ggplot2/

# world map 
worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')

# generate a spatial point file containing estimated slopes (random effect)
oc_aoo10_re <- brm_oc_aoo10_coef %>% 
  dplyr::select(cent_long, cent_lat, estimate_slope, realm) %>% 
  mutate(estimate_slope = ifelse(estimate_slope <= -1, -1, ifelse(estimate_slope >= 1, 1, estimate_slope))) %>%
  st_as_sf(coords = c('cent_long', 'cent_lat'))  %>% 
  st_set_crs(4326) %>%
  st_transform(crs = "+proj=moll")

# bounding of north america
namerica <- st_sfc(st_point(c(-115, 21)), st_point(c(-65, 55)), crs = 4326) %>%
  st_transform(crs = "+proj=moll") %>% 
  st_bbox() %>%
  st_as_sfc()

# bounding of Europe
europe <- st_sfc(st_point(c(-17, 38)), st_point(c(30, 60)), crs = 4326) %>%
  st_transform(crs = "+proj=moll") %>% 
  st_bbox() %>%
  st_as_sfc()

# plot to generate legend of realm
plot_legend_realm <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray60", colour = "gray50", size = 0.2) + 
  geom_sf(data = namerica, fill = NA, colour = "gray50") + 
  geom_sf(data = europe, fill = NA, colour = "gray50") + 
  geom_sf(data = oc_aoo10_re, aes(fill = estimate_slope, shape = realm), 
          size = 1.5, alpha = 0.8) + 
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) + 
  labs(shape = NULL) +
  theme_bw() +
  theme(legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 8)) + 
  guides(fill = FALSE) + 
  guides(shape = guide_legend(ncol = 1))

# get legends
legend_realm <- get_legend(plot_legend_realm)

# global map
re_map_global <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray80", colour = "gray50", size = 0.2) + 
  geom_sf(data = namerica, fill = NA, colour = "gray50") + 
  geom_sf(data = europe, fill = NA, colour = "gray50") + 
  geom_sf(data = oc_aoo10_re, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 1.5, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), ylim = c(-8172663, 9020048), expand = FALSE) + 
  #scale_colour_gradient2(breaks = c(-0.8, -0.4, 0, 0.4, 0.8), low = "blue", mid = "white", high = "red") +
  #scale_colour_gradient2(breaks = c(-0.8, -0.4, 0, 0.4, 0.8), low = "#155f49", mid = "#f7ffda", high = "#d17538") +
  scale_fill_gradient2(breaks = c(-0.8, -0.4, 0, 0.4, 0.8), low = "#155f49", mid = "#f7ffda", high = "#d17538") +
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) + 
  annotate("text", x =1.065*-11046300, y= 0.9*6386580, label = "NA", size = 2.81) +
  annotate("text", x =1.4*-1469518, y= 0.95*6876759, label = "EU", size = 2.81) +
  theme_bw() +
  theme(plot.margin = unit(c(0,-1,0,-1), "cm"),
        legend.position = c(0.55, 0.15), 
        legend.direction = 'horizontal',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.background = element_rect(fill = 'transparent'),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  guides(shape = FALSE) + 
  guides(fill = guide_colourbar(title = "Slope", title.position = "top", barwidth = 5))


re_map_global <- ggdraw(re_map_global) +
  draw_plot(legend_realm, width = 0.1, height = 0.1, x = 0.15, y = 0.15)

ggsave(file = "results/Fig2a_global.png", 
       width = 140, height = 70, units = 'mm', dpi = 600)

# map of north america
re_map_namerica <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray80", colour = "gray50", size = 0.2) + 
  geom_sf(data = oc_aoo10_re, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 2, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), expand = FALSE, 
           xlim = st_coordinates(namerica)[1:2, "X"], ylim = st_coordinates(namerica)[2:3, "Y"]) + 
  annotate("text", x =0.95*-11046300, y= 0.9*6386580, label = "NA", size = 2.81) +
  #scale_colour_gradient2(breaks = c(-0.8, -0.4, 0, 0.4, 0.8), low = "blue", mid = "white", high = "red") +
  scale_fill_gradient2(breaks = c(-0.8, -0.4, 0, 0.4, 0.8), low = "#155f49", mid = "#f7ffda", high = "#d17538") +
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "mm"),
        legend.position = "n", 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# map of Europe
re_map_europe <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray80", colour = "gray50", size = 0.2) + 
  geom_sf(data = oc_aoo10_re, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 2, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), expand = TRUE, 
           xlim = st_coordinates(europe)[1:2, "X"], ylim = st_coordinates(europe)[2:3, "Y"]) + 
  annotate("text", x =0.85*-1469518, y= 0.95*6876759, label = "EU", size = 2.81) +
  #scale_colour_gradient2(breaks = c(-0.8, -0.4, 0, 0.4, 0.8), low = "blue", mid = "white", high = "red") + 
  scale_fill_gradient2(breaks = c(-0.8, -0.4, 0, 0.4, 0.8), low = "#155f49", mid = "#f7ffda", high = "#d17538") +
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "mm"),
        legend.position = "n", 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

right <- cowplot::plot_grid(re_map_namerica, re_map_europe, ncol =1, nrow = 2, rel_heights = c(1, 1.15))
plot_re_map <- cowplot::plot_grid(re_map_global, right, ncol =2, nrow = 1, rel_widths = c(2.52, 1))

ggsave(plot_re_map, file = "results/Fig2a.png", 
       width = 180, height = 64, units = 'mm', dpi = 600) # 2.8:1


## Figure 2b: relationships between occupancy change and range size fro all combined and three realms
# colors for realm
realm_col <- c("All" = "black", "Terrestrial" = "#D95F02", "Freshwater" = "#1B9E77", "Marine" = "#7570B3")

brm_oc_aoo10_realm_all_fitted <- brm_oc_aoo10_realm_fitted %>%
  bind_rows(brm_oc_aoo10_fitted %>% mutate(realm = "All")) %>%
  mutate(realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))

brm_oc_aoo10_realm_all_line <- brm_oc_aoo10_realm_line %>%
  bind_rows(brm_oc_aoo10_line %>% mutate(realm = "All")) %>% 
  mutate(realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))

# generate figure
plot_oc_aoo <- ggplot() + 
  geom_ribbon(data = brm_oc_aoo10_realm_all_fitted,
              aes(x = aoo10, ymin = occup_change_logit_Q2.5, ymax = occup_change_logit_Q97.5, fill = realm), 
              size =1, alpha = 0.3) + 
  geom_segment(data = brm_oc_aoo10_realm_all_line, 
               aes(x = xmin, xend = xmax, 
                   y = intercept + slope * cl.xmin, 
                   yend = intercept + slope * cl.xmax,
                   colour = realm), 
               size = 1, alpha = 1) + 
  geom_hline(yintercept = 0, lty = 2, size = 0.3) + 
  scale_x_log10() + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  labs(x = "Range size (number of 10-km grid-cells)", y = "Occupancy change", colour = "", fill = "") +
  theme_classic() +
  theme(legend.position = c(0.9, 0.25),
        legend.margin = margin(0, 0, 0, -45),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))

ggsave(plot_oc_aoo, file = "results/Fig2b.png", width = 82, height = 82, units = 'mm')


## Figure 2c: frequency distribution of study-level slopes
# all studies
plot_histogram_slopes <- ggplot(brm_oc_aoo10_coef) + 
  geom_histogram(aes(estimate_slope, fill = sig_slope) ) +
  geom_vline(xintercept = 0, linetype =2, colour = "black") +
  geom_vline(data = brm_oc_aoo10_fixed %>% filter(term == "slope"), 
             aes(xintercept = Estimate), linetype =1) + 
  geom_rect(data = brm_oc_aoo10_fixed %>% filter(term == "slope"), 
            aes(xmin = Q2.5, xmax = Q97.5, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "gray50") + 
  annotate(geom = "text", x =1, y = 25, label = "All", size = 2.81) + 
  labs(x = "", y = "Number of studies") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) +
  theme_classic()+
  theme(legend.position = "n", 
        axis.title = element_text(size = 8, face = "plain"),
        axis.text = element_text(size = 6, face = "plain")) + 
  scale_fill_manual(values = c("neutral" = "light gray", "negative" = "#155f49", "positive" = "#d17538"))
  #scale_fill_manual(values = c("neutral" = "light gray", "negative" = "orange", "positive" = "skyblue"))
  #scale_fill_manual(values = c("neutral" = "light gray", "negative" = "blue", "positive" = "red"))


# different realms separately
brm_oc_aoo10_realm_slope <- brm_oc_aoo10_realm_fixed %>% 
  filter(term == "slope") %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

realm_labs <- data.frame(realm = factor(c("Terrestrial", "Freshwater", "Marine")),
                         label = c("Terrestrial", "Freshwater", "Marine"),
                         x = c(1.5, 1.5, 1.5), y = c(14, 14, 14))

plot_histogram_slopes_realms <- ggplot(brm_oc_aoo10_coef) + 
  facet_wrap(~ realm, ncol =1, scale = "free") + 
  geom_histogram(aes(estimate_slope, fill = sig_slope), bins = 20 ) +
  geom_vline(xintercept = 0, linetype =2, colour = "black") +
  geom_vline(data =brm_oc_aoo10_realm_slope, 
             aes(xintercept = Estimate), linetype =1) + 
  geom_rect(data = brm_oc_aoo10_realm_slope, 
            aes(xmin = Q2.5, xmax = Q97.5, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "gray50") + 
  geom_text(data = realm_labs, aes(x, y, label = label), size = 2.46 ) +
  labs(x = "", y = NULL) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) + 
  scale_x_continuous(limits = 1.15*range(brm_oc_aoo10_coef$estimate_slope), breaks = c( -1, 0, 1, 2)) +
  theme_classic() + 
  theme(legend.position = "n", 
        #plot.margin = unit(c(2,0.5,2,0), "mm"),
        axis.title = element_text(size = 8, face = "plain"),
        axis.text = element_text(size = 6, face = "plain"),
        strip.text = element_text(size = 8),
        line = element_line(size = 0.3),
        strip.text.x  = element_blank(),
        strip.background = element_blank()) + 
  scale_fill_manual(values = c("neutral" = "light gray", "negative" = "#155f49", "positive" = "#d17538"))
  #scale_fill_manual(values = c("neutral" = "light gray", "negative" = "orange", "positive" = "skyblue"))
  #scale_fill_manual(values = c("neutral" = "light gray", "negative" = "blue", "positive" = "red"))


plot_slope_histogram <- cowplot::plot_grid(plot_histogram_slopes, plot_histogram_slopes_realms, 
                                           ncol = 2, rel_widths =  c(2, 1)) + 
  cowplot::draw_label("Effect of range size", y = 0.05, x= 0.6, size = 8)

ggsave(plot_slope_histogram, file = "results/Fig2c.png", width = 120, height = 82, units = 'mm')


# combine panels of figure 2
bottom <- cowplot::plot_grid(plot_oc_aoo, plot_slope_histogram, nrow = 1, 
                   rel_widths =  c(1, 1.3), 
                   labels = c('B', 'C'), label_size = 10)

cowplot::plot_grid(plot_re_map, bottom, nrow = 2, 
                   rel_heights =  c(1, 1.1), 
                   labels = c('A', ""), label_size = 10)

ggsave(file = "results/Fig2.pdf", 
       width = 180, height = 135, units = 'mm')
ggsave(file = "results/Fig2.png", 
       width = 180, height = 135, units = 'mm', dpi = 600) # 2.8:2.1





############
## figure 3: relationship between study-level slope and percentage of sites in protected area

## fit models firstly: include uncertainty of the study-level estimate in the model 
# use the proportion of sites in protected areas that established before the late period of sampling
brm_slope_pa_bflate_realm <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + psamp_inPA_bflate:realm + (1|study)),
                          data = brm_oc_aoo10_coef %>% filter(!is.na(psamp_inPA_bflate)),
                          control = list(adapt_delta = 0.9, max_treedepth = 10),
                          cores = 4, chains = 4, iter = 8000, thin = 4,
                          file = "models/brm_output/brm_slope_pa_bflate_realm")
brm_slope_pa_bflate_realm
pp_check(brm_slope_pa_bflate_realm)
fixef(brm_slope_pa_bflate_realm, prob = c(0.025, 0.975, 0.1, 0.9))


## use the proportion of sites in any protected areas
brm_slope_pa_wdpa_realm <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + psamp_inPA_wdpa:realm + (1|study)),
                               data = brm_oc_aoo10_coef %>% filter(!is.na(psamp_inPA_wdpa)),
                               control = list(adapt_delta = 0.9, max_treedepth = 10),
                               cores = 4, chains = 4, iter = 8000, thin = 4,
                               file = "models/brm_output/brm_slope_pa_wdpa_realm")
brm_slope_pa_wdpa_realm
pp_check(brm_slope_pa_wdpa_realm)
fixef(brm_slope_pa_wdpa_realm, prob = c(0.025, 0.975, 0.1, 0.9))

# get the fitted values for each level of PA
# protected areas established before the late of sampling years
slope_pa_bflate_realm_fitted <- cbind(brm_slope_pa_bflate_realm$data,
                                      fitted(brm_slope_pa_bflate_realm, re_formula = NA)) %>% 
  as_tibble() %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

# all protected areas 
slope_pa_wdpa_realm_fitted <- cbind(brm_slope_pa_wdpa_realm$data,
                                    fitted(brm_slope_pa_wdpa_realm, re_formula = NA)) %>% 
  as_tibble() %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


## generate figure
# colors for realms
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#1B9E77", "Marine" = "#7570B3")
realm_labs <- data.frame(realm = factor(c("Terrestrial", "Freshwater", "Marine")),
                         label = c("Terrestrial", "Freshwater", "Marine"),
                         x = c(0.5, 0.5, 0.5), y = c(2.2, 2.2, 2.2))

# protected areas established before the late of sampling years
ggplot(slope_pa_bflate_realm_fitted) + 
  facet_wrap(~ realm, scales = "free") + 
  geom_point(aes(x = psamp_inPA_bflate, y = estimate_slope, color = realm), 
             shape = 19, size = 0.8, alpha = 0.8) +
  geom_ribbon(aes(x = psamp_inPA_bflate, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.5) + 
  geom_line(aes(x = psamp_inPA_bflate, y = Estimate, color = realm), size = 0.6) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  scale_y_continuous(limits = c(-1.2, 2.55), breaks = c( -1, 0, 1, 2)) + 
  geom_text(data = realm_labs, aes(x, y, label = label), size = 2.81) +
  labs(y = "Effect of range size", x = "Proportion of sites in protected areas") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = "n", 
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        strip.text.x  = element_blank(),
        line = element_line(size = 0.4),
        strip.background = element_blank())

ggsave("results/Fig3.pdf", width = 125, height = 50, units = 'mm')
ggsave("results/Fig3.png", width = 125, height = 50, units = 'mm', dpi = 600)


# supplement figure using all protected areas 
ggplot(slope_pa_wdpa_realm_fitted) + 
  facet_wrap(~ realm, scales = "free") + 
  geom_point(aes(x = psamp_inPA_wdpa, y = estimate_slope, color = realm), 
             shape = 19, size = 0.8, alpha = 0.8) +
  geom_ribbon(aes(x = psamp_inPA_wdpa, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.5) + 
  geom_line(aes(x = psamp_inPA_wdpa, y = Estimate, color = realm), size = 0.6) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  scale_y_continuous(limits = c(-1.2, 2.55), breaks = c( -1, 0, 1, 2)) + 
  geom_text(data = realm_labs, aes(x, y, label = label), size = 2.81) +
  labs(y = "Effect of range size", x = "Proportion of sites in protected areas") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = "n", 
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        strip.text.x  = element_blank(),
        line = element_line(size = 0.4),
        strip.background = element_blank())

ggsave("results/Fig.S_slope_PA_wdpa.png", width = 125, height = 50, units = 'mm', dpi = 600)


####################
## regress study-level slopes against P. of sites in protected areas and some other study characteristics together

# center predictor variables
brm_oc_aoo10_coef <- brm_oc_aoo10_coef %>%
  mutate(pa_cent = psamp_inPA_bflate - mean(psamp_inPA_bflate, na.rm= TRUE),
         lat_cent = abs(cent_lat) - mean(abs(cent_lat)),
         sprich_cent = log10(sprich) - mean(log10(sprich)),       
         duration_cent = log10(duration) - mean(log10(duration)),
         startyr_cent = start_year - mean(start_year),
         extent_cent = log10(extent_km2) - mean(log10(extent_km2)),
         nsamp_cent = log10(nsamp_used) - mean(log10(nsamp_used)))

brm_slope_realm_meta <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + psamp_inPA_bflate:realm +
                                 lat_cent:realm + sprich_cent:realm + 
                                 duration_cent:realm + startyr_cent:realm + 
                                 extent_cent:realm + nsamp_cent:realm + (1|study)),
                            data = brm_oc_aoo10_coef,
                            control = list(adapt_delta = 0.9, max_treedepth = 10),
                            cores = 4, chains = 4, iter = 8000, thin = 4,
                            file = "models/brm_output/brm_slope_realm_meta")
brm_slope_realm_meta
pp_check(brm_slope_realm_meta)

# get the fixed effects and save as a table 
fixef_slope_study.character <- fixef(brm_slope_realm_meta, probs = c(0.025,  0.975)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:Q97.5, round, 3))

write.csv(fixef_slope_study.character, file = "results/Table.S_slope_study.characters.csv")



#################
## compare proportion of sites in any protected areas and those established before the later period

## relationship between proportion of sites in protected areas estimated using all PA and those established before the later period
compare_PA_wdpa_bflate <- ggplot(dat_meta %>% filter(!is.na(psamp_inPA_wdpa)) %>%
         mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))) + 
  facet_wrap(~ realm, scales = "free") + 
  geom_point(aes(y = psamp_inPA_wdpa, x = psamp_inPA_bflate, color = realm), 
             shape = 19, size = 0.8, alpha = 0.8) +
  geom_abline(intercept = 0, slope =1, lty =2, size = 0.3) + 
  labs(y = "Proportion of sites in protected areas \nestablished at any time points", 
       x = "Proportion of sites in protected areas established before later period",
       tag = "A") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  #theme_bw()+ 
  theme(legend.position = "n", 
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        strip.background = element_blank())

histogram_PA_bflate <- ggplot(dat_meta %>% filter(!is.na(psamp_inPA_wdpa))) + 
  facet_wrap(~ realm, scales = "free") + 
  geom_histogram(aes(psamp_inPA_bflate, fill = realm)) +
  labs(x = "Proportion of sites in protected areas established before later period", y = "Number of studies", tag = "B") +
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic()+
  theme(legend.position = "n", 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        strip.text.x  = element_blank()) 

histogram_PA_wdpa <- ggplot(dat_meta %>% filter(!is.na(psamp_inPA_wdpa))) + 
  facet_wrap(~ realm, scales = "free") + 
  geom_histogram(aes(psamp_inPA_wdpa, fill = realm)) +
  labs(x = "Proportion of sites in protected areas established at any time points", y = "Number of studies", tag = "C") +
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic()+
  theme(legend.position = "n", 
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        plot.tag = element_text(size = 8, face = 'bold'),
        strip.text.x  = element_blank()) 

cowplot::plot_grid(compare_PA_wdpa_bflate, histogram_PA_bflate, histogram_PA_wdpa,
                   nrow= 3, align = "v")
ggsave("results/Fig.S_comopare_PA_all_before.late.png", width = 140, height = 160, unit = "mm")
