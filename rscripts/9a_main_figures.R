## draw Fig. 2

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


dat_meta <- dat_meta %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

brm_oc_aoo10_coef <- brm_oc_aoo10_coef %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

brm_oc_aoo10_realm_coef <- brm_oc_aoo10_realm_coef %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


# summary of used studies for table S1 
study_used <- brm_oc_aoo10_coef %>% dplyr::select(study, database, studyID, study_name, start_year, end_year, n_years, realm, taxon, taxon_new, region, psamp_inPA_bflate, estimate_slope)


# save study-level coefficients
write.csv(study_used, file = "results/Table_datasets.csv", row.names=FALSE)
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
  mutate(estimate_slope = ifelse(estimate_slope <= -0.06, -0.06, ifelse(estimate_slope >= 0.06, 0.06, estimate_slope))) %>%
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
  labs(shape = "Realm") +
  theme_bw() +
  theme(legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, hjust = 0.2)) + 
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
  scale_fill_gradient2(breaks = c(-0.05, 0, 0.05), low = "#155f49", mid = "#f7ffda", high = "#d17538") +
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
  scale_fill_gradient2(low = "#155f49", mid = "#f7ffda", high = "#d17538") +
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
  scale_fill_gradient2(low = "#155f49", mid = "#f7ffda", high = "#d17538") +
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
realm_col <- c("All" = "#000000", "Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")

brm_oc_aoo10_realm_all_fitted <- brm_oc_aoo10_realm_fitted %>%
  bind_rows(brm_oc_aoo10_fitted %>% mutate(realm = "All")) %>%
  mutate(realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))

brm_oc_aoo10_realm_all_line <- brm_oc_aoo10_realm_line %>%
  bind_rows(brm_oc_aoo10_line %>% mutate(realm = "All")) %>% 
  mutate(realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))

# generate figure
plot_oc_aoo <- ggplot() + 
  geom_ribbon(data = brm_oc_aoo10_realm_all_fitted,
              aes(x = aoo10, ymin = oc_sqroot_Q2.5, ymax = oc_sqroot_Q97.5, fill = realm), 
              size =1, alpha = 0.3) + 
  geom_segment(data = brm_oc_aoo10_realm_all_line, 
               aes(x = xmin, xend = xmax, 
                   y = intercept + slope * cl.xmin, 
                   yend = intercept + slope * cl.xmax,
                   colour = realm), 
               size = 1, alpha = 1) + 
  geom_hline(yintercept = 0, lty = 2, size = 0.3) + 
  scale_x_log10() + 
  scale_y_continuous(breaks = c(-0.08, -0.04, 0, 0.04, 0.08)) +
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  labs(x = "Range size (number of 10-km grid-cells)", y = "Occupancy change", colour = "", fill = "") +
  theme_classic() +
  theme(legend.position = c(0.9, 0.2),
        legend.margin = margin(0, 0, 0, -45),
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))

ggsave(plot_oc_aoo, file = "results/Fig2b.png", width = 82, height = 82, units = 'mm')


## Figure 2c: frequency distribution of study-level slopes
# all studies 

brm_oc_aoo10_coef <- brm_oc_aoo10_coef %>% 
  mutate(signal_slope = ifelse(estimate_slope > 0, "+ ns", "- ns"),
         signal_slope = ifelse(Q2.5_slope > 0, "+ sig", ifelse(Q97.5_slope <= 0, "- sig", signal_slope)),
         signal_slope = factor(signal_slope, levels = c("- sig","- ns", "+ ns", "+ sig")))

plot_histogram_slopes <- ggplot(brm_oc_aoo10_coef) + 
  geom_histogram(aes(estimate_slope, fill = signal_slope)) +
  geom_vline(xintercept = 0, linetype =2, colour = "black") +
  geom_vline(data = brm_oc_aoo10_fixed %>% filter(term == "slope"), 
             aes(xintercept = Estimate), linetype =1) + 
  geom_rect(data = brm_oc_aoo10_fixed %>% filter(term == "slope"), 
            aes(xmin = Q2.5, xmax = Q97.5, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "gray50") + 
  annotate(geom = "text", x =0.05, y = 20, label = "All", size = 2.81) + 
  labs(x = "", y = "Number of studies", fill = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) +
  theme_classic()+
  theme(legend.position = c(0.2,0.8), 
        legend.key.size = unit(3, 'mm'),
        axis.title = element_text(size = 8, face = "plain"),
        axis.text = element_text(size = 6, face = "plain")) + 
  scale_fill_manual(values = c("- ns" = "#75baa5", "+ ns" = "#e8b797", "- sig" = "#155f49", "+ sig" = "#d17538"))
  # scale_fill_manual(values = c("neutral" = "dark gray", "negative" = "#155f49", "positive" = "#d17538"))
  # scale_fill_manual(values = c("negative" = "#75baa5", "positive" = "#e8b797", "negative_sig" = "#155f49", "positive_sig" = "#d17538"))
  #scale_fill_manual(values = c("neutral" = "light gray", "negative" = "orange", "positive" = "skyblue"))
  #scale_fill_manual(values = c("neutral" = "light gray", "negative" = "blue", "positive" = "red"))


# different realms separately
brm_oc_aoo10_realm_slope <- brm_oc_aoo10_realm_fixed %>% 
  filter(term == "slope") %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

brm_oc_aoo10_realm_coef <- brm_oc_aoo10_realm_coef %>% 
  mutate(signal_slope = ifelse(estimate_slope > 0, "+ ns", "- ns"),
         signal_slope = ifelse(Q2.5_slope > 0, "+ sig", ifelse(Q97.5_slope <= 0, "- sig", signal_slope)),
         signal_slope = factor(signal_slope, levels = c("- sig","- ns", "+ ns", "+ sig")))


realm_labs <- data.frame(realm = factor(c("Terrestrial", "Freshwater", "Marine")),
                         label = c("Terrestrial", "Freshwater", "Marine"),
                         x = c(0.07, 0.07, 0.07), y = c(13, 15, 14))

plot_histogram_slopes_realms <- ggplot(brm_oc_aoo10_realm_coef) + 
  facet_wrap(~ realm, ncol =1, scale = "free") + 
  geom_histogram(aes(estimate_slope, fill = signal_slope), bins = 20 ) +
  geom_vline(xintercept = 0, linetype =2, colour = "black") +
  geom_vline(data =brm_oc_aoo10_realm_slope, 
             aes(xintercept = Estimate), linetype =1) + 
  geom_rect(data = brm_oc_aoo10_realm_slope, 
            aes(xmin = Q2.5, xmax = Q97.5, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "gray50") + 
  geom_text(data = realm_labs, aes(x, y, label = label), size = 2.11 ) +
  labs(x = "", y = NULL) + 
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) + 
  scale_x_continuous(limits = 1.15*range(brm_oc_aoo10_coef$estimate_slope), breaks = c(-0.05, 0, 0.05, 0.10)) +
  theme_classic() + 
  theme(legend.position = "n", 
        #plot.margin = unit(c(2,0.5,2,0), "mm"),
        axis.title = element_text(size = 8, face = "plain"),
        axis.text = element_text(size = 6, face = "plain"),
        strip.text = element_text(size = 8),
        line = element_line(size = 0.3),
        strip.text.x  = element_blank(),
        strip.background = element_blank()) + 
  scale_fill_manual(values = c("- ns" = "#75baa5", "+ ns" = "#e8b797", "- sig" = "#155f49", "+ sig" = "#d17538"))
  #scale_fill_manual(values = c("neutral" = "light gray", "negative" = "orange", "positive" = "skyblue"))
  #scale_fill_manual(values = c("neutral" = "light gray", "negative" = "blue", "positive" = "red"))


plot_slope_histogram <- cowplot::plot_grid(plot_histogram_slopes, plot_histogram_slopes_realms, 
                                           ncol = 2, rel_widths =  c(2, 1)) + 
  cowplot::draw_label("Effect of range size", y = 0.05, x= 0.6, size = 8)

ggsave(plot_slope_histogram, file = "results/Fig2c.png", width = 120, height = 82, units = 'mm')


# combine panels of figure 2
bottom <- cowplot::plot_grid(plot_oc_aoo, plot_slope_histogram, nrow = 1, 
                   rel_widths =  c(1, 1.3), 
                   labels = c('b', 'c'), label_size = 10)

cowplot::plot_grid(plot_re_map, bottom, nrow = 2, 
                   rel_heights =  c(1, 1.1), 
                   labels = c('a', ""), label_size = 10)

ggsave(file = "results/Fig2.pdf", 
       width = 179, height = 135, units = 'mm')
ggsave(file = "results/Fig2.png", 
       width = 179, height = 135, units = 'mm', dpi = 600) # 2.8:2.1
