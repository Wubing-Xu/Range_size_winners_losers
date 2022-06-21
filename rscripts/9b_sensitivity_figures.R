## draw figures showing results from sensitivity analyses using different estimates of range sizes and different subsets of assemblages or species

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ggplot2", "ggridges", "brms", "tidybayes", "cowplot")

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


############
## Figure showing frequency distribution of study-level slopes between occupancy change and range size
# using different measures of range size

# combine study-level slopes based on different range size measures
brm_oc_range_coef <- rbind(brm_oc_aoo50_coef %>% mutate(rangetype = "AOO_50km"),
      brm_oc_aoo100_coef %>% mutate(rangetype = "AOO_100km"),
      brm_oc_ahull6_coef %>% mutate(rangetype = "Alpha hull")) %>%
  mutate(rangetype = factor(rangetype, levels = c("AOO_50km", "AOO_100km", "Alpha hull")))

# combine global fixed slopes based on different range size measures
brm_oc_range_fixed <- rbind(brm_oc_aoo50_fixed %>% mutate(rangetype = "AOO_50km"),
                           brm_oc_aoo100_fixed %>% mutate(rangetype = "AOO_100km"),
                           brm_oc_ahull6_fixed %>% mutate(rangetype = "Alpha hull")) %>%
  mutate(rangetype = factor(rangetype, levels = c("AOO_50km", "AOO_100km", "Alpha hull")))

# generate figures
ggplot(brm_oc_range_coef) + 
  facet_wrap(~ rangetype, scales = "free_x") + 
  geom_histogram(aes(estimate_slope, fill = sig_slope, colour = NULL), alpha = 1) +
  geom_vline(xintercept = 0, linetype =2, colour = "black") +
  geom_vline(data = brm_oc_range_fixed %>% filter(term == "slope"), aes(xintercept = Estimate), linetype =1) + 
  geom_rect(data = brm_oc_range_fixed %>% filter(term == "slope"), 
            aes(xmin = Q2.5, xmax = Q97.5, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "gray50") + 
  labs(x = "Effect of range size on occupancy change", y = "Number of studies") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) +
  theme_bw() +
  theme(legend.position = "n", 
        axis.title = element_text(size = 12, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        strip.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("neutral" = "light gray", "negative" = "orange", "positive" = "skyblue"))

ggsave("results/Fig.S_sensitivity_slopes_rangetype.png", width = 150, height = 65, units = 'mm')


############
## Figure comparing slopes shown in main text with slopes based on different measures of range size

# combine study-level slopes
brm_oc_range_coef_sens <- brm_oc_aoo10_coef %>%
  dplyr::select(study, database, realm, taxon_new, slope_aoo10 = estimate_slope, Q2.5_aoo10 = Q2.5_slope, Q97.5_aoo10 = Q97.5_slope) %>%
  left_join(brm_oc_aoo50_coef %>%
              dplyr::select(study, slope_aoo50 = estimate_slope, Q2.5_aoo50 = Q2.5_slope, Q97.5_aoo50 = Q97.5_slope)) %>%
  left_join(brm_oc_aoo100_coef %>%
              dplyr::select(study, slope_aoo100 = estimate_slope, Q2.5_aoo100 = Q2.5_slope, Q97.5_aoo100 = Q97.5_slope)) %>%
  left_join(brm_oc_ahull6_coef %>%
              dplyr::select(study, slope_ahull6 = estimate_slope, Q2.5_ahull6 = Q2.5_slope, Q97.5_ahull6 = Q97.5_slope)) 

# combine global slopes
brm_oc_range_fixed_sens <- brm_oc_aoo10_fixed %>%
  dplyr::select(term, slope_aoo10 = Estimate, Q2.5_aoo10 = Q2.5, Q97.5_aoo10 = Q97.5) %>%
  left_join(brm_oc_aoo50_fixed %>%
              dplyr::select(term, slope_aoo50 = Estimate, Q2.5_aoo50 = Q2.5, Q97.5_aoo50 = Q97.5)) %>%
  left_join(brm_oc_aoo100_fixed %>%
              dplyr::select(term, slope_aoo100 = Estimate, Q2.5_aoo100 = Q2.5, Q97.5_aoo100 = Q97.5)) %>%
  left_join(brm_oc_ahull6_fixed %>%
              dplyr::select(term, slope_ahull6 = Estimate, Q2.5_ahull6 = Q2.5, Q97.5_ahull6 = Q97.5)) %>%
  filter(term == "slope")


slope_aoo10_aoo50 <- ggplot(data = brm_oc_range_coef_sens) +
  geom_linerange(aes(x = slope_aoo10, ymin = Q2.5_aoo50, ymax = Q97.5_aoo50), color = "darkgray", size = 0.3) + 
  geom_linerange(aes(y = slope_aoo50, xmin = Q2.5_aoo10, xmax = Q97.5_aoo10), color = "darkgray", size = 0.3) + 
  geom_point(aes(x = slope_aoo10, y = slope_aoo50), shape = 16, size = 0.8, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  geom_point(data = brm_oc_range_fixed_sens, aes(x = slope_aoo10, y = slope_aoo50), 
             shape = 16, size = 0.8, alpha = 1, color = "red") + 
  geom_linerange(data = brm_oc_range_fixed_sens, aes(x = slope_aoo10, ymin = Q2.5_aoo50, ymax = Q97.5_aoo50), size = 0.3, color = "red") + 
  geom_linerange(data = brm_oc_range_fixed_sens, aes(y = slope_aoo50, xmin = Q2.5_aoo10, xmax = Q97.5_aoo10), size = 0.3, color = "red") + 
  geom_abline(intercept = 0, slope = 1, linetype =2, size = 0.3, colour = "blue") + 
  xlim(min(brm_oc_range_coef_sens[,c(6, 9)]), max(brm_oc_range_coef_sens[,c(7, 10)])) + 
  ylim(min(brm_oc_range_coef_sens[,c(6, 9)]), max(brm_oc_range_coef_sens[,c(7, 10)])) + 
  labs(x = "", y = "Slopes based on AOO in 50 km", tag = "A") + 
  theme_classic() +
  #theme_bw() + 
  theme(legend.position = "n", 
        #plot.margin = unit(c(0.02, 0.05, 0, 0.05), "cm"),
        text = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

slope_aoo10_aoo100 <- ggplot(data = brm_oc_range_coef_sens) +
  geom_linerange(aes(x = slope_aoo10, ymin = Q2.5_aoo100, ymax = Q97.5_aoo100), color = "darkgray", size = 0.3) + 
  geom_linerange(aes(y = slope_aoo100, xmin = Q2.5_aoo10, xmax = Q97.5_aoo10), color = "darkgray", size = 0.3) + 
  geom_point(aes(x = slope_aoo10, y = slope_aoo100), shape = 16, size = 0.8, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  geom_point(data = brm_oc_range_fixed_sens, aes(x = slope_aoo10, y = slope_aoo100), 
             shape = 16, size = 0.8, alpha = 1, color = "red") + 
  geom_linerange(data = brm_oc_range_fixed_sens, aes(x = slope_aoo10, ymin = Q2.5_aoo100, ymax = Q97.5_aoo100), size = 0.3, color = "red") + 
  geom_linerange(data = brm_oc_range_fixed_sens, aes(y = slope_aoo100, xmin = Q2.5_aoo10, xmax = Q97.5_aoo10), size = 0.3, color = "red") + 
  geom_abline(intercept = 0, slope = 1, linetype =2, size = 0.3, colour = "blue") + 
  xlim(min(brm_oc_range_coef_sens[,c(6, 12)]), max(brm_oc_range_coef_sens[,c(7, 13)])) + 
  ylim(min(brm_oc_range_coef_sens[,c(6, 12)]), max(brm_oc_range_coef_sens[,c(7, 13)])) + 
  labs(x = "", y = "Slopes based on AOO in 100 km", tag = "B") + 
  theme_classic() +
  #theme_bw() + 
  theme(legend.position = "n", 
        #plot.margin = unit(c(0.02, 0.05, 0, 0.05), "cm"),
        text = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

slope_aoo10_ahull6 <- ggplot(data = brm_oc_range_coef_sens) +
  geom_linerange(aes(x = slope_aoo10, ymin = Q2.5_ahull6, ymax = Q97.5_ahull6), color = "darkgray", size = 0.3) + 
  geom_linerange(aes(y = slope_ahull6, xmin = Q2.5_aoo10, xmax = Q97.5_aoo10), color = "darkgray", size = 0.3) + 
  geom_point(aes(x = slope_aoo10, y = slope_ahull6), shape = 16, size = 0.8, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  geom_point(data = brm_oc_range_fixed_sens, aes(x = slope_aoo10, y = slope_ahull6), 
             shape = 16, size = 0.8, alpha = 1, color = "red") + 
  geom_linerange(data = brm_oc_range_fixed_sens, aes(x = slope_aoo10, ymin = Q2.5_ahull6, ymax = Q97.5_ahull6), size = 0.3, color = "red") + 
  geom_linerange(data = brm_oc_range_fixed_sens, aes(y = slope_ahull6, xmin = Q2.5_aoo10, xmax = Q97.5_aoo10), size = 0.3, color = "red") + 
  geom_abline(intercept = 0, slope = 1, linetype =2, size = 0.3, colour = "blue") + 
  xlim(min(brm_oc_range_coef_sens[,c(6, 12)]), max(brm_oc_range_coef_sens[,c(7, 13)])) + 
  ylim(min(brm_oc_range_coef_sens[,c(6, 12)]), max(brm_oc_range_coef_sens[,c(7, 13)])) + 
  labs(x = "", y = "Slopes based on alpha hull", tag = "C") + 
  theme_classic() +
  #theme_bw() + 
  theme(legend.position = "n", 
        #plot.margin = unit(c(0.02, 0.05, 0, 0.05), "cm"),
        text = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


cowplot::plot_grid(slope_aoo10_aoo50, slope_aoo10_aoo100, slope_aoo10_ahull6, nrow = 1, align = "hv") +
  cowplot::draw_label("Slopes reported in main text", y = 0.05, size = 8)

ggsave("results/Fig.S_sensitivity_slopes_rangesize.png", width = 180, height = 60, units = 'mm', dpi = 600)


############
## Figure comparing slopes shown in main text with 
# slopes based on samples in same locations across year, slopes based on first and last year data,
# based on species with relatively more occurrences in GBIF, based on species with initial occupancy >0 & <1

# combine study-level slopes
brm_oc_aoo10_coef_sens <- brm_oc_aoo10_coef %>%
  dplyr::select(study, database, realm, taxon_new, slope_all = estimate_slope, Q2.5_all = Q2.5_slope, Q97.5_all = Q97.5_slope) %>%
  left_join(brm_oc_aoo10_sloc_coef %>%
              dplyr::select(study, slope_sloc = estimate_slope, Q2.5_sloc = Q2.5_slope, Q97.5_sloc = Q97.5_slope)) %>%
  left_join(brm_oc_aoo10_2yr_coef %>%
              dplyr::select(study, slope_2yr = estimate_slope, Q2.5_2yr = Q2.5_slope, Q97.5_2yr = Q97.5_slope)) %>%
  left_join(brm_oc_aoo10_rgbif_coef %>%
              dplyr::select(study, slope_rgbif = estimate_slope, Q2.5_rgbif = Q2.5_slope, Q97.5_rgbif = Q97.5_slope)) %>% 
  left_join(brm_oc_aoo10_no01_coef %>%
              dplyr::select(study, slope_no01 = estimate_slope, Q2.5_no01 = Q2.5_slope, Q97.5_no01 = Q97.5_slope))

# combine global slopes
brm_oc_aoo10_fixed_sens <- brm_oc_aoo10_fixed %>%
  dplyr::select(term, slope_all = Estimate, Q2.5_all = Q2.5, Q97.5_all = Q97.5) %>%
  left_join(brm_oc_aoo10_sloc_fixed %>%
              dplyr::select(term, slope_sloc = Estimate, Q2.5_sloc = Q2.5, Q97.5_sloc = Q97.5)) %>%
  left_join(brm_oc_aoo10_2yr_fixed %>%
              dplyr::select(term, slope_2yr = Estimate, Q2.5_2yr = Q2.5, Q97.5_2yr = Q97.5)) %>%
  left_join(brm_oc_aoo10_rgbif_fixed %>%
              dplyr::select(term, slope_rgbif = Estimate, Q2.5_rgbif = Q2.5, Q97.5_rgbif = Q97.5)) %>%
  left_join(brm_oc_aoo10_no01_fixed %>%
              dplyr::select(term, slope_no01 = Estimate, Q2.5_no01 = Q2.5, Q97.5_no01 = Q97.5)) %>%
  filter(term == "slope")

slope_all_sloc <- ggplot(data = brm_oc_aoo10_coef_sens %>% filter(!is.na(slope_sloc))) +
  geom_linerange(aes(x = slope_all, ymin = Q2.5_sloc, ymax = Q97.5_sloc), color = "darkgray", size = 0.3) + 
  geom_linerange(aes(y = slope_sloc, xmin = Q2.5_all, xmax = Q97.5_all), color = "darkgray", size = 0.3) + 
  geom_point(aes(x = slope_all, y = slope_sloc), shape = 16, size = 0.8, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  geom_point(data = brm_oc_aoo10_fixed_sens, aes(x = slope_all, y = slope_sloc), 
             shape = 16, size = 0.8, alpha = 1, color = "red") + 
  geom_linerange(data = brm_oc_aoo10_fixed_sens, aes(x = slope_all, ymin = Q2.5_sloc, ymax = Q97.5_sloc), size = 0.3, color = "red") + 
  geom_linerange(data = brm_oc_aoo10_fixed_sens, aes(y = slope_sloc, xmin = Q2.5_all, xmax = Q97.5_all), size = 0.3, color = "red") + 
  geom_abline(intercept = 0, slope = 1, linetype =2, size = 0.3, colour = "blue") + 
  xlim(min(brm_oc_aoo10_coef_sens[,c(6, 9)]), max(brm_oc_aoo10_coef_sens[,c(7, 10)])) + 
  ylim(min(brm_oc_aoo10_coef_sens[,c(6, 9)]), max(brm_oc_aoo10_coef_sens[,c(7, 10)])) + 
  labs(x = "", y = "Slopes based on samples\n in same locations across years", tag = "A") + 
  theme_classic() +
  #theme_bw() + 
  theme(legend.position = "n", 
        plot.margin = unit(c(0.02, 0.05, 0, 0.05), "cm"),
        text = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

slope_all_2yr <- ggplot(data = brm_oc_aoo10_coef_sens) +
  geom_linerange(aes(x = slope_all, ymin = Q2.5_2yr, ymax = Q97.5_2yr), color = "darkgray", size = 0.3) + 
  geom_linerange(aes(y = slope_2yr, xmin = Q2.5_all, xmax = Q97.5_all), color = "darkgray", size = 0.3) + 
  geom_point(aes(x = slope_all, y = slope_2yr), shape = 16, size = 0.8, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  geom_point(data = brm_oc_aoo10_fixed_sens, aes(x = slope_all, y = slope_2yr), 
             shape = 16, size = 0.8, alpha = 1, color = "red") + 
  geom_linerange(data = brm_oc_aoo10_fixed_sens, aes(x = slope_all, ymin = Q2.5_2yr, ymax = Q97.5_2yr), size = 0.3, color = "red") + 
  geom_linerange(data = brm_oc_aoo10_fixed_sens, aes(y = slope_2yr, xmin = Q2.5_all, xmax = Q97.5_all), size = 0.3, color = "red") + 
  geom_abline(intercept = 0, slope = 1, linetype =2, size = 0.3, colour = "blue") + 
  xlim(min(brm_oc_aoo10_coef_sens[,c(6, 12)]), max(brm_oc_aoo10_coef_sens[,c(7, 13)])) + 
  ylim(min(brm_oc_aoo10_coef_sens[,c(6, 12)]), max(brm_oc_aoo10_coef_sens[,c(7, 13)])) + 
  labs(x = "", y = "Slopes based on samples\n in first and last year", tag = "B") + 
  theme_classic() +
  #theme_bw() + 
  theme(legend.position = "n", 
        plot.margin = unit(c(0.02, 0.05, 0, 0.05), "cm"),
        text = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

slope_all_rgbif <- ggplot(data = brm_oc_aoo10_coef_sens %>% filter(!is.na(slope_rgbif))) +
  geom_linerange(aes(x = slope_all, ymin = Q2.5_rgbif, ymax = Q97.5_rgbif), color = "darkgray", size = 0.3) + 
  geom_linerange(aes(y = slope_rgbif, xmin = Q2.5_all, xmax = Q97.5_all), color = "darkgray", size = 0.3) + 
  geom_point(aes(x = slope_all, y = slope_rgbif), shape = 16, size = 0.8, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  geom_point(data = brm_oc_aoo10_fixed_sens, aes(x = slope_all, y = slope_rgbif), 
             shape = 16, size = 0.8, alpha = 1, color = "red") + 
  geom_linerange(data = brm_oc_aoo10_fixed_sens, aes(x = slope_all, ymin = Q2.5_rgbif, ymax = Q97.5_rgbif), size = 0.3, color = "red") + 
  geom_linerange(data = brm_oc_aoo10_fixed_sens, aes(y = slope_rgbif, xmin = Q2.5_all, xmax = Q97.5_all), size = 0.3, color = "red") + 
  geom_abline(intercept = 0, slope = 1, linetype =2, size = 0.3, colour = "blue") + 
  xlim(min(brm_oc_aoo10_coef_sens[,c(6, 15)], na.rm=TRUE), max(brm_oc_aoo10_coef_sens[,c(7, 16)], na.rm=TRUE)) + 
  ylim(min(brm_oc_aoo10_coef_sens[,c(6, 15)], na.rm=TRUE), max(brm_oc_aoo10_coef_sens[,c(7, 16)], na.rm=TRUE)) + 
  labs(x = "", y = "Slopes based on species with\n relatively more occurrences in GBIF", tag = "C") + 
  theme_classic() +
  #theme_bw() + 
  theme(legend.position = "n", 
        plot.margin = unit(c(0.02, 0.05, 0.08, 0.05), "cm"),
        text = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

slope_all_no01  <- ggplot(data = brm_oc_aoo10_coef_sens %>% filter(!is.na(slope_no01))) +
  geom_linerange(aes(x = slope_all, ymin = Q2.5_no01, ymax = Q97.5_no01), color = "darkgray", size = 0.3) + 
  geom_linerange(aes(y = slope_no01, xmin = Q2.5_all, xmax = Q97.5_all), color = "darkgray", size = 0.3) + 
  geom_point(aes(x = slope_all, y = slope_no01), shape = 16, size = 0.8, alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype =2, size = 0.3, colour = "black") +
  geom_vline(xintercept = 0, linetype =2, size = 0.3, colour = "black") + 
  geom_point(data = brm_oc_aoo10_fixed_sens, aes(x = slope_all, y = slope_no01), 
             shape = 16, size = 0.8, alpha = 1, color = "red") + 
  geom_linerange(data = brm_oc_aoo10_fixed_sens, aes(x = slope_all, ymin = Q2.5_no01, ymax = Q97.5_no01), size = 0.3, color = "red") + 
  geom_linerange(data = brm_oc_aoo10_fixed_sens, aes(y = slope_no01, xmin = Q2.5_all, xmax = Q97.5_all), size = 0.3, color = "red") + 
  geom_abline(intercept = 0, slope = 1, linetype =2, size = 0.3, colour = "blue") + 
  xlim(min(brm_oc_aoo10_coef_sens[,c(6, 18)], na.rm=TRUE), max(brm_oc_aoo10_coef_sens[,c(7, 19)], na.rm=TRUE)) + 
  ylim(min(brm_oc_aoo10_coef_sens[,c(6, 18)], na.rm=TRUE), max(brm_oc_aoo10_coef_sens[,c(7, 19)], na.rm=TRUE)) + 
  labs(x = "", y = "Slopes based on species with\n initial occupancy >0 & <1", tag = "D") + 
  theme_classic() +
  #theme_bw() + 
  theme(legend.position = "n", 
        plot.margin = unit(c(0.02, 0.05, 0.08, 0.05), "cm"),
        text = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


cowplot::plot_grid(slope_all_sloc, slope_all_2yr, slope_all_rgbif, slope_all_no01, 
                   nrow = 2, rel_heights = c(1, 1.02)) +
  cowplot::draw_label("Slopes reported in main text", y = 0.015, size = 8)

ggsave("results/Fig.S_sensitivity_slopes_compare.png", width = 120, height = 110, units = 'mm', dpi = 600)

