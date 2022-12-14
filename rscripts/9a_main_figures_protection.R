## draw Fig. 3, 
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
load("results/fixef_prediction_protection.RDATA")


## fig. 3: prediction of occupancy change across range size for totally unprotected and protected regions
# the results came from the model using proportion of sites in early-established PA
oc_aoo10_psampbflate_fitted <- oc_aoo10_psampbflate_fitted %>%
  mutate(protection = factor(psamp_inPA_bflate, levels = c(0, 1), labels = c("Unprotected", "Protected")))

beta_psampbflate_text <- fixef_aoo10_psampbflate %>% 
  filter(term == "slope") %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         label = paste0('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')')) # beta(range %*% protection)

protection_col <- c("Unprotected" = "purple", "Protected" = "#1F78B4")

ggplot(oc_aoo10_psampbflate_fitted) +
  facet_wrap(~ realm, scales = "fixed") + 
  geom_ribbon(aes(x = aoo10, ymin = Q2.5, ymax = Q97.5, fill = protection), alpha = 0.3) +
  geom_line(aes(x = aoo10, y = Estimate, color = protection), size =1) + 
  geom_hline(yintercept = 0, lty = 2, size = 0.3) + 
  geom_text(data = beta_psampbflate_text,  aes(x= 1, y = Inf, label = label), 
            hjust = 0, vjust = 2.0,  parse = TRUE, size = 2.46) + #2.11
  scale_x_log10() + 
  labs(x = "Range size (number of 10-km grid-cells)", y = "Occupancy change", colour = "", fill = "") + 
  scale_colour_manual(name = NULL, values = protection_col) + 
  scale_fill_manual(name = NULL, values = protection_col) + 
  theme_bw() +
  theme(legend.position = c(0.20, 0.2),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(5, 'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8))

ggsave("results/Fig3.pdf", width = 170, height = 65, units = 'mm')
ggsave("results/Fig3.png", width = 170, height = 65, units = 'mm', dpi = 600)



## supplement figure from the model using proportion of sites in all PA
oc_aoo10_psampwdpa_fitted <- oc_aoo10_psampwdpa_fitted %>%
  mutate(protection = factor(psamp_inPA_bflate, levels = c(0, 1), labels = c("Unprotected", "Protected")))

beta_psampwdpa_text <- fixef_aoo10_psampwdpa %>% 
  filter(term == "slope") %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         label = paste0('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

ggplot(oc_aoo10_psampwdpa_fitted) +
  facet_wrap(~ realm, scales = "fixed") + 
  geom_ribbon(aes(x = aoo10, ymin = Q2.5, ymax = Q97.5, fill = protection), alpha = 0.3) +
  geom_line(aes(x = aoo10, y = Estimate, color = protection), size =1) + 
  geom_hline(yintercept = 0, lty = 2, size = 0.3) + 
  geom_text(data = beta_psampwdpa_text,  aes(x= 1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = TRUE, size = 2.46) +
  scale_x_log10() + 
  labs(x = "Range size (number of 10-km grid-cells)", y = "Occupancy change", colour = "", fill = "") + 
  scale_colour_manual(name = NULL, values = protection_col) + 
  scale_fill_manual(name = NULL, values = protection_col) + 
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(5, 'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8))

ggsave("results/Fig.S_OccupancyChange_aoo10_protection-psampwdpa.png", width = 170, height = 65, units = 'mm', dpi = 600)



## supplement figure from the model using proportion of regional area in early-established PA
oc_aoo10_pareabflate_fitted <- oc_aoo10_pareabflate_fitted %>%
  mutate(protection = factor(psamp_inPA_bflate, levels = c(0, 1), labels = c("Unprotected", "Protected")))

beta_pareabflate_text <- fixef_aoo10_pareabflate %>% 
  filter(term == "slope") %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         label = paste0('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

ggplot(oc_aoo10_pareabflate_fitted) +
  facet_wrap(~ realm, scales = "fixed") + 
  geom_ribbon(aes(x = aoo10, ymin = Q2.5, ymax = Q97.5, fill = protection), alpha = 0.3) +
  geom_line(aes(x = aoo10, y = Estimate, color = protection), size =1) + 
  geom_hline(yintercept = 0, lty = 2, size = 0.3) + 
  geom_text(data = beta_pareabflate_text,  aes(x= 1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = TRUE, size = 2.46) +
  scale_x_log10() + 
  labs(x = "Range size (number of 10-km grid-cells)", y = "Occupancy change", colour = "", fill = "") + 
  scale_colour_manual(name = NULL, values = protection_col) + 
  scale_fill_manual(name = NULL, values = protection_col) + 
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(5, 'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8))

ggsave("results/Fig.S_OccupancyChange_aoo10_protection-pareabflate.png", width = 170, height = 65, units = 'mm', dpi = 600)



## supplement figure from the model using proportion of regional area in all PA
oc_aoo10_pareawdpa_fitted <- oc_aoo10_pareawdpa_fitted %>%
  mutate(protection = factor(psamp_inPA_bflate, levels = c(0, 1), labels = c("Unprotected", "Protected")))

beta_pareawdpa_text <- fixef_aoo10_pareawdpa %>% 
  filter(term == "slope") %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         label = paste0('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

ggplot(oc_aoo10_pareawdpa_fitted) +
  facet_wrap(~ realm, scales = "fixed") + 
  geom_ribbon(aes(x = aoo10, ymin = Q2.5, ymax = Q97.5, fill = protection), alpha = 0.3) +
  geom_line(aes(x = aoo10, y = Estimate, color = protection), size =1) + 
  geom_hline(yintercept = 0, lty = 2, size = 0.3) + 
  geom_text(data = beta_pareawdpa_text,  aes(x= 1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = TRUE, size = 2.46) +
  scale_x_log10() + 
  labs(x = "Range size (number of 10-km grid-cells)", y = "Occupancy change", colour = "", fill = "") + 
  scale_colour_manual(name = NULL, values = protection_col) + 
  scale_fill_manual(name = NULL, values = protection_col) + 
  theme_bw() +
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(5, 'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8))

ggsave("results/Fig.S_OccupancyChange_aoo10_protection-pareawdpa.png", width = 170, height = 65, units = 'mm', dpi = 600)



#################
## compare proportion of sites in any protected areas and those established before the later period

# colors for realms
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")

dat_meta <- dat_meta %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## relationship between proportion of sites in protected areas estimated using all PA and those established before the later period
compare_PA_wdpa_bflate <- ggplot(dat_meta %>% filter(!is.na(psamp_inPA_wdpa))) + 
  facet_wrap(~ realm, scales = "fixed") + 
  geom_jitter(aes(y = psamp_inPA_wdpa, x = psamp_inPA_bflate, color = realm), 
              width = 0.01, height = 0.01,
             shape = 19, size = 1, alpha = 0.8) +
  geom_abline(intercept = 0, slope =1, lty =2, size = 0.3) + 
  labs(y = "Proportion of sites in protected areas \nestablished at any time points", 
       x = "Proportion of sites in protected areas established before later period",
       tag = "a") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  theme_bw()+ 
  theme(legend.position = "n", 
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## relationship between proportion of sites or region area in protected areas estimated before the later period
compare_PA_psamp_parea <- ggplot(dat_meta %>% filter(!is.na(psamp_inPA_bflate))) + 
  facet_wrap(~ realm, scales = "fixed") + 
  geom_jitter(aes(y = parea_inPA_bflate, x = psamp_inPA_bflate, color = realm), 
              width = 0.01, height = 0.01,
              shape = 19, size = 1, alpha = 0.8) +
  geom_abline(intercept = 0, slope =1, lty =2, size = 0.3) + 
  labs(y = "Proportion of regional area \nin protected areas", 
       x = "Proportion of sites in protected areas established before later period",
       tag = "b") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  theme_bw()+ 
  theme(legend.position = "n", 
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## histogram of proportion of sites in protected areas estimated before the later period
histogram_PA_bflate <- ggplot(dat_meta %>% filter(!is.na(psamp_inPA_wdpa))) + 
  facet_wrap(~ realm, scales = "free_y") + 
  geom_histogram(aes(psamp_inPA_bflate, fill = realm)) +
  labs(x = "Proportion of sites in protected areas established before later period",
       y = "Number of studies", 
       tag = "c") +
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()+
  theme(legend.position = "n", 
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        plot.tag = element_text(size = 8, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

cowplot::plot_grid(compare_PA_wdpa_bflate, compare_PA_psamp_parea, histogram_PA_bflate, 
                   nrow = 3, align = "v")
ggsave("results/Fig.S_comopare_protection_measures.png", width = 140, height = 160, unit = "mm")

