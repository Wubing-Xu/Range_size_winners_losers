## draw all other supplementary figures 

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
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) 

brm_oc_aoo10_coef <- brm_oc_aoo10_coef %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#1B9E77", "Marine" = "#7570B3")


###########################
## plot the proportion of species in different dynamics
# number of colonization, extinction species and species that increase or decrease occupancy
species_dynamic <- oc_period %>% 
  group_by(study, realm, studyID) %>%
  summarise(total = n(), 
            Gained = sum(dynamic == "colonization"),
            Lost = sum(dynamic == "extinction"),
            Persisted = sum(dynamic == "persistent"),
            Persisted_increased = sum(dynamic == "persistent" & occup_change > 0), 
            Persisted_decreased = sum(dynamic == "persistent" & occup_change < 0), 
            Persisted_stable = sum(dynamic == "persistent" & occup_change == 0)) %>%
  ungroup()

#get a summary of proportion of species in different dynamics
p_dynamic <- species_dynamic %>% 
  mutate(gained1 = Gained/total,
         lost1 = Lost/total,
         persisted1 = Persisted/total,
         persisted_increased1 = Persisted_increased/total,
         persisted_decreased1 = Persisted_decreased/total,
         persisted_stable1 = Persisted_stable/total,
         p_lost_gained = gained1  + lost1)

summary(p_dynamic)

species_dynamic_longer <- species_dynamic %>%
  mutate(p_persisted = Persisted/total) %>%
  arrange(p_persisted) %>% 
  mutate(study_id = row_number()) %>%
  dplyr::select(realm, study_id, Gained, Lost, Persisted_increased:Persisted_stable) %>%
  pivot_longer(cols = Gained:Persisted_stable, names_to = "dynamic", values_to = "number") %>%
  mutate(dynamic = factor(dynamic, levels = c("Gained", "Lost", "Persisted_increased", "Persisted_decreased", "Persisted_stable")))

# generate figures
plot_dynamic <- ggplot(species_dynamic_longer, aes( x = study_id, y = number, fill = dynamic)) + 
  geom_bar(position = "fill", stat = "identity",  width = 1) +
  labs(x = "Study rank", y = "Proportion of species") + 
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_continuous(expand = c(0.01,0.01), breaks = c(1, 40, 80, 120, 160, 203), labels = c(1, 40, 80, 120, 160, 204)) +
  theme_bw() +
  theme(legend.position = "top", 
        legend.key.size = unit(3, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(2, units = 'mm'),
        plot.margin = unit(c(2, 5, 2, 5), units = "mm"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +  
  scale_fill_manual(values = c("Gained" = "#AA4499", "Lost" = "#332288", 
                               "Persisted_increased" = "#0BFB27", "Persisted_decreased" = "#47C742", "Persisted_stable" = "#117733"))
# scale_fill_manual(values = c("Gained" = "#F8766D", "Lost" = "#00A5FF", 
#                            "Persisted_increased" = "green1", "Persisted_decreased" = "green3", "Persisted_stable" = "green4"))
# scale_fill_manual(values = c("Gained" = "#7570B3", "Lost" = "#D95F02", 
#                            "Persisted_increased" = "#1B9E77", "Persisted_decreased" = "#A1D99B", "Persisted_stable" = "#E5F5E0"))

ggsave(plot_dynamic, file = "results/Fig.S_species_dynamic.png", 
       width = 160, height = 80, units = 'mm')



####################
## Figure showing correlations among measures of range size
# four measures are used: aoo10, aoo50, aoo100, ahull6
library(GGally)

# get rang size of species used in analyses: 39,764 species
rangesize <- oc_period %>% 
  distinct(species, specieskey, taxon_new, taxon_final, aoo10, aoo50, aoo100, ahull6) %>% 
  mutate(aoo10 = log10(aoo10),
         aoo50 = log10(aoo50),
         aoo100 = log10(aoo100),
         ahull6 = log10(ahull6))

# generate figure
corplot_range <- ggpairs(rangesize, columns = c('aoo10', 'aoo50', 'aoo100', 'ahull6'),
        upper = list(continuous = wrap("cor")),
        lower = list(continuous = wrap("points", size = 0.5, alpha = 0.1, shape = 1)),
        columnLabels = c("AOO_10km", "AOO_50km", "AOO_100km", "Alpha hull")) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8),
        strip.text = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(corplot_range, file = "results/Fig.S_correlation_rangesize.png", width = 160, height = 140, units = 'mm')
ggsave(corplot_range, file = "results/Fig.S_correlation_rangesize.pdf", width = 160, height = 140, units = 'mm')


######################
# frequency distribution of range size (aoo-10km) for each taxon

rangesize_raw <- oc_period %>% 
  distinct(species, specieskey, taxon_new, taxon_final, aoo10, aoo50, aoo100, ahull6) 

ggplot(rangesize_raw) +
  facet_wrap(~ taxon_new, scale = "free_y") +
  geom_histogram(aes(aoo10, fill = taxon_new)) + 
  scale_x_log10() +
  labs(x = "Range size (number of 10-km grid cells)", y = "Number of species") + 
  theme_bw() +
  theme(legend.position = "no", 
        #strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        strip.text = element_text(size = 9, face = "plain"))

ggsave(file = "results/Fig.S_rangesize_distribution_aoo10.png", width = 159, height = 106, units = 'mm')



####################
## figure comparing number of occurrences from GBIF and community data

nocc <- oc_period %>% 
  distinct(species, specieskey, nocc_community, nocc_gbif, ratio_nocc_gbif) 

table(nocc$ratio_nocc_gbif<1)
length(which(nocc$ratio_nocc_gbif<1))/nrow(nocc)  # 4.9%
table(nocc$ratio_nocc_gbif<5) 
length(which(nocc$ratio_nocc_gbif<5))/nrow(nocc) # 18.3%

# compare number of occurrences from GBIF and community data
plot_nocc <- ggplot(nocc) + 
  geom_point(aes(nocc_community, nocc_gbif), size = 0.3, shape = 16, alpha = 0.3) + 
  geom_abline(intercept = 0, slope =1, lty = 2, col = "blue", size = 0.5) +
  scale_x_log10(limits = c(min(nocc[,c(3, 4)]), max(nocc[,c(3, 4)]))) + 
  scale_y_log10(limits = c(min(nocc[,c(3, 4)]), max(nocc[,c(3, 4)]))) +
  labs(x = "Number of occurrences from community data",
       y = "Number of occurrences from GBIF") +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# histogram of ratio of number of occurrnces from GBIF by those from community data
plot_rgbif <- ggplot(data = nocc %>% filter(ratio_nocc_gbif>0)) +
  geom_histogram(aes(ratio_nocc_gbif)) + 
  geom_vline(xintercept = 1, lty =2, size = 0.5) +
  scale_x_log10(breaks = c(0.1,  10,  1000,  100000)) +
  labs(x = "Ratio of occurrences (GBIF/community data)", y = "Number of species") +
  theme_classic() + 
  theme(axis.title = element_text(size = 8, face = "plain"),
        axis.text = element_text(size = 6, face = "plain"))


cowplot::plot_grid(plot_nocc, plot_rgbif, nrow = 1, align = "hv",
                   labels = c("A", "B"), label_size = 8)
ggsave("results/Fig.S_nocc_compare.png", width = 140, height = 70, units = 'mm', dpi = 600)



#######################
# histogram of proportion of species having estimates of range size, and its relation with study-level slopes 

# histogram of proportion of species having estimates of range size
hist_psprich_range <- ggplot(data = brm_oc_aoo10_coef) +
  geom_histogram(aes(Psprich_range, fill = realm)) + 
  labs(x = "Proportion of species", y = "Number of studies", tag = "A") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = c(0.28, 0.85),
        legend.text = element_text(size = 8),
        legend.key.size = unit(4, "mm"),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

# fit a model testing relation between study-level slopes and proportion of species having estimates of range size
brm_slope_PspRange <- brm(bf(estimate_slope| se(se_slope) ~ 0 + realm + Psprich_range:realm + (1|study)),
                                data = brm_oc_aoo10_coef,
                                control = list(adapt_delta = 0.9, max_treedepth = 10),
                                cores = 4, chains = 4, iter = 8000, thin = 4,
                                file = "models/brm_output/brm_slope_PspRange")
brm_slope_PspRange
pp_check(brm_slope_PspRange)
fixef(brm_slope_PspRange, probs = c(0.025, 0.1, 0.9, 0.975))

# relation with study-level slopes 
plot_slope_psprich_range <- ggplot(data = brm_oc_aoo10_coef, aes(x = Psprich_range, y = estimate_slope, colour = realm)) +
  geom_point(shape = 19, size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", aes(fill = realm)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3, colour = "black") + 
  labs(x = "Proportion of species", y = "Effect of range size", tag = "B") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() + 
  theme(legend.position = "no",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        plot.tag = element_text(size = 9, face = "bold"))

cowplot::plot_grid(hist_psprich_range, plot_slope_psprich_range, 
                   nrow = 1, align = "hv")
ggsave("results/Fig.S_psprich_range-1.png", width = 140, height = 70, units = 'mm', dpi = 600)


###############################
## relationship between occupancy in the early and late period

occupancy_first_last <- oc_period %>% 
  dplyr::select(study, realm, species, specieskey, occup_first, occup_last) 

ggplot(occupancy_first_last) +
  geom_bin2d(aes(occup_first, occup_last, fill = log10(..count..))) + 
  scale_fill_continuous(type = "viridis") + 
  labs(x = "Occupancy in the early period",
       y = "Occupany in the late period",
       fill = bquote(Log[10]*'(number of populations)')) + 
  theme_bw() +
  theme(legend.position= "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, units = 'mm'),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("results/Fig.S_occupancy_early_later.png", width = 80, height = 85, units = 'mm')



##########################
##  fixed effects of range size on occupancy change across taxa and regions

library(ggstance)
brm_oc_aoo10_realm_taxa_fixed <- brm_oc_aoo10_realm_taxa_fixed %>% 
  mutate(taxon_new =ifelse(taxon_new == "Amphibians and reptiles", "Amphibians/reptiles", as.character(taxon_new))) %>% 
  filter(term == "slope")

brm_oc_aoo10_realm_region_fixed <- brm_oc_aoo10_realm_region_fixed %>% 
  filter(term == "slope")

# limits of fixed effects across groups
xlimit_low <- min(brm_oc_aoo10_realm_taxa_fixed[,  "Q2.5"], brm_oc_aoo10_realm_region_fixed[,  "Q2.5"])
xlimit_high <- max(brm_oc_aoo10_realm_taxa_fixed[,  "Q97.5"], brm_oc_aoo10_realm_region_fixed[,  "Q97.5"])

# fixed effects across taxa
plot_taxon_fixef <- ggplot(brm_oc_aoo10_realm_taxa_fixed %>% filter(term == "slope")) +
  geom_linerangeh(aes(xmin=Q2.5, xmax=Q97.5, y = fct_rev(taxon_new), colour = realm), 
                 position = position_dodge(0.7), size =0.8) +
  geom_linerangeh(aes(xmin=Q10, xmax=Q90, y = fct_rev(taxon_new), colour = realm), 
                  position = position_dodge(0.7), size = 1.5) + 
  geom_point(aes(y = taxon_new, x = Estimate, colour = realm),
             position = position_dodge(0.7), size = 2) + 
  geom_text(aes(y = taxon_new, x = 1.4, label = paste('n == ', nstudy), group = realm), 
            position = position_dodge(0.7), parse = T, size = 2.1, hjust = 0) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  labs(x = "Effect of range size", ylab = NULL, colour = NULL) + 
  scale_x_continuous(limits = c(xlimit_low, 1.3*xlimit_high)) +
  scale_colour_manual(name = NULL, values = realm_col) +
  theme_classic() + 
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7),
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, units = 'mm'),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8)) 

# fixed effects across regions
plot_region_fixef <- ggplot(brm_oc_aoo10_realm_region_fixed %>% filter(term == "slope")) +
  geom_linerangeh(aes(xmin=Q2.5, xmax=Q97.5, y = fct_rev(region), colour = realm), 
                  position = position_dodge(0.7), size =0.8) +
  geom_linerangeh(aes(xmin=Q10, xmax=Q90, y = fct_rev(region), colour = realm), 
                  position = position_dodge(0.7), size = 1.5) + 
  geom_point(aes(y = region, x = Estimate, colour = realm),
             position = position_dodge(0.7), size = 2) + 
    geom_text(aes(y = region, x = 1.4, label = paste('n == ', nstudy), group = realm), 
            position = position_dodge(0.7), parse = T, size = 2.1, hjust = 0) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) + 
  scale_x_continuous(limits = c(xlimit_low, 1.3*xlimit_high)) +
  labs(x = "Effect of range size", ylab = NULL, colour = NULL) + 
  scale_colour_manual(name = NULL, values = realm_col) +
  theme_classic() + 
  theme(legend.position = "n",
        axis.title.y = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) 

cowplot::plot_grid(plot_taxon_fixef, plot_region_fixef, 
                   nrow = 2, rel_heights = c(1, 1.3), align = "v",
                   labels = c("A", "B"), label_size = 9)
ggsave("results/Fig.S_fixef_taxon_region.png", width = 100, height = 140, units = 'mm', dpi = 600)



##########################
## figure showing range-size occupancy change relations across taxon and region groups

# set factor levels
brm_oc_aoo10_realm_taxa_fitted <- brm_oc_aoo10_realm_taxa_fitted %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         taxon_new =ifelse(taxon_new == "Amphibians and reptiles", "Amphibians/reptiles", as.character(taxon_new)))

brm_oc_aoo10_realm_taxa_line <- brm_oc_aoo10_realm_taxa_line %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         taxon_new =ifelse(taxon_new == "Amphibians and reptiles", "Amphibians/reptiles", as.character(taxon_new)))

brm_oc_aoo10_realm_region_fitted <- brm_oc_aoo10_realm_region_fitted %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         region = factor(region, levels = c("Africa", "Asia", "Australia", "Europe", "North America", "South America",
                                          "Atlantic", "Pacific", "Arctic Ocean", "Indian Ocean", "Southern Ocean"),
                       labels = c("Africa", "Asia", "Australia", "Europe", "North America", "South America",
                                  "Atlantic Ocean", "Pacific Ocean", "Other Oceans", "Other Oceans", "Other Oceans"))) 

brm_oc_aoo10_realm_region_line <- brm_oc_aoo10_realm_region_line %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

# relationship across taxa
plot_oc_aoo_taxa <- ggplot() + 
  facet_wrap(~ realm, scale = "free") +
  geom_ribbon(data = brm_oc_aoo10_realm_taxa_fitted,
              aes(x = aoo10, ymin = occup_change_logit_Q2.5, ymax = occup_change_logit_Q97.5, fill = taxon_new), 
              size =1, alpha = 0.5) + 
  geom_segment(data = brm_oc_aoo10_realm_taxa_line, 
               aes(x = xmin, xend = xmax, 
                   y = intercept + slope * cl.xmin, 
                   yend = intercept + slope * cl.xmax,
                   colour = taxon_new), 
               size = 1, alpha = 1) + 
  geom_hline(yintercept = 0, lty = 2, size = 0.3) + 
  scale_y_continuous(limits = c(-4, 2), breaks = c( -4, -2, 0, 2)) + 
  scale_x_log10() + 
  labs(x = NULL, y = "Occupancy change", colour = "", fill = "") +
  theme_classic() +
  theme(legend.position = "right", 
        legend.key.size = unit(5, "mm"),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 9),
        strip.text = element_text(size = 9),
        strip.background = element_blank())

# relationship across regions
plot_oc_aoo_region <- ggplot() + 
  facet_wrap(~ realm, scale = "free") +
  geom_ribbon(data = brm_oc_aoo10_realm_region_fitted,
              aes(x = aoo10, ymin = occup_change_logit_Q2.5, ymax = occup_change_logit_Q97.5, fill = region), 
              size =1, alpha = 0.3) + 
  geom_segment(data = brm_oc_aoo10_realm_region_line, 
               aes(x = xmin, xend = xmax, 
                   y = intercept + slope * cl.xmin, 
                   yend = intercept + slope * cl.xmax,
                   colour = region), 
               size = 1, alpha = 1) + 
  geom_hline(yintercept = 0, lty = 2, size = 0.3) + 
  scale_y_continuous(limits = c(-4, 2), breaks = c( -4, -2, 0, 2)) + 
  scale_x_log10() + 
  labs(x = "Range size (number of 10-km grid-cells)", y = "Occupancy change", colour = NULL, fill = NULL) +
  theme_classic() +
  theme(legend.position = "right", 
        legend.key.size = unit(5, "mm"),
        legend.margin = margin(0, 0, 0, -10),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 9),
        strip.text = element_blank())

cowplot::plot_grid(plot_oc_aoo_taxa, plot_oc_aoo_region,
                   nrow = 2, rel_heights = c(1.04, 1), align = "v",
                   labels = c('A', "B"), label_size = 9)

ggsave(file = "results/Fig.S_oc_rangesize_taxa_region.png", 
       width = 180, height = 120, units = 'mm', dpi = 300)
