###############################################
## Calculate the number and percentage of of sites in protected areas for each study
# three values of protection status were used using protected areas that established any time points or those that 
# established before the start or late periods of assemblage sampling

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","ggplot2", "wdpar", "sf", "sp", "rgdal", "wdpar")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)



###########################
# prepared protected areas files (combine points and polygons and simplify) 
# the files were download from the World Database on Protected Areas (January 2022; https://www.protectedplanet.net/en)

## protected areas file 0
wdpa0_pts <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_0", layer = "WDPA_Jan2022_Public_shp-points")
wdpa0_pys <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_0", layer = "WDPA_Jan2022_Public_shp-polygons")

wdpa0_pts <- wdpa0_pts %>% 
  filter(REP_AREA >0 & (STATUS %in% c("Designated", "Inscribed", "Established")) & DESIG_ENG != "UNESCO-MAB Biosphere Reserve") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>% 
  st_transform("+proj=moll") %>% 
  st_make_valid()
wdpa0_pts <- wdpa0_pts %>% 
  st_buffer(dist = 1000*sqrt(wdpa0_pts$REP_AREA/pi)) %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 0.001) %>% 
  st_make_valid() %>%
  st_transform(crs = 4326) %>%
  mutate(STATUS_YR = ifelse(STATUS_YR==0, NA, STATUS_YR))

wdpa0_pys <- wdpa0_pys %>%
  filter((STATUS %in% c("Designated", "Inscribed", "Established")) & DESIG_ENG != "UNESCO-MAB Biosphere Reserve") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>% 
  st_simplify(preserveTopology = FALSE, dTolerance = 0.001) %>% 
  st_make_valid() %>%
  mutate(STATUS_YR = ifelse(STATUS_YR==0, NA, STATUS_YR))

save(wdpa0_pts, wdpa0_pys, file = "data/WDPA_Jan2022/wdpa0.RDATA")


## protected areas file 1
wdpa1_pts <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_1", layer = "WDPA_Jan2022_Public_shp-points")
wdpa1_pys <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_1", layer = "WDPA_Jan2022_Public_shp-polygons")

wdpa1_pts <- wdpa1_pts %>% 
  filter(REP_AREA >0 & (STATUS %in% c("Designated", "Inscribed", "Established")) & DESIG_ENG != "UNESCO-MAB Biosphere Reserve") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>% 
  st_transform("+proj=moll") %>% 
  st_make_valid()
wdpa1_pts <- wdpa1_pts %>% 
  st_buffer(dist = 1000*sqrt(wdpa1_pts$REP_AREA/pi)) %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 0.001) %>% 
  st_make_valid() %>%
  st_transform(crs = 4326) %>%
  mutate(STATUS_YR = ifelse(STATUS_YR==0, NA, STATUS_YR))

wdpa1_pys <- wdpa1_pys %>%
  filter((STATUS %in% c("Designated", "Inscribed", "Established")) & DESIG_ENG != "UNESCO-MAB Biosphere Reserve") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>% 
  st_simplify(preserveTopology = FALSE, dTolerance = 0.001) %>% 
  st_make_valid() %>%
  mutate(STATUS_YR = ifelse(STATUS_YR==0, NA, STATUS_YR))

save(wdpa1_pts, wdpa1_pys, file = "data/WDPA_Jan2022/wdpa1.RDATA")


## protected areas file 2
wdpa2_pts <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_2", layer = "WDPA_Jan2022_Public_shp-points")
wdpa2_pys <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_2", layer = "WDPA_Jan2022_Public_shp-polygons")

wdpa2_pts <- wdpa2_pts %>% 
  filter(REP_AREA >0 & (STATUS %in% c("Designated", "Inscribed", "Established")) & DESIG_ENG != "UNESCO-MAB Biosphere Reserve") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>% 
  st_transform("+proj=moll") %>% 
  st_make_valid()
wdpa2_pts <- wdpa2_pts %>% 
  st_buffer(dist = 1000*sqrt(wdpa2_pts$REP_AREA/pi)) %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 0.001) %>% 
  st_make_valid() %>%
  st_transform(crs = 4326) %>%
  mutate(STATUS_YR = ifelse(STATUS_YR==0, NA, STATUS_YR))

wdpa2_pys <- wdpa2_pys %>%
  filter((STATUS %in% c("Designated", "Inscribed", "Established")) & DESIG_ENG != "UNESCO-MAB Biosphere Reserve") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>% 
  st_simplify(preserveTopology = FALSE, dTolerance = 0.001) %>% 
  st_make_valid() %>%
  mutate(STATUS_YR = ifelse(STATUS_YR==0, NA, STATUS_YR))

save(wdpa2_pts, wdpa2_pys, file = "data/WDPA_Jan2022/wdpa2.RDATA")


# combine protected areas from different subsets for those originally repressed as points and polygons respectively
#  originally repressed as points
wdpa_pts <- bind_rows(wdpa0_pts, wdpa1_pts, wdpa2_pts) %>% 
  st_make_valid()

#  originally repressed as polygons
wdpa_pys <- bind_rows(wdpa0_pys, wdpa1_pys, wdpa2_pys) %>% 
  st_make_valid()

save(wdpa_pts, wdpa_pys, file = "data/WDPA_Jan2022/wdpa_combined.RDATA")


str(wdpa_pts)
plot(wdpa_pts[,1],axes=TRUE)
table(wdpa_pts$IUCN_CAT)
table(wdpa_pys$IUCN_CAT)
summary(wdpa_pys$REP_AREA)

table(wdpa_pys$STATUS_YR == 0)
table(wdpa_pys$MARINE)
summary(wdpa_pys$STATUS_YR[wdpa_pys$STATUS_YR >0])
summary(wdpa_pys$STATUS_YR[wdpa_pys$MARINE==1 & wdpa_pys$STATUS_YR >0])



################################
# Calculate the number and percentage of of sites in protected areas

## protected area and coordinates of samples
load("data/WDPA_Jan2022/wdpa_combined.RDATA")
load("data/Combined_assemblages_20220208.RDATA")
load("intermediate_results/occuapncy.RDATA")

# the start and end sampling and the year before the later period 
samp_year <- years_period %>%
  group_by(study) %>%
  summarise(middle_year = (max(year) + min(year))/2,
         start_year = min(year),
         end_year = max(year),
         start_late = min(year[period_combined == "late"]))

# the sites within regions
dat_samp <- dat %>% 
  dplyr::select(study, sample, latitude, longitude) %>% 
  distinct() %>%
  left_join(samp_year)

# transfer to spatial points
dat_pts <- dat_samp %>%  
  st_as_sf(coords = c('longitude', 'latitude'))  %>% 
  st_set_crs(4326)

# combine protected areas from points and polygons
wdpa <- bind_rows(wdpa_pts, wdpa_pys)


# over community samples with protected areas
samp_isinpa <- as_tibble(st_join(dat_pts, wdpa, join = st_intersects))

samp_inpa <- samp_isinpa %>% 
  dplyr::select(study, middle_year, start_year, end_year, start_late, sample, WDPAID, IUCN_CAT, STATUS_YR) %>% 
  filter(!is.na(WDPAID))

# number of sites in protected area
nsamp_inPA_bfstart <- samp_inpa %>% 
  filter(STATUS_YR <= start_year) %>%
  group_by(study) %>%
  summarise(nsamp_inPA_bfstart = n_distinct(sample))

nsamp_inPA_bflate <- samp_inpa %>% 
  filter(STATUS_YR < start_late) %>%
  group_by(study) %>%
  summarise(nsamp_inPA_bflate = n_distinct(sample))

nsamp_inPA_wdpa <- samp_inpa %>% 
  group_by(study) %>%
  summarise(nsamp_inPA_wdpa = n_distinct(sample))

nsamp_inPA <- nsamp_inPA_bfstart %>% 
  full_join(nsamp_inPA_bflate) %>%     
  full_join(nsamp_inPA_wdpa)

# percentage of samples in protected areas
dat_psamp_inPA <- dat_samp %>% 
  group_by(study) %>%
  # total number of sites and coordinates
  summarise(nsamp_total = n_distinct(sample),
            ncoord_total = n_distinct(latitude, longitude)) %>%
  ungroup() %>%
  # some studies only have one central coordinates, andt thus the values of psamp_inPA are possibly not accurate
  mutate(coordinate_local = ifelse(ncoord_total == 1, FALSE, TRUE)) %>% 
  left_join(nsamp_inPA) %>% 
  mutate(nsamp_inPA_bfstart = ifelse(is.na(nsamp_inPA_bfstart), 0, nsamp_inPA_bfstart),
         nsamp_inPA_bflate = ifelse(is.na(nsamp_inPA_bflate), 0, nsamp_inPA_bflate),
         nsamp_inPA_wdpa = ifelse(is.na(nsamp_inPA_wdpa), 0, nsamp_inPA_wdpa),
         psamp_inPA_bfstart = nsamp_inPA_bfstart/nsamp_total,
         psamp_inPA_bflate = nsamp_inPA_bflate/nsamp_total,
         psamp_inPA_wdpa = nsamp_inPA_wdpa/nsamp_total) %>%
  left_join(samp_year %>% dplyr::select(study, start_year, start_late)) %>%
  # add study-names and regional extent
  left_join(dat_meta %>% dplyr::select(study, study_name, realm, extent_km2))


ggplot(dat_psamp_inPA) +
  facet_wrap(~ realm, scales = "fixed") + 
  geom_point(aes(psamp_inPA_bflate, psamp_inPA_bfstart))


# for studies with only central coordinators
# output to check manually based on descriptions in original papers. 
# If the regions are largely protected, protection level was set as 1; If are not protected, set as 0. 
# For studies that can't determined based on original descriptions, 
# determine protection level based on central coordinates for regions with area <= 10 km2; and remove large regions
dat_psamp_inPA_manual <- dat_psamp_inPA %>%
  filter(!coordinate_local)
write.csv(dat_psamp_inPA_manual, file = "data/Assemblages_withCentralCoordinates_protection_manual.csv")

dat_psamp_inPA_manual <- read.csv( "data/Assemblages_withCentralCoordinates_protection_manual_filled.csv", row.names = 1)

# update the protection level using the manuuly corrected ones
id <- match(dat_psamp_inPA_manual$study, dat_psamp_inPA$study, nomatch = 0)
table(id == 0)
dat_psamp_inPA[id, "psamp_inPA_bfstart"] <- dat_psamp_inPA_manual$psamp_inPA_bfstart_new
dat_psamp_inPA[id, "psamp_inPA_bflate"] <- dat_psamp_inPA_manual$psamp_inPA_bflate_new
dat_psamp_inPA[id, "psamp_inPA_wdpa"] <- dat_psamp_inPA_manual$psamp_inPA_wdpa_new
dat_psamp_inPA[id, "comment"] <- dat_psamp_inPA_manual$comment


# save results
save(dat_psamp_inPA, file = "data/Assemblages_protection.RDATA")
