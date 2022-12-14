###############################################
## Calculate the number and percentage of of sites in protected areas for each study
# three values of protection status were used using protected areas that established any time points or those that 
# established before the start or late periods of assemblage sampling

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy",
                  "IDIVTS01" = "H:/wubing/iDiv/homogenization_occupancy")
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

sf_use_s2(FALSE)

###########################
# prepared protected areas files (combine points and polygons and simplify) 
# the files were download from the World Database on Protected Areas (January 2022; https://www.protectedplanet.net/en)

## protected areas file 0
wdpa0_pts <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Oct2022_Public_shp/WDPA_Oct2022_Public_shp_0", layer = "WDPA_Oct2022_Public_shp-points")
wdpa0_pys <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Oct2022_Public_shp/WDPA_Oct2022_Public_shp_0", layer = "WDPA_Oct2022_Public_shp-polygons")

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

save(wdpa0_pts, wdpa0_pys, file = "data/WDPA_Oct2022/wdpa0.RDATA")


## protected areas file 1
wdpa1_pts <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Oct2022_Public_shp/WDPA_Oct2022_Public_shp_1", layer = "WDPA_Oct2022_Public_shp-points")
wdpa1_pys <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Oct2022_Public_shp/WDPA_Oct2022_Public_shp_1", layer = "WDPA_Oct2022_Public_shp-polygons")

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

save(wdpa1_pts, wdpa1_pys, file = "data/WDPA_Oct2022/wdpa1.RDATA")


## protected areas file 2
wdpa2_pts <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Oct2022_Public_shp/WDPA_Oct2022_Public_shp_2", layer = "WDPA_Oct2022_Public_shp-points")
wdpa2_pys <- st_read(dsn = "C:/iDiv/Data/WDPA/WDPA_Oct2022_Public_shp/WDPA_Oct2022_Public_shp_2", layer = "WDPA_Oct2022_Public_shp-polygons")

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

save(wdpa2_pts, wdpa2_pys, file = "data/WDPA_Oct2022/wdpa2.RDATA")


# combine protected areas from different subsets for those originally repressed as points and polygons respectively
#  originally repressed as points
wdpa_pts <- bind_rows(wdpa0_pts, wdpa1_pts, wdpa2_pts) %>% 
  st_make_valid()

#  originally repressed as polygons
wdpa_pys <- bind_rows(wdpa0_pys, wdpa1_pys, wdpa2_pys) %>% 
  st_make_valid()

save(wdpa_pts, wdpa_pys, file = "data/WDPA_Oct2022/wdpa_combined.RDATA")


str(wdpa_pts)
plot(wdpa_pts[,1],axes=TRUE)
table(wdpa_pts$IUCN_CAT)
table(wdpa_pys$IUCN_CAT)
summary(wdpa_pys$REP_AREA)

table(wdpa_pys$STATUS_YR == 0)
table(is.na(wdpa_pys$STATUS_YR))
table(wdpa_pys$MARINE)
summary(wdpa_pys$STATUS_YR[wdpa_pys$STATUS_YR >0])
summary(wdpa_pys$STATUS_YR[wdpa_pys$MARINE==1 & wdpa_pys$STATUS_YR >0])



################################
# Calculate the number and proportion of of sites in protected areas
# calculate the proportion of region (convex hull of sampling sites) covered by protected area

## protected area and coordinates of samples
load("data/WDPA_Oct2022/wdpa_combined.RDATA")
load("data/Combined_assemblages.RDATA")
load("intermediate_results/occupancy.RDATA")

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

# combine protected areas from points and polygons
wdpa <- bind_rows(wdpa_pts, wdpa_pys) %>%
  st_make_valid()


###########
# Calculate the number and percentage of of sites in protected areas

# transfer to spatial points
dat_pts <- dat_samp %>%  
  st_as_sf(coords = c('longitude', 'latitude'))  %>% 
  st_set_crs(4326)

# over community samples with protected areas
samp_isinpa <- as_tibble(st_join(dat_pts, wdpa, join = st_intersects))

samp_inpa <- samp_isinpa %>% 
  dplyr::select(study, middle_year, start_year, end_year, start_late, sample, WDPAID, IUCN_CAT, STATUS_YR) %>% 
  filter(!is.na(WDPAID))

# number of sites in protected area
nsamp_inPA_bfstart <- samp_inpa %>% 
  filter(STATUS_YR < start_year) %>%
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


ggplot(dat_psamp_inPA, aes(x = psamp_inPA_wdpa, y = psamp_inPA_bflate)) +
  facet_wrap(~ realm) +
  geom_point() +
  geom_jitter(width = 0.1, height = 0.1)



################################
# calculate the proportion of region (convex hull of sampling sites) covered by protected area

# add metadata and divided 6 study which span 180/-180 longitude to 2 areas (east/west) (not perfect, but OK here to avoided following errors ) 
dat_substudy_samp<- dat_samp %>%
  left_join(dat_meta %>% dplyr::select(study, realm, extent_km2)) %>% 
  group_by(study) %>%
  mutate(long_span = max(longitude) - min(longitude)) %>%
  ungroup() %>%
  mutate(subid = ifelse(long_span> 180 & longitude >1, 2, 1)) %>%
  unite(col = study_sub, study, subid, remove = FALSE) %>%
  group_by(study_sub) %>%
  mutate(ncoord = n_distinct(latitude, longitude)) %>%
  ungroup()

dat_substudy_samp %>% distinct(study_sub, long_span, ncoord) %>% filter(ncoord ==1) # number of regions have only 1 coordinate
dat_substudy_samp %>% distinct(study_sub, long_span, ncoord) %>% filter(ncoord ==2)
dat_substudy_samp %>% distinct(study_sub, long_span, ncoord) %>% filter(long_span >180 & ncoord < 3)
dat_substudy_samp %>% distinct(study_sub, long_span, ncoord) %>% filter(ncoord >= 3)

# for regions with at least 3 points, generate a convex hull represting the region
dat_geom <- dat_substudy_samp %>% 
  filter(ncoord > 2) %>%
  distinct(study_sub, study, realm, start_late, latitude, longitude) %>%
  st_as_sf(coords = c('longitude', 'latitude'))  %>% 
  st_set_crs(4326) %>%
  group_by(study_sub, study, realm, start_late) %>% 
  summarise(geometry = st_combine(geometry)) %>% 
  st_convex_hull() %>% 
  st_make_valid()


# get the polygons of land and ocean, which will be used to intersect the regional polygons
land  <- rnaturalearth::ne_download(scale=10, type="land", category='physical', load=TRUE)
ocean  <- rnaturalearth::ne_download(scale=10, type="ocean", category='physical', load=TRUE)
land <- st_as_sf(land)
ocean <- st_as_sf(ocean)
ocean <- st_make_valid(ocean)

# calculate the intersection between regional polygons and land or ocean extent
# remove marine areas for terrestrial/freshwater regions, and remove land areas for marine regions
dat_geom_realm <- bind_rows(
  st_intersection(dat_geom %>% filter(realm != "Marine"), land),
  st_intersection(dat_geom %>% filter(realm == "Marine"), ocean)
)


# union protected ares, avoiding the overlapping areas 
wdpa_union <- st_union(wdpa)

# calculate the intersection between regional polygons and all protected areas
dat_geom_pa <- st_intersection(dat_geom_realm, wdpa_union)

# the areas covered by protected areas
region_pa_wdpa <- st_drop_geometry(dat_geom_pa) %>%
  mutate(area_pa_wdpa = as.numeric(st_area(dat_geom_pa))/10^6)

# the total areas of regional polygons and proportion covered by protected areas
region_pa_wdpa <- st_drop_geometry(dat_geom_realm) %>%
  mutate(area_total = as.numeric(st_area(dat_geom_realm))/10^6) %>% 
  left_join(region_pa_wdpa) %>%
  mutate(area_pa_wdpa = ifelse(is.na(area_pa_wdpa), 0, area_pa_wdpa)) %>%
  group_by(study, realm, extent_km2, start_late) %>%
  summarise(area_total = sum(area_total),
            area_pa_wdpa = sum(area_pa_wdpa)) %>%
  ungroup() %>%
  mutate(parea_inPA_wdpa = area_pa_wdpa/area_total)

save(wdpa_union, file = "data/WDPA_Oct2022/wdpa_combined_union.RDATA")
save(dat_geom_pa, file = "data/WDPA_Oct2022/Assemblages_polygons_inPA.RDATA")


# calculate the intersection between regional polygons and protected areas that was established before the late periods of sampling
# for each study, appropriate protected areas was combined firstly and then intersected with regional polygons
studies <- pull(dat_geom_realm, study) %>% unique()
region_pa_bflate <- NULL
for(i in 1:length(studies)){
  dat_geom_sub <- dat_geom_realm %>% filter(study == studies[i])
  print(i)
  
  wdpa_sub <- wdpa %>% filter(STATUS_YR < dat_geom_sub$start_late[1])
  wdpa_sub_union <- st_union(wdpa_sub)
  dat_geom_sub_pa <- st_intersection(dat_geom_sub, wdpa_sub_union)
  dat_sub_pa_area <- as.numeric(sum(st_area(dat_geom_sub_pa)))/10^6
  region_pa_bflate <- c(region_pa_bflate, dat_sub_pa_area)
}
save(region_pa_bflate, file = "data/WDPA_Oct2022/Assemblages_areas_inPA.RDATA")

# areas and proportion of regional polygons covered by protected areas
region_pa_bflate1 <- tibble(study = studies, area_pa_bflate = region_pa_bflate)

dat_parea_inPA <- left_join(region_pa_wdpa, region_pa_bflate1) %>%
  mutate(parea_inPA_bflate = area_pa_bflate/area_total)

ggplot(dat_parea_inPA, aes(x = parea_inPA_wdpa, y = parea_inPA_bflate)) +
  facet_wrap(~ realm) +
  geom_point() +
  geom_jitter(width = 0.05, height = 0.05)


######################
# combine protection data based on proportion of sites or regional polygons in protected areas


dat_coverage_PA <- dat_psamp_inPA %>%
  left_join(dat_parea_inPA %>% dplyr::select(- extent_km2)) %>%
  dplyr::select(study, study_name, realm, nsamp_total, ncoord_total, coordinate_local, extent_km2, start_year, start_late,
                psamp_inPA_bfstart, psamp_inPA_bflate, psamp_inPA_wdpa, area_total, parea_inPA_wdpa, parea_inPA_bflate) 

# for studies with only central coordinators
# output to check manually based on descriptions in original papers. 
# If the regions are largely protected, protection level was set as 1; If are not protected, set as 0. 
# For studies that can't determined based on original descriptions, 
# determine protection level based on central coordinates for regions with area <= 10 km2; and remove large regions
dat_coverage_PA_manual <- dat_coverage_PA %>%
  filter(!coordinate_local | is.na(parea_inPA_wdpa))

write.csv(dat_coverage_PA_manual, file = "data/Assemblages_withCentralCoordinates_protection_manual.csv")

dat_coverage_PA_manual <- read_csv( "data/Assemblages_withCentralCoordinates_protection_manual_filled.csv")[, -1]

# update the protection level using the manually corrected ones
id <- match(dat_coverage_PA_manual$study, dat_coverage_PA$study, nomatch = 0)
table(id == 0)
dat_coverage_PA[id, "psamp_inPA_bfstart"] <- dat_coverage_PA_manual$psamp_inPA_bfstart_new
dat_coverage_PA[id, "psamp_inPA_bflate"] <- dat_coverage_PA_manual$psamp_inPA_bflate_new
dat_coverage_PA[id, "psamp_inPA_wdpa"] <- dat_coverage_PA_manual$psamp_inPA_wdpa_new
dat_coverage_PA[id, "parea_inPA_bflate"] <- dat_coverage_PA_manual$parea_inPA_bflate_new
dat_coverage_PA[id, "parea_inPA_wdpa"] <- dat_coverage_PA_manual$parea_inPA_wdpa_new
dat_coverage_PA[id, "comment"] <- dat_coverage_PA_manual$comment


# save results
save(dat_coverage_PA, file = "data/Assemblages_protection.RDATA")

