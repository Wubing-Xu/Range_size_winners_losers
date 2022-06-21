####################################
## A example showing how we calculate range sizes

rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(dplyr)
library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(sf)

load("data/gbif/ranges/range_details_0127134-210914110416597.RDATA")
load("data/gbif/distribution_records/occ_0127134-210914110416597.RDATA")

land  <- rnaturalearth::ne_download(scale=10, type="land", category='physical', returnclass = 'sf')
behr <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
land_behr <- st_transform(land, crs = behr)


## prepare spatial polygons to show grid-cells that occupied by a species
values(ras100) <- 1
shp100 <- rasterToPolygons(ras100) 
colnames(shp100@data)[1] <- "ID"
shp100$ID <- which(as.vector(!is.na(ras100)))
shp100 <- st_as_sf(shp100)

# species list 
splist <- distinct(occ_rasid, species, specieskey)


## An example showing area of occupancy (AOO)
#  the id of grid-cells that occupied by an example species
occ_rasid_examp <- occ_rasid %>% 
  filter(specieskey == splist$specieskey[23]) %>%
  dplyr::select(species, specieskey, ras100id) %>%
  distinct()

AOO <- ggplot() +
  geom_sf(data = land_behr, fill = "gray60", colour = "gray60", size = 0.2) + 
  geom_sf(data = shp100 %>% filter(ID %in% occ_rasid_examp$ras100id), fill = "red", colour = "red", size = 0.01) + 
  labs(tag = "A", title  = "Area of occupancy (AOO): number of grid-cells") +
  theme_classic() + 
  theme(plot.tag.position = c(0, 0.98),
        plot.title = element_text(size = 10, hjust = 0.05),
        plot.tag = element_text(size = 10, face = "bold"))


## An example showing extent of occurrences (EOO) defined as the area of alpha hulls
# distribution occurrences of an example species
occ_examp <- occ %>% 
  filter(specieskey == splist$specieskey[23]) %>%
  distinct(specieskey, species, decimallongitude, decimallatitude) %>%
  st_as_sf(coords = c('decimallongitude', 'decimallatitude'))  %>% 
  st_set_crs(4326) %>%
  st_transform(crs = behr)

# the alpha hull of an example spcies
ahull_examp <- ahull6[ahull6@data$species == splist$specieskey[23], ] %>%
  st_as_sf() %>%
  st_transform(crs = behr)

EOO <- ggplot() +
  geom_sf(data = land_behr, fill = "gray60", colour = "gray50", size = 0.2) + 
  geom_sf(data = ahull_examp, fill = "blue", colour = "blue", alpha = 1) + 
  geom_sf(data = occ_examp, colour = "red", alpha = 0.5, size = 0.01) + 
  labs(tag = "B", title  = "Extent of occurrences (EOO): area of alpha hulls") +
  theme_classic() + 
  theme(plot.tag.position = c(0, 0.98),
        plot.title = element_text(size = 10, hjust = 0.05),
        plot.tag = element_text(size = 10, face = "bold"))

AOO_EOO_<- cowplot::plot_grid(AOO, EOO, nrow = 2)
ggsave(AOO_EOO_, file = "results/Fig.S_diagram_RangeSize.png", 
       width = 150, height = 125, units = 'mm', dpi = 600)
