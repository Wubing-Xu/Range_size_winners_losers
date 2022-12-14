## find the continents or oceans that each region (meta-community) located

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse", "sp", "rworldmap", "rnaturalearth", "oceanmap")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


## load coordinates of studies
load("data/Combined_assemblages.RDATA")

# download and prepare spatial data
ocean <- rnaturalearth::ne_download(scale=10, type="geography_marine_polys", category='physical', load=TRUE)
country <- getMap(resolution='low')
study_points <-  SpatialPoints(data.frame(dat_meta$cent_long, dat_meta$cent_lat), proj4string=CRS(proj4string(country)))  

# use over to get names of polygons containing each point
region <- dat_meta %>% 
  mutate(continent = over(study_points, country)$REGION, # return the continent (7 continent model)
         ocean = over(study_points, ocean)$name_en) # return the name of  ocean or sea 

region %>% filter(study_name == "sperandii_2020_Coastal dunes, Italy") %>% as.data.frame()

# visual check distributions of marine and non-marine studies 
world = map_data('world')
# oceans
ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray", size=0) +
  geom_point(data = filter(region, realm == "Marine"),
             aes(x = cent_long, y = cent_lat, colour = ocean), size = 1, stroke = 1.2)

ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray", size=0) +
  geom_point(data = filter(region, realm == "Marine" & is.na(ocean)),
             aes(x = cent_long, y = cent_lat), size = 0.5, stroke = 1.2)

# continents
ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray", size=0) +
  geom_point(data = filter(region, realm != "Marine"),
             aes(x = cent_long, y = cent_lat, colour = continent), size = 1, stroke = 1.2, alpha = 0.9) 

ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray", size=0) +
  geom_point(data = filter(region, realm != "Marine" & is.na(continent)),
             aes(x = cent_long, y = cent_lat), size = 0.5, stroke = 1.2) 

# marine studies to be checked
region %>% 
  filter(realm == "Marine" & is.na(ocean)) %>%
  dplyr::select(study:climate, continent, ocean, cent_long, cent_lat)

# terrestrial or freshwater studies to be checked. Most are located in islands
region %>% 
  filter(realm != "Marine" & is.na(continent)) %>%
  dplyr::select(study:climate, continent, ocean, cent_long, cent_lat)

# output to manually check, and than input
write_csv(region, file = "data/manual_continents_ocean.csv")
region <- read_csv("data/manual_continents_ocean_filled.csv")


# match names of sea to oceans in five models
table(region$ocean)
sort(unique(region$ocean))
ocean_manual <- data.frame(raw = c("Arctic", "Atlantic Ocean", "Baltic Sea", "Bay of Bengal", "Bay of Biscay",    
                                   "Bering Sea", "Black Sea", "Caribbean Sea", "Coral Sea", "Gulf of Bothnia",
                                   "Gulf of Maine", "Gulf of Saint Lawrence", "Gulf of Thailand", "Indian Ocean", "Irish Sea",         
                                   "Mediterranean Sea",  "Mozambique Channel", "North Sea",  "Pacific Ocean", "Philippine Sea",  
                                   "Salish Sea", "Sea of Japan", "Sea of the Hebrides", "Skagerrak",  "Southern Ocean",  
                                   "Tasman Sea", "Tyrrhenian Sea", "White Sea", "Cook Strait", "Arctic Ocean", 
                                   "Great Barrier Reef", "Mackenzie Bay"),
                           new = c("Arctic Ocean", "Atlantic", "Atlantic", "Indian Ocean", "Atlantic", 
                                   "Pacific", "Atlantic", "Atlantic", "Pacific", "Atlantic", 
                                   "Atlantic", "Atlantic", "Pacific", "Indian Ocean", "Atlantic", 
                                   "Atlantic", "Indian Ocean", "Atlantic", "Pacific", "Pacific", 
                                   "Pacific", "Pacific", "Atlantic", "Atlantic", "Southern Ocean", 
                                   "Pacific", "Atlantic", "Arctic Ocean", "Pacific", "Arctic Ocean", 
                                   "Pacific", "Arctic Ocean"))
# all is matched
sort(unique(region$ocean))[!sort(unique(region$ocean)) %in% ocean_manual$raw]

# update ocean names
id <- match(pull(region, ocean), ocean_manual[,1], nomatch = NA)
region$ocean[!is.na(id)] <- ocean_manual[id[!is.na(id)], 2]

region <- region %>% 
  mutate(region =  ifelse(realm == "Marine", ocean, continent))

save(region, file = "data/Assemblages_regions.RDATA")

