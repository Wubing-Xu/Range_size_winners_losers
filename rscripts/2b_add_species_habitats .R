###########################
# Add bird habitat data to gbif names and determine which habitat (land or ocean) is the main habitat for each species

rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(rgbif)

load("data/combined_checklists/splist_gbif.RDATA")
habitat <- read_csv("data/bird_habitats/bird_IUCN_dat.csv")

# checklist with only names from GBIF. removed duplicated
spgbif <- splist.gbif %>% 
  filter(!is.na(species)) %>% 
  dplyr::select(species, specieskey, realm, isbird) %>% 
  distinct()


# update "realm" for those species with inconsistent information of the same species 
check.realm <- spgbif %>% 
  distinct(specieskey, realm) %>% 
  filter(specieskey %in% pull(.[duplicated(specieskey),], specieskey)) %>%
  mutate(realm.new = NA)

# combine realms for the same species
check.realm.splist <- distinct(check.realm, specieskey) %>% unlist()
for(i in 1:length(check.realm.splist)){
  id <- check.realm[,1] == check.realm.splist[i]
  realm.new <- paste(sort(unique(unlist(strsplit(unlist(check.realm[id,2]), "_")))),collapse="_")
  check.realm$realm.new[id]<- realm.new
}
check.realm <- check.realm %>% distinct(specieskey, realm.new)

# use the updated to replace the raw
spgbif <- spgbif %>% 
  left_join(check.realm) %>% 
  mutate(realm = ifelse(is.na(realm.new), realm, realm.new)) %>% 
  dplyr::select(!realm.new) %>% distinct()


# update "isbird" for those species with inconsistent information of the same species 
check.isbird <- spgbif %>% 
  distinct(specieskey, isbird) %>% 
  filter(specieskey %in% 
           pull(.[duplicated(specieskey),], specieskey))
check.isbird  # no species have inconsistent information


## match bird names from gbif with names in the iucn habitat
# find gbif names for the names in the iucn habitat
hbtgbif <- lapply(1:nrow(habitat),function(x) name_backbone(name=habitat[x,"Sciname"],class = "Aves"))
hbtgbif <- bind_rows(hbtgbif)
habitat <- hbtgbif %>% dplyr::select(species, speciesKey) %>% bind_cols(habitat)
# write_csv(habitat, file = "Data/bird_habitats/bird_IUCN_dat_AddGbifSpeciesNames.csv")
# habitat <- read_csv("Data/bird_habitats/bird_IUCN_dat_AddGbifSpeciesNames.csv")

# match bird names with bird habitat data
birdgbif <- spgbif %>% filter(isbird == "bird")
id <- pull(birdgbif, specieskey) %in% pull(habitat, speciesKey)
table(id) #24 bird names can't be matched with bird habitat data
birdgbif <- birdgbif %>% 
  left_join(habitat %>% distinct(species, speciesKey, .keep_all = TRUE), by=c("species", "specieskey"="speciesKey"))

# determine haibitat (land or ocean)
birdgbif <- birdgbif %>% 
  mutate(habitat = ifelse((Habitat_Caves == 1 | Habitat_Desert == 1| Habitat_Forest == 1 | Habitat_Grassland == 1 | 
                             Habitat_Savanna == 1 |Habitat_Shrubland == 1) & (Habitat_Marine_Neritic == 0 | Habitat_Marine_Oceanic == 0), "land",
                          ifelse((Habitat_Caves == 0 | Habitat_Desert == 0| Habitat_Forest == 0 | Habitat_Grassland == 0 | 
                                    Habitat_Savanna == 0 |Habitat_Shrubland == 0) & (Habitat_Marine_Neritic == 1 | Habitat_Marine_Oceanic == 1), 
                                 "ocean", "land_ocean")))

# manual classification and fill missing data
write_csv(birdgbif, file="Data/bird_habitats/manual_birdgbif_habitat.csv")
birdgbif_checked <- read_csv("Data/bird_habitats/manual_birdgbif_habitat_filled.csv")

birdgbif <- birdgbif %>% 
  dplyr::select(-habitat) %>%
  left_join(birdgbif_checked %>% dplyr::select(species, specieskey, habitat))


# check whether species lack habitat data 
birdgbif %>% filter(is.na(habitat)) %>% dplyr::select(species, habitat)

# add bird habitat data to checklist
spgbif <- spgbif %>% 
  left_join(birdgbif %>% dplyr::select(c(species:isbird, habitat)) %>% distinct()) %>% 
  rename(bird.habitat = habitat)


# determine which areas (land or ocean) are the main habitat. It used to filter GBIF occurrences and clip ranges
spgbif <- spgbif %>%
  mutate(keep = ifelse((isbird == "notbird" & realm %in% c("terrestrial", "freshwater", "freshwater_terrestrial")) | 
                            (isbird == "bird" & bird.habitat == "land"), "land", NA)) %>%
  mutate(keep = ifelse((isbird == "notbird" & realm == "marine") | 
                         (isbird == "bird" & bird.habitat == "ocean"), "ocean", keep)) %>%
  mutate(keep = ifelse(isbird == "bird" & bird.habitat == "land_ocean", "both", keep))

table(spgbif$keep, useNA = "always") #75 species missing data; 45 species were recorded as "both" 


save(spgbif, file="Data/combined_checklists/spgbif_habitat.RDATA")

