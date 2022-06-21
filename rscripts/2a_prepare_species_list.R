## Combine species from four assemblage database, add and unify realm, kingdom, is_bird, and match to GBIF taxonomy backbone
# the matched GBIF species names were used to extact GBIF occurrences

rm(list = ls())
setwd("C:/Dropbox/iDiv/homogenization_occupancy")

library(tidyverse)
library(rgbif)

####################
#### combine all species from four database

#### BioTIME species list

load("data/BioTIME/BioTIME_4locations_10years_filtered_final.RDATA")

# species list in BioTIME
bt_species <- bt_filtered %>% 
  distinct(STUDY_ID, GENUS_SPECIES, REALM, HABITAT, TAXA, ORGANISMS)

# add kingdom according to descriptions in TAXA and ORGANISMS
id_plantae <- bt_species$TAXA %in% c("Terrestrial plants", "Marine plants", "Freshwater plants") &
  !bt_species$ORGANISMS %in% c("Macroalgae", "fungi", "moss", "lichen", "tropical algae")
table(id_plantae)

id_Animalia <- (bt_species$TAXA %in% c("Birds", "Fish", " Mammals", "Terrestrial invertebrates", "Marine invertebrates", 
                                       "Freshwater invertebrates", "Reptiles", "Amphibians") &
                  !bt_species$ORGANISMS %in% c("intertidal", "marine invertebrates & marine plants", "Microplankton", "Phytoplankton",
                                               "Phytoplankton + zooplankton (and some plants - changed taxa from /plants)", "turf algae")) |
  (bt_species$TAXA %in% c("Benthos") & bt_species$ORGANISMS %in% c("benthic animals", "Echinodermata",  "Epibenthic megafauna", 
                                                                   "macrozoobenthos", "Polychaeta","zoobenthos")) |
  (bt_species$TAXA %in% c("All") & bt_species$ORGANISMS %in% c("All animals observed around waterholes", "birds and some marine mammals", "Birds carnivores and ungulates",                                                             
                                                               "Cetaceans. seabirds + turtles", "fish and marine invertebrates", "Fish and marine invertebrates",  
                                                               "Herbivores carnivores and eagles", "Herpetofauna", "mostly seabirds + some marine mammals", "pelagic seabirds"))
table(id_Animalia)

bt_species[id_plantae, "kingdom"] <- "Plantae"
bt_species[id_Animalia, "kingdom"] <- "Animalia"

# the BioTIME species list
btsplist <- bt_species %>% 
  distinct(GENUS_SPECIES, REALM, TAXA, ORGANISMS, kingdom) %>%
  setNames(tolower(names(.))) %>% 
  rename(species=genus_species, taxon=taxa) %>% 
  mutate(realm = tolower(realm))

# For birds, we need their habitats to determine their ranges. Thus this taxonomic groups are indicated particularly
# check whether the species are birds for the taxon group named as "all"
check.bird <- btsplist %>% filter(taxon=="All") # & grepl(pattern="birds", unlist(.[,"organisms"]))
check.bird_gbif <- lapply(1:nrow(check.bird),function(x) name_backbone(name = check.bird[x,1], kingdom="Animalia"))
check.bird_gbif <- bind_rows(check.bird_gbif)
check.bird <- check.bird %>% mutate(class = check.bird_gbif$class) %>% 
  mutate(isbird = ifelse(class!="Aves" | is.na(class),"notbird","bird"))

# Add the column "isbird" for the checklist
btsplist <- btsplist %>% left_join(check.bird) %>% 
  mutate(isbird = ifelse(taxon=="Birds" | (isbird=="bird" & !is.na(isbird)), "bird","notbird")) %>%
  dplyr::select(species, realm, taxon, isbird, kingdom) %>% 
  distinct()



#### checklist of biotic homogenization database

load("data/Homogenization/homogenization_filtered.RDATA")

# one dataset in the meta data have multiple taxon types in the same region. Name it as "multiple"
bh_meta %>% distinct(dataset_id, regional, realm, taxon) %>% filter(duplicated(.[,c("dataset_id", "regional")]))
bh_meta %>% distinct(dataset_id, realm, taxon) %>% filter(duplicated(.[,c("dataset_id")])) 
bh_meta %>% distinct(dataset_id, regional, realm, taxon) %>% filter(dataset_id == "sorte_2018")
# bh_meta <- bh_meta %>% mutate(taxon = ifelse(dataset_id == "sorte_2018", "multiple", taxon))

# species list
bh_species <- bh_filtered %>% distinct(dataset_id, species, realm, taxon)

# check whether all species have realm and taxon information
bh_species %>% filter(is.na(realm))

# add kingdom to species list
unique(bh_species$taxon)
taxon_kingdom <- data.frame(taxon = c("Fish", "Invertebrates", "Marine plants", "Herpetofauna",  "Birds", "Plants", "Mammals", "Multiple taxa"),
                            kingdom = c("Animalia", "Animalia", "Plantae", "Animalia", "Animalia", "Plantae", "Animalia", NA ))

bh_species <- bh_species %>% left_join(taxon_kingdom)

# the prepared species list for biotic homogenization database
bhsplist <- bh_species %>% 
  mutate(isbird = ifelse(taxon=="Birds","bird","notbird")) %>%
  dplyr::select(species, realm, taxon, isbird, kingdom) %>% 
  distinct()


#### RivFishTIME species checklist
load("data/RivFishTIME/RivFishTIME_4sites_10years_filtered.RDATA")

ftsplist <- ft_filtered %>% distinct(species = Species) %>%
  mutate(realm = "freshwater", taxon = "fish", isbird = "notbird", kingdom = "Animalia")


#### species checklist of insect dataset
load("data/Role_insects/inset_metacommunities_filtered.RDATA")

itsplist <- it_filtered %>% 
  distinct(species, realm = Realm) %>%
  mutate(taxon = "invertebrate", isbird = "notbird", kingdom = "Animalia")



##combine checklist
splist <- bind_rows(bhsplist, btsplist, ftsplist, itsplist) %>% distinct()


# update "isbird" for those species with inconsistent information of the same species 
check.isbird <- splist %>% 
  distinct(species, isbird) %>% 
  filter(species %in% pull(.[duplicated(species),], species))

check.isbird_gbif <- lapply(1:nrow(check.isbird ),function(x) name_backbone(name=check.isbird [x, 1], kingdom="Animalia"))
check.isbird_gbif <- bind_rows(check.isbird_gbif)

check.isbird <- check.isbird %>% 
  mutate(class = check.isbird_gbif$class) %>% 
   mutate(isbird.new = ifelse(class=="Aves", "bird", "notbird")) %>% 
   mutate(isbird.new = ifelse(is.na(isbird.new), "unknown", isbird.new)) %>% 
   distinct(species, isbird.new)

# use the updated to replace the raw
splist <- splist %>% left_join(check.isbird) %>% 
  mutate(isbird = ifelse(is.na(isbird.new), isbird, isbird.new)) %>% 
  dplyr::select(!isbird.new) %>% 
  distinct()


# update "realm" for those species with inconsistent information of the same species 
splist <- splist %>%
  mutate(realm = tolower(realm))

check.realm <- splist %>% 
  distinct(species, realm) %>% 
  filter(species %in% pull(.[duplicated(species),], species)) %>%
  mutate(realm.new = NA)

# combine realms for the same species
check.realm.splist <- distinct(check.realm, species) %>% unlist()
for(i in 1:length(check.realm.splist)){
  id <- check.realm[,1] == check.realm.splist[i]
  check.realm$realm.new[id]<- paste(sort(unlist(check.realm[id, 2])), collapse="_")
}
check.realm <- check.realm %>% distinct(species, realm.new)

# use the updated to replace the raw
splist <- splist %>% left_join(check.realm) %>% 
  mutate(realm = ifelse(is.na(realm.new), realm, realm.new)) %>% 
  dplyr::select(!realm.new) %>% distinct()


# update "kingdom" for those species with inconsistent information of the same species 
check.kingdom <- splist %>% 
  distinct(species, kingdom) %>% 
  filter(species %in% pull(.[duplicated(species),],species)) %>%
  mutate(kingdom.new = NA)

check.kingdom.splist <- distinct(check.kingdom, species) %>% unlist()
for(i in 1:length(check.kingdom.splist)){
  id <- check.kingdom[,1] == check.kingdom.splist[i]
  kingdom <- unlist(check.kingdom[id, 2])
  kingdom <- kingdom[!is.na(kingdom)]
  check.kingdom$kingdom.new[id] <- ifelse(length(kingdom) ==1, kingdom, NA)
}

check.kingdom <- check.kingdom %>% 
  distinct(species, kingdom.new) %>%
  mutate(kingdom.new = ifelse(is.na(kingdom.new), "unknown", kingdom.new))

# use the updated to replace the raw
splist <- splist %>% left_join(check.kingdom) %>% 
  mutate(kingdom = ifelse(is.na(kingdom.new), kingdom, kingdom.new)) %>% 
  dplyr::select(!kingdom.new) %>% distinct()


# the final splist
splist <- distinct(splist, species, realm, isbird, kingdom)



#########
# find names in gbif backbone
spgbif <- lapply(1:nrow(splist),function(x) name_backbone(name=splist[x,1], kingdom=splist[x, 4]))
spgbif <- bind_rows(spgbif)

# Combine species checklist with information from gbif backbone
splist.gbif.all <- splist %>% 
  rename(input.species = species, input.kingdom = kingdom) %>% 
  bind_cols(spgbif %>% setNames(tolower(names(.))))

# use the class from gbif to update the raw "isbird". 16 species are updated
id <- ((splist.gbif.all[,"isbird"] == "bird" & splist.gbif.all[,"class"] != "Aves") | 
      (splist.gbif.all[,"isbird"] != "bird" & splist.gbif.all[,"class"] == "Aves")) &
      !is.na(splist.gbif.all[,"class"]) & !is.na(splist.gbif.all[,"species"])
table(id)
splist.gbif.all[id,"isbird"] <- ifelse(splist.gbif.all[id, "isbird"] == "bird","notbird", "bird")

# the simplified species list with names in gbif backbone
splist.gbif <- splist.gbif.all %>% select(c(input.species:input.kingdom, species, specieskey))

splist.gbif.specieskey <- splist.gbif %>% filter(!is.na(specieskey)) %>% distinct(specieskey) %>% pull()
length(splist.gbif.specieskey) # 17,747 species keys


# output
save(splist.gbif.all, file="data/combined_checklists/splist_gbif_all.columns.RDATA")
save(splist.gbif, splist.gbif.specieskey, file="data/combined_checklists/splist_gbif.RDATA")
write_csv(splist.gbif.all, file="data/combined_checklists/splist_gbif_all.columns.csv") #save as a csv to inspect

