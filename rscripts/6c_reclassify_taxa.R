##################################
## reclassify the taxonomic groups for those that are listed as "multiple taxa" and "Benthos".
# check the taxonomic rank of species in my filtered data using the GBIF backbone
# the taxonomic group with most species within a study will be used  

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/homogenization_occupancy",
                  "IDIVTS02" = "D:/ya90meri/homogenization_occupancy")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","vegan", "reshape2", "sf")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


load("data/Combined_assemblages_20220208.RDATA")
load("data/combined_checklists/splist_gbif_all.columns.RDATA")


#####################
## the studies are listed as "multiple taxa".
meta_mtaxa <-  dat_meta %>% 
    filter(taxon_new == "Multiple taxa")
meta_mtaxa$study

species_dat_mtaxa <- dat %>%
  filter(study %in% meta_mtaxa$study) %>%
  distinct(study, species, specieskey) %>%
  left_join(splist.gbif.all %>% 
              dplyr::select(specieskey, kingdom:family) %>%
              distinct())


# study bt_166: birds
species_dat_mtaxa %>% filter(study == "bt_166") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_166") %>% 
  count(order)
# order: 78 birds (47 Charadriiformes, 14 Procellariiformes, 10 Anseriformes, 7 others), 
# 23 mammals (18 Cetacea; 5 Carnivora)


# study bt_169: birds
species_dat_mtaxa %>% filter(study == "bt_169") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_169") %>% 
  count(order)
# order: 125 birds (52 Charadriiformes, 26 Procellariiformes, and others), 
# 27 mammals (20 Cetacea, 6 Carnivora, 1 Carcharhiniformes)


# study bt_172: mammals
species_dat_mtaxa %>% filter(study == "bt_172") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_172") %>% 
  count(order)
# order: 20 mammals (20 Cetacea), 
# 9 birds (4 Charadriiformes, 4 Procellariiformes, 1 Columbiformes); 
# 4 turtles (Testudines)


# study bt_182: fishes
species_dat_mtaxa %>% filter(study == "bt_182") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_182") %>% 
  count(phylum)

species_dat_mtaxa %>% filter(study == "bt_182") %>% 
  count(order)
# order: fishes: 18 fishes (6 Pleuronectiformes,5 Scorpaeniformes, 3 Perciformes, 2 Gadiformes, Lophiiformes, Rajiformes)
# invertebrates: 13 invertebrates (4 Decapoda, Amphilepidida, Camarodonta, Clypeasteroida, Dendrochirotida, Euryalida, Mytilida, Myxiniformes, Pectinida, Phyllodocida)
# Chordariales: 1


# study bt_274: invertebrates
species_dat_mtaxa %>% filter(study == "bt_274") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_274") %>% 
  count(phylum)
# phylum: invertebrate: 49 species ( 12 Cnidaria, 10 Bryozoa, 7 Annelida,  6 Chordata, 5 Echinodermata, 5 Mollusca, 4 Porifera)
# plants (alga): 31 species (17 Rhodophyta,  10 Ochrophyta, 2 Chlorophyta, 2 Tracheophyta)


# study bt_505: fishes
species_dat_mtaxa %>% filter(study == "bt_505") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_505") %>% 
  count(phylum)

species_dat_mtaxa %>% filter(study == "bt_505") %>% 
  count(order) %>% as.data.frame()
# 118 Chordata (most are fishes, 66 Perciformes,11 Scorpaeniformes, 7 Clupeiformes), 
# 16 Decapoda, 12 Mollusca, 4 Echinodermata, 2 Cnidaria


# study bt_527: birds
species_dat_mtaxa %>% filter(study == "bt_527") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_527") %>% 
  count(phylum)

species_dat_mtaxa %>% filter(study == "bt_527") %>% 
  count(order) %>% as.data.frame()
# phylum: 66 birds (32 Charadriiformes, 14 Anseriformes, 11 Procellariiformes, 4 Suliformes, 3 Gaviiformes, 1 Pelecaniformes, 1 Podicipediformes)
# 28 mammals (20 Cetacea; 8 Carnivora)



#########################
## the studies are listed as "Benthos".

meta_benthos <-  dat_meta %>% 
  filter(taxon_new == "Benthos")
meta_benthos$study
# "bt_78"  "bt_92"  "bt_110" "bt_162" "bt_163" "bt_187" "bt_196" "bt_204" "bt_213" "bt_468" "bt_469"

species_dat_benthos <- dat %>%
  filter(study %in% meta_benthos$study) %>%
  distinct(study, species, specieskey) %>%
  left_join(splist.gbif.all %>% 
              dplyr::select(specieskey, kingdom:family) %>%
              distinct())

# study bt_78: most are invertebrates
species_dat_benthos %>% filter(study == "bt_78") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_78") %>% 
  count(phylum)


# study bt_92: most are invertebrates
species_dat_benthos %>% filter(study == "bt_92") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_92") %>% 
  count(phylum)


# study bt_110: most are invertebrates
species_dat_benthos %>% filter(study == "bt_110") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_110") %>% 
  count(phylum)


# study bt_162: most are invertebrates
species_dat_benthos %>% filter(study == "bt_162") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_162") %>% 
  count(phylum)


# study bt_163: most are fish
species_dat_benthos %>% filter(study == "bt_163") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_163") %>% 
  count(phylum)

species_dat_benthos %>% filter(study == "bt_163") %>% 
  count(order) %>% as.data.frame()


# study bt_187: most are invertebrates
species_dat_benthos %>% filter(study == "bt_187") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_187") %>% 
  count(phylum)


# study bt_196: about half are invertebrates
species_dat_benthos %>% filter(study == "bt_196") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_196") %>% 
  count(phylum)


# study bt_204: all are invertebrates
species_dat_benthos %>% filter(study == "bt_204") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_204") %>% 
  count(phylum)


# study bt_213: most are fishes
species_dat_benthos %>% filter(study == "bt_213") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_213") %>% 
  count(phylum)

species_dat_benthos %>% filter(study == "bt_213") %>% 
  count(order)  %>% as.data.frame()


# study bt_468: about half are invertebrates
species_dat_benthos %>% filter(study == "bt_468") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_468") %>% 
  count(phylum)


# study bt_469: most are invertebrates
species_dat_benthos %>% filter(study == "bt_469") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_469") %>% 
  count(phylum)



#############
# to update taxonomic groups
dat_meta_mtaxa_manul <- data.frame(study= c("bt_166", "bt_169", "bt_172", "bt_182", "bt_274", "bt_505", "bt_527"),
                                  taxon_new = c("Birds", "Birds", "Mammals", "Fish", "Invertebrates", "Fish", "Birds"))


dat_meta_benthos_manual <- data.frame(study= c("bt_78","bt_92","bt_110", "bt_162", "bt_163", "bt_187", 
                                             "bt_196", "bt_204", "bt_213", "bt_468", "bt_469"),
                                  taxon_new = c("Invertebrates", "Invertebrates", "Invertebrates", "Invertebrates", "Fish", "Invertebrates", 
                                                "Invertebrates", "Invertebrates", "Fish", "Invertebrates", "Invertebrates"))

# taxon information
dat_meta_taxa <- dat_meta %>%
  dplyr::select(study, database, studyID, study_name, realm, taxon, taxon_new)

# update the "multiple taxa"
id <- match(dat_meta_mtaxa_manul[, "study"], dat_meta_taxa$study)
dat_meta_taxa[id, c("taxon_new")] <- dat_meta_mtaxa_manul[, c("taxon_new")]

# update the "benthos"
id <- match(dat_meta_benthos_manual[, "study"], dat_meta_taxa$study)
dat_meta_taxa[id, c("taxon_new")] <- dat_meta_benthos_manual[, c("taxon_new")]


# distinguish plants, invertebrates and fish in different realms (terrestrial, freshwater and marine)
dat_meta_taxa <- dat_meta_taxa %>% 
  mutate(taxon_final = ifelse(taxon_new %in% c("Fish", "Invertebrates", "Plants"), paste(realm, tolower(taxon_new), sep = " "), taxon_new))
table(dat_meta_taxa$taxon_final)

save(dat_meta_taxa, file = "data/Assemblages_taxa.RDATA")

