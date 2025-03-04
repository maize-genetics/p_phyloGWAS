#'##################################################################################################'
#'                                                                                                 #'
#' Project........Environmental GWAS (eGWAS) for PandAnd -- Hackathon                              #'
#' Title..........Geographic data from BIEN                                                        #'
#' Created at..... 11-06-2023                                                                      #'
#' Updated at..... 06-24-2024                                                                      #'
#' Author: G.Costa-Neto, germano.cneto<at>gmail.com                                                #'
#' Sheng-Kai Hsu re-write into function                                                            #'
#'##################################################################################################'
rm(list = ls())

require(tidyverse)
require(plyr)
require(reshape2)
require(foreach)

require(BIEN) # to access BIEN data base
require(rgbif) # to access GBIF data base

args = commandArgs(trailingOnly=TRUE)
# arg 1: path to species list (long format); arg 2: output directory

# dir.data = "/workdir/sh2246/p_phyloGWAS/output/testMetadata.txt"
# dir.output = "/workdir/sh2246/p_phyloGWAS/output/tMetadata_testout"
dir.data = args[1]
dir.output = args[2]


#'------------------------------------------------------------------------------------------------------------
# (1) importing the species names / ids ######
#'------------------------------------------------------------------------------------------------------------

# this is the reference file with the species names
androMetadata  <- read.delim(dir.data,header = T)
# androMetadata = androMetadata[1:20,]
# androMetadata %>% View

#'------------------------------------------------------------------------------------------------------------
# (2) using rgbif to access Global Biodiversity Information Facility, GBIF https://www.gbif.org/ ######
#'------------------------------------------------------------------------------------------------------------

missing_taxa = unique(androMetadata$names) %>% na.omit() %>% as.character() %>% sort()
missing_taxa # vector of species names. You can check it and edit it before running the codes below

# let's create a directory for saving .csv files from GBIF

path_species = paste0(dir.output,'/species_metadata_gbif') # where to save it

if(!dir.exists(path_species))dir.create(path_species)

gbif_data <- foreach::foreach(j = 1:length(missing_taxa),  .packages = c('rgbif','tidyverse'),.combine = 'rbind') %dopar%
  { print(paste("GBIF:",j, "out of",length(missing_taxa)))
    # rgbif::occ_search returns a list of features. the $data is what we are interested in.
    
    data = rgbif::occ_search(scientificName = missing_taxa[j],hasCoordinate = T)$data %>% as.data.frame()
    data_name = gsub(missing_taxa[j],pattern = ' ',replacement='_')
    
    
    # pulling only species ID, data key and coordinates
    if(nrow(data) == 0)
    {
      output <-  data.frame(scientificName=missing_taxa[j],decimalLatitude=NA,decimalLongitude=NA)
      data.table::fwrite(file = paste0(path_species,'/missing_species.csv'), data.frame(species=missing_taxa[j]),append = T)
    }
    
    if(nrow(data) >0)
    {
      # saving all metadata
      data.table::fwrite(file = paste0(path_species,'/',data_name,'.tsv'), data,sep ="\t")
      
      output <-   data.frame(scientificName=missing_taxa[j],decimalLatitude=data$decimalLatitude,decimalLongitude=data$decimalLongitude)
      #data[,c('key','scientificName','decimalLatitude','decimalLongitude')]
    }
    
    
    # merging everything ;)
    return(output)
  }

# to register as a derived dataset
gbifData = data.table::fread("/workdir/sh2246/p_phyloGWAS/output/metadataFormalOut/species_metadata_gbif_merged.tsv",fill=Inf,sep = "\t",na.strings = "")
derived_data <- gbifData[,6] %>%
  group_by(datasetKey) %>% 
  count()

# derived_dataset(user = "shengkaihsu",pwd = "***!",
#                 citation_data = derived_data,
#                 title = "occurrence data for 614 Poaceae species with associated genome assemblies",
#                 description = "This dataset is obtain using rgbif::occ_search() function for 614 Poaceae species.\n
#                 In the subsequent analysis, this dataset is combined with occurrence data from BIEN database and furhter cleaned with CoordinateCleaner package.",
#                 source_url = "https://doi.org/10.5281/zenodo.14967967")

#'------------------------------------------------------------------------------------------------------------
# (3) using BIEN to access Botanical Information and Ecology Network , https://bien.nceas.ucsb.edu/bien/ #####
#'------------------------------------------------------------------------------------------------------------

# path where to save it
path_species2 = paste0(dir.output,'/species_metadata_bien') # where to save it

if(!dir.exists(path_species2)) dir.create(path_species2)

# here I make an example without foreach. using conventional "for()"
bien_data <- c()
for(j in 1:length(missing_taxa)){ 
  print(paste("BIEN:",j, "out of",length(missing_taxa)))
  
  data = try(BIEN::BIEN_occurrence_species(species = missing_taxa[j],cultivated = FALSE,only.geovalid = TRUE) )
  
  data_name = gsub(missing_taxa[j],pattern = ' ',replacement='_')
  
  # pulling only species ID, data key and coordinates
  if(nrow(data) == 0)
  {
    output <-  data.frame(scientificName=missing_taxa[j],decimalLatitude=NA,decimalLongitude=NA)
    data.table::fwrite(file = paste0(path_species2,'/missing_species.csv'), data.frame(species=missing_taxa[j]),append = T)
  }
  
  if(nrow(data) >0)
  {
    # saving all metadata
    data.table::fwrite(file = paste0(path_species2,'/',data_name,'.tsv'), data,sep = "\t")
    
    output <-   data.frame(scientificName=missing_taxa[j],decimalLatitude=data$latitude,decimalLongitude=data$longitude)
    #data[,c('key','scientificName','decimalLatitude','decimalLongitude')]
  }
  bien_data <- rbind(bien_data,    output)
}

#'------------------------------------------------------------------------------------------------------------
# (4)Merging & Coordinate Cleaning ######
#'------------------------------------------------------------------------------------------------------------
# merging
bien_data$package = "BIEN"
gbif_data$package = "GBIF"

merge_data = rbind(bien_data,gbif_data)

# filtering -- there is other options and packages for it, but this is the main QC you can do
merge_data <-
  merge_data %>% 
  mutate(sum = decimalLatitude+decimalLongitude  ) %>%  
  mutate(diff = decimalLatitude-decimalLongitude  ) %>%  
  filter(!is.na(decimalLatitude)) %>% 
  filter(!sum == 0 |!diff == 0 ) 

merge_data <- merge_data [,1:4]   


# generating "keys" for the samples of a given specie AND coordiante
# for instance: specie A has 10 samples, but 3 of those 10 has the same coordiante. 
# So we will chose only 1 from those 3 (key 1, Key 2, key 3, we choose only key 1)
## comment: I'm doing this because there is many samples with the same coordinate
#  This is because they have diverse plants for a same coordinate
# as our purpose is just to sample coordinates, we can summarise it.

merge_data$source_key_specie <- NA
merge_data$id_key = NA
id_key = 1
merge_data$source_key_specie[1] <- paste0(merge_data$scientificName[1],'_',id_key )
merge_data$id_key[1] <- 1

for(i in 2:nrow(merge_data))
{
  if(merge_data$scientificName[i] == merge_data$scientificName[i-1])
  {
    coord_comp <- merge_data$decimalLatitude[i]==merge_data$decimalLatitude[i-1]
    
    if(isTRUE(coord_comp))
    {
      merge_data$source_key_specie[i] <-merge_data$source_key_specie[i-1] 
      merge_data$id_key[i] =id_key
    }
    else
    {
      id_key <- id_key+1
      merge_data$source_key_specie[i] <-paste0(merge_data$scientificName[i],'_',  id_key)
      merge_data$id_key[i] =id_key
    }
    
  }
  else
  {
    id_key = 1
    merge_data$source_key_specie[i] <- paste0(merge_data$scientificName[i],'_',id_key )
    merge_data$id_key[i] =id_key
  }
  
}

# removing overlapping coordinates in the same species
merge_data <-  merge_data %>% ddply(.(scientificName,source_key_specie),
                                  summarise,
                                  scientificName = scientificName[which.min(id_key)],
                                  source_key_specie = source_key_specie[which.min(id_key)], # key 1 only
                                  decimalLatitude = decimalLatitude[which.min(id_key)],
                                  decimalLongitude = decimalLongitude[which.min(id_key)])


# now we run a second QC. At this step we use CoordinateCleaner to find weird geographic points

species = unique(merge_data $scientificName)

merge_data_clean = c()

for(i in 1:length(species))
{
  data <- merge_data %>% filter(scientificName %in% species[i])
  
  if(nrow(data) < 6) # less than 6 points it is not necessary. Subjective criteria. I use 6 because it was a minimum number of samples for one of those species
  {
    merge_data_clean  <- rbind(merge_data_clean ,data)
  }
  
  if(nrow(data) >=6) # if higher than 6 points I did the cleaning.
  {
    data2<- try(CoordinateCleaner::cc_outl(lon = 'decimalLongitude',
                                           lat = 'decimalLatitude',
                                           species = 'scientificName',
                                           method = 'mad',mltpl = 30,
                                           #    method = 'quantile',
                                           x = data))
    if(isTRUE(nrow(data2) >=1))
    {
      merge_data_clean  <- rbind(merge_data_clean ,data2)
    } 
    else
    {
      merge_data_clean  <- rbind(merge_data_clean ,data)
    }
  }
}

all.equal(species, unique(merge_data_clean$scientificName))


#'------------------------------------------------------------------------------------------------------------
# (4)Output ######
#'------------------------------------------------------------------------------------------------------------

merge_data_final <- 
  merge_data_clean %>% 
  merge(androMetadata,by.x='scientificName',by.y='names')
data.table::fwrite(file = paste0(dir.output,'/coordinates_clean.csv'),merge_data_final)


# ######################## main analysis end here; some useful stuff down there ####################################
# 
# ## MANUALLY CORRECTION OF COORDINATES #####
# # below there is some examples of how I modified the species names manually
# #  "Andropogon aridus" No info
# 
# missing_taxa[1]
# 
# data = try(rgbif::occ_search(scientificName = 'Andropogon',hasCoordinate = T)$data %>% as.data.frame())
# unique(data$scientificName)
# 
# #specie_i_need = name_to_save_it = 'Zea perennis'
# 
# 
# # Andropogon gerardii -> Andropogon gerardi
# specie_i_need = 'Andropogon gerardi'
# name_to_save_it = 'Andropogon gerardii'
# 
# path_species = paste0(getwd(),'/species_metadata_bien')
# data = try(BIEN::BIEN_occurrence_species(species = specie_i_need,cultivated = FALSE,only.geovalid = TRUE) )
# data$NamePanAnd = name_to_save_it
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# 
# path_species = paste0(getwd(),'/species_metadata')
# data = try(rgbif::occ_search(scientificName = specie_i_need,hasCoordinate = T)$data %>% as.data.frame())
# data$NamePanAnd = name_to_save_it
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# 
# # Bothriochloa woodrovii # found it at tropics.org
# # Bothriochloa woodrovii -> Andropogon woodrovii or Amphilophis woodrovii 
# specie_i_need = 'Amphilophis woodrovii'
# name_to_save_it = 'Bothriochloa woodrovii'
# 
# 
# # Dimeria mooneyi did not found it
# data = try(rgbif::occ_search(scientificName = 'Dimeria',hasCoordinate = T)$data %>% as.data.frame())
# unique(data$scientificName) %>% sort()
# 
# # from literature I found this. Or I asked Toby (don't remember now)
# # Ischaemum koenigii -> Ischaemum aristatum
# 
# data = try(rgbif::occ_search(scientificName = 'Ischaemum',hasCoordinate = T)$data %>% as.data.frame())
# unique(data $scientificName) %>% sort()
# 
# 
# 
# # Ischaemum koenigii -> Ischaemum aristatum
# 
# specie_i_need = 'Ischaemum aristatum'
# name_to_save_it = 'Ischaemum koenigii'
# 
# path_species = paste0(getwd(),'/species_metadata_bien')
# data = try(BIEN::BIEN_occurrence_species(species = specie_i_need,cultivated = FALSE,only.geovalid = TRUE) )
# data$NameID = data$scrubbed_species_binomial
# data$scrubbed_species_binomial <-name_to_save_it 
# 
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# 
# path_species = paste0(getwd(),'/species_metadata')
# data = try(rgbif::occ_search(scientificName = specie_i_need,hasCoordinate = T)$data %>% as.data.frame())
# data$NamePanAnd = name_to_save_it
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# # Iseilema venkateswarlui # tropics/ kellog
# 
# data = try(rgbif::occ_search(scientificName = 'Iseilema',hasCoordinate = T)$data %>% as.data.frame())
# unique(data $scientificName) %>% sort()
# 
# missing_taxa
# 
# # Bothriochloa grahamii
# 
# data = try(rgbif::occ_search(scientificName = 'Bothriochloa',hasCoordinate = T)$data %>% as.data.frame())
# unique(data$scientificName)
# 
# data = try(BIEN::BIEN_occurrence_species(species ='Bothriochloa grahamii',cultivated = FALSE,only.geovalid = TRUE) )
# 
# 
# # Saccharum strictum -> Erianthus strictus
# 
# specie_i_need = 'Saccharum baldwinii'
# name_to_save_it = 'Saccharum strictumi'
# 
# path_species = paste0(getwd(),'/species_metadata_bien')
# data = try(BIEN::BIEN_occurrence_species(species = specie_i_need,cultivated = FALSE,only.geovalid = TRUE) )
# data$NameID = data$scrubbed_species_binomial
# data$scrubbed_species_binomial <-name_to_save_it 
# 
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# specie_i_need = 'Erianthus strictus'
# 
# 
# path_species = paste0(getwd(),'/species_metadata')
# data = try(rgbif::occ_search(scientificName = specie_i_need,hasCoordinate = T)$data %>% as.data.frame())
# data$NamePanAnd = name_to_save_it
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# 
# 
# # Sorghum controversum
# 
# # Sorghum controversum -> Andropogon controversus # nop
# # Sorghum controversum -> Andropogon laxum # nop
# # Sorghum controversum -> Sorghum halepense # yay!
# 
# 
# 
# data = try(BIEN::BIEN_occurrence_species(species ='Sorghum halepense',cultivated = FALSE,only.geovalid = TRUE) )
# 
# data = try(rgbif::occ_search(scientificName = 'Sorghum halepense',hasCoordinate = T)$data %>% as.data.frame())
# 
# 
# 
# specie_i_need = 'Sorghum halepense'
# name_to_save_it = 'Saccharum strictumi'
# 
# path_species = paste0(getwd(),'/species_metadata_bien')
# data = try(BIEN::BIEN_occurrence_species(species = specie_i_need,cultivated = FALSE,only.geovalid = TRUE) )
# data$NameID = data$scrubbed_species_binomial
# data$scrubbed_species_binomial <-name_to_save_it 
# 
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# 
# 
# 
# path_species = paste0(getwd(),'/species_metadata')
# data = try(rgbif::occ_search(scientificName = specie_i_need,hasCoordinate = T)$data %>% as.data.frame())
# data$NamePanAnd = name_to_save_it
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# 
# 
# # Sorghum controversum -> !Sorghum halepense (L.) Pers???
# 
# 
# specie_i_need = 'Andropogon controversus'
# name_to_save_it = 'Sorghum controversum'
# 
# path_species = paste0(getwd(),'/species_metadata_bien')
# data = try(BIEN::BIEN_occurrence_species(species = specie_i_need,cultivated = FALSE,only.geovalid = TRUE) )
# data$NameID = data$scrubbed_species_binomial
# data$scrubbed_species_binomial <-name_to_save_it 
# 
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# 
# 
# # Themeda saxicola ??
# data = try(rgbif::occ_search(scientificName = 'Themeda',hasCoordinate = T)$data %>% as.data.frame())
# unique(data$scientificName)
# 
# path_species = paste0(getwd(),'/species_metadata')
# data = try(rgbif::occ_search(scientificName = specie_i_need,hasCoordinate = T)$data %>% as.data.frame())
# data$NamePanAnd = name_to_save_it
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# # not able to: Themeda saxicola, Tripsacum cundinamarce
# 
# data = try(rgbif::occ_search(scientificName = 'Tripsacum',hasCoordinate = T)$data %>% as.data.frame())
# data$scientificName %>% unique() %>% sort()
# 
# # Tripsacum zopilotense
# data = try(BIEN::BIEN_occurrence_species(species ='Tripsacum zopilotense',cultivated = FALSE,only.geovalid = TRUE) )
# 
# data = try(rgbif::occ_search(scientificName = 'Tripsacum zopilotense',hasCoordinate = T)$data %>% as.data.frame())
# 
# 
# 
# specie_i_need = 'Tripsacum zopilotense'
# name_to_save_it = 'Tripsacum zopolitense'
# 
# path_species = paste0(getwd(),'/species_metadata_bien')
# data = try(BIEN::BIEN_occurrence_species(species = specie_i_need,cultivated = FALSE,only.geovalid = TRUE) )
# data$NameID = data$scrubbed_species_binomial
# data$scrubbed_species_binomial <-name_to_save_it 
# 
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# 
# 
# 
# path_species = paste0(getwd(),'/species_metadata')
# data = try(rgbif::occ_search(scientificName = specie_i_need,hasCoordinate = T)$data %>% as.data.frame())
# data$NamePanAnd = name_to_save_it
# data.table::fwrite(file = paste0(path_species,'/',name_to_save_it ,'.csv'), data)
# 
# ### end ####
# 
# 
# #{
# # as a sack of exemplification, here was a situation where I had the specie but with a different name.
# # I was just checking it manually
# # PanAndTracker$species.name[c(2,26)] = c('Andropogon gerardii','Zea mays subsp. huehuetenangensis')
# # }
# 
# # then I made some venn diagrams (check #pandand_gwas channel)
# 
# # some species / taxa we already have some coordianates. So I used it:
# PanAndTracker <- 
#   read_excel("~/Library/CloudStorage/OneDrive-CornellUniversity/Buckler Lab/PanAnd/Envirotyping/Checklist_ENA-PlantSample_PanAndGenomesMetadata .xlsx", 
#              col_types = c("numeric", "text", "text", 
#                            "text", "text", "text", "text", "numeric", 
#                            "numeric", "text", "text", "text", 
#                            "text", "text", "text", "text", "text", 
#                            "numeric", "text", "text", "text", 
#                            "text"))
# 
# 
# head(PanAndTracker)
# names(PanAndTracker)[c(2,8,9)] = c('species.name','latitude','longitude')
# 
# dim(PanAndTracker)
# 
# PanAndTracker$source = 'PandAndTracker'
# 
# new=unique(androMetadata2023_08_24$speciesName)
# 
# new[!new %in% unique(bien_data_final$scientificName)] # comparting with bien
# 
# 
# we_have = c(unique(bien_data_final$scientificName),unique(gbif_data_final$scientificName))
# 
# missing_taxa = new[!new %in% we_have ] %>% sort()
# 
# # {
# # so let's say that after this i decided again to run bien and pull more missing entries.
# # this is a manual process, unfortunetly, and I'm only trying to give a big picture of it
# 
# # path_species = paste0(getwd(),'/species_metadata_bien')
# 
# #bien_data <- c()
# 
# 
# #for(j in 1:length(missing_taxa))
# #{
# #  
# #  cat(paste0('doing.... ',j,' out of ',length(missing_taxa),'\n'))
# #  data = try(BIEN::BIEN_occurrence_species(species = missing_taxa[j],cultivated = FALSE,only.geovalid = TRUE) )
# #  
# #  data_name = gsub(missing_taxa[j],pattern = ' ',replacement='_')
# #  
# #  # pulling only species ID, data key and coordinates
# #  if(nrow(data) == 0)
# #  {
# #    output <-  data.frame(date_collected=NA,scientificName=missing_taxa[j],decimalLatitude=NA,decimalLongitude=NA)
# #    data.table::fwrite(file = paste0(path_species,'/missing_species_v2.csv'), data.frame(species=missing_taxa[j]),append = T)
# #  }
# #  
# #  if(nrow(data) >0)
# #  {
# #   # saving all metadata
# #    data.table::fwrite(file = paste0(path_species,'/',data_name,'.csv'), data)
# #    
# #    output <-   data.frame(date_collected=data$date_collected,scientificName=missing_taxa[j],decimalLatitude=data$latitude,decimalLongitude=data$longitude)
# #    #data[,c('key','scientificName','decimalLatitude','decimalLongitude')]
# #  }
# #  bien_data <- rbind(bien_data,    output)
# #  
# #}
# 
# ##bien_data <-   bien_data %>% filter(!is.na(decimalLatitude)) %>% droplevels()
# 
# #missing_taxa[!missing_taxa %in% unique(bien_data$scientificName)]
# 
# #gbif_data <-     gbif_data %>% filter(!is.na(decimalLatitude)) %>% droplevels()
# 
# # missing_taxa[!missing_taxa %in% unique(gbif_data$scientificName)]
# #}
# 
# 
# 
# 
# 
