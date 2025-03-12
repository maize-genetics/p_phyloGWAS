rm(list = ls())

require(tidyverse)
require(plyr)
require(reshape2)
require(foreach)

require(BIEN) # to access BIEN data base
require(rgbif) 
library(parallel)
dir.data = "/workdir/sh2246/p_phyloGWAS/data/spNameMetadata_20240819.txt"
androMetadata  <- read.delim(dir.data,header = T)
ass_taxa = unique(androMetadata$names) %>% na.omit() %>% as.character() %>% sort()
grassSp = BIEN_taxonomy_family("Poaceae")

SpTw = BIEN_list_country("Taiwan",cultivated = F)
grassSpTw = intersect(grassSp[,7],SpTw[,2])

SpALL = mclapply(unique(BIEN_metadata_list_political_names()[,1]),function(x) BIEN_list_country(x,cultivated = F)[,2],mc.cores = 40)
names(SpALL) = unique(BIEN_metadata_list_political_names()[,1])
grassSpALL = mclapply(SpALL,function(x) intersect(grassSp[,7],unique(x)))


setdiff(setdiff(grassSpTw,grassSpALL),ass_taxa)
