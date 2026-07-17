# Rscript to plot Figure S1
# author: Sheng-Kai Hsu
# date created: 2024.04.30
# date last edit: 2024.04.30

library(raster)
library(terra)
require(tidyverse)
require(plyr)
require(reshape2)

bien_data_clean = data.table::fread('/workdir/sh2246/p_phyloGWAS/output/envData/speciesRange/bien_coordinates_clean_2023.12.04.csv')
head(bien_data_clean)
dim(bien_data_clean)

# creating a fake "environmental unit": species - sample
bien_data_clean <-
  bien_data_clean %>% 
  ddply(.(scientificName),mutate,envScientificName = paste0('env_',scientificName,'_',1:length(decimalLatitude)))
bien_data_clean <- bien_data_clean %>% na.omit()


test = subset(brick("/workdir/sh2246/p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif"),1)
set.seed(1234)
coordinates<-bien_data_clean[grep("Andropogon gerardii",bien_data_clean[,3]),]
coordinates <- coordinates %>%
  group_by(scientificName) %>%
  sample_n(size = min(500, n()), replace = FALSE)


v <- vect(coordinates, c("decimalLongitude", "decimalLatitude"), crs="+proj=longlat")
vv <- project(v, crs(test))

png("/workdir/sh2246/p_phyloGWAS/output/figure/suppFig/suppFig_envPCPipeline_a1.png",
    width = 6,height = 3,units = "cm",res = 600,pointsize = 6)
par(mar =c (2,2,1,1))
plot(test,col = colorRampPalette(c("blue", "lightblue", "yellow","red"))(255),
     ylim = c(-60,70),cex.axis = .75,asp = NA,legend = F,mgp = c(3,.5,0))
points(vv, cex=.1,col=alpha("forestgreen",.5),pch = 1)
dev.off()
