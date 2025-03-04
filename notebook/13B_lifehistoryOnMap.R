require(tidyverse)
require(plyr)
require(reshape2)
library(raster)
library(terra)

#'------------------------------------------------------------------------------------------------------------
# (0) merge the sampleTracker name to assemblyID
#'------------------------------------------------------------------------------------------------------------
metadata = read.csv("/workdir/sh2246/p_phyloGWAS/data/Poaceae_metadata_2024.08.21.csv",header = T)
specimenCoordinate = read.table("/workdir/sh2246/p_phyloGWAS/output/panAnd_sample_coordinate.tsv",header =T)
specimenCoordinate = specimenCoordinate[!duplicated(specimenCoordinate$sample),]
specimenCoordinate_merged = merge(metadata,specimenCoordinate,by.x = "tracker_sample_name",by.y = "sample")


test = brick("/workdir/sh2246/p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_bio/wc2.1_2.5m_bio_18.tif")

tempannualCoord = na.omit(specimenCoordinate_merged[specimenCoordinate_merged$lifeHistory=="annual"&abs(specimenCoordinate_merged$Lat)>23.5,c(3,12,13)])
tropannualCoord = na.omit(specimenCoordinate_merged[specimenCoordinate_merged$lifeHistory=="annual"&abs(specimenCoordinate_merged$Lat)<23.5,c(3,12,13)])
tempperennialCoord = na.omit(specimenCoordinate_merged[specimenCoordinate_merged$lifeHistory=="perennial"&abs(specimenCoordinate_merged$Lat)>23.5,c(3,12,13)])
tropperennialCoord = na.omit(specimenCoordinate_merged[specimenCoordinate_merged$lifeHistory=="perennial"&abs(specimenCoordinate_merged$Lat)<23.5,c(3,12,13)])

colnames(tempannualCoord)[c(3,2)] = c("decimalLongitude", "decimalLatitude")
colnames(tempperennialCoord)[c(3,2)] = c("decimalLongitude", "decimalLatitude")
colnames(tropannualCoord)[c(3,2)] = c("decimalLongitude", "decimalLatitude")
colnames(tropperennialCoord)[c(3,2)] = c("decimalLongitude", "decimalLatitude")

v <- vect(tempannualCoord, c("decimalLongitude", "decimalLatitude"), crs="+proj=longlat")
v2 <- vect(tempperennialCoord, c("decimalLongitude", "decimalLatitude"), crs="+proj=longlat")
v3 <- vect(tropannualCoord, c("decimalLongitude", "decimalLatitude"), crs="+proj=longlat")
v4 <- vect(tropperennialCoord, c("decimalLongitude", "decimalLatitude"), crs="+proj=longlat")
vv <- project(v, crs(test))
vv2 <- project(v2, crs(test))
vv3 <- project(v3, crs(test))
vv4 <- project(v4, crs(test))

# png("/workdir/sh2246/p_evolBNI/output/publication_plot/Figure3J.png",width = 17.4,height = 8.7,units = "cm",res = 600,pointsize = 8)
par(mar =c (3,3,1,1))
plot(log10(test),col = colorRampPalette(c("peru", "lightyellow", "forestgreen"))(255))
points(vv, cex=.9,col="darkblue",pch = 1)
points(vv2, cex=.9, col="darkblue",pch = 19)
points(vv3, cex=.9, col="darkred",pch = 1)
points(vv4, cex=.9, col="darkred",pch = 19)
# dev.off()
