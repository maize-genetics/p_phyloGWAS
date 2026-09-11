#'##################################################################################################'
#'                                                                                                 #'
#' Project........ Environmental GWAS (eGWAS) for PandAnd -- Hackathon                             #'
#' Title.......... Environmental Data extraction for each occurence record                         #'
#' Created at..... 11-06-2023                                                                      #'
#' Updated at..... 04-30-2024 (Sheng-Kai Hsu)                                                      #'
#' Author: G.Costa-Neto, germano.cneto<at>gmail.com                                                #'
#'                                                                                                 #'
#'##################################################################################################'

require(tidyverse)
require(plyr)
require(reshape2)

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")

#'------------------------------------------------------------------------------------------------------------
# (0) merge the sampleTracker name to assemblyID
#'------------------------------------------------------------------------------------------------------------
metadata = read.csv(file.path(PHYLOGWAS_ROOT, "data/Poaceae_metadata_2024.08.21.csv"),header = T)
specimenCoordinate = read.table(file.path(PHYLOGWAS_ROOT, "output/panAnd_sample_coordinate.tsv"),header =T)
specimenCoordinate = specimenCoordinate[!duplicated(specimenCoordinate$sample),]
specimenCoordinate_merged = merge(metadata,specimenCoordinate,by.x = "tracker_sample_name",by.y = "sample")

#'------------------------------------------------------------------------------------------------------------
# (1) load geo data 
#'------------------------------------------------------------------------------------------------------------
# data_clean = data.table::fread(file.path(PHYLOGWAS_ROOT, 'output/metadataFormalOut/coordinates_clean.csv'))
data_clean = specimenCoordinate_merged[,c(2,3,12,13)]
colnames(data_clean)[2:4] = c("latest_name","decimalLatitude","decimalLongitude")
head(data_clean)
dim(data_clean)

#'------------------------------------------------------------------------------------------------------------
# (2) extract env. features for each coordinate
#'------------------------------------------------------------------------------------------------------------
# creating a fake "environmental unit": species - sample
data_clean <-
  data_clean %>% 
  ddply(.(latest_name),mutate,envScientificName = paste0('env_',latest_name,'_',1:length(decimalLatitude)))
dim(data_clean)
head(data_clean)


data_clean <- data_clean %>% na.omit()

# table(data_clean$scientificName) %>% sort()
# table(data_clean$scientificName) %>% hist(breaks = seq(0,3500,10))

########### Bioclimatic variables
# check src_generating_FAO_GAEZ.R to see how to generate enviromeDB::WC_Bioclimate since the package is broken
source(file.path(PHYLOGWAS_ROOT, "src/get_spatial_fun.R"))

url = file.path(PHYLOGWAS_ROOT, 'output/envData/GIS_raster/WC_Bioclim.rds')
tmp = readRDS(url)
geographic_ranges_bien  = 
  get_spatial( env.dataframe =data_clean,
                            lat = 'decimalLatitude',
                            lng = 'decimalLongitude',
                            env.id = 'envScientificName',
                            digital.raster = readRDS(url), # using a certain url)
  )

########### Elevation
url = '/workdir/sh2246/p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif' 
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
               lat = 'decimalLatitude',
               lng = 'decimalLongitude',
               env.id = 'envScientificName',
               name.feature = 'Elevation_m',
               digital.raster = terra::rast(url))#envirotypeR::SRTM_elevation #terra::rast(url), # using a certain url)

########### Global Hydrologic Soil Groups
url = '/workdir/sh2246/p_evolBNI/data/GIS_env_data/Global_Hydrologic_Soil_Group_1566/Global_Hydrologic_Soil_Group_1566/data/HYSOGs250m.tif'
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
               lat = 'decimalLatitude',
               lng = 'decimalLongitude',
               env.id = 'envScientificName',
               name.feature = 'HYSOGs',
               digital.raster = terra::rast(url) #terra::rast(url), # using a certain url)
  ) 


########### FAO-GAEZ 
url = file.path(PHYLOGWAS_ROOT, 'output/envData/GIS_raster/GAEZ_AEZ.rds')
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
                            lat = 'decimalLatitude',
                            lng = 'decimalLongitude',
                            env.id = 'envScientificName',
                            digital.raster = readRDS(url))#

########### Soil Temperature 
url = file.path(PHYLOGWAS_ROOT, 'output/envData/GIS_raster/TEMP_soil.rds')
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
                            lat = 'decimalLatitude',
                            lng = 'decimalLongitude',
                            env.id = 'envScientificName',
                            digital.raster = readRDS(url))#


########### Soil Features from GSDE
# the file needs to be converted (reading as brick)
# and this is a too big .rds file.
# so let's pull each layer per time
urlList = list.files('/workdir/sh2246/p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files/',pattern = "*.nc",recursive = T,full.names = T)


for(i in 1:length(urlList))
{
    # this takes a long time in relation to the previous raster files. Don't worry.
  geographic_ranges_bien = 
    get_spatial( env.dataframe = geographic_ranges_bien,
                 lat = 'decimalLatitude',
                 lng = 'decimalLongitude',
                 env.id = 'envScientificName',#which.raster.number = 1,
                 digital.raster = raster::brick(urlList[i]))#
}

# name the GSDE variablees
GSDE_varNames = limma::strsplit2(list.dirs("/workdir/sh2246/p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files",recursive = F),"/")[,8]

colnames(geographic_ranges_bien)[-c(1:71)] = paste(rep(GSDE_varNames,each = 4),c(5,15,30,200),sep = "_GSDE_")

# take the first layer only
rmIdx = grep('GSDE_15|_30|_200',colnames(geographic_ranges_bien))
geographic_ranges_bien_filtered = geographic_ranges_bien[,-rmIdx]

# geographic_ranges_bien_filtered = geographic_ranges_bien_filtered[,c(6,4,5,7:104)]

geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==Inf] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-Inf] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-99] = NA
geographic_ranges_bien_filtered[geographic_ranges_bien_filtered==-999] = NA
noNAIdx = apply(geographic_ranges_bien_filtered,2,function(x) !any(is.na(x)))
geographic_ranges_bien_filtered[,noNAIdx][geographic_ranges_bien_filtered[,noNAIdx]==156] = NA # for % data in GSDE, NA -> 156...


# put NA on Weird results
# .ControlData <- function(x)
# {
#   if(isTRUE(x ==  Inf) ) x <- NA
#   if(isTRUE(is.nan(x)) ) x <- NA
#   if(isTRUE(x == -99)  ) x <- NA
#   if(isTRUE(x == -999)) x <- NA
#   if(isTRUE(x == -Inf) ) x <- NA
#   return(x)
# }


write.table(geographic_ranges_bien_filtered,
            file.path(PHYLOGWAS_ROOT, "output/panand_specimen_envData_20241126.txt"),quote = F,sep = "\t")


####comparison to perSpecies estimates####
rownames(geographic_ranges_bien_filtered) = geographic_ranges_bien_filtered$assemblyID
envData = readRDS(file.path(PHYLOGWAS_ROOT, 'output/ePC_20240827.rds'))
commonID = intersect(geographic_ranges_bien_filtered$assemblyID,rownames(envData$environmental.features))


r=c()
for (i in 1:95){
  tmpCol = colnames(envData$environmental.features)[2+(i-1)*3]
  r = c(r,cor(envData$environmental.features[commonID,tmpCol],
              geographic_ranges_bien_filtered[commonID,gsub("_quan50","",tmpCol)],use = "complete",method = "spearman"))
}
names(r) = gsub("_quan50","",colnames(envData$environmental.features)[seq(2,284,3)])
r[c("bio01_Annual_Mean_Temperature","bio12_Annual_Precipitation","eta",'PMEH1_GSDE_5')]
par(mfrow = c(2,2),mar =c(5,5,2,2))
plot(envData$environmental.features[commonID,2],geographic_ranges_bien_filtered[commonID,6],
     xlab = "Per Species Median",ylab= "Per Sample",main = "Bio01 - Annual Mean Temperature (rho = 0.742)",asp = 1)
plot(envData$environmental.features[commonID,2+3*3],geographic_ranges_bien_filtered[commonID,5+4],
     xlab = "Per Species Median",ylab= "Per Sample",main = "Bio01 - Annual Precipitation  (rho = 0.606)",asp = 1)
plot(envData$environmental.features[commonID,"eta_quan50"],geographic_ranges_bien_filtered[commonID,"eta"],
     xlab = "Per Species Median",ylab= "Per Sample",main = "Eta: evaportranspiration rate  (rho = 0.712)",asp = 1)
plot(envData$environmental.features[commonID,"PMEH1_GSDE_5_quan50"],geographic_ranges_bien_filtered[commonID,"PMEH1_GSDE_5"],
     xlab = "Per Species Median",ylab= "Per Sample",main = "Soil exchangeable P content  (rho = 0.242)",asp = 1)

par(mfrow = c(1,1),mar =c(12,5,1,1))
barplot(sort(r),las =2, cex.names=.6,ylab = "Spearman's correlation coefficient (rho)")

allSpeciesCoordinate = read.delim(file.path(PHYLOGWAS_ROOT, "output/metadataFormalOut/formal_envData_20240820.txt"),
                                        header = T)
par(mfrow = c(1,1),mar =c(12,5,1,1))
tmp = allSpeciesCoordinate[grep("Zea",allSpeciesCoordinate$latest_name),]
boxplot(tmp[,23]~tmp[,1],las =2,xlab = "",ylab = "elevation",cex.axis = .8)
