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


#'------------------------------------------------------------------------------------------------------------
# (1) load geo data 
#'------------------------------------------------------------------------------------------------------------
data_clean = data.table::fread('/workdir/sh2246/p_phyloGWAS/output/metadataFormalOut/coordinates_clean.csv')
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
source('https://raw.githubusercontent.com/gcostaneto/envirotypeR/main/R/get_spatial_fun.R')

url = '/workdir/sh2246/p_phyloGWAS/output/envData/GIS_raster/WC_Bioclim.rds'
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
url = '/workdir/sh2246/p_phyloGWAS/output/envData/GIS_raster/GAEZ_AEZ.rds'
geographic_ranges_bien = 
  get_spatial( env.dataframe = geographic_ranges_bien,
                            lat = 'decimalLatitude',
                            lng = 'decimalLongitude',
                            env.id = 'envScientificName',
                            digital.raster = readRDS(url))#

########### Soil Temperature 
url = '/workdir/sh2246/p_phyloGWAS/output/envData/GIS_raster/TEMP_soil.rds'
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

geographic_ranges_bien_filtered = geographic_ranges_bien_filtered[,c(6,4,5,7:104)]

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
            "/workdir/sh2246/p_phyloGWAS/output/metadataFormalOut/formal_envData_20240820.txt",quote = F,sep = "\t")


