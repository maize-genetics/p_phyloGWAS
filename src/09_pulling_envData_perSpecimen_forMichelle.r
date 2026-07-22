#'##################################################################################################'
#'                                                                                                 #'
#' Project........ Environmental GWAS (eGWAS) for PandAnd                                          #'
#' Title.......... Environmental Data extraction for each occurrence record                        #'
#' Created at..... 11-06-2023                                                                      #'
#' Updated at..... 10-01-2025 (Sheng-Kai Hsu)                                                      #'
#' Author: G.Costa-Neto, germano.cneto<at>gmail.com, Sheng-Kai Hsu (sh2246<at>cornell.edu)         #'
#'                                                                                                 #'
#'##################################################################################################'

require(tidyverse)
require(plyr)
require(reshape2)

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")

#'------------------------------------------------------------------------------------------------------------
# (1) load geo data
#'------------------------------------------------------------------------------------------------------------
specimenCoordinate = read.table(file.path(PHYLOGWAS_ROOT, "data/combined_latlong_for_shengkai.txt"),header =T)
colnames(specimenCoordinate)[1] = c("sample")
head(specimenCoordinate)
dim(specimenCoordinate)

#'------------------------------------------------------------------------------------------------------------
# (2) extract env. features for each coordinate
#'------------------------------------------------------------------------------------------------------------
# creating a fake "environmental unit": species - sample
specimenCoordinate <-
  specimenCoordinate %>% 
  ddply(.(sample),mutate,envSampleName = paste0('env_',sample,'_',1:length(sample)))
dim(specimenCoordinate)
head(specimenCoordinate)

specimenCoordinateNoNA <- specimenCoordinate %>% na.omit()

########### Bioclimatic variables
# check src_generating_FAO_GAEZ.R to see how to generate enviromeDB::WC_Bioclimate since the package is broken
source('https://raw.githubusercontent.com/gcostaneto/envirotypeR/main/R/get_spatial_fun.R')

url = file.path(PHYLOGWAS_ROOT, 'output/envData/GIS_raster/WC_Bioclim.rds')
tmp = readRDS(url)
geographic_ranges  = 
  get_spatial( env.dataframe =specimenCoordinateNoNA,
                            lat = 'approxlat',
                            lng = 'approxlong',
                            env.id = 'envSampleName',
                            digital.raster = readRDS(url), # using a certain url)
  )

########### Elevation
url = '/workdir/sh2246/p_evolBNI/data/GIS_env_data/WorldClim_raw_2.5m_files/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif' 
geographic_ranges = 
  get_spatial( env.dataframe = geographic_ranges,
               lat = 'approxlat',
               lng = 'approxlong',
               env.id = 'envSampleName',
               name.feature = 'Elevation_m',
               digital.raster = terra::rast(url))#envirotypeR::SRTM_elevation #terra::rast(url), # using a certain url)

########### Global Hydrologic Soil Groups
url = '/workdir/sh2246/p_evolBNI/data/GIS_env_data/Global_Hydrologic_Soil_Group_1566/Global_Hydrologic_Soil_Group_1566/data/HYSOGs250m.tif'
geographic_ranges = 
  get_spatial( env.dataframe = geographic_ranges,
               lat = 'approxlat',
               lng = 'approxlong',
               env.id = 'envSampleName',
               name.feature = 'HYSOGs',
               digital.raster = terra::rast(url) #terra::rast(url), # using a certain url)
  ) 


########### FAO-GAEZ 
url = file.path(PHYLOGWAS_ROOT, 'output/envData/GIS_raster/GAEZ_AEZ.rds')
geographic_ranges = 
  get_spatial( env.dataframe = geographic_ranges,
               lat = 'approxlat',
               lng = 'approxlong',
               env.id = 'envSampleName',
               digital.raster = readRDS(url))#

########### Soil Temperature 
url = file.path(PHYLOGWAS_ROOT, 'output/envData/GIS_raster/TEMP_soil.rds')
geographic_ranges = 
  get_spatial( env.dataframe = geographic_ranges,
               lat = 'approxlat',
               lng = 'approxlong',
               env.id = 'envSampleName',
               digital.raster = readRDS(url))#


########### Soil Features from GSDE
# the file needs to be converted (reading as brick)
# and this is a too big .rds file.
# so let's pull each layer per time
urlList = list.files('/workdir/sh2246/p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files/',pattern = "*.nc",recursive = T,full.names = T)


for(i in 1:length(urlList))
{
    # this takes a long time in relation to the previous raster files. Don't worry.
  geographic_ranges = 
    get_spatial( env.dataframe = geographic_ranges,
                 lat = 'approxlat',
                 lng = 'approxlong',
                 env.id = 'envSampleName',
                 digital.raster = raster::brick(urlList[i]))#
}

# name the GSDE variablees
GSDE_varNames = limma::strsplit2(list.dirs("/workdir/sh2246/p_evolBNI/data/GIS_env_data/GSDE_raw_nc_files",recursive = F),"/")[,8]

colnames(geographic_ranges)[-c(1:71)] = paste(rep(GSDE_varNames,each = 4),c(5,15,30,200),sep = "_GSDE_")

# take the two layers only
rmIdx = grep('GSDE_30|_200',colnames(geographic_ranges))
geographic_ranges_filtered = geographic_ranges[,-rmIdx]

geographic_ranges_filtered[geographic_ranges_filtered==Inf] = NA
geographic_ranges_filtered[geographic_ranges_filtered==-Inf] = NA
geographic_ranges_filtered[geographic_ranges_filtered==-99] = NA
geographic_ranges_filtered[geographic_ranges_filtered==-999] = NA
noNAIdx = apply(geographic_ranges_filtered,2,function(x) !any(is.na(x)))
geographic_ranges_filtered[,noNAIdx][geographic_ranges_filtered[,noNAIdx]==156] = NA # for % data in GSDE, NA -> 156...



write.table(geographic_ranges_filtered,
            file.path(PHYLOGWAS_ROOT, "output/panand_specimen_envDataToMichelle_20251001.txt"),quote = F,sep = "\t")

