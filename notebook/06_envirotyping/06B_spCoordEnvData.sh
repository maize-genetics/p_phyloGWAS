#!/bin/bash
# env. data
# run 12_metadataProcessing.ipynb
# species name -> coordinate
mkdir output/metadataFormalOut
Rscript src/08_pulling_geo_data.R data/spNameMetadata_20240819.txt output/metadataFormalOut

#run 09_pulling_envData.r