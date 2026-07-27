#!/bin/bash
export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"
cd "${PHYLOGWAS_ROOT}"

# env. data
# run 02A_metadataProcessing.ipynb first (species name -> metadata)
# species name -> coordinate
mkdir output/metadataFormalOut
Rscript src/08_pulling_geo_data.R data/spNameMetadata_20240819.txt output/metadataFormalOut

# species coordinates -> environmental features (GIS raster extraction)
Rscript src/09_pulling_envData.r output/metadataFormalOut/coordinates_clean.csv output/metadataFormalOut/formal_envData_20240820.txt

# environmental features -> envPCs (PCA) + derived pipeline outputs (envData_707Poaceae*, HyPhy adaptation lists, loading diagnostics)
Rscript src/S05_envPC_analysis.R
