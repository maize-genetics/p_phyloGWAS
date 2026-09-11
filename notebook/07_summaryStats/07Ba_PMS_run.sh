#!/bin/bash
# Premature stop codon presence/absence per sequence, across all OGs.
# (replaces the SCINET-only 07Ba_PMS_SCINET.sh - output/orthofinderProteinMSAs_fullset_20250710/
#  already has the full, current OG set locally, so this runs as a plain local GNU-parallel job
#  matching the convention established in stage 03; see archived/07Ba_PMS_SCINET.sh for the
#  original SCINET/SLURM-array version)
export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"
cd "${PHYLOGWAS_ROOT}"

mkdir -p output/prematureStop
find output/orthofinderProteinMSAs_fullset_20250710/ -type f -name '*.fa' | \
  parallel -j 40 "python src/11_find_premature_stops.py {} output/prematureStop/{/.}_prematureStops.txt"

awk 'FNR!=1' output/prematureStop/*.txt > output/combined_PMS_20250421.txt
