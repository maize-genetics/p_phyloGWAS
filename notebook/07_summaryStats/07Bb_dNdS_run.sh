#!/bin/bash
# dN/dS to a reference outgroup (ASM1935983v1), per sequence, across all OGs.
# src/07B_getOmega2Ref.R was already correct/portable (CLI args) but never wired up to a
# driver in this repo - this reproduces data/fullSetOGs_240903.txt, the real combined table
# 08A_masterDataTableGeneration.ipynb already reads.
export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"
cd "${PHYLOGWAS_ROOT}"

REF=ASM1935983v1

mkdir -p output/omega2Ref
find output/CDSMSAPerOG_gs/ -type f -name '*.gs.fa' | \
  parallel -j 40 "Rscript src/07B_getOmega2Ref.R --input {} --output output/omega2Ref/{/.}.dnds.txt --ref ${REF}"

cat output/omega2Ref/*.dnds.txt > output/fullSetOGs.txt
