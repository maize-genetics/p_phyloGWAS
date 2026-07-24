#!/bin/bash
# Ancestral sequence reconstruction for the OGs filtered in 03B_OGFilter.ipynb
export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"

mkdir output/poaceaeHelixOGMSA_plusAnc
cat "${PHYLOGWAS_ROOT}/output/poaceaeHelixerOG_filtered_20250331.txt" | parallel -j 35 'Rscript src/S02_phangornAncestralSeqReconstruct_AA.R --input output/OrthoFinder/Results_Jun06/MultipleSequenceAlignments/{}.fa --output output/poaceaeHelixOGMSA_plusAnc/{}'

# extract the ancestral sequence
mkdir output/poaceaeHelixOG_ancSeq/
find output/poaceaeHelixOGMSA_plusAnc/ -type f| cut -d / -f 3| sort |parallel -j 40 'tail -n 2 output/poaceaeHelixOGMSA_plusAnc/{} > output/poaceaeHelixOG_ancSeq/{}'
# rename
mkdir output/poaceaeHelixOG_ancSeq_renamed
find output/poaceaeHelixOG_ancSeq/ -type f| cut -d / -f 3|sort|parallel -j 40 "cat output/poaceaeHelixOG_ancSeq/{} | sed 's/MRCA/{= s/.fa// =}/g' > output/poaceaeHelixOG_ancSeq_renamed/{}"
# add back
find output/poaceaeHelixOG_ancSeq_renamed/ -type f| sort | xargs cat > output/poaceaeHelixerOG_ancSeq_filtered.fa

# remove gap
cat output/poaceaeHelixerOG_ancSeq_filtered.fa |sed 's/-//g' > output/poaceaeHelixerOG_ancSeq_filtered_gapRemoved.fa

/programs/seqkit-0.15.0/seqkit seq -m 1 output/poaceaeHelixerOG_ancSeq_filtered_gapRemoved.fa > output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa
