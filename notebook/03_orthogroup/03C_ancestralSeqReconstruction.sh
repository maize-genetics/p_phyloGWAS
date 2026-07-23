#!/bin/bash
# Ancestral sequence reconstruction for the OGs filtered in 03B_OGFilter.ipynb
export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"

mkdir output/poaceaeHelixOGMSA_plusAnc_toAdd
cat "${PHYLOGWAS_ROOT}/output/poaceaeHelixerOG_filtered_20250331.txt" | parallel -j 35 'Rscript src/S02_phangornAncestralSeqReconstruct_AA.R --input output/OrthoFinder/Results_Jun06/MultipleSequenceAlignments/{}.fa --output output/poaceaeHelixOGMSA_plusAnc_toAdd/{}'

# extract the ancestral sequence
mkdir output/poaceaeHelixOG_ancSeq_toAdd/
find output/poaceaeHelixOGMSA_plusAnc_toAdd/ -type f| cut -d / -f 3| sort |parallel -j 40 'tail -n 2 output/poaceaeHelixOGMSA_plusAnc_toAdd/{} > output/poaceaeHelixOG_ancSeq_toAdd/{}'
# rename
mkdir output/poaceaeHelixOG_ancSeq_renamed_toAdd
find output/poaceaeHelixOG_ancSeq_toAdd/ -type f| cut -d / -f 3|sort|parallel -j 40 "cat output/poaceaeHelixOG_ancSeq_toAdd/{} | sed 's/MRCA/{= s/.fa// =}/g' > output/poaceaeHelixOG_ancSeq_renamed_toAdd/{}"
# add back
find output/poaceaeHelixOG_ancSeq_renamed_toAdd/ -type f| sort | xargs cat > output/poaceaeHelixerOG_ancSeq_filtered_toAdd.fa

# remove gap
cat output/poaceaeHelixerOG_ancSeq_filtered_toAdd.fa |sed 's/-//g' > output/poaceaeHelixerOG_ancSeq_filtered_toAdd_gapRemoved.fa

/programs/seqkit-0.15.0/seqkit seq -m 1 output/poaceaeHelixerOG_ancSeq_filtered_toAdd_gapRemoved.fa > output/poaceaeHelixerOG_ancSeq_filtered_toAdd_gapRemoved_noempty.fa

# NOTE (author-confirmed, not scripted): this "_toAdd" batch was manually merged into
# output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa (the file 03D_TableS5Generation
# actually reads) alongside an earlier base ancestral-sequence set. That merge step itself
# was never captured in a script -- a fresh re-run reproduces this "_toAdd" batch but would
# need that merge redone by hand to reproduce the final v2 file exactly.
