#!/bin/bash
# Ancestral sequence reconstruction for the OGs filtered in 03B_OGFilter.ipynb
export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"
cd ${PHYLOGWAS_ROOT}

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


# miniprot ID matching
miniprot --gff-only -t 30 data/Angiosperms353_orysaSequences.fasta output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa > output/angiosperm353ToOG.gff
grep "Rank=1" output/angiosperm353ToOG.gff|grep "mRNA" |awk '{print $1"\t"$9}'|sed 's/;Frameshift=*//g'|sed 's/;StopCodon=*//g'|sed 's/\t/;/g'|cut -d ";" -f 1,6|sed 's/;Target=/\t/g' > output/angiosperm353ToOG_mapping.txt

miniprot --gff-only -t 30 data/Zea_mays_v5_mrna.fa output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa > output/OGToZm_v2.gff
grep "Rank=1;" output/OGToZm_v2.gff|grep "mRNA" |awk '{print $1"\t"$9}'|sed 's/;Frameshift=*//g'|sed 's/;StopCodon=*//g'|sed 's/\t/;/g'|cut -d ";" -f 1,6|sed 's/;Target=/\t/g' > data/OGToZm_mapping_v2.txt

miniprot --gff-only -t 30 output/Pv-672_v3.0_cds.fa output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa > output/OGToPv_v2.gff
grep "Rank=1;" output/OGToPv_v2.gff|grep "mRNA" |awk '{print $1"\t"$9}'|sed 's/;Frameshift=*//g'|sed 's/;StopCodon=*//g'|sed 's/\t/;/g'|cut -d ";" -f 1,6|sed 's/;Target=/\t/g' > output/OGToPv_mapping_v2.txt


miniprot --gff-only -t 30 data/Araport11_seq_20220914_representative_gene_model.fa.gz output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa > output/OGToAt_v2.gff
grep "Rank=1" output/OGToAt_v2.gff|grep "mRNA" |awk '{print $1"\t"$9}'|sed 's/;Frameshift=*//g'|sed 's/;StopCodon=*//g'|sed 's/\t/;/g'|cut -d ";" -f 1,6|sed 's/;Target=/\t/g' > output/OGToAt_mapping_v2.txt

