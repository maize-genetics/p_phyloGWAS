#!/bin/bash

cd /workdir/sh2246/p_phyloGWAS/


## step 1: MSA cleaning
mkdir output/CDSMSAPerOG_HyPhy_20250203
find output/CDSMSAPerOG_gs/ -name "*.l0.fa" -type f| cut -d / -f 3|sed 's/.gs.l0.fa//g' | parallel -j 40 " /programs/hyphy-2.5.49/bin/hyphy cln Universal output/CDSMSAPerOG_gs/{}.gs.l0.fa No/No output/CDSMSAPerOG_HyPhy_20250203/{}.fa" > output/MSAcleaning.log 2>&1

## step 2: get target tip names (based on data/annual_assemblies_20250203.txt)
mkdir output/target_PAML_20250203
find output/CDSMSAPerOG_HyPhy_20250203/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' |parallel -j 40 "grep -f data/annual_assemblies_20250203.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/target_PAML_20250203/{}.target"

## step 3: RAxML gene tree generation
mkdir cd output/geneTree_allOGs_20250203/
find /workdir/sh2246/p_phyloGWAS/output/CDSMSAPerOG_HyPhy_20250203/ -name '*.fa' -type f| cut -d / -f 7|sed 's/.fa//g' | parallel -j 35 'src/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s /workdir/sh2246/p_phyloGWAS/output/CDSMSAPerOG_HyPhy_20250203/{}.fa -# 1 -w /workdir/sh2246/p_phyloGWAS/output/geneTree_allOGs_20250203/ -n {}.tree'

## step 4: tree labeling
find output/geneTree_allOGs/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' | parallel -j 10 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/target_PAML_20250203/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/geneTree_allOGs_20250203/RAxML_Labeled_bestTree.{}_Relax.tree"

## step 5: HyPhy RELAX test
mkdir output/HyPhyResult
find output/geneTree_allOGs/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' | parallel -j 10 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_20250203/{}.fa --tree output/geneTree_allOGs_20250203/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/{}.RELAX.json"

## to extract test statistics from json files into final output
