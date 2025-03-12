#!/bin/bash

# 0) count missing taxa
find output/CDSPerOG_noATGFiltered/ -type f|sort |parallel -k -j 20 "grep '^-$' {} |wc -l >> output/numMissingTaxa.txt"
find output/CDSPerOG_noATGFiltered/ -type f | sort | parallel -j 20 --linebuffer "echo -n '{} '; grep '^-$' {} | wc -l" > output/numMissingTaxa.txt

# 1) MSA to dist matrix
mkdir output/distMatrix/
find output/MSAPerOG_noATGFiltered/ -type f | parallel -j 30 'Rscript src/01_MSA2Dist.R --input {} --output output/distMatrix/'
# Rscript src/01_MSA2Dist.R --input output/MSAPerOG_noATGFiltered --output output/distMatrix/ --thread 20

# 2) dist matrix to pseudogenization approximates
Rscript src/02_getDistPAV.R --input output/distMatrix --reference 473 --output output/

# 3) gene tree contruction
mkdir output/ref/
find output/MSAPerOG_noATGFiltered/ -type f |parallel -j 30 "grep 'Pavag' {} |sed 's/>//g'|shuf -n 1 > output/ref/{/.}.ref"

mkdir output/MSAPerOG_gs
find output/MSAPerOG_noATGFiltered/ -type f |parallel -j 30 "Rscript /workdir/sh2246/p_panAndOGASR/paml_pipeline/Ungap_MSA_SKH.R {} output/ref/{/.}.ref output/MSAPerOG_gs/{/.}.gs.fa"

mkdir output/geneTree/
cd output/geneTree/
find /workdir/sh2246/p_phyloGWAS/output/MSAPerOG_gs -name *.gs.fa -type f | parallel -j 30 '/workdir/sh2246/p_panAndOGASR/paml_pipeline/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s {} -# 1 -n {/.}.tree'
cd /workdir/sh2246/p_phyloGWAS

# 4) neutral species and genic neutral tree
export PATH=/programs/seqtk:$PATH
find output/MSAPerOG_gs -type f | parallel -j 30 "cat {} | seqtk seq -l0 > {= s/.gs.fa/.gs.l0.fa/  =}"

# generate dummy features
mkdir output/featureFiles/
find output/MSAPerOG_gs -name *.gs.l0.fa -type f| cut -d / -f 3|sed 's/.gs.l0.fa//g' |parallel -j 30 "head -n 2 output/MSAPerOG_gs/{}.gs.l0.fa |tail -n 1 |wc -c| awk '{print \"dummy\tCDS\t\"1\"\t\"\$0-1\"\t.\t+\t.\"\"\tID = gene1\"}' > output/featureFiles/{}.tmp"

find output/MSAPerOG_gs/ -name *.gs.l0.fa -type f| cut -d / -f 3|sed 's/.gs.l0.fa//g' | parallel -j 30 'paste -d "\t" output/ref/{}.ref output/featureFiles/{}.tmp > output/featureFiles/{}.gff3'

# get 4d MSA
mkdir output/MSAPerOG_4d
find output/MSAPerOG_gs -name *.gs.l0.fa -type f| cut -d / -f 3|sed 's/.gs.l0.fa//g' | parallel -j 30 'Rscript src/03_run_get4d.msa.R --feature output/featureFiles/{}.gff3 --input output/MSAPerOG_gs/{}.gs.l0.fa --output output/MSAPerOG_4d/{}.4d.fa' 

# construct RAxML tree
# mkdir output/localNeutralTree/
# cd output/localNeutralTree/
# find /workdir/sh2246/p_evolBNI/output/local_neutral_aln -name *.fa -type f| cut -d / -f 7|sed 's/.fa//g' | parallel -j 30 '/workdir/sh2246/p_panAndOGASR/paml_pipeline/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s /workdir/sh2246/p_evolBNI/output/local_neutral_aln/{}.fa -# 1 -n {}.tree'
# cd /workdir/sh2246/p_evolBNI/

mkdir output/geneNeutralTree/
cd output/geneNeutralTree/
find output/MSAPerOG_4d/ -name *.4d.fa -type f| cut -d / -f 3|sed 's/.4d.fa//g' | parallel -j 30 '/workdir/sh2246/p_panAndOGASR/paml_pipeline/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s /workdir/sh2246/p_phyloGWAS/output/MSAPerOG_4d/{}.4d.fa -o {} -# 1 -n {}.tree'
cd /workdir/sh2246/p_phyloGWAS/

mkdir output/neutralTree/
cd output/neutralTree/
/workdir/sh2246/p_panAndOGASR/paml_pipeline/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s /workdir/sh2246/p_phyloGWAS/output/aln_4d_random5k_20bp.fasta -# 1 -n aln_4d_random5k.fasta.tree
cd /workdir/sh2246/p_phyloGWAS/
