#!/bin/bash

#miniProt results to mafft alignment

# filter for primary alignment only
mkdir data/miniprotPrimaryAlignmentVsPaspalum/
find data/miniprotAlignmentVsPaspalum_unfiltered/ -name "*.gff" -type f |parallel -j 30 "cat {}|grep 'Rank=1' > data/miniprotPrimaryAlignmentVsPaspalum/{/.}_primary.gff"

# gffread to get CDS
# to get the assemblies...
mkdir output/CDSPerTaxa
find /workdir/coh22/androMotifs/genomes/fastas/ -name "*.fa" -type f | parallel -j 20 "gffread -g {} -x output/CDSPerTaxa/{/.}.fa data/miniprotAlignmentVsPaspalum_unfiltered/{/.}.gff"

# filter for 70% CDS coverage
mkdir output/miniProtFilter_outTab
mkdir output/miniProtFilter_outFasta
mkdir output/miniProtFilter_outGFF
find /workdir/coh22/androMotifs/genomes/fastas/ -name "*.fa" -type f | parallel -j 20 "Rscript src/06_miniProtFilter.R --gffDir data/miniprotAlignmentVsPaspalum_unfiltered/{/.}.gff --queryFa output/Pv-672_v3.0_peptide.fa --inFa output/CDSPerTaxa/{/.}.fa --method p --cov 0.7 --outTab output/miniProtFilter_outTab/{/.}.txt --outFa output/miniProtFilter_outFasta/{/.}.fa --outGff output/miniProtFilter_outGFF/{/.}.gff"

# run p_phyloGWAS/notebook/05B_seqFilesRearrange

# move to SCINET and perform MAFFT alignment there
# see 03B_MAFFT_SCINET.sh

# mkdir output/CDSMSAPerOG
# find output/CDSPerOG -type f |parallel -j 10 "/programs/mafft/bin/mafft --ep 0 --genafpair --maxiterate 1000 {} > output/CDSMSAPerOG/{/.}.fa" #mem intensive
