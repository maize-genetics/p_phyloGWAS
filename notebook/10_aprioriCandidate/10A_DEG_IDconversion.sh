#!/bin/bash

# gather all the genomes -> cds.fa

# miniprot OG protein against cds.fa

find data/DEG_CDS_fasta/ -type f |cut -d '/' -f3 |sed 's/.fa.gz//g'|parallel -j 10 'miniprot --gff-only -t 4 data/DEG_CDS_fasta/{}.fa.gz output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa > data/DEG_mappingFiles/OGTo{}.gff'

find data/DEG_CDS_fasta/ -type f |cut -d '/' -f3 |sed 's/.fa.gz//g'|parallel -j 23 "grep \"Rank=1;\" data/DEG_mappingFiles/OGTo{}.gff|grep \"mRNA\" |awk '{print \$1\"\t\"\$9}'|sed 's/;Frameshift=*//g'|sed 's/;StopCodon=*//g'|sed 's/\t/;/g'|cut -d \";\" -f 1,6|sed 's/;Target=/\t/g' > data/DEG_mappingFiles/OG_mapping_{}.txt"