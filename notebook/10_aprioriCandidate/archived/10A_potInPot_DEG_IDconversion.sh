#!/bin/bash
# Pot-in-pot dataset ID conversion, extracted from 10A_DEG_IDconversion.sh. Author had
# already commented this block out there (a Helixer-annotation-based variant of the same
# miniprot ID-conversion step, for a different assembly set under
# output/candidateGenes_remob/assemblies_pp/) - moved here rather than left disabled in
# place, matching the archived-code convention. Kept commented out exactly as found.

# #pot-in-pot dataset
# mkdir data/DEG_helixer_CDS_fasta
# find output/candidateGenes_remob/assemblies_pp/ -type f -name '*.fasta'|cut -d '/' -f4|sed 's/.fasta//g'|parallel -j 13 'gffread output/candidateGenes_remob/assemblies_pp/{}_hxr.gff -g output/candidateGenes_remob/assemblies_pp/{}.fasta -x data/DEG_helixer_CDS_fasta/{}_hxrCDS.fa'

# mkdir data/DEG_mappingFiles_helixer
# find data/DEG_helixer_CDS_fasta/ -type f |cut -d '/' -f3 |sed 's/.fa//g'|parallel -j 13 'miniprot --gff-only -t 4 data/DEG_helixer_CDS_fasta/{}.fa output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa > data/DEG_mappingFiles_helixer/OGTo{}.gff'

# find data/DEG_helixer_CDS_fasta/ -type f |cut -d '/' -f3 |sed 's/.fa//g'|parallel -j 13 "grep \"Rank=1;\" data/DEG_mappingFiles_helixer/OGTo{}.gff|grep \"mRNA\" |awk '{print \$1\"\t\"\$9}'|sed 's/;Frameshift=*//g'|sed 's/;StopCodon=*//g'|sed 's/\t/;/g'|cut -d \";\" -f 1,6|sed 's/;Target=/\t/g' > data/DEG_mappingFiles_helixer/OG_mapping_{}.txt"
