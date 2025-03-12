#!/bin/bash

ragtag.py scaffold -q 60 -f 100 -i 0.2 --remove-small -o output/scaffoldingTest/AN20T005 ../p_panAndOGASR/data/assemblies/final_8_19_22/Sb-JGI-v3.fasta data/testAssembly/AN20T005.final.contigs.fa 

ragtag.py scaffold -q 60 -f 100 -i 0.2 -o output/scaffoldingTest/AN20T005_v2 ../p_panAndOGASR/data/assemblies/final_8_19_22/Sb-JGI-v3.fasta data/testAssembly/AN20T005.final.contigs.fa 


mkdir output/scaffoldingTest/miniProt/
miniprot --gff-only -t 30 output/scaffoldingTest/AN20T005/ragtag.scaffold.fasta output/Pv-672_v3.0_peptide.fa > output/scaffoldingTest/miniProt/AN20T005.scaffolded.gff

miniprot --gff-only -t 30 output/scaffoldingTest/AN20T005_v2/ragtag.scaffold.fasta output/Pv-672_v3.0_peptide.fa > output/scaffoldingTest/miniProt/AN20T005.scaffolded_v2.gff


for i in `cat data/testSamples.txt |cut -f1`
do
scp cbsublfs1://data1/users/ajs692/panand_data/megahit_shortread_assemblies/final_contigs/${i}* data/testAssembly/
done

parallel --link -j 10 "ragtag.py scaffold -q 60 -f 100 -i 0.2 --remove-small -o output/scaffoldingTest/{1} data/final_8_19_22/{2} data/testAssembly/{1}.final.contigs.fa" ::: `cat data/testSamples.txt |cut -f1` ::: `cat data/testSamples.txt |cut -f2`

ragtag.py scaffold -q 60 -f 100 -i 0.2 --remove-small -o output/scaffoldingTest/AN21TSTL0252 data/final_8_19_22/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta data/testAssembly/AN21TSTL0252.final.contigs.fa 

#fai step takes longer than expected

cat data/testSamples.txt |cut -f1 |parallel -j 10 "miniprot --gff-only -t 2 data/testAssembly/{}.final.contigs.fa output/Pv-672_v3.0_peptide.fa > output/scaffoldingTest/miniProt/{}.gff"

cat data/testSamples.txt |cut -f1 |parallel -j 10 "miniprot --gff-only -t 2 output/scaffoldingTest/{}/ragtag.scaffold.fasta output/Pv-672_v3.0_peptide.fa > output/scaffoldingTest/miniProt/{}.scaffolded.gff"



#ploidy test

mkdir output/scaffoldingTest/ploidy/

find data/contigged_PanandAssemblies/ -type f|cut -d "/" -f 3|sed 's/_710bpContigs.fasta//g'|parallel -j 5 "ragtag.py scaffold -q 60 -f 100 -i 0.2 --remove-small -o output/scaffoldingTest/ploidy/{} data/final_8_19_22/{}.fasta data/contigged_PanandAssemblies/{}_710bpContigs.fasta"


miniprot --gff-only -t 30 output/scaffoldingTest/ploidy/Cs-KelloggPI219580-DRAFT-PanAnd-1.0/ragtag.scaffold.fasta output/Pv-672_v3.0_peptide.fa > output/scaffoldingTest/ploidy/Cs-KelloggPI219580-DRAFT-PanAnd-1.0/scaffolded.gff
