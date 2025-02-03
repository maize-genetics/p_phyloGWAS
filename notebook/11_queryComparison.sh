#!/bin/bash
# consensus query test
# miniProt different set of concensus to the long read genomes
# genomes: /data1/users/mcs368/panand_genomes/genomes

mkdir output/queryComparison
mkdir output/queryComparison/miniProt_Pv/
mkdir output/queryComparison/miniProt_helixerOG/
mkdir output/queryComparison/miniProt_poaceaeHelixerOG/
mkdir output/queryComparison/miniProt_eggNOG/

#Pv protein
for i in /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta
do 
f="$(basename -s .fasta $i)"
# echo $f
miniprot --gff-only -t 30 $i output/Pv-672_v3.0_peptide.fa > output/queryComparison/miniProt_Pv/${f}.gff
done

#HelixerOG
for i in /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta
do 
f="$(basename -s .fasta $i)"
# echo $f
miniprot --gff-only -t 30 $i output/HelixerOG_ancSeq_gapRemoved.fa > output/queryComparison/miniProt_helixerOG/${f}.gff
done

#HelixerOG - poaceae
for i in /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta
do 
f="$(basename -s .fasta $i)"
# echo $f
miniprot --gff-only -t 30 $i output/poaceaeHelixerOG_ancSeq_gapRemoved.fa > output/queryComparison/miniProt_poaceaeHelixerOG/${f}.gff
done


# eggNOG
find data/consensus/ -type f| sort | xargs cat > output/eggNOGv5_consensus.fa

for i in /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta
do 
f="$(basename -s .fasta $i)"
# echo $f
miniprot --gff-only -t 30 $i output/eggNOGv5_consensus.fa > output/queryComparison/miniProt_eggNOG/${f}.gff
done

## extract CDS
mkdir output/queryComparison/CDSPerTaxa_Pv
for i in /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta
do 
f="$(basename -s .fasta $i)"
# echo $f
gffread -g $i -x output/queryComparison/CDSPerTaxa_Pv/$f.fa output/queryComparison/miniProt_Pv/$f.gff
done

mkdir output/queryComparison/CDSPerTaxa_helixerOG
for i in /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta
do 
f="$(basename -s .fasta $i)"
# echo $f
gffread -g $i -x output/queryComparison/CDSPerTaxa_helixerOG/$f.fa output/queryComparison/miniProt_helixerOG/$f.gff
done

mkdir output/queryComparison/CDSPerTaxa_poaceaeHelixerOG
for i in /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta
do 
f="$(basename -s .fasta $i)"
# echo $f
gffread -g $i -x output/queryComparison/CDSPerTaxa_poaceaeHelixerOG/$f.fa output/queryComparison/miniProt_poaceaeHelixerOG/$f.gff
done

mkdir output/queryComparison/CDSPerTaxa_eggNOG
for i in /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta
do 
f="$(basename -s .fasta $i)"
# echo $f
gffread -g $i -x output/queryComparison/CDSPerTaxa_eggNOG/$f.fa output/queryComparison/miniProt_eggNOG/$f.gff
done

# filter miniProt gff

mkdir output/queryComparison/miniProtFilter_outTab_Pv
mkdir output/queryComparison/miniProtFilter_outFasta_Pv
mkdir output/queryComparison/miniProtFilter_outGFF_Pv
find /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta -type f |cut -d / -f 7|sed 's/.fasta//g'|parallel -j 20 "Rscript src/06_miniProtFilter.R --gffDir output/queryComparison/miniProt_Pv/{}.gff --queryFa output/Pv-672_v3.0_peptide.fa --inFa output/queryComparison/CDSPerTaxa_Pv/{}.fa --method p --cov 0 --outTab output/queryComparison/miniProtFilter_outTab_Pv/{}.txt --outFa output/queryComparison/miniProtFilter_outFasta_Pv/{}.fa --outGff output/queryComparison/miniProtFilter_outGFF_Pv/{}.gff"

mkdir output/queryComparison/miniProtFilter_outTab_helixerOG
mkdir output/queryComparison/miniProtFilter_outFasta_helixerOG
mkdir output/queryComparison/miniProtFilter_outGFF_helixerOG
find /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta -type f |cut -d / -f 7|sed 's/.fasta//g'|parallel -j 20 "Rscript src/06_miniProtFilter.R --gffDir output/queryComparison/miniProt_helixerOG/{}.gff --queryFa output/HelixerOG_ancSeq_gapRemoved.fa --inFa output/queryComparison/CDSPerTaxa_helixerOG/{}.fa --method p --cov 0 --outTab output/queryComparison/miniProtFilter_outTab_helixerOG/{}.txt --outFa output/queryComparison/miniProtFilter_outFasta_helixerOG/{}.fa --outGff output/queryComparison/miniProtFilter_outGFF_helixerOG/{}.gff"

mkdir output/queryComparison/miniProtFilter_outTab_eggNOG
mkdir output/queryComparison/miniProtFilter_outFasta_eggNOG
mkdir output/queryComparison/miniProtFilter_outGFF_eggNOG
find /workdir/sh2246/p_phyloGWAS/data/genomes/*.fasta -type f |cut -d / -f 7|sed 's/.fasta//g'|parallel -j 20 "Rscript src/06_miniProtFilter.R --gffDir output/queryComparison/miniProt_eggNOG/{}.gff --queryFa output/eggNOGv5_consensus.fa --inFa output/queryComparison/CDSPerTaxa_eggNOG/{}.fa --method p --cov 0 --outTab output/queryComparison/miniProtFilter_outTab_eggNOG/{}.txt --outFa output/queryComparison/miniProtFilter_outFasta_eggNOG/{}.fa --outGff output/queryComparison/miniProtFilter_outGFF_eggNOG/{}.gff"

# MAFFT MSA generation
Rscript src/07_seqRearrange.R output/queryComparison/miniProtFilter_outFasta_Pv 10 output/Pv-672_v3.0_peptide.fa output/queryComparison/CDSPerOG_Pv/ &
Rscript src/07_seqRearrange.R output/queryComparison/miniProtFilter_outFasta_helixerOG 10 output/HelixerOG_ancSeq_gapRemoved.fa output/queryComparison/CDSPerOG_helixerOG/ &
Rscript src/07_seqRearrange.R output/queryComparison/miniProtFilter_outFasta_eggNOG 10 output/eggNOGv5_consensus.fa output/queryComparison/CDSPerOG_eggNOG/


mkdir output/queryComparison/CDSMSAPerOG_Pv
mkdir output/queryComparison/CDSMSAPerOG_helixerOG
mkdir output/queryComparison/CDSMSAPerOG_eggNOG
find output/queryComparison/CDSPerOG_Pv -type f |parallel -j 20 "/programs/mafft/bin/mafft --ep 0 --genafpair --maxiterate 1000 {} > output/queryComparison/CDSMSAPerOG_Pv/{/.}.fa" #mem intensive
find output/queryComparison/CDSPerOG_helixerOG -type f |parallel -j 20 "/programs/mafft/bin/mafft --ep 0 --genafpair --maxiterate 1000 {} > output/queryComparison/CDSMSAPerOG_helixerOG/{/.}.fa" #mem intensive
find output/queryComparison/CDSPerOG_eggNOG -type f |parallel -j 20 "/programs/mafft/bin/mafft --ep 0 --genafpair --maxiterate 1000 {} > output/queryComparison/CDSMSAPerOG_eggNOG/{/.}.fa" #mem intensive

#diagnosis of variable 'ATG'
mkdir output/queryComparison/ATGProfile_Pv
mkdir output/queryComparison/ATGProfile_helixerOG
mkdir output/queryComparison/ATGProfile_eggNOG
find /workdir/sh2246/p_phyloGWAS/output/queryComparison/CDSMSAPerOG_Pv/ -type f| parallel -j 20 "python ./src/04_ATGdiagnosis.py {} output/queryComparison/ATGProfile_Pv/{/.}.atgProfile.txt"
find /workdir/sh2246/p_phyloGWAS/output/queryComparison/CDSMSAPerOG_helixerOG/ -type f| parallel -j 20 "python ./src/04_ATGdiagnosis.py {} output/queryComparison/ATGProfile_helixerOG/{/.}.atgProfile.txt"
find /workdir/sh2246/p_phyloGWAS/output/queryComparison/CDSMSAPerOG_eggNOG/ -type f| parallel -j 20 "python ./src/04_ATGdiagnosis.py {} output/queryComparison/ATGProfile_eggNOG/{/.}.atgProfile.txt"

