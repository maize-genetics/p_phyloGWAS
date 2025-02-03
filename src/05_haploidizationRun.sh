#!/bin/bash
# script to purge allelic assemblies to haploid
# author: Sheng-Kai Hsu
# date created: 2024.03.18
# date last edited: 2024.03.20
# mkdir output/haploidizationRun/
# mkdir output/haploidizationRun/outFasta
# usage bash src/05_haploidizationRun.sh <path_to_fasta>
# find /workdir/sh2246/p_phyloGWAS/data/final_8_19_22/ -name "*.fasta" -type f|parallel -j 11 "bash src/05_haploidizationRun.sh {}"

# sampleID=`basename $1 .fasta`
sampleID=`basename $1 .fa`

mkdir output/haploidizationRun/$sampleID
mkdir output/haploidizationRun/$sampleID/gffPerScaf
mkdir output/haploidizationRun/$sampleID/CDSPerScaf
mkdir output/haploidizationRun/$sampleID/tmp
mkdir output/haploidizationRun/$sampleID/log
cd output/haploidizationRun/$sampleID

/programs/seqkit-0.15.0/seqkit fx2tab -nl $1 |sort -k2,2nr |awk '{if($2>=1000000) print $1}' > goodScafName.txt 
# /programs/seqkit-0.15.0/seqkit fx2tab -nl $1 |sort -k2,2nr |awk '{if($2<1000000) print $1}' > badScafName.txt
# helixer annotation -> gff per scaffold -> CDS per scaffold
cat goodScafName.txt | parallel -j 2 "grep -w '^{}' /workdir/sh2246/p_phyloGWAS/data/annotations/\${sampleID}_helixer.gff > ./gffPerScaf/{}.gff3"

cat goodScafName.txt | parallel -j 2 "anchorwave gff2seq -i gffPerScaf/{}.gff3 -r \${1} -o CDSPerScaf/{}.cds.fa"

# iteration:
# find the smallest scaffold (scaffold number ordered by length already; good!)
/programs/seqtk/seqtk subseq $1 goodScafName.txt | /programs/seqkit-0.15.0/seqkit sort -l -r > tmp/remain.fa

for i in `tac goodScafName.txt | sed '$ d'`
do
	# i="scaf_29"
	echo iteration ${i}
	# get sequences
	cp tmp/remain.fa tmp/remain0.fa
	samtools faidx tmp/remain0.fa `grep -v ${i} goodScafName.txt` > tmp/remain.fa
	samtools faidx tmp/remain0.fa ${i} > tmp/${i}.fa
	
	# anchorwave anchor minimap (R1Q1) to the rest
	minimap2 --eqx -x splice -t 2 -k 12 -a -p 0.4 -N 20 tmp/remain.fa CDSPerScaf/${i}.cds.fa > tmp/remain.sam
	minimap2 --eqx -x splice -t 2 -k 12 -a -p 0.4 -N 20 tmp/${i}.fa CDSPerScaf/${i}.cds.fa > tmp/${i}.sam
	anchorwave proali -i gffPerScaf/${i}.gff3 \
	-as CDSPerScaf/${i}.cds.fa \
	-r tmp/${i}.fa \
	-a tmp/remain.sam \
	-ar tmp/${i}.sam \
	-s tmp/remain.fa \
	-n tmp/${i}.anchors \
	-R 1 -Q 1 -ns
	# -f output/haploidizationTest/tmp/${i}.f.maf #> output/haploidizationTest/log/${i}.log
	grep -v "#" tmp/${i}.anchors|tail -n +2 |grep -v "interanchor"|awk '{print $8}'> tmp/${i}.anchorsNames.txt

	# distance calculation -> decision
	samtools view -F 4 -F 256 tmp/remain.sam |grep -f tmp/${i}.anchorsNames.txt | awk '{
    cigar = $6;
    while (match(cigar, /[0-9]+[=X]/)) {
        val = substr(cigar, RSTART, RLENGTH - 1);
        if (substr(cigar, RSTART + RLENGTH - 1, 1) == "=") {
            total_matches += val;
        } else if (substr(cigar, RSTART + RLENGTH - 1, 1) == "X") {
            total_mismatches += val;
        }
        cigar = substr(cigar, RSTART + RLENGTH);
    }
}
END {
    if (total_matches + total_mismatches > 0) {
        overall_divergence = total_mismatches / (total_matches + total_mismatches);
        print overall_divergence;
    } else {
        print 1;
    }
}' >> tmp/dist.tab.txt

	# ~/MohamedsBioinformaticsTools-main/MAFDivergence/MAFDivergence -t 8 output/haploidizationTest/tmp/${i}.f.maf > output/haploidizationTest/tmp/${i}.dist.txt
	# awk '{sumBP+=$4; sumSNP+=$5} END {if(sumBP!=0) print sumSNP/sumBP}' output/haploidizationTest/tmp/${i}.dist.txt >> output/haploidizationTest/tmp/dist.tab.txt
done

