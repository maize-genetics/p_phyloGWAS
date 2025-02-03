#!/bin/bash

cd /workdir/sh2246/p_panAndOGASR/


# get target tip names
mkdir output/target_PAML
find output/geneTree_allOGs/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' |parallel -j 20 "grep -f data/annual_assemblies.txt output/geneTree_allOGs/{}.fa|sed 's/>//g' > output/target_PAML/{}.target"

## test
mkdir output/paml_testRun
find output/geneTree_allOGs/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' |head | parallel -j 6 "cp output/geneTree_allOGs/{}.fa output/paml_testRun/"
find output/geneTree_allOGs/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' |head | parallel --dryrun -j 6 "bash /workdir/sh2246/p_panAndOGASR/paml_pipeline/Run_PAML_scripts_SKH_v4.sh output/paml_testRun/{}.fa 'ASM1935983v1' output/target_PAML/{}.target 'ASM1935983v1'"

output/paml_testRun/OG0000928_cleaned.fa
output/paml_testRun/RAxML_Labeled_bestTree.OG0000928.tree

## HyPhy
/programs/hyphy-2.5.49/bin/hyphy busted --alignment output/paml_testRun/OG0000928_cleaned.fa --tree output/paml_testRun/RAxML_Labeled_bestTree.OG0000928.tree --branches Foreground

/programs/hyphy-2.5.49/bin/hyphy absrel --alignment output/paml_testRun/OG0000928_cleaned.fa --tree output/paml_testRun/RAxML_Labeled_bestTree.OG0000928.tree --branches Foreground

/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/paml_testRun/OG0000928_cleaned.fa --tree output/paml_testRun/RAxML_Labeled_bestTree.OG0000928_Relax.tree --test Foreground --reference Background --models Minimal
/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/paml_testRun/subset.fa --tree output/paml_testRun/RAxML_Labeled_bestTree.OG0000928_Relax.tree --test Foreground --reference Background --models Minimal

cat output/PAML_testRun_geneID.txt | parallel -j 35 "bash paml_pipeline/Run_PAML_scripts_SKH_v3.sh output/paml_testRun3/{}.aln.gs.fw.l0.fa output/ref/{}.ref output/target/{}.target output/root/{}.root"

mkdir output/paml_testRun2
cat output/PAML_testRun_geneID.txt | parallel -j 35 "cp output/aln_with_paspalumCDS_forward/{}.aln.gs.fw.l0.fa output/paml_testRun2/"
cat output/PAML_testRun_geneID.txt | parallel -j 35 "bash paml_pipeline/Run_PAML_scripts_SKH_v2.sh output/paml_testRun2/{}.aln.gs.fw.l0.fa output/ref/{}.ref output/target/{}.target output/root/{}.root"

# transporter test
find output/transporter_candidate_testrun1/ -type f|cut -d / -f 3 |awk '{print substr($1,1,14)}'sort|uniq| parallel -j 2 "bash paml_pipeline/Run_PAML_scripts_SKH_v3.sh output/transporter_candidate_testrun1/{}.aln.gs.fw.l0.fa output/ref/{}.ref output/target/{}.target output/root/{}.root"

# reproduction test
find output/reproduction_candidate_testrun/ -type f|cut -d / -f 3 |awk '{print substr($1,1,14)}'|sort|uniq| parallel -j 2 "bash paml_pipeline/Run_PAML_scripts_SKH_v3.sh output/reproduction_candidate_testrun/{}.aln.gs.fw.l0.fa output/ref/{}.ref output/target/{}.target output/root/{}.root"


bash paml_pipeline/Run_PAML_scripts_SKH.sh output/aln_with_paspalumCDS_forward/Pavag03G194200.aln.gs.fw.l0.fa output/ref/Pavag03G194200.ref output/target/Pavag03G194200.target output/root/Pavag03G194200.root

bash paml_pipeline/Run_PAML_scripts_SKH.sh output/aln_with_paspalumCDS_forward/Pavag02G345800.aln.gs.fw.l0.fa output/ref/Pavag02G345800.ref output/target/Pavag02G345800.target output/root/Pavag02G345800.root


##############################################################################################

find msa/ -type f| awk '{print substr($0,5,9)}'|parallel -j 11 "grep 'Pavag' msa/{}_mafft.aln|sed 's/>//g'|shuf -n 1 > ref/{}"

find msa/ -type f| awk '{print substr($0,5,9)}'|parallel -j 11 "grep 'Zm00001eb' msa/{}_mafft.aln|sed 's/>//g' > target/{}"

find msa/ -type f| awk '{print substr($0,5,9)}'|parallel -j 11 "bash Run"

bash ../paml_pipeline/Run_PAML_scripts_SKH.sh msa/OG0016935_mafft.aln ref/OG0016935 target/OG0016935 ref/OG0016935 
bash ../paml_pipeline/Run_PAML_scripts_SKH.sh msa/OG0004100_mafft.aln ref/OG0004100 target/OG0004100 ref/OG0004100 

bash ../paml_pipeline/Run_PAML_scripts_SKH.sh msa/OG0009148_mafft.aln ref/OG0009148 target/OG0009148 ref/OG0009148 

bash ../paml_pipeline/Run_PAML_scripts_SKH.sh msa/OG0005207_mafft.aln ref/OG0005207 target/OG0005207 ref/OG0005207 &

bash ../paml_pipeline/Run_PAML_scripts_SKH.sh msa/ ref/OG0005207 target/OG0005207 ref/OG0005207 &


grep 'zmB73\|zmhuet\|zluxur\|znicar\|zTIL\|irugos\|sbicol\|telega' msa/Pavag03G194200.aln.fa|sed 's/>//g' > target/Pavag03G194200

# split -l1 -d target/Pavag03G194200 target_split/Pavag03G194200_

bash ../paml_pipeline/Run_PAML_scripts_SKH.sh msa/Pavag03G194200.aln.fa ref/Pavag03G194200 target/Pavag03G194200 root/Pavag03G194200


parallel -j 20 'bash ../paml_pipeline/Run_PAML_scripts_SKH.sh msa/Pavag03G194200.aln.fa ref/Pavag03G194200 target_split/Pavag03G194200_{} root/Pavag03G194200' ::: {00..19}

grep 'Zm00001eb\|Zh\|Zn\|Zv\|Zx\|Ir\|Sobic\|Te' msa/OG0005207_mafft.aln|sed 's/>//g' > target/OG0005207

parallel -j 18 'bash ../paml_pipeline/Run_PAML_scripts_SKH.sh msa/OG0005207_mafft.aln ref/OG0005207 target_split/OG0005207_{} ref/OG0005207' ::: {00..17}


