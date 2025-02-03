mkdir output/distMatrix/
Rscript src/01_MSA2Dist.R --input data/MSAs --output output/distMatrix/ --thread 10

Rscript src/02_getDistPAV.R --input output/distMatrix --reference B73_REF_FULL --output output/distPAV/

mkdir output/MSAs_gs
find data/MSAs -type f|cut -d / -f 3|sed 's/.fa//g' |parallel -j 30 "Rscript /workdir/sh2246/p_panAndOGASR/paml_pipeline/Ungap_MSA_SKH_v2.R data/MSAs/{}.fa B73_REF_FULL output/MSAs_gs/{}.gs.fa"

mkdir output/geneTree/
cd output/geneTree/
find /workdir/hackathon_202309/output/MSAs_gs -name *.gs.fa -type f| cut -d / -f 6|sed 's/_msa_padded10N_rc.gs.fa//g' | parallel -j 30 '/workdir/sh2246/p_panAndOGASR/paml_pipeline/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s /workdir/hackathon_202309/output/MSAs_gs/{}_msa_padded10N_rc.gs.fa -# 1 -n {}.tree'

mkdir output/simPAV/ePC1
parallel -j 8 'mkdir output/simPAV/ePC1/freq_{}' ::: 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8
parallel -j 8 'Rscript src/03_simulatePAV.R --treeDir data/speciesTrees_finaltaxa_median4dDz_tabasco4000.newick --pMatrix data/ePCs.09.14.23.csv --traitID ePC_1 --freq {} --replicate 100 --output output/simPAV/ePC1/freq_{}/simPAV_freq_{}' ::: 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8

mkdir output/simPAV/lifeHistory
parallel -j 8 'mkdir output/simPAV/lifeHistory/freq_{}' ::: 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8
parallel -j 8 'Rscript src/03B_simulatePAV_cat.R --treeDir data/speciesTrees_finaltaxa_median4dDz_tabasco4000.newick --pMatrix data/Perennitality_Numeric_Phenotype_325Taxa.csv --traitID lifeHistory_01 --freq {} --replicate 100 --output output/simPAV/lifeHistory/freq_{}/simPAV_freq_{}' ::: 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8

