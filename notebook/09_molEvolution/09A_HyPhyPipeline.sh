#!/bin/bash

export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"
cd "${PHYLOGWAS_ROOT}/"


## step 1: MSA cleaning
mkdir output/CDSMSAPerOG_HyPhy_20250203
find output/CDSMSAPerOG_gs/ -name "*.l0.fa" -type f| cut -d / -f 3|sed 's/.gs.l0.fa//g' | parallel -j 40 " /programs/hyphy-2.5.49/bin/hyphy cln Universal output/CDSMSAPerOG_gs/{}.gs.l0.fa No/No output/CDSMSAPerOG_HyPhy_20250203/{}.fa" > output/MSAcleaning.log 2>&1

## step 2: RAxML gene tree generation
mkdir output/geneTree_allOGs_20250203/
find "${PHYLOGWAS_ROOT}/output/CDSMSAPerOG_HyPhy_20250203/" -name '*.fa' -type f| cut -d / -f 7|sed 's/.fa//g' | parallel -j 35 'src/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s ${PHYLOGWAS_ROOT}/output/CDSMSAPerOG_HyPhy_20250203/{}.fa -# 1 -w ${PHYLOGWAS_ROOT}/output/geneTree_allOGs_20250203/ -n {}.tree'

# clean up MSA and rebuild tree for rhizome analyses (keep only perennial)
mkdir output/CDSMSAPerOG_HyPhy_rhizome_20260303
cat output/candidateOG_rhizome.txt | parallel -j 35 'PAT=$(paste -sd"|" data/perennial_assemblies_20260213.txt)
/programs/seqkit-0.15.0/seqkit grep -n -r -p ":(${PAT}):" output/CDSMSAPerOG_gs/{}.gs.l0.fa -o output/CDSMSAPerOG_HyPhy_rhizome_20260303/{}.gs.l0.fa'

find output/CDSMSAPerOG_HyPhy_rhizome_20260303/ -name "*.l0.fa" -type f| cut -d / -f 3|sed 's/.gs.l0.fa//g' | parallel -j 40 " /programs/hyphy-2.5.49/bin/hyphy cln Universal output/CDSMSAPerOG_HyPhy_rhizome_20260303/{}.gs.l0.fa No/No output/CDSMSAPerOG_HyPhy_rhizome_20260303/{}.fa" > output/MSAcleaning_perennial.log 2>&1
rm output/CDSMSAPerOG_HyPhy_rhizome_20260303/*.gs.l0.fa
mkdir output/geneTree_perennialOnly_rhizomeOGs_20260303/
find "${PHYLOGWAS_ROOT}/output/CDSMSAPerOG_HyPhy_rhizome_20260303/" -name '*.fa' -type f| cut -d / -f 7|sed 's/.fa//g' | parallel -j 35 'src/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s ${PHYLOGWAS_ROOT}/output/CDSMSAPerOG_HyPhy_rhizome_20260303/{}.fa -# 1 -w ${PHYLOGWAS_ROOT}/output/geneTree_perennialOnly_rhizomeOGs_20260303/ -n {}.tree'

# clean up MSA and rebuild tree for annual analyses (keep only nonrhizomatous)
mkdir output/CDSMSAPerOG_HyPhy_LH_20260415
cat output/candidateOG_lifeHistory2.txt | parallel -j 35 'PAT=$(paste -sd"|" data/nonrhizomatous_assemblies_20260415.txt)
/programs/seqkit-0.15.0/seqkit grep -n -r -p ":(${PAT}):" output/CDSMSAPerOG_gs/{}.gs.l0.fa -o output/CDSMSAPerOG_HyPhy_LH_20260415/{}.gs.l0.fa'

find output/CDSMSAPerOG_HyPhy_LH_20260415/ -name "*.l0.fa" -type f| cut -d / -f 3|sed 's/.gs.l0.fa//g' | parallel -j 40 " /programs/hyphy-2.5.49/bin/hyphy cln Universal output/CDSMSAPerOG_HyPhy_LH_20260415/{}.gs.l0.fa No/No output/CDSMSAPerOG_HyPhy_LH_20260415/{}.fa" > output/MSAcleaning_nonrhizomatous.log 2>&1
rm output/CDSMSAPerOG_HyPhy_LH_20260415/*.gs.l0.fa
mkdir output/geneTree_nonrhizomatousOnly_LHOGs_20260415/
find "${PHYLOGWAS_ROOT}/output/CDSMSAPerOG_HyPhy_LH_20260415/" -name '*.fa' -type f| cut -d / -f 7|sed 's/.fa//g' | parallel -j 35 'src/standard-RAxML/raxmlHPC -m GTRGAMMA -p 12345 -s ${PHYLOGWAS_ROOT}/output/CDSMSAPerOG_HyPhy_LH_20260415/{}.fa -# 1 -w ${PHYLOGWAS_ROOT}/output/geneTree_nonrhizomatousOnly_LHOGs_20260415/ -n {}.tree'


## step 3: get target tip names
# life history
mkdir output/targetHyPhy_lifeHistory_20260213
cat output/candidateOG_lifeHistory.txt |parallel -j 40 "grep -f data/annual_assemblies_20260213.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_lifeHistory_20260213/{}.target"

mkdir output/targetHyPhy_lifeHistory_20260416
cat output/candidateOG_lifeHistory2.txt |parallel -j 40 "grep -f data/annual_assemblies_20260213.txt output/CDSMSAPerOG_HyPhy_LH_20260415/{}.fa|sed 's/>//g' > output/targetHyPhy_lifeHistory_20260416/{}.target"

mkdir output/targetHyPhy_lifeHistory_20260515
cat output/candidateOG_LH3ClassModel.txt |parallel -j 40 "grep -f data/annual_assemblies_20260213.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_lifeHistory_20260515/{}.target"

# rhizome
mkdir output/targetHyPhy_rhizome_20260303
cat output/candidateOG_rhizome.txt |parallel -j 40 "grep -f data/perennialRhizome_assemblies_20260303.txt output/CDSMSAPerOG_HyPhy_rhizome_20260303/{}.fa|sed 's/>//g' > output/targetHyPhy_rhizome_20260303/{}.target"

mkdir output/targetHyPhy_rhizome_20260515
cat output/candidateOG_LH3ClassModel.txt |parallel -j 40 "grep -f data/perennialRhizome_assemblies_20260303.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_rhizome_20260515/{}.target"

#warm/cold adapted
mkdir output/targetHyPhy_cold_20250415
find output/CDSMSAPerOG_HyPhy_20250203/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' |parallel -j 40 "grep -f output/coldAdaptedAssemblies.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_cold_20250415/{}.target"
mkdir output/targetHyPhy_warm_20250415
find output/CDSMSAPerOG_HyPhy_20250203/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' |parallel -j 40 "grep -f output/warmAdaptedAssemblies.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_warm_20250415/{}.target"
#wet/dry adapted
mkdir output/targetHyPhy_drought_20251204
cat output/candidateOG_envPC2.txt|parallel -j 40 "grep -f output/droughtAdaptedAssemblies.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_drought_20251204/{}.target"
mkdir output/targetHyPhy_wet_20251204
cat output/candidateOG_envPC2.txt|parallel -j 40 "grep -f output/wetAdaptedAssemblies.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_wet_20251204/{}.target"

#sand/clay adapted
mkdir output/targetHyPhy_sand_20251204
cat output/candidateOG_envPC3.txt|parallel -j 40 "grep -f output/sandAdaptedAssemblies.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_sand_20251204/{}.target"
mkdir output/targetHyPhy_clay_20251204
cat output/candidateOG_envPC3.txt|parallel -j 40 "grep -f output/clayAdaptedAssemblies.txt output/CDSMSAPerOG_HyPhy_20250203/{}.fa|sed 's/>//g' > output/targetHyPhy_clay_20251204/{}.target"

## step 4: tree labeling
# life history
mkdir output/labeledGeneTree_lifeHistory
cat output/candidateOG_lifeHistory.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_lifeHistory_20260213/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_lifeHistory/RAxML_Labeled_bestTree.{}_Relax.tree"

mkdir output/labeledGeneTree_lifeHistory_20260416
cat output/candidateOG_lifeHistory2.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_nonrhizomatousOnly_LHOGs_20260415/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_lifeHistory_20260416/{}.target output/CDSMSAPerOG_HyPhy_LH_20260415/{}.fa output/labeledGeneTree_lifeHistory_20260416/RAxML_Labeled_bestTree.{}_Relax.tree"

mkdir output/labeledGeneTree_lifeHistory_20260515
cat output/candidateOG_LH3ClassModel.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_lifeHistory_20260515/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_lifeHistory_20260515/RAxML_Labeled_bestTree.{}_Relax.tree"

# rhizome
mkdir output/labeledGeneTree_rhizome
cat output/candidateOG_rhizome.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_perennialOnly_rhizomeOGs_20260303/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_rhizome_20260303/{}.target output/CDSMSAPerOG_HyPhy_rhizome_20260303/{}.fa output/labeledGeneTree_rhizome/RAxML_Labeled_bestTree.{}_Relax.tree"

mkdir output/labeledGeneTree_rhizome_20260515
cat output/candidateOG_LH3ClassModel.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_rhizome_20260515/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_rhizome_20260515/RAxML_Labeled_bestTree.{}_Relax.tree"


#warm/cold adapted
mkdir output/labeledGeneTree_cold
find output/CDSMSAPerOG_HyPhy_20250203/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_cold_20250415/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_cold/RAxML_Labeled_bestTree.{}_Relax.tree"
mkdir output/labeledGeneTree_warm
find output/CDSMSAPerOG_HyPhy_20250203/ -name "*.fa" -type f| cut -d / -f 3|sed 's/.fa//g' | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_warm_20250415/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_warm/RAxML_Labeled_bestTree.{}_Relax.tree"

#wet/dry adapted
mkdir output/labeledGeneTree_drought
cat output/candidateOG_envPC2.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_drought_20251204/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_drought/RAxML_Labeled_bestTree.{}_Relax.tree"
mkdir output/labeledGeneTree_wet
cat output/candidateOG_envPC2.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_wet_20251204/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_wet/RAxML_Labeled_bestTree.{}_Relax.tree"

#sand/clay adapted
mkdir output/labeledGeneTree_sand
cat output/candidateOG_envPC3.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_sand_20251204/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_sand/RAxML_Labeled_bestTree.{}_Relax.tree"
mkdir output/labeledGeneTree_clay
cat output/candidateOG_envPC3.txt | parallel -j 30 "Rscript src/LabelNodes_SKH_v7_HyPhyRelax.R output/geneTree_allOGs_20250203/RAxML_bestTree.{}.tree ASM1935983v1 output/targetHyPhy_clay_20251204/{}.target output/CDSMSAPerOG_HyPhy_20250203/{}.fa output/labeledGeneTree_clay/RAxML_Labeled_bestTree.{}_Relax.tree"


## step 5: HyPhy RELAX test
mkdir output/HyPhyResult
mkdir output/HyPhyResult/lifeHistory
cat output/candidateOG_lifeHistory.txt | parallel -j 30 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_20250203/{}.fa --tree output/labeledGeneTree_lifeHistory/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/lifeHistory/{}.RELAX.json"

mkdir output/HyPhyResult/lifeHistory_20260416
cat output/candidateOG_lifeHistory2.txt | parallel -j 30 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_LH_20260415/{}.fa --tree output/labeledGeneTree_lifeHistory_20260416/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/lifeHistory_20260416/{}.RELAX.json"

mkdir output/HyPhyResult/rhizome
cat output/candidateOG_rhizome.txt | parallel -j 30 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_rhizome_20260303/{}.fa --tree output/labeledGeneTree_rhizome/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/rhizome/{}.RELAX.json"

mkdir output/HyPhyResult/lifeHistory_20260515
cat output/candidateOG_LH3ClassModel.txt | parallel -j 30 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_20250203/{}.fa --tree output/labeledGeneTree_lifeHistory_20260515/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/lifeHistory_20260515/{}.RELAX.json"

mkdir output/HyPhyResult/rhizome_20260515
cat output/candidateOG_LH3ClassModel.txt | parallel -j 30 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_20250203/{}.fa --tree output/labeledGeneTree_rhizome_20260515/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/rhizome_20260515/{}.RELAX.json"

mkdir output/HyPhyResult/envPC2_drought
mkdir output/HyPhyResult/envPC2_wet
cat output/candidateOG_envPC2.txt | parallel -j 30 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_20250203/{}.fa --tree output/labeledGeneTree_drought/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/envPC2_drought/{}.RELAX.json"
cat output/candidateOG_envPC2.txt | parallel -j 30 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_20250203/{}.fa --tree output/labeledGeneTree_wet/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/envPC2_wet/{}.RELAX.json"

mkdir output/HyPhyResult/envPC3_sand
mkdir output/HyPhyResult/envPC3_clay
cat output/candidateOG_envPC3.txt | parallel -j 35 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_20250203/{}.fa --tree output/labeledGeneTree_sand/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/envPC3_sand/{}.RELAX.json"
cat output/candidateOG_envPC3.txt | parallel -j 35 "/programs/hyphy-2.5.49/bin/hyphy relax --alignment output/CDSMSAPerOG_HyPhy_20250203/{}.fa --tree output/labeledGeneTree_clay/RAxML_Labeled_bestTree.{}_Relax.tree --test Foreground --reference Background --models Minimal --output output/HyPhyResult/envPC3_clay/{}.RELAX.json"

## to extract test statistics from json files into final output



