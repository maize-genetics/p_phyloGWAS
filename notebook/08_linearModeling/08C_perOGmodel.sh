#per OG envPC-activity model

export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"

for i in {1..3}; do time Rscript src/12_runPermulation_perOGModel.R "${PHYLOGWAS_ROOT}/output/envData_707Poaceae_20250804.txt" "${PHYLOGWAS_ROOT}/output/masterDataTable_PAVFill_20251001.txt" "${PHYLOGWAS_ROOT}/output/phyloK_728Poaceae_astral_20250407.txt" "${PHYLOGWAS_ROOT}/output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk" envPC_${i} "${PHYLOGWAS_ROOT}/output/finalModels_20251002" FALSE TRUE > log/perOGModeling_GLM_envPC${i}_NOpermulation.log 2>&1; done
