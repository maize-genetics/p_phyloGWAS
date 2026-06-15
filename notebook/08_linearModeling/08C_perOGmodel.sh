#per OG envPC-activity model

for i in {1..3}; do time Rscript src/12_runPermulation_perOGModel.R /workdir/sh2246/p_phyloGWAS/output/envData_707Poaceae_20250804.txt /workdir/sh2246/p_phyloGWAS/output/masterDataTable_PAVFill_20251001.txt /workdir/sh2246/p_phyloGWAS/output/phyloK_728Poaceae_astral_20250407.txt /workdir/sh2246/p_phyloGWAS/output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk envPC_${i} /workdir/sh2246/p_phyloGWAS/output/finalModels_20251002 FALSE TRUE > log/perOGModeling_GLM_envPC${i}_NOpermulation.log 2>&1; done
