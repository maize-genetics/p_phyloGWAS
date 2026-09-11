#per OG envPC-activity model
#
# permulation=TRUE, GLM=FALSE is the real pipeline: this is what actually produced
# output/finalModels_20251002/ (ASREML_res.txt -> permulation/ -> ASREML_res_empPadded.txt,
# confirmed by file timestamps Oct 2 -> Oct 4 -> Oct 15 2025). A prior committed version of
# this script ran permulation=FALSE, GLM=TRUE -- leftover from GLM-variant testing, and
# produced only the separate, never-emp_p-padded output/finalModels_20251002/envPC_*/ASREML_GLMres.txt
# (Dec 3 2025), unrelated to the candidate-defining results. Fixed here; not rerun.

export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"

for i in {1..3}; do time Rscript src/12_runPermulation_perOGModel.R "${PHYLOGWAS_ROOT}/output/envData_707Poaceae_20250804.txt" "${PHYLOGWAS_ROOT}/output/masterDataTable_PAVFill_20251001.txt" "${PHYLOGWAS_ROOT}/output/phyloK_728Poaceae_astral_20250407.txt" "${PHYLOGWAS_ROOT}/output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk" envPC_${i} "${PHYLOGWAS_ROOT}/output/finalModels_20251002" TRUE FALSE > log/perOGModeling_envPC${i}_permulation.log 2>&1; done
