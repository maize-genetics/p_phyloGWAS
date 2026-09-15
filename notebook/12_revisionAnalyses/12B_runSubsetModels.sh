#!/bin/bash
# Reruns 12_runPermulation_perOGModel.R (permulation=TRUE, GLM=FALSE by default -- the real
# pipeline, see 08D_perOGmodel.sh's fix commit) on one or more taxon subsets built
# by 12A_buildTraitSubsets.R, for envPC_1/2/3 each.
#
# Usage:
#   MC_CORES=8 notebook/12_revisionAnalyses/12B_runSubsetModels.sh <subsetName> [<subsetName> ...]
#   MC_CORES=8 notebook/12_revisionAnalyses/12B_runSubsetModels.sh --all   # every subset in manifest.csv
#   PERMULATION=FALSE MC_CORES=8 notebook/12_revisionAnalyses/12B_runSubsetModels.sh <subsetName>
#     -- skip the 1000-iteration permulation/emp_p step, base-model fit only. Use this when
#     only the rank-correlation comparison is needed (not the 330-candidate emp_p concordance
#     check) -- permulation convergence can be too low/slow to be worth running on a taxon-
#     reduced subset (observed as low as ~3% per iteration on the Andropogoneae subset,
#     implying ~35+ hours per envPC and still mostly NA emp_p afterward; see revision-analyses.md).
#
# subsetName must match a "subset" value in output/revisionAnalyses/traitTables/manifest.csv
# (e.g. andropogoneae, oneTipPerSpecies, drop_Andropogoninae, ...). Each subset x envPC
# combination runs sequentially (not in parallel with each other) -- only the within-run
# mclapply parallelism (MC_CORES, default 8) is used, to stay within a conservative core
# budget on a shared machine. Skips (does not rerun) any subset/envPC pair whose output
# already has what this invocation would produce (ASREML_res_empPadded.txt when permulation
# is on, plain ASREML_res.txt when it's off), so this script is safe to re-invoke to pick up
# where an earlier run left off.

set -euo pipefail
export PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"
export MC_CORES="${MC_CORES:-8}"
PERMULATION="${PERMULATION:-TRUE}"
RSCRIPT="${RSCRIPT:-/home/sh2246/miniconda3/envs/phyloGWAS-r/bin/Rscript}"

TRAIT_DIR="${PHYLOGWAS_ROOT}/output/revisionAnalyses/traitTables"
OUT_ROOT="${PHYLOGWAS_ROOT}/output/revisionAnalyses/perSubsetModels"
MASTER_TAB="${PHYLOGWAS_ROOT}/output/masterDataTable_PAVFill_20251001.txt"
PHYLOK="${PHYLOGWAS_ROOT}/output/phyloK_728Poaceae_astral_20250407.txt"
TREE="${PHYLOGWAS_ROOT}/output/PoaceaeTree_angiosperm353_astral_filtered_20250407.nwk"

manifest="${TRAIT_DIR}/manifest.csv"
[ -f "$manifest" ] || { echo "Missing $manifest -- run 12A_buildTraitSubsets.R first" >&2; exit 1; }

if [ "${1:-}" = "--all" ]; then
  subsets=$(tail -n +2 "$manifest" | cut -d',' -f1 | tr -d '"')
else
  [ $# -ge 1 ] || { echo "Usage: $0 <subsetName> [...] | --all" >&2; exit 1; }
  subsets="$*"
fi

traitFileFor() {
  case "$1" in
    andropogoneae) echo "envData_andropogoneae.txt" ;;
    oneTipPerSpecies) echo "envData_oneTipPerSpecies.txt" ;;
    drop_*) echo "envData_dropSubtribe_${1#drop_}.txt" ;;
    *) echo ""; return 1 ;;
  esac
}

mkdir -p "${PHYLOGWAS_ROOT}/log"

for subset in $subsets; do
  traitFile=$(traitFileFor "$subset") || { echo "Unknown subset: $subset" >&2; exit 1; }
  traitPath="${TRAIT_DIR}/${traitFile}"
  [ -f "$traitPath" ] || { echo "Missing trait table: $traitPath" >&2; exit 1; }
  outDir="${OUT_ROOT}/${subset}"
  mkdir -p "$outDir"

  for i in 1 2 3; do
    pcOut="${outDir}/envPC_${i}"
    if [ "$PERMULATION" = "TRUE" ] && [ -f "${pcOut}/ASREML_res_empPadded.txt" ]; then
      echo "[skip] ${subset} envPC_${i} -- already has ASREML_res_empPadded.txt"
      continue
    fi
    if [ "$PERMULATION" != "TRUE" ] && [ -f "${pcOut}/ASREML_res.txt" ]; then
      echo "[skip] ${subset} envPC_${i} -- already has ASREML_res.txt (permulation off)"
      continue
    fi
    logFile="${PHYLOGWAS_ROOT}/log/revisionAnalyses_${subset}_envPC${i}.log"
    echo "[run]  ${subset} envPC_${i} -> ${pcOut} (log: ${logFile}, permulation=${PERMULATION})"
    time "$RSCRIPT" "${PHYLOGWAS_ROOT}/src/12_runPermulation_perOGModel.R" \
      "$traitPath" "$MASTER_TAB" "$PHYLOK" "$TREE" \
      "envPC_${i}" "$outDir" "$PERMULATION" FALSE \
      > "$logFile" 2>&1
  done
done
