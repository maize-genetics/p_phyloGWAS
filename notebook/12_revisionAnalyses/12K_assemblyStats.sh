#!/bin/bash
# Standard assembly statistics for every assembly in data/assemblies/ (727 files).
# One row per assembly, keyed by assemblyID (filename minus .fa.gz), so it joins
# to any version of the Poaceae metadata. See 12K_assemblyStats.py for the stat
# definitions (as-deposited, contig-level split at >=10 N, and >=1/10 kb bins).
#
# Usage: JOBS=8 notebook/12_revisionAnalyses/12K_assemblyStats.sh
# Each job uses ~2 seqkit threads plus one python process; memory peaks at a few
# GB on the largest chromosome-level assemblies (A17, 14 Gb: ~2.5 GB).

set -euo pipefail
PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"
JOBS="${JOBS:-8}"
PY="${PY:-/home/sh2246/miniconda3/envs/daily/bin/python}"
SCRIPT="${PHYLOGWAS_ROOT}/notebook/12_revisionAnalyses/12K_assemblyStats.py"
OUTDIR="${PHYLOGWAS_ROOT}/output/assemblyStats"
OUT="${OUTDIR}/assemblyStats_727.tsv"
mkdir -p "$OUTDIR" "${PHYLOGWAS_ROOT}/log"

# largest first, so the long tail doesn't end on one big file running alone
ls -S "${PHYLOGWAS_ROOT}"/data/assemblies/*.fa.gz \
  | parallel -j "$JOBS" --joblog "${PHYLOGWAS_ROOT}/log/12K_assemblyStats.joblog" \
      "$PY" "$SCRIPT" {} \
  > "${OUTDIR}/rows.tmp"

{ "$PY" "$SCRIPT" --header; sort -k1,1 "${OUTDIR}/rows.tmp"; } > "$OUT"
rm "${OUTDIR}/rows.tmp"
echo "wrote $OUT: $(($(wc -l < "$OUT") - 1)) assemblies"
