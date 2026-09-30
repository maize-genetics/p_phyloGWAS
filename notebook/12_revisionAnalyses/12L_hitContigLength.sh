#!/usr/bin/env bash
# Revision analysis 12L: length of the contigs that carry miniprot gene hits.
#
# 12K's contig N50 is low for short-read assemblies (median ~1.2 kb), but N50 is
# dominated by intergenic/repeat fragments. This asks the gene-space question
# directly: for every miniprot mRNA hit, how long is its host contig (gap-split
# at >= 10 N, same rule as 12K), how close is the hit to a contig end, and how
# much of the query does it cover. Per-hit tables go to output/hitContigLength/
# perHit/ (gitignored); 12L_hitContigLength_summary.R turns them into the
# tracked per-assembly summary and figure.
#
# ~1 min and < 1 GB per assembly for typical sizes.
set -euo pipefail

PHYLOGWAS_ROOT="${PHYLOGWAS_ROOT:-/workdir/sh2246/p_phyloGWAS}"
JOBS="${JOBS:-8}"
PY="${PY:-/home/sh2246/miniconda3/envs/daily/bin/python}"
SEQKIT="${SEQKIT:-/home/sh2246/miniconda3/envs/daily/bin/seqkit}"
SCRIPT="${PHYLOGWAS_ROOT}/notebook/12_revisionAnalyses/12L_hitContigLength.py"
OUTDIR="${PHYLOGWAS_ROOT}/output/hitContigLength"
MP="${PHYLOGWAS_ROOT}/output/orthofinderMiniProt"
MP_ADD="${PHYLOGWAS_ROOT}/output/orthofinderMiniProt_additionalOGs"
QLEN="${OUTDIR}/queryLengths.tsv"

mkdir -p "${OUTDIR}/perHit" "${PHYLOGWAS_ROOT}/log"

# Query (ancestral OG protein) lengths, for query coverage.
"$SEQKIT" fx2tab -n -i -l "${PHYLOGWAS_ROOT}/output/poaceaeHelixerOG_ancSeq_gapRemoved_v2_20240909.fa" \
  | cut -f1,2 > "$QLEN"

# Largest first so the long tail doesn't end up running alone.
# {= ... =} strips the directory and the .fa.gz suffix to get the assemblyID.
ls -S "${PHYLOGWAS_ROOT}"/data/assemblies/*.fa.gz \
  | parallel -j "$JOBS" --joblog "${PHYLOGWAS_ROOT}/log/12L_hitContigLength.joblog" \
      "$PY" "$SCRIPT" {} "$QLEN" "${OUTDIR}/perHit" \
      "${MP}/{= s:.*/::; s:\.fa\.gz$:: =}.gff" "${MP_ADD}/{= s:.*/::; s:\.fa\.gz$:: =}.gff"

echo "wrote $(ls "${OUTDIR}"/perHit/*.tsv.gz | wc -l) per-assembly hit tables to ${OUTDIR}/perHit"
