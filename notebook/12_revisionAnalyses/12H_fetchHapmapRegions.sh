#!/usr/bin/env bash
# Download the Panzea HapMap3.2.1 282-panel VCFs (B73v5, SnpEff-annotated,
# 277 samples) from MaizeGDB's Box share, then extract the target-gene regions
# into one small VCF.
#
# Remote tabix over HTTPS against Box was tried first but is unreliable (Box's
# tokenised redirect URLs break htslib range requests), so the per-chromosome
# files are downloaded in full (5.75 GB, ~2 min at ~70 MB/s) and kept under
# data/hapmap282/full/ for reuse. Files already present with the expected size
# are not re-downloaded.
#
# Usage: bash src/maize_02a_fetch_hapmap_regions.sh <target_genes.bed> <out.vcf.gz> [download_dir]
#   target_genes.bed: chrom  start  end  strand  geneID  (1-based, from GFF3)
#
# Box share: https://iastate.app.box.com/s/9edkm5f51q88ltzrb0hxaz4p36fdsr7q/folder/364949363629
set -euo pipefail
BED=$1
OUT=$2
DL=${3:-data/hapmap282/full}
mkdir -p "$DL"

BASE="https://iastate.app.box.com/index.php?rm=box_download_shared_file&shared_name=9edkm5f51q88ltzrb0hxaz4p36fdsr7q&file_id="
# file_id and expected byte size, taken from the Box folder listing (2026-02 upload by MaizeGDB)
declare -A VCF=( [chr1]=f_2130578606035:879752128 [chr2]=f_2130575682763:686351906 [chr3]=f_2130565447010:664030955
                 [chr4]=f_2130570899115:685153194 [chr5]=f_2130579254281:568338711 [chr6]=f_2130565599717:460079204
                 [chr7]=f_2130570462632:486018751 [chr8]=f_2130569793072:469949139 [chr9]=f_2130571207139:439869790
                 [chr10]=f_2130575142039:413000001 )
declare -A TBI=( [chr1]=f_2130559398170 [chr2]=f_2130566946970 [chr3]=f_2130567490839 [chr4]=f_2130573639646 [chr5]=f_2130562468058
                 [chr6]=f_2130576892054 [chr7]=f_2130570215608 [chr8]=f_2130578449406 [chr9]=f_2130571957694 [chr10]=f_2130565858948 )

fetch() { # url  dest  [expected_size]
  local url=$1 dest=$2 exp=${3:-}
  if [ -s "$dest" ] && { [ -z "$exp" ] || [ "$(stat -c %s "$dest")" -eq "$exp" ]; }; then
    echo "  have $dest"; return
  fi
  curl -sL --retry 3 -o "$dest" "$url"
  if [ -n "$exp" ] && [ "$(stat -c %s "$dest")" -ne "$exp" ]; then
    echo "size mismatch for $dest (got $(stat -c %s "$dest"), expected $exp)" >&2; exit 1
  fi
  echo "  fetched $dest"
}

TMP=$(mktemp -d)
first=1
for chr in $(cut -f1 "$BED" | sort -V | uniq); do
  IFS=: read -r fid size <<< "${VCF[$chr]}"
  vcf="$DL/panzea_hapmap_282_${chr}.vcf.gz"
  echo "[$chr]"
  fetch "${BASE}${fid}" "$vcf" "$size"
  fetch "${BASE}${TBI[$chr]}" "$vcf.tbi"
  regions=$(awk -v c="$chr" '$1==c {printf "%s:%s-%s ", $1, $2, $3}' "$BED")
  if [ $first -eq 1 ]; then tabix -h "$vcf" $regions > "$TMP/all.vcf"; first=0
  else tabix "$vcf" $regions >> "$TMP/all.vcf"; fi
  echo "  regions: $regions -> $(grep -vc '^#' "$TMP/all.vcf") records so far"
done

bgzip -c "$TMP/all.vcf" > "$OUT"
tabix -p vcf "$OUT"
rm -rf "$TMP"
echo "wrote $OUT ($(zcat "$OUT" | grep -vc '^#') records)"
