#!/usr/bin/env bash
# ================== Find atac peaks =====================
# Find peaks in atac data using macs2
# deps:
#   - macs2
# usage:
#   - chmod +x find-atac-peaks.sh && ./find-atac-peaks.sh
# ========================================================

set -oue pipefail

PROJECT_DIR=".."
BAM_DIR="$PROJECT_DIR/3_align/bam"
pval="0.5" # be permissive here, can filter peaks later

shopt -s nullglob
bam_files=("$BAM_DIR"/*.filtered.bam)
shopt -u nullglob

for bam in "${bam_files[@]}"; do
  base="$(basename "$bam")"
  sample="${base%.filtered.bam}"
  sample_out="$sample"
  mkdir -p "$sample_out"

  macs2 callpeak \
    -t "$bam" \
    -f BAMPE \
    -g hs \
    --keep-dup all \
    --nomodel \
    -p "$pval" \
    -n "$sample" \
    --nolambda \
    --shift -100 \
    --extsize 200 \
    --outdir "$sample_out"
  done
