#!/usr/bin/env bash
set -euo pipefail

# ATAC-seq samples from GSE159654
#   SRR12849008 - ATAC_hESC_r1
#   SRR12849009 - ATAC_hESC_r2
#   SRR12849010 - ATAC_hPGC_r1
#   SRR12849011 - ATAC_hPGC_r2
#   SRR12849012 - ATAC_hPGCLC_d2_r1
#   SRR12849013 - ATAC_hPGCLC_d2_r2
#   SRR12849014 - ATAC_hPGCLC_d4_r1
#   SRR12849015 - ATAC_hPGCLC_d4_r2

declare -A SAMPLES=(
    ["SRR12849008"]="ATAC_hESC_r1"
    ["SRR12849009"]="ATAC_hESC_r2"
    ["SRR12849010"]="ATAC_hPGC_r1"
    ["SRR12849011"]="ATAC_hPGC_r2"
    ["SRR12849012"]="ATAC_hPGCLC_d2_r1"
    ["SRR12849013"]="ATAC_hPGCLC_d2_r2"
    ["SRR12849014"]="ATAC_hPGCLC_d4_r1"
    ["SRR12849015"]="ATAC_hPGCLC_d4_r2"
)

OUTDIR="fastqs"
mkdir -p "$OUTDIR"

for srr in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$srr]}"

    prefetch "$srr" --output-directory .

    fasterq-dump "$srr" --split-files --outdir "$OUTDIR" --temp .

    mv "$OUTDIR/${srr}_1.fastq" "$OUTDIR/${sample}_1.fastq"
    mv "$OUTDIR/${srr}_2.fastq" "$OUTDIR/${sample}_2.fastq"
    gzip "$OUTDIR/${sample}_1.fastq" &
    gzip "$OUTDIR/${sample}_2.fastq" &
    wait

    rm -rf "$srr"
done
