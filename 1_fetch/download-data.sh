#!/usr/bin/env bash
# ================== Data download ==================
# download fastq files
# deps:
#   - fasterq-dump
#   - prefetch
# usage:
#   - chmod +x download-data.sh && ./download-data.sh
#   - add more datasets to <SAMPLES> if needed
# ===================================================

set -oue pipefail

# ATAC-seq samples from GSE159654
#   SRR12849008 - ATAC_hESC_r1
#   SRR12849009 - ATAC_hESC_r2
#   SRR12849010 - ATAC_hPGC_r1
#   SRR12849011 - ATAC_hPGC_r2
#   SRR12849012 - ATAC_hPGCLC_d2_r1
#   SRR12849013 - ATAC_hPGCLC_d2_r2
#   SRR12849014 - ATAC_hPGCLC_d4_r1
#   SRR12849015 - ATAC_hPGCLC_d4_r2

# uncomment other datasets if needed
declare -A SAMPLES=(
	# ["SRR12849008"]="ATAC_hESC_r1"
	# ["SRR12849009"]="ATAC_hESC_r2"
	# ["SRR12849010"]="ATAC_hPGC_r1"
	# ["SRR12849011"]="ATAC_hPGC_r2"
	["SRR12849012"]="ATAC_hPGCLC_d2_r1"
	["SRR12849013"]="ATAC_hPGCLC_d2_r2"
	# ["SRR12849014"]="ATAC_hPGCLC_d4_r1"
	# ["SRR12849015"]="ATAC_hPGCLC_d4_r2"
)

for srr in "${!SAMPLES[@]}"; do
	sample="${SAMPLES[$srr]}"
	if [[ -d "$srr" ]]; then
		echo "Skipping prefetch: ${sample} (already exists)"
	else
		echo "Prefetching: ${sample}"
		prefetch "$srr" --output-directory . &
	fi
done
wait

for srr in "${!SAMPLES[@]}"; do
	(
		echo "Generating fastqs: ${SAMPLES[$srr]}"
		sample="${SAMPLES[$srr]}"
		tmpdir="tmp_${srr}"
		mkdir -p "$tmpdir"
		fasterq-dump "$srr" --split-files --outdir . --temp "$tmpdir"

		mv "${srr}_1.fastq" "${sample}_1.fastq"
		mv "${srr}_2.fastq" "${sample}_2.fastq"
		gzip "${sample}_1.fastq" &
		gzip "${sample}_2.fastq" &
		wait
		rm -rf "$srr" "$tmpdir"
	) &
done
wait
echo "Data download finished."
