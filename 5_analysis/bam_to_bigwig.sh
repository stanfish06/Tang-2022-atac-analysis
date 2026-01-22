#!/usr/bin/env bash
# ================== Create bigwigs for atac datasets ==================
# Generate bigwigs files for genome track visualization
# deps:
# usage:
#   - chmod +x bam_to_bigwigs.sh && ./bam_to_bigwigs.sh
# ======================================================================

set -oue pipefail

PROJECT_DIR=".."
BAM_DIR="$PROJECT_DIR/3_align/bam"
BW_DIR="$PROJECT_DIR/5_analysis/bigwigs"
mkdir -p "$BW_DIR"

find "$BAM_DIR" -name "*.filtered.bam" | while read -r bam_file; do
	filename=$(basename "$bam_file")
	sample="${filename%.filtered.bam}"
	out_bw="$BW_DIR/${sample}.bw"

	if [[ ! -f "${bam_file}.bai" ]]; then
		echo "Indexing $bam_file..."
		samtools index "$bam_file"
	fi

	echo "Processing $sample..."
	echo "Input: $bam_file"
	echo "Output: $out_bw"

	bamCoverage \
		--bam "$bam_file" \
		--outFileName "$out_bw" \
		--outFileFormat bigwig \
		--normalizeUsing CPM \
		--binSize 10 \
		--smoothLength 30 \
		--numberOfProcessors 8

	echo "Done: $out_bw"
done
echo "All done"
