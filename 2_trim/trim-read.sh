#!/usr/bin/env bash
# ================== Trim reads =====================
# Trim reads (using standard parameters for atac)
# deps:
#   - trim galore
# usage:
#   - chmod +x trim-read.sh && ./trim-read.sh
# ===================================================

set -oue pipefail

PROJECT_DIR=".."
FASTQ_DIR="$PROJECT_DIR/1_fetch/"

find "$FASTQ_DIR" -name "*_1.fastq*" | sort | while read -r r1; do
	r2="${r1/_1.fastq/_2.fastq}"
	base=$(basename "$r1" | sed 's/_1\.fastq.*//')
	out1="${base}_1_val_1.fq.gz"
	out2="${base}_2_val_2.fq.gz"

	if [[ -f "$out1" && -f "$out2" ]]; then
		echo "Skipping: ${base} (already trimmed)"
		continue
	fi

	echo "Processing: $(basename "$r1")"
	echo "            $(basename "$r2")"

	trim_galore \
		--nextera \
		--paired \
		--fastqc \
		--length 20 \
		-j 4 \
		--output_dir . \
		"$r1" "$r2"

	echo "Done: $(basename "$r1")"
	echo "      $(basename "$r2")"
done
