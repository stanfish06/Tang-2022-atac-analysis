#!/usr/bin/env bash
# ============ Fetch aligned data from Zenodo ============
# Download pre-aligned BAM files from Zenodo
# deps:
#   - curl or wget
# usage:
#   - chmod +x fetch-from-s3.sh && ./fetch-from-s3.sh
# ========================================================

set -euo pipefail

ZENODO_RECORD="XXXXXXX"
ZENODO_URL="https://zenodo.org/records/${ZENODO_RECORD}/files"

SAMPLES=(
	# "ATAC_hESC_r1"
	# "ATAC_hESC_r2"
	# "ATAC_hPGC_r1"
	# "ATAC_hPGC_r2"
	"ATAC_hPGCLC_d2_r1"
	"ATAC_hPGCLC_d2_r2"
	# "ATAC_hPGCLC_d4_r1"
	# "ATAC_hPGCLC_d4_r2"
)

mkdir -p bam

for sample in "${SAMPLES[@]}"; do
	if [[ -f "bam/${sample}.filtered.bam" ]]; then
		echo "[$sample] Already exists, skipping"
		continue
	fi

	echo "[$sample] Downloading..."
	curl -L -o "bam/${sample}.filtered.bam" "$ZENODO_URL/${sample}.filtered.bam"
	curl -L -o "bam/${sample}.filtered.bam.bai" "$ZENODO_URL/${sample}.filtered.bam.bai"
done

echo "All done"
