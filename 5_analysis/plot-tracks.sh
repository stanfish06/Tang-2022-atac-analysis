#!/usr/bin/env bash
# ================== plot tracks  ==================
# plot ATAC and motif tracks
# deps:
#   - pyGenomeTracks
# usage:
#   - chmod +x plot_tracks.sh && ./plot_tracks.sh
# ==================================================

set -oue pipefail

PROJECT_DIR=".."

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

TARGETS_FILE="targets.bed"
GENES_FILE="$PROJECT_DIR/0_ref/gencode.v19.annotation.gtf.gz"
TOP_MOTIFS=("ERG" "ESR1" "ETS1" "RUNX2" "SNAI2" "TAL1")
COLORS=("#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30" "#4DBEEE" "#A2142F")

for sample in "${SAMPLES[@]}"; do
	echo "=== Processing sample: $sample ==="

	BW_FILE="bigwigs/${sample}.bw"
	SAMPLE_DIR="$sample"
	OUTPUT_DIR="${sample}/tracks"
	mkdir -p "$OUTPUT_DIR"

	while read -r chrom start end name; do

		echo "Processing $name..."

		gene_dir="${SAMPLE_DIR}/fimo_out/erk_viz_${name}"
		fimo_gff="${gene_dir}/fimo.gff"
		bed_dir="${gene_dir}/beds"

		if [[ -f "$fimo_gff" ]]; then
			mkdir -p "$bed_dir"
			awk -v out_dir="$bed_dir" '
        BEGIN {OFS="\t"}
        !/^#/ {
            split($9, a, ";");
            motif_name = "";
            for (i in a) {
                if (a[i] ~ /^Name=/) {
                    name_attr = substr(a[i], 6);
                    split(name_attr, n, "_");
                    motif_name = n[1];
                }
            }
            if (motif_name != "") {
                # BED: chrom start end name score strand
                print $1, $4-1, $5, motif_name, $6, $7 >> (out_dir "/" motif_name ".bed");
            }
        }' "$fimo_gff"
		else
			echo "  > Warning: FIMO GFF not found ($fimo_gff). Skipping."
			continue
		fi

		ini_file="${OUTPUT_DIR}/config_${name}.ini"

		cat <<EOF >"$ini_file"
[x-axis]
fontsize = 16

[ATAC-seq Signal]
file = $BW_FILE
title = ATAC-seq
height = 4
min_value = 0
max_value = 0.5
file_type = bigwig

EOF

		c_idx=0

		for motif in "${TOP_MOTIFS[@]}"; do
			bed_file="$bed_dir/${motif}.bed"
			color="${COLORS[$c_idx]}"

			if [[ -f "$bed_file" ]]; then
				cat <<EOF >>"$ini_file"
[$motif]
file = $bed_file
title = $motif
height = 0.5
color = $color
labels = false
display = collapsed
border_color = none
file_type = bed

EOF
			else
				:
			fi
			c_idx=$(((c_idx + 1) % 7))
		done

		cat <<EOF >>"$ini_file"
[spacer]
height = 1

[Genes]
file = $GENES_FILE
title = Genes
height = 3
color = darkblue
height_utr = 1
style = UCSC
file_type = gtf
gene_rows = 4
fontsize = 14
EOF

		output="${OUTPUT_DIR}/${name}.png"
		region="${chrom}:${start}-${end}"

		echo "  > Plotting to $output..."
		pyGenomeTracks \
			--tracks "$ini_file" \
			--region "$region" \
			--outFileName "$output" \
			--dpi 600 \
			--width 40 \
			--fontSize 10

	done <"$TARGETS_FILE"

	echo "=== Done: $sample ==="
done

echo "All done"
