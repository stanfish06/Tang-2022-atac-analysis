#!/usr/bin/env bash
# ============ Read alignemnt ============
# Align trimmed fastq data to genome
# deps:
#   - bowtie2
# usage:
#   - chmod +x align.sh && ./align.sh
# =======================================

#SBATCH --job-name=atac_align
#SBATCH --output=logs/align_%A_%a.out
#SBATCH --error=logs/align_%A_%a.err
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=zyyu@umich.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --time=8:00:00
#SBATCH --array=0-7

set -euo pipefail

PROJECT_DIR=".."
INDEX="$PROJECT_DIR/0_ref/bowtie2_index/GRCh37"
INPUT_DIR="$PROJECT_DIR/2_trim"
THREADS=12

SAMPLES=(
	"ATAC_hESC_r1"
	"ATAC_hESC_r2"
	"ATAC_hPGC_r1"
	"ATAC_hPGC_r2"
	"ATAC_hPGCLC_d2_r1"
	"ATAC_hPGCLC_d2_r2"
	"ATAC_hPGCLC_d4_r1"
	"ATAC_hPGCLC_d4_r2"
)

mkdir -p bam metrics logs

run_sample() {
	local sample="$1"
	local r1="$INPUT_DIR/${sample}_1_val_1.fq.gz"
	local r2="$INPUT_DIR/${sample}_2_val_2.fq.gz"

	if [[ -f "bam/${sample}.filtered.bam" ]]; then
		echo "[$sample] Already done, skipping"
		return 0
	fi

	echo "[$sample] Aligning..."
	bowtie2 --very-sensitive -X 2000 -k 10 -p "$THREADS" -x "$INDEX" \
		--rg-id "$sample" --rg "SM:$sample" \
		-1 "$r1" -2 "$r2" 2>"logs/${sample}.bowtie2.log" |
		samtools view -bS -@ 2 - >"bam/${sample}.bam"

	echo "[$sample] Sorting..."
	samtools sort -@ "$THREADS" -m 2G -o "bam/${sample}.sorted.bam" "bam/${sample}.bam"
	samtools index "bam/${sample}.sorted.bam"

	echo "[$sample] Marking duplicates..."
	picard MarkDuplicates \
		I="bam/${sample}.sorted.bam" \
		O="bam/${sample}.dedup.bam" \
		M="metrics/${sample}.dup_metrics.txt" \
		VALIDATION_STRINGENCY=LENIENT QUIET=true

	echo "[$sample] Filtering..."
	samtools view -h -q 30 -F 1024 "bam/${sample}.dedup.bam" |
		awk 'BEGIN{OFS="\t"} $1 ~ /^@/ || $3 != "chrM"' |
		samtools view -b -@ 2 -o "bam/${sample}.filtered.bam"
	samtools index "bam/${sample}.filtered.bam"

	echo "[$sample] Done"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
	run_sample "${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
else
	for sample in "${SAMPLES[@]}"; do
		run_sample "$sample"
	done
	echo "All done"
fi
