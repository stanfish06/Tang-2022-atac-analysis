#!/usr/bin/env bash
set -euo pipefail

INDEX="../ref/bowtie2_index/GRCh37"
INPUT_DIR="../2_trim"
THREADS=12

SAMPLES=(
  "ATAC_hPGCLC_d2_r1"
  "ATAC_hESC_r1"
)

mkdir -p bam metrics logs

for sample in "${SAMPLES[@]}"; do
  r1="$INPUT_DIR/${sample}_1_val_1.fq.gz"
  r2="$INPUT_DIR/${sample}_2_val_2.fq.gz"

  if [[ -f "bam/${sample}.filtered.bam" ]]; then
    echo "[$sample] Already done, skipping"
    continue
  fi

  echo "[$sample] Aligning..."
  bowtie2 --very-sensitive -X 2000 -k 10 -p "$THREADS" -x "$INDEX" \
    --rg-id "$sample" --rg "SM:$sample" \
    -1 "$r1" -2 "$r2" 2>"logs/${sample}.bowtie2.log" \
    | samtools view -bS -@ 2 - > "bam/${sample}.bam"

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
  samtools view -h -q 30 -F 1024 "bam/${sample}.dedup.bam" \
    | awk 'BEGIN{OFS="\t"} $1 ~ /^@/ || $3 != "chrM"' \
    | samtools view -b -@ 2 -o "bam/${sample}.filtered.bam"
  samtools index "bam/${sample}.filtered.bam"

  echo "[$sample] Done"
done

echo "All done"
