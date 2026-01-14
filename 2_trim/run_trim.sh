#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")" 

FASTQ_DIR="$PROJECT_DIR/1_fetch/fastqs"
TRIM_DIR="$PROJECT_DIR/2_trim"

find "$FASTQ_DIR" -name "*_1.fastq*" | sort | while read -r r1; do
    r2="${r1/_1.fastq/_2.fastq}"

    echo "Processing: $(basename "$r1")"
    echo "            $(basename "$r2")"

    ./TrimGalore-0.6.10/trim_galore \
        --nextera \
        --paired \
        --fastqc \
        --length 20 \
        -j 4 \
        --output_dir . \
        "$r1" "$r2"

done
