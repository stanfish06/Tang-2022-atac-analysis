#!/usr/bin/env bash

find ../1_fetch/fastqs -name "*_1.fastq*" | sort | while read -r r1; do
    r2="${r1/_1.fastq/_2.fastq}"
    trim_galore \
        --nextera \
        --paired \
        --fastqc \
        --length 20 \
        -j 4 \
        --output_dir . \
        "$r1" "$r2"
done

echo "done"
