#!/usr/bin/env bash
# ================== Filter atac peaks ==================
# Find peaks in atac data using macs2
# deps:
# usage:
#   - chmod +x filter-peaks.sh && ./filter_peaks.sh
# ======================================================
set -oue pipefail

# Filtering parameters
Q_THRESH=-1         # -log10(q-value) > 4 (q < 0.0001) [q is FDR corrected p value]
MIN_WIDTH=200       # Minimum peak width (bp)
MAX_WIDTH=50000000  # Maximum peak width (bp)

PROJECT_DIR=".."
PEAKS_DIR="$PROJECT_DIR"

find "$PEAKS_DIR" -name "*_peaks.narrowPeak" ! -name "*.filtered.*" | while read -r peak_file; do
    dir=$(dirname "$peak_file")
    base=$(basename "$peak_file" .narrowPeak)
    out_file="$dir/${base}.filtered.narrowPeak"
    
    echo "Processing: $base"
    
    awk -v q="$Q_THRESH" -v minw="$MIN_WIDTH" -v maxw="$MAX_WIDTH" -v sig="0" \
    '\
    BEGIN {OFS="\t"}
    {
        width = $3 - $2;
        if ($9 > q && width >= minw && width <= maxw && $7 > sig) {
            start = $2;
            if (start < 0) start = 0;
            end = $3;
            print $1, start, end, $4, $5, $6, $7, $8, $9, $10
        }
    }' "$peak_file" > "$out_file"
    
    original_count=$(wc -l < "$peak_file")
    filtered_count=$(wc -l < "$out_file")
    echo "Original: $original_count"
    echo "Filtered: $filtered_count"

done

echo "Done"
