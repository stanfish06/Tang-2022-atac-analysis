#!/usr/bin/env bash
# ================== Motif analysis  ==================
# Scan motifs of ERK targets in ATAC peak regions
# deps:
#   - meme
# usage:
#   - chmod +x motif-analysis.sh && ./motif-analysis.sh
# =====================================================

set -euo pipefail

PROJECT_DIR=".."
WORK_DIR="$PROJECT_DIR/5_analysis"
PEAKS_DIR="$PROJECT_DIR/4_peaks"
REF_FASTA="$PROJECT_DIR/0_ref/GRCh37.p13.genome.fa"
TARGETS_BED="$WORK_DIR/targets.bed"

DB_HOCOMOCO="$WORK_DIR/meme_databases/HOCOMOCOv11_core_HUMAN_mono.meme"
DB_JASPAR="$WORK_DIR/meme_databases/JASPAR2024_CORE_vertebrates_nr.meme"
DB_ERK="$WORK_DIR/meme_databases/erk_targets.meme"

# Find assigned peaks (+- some distances to target gene TSS)
find "$PEAKS_DIR" -name "*.filtered.narrowPeak" | sort | while read -r peak_file; do
    
    filename=$(basename "$peak_file")
    sample="${filename%_peaks.filtered.narrowPeak}"
    sample="${sample%.filtered.narrowPeak}"
    
    echo "Processing Sample: $sample"
    
    out_dir="$WORK_DIR/$sample"
    mkdir -p "$out_dir"
    
    echo "Intersecting with targets..."
    temp_intersect="$out_dir/temp_intersect.bed"
    target_peaks="$out_dir/target_peaks.bed"
    
    bedtools intersect -a "$peak_file" -b "$TARGETS_BED" -wa -wb > "$temp_intersect"
    
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4 "::" $14, $5, $6}' "$temp_intersect" > "$target_peaks"
    
    peak_count=$(wc -l < "$target_peaks")
    echo "  > Found $peak_count overlaps."
    
    echo "Extracting sequences..."
    target_seqs_genomic="$out_dir/target_genomic.fa"
    target_seqs_labeled="$out_dir/target_labeled.fa"
    bedtools getfasta -fi "$REF_FASTA" -bed "$target_peaks" -name -fo "$target_seqs_genomic"
    cp "$target_seqs_genomic" "$target_seqs_labeled"
    sed -i 's/::chr.*//' "$target_seqs_labeled"
    sed -i 's/::/__/g' "$target_seqs_labeled"
    
    fimo_out="$out_dir/fimo_out"
    mkdir -p "$fimo_out"
    
    echo "Running full FIMO (for quantification)..."
    fimo --oc "$fimo_out/erk_stats" --thresh 1e-4 "$DB_ERK" "$target_seqs_labeled"
    
    echo "Running Per-Gene FIMO (for visualization)..."
    cut -f4 "$TARGETS_BED" | sort | uniq | while read -r gene_name; do
        echo "Processing Gene: $gene_name"
        gene_seqs="$out_dir/target_genomic_${gene_name}.fa"
        
        awk -v pattern="::${gene_name}::" '
            /^>/ {
                if ($0 ~ pattern) {print_flag=1; print $0} else {print_flag=0}
                next
            }
            print_flag {print $0}
        ' "$target_seqs_genomic" > "$gene_seqs"
        
        if [[ ! -s "$gene_seqs" ]]; then
            echo "No sequences found for $gene_name. Skipping."
            rm "$gene_seqs"
            continue
        fi
        
        fimo --oc "$fimo_out/erk_viz_${gene_name}" --thresh 1e-4 "$DB_ERK" "$gene_seqs"
        rm "$gene_seqs"
    done
    
    echo "Done: $out_dir"
done

echo "All done"
