#!/usr/bin/env bash
# ================== Create Motif database for ERK ==================
# Based on direct targets of ERK, generate a database of relevant
# motifs for those direct targets
# deps:
# usage:
#   - chmod +x create_erk_meme.sh && ./create_erk_meme.sh
# ===================================================================

set -oue pipefail

PROJECT_DIR=".."
TARGETS_FILE="$PROJECT_DIR/5_analysis/direct_erk_targets.txt"
HOCOMOCO_DB="$PROJECT_DIR/5_analysis/meme_databases/HOCOMOCOv11_core_HUMAN_mono.meme"
OUTPUT_DB="$PROJECT_DIR/5_analysis/meme_databases/erk_targets.meme"

echo "Creating ERK Targets MEME database..."

# MEME Header
cat <<EOF > "$OUTPUT_DB"
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

EOF

TARGETS_REGEX=$(sed 's/$/_HUMAN/' "$TARGETS_FILE" | tr '\n' '|' | sed 's/|$//')

awk -v targets="$TARGETS_REGEX" '
    $1 == "MOTIF" && $2 ~ "^(" targets ")" {
        print_flag = 1
        print ""
        print "# Direct ERK Target: " $2
        print $0
        next
    }
    
    $1 == "MOTIF" {
        print_flag = 0
    }
    
    print_flag == 1 {
        print $0
    }
' "$HOCOMOCO_DB" >> "$OUTPUT_DB"

count=$(grep "MOTIF" "$OUTPUT_DB" | wc -l)
echo "Done. Extracted $count motifs to $OUTPUT_DB"
