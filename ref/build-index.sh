#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
threads=8

cd "$script_dir"

genome_gz="GRCh37.p13.genome.fa.gz"
genome_fa="GRCh37.p13.genome.fa"
index_prefix="bowtie2_index/GRCh37"

if [[ ! -f "$genome_fa" ]]; then
  echo "Decompressing $genome_gz..."
  gunzip -k "$genome_gz"
fi

mkdir -p "$(dirname "$index_prefix")"

bowtie2-build \
  --threads "$threads" \
  "$genome_fa" \
  "$index_prefix"
