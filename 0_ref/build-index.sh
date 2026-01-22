#!/usr/bin/env bash
# ================== Download reference genome ==================
# download fastq files
# deps:
#   - bowtie2
# usage:
#   - chmod +x build-index.sh && ./build-index.sh
# ===============================================================

set -oue pipefail

echo "Download reference genome from gencode"
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

threads=8
genome_gz="GRCh37.p13.genome.fa.gz"
genome_fa="GRCh37.p13.genome.fa"
index_prefix="bowtie2_index/GRCh37"

if [[ ! -f "$genome_fa" ]]; then
	echo "Decompressing $genome_gz..."
	gunzip -k "$genome_gz"
fi

mkdir -p "$(dirname "$index_prefix")"

echo "Build bowtie2 index for alignment"
bowtie2-build \
	--threads "$threads" \
	"$genome_fa" \
	"$index_prefix"
