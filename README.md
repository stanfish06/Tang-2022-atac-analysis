# Scripts for motif analysis on ATAC data
## Basic workflow
Create environment
------------------
```bash
conda env create -f conda_env.yaml
```

Folder intro
------------
| Name       | Description                                                              |
| ---------- | ------------------------------------------------------------------------ |
| 0_ref      | Download and build bowtie2 index with `build-index.sh`                   |
| 1_fetch    | Fetch fastq data from GEO with `download-data.sh`                        |
| 2_trim     | Trim fastq data with `trim-read.sh`                                      |
| 3_align    | Align trimmed data with `align.sh`                                       |
| 4_peaks    | Run MACS2 to find peaks with `filter-peaks`                              |
| 5_analysis | Generate bigwigs `bam-to-bigwigs.sh` and scan motifs `motif-analysis.sh` |

General instruction
-------------------
- Create and activate conda env for all steps
- Go into each folder and run its scripts inside that folder
- Follow the index order specified in folder names
