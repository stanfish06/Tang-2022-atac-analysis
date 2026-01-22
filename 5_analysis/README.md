# Folder for motif analysis
Scripts intro
-------------
| Name                 | Description                                                                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `bam-to-bigwig.sh`   | Generate bigwigs (for visualization) from the aligned bam files after dedup                                                              |
| `get-erk-targets.py` | Find direct erk targets from the file [A_Compendium_of_ERK_targets.csv](A_Compendium_of_ERK_targets.csv)                                 |
| `create-erk-meme.sh` | Create `.meme` database (store motifs) file for erk targets, used for scan                                                               |
| `motif-analysis.sh`  | Do motif scan in three steps: 1. find +- n kb window for each gene 2. find atac peaks in this window 3. scan motifs in those regions     |
| `plot-heatmap.py`    | After counting number of each motif for each target gene, plot a summary heatmap                                                         |
| `plot-tracks.py`     | plot genome tracks for atac signal and motif positions                                                                                   |


Instruction
-----------
- First, run `bam-to-bigwig.sh` first to generate bigwig files
- Next, run `motif-analysis.sh` to scan motifs
- Finally, run `plot-heatmap.py` and `plot-tracks.py` to generate plots (order does not matter)
- Optionally, rerun `get-erk-targets.py` and/or `create-erk-mem.sh` if motif list needs update
