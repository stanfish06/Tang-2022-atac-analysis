#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob

SAMPLE = "ATAC_hPGCLC_d2_r1"
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIMO_BASE = os.path.join(BASE_DIR, SAMPLE, "fimo_out")
OUTPUT_FILE = os.path.join(BASE_DIR, f"{SAMPLE}_motif_heatmap.png")

Q_THRESHOLD = 1 

gene_dirs = glob.glob(os.path.join(FIMO_BASE, "erk_viz_*"))

data = {}

for d in gene_dirs:
    dirname = os.path.basename(d)
    gene_name = dirname.replace("erk_viz_", "")
    
    fimo_file = os.path.join(d, "fimo.tsv")
    
    if not os.path.exists(fimo_file):
        print(f"Warning: No fimo.tsv in {d}")
        continue
        
    df = pd.read_csv(fimo_file, sep="\t", comment="#")
        
    if df.empty:
        continue
            
    if 'q-value' in df.columns:
        df = df[df['q-value'] <= Q_THRESHOLD]
            
    counts = df['motif_id'].value_counts()
        
    data[gene_name] = counts

matrix = pd.DataFrame(data).fillna(0).T

gene_map = {name: name.replace("_25kb", "") for name in matrix.index}
matrix.rename(index=gene_map, inplace=True)

endo_genes = ["EOMES", "FOXA2", "OTX2"]
pgc_genes = ["NANOG", "PRDM1", "SOX17", "TFAP2C", "KLF4"]
control_genes = ["SPRY1", "SPRY4", "DUSP5", "DUSP6"]

def get_group_label(gene):
    if gene in endo_genes:
        return f"{gene}", "A_Endo"
    elif gene in pgc_genes:
        return f"{gene}", "B_PGC"
    elif gene in control_genes:
        return f"{gene} (Ctrl)", "C_Control"
    else:
        return f"{gene}", "D_Other"

new_index = []
sort_keys = []
for gene in matrix.index:
    label, sort_key = get_group_label(gene)
    new_index.append(label)
    sort_keys.append(sort_key)

matrix.index = new_index
matrix["_sort"] = sort_keys
matrix = matrix.sort_values("_sort")
matrix = matrix.drop(columns=["_sort"])

matrix.columns = [c.split("_")[0] for c in matrix.columns]

matrix = matrix.reindex(sorted(matrix.columns), axis=1)

norm_matrix = matrix.apply(lambda x: x / x.max() if x.max() > 0 else x, axis=0)

plt.figure(figsize=(18, 6)) 
sns.set_theme(style="whitegrid")
sns.heatmap(norm_matrix, annot=matrix, fmt="g", cmap="viridis", linewidths=1, cbar=False)
plt.title(f"Motif Frequency: ±25kb from TSS (colors normalized by columns, unclustered)", fontsize=16)
plt.ylabel("Target Locus", fontsize=12)
plt.xlabel("Motif", fontsize=12)
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=300)

import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch

OUTPUT_FILE_CLUSTERED = os.path.join(BASE_DIR, f"{SAMPLE}_motif_heatmap_control_clustered.png")

control_rows = [r for r in matrix.index if "(Ctrl)" in r]
if control_rows:
    ctrl_data = norm_matrix.loc[control_rows]
    
    ctrl_data = ctrl_data.fillna(0)
    valid_cols = ctrl_data.columns[ctrl_data.sum() > 0]
    
    if len(valid_cols) > 1:
        ctrl_subset = ctrl_data[valid_cols]
        
        d = ssd.pdist(ctrl_subset.T, metric='euclidean')
        L = sch.linkage(d, method='average')
        dendro = sch.dendrogram(L, no_plot=True)
        leaves = dendro['leaves']
        
        sorted_valid = [valid_cols[i] for i in leaves]
        unused = [c for c in matrix.columns if c not in sorted_valid]
        final_order = sorted_valid + unused
        
        norm_matrix_sorted = norm_matrix[final_order]
        matrix_sorted = matrix[final_order]
        
        plt.figure(figsize=(18, 6))
        sns.heatmap(norm_matrix_sorted, annot=matrix_sorted, fmt="g", cmap="viridis", linewidths=1, cbar=False)
        plt.title(f"Motif Frequency (colors normalized by columns, columns clustered using control genes)", fontsize=16)
        plt.ylabel("Target Locus", fontsize=12)
        plt.xlabel("Motif (Sorted by Similarity in Controls)", fontsize=12)
        plt.xticks(rotation=30, ha='right')
        plt.tight_layout()
        plt.savefig(OUTPUT_FILE_CLUSTERED, dpi=300)
