#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch
import os
import glob

SAMPLES = [
    # "ATAC_hESC_r1",
    # "ATAC_hESC_r2",
    # "ATAC_hPGC_r1",
    # "ATAC_hPGC_r2",
    "ATAC_hPGCLC_d2_r1",
    "ATAC_hPGCLC_d2_r2",
    # "ATAC_hPGCLC_d4_r1",
    # "ATAC_hPGCLC_d4_r2",
]

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
Q_THRESHOLD = 1

ENDO_GENES = ["EOMES", "FOXA2", "OTX2"]
PGC_GENES = ["NANOG", "PRDM1", "SOX17", "TFAP2C", "KLF4"]
CONTROL_GENES = ["SPRY1", "SPRY4", "DUSP5", "DUSP6"]


def get_group_label(gene: str) -> tuple[str, str]:
    if gene in ENDO_GENES:
        return f"{gene}", "A_Endo"
    elif gene in PGC_GENES:
        return f"{gene}", "B_PGC"
    elif gene in CONTROL_GENES:
        return f"{gene} (Ctrl)", "C_Control"
    else:
        return f"{gene}", "D_Other"


def read_fimo_table(path: str) -> pd.DataFrame:
    header = None
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                possible = stripped.lstrip("#").strip()
                if "\t" in possible:
                    header = possible.split("\t")
                continue
            break

    if header:
        df = pd.read_csv(path, sep="\t", comment="#", header=None, names=header)
    else:
        df = pd.read_csv(path, sep="\t")

    df.columns = [c.lstrip("#").strip() for c in df.columns]
    return df


def process_sample(sample: str) -> None:
    print(f"=== Processing sample: {sample} ===")

    fimo_base = os.path.join(BASE_DIR, sample, "fimo_out")
    output_file = os.path.join(BASE_DIR, f"{sample}_motif_heatmap.png")
    output_file_clustered = os.path.join(
        BASE_DIR, f"{sample}_motif_heatmap_control_clustered.png"
    )

    gene_dirs = glob.glob(os.path.join(fimo_base, "erk_viz_*"))

    if not gene_dirs:
        print(f"  > No FIMO results found for {sample}. Skipping.")
        return

    data = {}

    for d in gene_dirs:
        dirname = os.path.basename(d)
        gene_name = dirname.replace("erk_viz_", "")

        fimo_file = os.path.join(d, "fimo.txt")
        if not os.path.exists(fimo_file):
            fimo_alt = os.path.join(d, "fimo.tsv")
            if os.path.exists(fimo_alt):
                fimo_file = fimo_alt

        if not os.path.exists(fimo_file):
            print(f"  > Warning: No FIMO table in {d}")
            continue

        df = read_fimo_table(fimo_file)

        if df.empty:
            continue

        q_col = None
        for col in ("q-value", "q_value", "qvalue"):
            if col in df.columns:
                q_col = col
                break
        if q_col is not None:
            df = df[df[q_col] <= Q_THRESHOLD]

        motif_col = None
        for col in ("motif_id", "pattern name", "pattern_name"):
            if col in df.columns:
                motif_col = col
                break
        if motif_col is None:
            print(f"  > Warning: No motif column in {fimo_file}")
            continue

        counts = df[motif_col].value_counts()

        data[gene_name] = counts

    if not data:
        print(f"  > No data found for {sample}. Skipping.")
        return

    matrix = pd.DataFrame(data).fillna(0).T

    gene_map = {name: name.replace("_25kb", "") for name in matrix.index}
    matrix.rename(index=gene_map, inplace=True)

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
    sns.heatmap(
        norm_matrix, annot=matrix, fmt="g", cmap="viridis", linewidths=1, cbar=False
    )
    plt.title(
        f"{sample}: Motif Frequency ±25kb from TSS (normalized by columns)", fontsize=16
    )
    plt.ylabel("Target Locus", fontsize=12)
    plt.xlabel("Motif", fontsize=12)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"  > Saved: {output_file}")

    control_rows = [r for r in matrix.index if "(Ctrl)" in r]
    if control_rows:
        ctrl_data = norm_matrix.loc[control_rows]

        ctrl_data = ctrl_data.fillna(0)
        valid_cols = ctrl_data.columns[ctrl_data.sum() > 0]

        if len(valid_cols) > 1:
            ctrl_subset = ctrl_data[valid_cols]

            d = ssd.pdist(ctrl_subset.T, metric="euclidean")
            L = sch.linkage(d, method="average")
            dendro = sch.dendrogram(L, no_plot=True)
            leaves = dendro["leaves"]

            sorted_valid = [valid_cols[i] for i in leaves]
            unused = [c for c in matrix.columns if c not in sorted_valid]
            final_order = sorted_valid + unused

            norm_matrix_sorted = norm_matrix[final_order]
            matrix_sorted = matrix[final_order]

            plt.figure(figsize=(18, 6))
            sns.heatmap(
                norm_matrix_sorted,
                annot=matrix_sorted,
                fmt="g",
                cmap="viridis",
                linewidths=1,
                cbar=False,
            )
            plt.title(
                f"{sample}: Motif Frequency (clustered by control genes)", fontsize=16
            )
            plt.ylabel("Target Locus", fontsize=12)
            plt.xlabel("Motif (Sorted by Similarity in Controls)", fontsize=12)
            plt.xticks(rotation=30, ha="right")
            plt.tight_layout()
            plt.savefig(output_file_clustered, dpi=300)
            plt.close()
            print(f"  > Saved: {output_file_clustered}")

    print(f"=== Done: {sample} ===")


if __name__ == "__main__":
    for sample in SAMPLES:
        process_sample(sample)
    print("All done")
