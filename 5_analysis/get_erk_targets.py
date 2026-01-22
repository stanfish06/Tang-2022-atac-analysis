import pandas as pd
import numpy as np

# Fetched from https://sys-bio.net/erk_targets/ on 2026-01-22
dat = pd.read_csv("./A_Compendium_of_ERK_targets.csv", index_col=0)

mask = np.logical_and(dat["class"] == "direct", dat["organism"].str.contains("human"))
targets = dat[mask]["Gene Name"].unique()

with open("direct_human_erk_targets.txt", "w") as f:
    for gene in targets:
        f.write(gene + "\n")
