#!/usr/bin/env python3
"""
Quick differential expression for GSE148405 microglia (5xFAD vs WT).
Assumptions:
- All cells are microglia (sorted in the study).
- Age ~30 weeks (middle-aged; used as menopause proxy per user request).
- Sex metadata is absent in GEO; marked as 'unknown'.
- Genotype/region derived from the Sample_title line of the series matrix.
Outputs: DE table + up/down gene lists + quick volcano plot.
"""
import argparse
import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np


def parse_sample_titles(series_matrix_path):
    """Return list of sample titles in column order."""
    with gzip.open(series_matrix_path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                parts = line.strip().split("\t")[1:]  # drop the key
                # strip surrounding quotes
                return [p.strip('"') for p in parts]
    raise RuntimeError("!Sample_title line not found in series matrix")


def build_metadata(sample_titles):
    rows = []
    for title in sample_titles:
        # Example: P211_A1_5xFAD_C57BL/6J_hippocampus_rep1
        toks = title.split("_")
        sample_id = "_".join(toks[:2])
        if "blank" in title.lower():
            genotype = None
            region = None
        else:
            genotype = "5xFAD" if "5xFAD" in toks else "WT"
            region = (
                "hippocampus"
                if "hippocampus" in toks
                else "cortex"
                if "cortex" in toks
                else "cerebellum"
                if "cerebellum" in toks
                else "unknown"
            )
        rows.append(
            {
                "cell_id": sample_id,
                "sample_title": title,
                "genotype": genotype,
                "region": region,
                "age_weeks": 30,
                "sex": "unknown",
            }
        )
    meta = pd.DataFrame(rows).set_index("cell_id")
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, help="GSE148405_counts.csv.gz")
    ap.add_argument("--series", required=True, help="GSE148405_series_matrix.txt.gz")
    ap.add_argument(
        "--outdir",
        default="projects/mouse_microglia/results/gse148405_proxy",
        help="Output directory",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Reading series matrix for metadata...")
    titles = parse_sample_titles(args.series)
    meta = build_metadata(titles)

    print("Reading counts matrix...")
    counts = pd.read_csv(args.counts, index_col=0)
    if counts.shape[1] != len(titles):
        raise ValueError(f"Counts columns {counts.shape[1]} != titles {len(titles)}")

    counts.columns = meta.index  # ensure names align

    # Drop blanks / missing genotype
    keep_cols = meta[meta["genotype"].notna()].index
    counts = counts[keep_cols]
    meta = meta.loc[keep_cols]

    # Drop ERCC spike-ins
    counts = counts[~counts.index.str.startswith("ERCC")]

    # Build AnnData
    adata = sc.AnnData(counts.T)  # cells x genes
    adata.obs = meta.copy()
    adata.var_names = counts.index

    # Basic QC
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_genes(adata, min_cells=20)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # DE: 5xFAD vs WT
    sc.tl.rank_genes_groups(
        adata, "genotype", groups=["5xFAD"], reference="WT", method="wilcoxon"
    )
    de = sc.get.rank_genes_groups_df(adata, group="5xFAD")
    de.to_csv(outdir / "DE_5xFAD_vs_WT.tsv", sep="\t", index=False)

    # Up / down lists
    de[(de["pvals_adj"] < 0.05) & (de["logfoldchanges"] > 0.25)]["names"].to_csv(
        outdir / "up_genes.txt", index=False, header=False
    )
    de[(de["pvals_adj"] < 0.05) & (de["logfoldchanges"] < -0.25)]["names"].to_csv(
        outdir / "down_genes.txt", index=False, header=False
    )

    # Quick volcano
    plt.figure(figsize=(6, 5))
    plt.scatter(
        de["logfoldchanges"],
        -np.log10(de["pvals_adj"] + 1e-300),
        s=6,
        c="gray",
        alpha=0.5,
    )
    plt.xlabel("log2 fold change (5xFAD vs WT)")
    plt.ylabel("-log10 adj p")
    plt.axvline(0.25, color="orange", linestyle="--", linewidth=0.8)
    plt.axvline(-0.25, color="orange", linestyle="--", linewidth=0.8)
    plt.axhline(-np.log10(0.05), color="blue", linestyle="--", linewidth=0.8)
    plt.title("GSE148405 microglia DE (proxy female, 30w)")
    plt.tight_layout()
    plt.savefig(outdir / "volcano.png", dpi=200)
    plt.close()

    print("Done. Results at", outdir)


if __name__ == "__main__":
    main()
