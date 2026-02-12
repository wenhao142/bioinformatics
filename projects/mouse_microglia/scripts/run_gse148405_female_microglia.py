#!/usr/bin/env python3
"""
Draft pipeline steps (placeholders) for GSE148405 female microglia DE.
This file is a skeleton; fill download paths once data is available.
"""
import argparse
import scanpy as sc
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Input h5ad (counts + metadata)")
    parser.add_argument("--sex_col", default="sex")
    parser.add_argument("--genotype_col", default="genotype")
    parser.add_argument("--case", default="5xFAD")
    parser.add_argument("--ctrl", default="WT")
    parser.add_argument("--out_prefix", default="results/microglia_5xFAD_female")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)

    # Filter female
    adata = adata[adata.obs[args.sex_col] == "female"].copy()

    # Basic QC (placeholder)
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # If already microglia only, skip; otherwise subset by annotation if available.
    if "cell_type" in adata.obs:
        adata = adata[adata.obs["cell_type"].str.contains("micro", case=False, na=False)].copy()

    # DE case vs ctrl
    adata.obs["group"] = adata.obs[args.genotype_col]
    sc.tl.rank_genes_groups(adata, "group", groups=[args.case], reference=args.ctrl, method="wilcoxon")
    de = sc.get.rank_genes_groups_df(adata, group=args.case)
    de.to_csv(f"{args.out_prefix}_DE.tsv", sep="\t", index=False)

    # Save top genes lists
    de[(de["pvals_adj"] < 0.05) & (de["logfoldchanges"] > 1)]["names"].to_csv(
        f"{args.out_prefix}_up.txt", index=False, header=False
    )
    de[(de["pvals_adj"] < 0.05) & (de["logfoldchanges"] < -1)]["names"].to_csv(
        f"{args.out_prefix}_down.txt", index=False, header=False
    )

    print("DE written to", f"{args.out_prefix}_DE.tsv")


if __name__ == "__main__":
    main()
