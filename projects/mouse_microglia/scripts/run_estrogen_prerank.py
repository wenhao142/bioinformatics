#!/usr/bin/env python3
"""
Preranked GSEA for estrogen response (Hallmark EARLY/LATE) on 5xFAD vs WT DE.
Input: DE table from run_gse148405_proxy.py (columns: names, logfoldchanges, pvals_adj, etc.)
Uses Enrichr Hallmark gene sets via gseapy (no local MSigDB download needed).
"""
import argparse
from pathlib import Path

import pandas as pd
import gseapy as gp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--de",
        default="projects/mouse_microglia/results/gse148405_proxy/DE_5xFAD_vs_WT.tsv",
        help="DE table with columns: names, logfoldchanges",
    )
    ap.add_argument(
        "--outdir",
        default="projects/mouse_microglia/results/gse148405_proxy/gsea_estrogen",
        help="Output directory",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.de, sep="\t")
    rank_series = df.set_index("names")["logfoldchanges"]
    # sort descending as required by prerank
    rank_series = rank_series.sort_values(ascending=False)
    rank_file = outdir / "prerank.rnk"
    rank_series.to_csv(rank_file, sep="\t", header=False)

    # Run prerank against Hallmark (will include ESTROGEN_RESPONSE_EARLY/LATE)
    prerank_res = gp.prerank(
        rnk=rank_series,
        gene_sets="MSigDB_Hallmark_2020",
        outdir=str(outdir),
        min_size=10,
        max_size=500,
        permutation_num=1000,
        seed=42,
        processes=4,
        format="png",
    )

    # Extract estrogen subsets from the main report CSV that gseapy writes.
    report_path = outdir / "gseapy.gene_set.prerank.report.csv"
    if report_path.exists():
        enr = pd.read_csv(report_path)
        estro = enr[enr["Term"].str.contains("Estrogen Response", case=False)]
        estro.to_csv(outdir / "estrogen_only.tsv", sep="\t", index=False)
    else:
        print("Warning: report file not found; skipping estrogen subset.")

    print("GSEA finished. Key files:")
    print("  Ranked list:", rank_file)
    print("  Full results:", outdir / "gseapy.prerank.gene_sets.report.csv")
    print("  Estrogen subset:", outdir / "estrogen_only.tsv")


if __name__ == "__main__":
    main()
