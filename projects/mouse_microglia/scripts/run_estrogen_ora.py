#!/usr/bin/env python3
"""
Simple over-representation analysis (ORA) for estrogen Hallmark gene sets.
Uses DE table from GSE148405 5xFAD vs WT.
Significant genes: adj p < 0.05 and |logFC| > 0.25 (same as volcano thresholds).
Tests enrichment with Fisher's exact test.
"""
import argparse
from pathlib import Path
import pandas as pd
import gseapy as gp
from scipy.stats import fisher_exact


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--de",
        default="projects/mouse_microglia/results/gse148405_proxy/DE_5xFAD_vs_WT.tsv",
        help="DE table with columns: names, logfoldchanges, pvals_adj",
    )
    ap.add_argument(
        "--out",
        default="projects/mouse_microglia/results/gse148405_proxy/gsea_estrogen/ora_estrogen.tsv",
        help="Output TSV",
    )
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.de, sep="\t")

    # define significant set
    sig = df[(df.pvals_adj < 0.05) & (df.logfoldchanges.abs() > 0.25)]
    gene_universe = set(df.names)
    sig_genes = set(sig.names)

    # pull hallmark gene sets
    # gseapy Hallmark name (mouse) is 'MSigDB_Hallmark_2020' and works for mouse symbols.
    hallmark = gp.get_library(name="MSigDB_Hallmark_2020", organism="Mouse")
    # filter estrogen sets
    estro_sets = {k: v for k, v in hallmark.items() if "ESTROGEN RESPONSE" in k.upper()}

    rows = []
    for term, genes in estro_sets.items():
        gs = set(genes)
        a = len(sig_genes & gs)  # in set & sig
        b = len(gs - sig_genes)  # in set not sig
        c = len(sig_genes - gs)  # sig not in set
        d = len(gene_universe - sig_genes - gs)  # neither
        table = [[a, b], [c, d]]
        _, p = fisher_exact(table, alternative="greater")
        rows.append({"Term": term, "overlap": a, "set_size": len(gs), "pval": p})

    if rows:
        res = pd.DataFrame(rows).sort_values("pval")
        res.to_csv(out_path, sep="\t", index=False)
        print("ORA written to", out_path)
    else:
        print("No estrogen sets found; nothing written.")


if __name__ == "__main__":
    main()
