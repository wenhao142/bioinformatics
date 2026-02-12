#!/usr/bin/env python3
"""
Create up/down gene lists from DE table.
Usage: python py_build_signature.py de.tsv logfc_cut fdr_cut out_prefix
"""
import sys
import pandas as pd


def main():
    if len(sys.argv) < 5:
        print("Usage: py_build_signature.py de.tsv logfc_cut fdr_cut out_prefix", file=sys.stderr)
        sys.exit(1)
    de_path, logfc_cut, fdr_cut, out_prefix = sys.argv[1:5]
    logfc_cut = float(logfc_cut)
    fdr_cut = float(fdr_cut)

    de = pd.read_csv(de_path, sep="\t")
    if "feature" not in de.columns:
        raise ValueError("DE table must contain 'feature' column")

    up = de[(de["logFC"] >= logfc_cut) & (de["adj.P.Val"] <= fdr_cut)]["feature"]
    down = de[(de["logFC"] <= -logfc_cut) & (de["adj.P.Val"] <= fdr_cut)]["feature"]

    up.to_csv(f"{out_prefix}.up.txt", index=False, header=False)
    down.to_csv(f"{out_prefix}.down.txt", index=False, header=False)
    print(f"[signature] up: {len(up)} down: {len(down)}", file=sys.stderr)


if __name__ == "__main__":
    main()
