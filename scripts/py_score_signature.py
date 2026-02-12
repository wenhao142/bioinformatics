#!/usr/bin/env python3
"""
Rank-based signature score (up mean rank - down mean rank).
Usage: python py_score_signature.py expr.tsv.gz up.txt down.txt out.tsv
"""
import sys
import pandas as pd
import numpy as np


def main():
    if len(sys.argv) < 5:
        print("Usage: py_score_signature.py expr.tsv.gz up.txt down.txt out.tsv", file=sys.stderr)
        sys.exit(1)
    expr_path, up_path, down_path, out_path = sys.argv[1:5]

    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    up = set(pd.read_csv(up_path, header=None)[0].tolist())
    down = set(pd.read_csv(down_path, header=None)[0].tolist())

    common_up = [g for g in expr.index if g in up]
    common_down = [g for g in expr.index if g in down]
    n_genes = expr.shape[0]

    ranks = expr.rank(axis=0, method="average")

    def score(col):
        up_score = ranks.loc[common_up, col].mean() / n_genes if common_up else 0
        down_score = ranks.loc[common_down, col].mean() / n_genes if common_down else 0
        return up_score - down_score

    scores = pd.Series({s: score(s) for s in expr.columns}, name="signature_score")
    scores.index.name = "sample"
    scores.to_csv(out_path, sep="\t")


if __name__ == "__main__":
    main()
