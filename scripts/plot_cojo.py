#!/usr/bin/env python3
"""
Quick visualizations for GCTA COJO outputs.

Given a *.cojo.cma file, produce:
1) Regional Manhattan: bp vs -log10(pC) (falls back to p).
2) Marginal vs conditional effect scatter: b vs bC (if columns exist).

Usage:
    python scripts/plot_cojo.py --cma projects/human_esr/results/gcta/NA12878.ESR1.cojo.cma \
        --out results/gcta/plots/NA12878.ESR1

Outputs:
    <out>.manhattan.png
    <out>.beta_compare.png

The script is defensive: if the file still contains the placeholder text or
columns are missing, it will print a clear message and exit with code 1.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_cma(path: Path) -> pd.DataFrame:
    txt = path.read_text().strip()
    if not txt or txt.lower().startswith(("placeholder", "todo")):
        raise ValueError("COJO output is a placeholder; run GCTA COJO to generate real results.")
    # Try tab/space separated with header; pandas will infer.
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    return df


def pick_columns(df: pd.DataFrame):
    """
    Map a variety of possible COJO column spellings to canonical keys.
    """
    col_map = {}
    lower_cols = {c.lower(): c for c in df.columns}

    def grab(options, required=False):
        for opt in options:
            if opt in lower_cols:
                return lower_cols[opt]
        if required:
            raise KeyError(f"Missing required column; tried {options}")
        return None

    col_map["bp"] = grab(["bp", "pos"])
    col_map["p_cond"] = grab(["pc", "p_cond", "p_conditional", "pc_cond"])
    col_map["p_marg"] = grab(["p", "pval", "p_marg"])
    col_map["b_cond"] = grab(["bc", "b_cond"])
    col_map["b_marg"] = grab(["b", "beta", "b_marg"])
    return col_map


def plot_manhattan(df: pd.DataFrame, cols: dict, out_png: Path):
    p_col = cols["p_cond"] or cols["p_marg"]
    if p_col is None or cols["bp"] is None:
        raise KeyError("Need bp and p/pC columns to draw the regional Manhattan plot.")
    x = df[cols["bp"]].astype(float)
    pvals = df[p_col].astype(float).replace(0, np.nan)
    y = -np.log10(pvals)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.scatter(x, y, s=10, alpha=0.7, edgecolor="none")
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel(f"-log10({p_col})")
    ax.set_title("COJO conditional signals")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def plot_beta_compare(df: pd.DataFrame, cols: dict, out_png: Path):
    if cols["b_cond"] is None or cols["b_marg"] is None:
        raise KeyError("Need b and bC columns to draw marginal vs conditional effects.")
    x = df[cols["b_marg"]].astype(float)
    y = df[cols["b_cond"]].astype(float)
    lim = np.nanmax(np.abs(pd.concat([x, y]))) * 1.1 or 1

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.scatter(x, y, s=12, alpha=0.7, edgecolor="none")
    ax.plot([-lim, lim], [-lim, lim], color="red", linestyle="--", linewidth=1)
    ax.set_xlabel(f"Marginal effect ({cols['b_marg']})")
    ax.set_ylabel(f"Conditional effect ({cols['b_cond']})")
    ax.set_title("Marginal vs conditional β")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot GCTA COJO .cojo.cma results.")
    parser.add_argument("--cma", required=True, type=Path, help="Path to *.cojo.cma file")
    parser.add_argument("--out", required=True, type=Path, help="Output prefix (pngs)")
    args = parser.parse_args()

    try:
        df = load_cma(args.cma)
        cols = pick_columns(df)
        plot_manhattan(df, cols, args.out.with_suffix(".manhattan.png"))
        try:
            plot_beta_compare(df, cols, args.out.with_suffix(".beta_compare.png"))
        except KeyError as e:
            print(f"[warn] Skip beta comparison: {e}", file=sys.stderr)
        print(f"Saved plots to {args.out.with_suffix('.manhattan.png')}")
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
