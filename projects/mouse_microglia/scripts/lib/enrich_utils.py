"""
Reusable enrichment helpers (Enrichr / prerank GSEA) for transcriptome projects.

Key functions:
  - load_gene_list(path)
  - run_enrichr(gene_list, library, out_prefix, fdr_cutoff=0.05, topn=10, title=None)
  - run_prerank(rnk_series, gene_sets, outdir, **kwargs)
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Tuple

import gseapy as gp
import pandas as pd


def load_gene_list(path: str | Path) -> list[str]:
    """Read one gene symbol per line, strip blanks."""
    with open(path) as fh:
        return [g.strip() for g in fh if g.strip()]


def run_enrichr(
    gene_list: Iterable[str],
    library: str,
    out_prefix: Path,
    fdr_cutoff: float = 0.05,
    topn: int = 10,
    title: str | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run Enrichr on a gene list.
    Saves raw and significant tables with suffixes _raw.tsv / _sig.tsv.
    Optionally writes a barplot (topn by adjusted p) when sig exists.
    """
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    title = title or out_prefix.name

    enr = gp.enrichr(
        gene_list=list(gene_list),
        gene_sets=library,
        outdir=str(out_prefix.parent),
        cutoff=1.0,  # handle filtering ourselves
        no_plot=True,
    )
    raw = enr.results

    raw_path = Path(f"{out_prefix}_raw.tsv")
    sig_path = Path(f"{out_prefix}_sig.tsv")
    raw.to_csv(raw_path, sep="\t", index=False)

    sig = raw[raw["Adjusted P-value"] < fdr_cutoff]
    sig.to_csv(sig_path, sep="\t", index=False)

    if not sig.empty:
        top = sig.nsmallest(topn, "Adjusted P-value")
        gp.barplot(df=top, title=title, ofname=f"{out_prefix}_barplot.png", cutoff=fdr_cutoff)

    return raw, sig


def run_prerank(
    rnk_series: pd.Series,
    gene_sets: str,
    outdir: Path,
    min_size: int = 10,
    max_size: int = 500,
    permutations: int = 1000,
    threads: int = 4,
    seed: int = 42,
    tag: str = "prerank",
) -> gp.prerank:
    """
    Run gseapy prerank and return the result object.
    rnk_series: indexed by gene, values = statistic (e.g., logFC).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    rnk_series = rnk_series.sort_values(ascending=False)
    rnk_file = outdir / f"{tag}.rnk"
    rnk_series.to_csv(rnk_file, sep="\t", header=False)

    res = gp.prerank(
        rnk=rnk_series,
        gene_sets=gene_sets,
        outdir=str(outdir),
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutations,
        seed=seed,
        processes=threads,
        format="png",
    )
    return res

