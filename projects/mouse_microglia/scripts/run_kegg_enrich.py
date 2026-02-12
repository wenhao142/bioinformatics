#!/usr/bin/env python3
"""
KEGG enrichment (Enrichr) for a gene list; reusable across projects.
Default gene list: up-regulated genes from GSE148405 proxy run.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.append(str(ROOT / "lib"))

from enrich_utils import load_gene_list, run_enrichr  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default="projects/mouse_microglia/results/gse148405_proxy/up_genes.txt")
    ap.add_argument("--outdir", default="projects/mouse_microglia/results/gse148405_proxy/kegg")
    ap.add_argument("--library", default="KEGG_2019_Mouse")
    ap.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff")
    ap.add_argument("--topn", type=int, default=10, help="Top N terms for barplot")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    genes = load_gene_list(args.genes)

    raw, sig = run_enrichr(
        genes,
        args.library,
        outdir / "kegg_enrichr",
        fdr_cutoff=args.fdr,
        topn=args.topn,
        title=f"KEGG ({args.library})",
    )

    print(f"KEGG Enrichr done. Significant terms (FDR<{args.fdr}): {len(sig)}")


if __name__ == "__main__":
    main()
