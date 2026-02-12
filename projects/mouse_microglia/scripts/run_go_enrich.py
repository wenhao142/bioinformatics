#!/usr/bin/env python3
"""
GO Biological Process enrichment (Enrichr) for up/down DE gene lists.
Uses shared helpers in scripts/lib/enrich_utils.py for easier reuse across projects.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.append(str(ROOT / "lib"))

from enrich_utils import load_gene_list, run_enrichr  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--up", default="projects/mouse_microglia/results/gse148405_proxy/up_genes.txt")
    ap.add_argument("--down", default="projects/mouse_microglia/results/gse148405_proxy/down_genes.txt")
    ap.add_argument("--outdir", default="projects/mouse_microglia/results/gse148405_proxy/go_bp")
    ap.add_argument("--library", default="GO_Biological_Process_2023")
    ap.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff")
    ap.add_argument("--topn", type=int, default=10, help="Top N terms for barplot")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    up_genes = load_gene_list(args.up)
    down_genes = load_gene_list(args.down)

    _, sig_up = run_enrichr(
        up_genes,
        args.library,
        outdir / "up_go_bp",
        fdr_cutoff=args.fdr,
        topn=args.topn,
        title="GO BP (Up)",
    )
    _, sig_down = run_enrichr(
        down_genes,
        args.library,
        outdir / "down_go_bp",
        fdr_cutoff=args.fdr,
        topn=args.topn,
        title="GO BP (Down)",
    )

    print(f"GO BP enrichment done. Significant terms (FDR<{args.fdr}): up={len(sig_up)}, down={len(sig_down)}")


if __name__ == "__main__":
    main()
