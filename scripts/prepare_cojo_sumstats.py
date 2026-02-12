#!/usr/bin/env python3
"""
Filter a large GWAS summary statistics file to a target region and
reshape columns for GCTA COJO (--cojo-file).

Input (GWAS catalog pre-SSF, e.g. GCST90027158_buildGRCh38.tsv.gz):
    variant_id, p_value, chromosome, base_pair_location,
    effect_allele, other_allele, effect_allele_frequency,
    beta, standard_error, n_cases, n_controls, ...

Output columns (tab-delimited):
    SNP   A1   A2   freq   b   se   p   N   chr   bp
SNP is set to "chr:bp" to match the plink IDs in this repo (e.g. 6:146811054).

Example:
    python scripts/prepare_cojo_sumstats.py \
        --sumstats projects/human_esr/gwas/GCST90027158_buildGRCh38.tsv.gz \
        --bim projects/human_esr/results/plink/NA12878.ESR1.bim \
        --chr 6 --start 146809107 --end 157129619 \
        --out projects/human_esr/results/gcta/GCST90027158.ESR1.cojo_input.txt
"""

import argparse
import gzip
from pathlib import Path

import pandas as pd


def load_positions_from_bim(bim_path: Path) -> set[int]:
    """Return set of base-pair positions present in the LD reference (.bim)."""
    positions = set()
    with bim_path.open() as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                positions.add(int(parts[3]))
    return positions


def main():
    parser = argparse.ArgumentParser(description="Subset summary stats to region and format for GCTA COJO.")
    parser.add_argument("--sumstats", required=True, type=Path, help="GWAS summary stats TSV[.gz]")
    parser.add_argument("--bim", required=True, type=Path, help="plink .bim file for LD reference")
    parser.add_argument("--chr", required=True, type=int, help="chromosome number (int)")
    parser.add_argument("--start", required=True, type=int, help="start bp (inclusive)")
    parser.add_argument("--end", required=True, type=int, help="end bp (inclusive)")
    parser.add_argument("--out", required=True, type=Path, help="output path (tsv)")
    parser.add_argument("--chunksize", type=int, default=250_000, help="rows per chunk to stream")
    args = parser.parse_args()

    positions = load_positions_from_bim(args.bim)
    if not positions:
        raise SystemExit(f"No positions read from {args.bim}")

    cols_needed = [
        "chromosome",
        "base_pair_location",
        "effect_allele",
        "other_allele",
        "effect_allele_frequency",
        "beta",
        "standard_error",
        "p_value",
        "n_cases",
        "n_controls",
    ]

    out_path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0

    # Write header
    with out_path.open("w") as out_f:
        out_f.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tN\tchr\tbp\n")

    # Stream through chunks
    reader = pd.read_csv(
        args.sumstats,
        sep="\t",
        compression="gzip" if args.sumstats.suffix == ".gz" else "infer",
        usecols=cols_needed,
        chunksize=args.chunksize,
        dtype={
            "chromosome": str,
            "base_pair_location": int,
            "effect_allele": str,
            "other_allele": str,
            "effect_allele_frequency": float,
            "beta": float,
            "standard_error": float,
            "p_value": float,
            "n_cases": float,
            "n_controls": float,
        },
    )

    for chunk in reader:
        # Filter chromosome and interval
        chunk = chunk[chunk["chromosome"].astype(int) == args.chr]
        if chunk.empty:
            continue
        chunk = chunk[
            (chunk["base_pair_location"] >= args.start)
            & (chunk["base_pair_location"] <= args.end)
        ]
        if chunk.empty:
            continue

        # Keep only SNPs present in LD reference positions
        chunk = chunk[chunk["base_pair_location"].isin(positions)]
        if chunk.empty:
            continue

        # Compute N = n_cases + n_controls (fall back to max if missing)
        n_cases = chunk["n_cases"].fillna(0)
        n_controls = chunk["n_controls"].fillna(0)
        N = n_cases + n_controls
        # If all zeros (rare), set to sample size from metadata (487511) as fallback.
        N = N.replace(0, 487511)

        out_df = pd.DataFrame(
            {
                "SNP": chunk["chromosome"].astype(str) + ":" + chunk["base_pair_location"].astype(str),
                "A1": chunk["effect_allele"].str.upper(),
                "A2": chunk["other_allele"].str.upper(),
                "freq": chunk["effect_allele_frequency"],
                "b": chunk["beta"],
                "se": chunk["standard_error"],
                "p": chunk["p_value"],
                "N": N.astype(int),
                "chr": chunk["chromosome"].astype(int),
                "bp": chunk["base_pair_location"].astype(int),
            }
        )

        out_df.to_csv(out_path, sep="\t", header=False, index=False, mode="a")
        written += len(out_df)

    if written == 0:
        raise SystemExit("No overlapping variants found for the specified region and LD reference.")
    print(f"Wrote {written} variants to {out_path}")


if __name__ == "__main__":
    main()
