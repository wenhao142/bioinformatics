#!/usr/bin/env python3
"""
Download a GEO series with GEOparse and export:
- expression matrix (probes x samples) as TSV.gz
- phenotype table (sample metadata) as TSV

Usage:
  python scripts/py_fetch_geo.py GSE_ID out_dir
"""
import sys
import os
import gzip
import pandas as pd
import GEOparse


def main():
    if len(sys.argv) < 3:
        print("Usage: py_fetch_geo.py GSE_ID out_dir", file=sys.stderr)
        sys.exit(1)
    gse_id, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    gse = GEOparse.get_GEO(geo=gse_id, destdir=out_dir, annotate_gpl=True, how="full")

    # Expression matrix: probes x samples (VALUE)
    expr = gse.pivot_samples("VALUE")
    expr.to_csv(os.path.join(out_dir, "expression.tsv.gz"), sep="\t", compression="gzip")

    # Phenotype metadata: flatten characteristics
    rows = []
    for gsm_name, gsm in gse.gsms.items():
        row = {"sample": gsm_name}
        md = gsm.metadata
        # keep title and source_name_ch1 if available
        for key in ["title", "source_name_ch1"]:
            if key in md:
                row[key] = ";".join(md[key])
        # characteristics_ch1 fields
        for key in md:
            if key.startswith("characteristics_ch1"):
                vals = md[key]
                # try to split "label: value"
                parsed = []
                for v in vals:
                    if ":" in v:
                        k2, v2 = v.split(":", 1)
                        row[k2.strip()] = v2.strip()
                    else:
                        parsed.append(v)
                if parsed:
                    row[key] = ";".join(parsed)
        rows.append(row)
    pheno = pd.DataFrame(rows)
    pheno.to_csv(os.path.join(out_dir, "pheno.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
