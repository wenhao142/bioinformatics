#!/usr/bin/env python3
"""
Two-group differential expression using Welch's t-test.
Inputs:
  expr.tsv.gz : genes/probes x samples
  pheno.tsv   : must contain 'sample' and group_col

Usage:
  python scripts/py_de_ttest.py expr.tsv.gz pheno.tsv group_col case_label ctrl_label out_prefix
"""
import sys
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy import stats


def main():
    if len(sys.argv) < 7:
        print("Usage: py_de_ttest.py expr.tsv.gz pheno.tsv group_col case_label ctrl_label out_prefix", file=sys.stderr)
        sys.exit(1)
    expr_path, pheno_path, group_col, case, ctrl, out_prefix = sys.argv[1:7]

    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    pheno = pd.read_csv(pheno_path, sep="\t")
    if "sample" not in pheno.columns:
        raise ValueError("pheno.tsv must contain 'sample' column")
    # intersect samples
    common = [s for s in expr.columns if s in set(pheno["sample"])]
    expr = expr[common]
    pheno = pheno.set_index("sample").loc[common]
    if group_col not in pheno.columns:
        raise ValueError(f"group_col {group_col} not found in pheno")

    labels = pheno[group_col].astype(str)
    mask_case = labels == case
    mask_ctrl = labels == ctrl
    if mask_case.sum() == 0 or mask_ctrl.sum() == 0:
        raise ValueError("case or ctrl group has zero samples")

    case_expr = expr.loc[:, mask_case.values]
    ctrl_expr = expr.loc[:, mask_ctrl.values]

    # Welch's t-test per row
    tstat, pvals = stats.ttest_ind(case_expr.values, ctrl_expr.values, axis=1, equal_var=False, nan_policy="omit")
    case_mean = np.nanmean(case_expr.values, axis=1)
    ctrl_mean = np.nanmean(ctrl_expr.values, axis=1)
    logfc = np.log2((case_mean + 1e-8) / (ctrl_mean + 1e-8))

    # FDR
    adj = multipletests(pvals, method="fdr_bh")[1]

    out = pd.DataFrame({
        "feature": expr.index,
        "logFC": logfc,
        "P.Value": pvals,
        "adj.P.Val": adj,
        "case_mean": case_mean,
        "ctrl_mean": ctrl_mean,
        "n_case": mask_case.sum(),
        "n_ctrl": mask_ctrl.sum(),
    })
    out = out.sort_values("adj.P.Val")
    out.to_csv(f"{out_prefix}.de.tsv", sep="\t", index=False)

    # simple volcano data (for plotting later if needed)
    out.to_csv(f"{out_prefix}.volcano.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
