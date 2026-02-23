import math
from collections import defaultdict
from typing import Any

from fastapi import APIRouter, Depends, Query

from api import omics, variants
from api.auth import current_user

router = APIRouter(prefix="/evidence", tags=["evidence"])


def _effective_pvalue(expr_row: dict[str, Any] | None) -> float:
    if not expr_row:
        return 1.0
    adj_pval = expr_row.get("adj_pval")
    if isinstance(adj_pval, (int, float)) and adj_pval > 0:
        return float(adj_pval)
    pval = expr_row.get("pval")
    if isinstance(pval, (int, float)) and pval > 0:
        return float(pval)
    return 1.0


def _expr_index_by_gene() -> dict[str, dict[str, Any]]:
    by_gene: dict[str, dict[str, Any]] = {}
    for row in omics._EXPR_TABLE:
        gene = row["gene"].strip()
        if not gene:
            continue
        key = gene.upper()
        existing = by_gene.get(key)
        if existing is None:
            by_gene[key] = dict(row)
            continue
        existing_p = _effective_pvalue(existing)
        new_p = _effective_pvalue(row)
        if new_p < existing_p:
            by_gene[key] = dict(row)
            continue
        if new_p == existing_p and abs(row["logfc"]) > abs(existing["logfc"]):
            by_gene[key] = dict(row)
    return by_gene


def _score_gene_bucket(
    gene: str,
    bucket: list[dict[str, Any]],
    expr_row: dict[str, Any] | None,
) -> dict[str, Any]:
    qualities = [float(v["qual"]) for v in bucket if v.get("qual") is not None]
    positions = [int(v["pos"]) for v in bucket]
    max_qual = max(qualities) if qualities else 0.0
    mean_qual = (sum(qualities) / len(qualities)) if qualities else 0.0

    logfc = expr_row["logfc"] if expr_row else None
    pval = expr_row["pval"] if expr_row else None
    adj_pval = expr_row["adj_pval"] if expr_row else None
    effective_p = _effective_pvalue(expr_row) if expr_row else None

    variant_burden_component = min(3.0, math.log10(len(bucket) + 1) * 2.0)
    variant_quality_component = min(2.0, max_qual / 30.0) if max_qual > 0 else 0.0
    omics_effect_component = min(3.0, abs(float(logfc)) * 1.5) if logfc is not None else 0.0
    omics_significance_component = (
        min(4.0, -math.log10(max(effective_p, 1e-12))) if effective_p is not None else 0.0
    )
    total_score = (
        variant_burden_component
        + variant_quality_component
        + omics_effect_component
        + omics_significance_component
    )

    return {
        "gene": gene,
        "score": round(total_score, 4),
        "locus": {
            "chr": bucket[0]["chr"],
            "start": min(positions),
            "end": max(positions),
        },
        "feature_breakdown": {
            "variant_count": len(bucket),
            "variant_quality": {
                "mean_qual": round(mean_qual, 4) if qualities else None,
                "max_qual": round(max_qual, 4) if qualities else None,
            },
            "omics": {
                "logfc": logfc,
                "pval": pval,
                "adj_pval": adj_pval,
                "effective_pval": effective_p,
            },
            "components": {
                "variant_burden": round(variant_burden_component, 4),
                "variant_quality": round(variant_quality_component, 4),
                "omics_effect": round(omics_effect_component, 4),
                "omics_significance": round(omics_significance_component, 4),
                "total": round(total_score, 4),
            },
        },
    }


def compute_ranked_evidence(chr: str, start: int, end: int, top_n: int) -> dict[str, Any]:
    region_variants = variants.STORE.query(chr, start, end)
    expr_index = _expr_index_by_gene()

    buckets: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in region_variants:
        gene = variants.nearest_gene(record["chr"], record["pos"])
        if not gene:
            continue
        item = dict(record)
        item["gene"] = gene
        buckets[gene].append(item)

    ranked_genes: list[dict[str, Any]] = []
    for gene, bucket in buckets.items():
        ranked_genes.append(_score_gene_bucket(gene, bucket, expr_index.get(gene.upper())))

    ranked_genes.sort(
        key=lambda row: (
            -row["score"],
            -row["feature_breakdown"]["variant_count"],
            row["gene"],
        )
    )
    ranked_genes = ranked_genes[:top_n]

    for idx, row in enumerate(ranked_genes, start=1):
        row["rank"] = idx

    ranked_loci = [
        {
            "rank": row["rank"],
            "gene": row["gene"],
            "chr": row["locus"]["chr"],
            "start": row["locus"]["start"],
            "end": row["locus"]["end"],
            "score": row["score"],
            "feature_breakdown": row["feature_breakdown"],
        }
        for row in ranked_genes
    ]

    return {
        "region": {"chr": chr, "start": start, "end": end},
        "ranked_genes": ranked_genes,
        "ranked_loci": ranked_loci,
        "meta": {
            "variants_considered": len(region_variants),
            "genes_ranked": len(ranked_genes),
            "omics_rows_loaded": len(omics._EXPR_TABLE),
        },
    }


def _aggregate_genes() -> list[dict[str, Any]]:
    expr_index = _expr_index_by_gene()
    buckets: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for record in variants.STORE.all_variants():
        gene = variants.nearest_gene(record["chr"], record["pos"])
        if not gene:
            continue
        buckets[gene].append(record)

    rows: list[dict[str, Any]] = []
    for gene, bucket in buckets.items():
        scored = _score_gene_bucket(gene, bucket, expr_index.get(gene.upper()))
        rows.append(
            {
                "gene": gene,
                "variant_count": len(bucket),
                "score": scored["score"],
                "mean_abs_logfc": (
                    abs(float(scored["feature_breakdown"]["omics"]["logfc"]))
                    if scored["feature_breakdown"]["omics"]["logfc"] is not None
                    else None
                ),
            }
        )

    rows.sort(key=lambda row: (-row["score"], -row["variant_count"], row["gene"]))
    return rows


@router.get("/rank")
def rank_evidence(
    chr: str = Query(..., description="chromosome, e.g., chr1"),
    start: int = Query(1),
    end: int = Query(250_000_000),
    top_n: int = Query(20, ge=1, le=200),
    user=Depends(current_user),
):
    return compute_ranked_evidence(chr=chr, start=start, end=end, top_n=top_n)


@router.get("/aggregate")
def aggregate_evidence(user=Depends(current_user)):
    return {"genes": _aggregate_genes()}
