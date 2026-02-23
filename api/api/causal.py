import hashlib
import json
import math
import platform
import time
from importlib import metadata
from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query

from api import omics, variants
from api.auth import current_user

router = APIRouter(prefix="/causal", tags=["causal"])

SCORING_VERSION = "causal-v1"
_CAUSAL_RUNS: list[dict[str, Any]] = []
_CAUSAL_RESULTS: dict[str, dict[str, Any]] = {}


def _sha256_json(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _safe_package_version(package: str) -> str:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return "unknown"


def _tool_versions() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "fastapi": _safe_package_version("fastapi"),
        "pydantic": _safe_package_version("pydantic"),
        "pyjwt": _safe_package_version("PyJWT"),
        "scoring": SCORING_VERSION,
    }


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
    index: dict[str, dict[str, Any]] = {}
    for row in omics.get_expr_table():
        gene = row["gene"].strip().upper()
        if not gene:
            continue
        existing = index.get(gene)
        if existing is None:
            index[gene] = row
            continue
        if _effective_pvalue(row) < _effective_pvalue(existing):
            index[gene] = row
    return index


def _expr_component(expr_row: dict[str, Any] | None) -> float:
    if not expr_row:
        return 0.0
    logfc = float(expr_row.get("logfc", 0.0))
    p_signal = min(1.0, -math.log10(max(_effective_pvalue(expr_row), 1e-12)) / 8.0)
    effect_signal = min(1.0, abs(logfc) / 3.0)
    return 0.6 * effect_signal + 0.4 * p_signal


def _annotation_weight(record: dict[str, Any], gene: str | None) -> float:
    weight = 1.0
    if str(record.get("filter", "")).upper() != "PASS":
        weight -= 0.15
    if gene:
        weight += 0.05
    ref = str(record.get("ref", ""))
    alt = str(record.get("alt", ""))
    if len(ref) == 1 and len(alt) == 1:
        weight += 0.05
    return max(0.6, min(1.2, weight))


def _compute_scores(
    chr: str,
    start: int,
    end: int,
    top_n: int,
    ld_window_bp: int,
) -> dict[str, Any]:
    region_variants = variants.STORE.query(chr, start, end)
    if not region_variants:
        raise HTTPException(status_code=400, detail="No variants in region")

    expr_index = _expr_index_by_gene()
    lead = max(
        region_variants,
        key=lambda row: float(row["qual"]) if row.get("qual") is not None else 0.0,
    )
    lead_pos = int(lead["pos"])

    variant_scores: list[dict[str, Any]] = []
    gene_buckets: dict[str, list[dict[str, Any]]] = {}
    for row in region_variants:
        pos = int(row["pos"])
        qual = float(row["qual"]) if row.get("qual") is not None else 0.0
        qual_signal = min(1.0, max(0.0, qual / 60.0))
        distance = abs(pos - lead_pos)
        ld_signal = max(0.0, 1.0 - (distance / max(1, ld_window_bp)))
        gene = variants.nearest_gene(row["chr"], pos)
        expr_signal = _expr_component(expr_index.get(gene.upper())) if gene else 0.0
        ann_weight = _annotation_weight(row, gene)
        raw_score = (0.45 * qual_signal + 0.35 * ld_signal + 0.2 * expr_signal) * ann_weight
        score = round(raw_score, 6)

        variant_item = {
            "chr": row["chr"],
            "pos": pos,
            "ref": row["ref"],
            "alt": row["alt"],
            "filter": row["filter"],
            "qual": row["qual"],
            "gene": gene,
            "ld_proxy_signal": round(ld_signal, 6),
            "annotation_weight": round(ann_weight, 6),
            "expr_signal": round(expr_signal, 6),
            "score": score,
        }
        variant_scores.append(variant_item)
        if gene:
            gene_buckets.setdefault(gene, []).append(variant_item)

    variant_scores.sort(key=lambda item: (-item["score"], item["pos"], item["alt"]))
    variant_scores = variant_scores[:top_n]

    gene_scores: list[dict[str, Any]] = []
    for gene, bucket in gene_buckets.items():
        bucket_sorted = sorted(bucket, key=lambda item: item["score"], reverse=True)
        top_scores = [item["score"] for item in bucket_sorted[:3]]
        mean_top3 = sum(top_scores) / len(top_scores)
        expr_row = expr_index.get(gene.upper())
        expr_signal = _expr_component(expr_row)
        burden_signal = min(0.2, 0.05 * len(bucket))
        gene_score = 0.7 * mean_top3 + 0.3 * expr_signal + burden_signal
        gene_scores.append(
            {
                "gene": gene,
                "score": round(gene_score, 6),
                "variant_count": len(bucket),
                "top_variant_score": round(bucket_sorted[0]["score"], 6),
                "mean_top3_variant_score": round(mean_top3, 6),
                "expr_signal": round(expr_signal, 6),
            }
        )

    gene_scores.sort(key=lambda item: (-item["score"], -item["variant_count"], item["gene"]))
    gene_scores = gene_scores[:top_n]
    for index, row in enumerate(gene_scores, start=1):
        row["rank"] = index

    return {
        "region": {"chr": chr, "start": start, "end": end},
        "lead_variant": {"chr": lead["chr"], "pos": int(lead["pos"]), "qual": lead["qual"]},
        "variant_scores": variant_scores,
        "gene_scores": gene_scores,
        "meta": {
            "ld_window_bp": ld_window_bp,
            "variants_considered": len(region_variants),
            "genes_scored": len(gene_scores),
            "omics_rows_loaded": len(omics.get_expr_table()),
            "inference_label": "inference",
        },
    }


def _latest_same_signature(signature_hash: str) -> dict[str, Any] | None:
    for run in reversed(_CAUSAL_RUNS):
        if run["signature_hash"] == signature_hash:
            return run
    return None


@router.post("/score")
def run_causal_scoring(
    chr: str = Query(..., description="chromosome, e.g., chr1"),
    start: int = Query(1),
    end: int = Query(250_000_000),
    top_n: int = Query(20, ge=1, le=200),
    ld_window_bp: int = Query(50_000, ge=100, le=1_000_000),
    method: str | None = Query(default=None, max_length=120),
    project_id: str | None = Query(default=None),
    user=Depends(current_user),
):
    params = {
        "chr": chr,
        "start": start,
        "end": end,
        "top_n": top_n,
        "ld_window_bp": ld_window_bp,
        "project_id": project_id,
    }
    region_variants = variants.STORE.query(chr, start, end)
    if not region_variants:
        raise HTTPException(status_code=400, detail="No variants in region")
    expr_rows = omics.get_expr_table()
    input_hashes = {
        "variants_region_sha256": _sha256_json(region_variants),
        "omics_table_sha256": _sha256_json(expr_rows),
    }
    signature_hash = _sha256_json({"params": params, "input_hashes": input_hashes, "scoring": SCORING_VERSION})

    result = _compute_scores(chr=chr, start=start, end=end, top_n=top_n, ld_window_bp=ld_window_bp)
    result_hash = _sha256_json(result)

    previous = _latest_same_signature(signature_hash)
    stable_with_previous = bool(previous and previous["result_hash"] == result_hash)

    run_id = f"causal-{len(_CAUSAL_RUNS) + 1}"
    run_record = {
        "run_id": run_id,
        "kind": "causal_score",
        "created_at": time.time(),
        "created_by": user["email"],
        "project_id": project_id,
        "selected_method": method or "causal-score",
        "params": params,
        "input_hashes": input_hashes,
        "tool_versions": _tool_versions(),
        "signature_hash": signature_hash,
        "result_hash": result_hash,
        "stable_with_previous": stable_with_previous,
        "rerun_of": previous["run_id"] if stable_with_previous else None,
        "counts": {
            "variants_in_region": len(region_variants),
            "omics_rows": len(expr_rows),
            "variant_scores": len(result["variant_scores"]),
            "gene_scores": len(result["gene_scores"]),
        },
        "inference_label": "inference",
    }
    _CAUSAL_RUNS.append(run_record)
    _CAUSAL_RESULTS[run_id] = result

    return {"run": run_record, "result": result}


@router.get("/runs")
def list_causal_runs(limit: int = Query(50, ge=1, le=500), user=Depends(current_user)):
    visible = _CAUSAL_RUNS if user["role"] == "admin" else [r for r in _CAUSAL_RUNS if r["created_by"] == user["email"]]
    return {"runs": list(reversed(visible))[:limit]}


@router.get("/runs/{run_id}")
def get_causal_run(run_id: str, user=Depends(current_user)):
    for run in _CAUSAL_RUNS:
        if run["run_id"] != run_id:
            continue
        if user["role"] == "admin" or run["created_by"] == user["email"]:
            return {"run": run}
        raise HTTPException(status_code=403, detail="Forbidden")
    raise HTTPException(status_code=404, detail="Run not found")


@router.get("/runs/{run_id}/result")
def get_causal_run_result(run_id: str, user=Depends(current_user)):
    for run in _CAUSAL_RUNS:
        if run["run_id"] != run_id:
            continue
        if user["role"] != "admin" and run["created_by"] != user["email"]:
            raise HTTPException(status_code=403, detail="Forbidden")
        result = _CAUSAL_RESULTS.get(run_id)
        if result is None:
            raise HTTPException(status_code=404, detail="Run result not found")
        return {"run_id": run_id, "result": result}
    raise HTTPException(status_code=404, detail="Run not found")
