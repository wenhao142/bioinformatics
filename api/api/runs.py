import hashlib
import json
import platform
import time
from importlib import metadata
from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query

from api import evidence, omics, variants
from api.auth import current_user

router = APIRouter(prefix="/runs", tags=["runs"])

SCORING_VERSION = "evidence-v1"
_RUNS: list[dict[str, Any]] = []
_RUN_RESULTS: dict[str, dict[str, Any]] = {}


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


def _normalized_region_variants(chr: str, start: int, end: int) -> list[dict[str, Any]]:
    rows = variants.STORE.query(chr, start, end)
    normalized: list[dict[str, Any]] = []
    for row in rows:
        normalized.append(
            {
                "chr": row["chr"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "qual": row["qual"],
                "filter": row["filter"],
                "gene": variants.nearest_gene(row["chr"], row["pos"]),
            }
        )
    normalized.sort(key=lambda item: (item["chr"], item["pos"], item["ref"], item["alt"]))
    return normalized


def _normalized_omics_rows() -> list[dict[str, Any]]:
    rows = [
        {
            "gene": row["gene"],
            "logfc": row["logfc"],
            "pval": row["pval"],
            "adj_pval": row["adj_pval"],
        }
        for row in omics._EXPR_TABLE
    ]
    rows.sort(key=lambda item: (item["gene"].upper(), item["pval"], item["logfc"]))
    return rows


def _find_latest_same_signature(signature_hash: str) -> dict[str, Any] | None:
    for run in reversed(_RUNS):
        if run["signature_hash"] == signature_hash:
            return run
    return None


@router.post("/evidence")
def run_evidence(
    chr: str = Query(..., description="chromosome, e.g., chr1"),
    start: int = Query(1),
    end: int = Query(250_000_000),
    top_n: int = Query(20, ge=1, le=200),
    method: str | None = Query(default=None, max_length=120),
    project_id: str | None = Query(default=None),
    user=Depends(current_user),
):
    params = {"chr": chr, "start": start, "end": end, "top_n": top_n}
    region_variants = _normalized_region_variants(chr, start, end)
    omics_rows = _normalized_omics_rows()

    input_hashes = {
        "variants_region_sha256": _sha256_json(region_variants),
        "omics_table_sha256": _sha256_json(omics_rows),
    }
    signature_hash = _sha256_json({"params": params, "input_hashes": input_hashes, "scoring": SCORING_VERSION})

    result = evidence.compute_ranked_evidence(chr=chr, start=start, end=end, top_n=top_n)
    result_hash = _sha256_json({"ranked_genes": result["ranked_genes"], "ranked_loci": result["ranked_loci"]})

    previous = _find_latest_same_signature(signature_hash)
    stable_with_previous = bool(previous and previous["result_hash"] == result_hash)

    run_id = f"run-{len(_RUNS) + 1}"
    run_record = {
        "run_id": run_id,
        "created_at": time.time(),
        "created_by": user["email"],
        "project_id": project_id,
        "kind": "evidence_rank",
        "selected_method": method or "evidence-rank",
        "params": params,
        "input_hashes": input_hashes,
        "tool_versions": _tool_versions(),
        "signature_hash": signature_hash,
        "result_hash": result_hash,
        "stable_with_previous": stable_with_previous,
        "rerun_of": previous["run_id"] if stable_with_previous else None,
        "counts": {
            "variants_in_region": len(region_variants),
            "omics_rows": len(omics_rows),
        },
    }
    _RUNS.append(run_record)
    _RUN_RESULTS[run_id] = result

    return {
        "run": run_record,
        "result": result,
    }


@router.get("")
def list_runs(limit: int = Query(50, ge=1, le=500), user=Depends(current_user)):
    visible = _RUNS if user["role"] == "admin" else [run for run in _RUNS if run["created_by"] == user["email"]]
    return {"runs": list(reversed(visible))[:limit]}


@router.get("/{run_id}")
def get_run(run_id: str, user=Depends(current_user)):
    for run in _RUNS:
        if run["run_id"] != run_id:
            continue
        if user["role"] == "admin" or run["created_by"] == user["email"]:
            return {"run": run}
        raise HTTPException(status_code=403, detail="Forbidden")
    raise HTTPException(status_code=404, detail="Run not found")


@router.get("/{run_id}/result")
def get_run_result(run_id: str, user=Depends(current_user)):
    for run in _RUNS:
        if run["run_id"] != run_id:
            continue
        if user["role"] != "admin" and run["created_by"] != user["email"]:
            raise HTTPException(status_code=403, detail="Forbidden")
        result = _RUN_RESULTS.get(run_id)
        if result is None:
            raise HTTPException(status_code=404, detail="Run result not found")
        return {"run_id": run_id, "result": result}
    raise HTTPException(status_code=404, detail="Run not found")
