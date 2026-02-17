import csv
import io
import statistics
from typing import List, TypedDict

from fastapi import APIRouter, Depends, File, HTTPException, UploadFile

from api.auth import current_user


class ExprRecord(TypedDict):
    gene: str
    logfc: float
    pval: float
    adj_pval: float | None


router = APIRouter(prefix="/omics", tags=["omics"])

_EXPR_TABLE: List[ExprRecord] = []


def _parse_expr_file(file_bytes: bytes) -> List[ExprRecord]:
    try:
        text = file_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:  # pragma: no cover - unlikely
        raise HTTPException(status_code=400, detail=f"Decode error: {exc}")

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    required = {"gene", "logfc", "pval"}
    if not required.issubset(set(h.lower() for h in reader.fieldnames or [])):
        raise HTTPException(status_code=400, detail="Missing headers: gene, logfc, pval")

    records: List[ExprRecord] = []
    for row in reader:
        try:
            gene = row.get("gene") or row.get("Gene") or row.get("GENE")
            logfc = float(row.get("logfc") or row.get("logFC") or row.get("LOGFC"))
            pval = float(row.get("pval") or row.get("p") or row.get("P"))
            adj_raw = row.get("adj_pval") or row.get("padj") or row.get("FDR") or None
            adj_pval = float(adj_raw) if adj_raw not in (None, "", "NA", "NaN") else None
        except (TypeError, ValueError):
            continue
        if gene:
            records.append({"gene": gene, "logfc": logfc, "pval": pval, "adj_pval": adj_pval})
    return records


def get_expr_table() -> List[ExprRecord]:
    # return a shallow copy to avoid external mutation
    return list(_EXPR_TABLE)


@router.post("/expr/upload")
def upload_expr(file: UploadFile = File(...), user=Depends(current_user)):
    data = file.file.read()
    records = _parse_expr_file(data)
    if not records:
        raise HTTPException(status_code=400, detail="No rows parsed")
    _EXPR_TABLE.clear()
    _EXPR_TABLE.extend(records)
    return {"ingested": len(records), "example": _EXPR_TABLE[:3]}


@router.get("/expr")
def list_expr(limit: int = 100, user=Depends(current_user)):
    return {"rows": _EXPR_TABLE[:limit]}


@router.get("/expr/stats")
def expr_stats(user=Depends(current_user)):
    if not _EXPR_TABLE:
        return {"count": 0, "logfc_mean": None, "sig_genes": 0}
    logfcs = [r["logfc"] for r in _EXPR_TABLE]
    sig = [r for r in _EXPR_TABLE if r["adj_pval"] is not None and r["adj_pval"] < 0.05]
    return {
        "count": len(_EXPR_TABLE),
        "logfc_mean": statistics.fmean(logfcs),
        "sig_genes": len(sig),
    }
