import json
import os
import re
import time
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from fastapi import APIRouter, Depends, HTTPException, Query

from api.auth import current_user

router = APIRouter(prefix="/literature", tags=["literature"])

TRUE_VALUES = {"1", "true", "yes", "on"}
_PUBMED_RECORDS: list[dict[str, Any]] = []
_PUBMED_KEYS: set[str] = set()


def _online_mode_enabled() -> bool:
    return (
        os.getenv("ONLINE_MODE", "false").lower() in TRUE_VALUES
        or os.getenv("ENABLE_PUBMED", "false").lower() in TRUE_VALUES
    )


def _http_get_json(url: str) -> dict[str, Any]:
    req = Request(url, headers={"User-Agent": "AD-Evidence-Platform/1.0"})
    with urlopen(req, timeout=12) as resp:
        body = resp.read().decode("utf-8")
    return json.loads(body)


def _extract_year(pubdate: str | None) -> int | None:
    if not pubdate:
        return None
    m = re.search(r"(19|20)\d{2}", pubdate)
    if not m:
        return None
    return int(m.group(0))


def _pubmed_search_ids(term: str, retmax: int) -> list[str]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "pubmed", "retmode": "json", "retmax": str(retmax), "term": term}
    payload = _http_get_json(f"{base}?{urlencode(params)}")
    return payload.get("esearchresult", {}).get("idlist", []) or []


def _pubmed_summaries(pmids: list[str]) -> list[dict[str, Any]]:
    if not pmids:
        return []
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {"db": "pubmed", "retmode": "json", "id": ",".join(pmids)}
    payload = _http_get_json(f"{base}?{urlencode(params)}")
    result = payload.get("result", {})
    uids = result.get("uids", []) or []
    out: list[dict[str, Any]] = []
    for uid in uids:
        row = result.get(uid) or {}
        out.append(
            {
                "pmid": str(uid),
                "title": row.get("title") or "",
                "year": _extract_year(row.get("pubdate")),
            }
        )
    return out


def _store_records(records: list[dict[str, Any]]):
    for row in records:
        key = f"{row['gene']}|{row['pmid']}"
        if key in _PUBMED_KEYS:
            continue
        _PUBMED_KEYS.add(key)
        _PUBMED_RECORDS.append(row)


@router.get("/pubmed/fetch")
def fetch_pubmed_metadata(
    genes: list[str] = Query(..., description="Repeat query param: genes=APP&genes=APOE"),
    disease: str = Query("Alzheimer disease"),
    max_per_gene: int = Query(3, ge=1, le=20),
    user=Depends(current_user),
):
    if not _online_mode_enabled():
        raise HTTPException(status_code=503, detail="PubMed provider disabled in offline mode")

    rows: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    for gene in genes:
        query = f"({gene}[Title/Abstract]) AND ({disease}[Title/Abstract])"
        try:
            pmids = _pubmed_search_ids(query, max_per_gene)
            papers = _pubmed_summaries(pmids)
            for paper in papers:
                rows.append(
                    {
                        "gene": gene,
                        "pmid": paper["pmid"],
                        "title": paper["title"],
                        "year": paper["year"],
                        "source": "pubmed",
                        "fetched_at": time.time(),
                    }
                )
        except Exception as exc:
            errors.append({"gene": gene, "error": str(exc)})

    _store_records(rows)
    return {
        "mode": "online",
        "query": {"genes": genes, "disease": disease, "max_per_gene": max_per_gene},
        "records": rows,
        "errors": errors,
        "stored_total": len(_PUBMED_RECORDS),
    }


@router.get("/pubmed")
def list_pubmed_records(
    gene: str | None = None,
    limit: int = Query(100, ge=1, le=1000),
    user=Depends(current_user),
):
    rows = _PUBMED_RECORDS
    if gene:
        rows = [row for row in rows if row["gene"].lower() == gene.lower()]
    return {"records": rows[:limit], "total": len(rows)}
