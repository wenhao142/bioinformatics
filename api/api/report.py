import json
from datetime import datetime, timezone
from html import escape
from typing import Any

import jwt
from fastapi import APIRouter, Depends, HTTPException, Query, Response
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

from api import causal, research, runs
from api.auth import JWT_SECRET, current_user

router = APIRouter(prefix="/report", tags=["report"])
bearer = HTTPBearer(auto_error=False)


def _format_ts(value: float | int | None) -> str:
    if value is None:
        return "unknown"
    return datetime.fromtimestamp(float(value), tz=timezone.utc).isoformat()


def _ensure_access(run_record: dict[str, Any], user: dict[str, Any]) -> None:
    if user["role"] == "admin":
        return
    if run_record.get("created_by") == user["email"]:
        return
    raise HTTPException(status_code=403, detail="Forbidden")


def _auth_user_from_header_or_query(
    token: str | None,
    credentials: HTTPAuthorizationCredentials | None,
) -> dict[str, Any]:
    if credentials and credentials.credentials:
        return current_user(credentials)
    if token:
        try:
            payload = jwt.decode(token, JWT_SECRET, algorithms=["HS256"])
            return {"email": payload.get("sub"), "role": payload.get("role")}
        except jwt.PyJWTError:
            raise HTTPException(status_code=401, detail="Invalid token")
    raise HTTPException(status_code=401, detail="Not authenticated")


def _resolve_run(project_id: str, run_id: str, user: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    record = next((r for r in runs._RUNS if r["run_id"] == run_id), None)
    result = runs._RUN_RESULTS.get(run_id)
    if record is None:
        record = next((r for r in causal._CAUSAL_RUNS if r["run_id"] == run_id), None)
        result = causal._CAUSAL_RESULTS.get(run_id)
    if record is None or result is None:
        raise HTTPException(status_code=404, detail="Run not found")

    _ensure_access(record, user)
    expected_project = record.get("project_id") or "default"
    if project_id != expected_project:
        raise HTTPException(status_code=404, detail="Project/run pair not found")
    return record, result


def _rows_from_result(run_record: dict[str, Any], result: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if run_record["kind"] == "evidence_rank":
        gene_rows = result.get("ranked_genes", [])
        locus_rows = result.get("ranked_loci", [])
        return gene_rows, locus_rows
    gene_rows = result.get("gene_scores", [])
    variant_rows = result.get("variant_scores", [])
    return gene_rows, variant_rows


def _top_genes(gene_rows: list[dict[str, Any]], limit: int = 5) -> list[str]:
    seen: set[str] = set()
    genes: list[str] = []
    for row in gene_rows:
        gene = str(row.get("gene", "")).strip()
        if not gene:
            continue
        key = gene.upper()
        if key in seen:
            continue
        seen.add(key)
        genes.append(gene)
        if len(genes) >= limit:
            break
    return genes


def _research_section(gene_rows: list[dict[str, Any]]) -> dict[str, Any]:
    genes = _top_genes(gene_rows, limit=5)
    if not genes:
        return {
            "disease": "Alzheimer disease",
            "summary": "No ranked genes available for research-direction synthesis.",
            "hypotheses": [],
            "experiments": [],
            "citations": [],
            "inference_label": "inference",
        }

    disease = "Alzheimer disease"
    top_n = min(5, len(genes))
    hypotheses = [research._build_hypothesis(gene, disease, idx + 1) for idx, gene in enumerate(genes[:top_n])]
    experiments = [research._build_experiment(gene, disease, idx + 1) for idx, gene in enumerate(genes[:top_n])]
    summary_request = research.ResearchSummaryRequest(
        disease=disease,
        genes=genes[:top_n],
        top_n=top_n,
        include_pubmed=True,
    )
    citation_bundle = research._build_citation_bundle(summary_request, genes[:top_n])
    summary_text, cited_ids = research._offline_summary(summary_request, citation_bundle)
    citations = [citation_bundle["citation_lookup"][citation_id] for citation_id in cited_ids]
    return {
        "disease": disease,
        "summary": summary_text,
        "hypotheses": hypotheses,
        "experiments": experiments,
        "citations": citations,
        "inference_label": "inference",
    }


def _markdown_report(project_id: str, run_record: dict[str, Any], result: dict[str, Any]) -> str:
    genes, evidence_rows = _rows_from_result(run_record, result)
    research_block = _research_section(genes)

    meta_payload = {
        "selected_method": run_record.get("selected_method", run_record.get("kind")),
        "params": run_record.get("params", {}),
        "input_hashes": run_record.get("input_hashes", {}),
        "tool_versions": run_record.get("tool_versions", {}),
        "signature_hash": run_record.get("signature_hash", "NA"),
        "result_hash": run_record.get("result_hash", "NA"),
        "stable_with_previous": run_record.get("stable_with_previous", False),
        "rerun_of": run_record.get("rerun_of"),
        "counts": run_record.get("counts", {}),
    }

    lines = [
        "# Analysis Report",
        "",
        "## Summary",
        f"- Project: `{project_id}`",
        f"- Run ID: `{run_record['run_id']}`",
        f"- Kind: `{run_record['kind']}`",
        f"- Created by: `{run_record.get('created_by', 'unknown')}`",
        f"- Created at (UTC): `{_format_ts(run_record.get('created_at'))}`",
        "",
        "## Top Genes",
        "| Rank | Gene | Score |",
        "| --- | --- | ---: |",
    ]

    for index, row in enumerate(genes[:10], start=1):
        rank = row.get("rank", index)
        gene = row.get("gene", "NA")
        score = row.get("score", 0)
        lines.append(f"| {rank} | {gene} | {float(score):.4f} |")

    evidence_heading = "## Top Loci Evidence"
    evidence_schema = "| Item | Chr | Start/Pos | End/Alt | Score |"
    if run_record["kind"] != "evidence_rank":
        evidence_heading = "## Top Variant Evidence"
    lines.extend(["", evidence_heading, evidence_schema, "| --- | --- | --- | --- | ---: |"])

    for row in evidence_rows[:15]:
        if run_record["kind"] == "evidence_rank":
            item = row.get("gene", "NA")
            chr_val = row.get("chr", "NA")
            start_val = row.get("start", "NA")
            end_val = row.get("end", "NA")
            score = row.get("score", 0)
            lines.append(f"| {item} | {chr_val} | {start_val} | {end_val} | {float(score):.4f} |")
        else:
            item = row.get("gene", "NA")
            chr_val = row.get("chr", "NA")
            pos = row.get("pos", "NA")
            alt = row.get("alt", "NA")
            score = row.get("score", 0)
            lines.append(f"| {item} | {chr_val} | {pos} | {alt} | {float(score):.4f} |")

    lines.extend(
        [
            "",
            "## Research Direction (Inference)",
            f"- Disease: `{research_block['disease']}`",
            f"- Inference label: `{research_block['inference_label']}`",
            f"- Summary: {research_block['summary']}",
            "",
            "### Hypotheses",
        ]
    )
    if research_block["hypotheses"]:
        for row in research_block["hypotheses"]:
            lines.append(f"- {row['id']} ({row['gene']}): {row['statement']}")
    else:
        lines.append("- No hypotheses available.")

    lines.extend(["", "### Experiment Plans"])
    if research_block["experiments"]:
        for row in research_block["experiments"]:
            lines.append(f"- {row['id']} ({row['gene']}): {row['title']}")
    else:
        lines.append("- No experiment plans available.")

    lines.extend(["", "### Evidence Citations", "| ID | Type | Label |", "| --- | --- | --- |"])
    if research_block["citations"]:
        for citation in research_block["citations"]:
            lines.append(f"| {citation['id']} | {citation['source_type']} | {citation['label']} |")
    else:
        lines.append("| NA | NA | No citations available |")

    lines.extend(
        [
            "",
            "## Run Metadata",
            "```json",
            json.dumps(meta_payload, indent=2, sort_keys=True),
            "```",
        ]
    )
    return "\n".join(lines) + "\n"


def _html_report(project_id: str, run_record: dict[str, Any], result: dict[str, Any]) -> str:
    markdown = _markdown_report(project_id, run_record, result)
    escaped = escape(markdown)
    return (
        "<html><head><meta charset='utf-8'><title>Analysis Report</title>"
        "<style>body{font-family:Segoe UI,Arial,sans-serif;padding:24px;line-height:1.4}"
        "pre{white-space:pre-wrap;background:#f8fafc;border:1px solid #dbe2ea;padding:16px;border-radius:8px}</style>"
        "</head><body><h1>Analysis Report</h1><pre>"
        + escaped
        + "</pre></body></html>"
    )


@router.get("/{project_id}/{run_id}")
def export_report(
    project_id: str,
    run_id: str,
    format: str = Query("markdown", pattern="^(markdown|html)$"),
    token: str | None = Query(default=None, description="JWT token fallback for direct-open links"),
    credentials: HTTPAuthorizationCredentials | None = Depends(bearer),
):
    user = _auth_user_from_header_or_query(token, credentials)
    run_record, result = _resolve_run(project_id, run_id, user)
    if format == "html":
        return Response(content=_html_report(project_id, run_record, result), media_type="text/html")
    return Response(content=_markdown_report(project_id, run_record, result), media_type="text/markdown")
