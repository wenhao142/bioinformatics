from datetime import datetime, timezone
from html import escape
from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query, Response

from api import causal, runs
from api.auth import current_user

router = APIRouter(prefix="/report", tags=["report"])


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


def _markdown_report(project_id: str, run_record: dict[str, Any], result: dict[str, Any]) -> str:
    genes, evidence_rows = _rows_from_result(run_record, result)
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

    lines.extend(
        [
            "",
            "## Evidence Table",
            "| Item | Chr | Start/Pos | End/Alt | Score |",
            "| --- | --- | --- | --- | ---: |",
        ]
    )

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
            "## Method",
            f"- Scoring engine: `{run_record.get('tool_versions', {}).get('scoring', 'unknown')}`",
            "- Inference label: `inference`",
            "",
            "## Reproducibility Metadata",
            f"- Params: `{run_record.get('params', {})}`",
            f"- Input hashes: `{run_record.get('input_hashes', {})}`",
            f"- Tool versions: `{run_record.get('tool_versions', {})}`",
            f"- Signature hash: `{run_record.get('signature_hash', 'NA')}`",
            f"- Result hash: `{run_record.get('result_hash', 'NA')}`",
            f"- Stable with previous: `{run_record.get('stable_with_previous', False)}`",
            f"- Rerun of: `{run_record.get('rerun_of')}`",
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
    user=Depends(current_user),
):
    run_record, result = _resolve_run(project_id, run_id, user)
    if format == "html":
        return Response(content=_html_report(project_id, run_record, result), media_type="text/html")
    return Response(content=_markdown_report(project_id, run_record, result), media_type="text/markdown")
