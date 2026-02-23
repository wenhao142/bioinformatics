import json
import os
from typing import Any
from urllib.request import Request, urlopen

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field

from api import evidence, literature
from api.auth import current_user

router = APIRouter(prefix="/research", tags=["research"])
TRUE_VALUES = {"1", "true", "yes", "on"}


class ResearchDirectionRequest(BaseModel):
    disease: str = Field(default="Alzheimer disease")
    genes: list[str] | None = None
    chr: str | None = None
    start: int | None = None
    end: int | None = None
    top_n: int = Field(default=5, ge=1, le=20)


class ResearchSummaryRequest(BaseModel):
    disease: str = Field(default="Alzheimer disease")
    genes: list[str] | None = None
    chr: str | None = None
    start: int | None = None
    end: int | None = None
    top_n: int = Field(default=5, ge=1, le=20)
    include_pubmed: bool = True
    llm_mode: str = Field(default="auto", pattern="^(auto|offline)$")
    llm_model: str | None = Field(default=None, max_length=120)


def _normalized_genes(genes: list[str] | None) -> list[str]:
    if not genes:
        return []
    seen: set[str] = set()
    out: list[str] = []
    for raw in genes:
        gene = raw.strip()
        if not gene:
            continue
        key = gene.upper()
        if key in seen:
            continue
        seen.add(key)
        out.append(gene)
    return out


def _genes_from_region(request: ResearchDirectionRequest) -> list[str]:
    if not request.chr:
        return []
    start = request.start if request.start is not None else 1
    end = request.end if request.end is not None else 250_000_000
    ranked = evidence.compute_ranked_evidence(chr=request.chr, start=start, end=end, top_n=request.top_n)
    return [row["gene"] for row in ranked["ranked_genes"]]


def _build_hypothesis(gene: str, disease: str, index: int) -> dict[str, Any]:
    return {
        "id": f"H{index}",
        "gene": gene,
        "statement": f"{gene} may modulate {disease} risk through transcriptomic and variant-linked mechanisms.",
        "inference_label": "inference",
        "evidence_basis": "Template/rule-based offline synthesis from ranked locus evidence.",
    }


def _build_experiment(gene: str, disease: str, index: int) -> dict[str, Any]:
    return {
        "id": f"E{index}",
        "gene": gene,
        "title": f"Targeted validation plan for {gene}",
        "rationale": f"Validate whether {gene} perturbation changes AD-relevant pathways and molecular readouts.",
        "steps": [
            f"Replicate {gene} differential expression in an independent {disease} cohort.",
            f"Perform genotype-stratified expression analysis at the {gene} locus.",
            f"Run CRISPR/siRNA perturbation for {gene} in a neural or microglial model.",
        ],
        "primary_readouts": [
            "Expression shift of disease-relevant markers",
            "Amyloid/tau-related pathway activity",
            "Effect size consistency across cohorts",
        ],
        "inference_label": "inference",
    }


def _llm_summary_enabled() -> bool:
    return (
        os.getenv("ENABLE_LLM_SUMMARY", "false").lower() in TRUE_VALUES
        and os.getenv("ONLINE_MODE", "false").lower() in TRUE_VALUES
        and bool(os.getenv("OPENAI_API_KEY"))
    )


def _call_cloud_llm(prompt: str, model_override: str | None = None) -> dict[str, Any]:
    api_key = os.getenv("OPENAI_API_KEY")
    model = model_override or os.getenv("OPENAI_MODEL", "gpt-4o-mini")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY missing")
    payload = {
        "model": model,
        "temperature": 0.2,
        "response_format": {"type": "json_object"},
        "messages": [
            {
                "role": "system",
                "content": (
                    "You generate concise research summaries for biomedical evidence bundles. "
                    "Use only provided citation IDs. Return JSON with keys: summary, citations."
                ),
            },
            {"role": "user", "content": prompt},
        ],
    }
    req = Request(
        "https://api.openai.com/v1/chat/completions",
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        },
    )
    with urlopen(req, timeout=20) as response:
        raw = json.loads(response.read().decode("utf-8"))
    content = raw["choices"][0]["message"]["content"]
    return json.loads(content)


def _build_citation_bundle(request: ResearchSummaryRequest, genes: list[str]) -> dict[str, Any]:
    citation_rows: list[dict[str, Any]] = []
    citation_lookup: dict[str, dict[str, Any]] = {}
    next_id = 1

    def add(source_type: str, label: str, payload: dict[str, Any]) -> str:
        nonlocal next_id
        citation_id = f"C{next_id}"
        next_id += 1
        row = {
            "id": citation_id,
            "source_type": source_type,
            "label": label,
            "payload": payload,
        }
        citation_rows.append(row)
        citation_lookup[citation_id] = row
        return citation_id

    for gene in genes:
        add("input-gene", f"Requested gene {gene}", {"gene": gene})

    ranked: dict[str, Any] | None = None
    if request.chr:
        start = request.start if request.start is not None else 1
        end = request.end if request.end is not None else 250_000_000
        ranked = evidence.compute_ranked_evidence(chr=request.chr, start=start, end=end, top_n=request.top_n)
        for row in ranked["ranked_genes"][: request.top_n]:
            add(
                "ranked-evidence",
                f"{row['gene']} score={row['score']} rank={row['rank']}",
                {"gene": row["gene"], "score": row["score"], "rank": row["rank"]},
            )

    if request.include_pubmed:
        for gene in genes:
            records = [r for r in literature._PUBMED_RECORDS if r["gene"].lower() == gene.lower()]
            for record in records[:3]:
                add(
                    "pubmed",
                    f"PMID {record['pmid']} ({record.get('year') or 'n/a'})",
                    {
                        "gene": gene,
                        "pmid": record["pmid"],
                        "title": record["title"],
                        "year": record.get("year"),
                    },
                )

    if not citation_rows:
        raise HTTPException(status_code=400, detail="No evidence available to summarize")

    return {
        "genes": genes,
        "citations": citation_rows,
        "citation_lookup": citation_lookup,
        "ranked": ranked,
    }


def _offline_summary(request: ResearchSummaryRequest, bundle: dict[str, Any]) -> tuple[str, list[str]]:
    genes = bundle["genes"][: request.top_n]
    ranked = bundle["ranked"]
    if ranked and ranked["ranked_genes"]:
        top = ", ".join(f"{row['gene']} ({row['score']:.2f})" for row in ranked["ranked_genes"][:3])
        summary = (
            f"Offline template summary for {request.disease}: prioritized genes are {top}. "
            "Signals combine locus-level variant evidence with available omics effect sizes."
        )
    else:
        summary = (
            f"Offline template summary for {request.disease}: candidate genes are {', '.join(genes)}. "
            "Interpret as inference pending literature and lab validation."
        )
    chosen = [row["id"] for row in bundle["citations"][: min(6, len(bundle["citations"]))]]
    return summary, chosen


@router.post("/directions")
def generate_research_directions(request: ResearchDirectionRequest, user=Depends(current_user)):
    genes = _normalized_genes(request.genes)
    if not genes:
        genes = _genes_from_region(request)

    if not genes:
        raise HTTPException(
            status_code=400,
            detail="No genes provided and no ranked genes found for the requested region",
        )

    hypotheses = [_build_hypothesis(gene, request.disease, i + 1) for i, gene in enumerate(genes[: request.top_n])]
    experiments = [_build_experiment(gene, request.disease, i + 1) for i, gene in enumerate(genes[: request.top_n])]

    return {
        "mode": "offline-template",
        "disease": request.disease,
        "genes_considered": genes[: request.top_n],
        "hypotheses": hypotheses,
        "experiments": experiments,
        "uncertainty_notes": [
            "These outputs are template/rule-based inferences, not literature-verified causal claims.",
            "Prioritize wet-lab validation and replication before decision-making.",
        ],
    }


@router.post("/summary")
def summarize_research(request: ResearchSummaryRequest, user=Depends(current_user)):
    genes = _normalized_genes(request.genes)
    if not genes:
        genes = _genes_from_region(
            ResearchDirectionRequest(
                disease=request.disease,
                genes=request.genes,
                chr=request.chr,
                start=request.start,
                end=request.end,
                top_n=request.top_n,
            )
        )
    if not genes:
        raise HTTPException(status_code=400, detail="No genes provided and no ranked genes found for the requested region")

    bundle = _build_citation_bundle(request, genes[: request.top_n])
    citation_lookup: dict[str, dict[str, Any]] = bundle["citation_lookup"]
    mode = "offline-template"
    warnings: list[str] = []
    llm_mode_requested = request.llm_mode
    llm_model_used: str | None = None

    summary_text: str
    cited_ids: list[str]
    if request.llm_mode == "offline":
        summary_text, cited_ids = _offline_summary(request, bundle)
    elif _llm_summary_enabled():
        prompt = json.dumps(
            {
                "disease": request.disease,
                "genes": bundle["genes"][: request.top_n],
                "citations": bundle["citations"],
                "instructions": "Use only citation IDs from citations[]. Return JSON {summary, citations:[id,...]}",
            },
            ensure_ascii=True,
        )
        try:
            llm_out = _call_cloud_llm(prompt, request.llm_model)
            summary_text = str(llm_out.get("summary", "")).strip()
            raw_ids = llm_out.get("citations") if isinstance(llm_out.get("citations"), list) else []
            cleaned: list[str] = []
            for citation_id in raw_ids:
                if isinstance(citation_id, str) and citation_id in citation_lookup and citation_id not in cleaned:
                    cleaned.append(citation_id)
            cited_ids = cleaned[:8]
            if not summary_text:
                raise ValueError("Empty summary from LLM")
            if not cited_ids:
                raise ValueError("No valid citations from LLM")
            mode = "online-llm"
            llm_model_used = request.llm_model or os.getenv("OPENAI_MODEL", "gpt-4o-mini")
        except Exception as exc:
            warnings.append(f"LLM summary unavailable, fallback to offline template: {exc}")
            summary_text, cited_ids = _offline_summary(request, bundle)
            mode = "offline-template"
    else:
        if request.llm_mode == "auto":
            warnings.append("LLM summary is disabled or unavailable, using offline template.")
        summary_text, cited_ids = _offline_summary(request, bundle)

    citations = [citation_lookup[citation_id] for citation_id in cited_ids]
    return {
        "mode": mode,
        "disease": request.disease,
        "genes_considered": bundle["genes"][: request.top_n],
        "summary": summary_text,
        "inference_label": "inference",
        "citations": citations,
        "citation_ids": cited_ids,
        "warnings": warnings,
        "llm_mode_requested": llm_mode_requested,
        "llm_model_used": llm_model_used,
    }
