from typing import Any

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field

from api import evidence
from api.auth import current_user

router = APIRouter(prefix="/research", tags=["research"])


class ResearchDirectionRequest(BaseModel):
    disease: str = Field(default="Alzheimer disease")
    genes: list[str] | None = None
    chr: str | None = None
    start: int | None = None
    end: int | None = None
    top_n: int = Field(default=5, ge=1, le=20)


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
