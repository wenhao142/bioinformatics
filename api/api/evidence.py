from collections import defaultdict
from typing import Dict, List

from fastapi import APIRouter, Depends

from api.auth import current_user
from api import variants, omics

router = APIRouter(prefix="/evidence", tags=["evidence"])


def _aggregate():
    # map gene -> data
    table: Dict[str, Dict] = defaultdict(
        lambda: {"gene": None, "variant_count": 0, "mean_abs_logfc": None, "score": 0.0}
    )

    expr_map: Dict[str, List[float]] = defaultdict(list)
    for rec in omics.get_expr_table():
        expr_map[rec["gene"]].append(abs(rec["logfc"]))

    for gene, vals in expr_map.items():
        table[gene]["gene"] = gene
        table[gene]["mean_abs_logfc"] = sum(vals) / len(vals)

    all_vars = variants.STORE.all_variants()
    for v in all_vars:
        gene = variants.nearest_gene(v["chr"], v["pos"])
        if not gene:
            continue
        slot = table[gene]
        slot["gene"] = gene
        slot["variant_count"] += 1

    # compute score: variant_count + 2 * mean_abs_logfc (if present)
    for g, data in table.items():
        score = data["variant_count"]
        if data["mean_abs_logfc"] is not None:
            score += 2 * data["mean_abs_logfc"]
        data["score"] = score

    results = list(table.values())
    results.sort(key=lambda x: x["score"], reverse=True)
    return results


@router.get("/aggregate")
def aggregate(user=Depends(current_user)):
    return {"genes": _aggregate()}
