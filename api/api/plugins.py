import hashlib
import json
import os
import subprocess
import sys
import time
from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query, status
from pydantic import BaseModel, Field

from api import canonical_types
from api import evidence
from api import manifest_schema
from api import method_registry
from api.auth import current_user

router = APIRouter(prefix="/plugins", tags=["plugins"])

_PLUGIN_RUNS: dict[str, list[dict[str, Any]]] = {}
_PLUGIN_RESULTS: dict[str, dict[str, Any]] = {}


class PluginResources(BaseModel):
    cpu_millicores: int = Field(default=1000, ge=100, le=64000)
    memory_mb: int = Field(default=1024, ge=128, le=262144)
    gpu: bool = False


class AdapterSpec(BaseModel):
    from_type: str = Field(min_length=1, max_length=120)
    to_type: str = Field(min_length=1, max_length=120)
    strategy: str = Field(default="transform", min_length=1, max_length=80)


class PluginManifest(BaseModel):
    plugin_id: str = Field(pattern=r"^[a-z0-9][a-z0-9._-]{1,63}$")
    name: str = Field(min_length=1, max_length=120)
    version: str = Field(pattern=r"^\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?$")
    image: str = Field(min_length=3, max_length=512)
    description: str | None = Field(default=None, max_length=2000)
    input_schema: dict[str, Any]
    output_schema: dict[str, Any]
    resources: PluginResources = PluginResources()
    enabled: bool = True
    tags: list[str] = Field(default_factory=list, max_length=30)
    adapter: AdapterSpec | None = None
    deprecated: bool = False
    deprecation_message: str | None = Field(default=None, max_length=400)
    superseded_by: str | None = Field(default=None, max_length=64)


class PluginRunRequest(BaseModel):
    project_id: str = Field(min_length=1, max_length=120)
    chr: str = Field(min_length=1, max_length=32)
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    top_n: int = Field(default=20, ge=1, le=200)
    parameters: dict[str, Any] | None = None


class PluginEnableRequest(BaseModel):
    enabled: bool


class PluginCompareRequest(BaseModel):
    left_run_id: str = Field(min_length=1, max_length=120)
    right_run_id: str = Field(min_length=1, max_length=120)


def _sha256_json(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _validate_json_schema(name: str, schema: dict[str, Any]) -> None:
    if not isinstance(schema, dict):
        raise HTTPException(status_code=400, detail=f"{name} must be a JSON object")
    if schema.get("type") != "object":
        raise HTTPException(status_code=400, detail=f"{name}.type must be 'object'")

    properties = schema.get("properties")
    if not isinstance(properties, dict) or not properties:
        raise HTTPException(status_code=400, detail=f"{name}.properties must be a non-empty object")

    allowed_types = {"string", "number", "integer", "boolean", "object", "array"}
    for field_name, field_spec in properties.items():
        if not isinstance(field_spec, dict):
            raise HTTPException(status_code=400, detail=f"{name}.properties.{field_name} must be an object")
        field_type = field_spec.get("type")
        if field_type not in allowed_types:
            raise HTTPException(
                status_code=400,
                detail=f"{name}.properties.{field_name}.type must be one of {sorted(allowed_types)}",
            )

    required = schema.get("required", [])
    if required is not None:
        if not isinstance(required, list) or any(not isinstance(item, str) for item in required):
            raise HTTPException(status_code=400, detail=f"{name}.required must be a list of strings")
        missing = [item for item in required if item not in properties]
        if missing:
            raise HTTPException(status_code=400, detail=f"{name}.required contains unknown fields: {missing}")

    additional = schema.get("additionalProperties", True)
    if not isinstance(additional, (bool, dict)):
        raise HTTPException(status_code=400, detail=f"{name}.additionalProperties must be bool or object")


def _admin_only(user: dict[str, Any]) -> None:
    if user["role"] != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin only")


def _validate_output_against_schema(schema: dict[str, Any], payload: dict[str, Any]) -> None:
    required = schema.get("required", [])
    for key in required:
        if key not in payload:
            raise HTTPException(status_code=500, detail=f"Plugin output missing required field: {key}")


def _run_builtin_baseline(request: PluginRunRequest) -> dict[str, Any]:
    return evidence.compute_ranked_evidence(
        chr=request.chr,
        start=request.start,
        end=request.end,
        top_n=request.top_n,
    )


def _run_workflow_wrapper(kind: str, payload: dict[str, Any]) -> dict[str, Any]:
    params = payload.get("parameters") or {}
    workflow_file = str(params.get("workflow_file", "")).strip() or "<unset>"
    dry_run = bool(params.get("dry_run", True))
    command_preview = {"engine": kind, "workflow_file": workflow_file, "dry_run": dry_run}
    return {
        "status": "planned",
        "workflow_engine": kind,
        "command_preview": command_preview,
        "ranked_genes": [],
        "ranked_loci": [],
    }


def _docker_image_digest(image: str) -> str | None:
    try:
        completed = subprocess.run(
            ["docker", "image", "inspect", image, "--format", "{{index .RepoDigests 0}}"],
            text=True,
            capture_output=True,
            timeout=20,
            check=False,
        )
    except Exception:
        return None
    if completed.returncode != 0:
        return None
    digest = (completed.stdout or "").strip()
    return digest or None


def _runner_tool_versions() -> dict[str, str]:
    return {
        "plugin_runner": "1.1.0",
        "python": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
    }


def _run_docker_plugin(image: str, payload: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    allow_docker = os.getenv("PLUGIN_RUNNER_ALLOW_DOCKER", "false").lower() in {"1", "true", "yes", "on"}
    if not allow_docker:
        raise HTTPException(status_code=400, detail="Docker plugin execution is disabled")

    started = time.time()
    try:
        completed = subprocess.run(
            ["docker", "run", "--rm", "--network", "none", "--read-only", "--cap-drop", "ALL", "-i", image],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            timeout=120,
            check=False,
        )
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Docker run failed: {exc}")

    if completed.returncode != 0:
        stderr = (completed.stderr or "").strip()
        raise HTTPException(status_code=500, detail=f"Docker plugin failed: {stderr or 'unknown error'}")
    try:
        parsed = json.loads(completed.stdout or "{}")
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=500, detail=f"Plugin output is not valid JSON: {exc}")
    return parsed, {
        "duration_ms": int((time.time() - started) * 1000),
        "stdout_tail": (completed.stdout or "")[-2000:],
        "stderr_tail": (completed.stderr or "")[-2000:],
    }


def _validate_image_reference(image: str) -> None:
    if image.startswith("builtin/"):
        return
    if image.endswith(":latest"):
        raise HTTPException(status_code=400, detail="image tag 'latest' is not allowed; use pinned tag or digest")
    last = image.rsplit("/", 1)[-1]
    has_digest = "@sha256:" in image
    has_tag = ":" in last
    if not has_digest and not has_tag:
        raise HTTPException(status_code=400, detail="image must include explicit tag or digest")
    if has_tag and last.endswith(":"):
        raise HTTPException(status_code=400, detail="image tag cannot be empty")


def _validate_manifest_policy(manifest: PluginManifest) -> None:
    _validate_image_reference(manifest.image)
    if manifest.deprecated and not manifest.deprecation_message:
        raise HTTPException(status_code=400, detail="deprecated plugins require deprecation_message")
    if "adapter" in manifest.tags and manifest.adapter is None:
        raise HTTPException(status_code=400, detail="adapter plugin requires adapter specification")
    if manifest.adapter is not None:
        if not canonical_types.is_known_canonical_type(manifest.adapter.from_type):
            raise HTTPException(status_code=400, detail=f"adapter.from_type unknown: {manifest.adapter.from_type}")
        if not canonical_types.is_known_canonical_type(manifest.adapter.to_type):
            raise HTTPException(status_code=400, detail=f"adapter.to_type unknown: {manifest.adapter.to_type}")


def _plugin_runs(plugin_id: str) -> list[dict[str, Any]]:
    return _PLUGIN_RUNS.setdefault(plugin_id, [])


def _plugin_or_404(plugin_id: str, version: str | None = None) -> dict[str, Any]:
    plugin = method_registry.get_method(plugin_id, version=version)
    if plugin is None:
        raise HTTPException(status_code=404, detail="Plugin not found")
    return plugin


def _ensure_run_access(run: dict[str, Any], user: dict[str, Any]) -> None:
    if user["role"] == "admin":
        return
    if run["created_by"] == user["email"]:
        return
    raise HTTPException(status_code=403, detail="Forbidden")


def _extract_gene_scores(result: dict[str, Any]) -> dict[str, float]:
    table: dict[str, float] = {}
    for row in result.get("gene_scores", []):
        gene = row.get("gene")
        score = row.get("score")
        if isinstance(gene, str) and isinstance(score, (int, float)):
            table[gene] = float(score)
    for row in result.get("ranked_genes", []):
        gene = row.get("gene")
        score = row.get("score")
        if isinstance(gene, str) and isinstance(score, (int, float)):
            table[gene] = float(score)
    return table


@router.post("/register", status_code=201)
def register_plugin(manifest: dict[str, Any], user=Depends(current_user)):
    _admin_only(user)
    try:
        manifest_schema.validate_plugin_manifest_or_raise(manifest)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    parsed = PluginManifest.model_validate(manifest)
    _validate_manifest_policy(parsed)
    if method_registry.get_method(parsed.plugin_id, version=parsed.version) is not None:
        raise HTTPException(status_code=409, detail="Plugin already registered")

    _validate_json_schema("input_schema", parsed.input_schema)
    _validate_json_schema("output_schema", parsed.output_schema)
    try:
        canonical_types.validate_schema_canonical_types("input_schema", parsed.input_schema)
        canonical_types.validate_schema_canonical_types("output_schema", parsed.output_schema)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    record = parsed.model_dump()
    record["registered_by"] = user["email"]
    created = method_registry.register_method(record)
    if created is None:
        raise HTTPException(status_code=409, detail="Plugin already registered")
    return {"plugin": created}


@router.post("/{plugin_id}/run")
def run_plugin(plugin_id: str, body: PluginRunRequest, user=Depends(current_user)):
    plugin = _plugin_or_404(plugin_id)
    if not plugin.get("enabled", True):
        raise HTTPException(status_code=400, detail="Plugin is disabled")
    if plugin.get("deprecated"):
        raise HTTPException(status_code=400, detail="Plugin is deprecated and cannot be executed")

    payload = {
        "project_id": body.project_id,
        "chr": body.chr,
        "start": body.start,
        "end": body.end,
        "top_n": body.top_n,
        "parameters": body.parameters or {},
    }
    input_hash = _sha256_json(payload)
    params_hash = _sha256_json(payload["parameters"])
    started_at = time.time()

    if "baseline" in (plugin.get("tags") or []) or plugin.get("image", "").startswith("builtin/"):
        tags = plugin.get("tags") or []
        if "snakemake-wrapper" in tags:
            engine = "snakemake-wrapper"
            result = _run_workflow_wrapper("snakemake", payload)
        elif "nextflow-wrapper" in tags:
            engine = "nextflow-wrapper"
            result = _run_workflow_wrapper("nextflow", payload)
        else:
            engine = "builtin"
            result = _run_builtin_baseline(body)
        execution_meta = {"duration_ms": int((time.time() - started_at) * 1000), "stdout_tail": "", "stderr_tail": ""}
    else:
        engine = "docker"
        result, execution_meta = _run_docker_plugin(plugin["image"], payload)

    if not isinstance(result, dict):
        raise HTTPException(status_code=500, detail="Plugin output must be a JSON object")
    _validate_output_against_schema(plugin["output_schema"], result)

    runs = _plugin_runs(plugin_id)
    run_id = f"{plugin_id}-run-{len(runs) + 1}"
    result_hash = _sha256_json(result)
    run_record = {
        "run_id": run_id,
        "plugin_id": plugin_id,
        "plugin_version": plugin["version"],
        "project_id": body.project_id,
        "created_at": time.time(),
        "created_by": user["email"],
        "engine": engine,
        "image": plugin["image"],
        "container_image_digest": _docker_image_digest(plugin["image"]) if engine == "docker" else None,
        "tool_versions": _runner_tool_versions(),
        "params_hash": params_hash,
        "parameters": payload["parameters"],
        "input_hash": input_hash,
        "result_hash": result_hash,
        "status": "completed",
        "execution": {
            "started_at": started_at,
            "finished_at": time.time(),
            "duration_ms": execution_meta["duration_ms"],
        },
        "logs": {
            "stdout_tail": execution_meta["stdout_tail"],
            "stderr_tail": execution_meta["stderr_tail"],
        },
        "counts": {
            "ranked_genes": len(result.get("ranked_genes", [])),
            "ranked_loci": len(result.get("ranked_loci", [])),
        },
    }
    runs.append(run_record)
    _PLUGIN_RESULTS[run_id] = result

    return {"run": run_record, "result": result}


@router.get("/{plugin_id}/runs")
def list_plugin_runs(plugin_id: str, user=Depends(current_user)):
    _plugin_or_404(plugin_id)
    runs = _plugin_runs(plugin_id)
    visible = runs if user["role"] == "admin" else [run for run in runs if run["created_by"] == user["email"]]
    return {"runs": list(reversed(visible))}


@router.get("/{plugin_id}/runs/{run_id}")
def get_plugin_run(plugin_id: str, run_id: str, user=Depends(current_user)):
    _plugin_or_404(plugin_id)
    for run in _plugin_runs(plugin_id):
        if run["run_id"] != run_id:
            continue
        _ensure_run_access(run, user)
        return {"run": run}
    raise HTTPException(status_code=404, detail="Run not found")


@router.get("/{plugin_id}/runs/{run_id}/result")
def get_plugin_run_result(plugin_id: str, run_id: str, user=Depends(current_user)):
    _plugin_or_404(plugin_id)
    run = next((r for r in _plugin_runs(plugin_id) if r["run_id"] == run_id), None)
    if run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    _ensure_run_access(run, user)
    result = _PLUGIN_RESULTS.get(run_id)
    if result is None:
        raise HTTPException(status_code=404, detail="Run result not found")
    return {"run_id": run_id, "result": result}


@router.patch("/{plugin_id}/enabled")
def set_plugin_enabled(
    plugin_id: str,
    body: PluginEnableRequest,
    version: str | None = Query(default=None, description="Optional plugin version; default is latest"),
    user=Depends(current_user),
):
    _admin_only(user)
    updated = method_registry.set_enabled(plugin_id, body.enabled, version=version)
    if updated is None:
        raise HTTPException(status_code=404, detail="Plugin not found")
    return {
        "plugin_id": plugin_id,
        "version": updated["version"],
        "enabled": body.enabled,
    }


@router.post("/{plugin_id}/compare")
def compare_plugin_runs(plugin_id: str, body: PluginCompareRequest, user=Depends(current_user)):
    _plugin_or_404(plugin_id)
    runs = _plugin_runs(plugin_id)
    left_run = next((r for r in runs if r["run_id"] == body.left_run_id), None)
    right_run = next((r for r in runs if r["run_id"] == body.right_run_id), None)
    if left_run is None or right_run is None:
        raise HTTPException(status_code=404, detail="Run not found")
    _ensure_run_access(left_run, user)
    _ensure_run_access(right_run, user)

    left_result = _PLUGIN_RESULTS.get(left_run["run_id"])
    right_result = _PLUGIN_RESULTS.get(right_run["run_id"])
    if left_result is None or right_result is None:
        raise HTTPException(status_code=404, detail="Run result not found")

    left_scores = _extract_gene_scores(left_result)
    right_scores = _extract_gene_scores(right_result)
    left_genes = set(left_scores.keys())
    right_genes = set(right_scores.keys())
    overlap = sorted(left_genes & right_genes)

    score_deltas = [
        {
            "gene": gene,
            "left_score": round(left_scores[gene], 6),
            "right_score": round(right_scores[gene], 6),
            "delta": round(left_scores[gene] - right_scores[gene], 6),
        }
        for gene in overlap
    ]
    score_deltas.sort(key=lambda row: abs(row["delta"]), reverse=True)

    return {
        "plugin_id": plugin_id,
        "left_run_id": left_run["run_id"],
        "right_run_id": right_run["run_id"],
        "overlap_genes": overlap,
        "only_left_genes": sorted(left_genes - right_genes),
        "only_right_genes": sorted(right_genes - left_genes),
        "score_deltas": score_deltas,
    }


@router.get("")
def list_plugins(
    all_versions: bool = Query(default=False, description="If true, return all versions; else latest per plugin_id"),
    user=Depends(current_user),
):
    return {"plugins": method_registry.list_methods(latest_only=not all_versions)}


@router.get("/canonical-types")
def list_canonical_types(user=Depends(current_user)):
    return {
        "registry_version": canonical_types.CANONICAL_TYPE_REGISTRY_VERSION,
        "types": canonical_types.list_canonical_types(),
    }


@router.get("/search")
def search_plugins(
    q: str | None = Query(default=None, description="Search by plugin_id or name substring"),
    tag: str | None = Query(default=None, description="Filter by tag"),
    enabled: bool | None = Query(default=None, description="Filter by enabled state"),
    all_versions: bool = Query(default=False, description="Search all versions when true"),
    user=Depends(current_user),
):
    rows = method_registry.list_methods(latest_only=not all_versions)
    filtered: list[dict[str, Any]] = []
    q_norm = (q or "").strip().lower()
    tag_norm = (tag or "").strip().lower()
    for row in rows:
        if q_norm and q_norm not in row.get("plugin_id", "").lower() and q_norm not in row.get("name", "").lower():
            continue
        if tag_norm and tag_norm not in [t.lower() for t in (row.get("tags") or [])]:
            continue
        if enabled is not None and bool(row.get("enabled", True)) != enabled:
            continue
        filtered.append(row)
    return {"plugins": filtered}


@router.get("/adapters")
def list_adapters(user=Depends(current_user)):
    rows = method_registry.list_methods(latest_only=False)
    adapters = []
    for row in rows:
        if "adapter" in (row.get("tags") or []) and row.get("adapter"):
            adapters.append(row)
    return {"adapters": adapters}


@router.get("/{plugin_id}")
def get_plugin(
    plugin_id: str,
    version: str | None = Query(default=None, description="Optional plugin version; default is latest"),
    user=Depends(current_user),
):
    plugin = _plugin_or_404(plugin_id, version=version)
    return {"plugin": plugin}


def reset_plugin_state_for_tests() -> None:
    method_registry.clear_registry()
    _PLUGIN_RUNS.clear()
    _PLUGIN_RESULTS.clear()
