import itertools
import json
import os
import time
import hashlib
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field
from minio import Minio

from api import canonical_types
from api.auth import current_user

router = APIRouter(prefix="/workflows", tags=["workflows"])

_WORKFLOWS: dict[str, dict[str, Any]] = {}
_WORKFLOW_RUNS: dict[str, dict[str, Any]] = {}


class WorkflowNode(BaseModel):
    node_id: str = Field(pattern=r"^[a-zA-Z0-9][a-zA-Z0-9._-]{0,63}$")
    plugin_id: str = Field(min_length=1, max_length=120)
    version: str | None = Field(default=None, max_length=64)
    input_types: dict[str, str] = Field(default_factory=dict)
    output_types: dict[str, str] = Field(default_factory=dict)
    parameters: dict[str, Any] = Field(default_factory=dict)


class WorkflowEdge(BaseModel):
    from_node: str = Field(min_length=1, max_length=64)
    from_output: str = Field(min_length=1, max_length=120)
    to_node: str = Field(min_length=1, max_length=64)
    to_input: str = Field(min_length=1, max_length=120)


class WorkflowDefinition(BaseModel):
    workflow_id: str = Field(pattern=r"^[a-zA-Z0-9][a-zA-Z0-9._-]{0,119}$")
    nodes: list[WorkflowNode] = Field(min_length=1, max_length=500)
    edges: list[WorkflowEdge] = Field(default_factory=list, max_length=5000)
    parameter_sweeps: dict[str, list[Any]] = Field(default_factory=dict)


class DistributedRunRequest(BaseModel):
    max_workers: int = Field(default=2, ge=1, le=32)
    limit_runs: int | None = Field(default=None, ge=1, le=1000)


def _node_map(defn: WorkflowDefinition) -> dict[str, WorkflowNode]:
    table: dict[str, WorkflowNode] = {}
    for node in defn.nodes:
        if node.node_id in table:
            raise ValueError(f"Duplicate node_id: {node.node_id}")
        table[node.node_id] = node
    return table


def _validate_node_type_registry(defn: WorkflowDefinition) -> None:
    for node in defn.nodes:
        for port, type_id in node.input_types.items():
            if not canonical_types.is_known_canonical_type(type_id):
                raise ValueError(f"Unknown canonical input type for {node.node_id}.{port}: {type_id}")
        for port, type_id in node.output_types.items():
            if not canonical_types.is_known_canonical_type(type_id):
                raise ValueError(f"Unknown canonical output type for {node.node_id}.{port}: {type_id}")


def _build_graph(defn: WorkflowDefinition, nodes: dict[str, WorkflowNode]) -> tuple[dict[str, list[str]], dict[str, int]]:
    graph: dict[str, list[str]] = {node_id: [] for node_id in nodes}
    indegree: dict[str, int] = {node_id: 0 for node_id in nodes}
    for edge in defn.edges:
        if edge.from_node not in nodes:
            raise ValueError(f"Edge references unknown from_node: {edge.from_node}")
        if edge.to_node not in nodes:
            raise ValueError(f"Edge references unknown to_node: {edge.to_node}")
        from_node = nodes[edge.from_node]
        to_node = nodes[edge.to_node]
        if edge.from_output not in from_node.output_types:
            raise ValueError(f"Unknown output port: {edge.from_node}.{edge.from_output}")
        if edge.to_input not in to_node.input_types:
            raise ValueError(f"Unknown input port: {edge.to_node}.{edge.to_input}")
        from_type = from_node.output_types[edge.from_output]
        to_type = to_node.input_types[edge.to_input]
        if from_type != to_type:
            raise ValueError(
                f"I/O type mismatch: {edge.from_node}.{edge.from_output}({from_type}) -> "
                f"{edge.to_node}.{edge.to_input}({to_type})"
            )
        graph[edge.from_node].append(edge.to_node)
        indegree[edge.to_node] += 1
    return graph, indegree


def _topological_order(graph: dict[str, list[str]], indegree: dict[str, int]) -> list[str]:
    q = deque(sorted([node_id for node_id, deg in indegree.items() if deg == 0]))
    order: list[str] = []
    local_indegree = dict(indegree)
    while q:
        node = q.popleft()
        order.append(node)
        for nxt in graph.get(node, []):
            local_indegree[nxt] -= 1
            if local_indegree[nxt] == 0:
                q.append(nxt)
    if len(order) != len(indegree):
        raise ValueError("Workflow graph contains a cycle")
    return order


def _branching_nodes(edges: list[WorkflowEdge]) -> list[str]:
    out_degree: dict[str, int] = defaultdict(int)
    for edge in edges:
        out_degree[edge.from_node] += 1
    return sorted([node_id for node_id, deg in out_degree.items() if deg > 1])


def _expand_parameter_sweeps(defn: WorkflowDefinition) -> list[dict[str, Any]]:
    if not defn.parameter_sweeps:
        return [{"index": 0, "overrides": {}}]

    node_ids = {node.node_id for node in defn.nodes}
    keys = sorted(defn.parameter_sweeps.keys())
    values_lists: list[list[Any]] = []
    for key in keys:
        values = defn.parameter_sweeps[key]
        if not isinstance(values, list) or not values:
            raise ValueError(f"parameter_sweeps.{key} must be a non-empty array")
        if "." not in key:
            raise ValueError(f"parameter_sweeps key must be node.param: {key}")
        node_id, param_name = key.split(".", 1)
        if node_id not in node_ids:
            raise ValueError(f"parameter_sweeps references unknown node: {node_id}")
        if not param_name:
            raise ValueError(f"parameter_sweeps key missing param name: {key}")
        values_lists.append(values)

    expanded: list[dict[str, Any]] = []
    for idx, combo in enumerate(itertools.product(*values_lists)):
        overrides: dict[str, dict[str, Any]] = {}
        for key, value in zip(keys, combo):
            node_id, param_name = key.split(".", 1)
            overrides.setdefault(node_id, {})[param_name] = value
        expanded.append({"index": idx, "overrides": overrides})
    return expanded


def _validation_report(defn: WorkflowDefinition) -> dict[str, Any]:
    nodes = _node_map(defn)
    _validate_node_type_registry(defn)
    graph, indegree = _build_graph(defn, nodes)
    order = _topological_order(graph, indegree)
    sweeps = _expand_parameter_sweeps(defn)
    return {
        "workflow_id": defn.workflow_id,
        "node_count": len(defn.nodes),
        "edge_count": len(defn.edges),
        "acyclic": True,
        "topological_order": order,
        "branching_nodes": _branching_nodes(defn.edges),
        "sweep_count": len(sweeps),
    }


def _sha256_json(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _get_workflow_or_404(workflow_id: str) -> dict[str, Any]:
    item = _WORKFLOWS.get(workflow_id)
    if item is None:
        raise HTTPException(status_code=404, detail="Workflow not found")
    return item


def _minio_client() -> Minio | None:
    endpoint = os.getenv("MINIO_ENDPOINT", "").strip()
    if not endpoint:
        return None
    host = endpoint.replace("http://", "").replace("https://", "").strip("/")
    secure = endpoint.startswith("https://")
    access_key = os.getenv("MINIO_ACCESS_KEY", os.getenv("MINIO_ROOT_USER", "minioadmin"))
    secret_key = os.getenv("MINIO_SECRET_KEY", os.getenv("MINIO_ROOT_PASSWORD", "minioadmin123"))
    return Minio(host, access_key=access_key, secret_key=secret_key, secure=secure)


def _cloud_export_path(workflow_id: str) -> Path:
    return Path("results") / "workflows" / "cloud" / f"{workflow_id}.json"


@router.post("/import", status_code=201)
def import_workflow(defn: WorkflowDefinition, user=Depends(current_user)):
    if defn.workflow_id in _WORKFLOWS:
        raise HTTPException(status_code=409, detail="Workflow already exists")
    try:
        report = _validation_report(defn)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    _WORKFLOWS[defn.workflow_id] = {
        "workflow": defn.model_dump(),
        "created_by": user["email"],
        "report": report,
    }
    return {"workflow_id": defn.workflow_id, "report": report}


@router.get("")
def list_workflows(user=Depends(current_user)):
    rows = []
    for workflow_id, item in sorted(_WORKFLOWS.items()):
        rows.append(
            {
                "workflow_id": workflow_id,
                "created_by": item["created_by"],
                "node_count": item["report"]["node_count"],
                "edge_count": item["report"]["edge_count"],
            }
        )
    return {"workflows": rows}


@router.get("/{workflow_id}/export")
def export_workflow(workflow_id: str, user=Depends(current_user)):
    item = _get_workflow_or_404(workflow_id)
    return {"workflow": item["workflow"]}


@router.post("/{workflow_id}/validate")
def validate_workflow(workflow_id: str, user=Depends(current_user)):
    item = _get_workflow_or_404(workflow_id)
    defn = WorkflowDefinition.model_validate(item["workflow"])
    try:
        report = _validation_report(defn)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    item["report"] = report
    return {"report": report}


@router.get("/{workflow_id}/plan")
def workflow_plan(workflow_id: str, user=Depends(current_user)):
    item = _get_workflow_or_404(workflow_id)
    return {
        "workflow_id": workflow_id,
        "topological_order": item["report"]["topological_order"],
        "branching_nodes": item["report"]["branching_nodes"],
    }


@router.post("/{workflow_id}/sweeps/expand")
def expand_sweeps(workflow_id: str, user=Depends(current_user)):
    item = _get_workflow_or_404(workflow_id)
    defn = WorkflowDefinition.model_validate(item["workflow"])
    try:
        expanded = _expand_parameter_sweeps(defn)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    return {"workflow_id": workflow_id, "expanded_runs": expanded, "count": len(expanded)}


@router.post("/{workflow_id}/execute/distributed")
def execute_distributed(workflow_id: str, body: DistributedRunRequest, user=Depends(current_user)):
    item = _get_workflow_or_404(workflow_id)
    defn = WorkflowDefinition.model_validate(item["workflow"])
    expanded = _expand_parameter_sweeps(defn)
    if body.limit_runs is not None:
        expanded = expanded[: body.limit_runs]

    plan = item["report"]["topological_order"]
    start = time.time()
    run_id = f"{workflow_id}-dist-{int(start * 1000)}"

    def _execute_one(expanded_run: dict[str, Any]) -> dict[str, Any]:
        payload = {
            "workflow_id": workflow_id,
            "index": expanded_run["index"],
            "overrides": expanded_run["overrides"],
            "plan": plan,
        }
        return {
            "index": expanded_run["index"],
            "status": "completed",
            "plan": plan,
            "overrides": expanded_run["overrides"],
            "signature": _sha256_json(payload),
        }

    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=body.max_workers) as executor:
        futures = [executor.submit(_execute_one, entry) for entry in expanded]
        for fut in as_completed(futures):
            results.append(fut.result())
    results.sort(key=lambda row: row["index"])

    summary = {
        "run_id": run_id,
        "workflow_id": workflow_id,
        "submitted_runs": len(expanded),
        "completed_runs": len(results),
        "max_workers": body.max_workers,
        "duration_ms": int((time.time() - start) * 1000),
    }
    _WORKFLOW_RUNS[run_id] = {"summary": summary, "results": results}
    return {"summary": summary, "results": results}


@router.get("/runs/{run_id}")
def get_distributed_run(run_id: str, user=Depends(current_user)):
    item = _WORKFLOW_RUNS.get(run_id)
    if item is None:
        raise HTTPException(status_code=404, detail="Workflow run not found")
    return item


@router.post("/{workflow_id}/cloud/export")
def cloud_export_workflow(workflow_id: str, user=Depends(current_user)):
    item = _get_workflow_or_404(workflow_id)
    payload = json.dumps(item["workflow"], ensure_ascii=True).encode("utf-8")
    object_name = f"{workflow_id}.json"

    client = _minio_client()
    if client is not None:
        try:
            bucket = os.getenv("WORKFLOW_BUCKET", "workflows")
            found = client.bucket_exists(bucket)
            if not found:
                client.make_bucket(bucket)
            from io import BytesIO

            client.put_object(bucket, object_name, BytesIO(payload), len(payload), content_type="application/json")
            return {"workflow_id": workflow_id, "storage": "minio", "uri": f"s3://{bucket}/{object_name}"}
        except Exception:
            pass

    target = _cloud_export_path(workflow_id)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(payload)
    return {"workflow_id": workflow_id, "storage": "local", "uri": str(target)}


@router.post("/cloud/import/{workflow_id}")
def cloud_import_workflow(workflow_id: str, user=Depends(current_user)):
    if workflow_id in _WORKFLOWS:
        raise HTTPException(status_code=409, detail="Workflow already exists")

    object_name = f"{workflow_id}.json"
    client = _minio_client()
    payload: dict[str, Any] | None = None
    if client is not None:
        try:
            bucket = os.getenv("WORKFLOW_BUCKET", "workflows")
            response = client.get_object(bucket, object_name)
            raw = response.read()
            payload = json.loads(raw.decode("utf-8"))
        except Exception:
            payload = None

    if payload is None:
        src = _cloud_export_path(workflow_id)
        if not src.exists():
            raise HTTPException(status_code=404, detail="Cloud workflow artifact not found")
        payload = json.loads(src.read_text(encoding="utf-8"))

    defn = WorkflowDefinition.model_validate(payload)
    report = _validation_report(defn)
    _WORKFLOWS[defn.workflow_id] = {
        "workflow": defn.model_dump(),
        "created_by": user["email"],
        "report": report,
    }
    return {"workflow_id": defn.workflow_id, "report": report}


def reset_workflow_store_for_tests() -> None:
    _WORKFLOWS.clear()
    _WORKFLOW_RUNS.clear()
