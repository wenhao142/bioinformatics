import itertools
import json
import os
import time
import hashlib
import shlex
import shutil
import subprocess
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel, Field
from minio import Minio
import psycopg2

from api import canonical_types
from api.auth import current_user

router = APIRouter(prefix="/workflows", tags=["workflows"])

_WORKFLOWS: dict[str, dict[str, Any]] = {}
_WORKFLOW_RUNS: dict[str, dict[str, Any]] = {}


def _pg_dsn() -> str:
    dsn = os.getenv("DATABASE_URL")
    if dsn:
        return dsn
    return (
        f"dbname={os.getenv('POSTGRES_DB', 'omics')} "
        f"user={os.getenv('POSTGRES_USER', 'postgres')} "
        f"password={os.getenv('POSTGRES_PASSWORD', 'postgres')} "
        f"host={os.getenv('POSTGRES_HOST', 'localhost')} "
        f"port={os.getenv('POSTGRES_PORT', '5432')}"
    )


def _workflow_definition_root() -> Path:
    return Path("results") / "workflows" / "definitions"


def _workflow_definition_path(workflow_id: str) -> Path:
    return _workflow_definition_root() / f"{workflow_id}.json"


def _workflow_run_record_path(run_id: str) -> Path:
    return _workflow_run_root() / run_id / "run-record.json"


def _workflow_project_id(workflow_id: str) -> str | None:
    if "__" not in workflow_id:
        return None
    project_id, _ = workflow_id.split("__", 1)
    return project_id or None


def _write_workflow_definition_file(
    workflow_id: str,
    workflow: dict[str, Any],
    report: dict[str, Any],
    created_by: str,
) -> str:
    path = _workflow_definition_path(workflow_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "workflow_id": workflow_id,
        "project_id": _workflow_project_id(workflow_id),
        "created_by": created_by,
        "workflow": workflow,
        "report": report,
        "saved_at": int(time.time()),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return str(path.resolve())


def _write_run_record_file(
    run_id: str,
    workflow_id: str,
    created_by: str,
    summary: dict[str, Any],
    results: list[dict[str, Any]],
) -> str:
    path = _workflow_run_record_path(run_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "run_id": run_id,
        "workflow_id": workflow_id,
        "project_id": _workflow_project_id(workflow_id),
        "created_by": created_by,
        "summary": summary,
        "results": results,
        "saved_at": int(time.time()),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return str(path.resolve())


def _host_visible_path(path: str | None) -> str | None:
    if not path:
        return path
    try:
        resolved = Path(path)
        container_root = Path("/app/results")
        if resolved.is_absolute() and str(resolved).startswith(str(container_root)):
            return str(Path("results") / resolved.relative_to(container_root))
    except Exception:
        pass
    return path


def _present_workflow_row(row: dict[str, Any]) -> dict[str, Any]:
    payload = dict(row)
    payload["local_export_path"] = _host_visible_path(payload.get("local_export_path"))
    return payload


def _present_run_summary(summary: dict[str, Any]) -> dict[str, Any]:
    payload = dict(summary)
    payload["local_record_path"] = _host_visible_path(payload.get("local_record_path"))
    return payload


def _present_run_result(result: dict[str, Any]) -> dict[str, Any]:
    payload = dict(result)
    workflow_artifacts = dict(payload.get("workflow_artifacts") or {})
    for key in ("run_dir", "snakefile", "configfile"):
        workflow_artifacts[key] = _host_visible_path(workflow_artifacts.get(key))
    payload["workflow_artifacts"] = workflow_artifacts
    payload["outputs"] = [_host_visible_path(item) or item for item in payload.get("outputs", [])]
    return payload


def _present_run_record(item: dict[str, Any]) -> dict[str, Any]:
    payload = dict(item)
    local_record_path = _host_visible_path(payload.get("local_record_path"))
    payload["summary"] = _present_run_summary({**payload["summary"], "local_record_path": local_record_path})
    payload["results"] = [_present_run_result(row) for row in payload["results"]]
    payload["local_record_path"] = local_record_path
    return payload


class MemoryWorkflowStore:
    def create_workflow(
        self,
        workflow_id: str,
        workflow: dict[str, Any],
        report: dict[str, Any],
        created_by: str,
        local_export_path: str,
    ) -> bool:
        if workflow_id in _WORKFLOWS:
            return False
        _WORKFLOWS[workflow_id] = {
            "workflow": workflow,
            "created_by": created_by,
            "report": report,
            "local_export_path": local_export_path,
        }
        return True

    def get_workflow(self, workflow_id: str) -> dict[str, Any] | None:
        item = _WORKFLOWS.get(workflow_id)
        return dict(item) if item else None

    def list_workflows(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for workflow_id, item in sorted(_WORKFLOWS.items()):
            rows.append(
                {
                    "workflow_id": workflow_id,
                    "created_by": item["created_by"],
                    "report": item["report"],
                    "local_export_path": item.get("local_export_path"),
                }
            )
        return rows

    def update_workflow_report(self, workflow_id: str, report: dict[str, Any]) -> None:
        if workflow_id not in _WORKFLOWS:
            raise KeyError(workflow_id)
        _WORKFLOWS[workflow_id]["report"] = report

    def save_run(self, run_id: str, workflow_id: str, created_by: str, summary: dict[str, Any], results: list[dict[str, Any]], local_record_path: str) -> None:
        _WORKFLOW_RUNS[run_id] = {
            "summary": summary,
            "results": results,
            "created_by": created_by,
            "workflow_id": workflow_id,
            "local_record_path": local_record_path,
        }

    def get_run(self, run_id: str) -> dict[str, Any] | None:
        item = _WORKFLOW_RUNS.get(run_id)
        return dict(item) if item else None

    def list_runs(self, limit: int = 50) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for run_id, item in sorted(_WORKFLOW_RUNS.items(), key=lambda pair: pair[1]["summary"]["run_id"], reverse=True):
            rows.append(
                {
                    **item["summary"],
                    "local_record_path": item.get("local_record_path"),
                }
            )
            if len(rows) >= limit:
                break
        return rows

    def reset(self) -> None:
        _WORKFLOWS.clear()
        _WORKFLOW_RUNS.clear()


class PgWorkflowStore:
    def __init__(self):
        self.conn = psycopg2.connect(_pg_dsn())
        self.conn.autocommit = True
        self._ensure_tables()

    def _ensure_tables(self) -> None:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                CREATE TABLE IF NOT EXISTS workflow_definitions (
                    workflow_id TEXT PRIMARY KEY,
                    project_id TEXT,
                    created_by TEXT NOT NULL,
                    workflow_json JSONB NOT NULL,
                    report_json JSONB NOT NULL,
                    local_export_path TEXT NOT NULL,
                    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
                );
                CREATE TABLE IF NOT EXISTS workflow_runs (
                    run_id TEXT PRIMARY KEY,
                    workflow_id TEXT NOT NULL,
                    project_id TEXT,
                    created_by TEXT NOT NULL,
                    summary_json JSONB NOT NULL,
                    results_json JSONB NOT NULL,
                    local_record_path TEXT NOT NULL,
                    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
                );
                ALTER TABLE workflow_definitions ADD COLUMN IF NOT EXISTS project_id TEXT;
                ALTER TABLE workflow_runs ADD COLUMN IF NOT EXISTS project_id TEXT;
                CREATE INDEX IF NOT EXISTS workflow_definitions_project_id_idx ON workflow_definitions(project_id);
                CREATE INDEX IF NOT EXISTS workflow_runs_workflow_id_idx ON workflow_runs(workflow_id);
                CREATE INDEX IF NOT EXISTS workflow_runs_project_id_idx ON workflow_runs(project_id);
                CREATE INDEX IF NOT EXISTS workflow_runs_created_by_idx ON workflow_runs(created_by);
                """
            )

    def create_workflow(
        self,
        workflow_id: str,
        workflow: dict[str, Any],
        report: dict[str, Any],
        created_by: str,
        local_export_path: str,
    ) -> bool:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                INSERT INTO workflow_definitions(workflow_id, project_id, created_by, workflow_json, report_json, local_export_path)
                VALUES (%s, %s, %s, %s::jsonb, %s::jsonb, %s)
                ON CONFLICT (workflow_id) DO NOTHING
                RETURNING workflow_id
                """,
                (
                    workflow_id,
                    _workflow_project_id(workflow_id),
                    created_by,
                    json.dumps(workflow, sort_keys=True),
                    json.dumps(report, sort_keys=True),
                    local_export_path,
                ),
            )
            created = cur.fetchone()
        return bool(created)

    def get_workflow(self, workflow_id: str) -> dict[str, Any] | None:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                SELECT project_id, created_by, workflow_json, report_json, local_export_path
                FROM workflow_definitions
                WHERE workflow_id=%s
                """,
                (workflow_id,),
            )
            row = cur.fetchone()
        if not row:
            return None
        return {
            "project_id": row[0],
            "workflow": row[2],
            "created_by": row[1],
            "report": row[3],
            "local_export_path": row[4],
        }

    def list_workflows(self) -> list[dict[str, Any]]:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                SELECT workflow_id, project_id, created_by, report_json, local_export_path
                FROM workflow_definitions
                ORDER BY workflow_id
                """
            )
            rows = cur.fetchall()
        return [
            {
                "workflow_id": row[0],
                "project_id": row[1],
                "created_by": row[2],
                "report": row[3],
                "local_export_path": row[4],
            }
            for row in rows
        ]

    def update_workflow_report(self, workflow_id: str, report: dict[str, Any]) -> None:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                UPDATE workflow_definitions
                SET report_json=%s::jsonb, updated_at=NOW()
                WHERE workflow_id=%s
                """,
                (json.dumps(report, sort_keys=True), workflow_id),
            )

    def save_run(self, run_id: str, workflow_id: str, created_by: str, summary: dict[str, Any], results: list[dict[str, Any]], local_record_path: str) -> None:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                INSERT INTO workflow_runs(run_id, workflow_id, project_id, created_by, summary_json, results_json, local_record_path)
                VALUES (%s, %s, %s, %s, %s::jsonb, %s::jsonb, %s)
                ON CONFLICT (run_id) DO UPDATE
                SET project_id=EXCLUDED.project_id,
                    summary_json=EXCLUDED.summary_json,
                    results_json=EXCLUDED.results_json,
                    local_record_path=EXCLUDED.local_record_path
                """,
                (
                    run_id,
                    workflow_id,
                    _workflow_project_id(workflow_id),
                    created_by,
                    json.dumps(summary, sort_keys=True),
                    json.dumps(results, sort_keys=True),
                    local_record_path,
                ),
            )

    def get_run(self, run_id: str) -> dict[str, Any] | None:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                SELECT workflow_id, project_id, created_by, summary_json, results_json, local_record_path
                FROM workflow_runs
                WHERE run_id=%s
                """,
                (run_id,),
            )
            row = cur.fetchone()
        if not row:
            return None
        return {
            "workflow_id": row[0],
            "project_id": row[1],
            "created_by": row[2],
            "summary": row[3],
            "results": row[4],
            "local_record_path": row[5],
        }

    def list_runs(self, limit: int = 50) -> list[dict[str, Any]]:
        with self.conn.cursor() as cur:
            cur.execute(
                """
                SELECT project_id, summary_json, local_record_path
                FROM workflow_runs
                ORDER BY created_at DESC
                LIMIT %s
                """,
                (limit,),
            )
            rows = cur.fetchall()
        return [{**row[1], "project_id": row[0], "local_record_path": row[2]} for row in rows]

    def reset(self) -> None:
        with self.conn.cursor() as cur:
            cur.execute("DELETE FROM workflow_runs")
            cur.execute("DELETE FROM workflow_definitions")


def _create_workflow_store() -> MemoryWorkflowStore | PgWorkflowStore:
    use_db = os.getenv("WORKFLOW_USE_DB", "true").lower() in {"1", "true", "yes", "on"}
    if use_db:
        try:
            return PgWorkflowStore()
        except Exception:
            pass
    return MemoryWorkflowStore()


_STORE = _create_workflow_store()


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


def _canonical_entry(type_id: str) -> dict[str, Any] | None:
    for entry in canonical_types.list_canonical_types():
        if entry["id"] == type_id:
            return entry
    return None


def _canonical_extension(type_id: str) -> str:
    entry = _canonical_entry(type_id)
    if entry is None:
        return ".txt"
    extensions = entry.get("extensions") or []
    if not extensions:
        return ".txt"
    return str(extensions[0])


def _canonical_sidecars(type_id: str) -> list[str]:
    entry = _canonical_entry(type_id)
    if entry is None:
        return []
    sidecars = entry.get("requires_sidecars") or []
    return [str(item) for item in sidecars]


def _workflow_run_root() -> Path:
    return Path("results") / "workflows" / "runs"


def _rule_name(node_id: str, index: int) -> str:
    safe = "".join(ch if ch.isalnum() else "_" for ch in node_id)
    return f"node_{index}_{safe}".strip("_")


def _port_output_path(node_id: str, port: str, type_id: str) -> str:
    extension = _canonical_extension(type_id)
    filename = f"{port}{extension}"
    return str(Path("artifacts") / node_id / filename)


def _workflow_edges_by_target(defn: WorkflowDefinition) -> dict[tuple[str, str], WorkflowEdge]:
    table: dict[tuple[str, str], WorkflowEdge] = {}
    for edge in defn.edges:
        table[(edge.to_node, edge.to_input)] = edge
    return table


def _project_datasets(project_id: str | None) -> list[dict[str, Any]]:
    if not project_id:
        return []
    from api import datasets

    rows: list[dict[str, Any]] = []
    for row in datasets._DATASETS:
        if row.get("project_id") == project_id:
            rows.append(dict(row))
    rows.sort(key=lambda item: (item.get("uploaded_at", 0), item.get("id", 0), item.get("filename", "")))
    return rows


def _resolve_external_inputs(project_id: str | None, defn: WorkflowDefinition) -> dict[str, Any]:
    datasets_by_type: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in _project_datasets(project_id):
        if row.get("validation_status") != "validated":
            continue
        canonical_type = row.get("canonical_type")
        if canonical_type:
            datasets_by_type[str(canonical_type)].append(row)

    edges_by_target = _workflow_edges_by_target(defn)
    bindings: dict[tuple[str, str], dict[str, Any]] = {}
    missing: list[dict[str, str]] = []
    for node in defn.nodes:
        for port, type_id in node.input_types.items():
            if edges_by_target.get((node.node_id, port)) is not None:
                continue
            candidates = datasets_by_type.get(type_id, [])
            if not candidates:
                missing.append({"node_id": node.node_id, "port": port, "type_id": type_id})
                continue
            selected = candidates[0]
            bindings[(node.node_id, port)] = {
                "dataset_id": selected["id"],
                "filename": selected["filename"],
                "canonical_type": selected["canonical_type"],
                "path": selected["uri"],
                "project_id": selected["project_id"],
            }

    return {
        "project_id": project_id,
        "bindings": bindings,
        "missing": missing,
        "validated_datasets": [row["filename"] for rows in datasets_by_type.values() for row in rows],
    }


def _resolved_parameters(node: WorkflowNode, overrides: dict[str, Any]) -> dict[str, Any]:
    merged = dict(node.parameters)
    merged.update(overrides)
    return merged


def _snakefile_rule_text(rule_name: str, node_payload: dict[str, Any]) -> str:
    inputs = node_payload["inputs"]
    outputs = node_payload["outputs"]
    lines = [f"rule {rule_name}:"]
    if inputs:
        lines.append("    input:")
        for port, input_path in inputs.items():
            lines.append(f"        {port}={input_path!r},")
    lines.append("    output:")
    for port, spec in outputs.items():
        lines.append(f"        {port}={spec['path']!r},")
    lines.append(f"    params: node_payload={repr(node_payload)}")
    lines.extend(
        [
            "    run:",
            "        import json",
            "        from pathlib import Path",
            "        payload = params.node_payload",
            "        for spec in payload['outputs'].values():",
            "            main_path = Path(spec['path'])",
            "            main_path.parent.mkdir(parents=True, exist_ok=True)",
            "            suffix = ''.join(main_path.suffixes)",
            "            body = {",
            "                'workflow_id': payload['workflow_id'],",
            "                'node_id': payload['node_id'],",
            "                'plugin_id': payload['plugin_id'],",
            "                'parameters': payload['parameters'],",
            "                'inputs': payload['inputs'],",
            "                'output_type': spec['type'],",
            "            }",
            "            if suffix == '.html':",
            "                title = f\"{payload['node_id']} report\"",
            "                content = '<html><body><h1>' + title + '</h1><pre>' + json.dumps(body, indent=2, sort_keys=True) + '</pre></body></html>'",
            "            else:",
            "                content = json.dumps(body, indent=2, sort_keys=True) + '\\n'",
            "            main_path.write_text(content, encoding='utf-8')",
            "            for sidecar in spec.get('sidecars', []):",
            "                sidecar_path = Path(str(main_path) + sidecar)",
            "                sidecar_path.write_text('placeholder\\n', encoding='utf-8')",
        ]
    )
    return "\n".join(lines)


def _render_snakefile(run_payload: dict[str, Any]) -> str:
    lines = [
        'configfile: "config.json"',
        "",
        "rule all:",
        "    input:",
    ]
    for target in run_payload["all_targets"]:
        lines.append(f"        {target!r},")
    lines.append("")
    for rule in run_payload["rules"]:
        lines.append(_snakefile_rule_text(rule["rule_name"], rule["payload"]))
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _materialize_workflow_run(
    workflow_id: str,
    defn: WorkflowDefinition,
    expanded_run: dict[str, Any],
    plan: list[str],
    run_id: str,
    external_inputs: dict[tuple[str, str], dict[str, Any]],
) -> dict[str, Any]:
    run_root = _workflow_run_root()
    run_dir = run_root / run_id / f"task-{expanded_run['index']:03d}"
    run_dir.mkdir(parents=True, exist_ok=True)

    nodes = _node_map(defn)
    edges_by_target = _workflow_edges_by_target(defn)
    output_paths: dict[tuple[str, str], str] = {}
    rules: list[dict[str, Any]] = []
    all_targets: list[str] = []

    for node_id in plan:
        node = nodes[node_id]
        outputs: dict[str, dict[str, Any]] = {}
        for port, type_id in node.output_types.items():
            output_path = _port_output_path(node.node_id, port, type_id)
            output_paths[(node.node_id, port)] = output_path
            outputs[port] = {
                "path": output_path,
                "type": type_id,
                "sidecars": _canonical_sidecars(type_id),
            }
            all_targets.append(output_path)

        inputs: dict[str, str] = {}
        external_input_meta: dict[str, dict[str, Any]] = {}
        for port, _type_id in node.input_types.items():
            edge = edges_by_target.get((node.node_id, port))
            if edge is not None:
                upstream_path = output_paths[(edge.from_node, edge.from_output)]
                inputs[port] = upstream_path
                continue
            external_binding = external_inputs.get((node.node_id, port))
            if external_binding is not None:
                inputs[port] = external_binding["path"]
                external_input_meta[port] = external_binding

        payload = {
            "workflow_id": workflow_id,
            "run_id": run_id,
            "index": expanded_run["index"],
            "node_id": node.node_id,
            "plugin_id": node.plugin_id,
            "version": node.version,
            "parameters": _resolved_parameters(node, expanded_run["overrides"].get(node.node_id, {})),
            "inputs": inputs,
            "external_inputs": external_input_meta,
            "outputs": outputs,
        }
        rules.append(
            {
                "rule_name": _rule_name(node.node_id, len(rules) + 1),
                "payload": payload,
            }
        )

    config_payload = {
        "workflow_id": workflow_id,
        "run_id": run_id,
        "index": expanded_run["index"],
        "plan": plan,
        "overrides": expanded_run["overrides"],
        "generated_at": int(time.time()),
        "rules": rules,
    }
    run_payload = {"all_targets": all_targets, "rules": rules}
    snakefile = run_dir / "Snakefile"
    configfile = run_dir / "config.json"
    snakefile.write_text(_render_snakefile(run_payload), encoding="utf-8")
    configfile.write_text(json.dumps(config_payload, indent=2, sort_keys=True), encoding="utf-8")

    return {
        "run_dir": str(run_dir.resolve()),
        "snakefile": str(snakefile.resolve()),
        "configfile": str(configfile.resolve()),
        "all_targets": [str((run_dir / target).resolve()) for target in all_targets],
        "rules": rules,
    }


def _resolve_snakemake_command() -> list[str] | None:
    raw = os.getenv("WORKFLOW_SNAKEMAKE_BIN", "snakemake").strip() or "snakemake"
    parts = shlex.split(raw)
    if not parts:
        return None
    command = shutil.which(parts[0])
    if command is None:
        return None
    parts[0] = command
    return parts


def _execute_fallback(materialized: dict[str, Any]) -> dict[str, Any]:
    start = time.time()
    for rule in materialized["rules"]:
        payload = rule["payload"]
        for spec in payload["outputs"].values():
            main_path = Path(materialized["run_dir"]) / spec["path"]
            main_path.parent.mkdir(parents=True, exist_ok=True)
            body = {
                "workflow_id": payload["workflow_id"],
                "node_id": payload["node_id"],
                "plugin_id": payload["plugin_id"],
                "parameters": payload["parameters"],
                "inputs": payload["inputs"],
                "output_type": spec["type"],
            }
            suffix = "".join(main_path.suffixes)
            if suffix == ".html":
                content = "<html><body><pre>" + json.dumps(body, indent=2, sort_keys=True) + "</pre></body></html>"
            else:
                content = json.dumps(body, indent=2, sort_keys=True) + "\n"
            main_path.write_text(content, encoding="utf-8")
            for sidecar in spec.get("sidecars", []):
                Path(str(main_path) + sidecar).write_text("placeholder\n", encoding="utf-8")
    return {
        "engine": "python-fallback",
        "snakemake_available": False,
        "status": "completed",
        "dry_run": {"status": "unavailable", "stdout": "", "stderr": "snakemake command not found"},
        "run": {"status": "completed", "stdout": "", "stderr": "", "duration_ms": int((time.time() - start) * 1000)},
    }


def _execute_snakemake(materialized: dict[str, Any]) -> dict[str, Any]:
    command = _resolve_snakemake_command()
    if command is None:
        return _execute_fallback(materialized)

    cwd = materialized["run_dir"]
    dry_start = time.time()
    dry_run = subprocess.run(
        command + ["--snakefile", "Snakefile", "--configfile", "config.json", "--cores", "1", "--dry-run"],
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
        timeout=60,
    )
    dry_meta = {
        "status": "completed" if dry_run.returncode == 0 else "failed",
        "stdout": (dry_run.stdout or "")[-4000:],
        "stderr": (dry_run.stderr or "")[-4000:],
        "duration_ms": int((time.time() - dry_start) * 1000),
    }
    if dry_run.returncode != 0:
        return {
            "engine": "snakemake",
            "snakemake_available": True,
            "status": "failed",
            "dry_run": dry_meta,
            "run": {"status": "skipped", "stdout": "", "stderr": "dry-run failed", "duration_ms": 0},
        }

    run_start = time.time()
    completed = subprocess.run(
        command + ["--snakefile", "Snakefile", "--configfile", "config.json", "--cores", "1"],
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    return {
        "engine": "snakemake",
        "snakemake_available": True,
        "status": "completed" if completed.returncode == 0 else "failed",
        "dry_run": dry_meta,
        "run": {
            "status": "completed" if completed.returncode == 0 else "failed",
            "stdout": (completed.stdout or "")[-4000:],
            "stderr": (completed.stderr or "")[-4000:],
            "duration_ms": int((time.time() - run_start) * 1000),
        },
    }


def _get_workflow_or_404(workflow_id: str) -> dict[str, Any]:
    item = _STORE.get_workflow(workflow_id)
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
    if _STORE.get_workflow(defn.workflow_id) is not None:
        raise HTTPException(status_code=409, detail="Workflow already exists")
    try:
        report = _validation_report(defn)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    workflow_payload = defn.model_dump()
    local_export_path = _write_workflow_definition_file(
        workflow_id=defn.workflow_id,
        workflow=workflow_payload,
        report=report,
        created_by=user["email"],
    )
    created = _STORE.create_workflow(
        workflow_id=defn.workflow_id,
        workflow=workflow_payload,
        report=report,
        created_by=user["email"],
        local_export_path=local_export_path,
    )
    if not created:
        raise HTTPException(status_code=409, detail="Workflow already exists")
    return {"workflow_id": defn.workflow_id, "report": report, "local_export_path": _host_visible_path(local_export_path)}


@router.get("")
def list_workflows(user=Depends(current_user)):
    rows = []
    for item in _STORE.list_workflows():
        item = _present_workflow_row(item)
        report = item["report"]
        rows.append(
            {
                "workflow_id": item["workflow_id"],
                "created_by": item["created_by"],
                "node_count": report["node_count"],
                "edge_count": report["edge_count"],
                "local_export_path": item.get("local_export_path"),
            }
        )
    return {"workflows": rows}


@router.get("/{workflow_id}/export")
def export_workflow(workflow_id: str, user=Depends(current_user)):
    item = _present_workflow_row(_get_workflow_or_404(workflow_id))
    return {"workflow": item["workflow"], "local_export_path": item.get("local_export_path")}


@router.post("/{workflow_id}/validate")
def validate_workflow(workflow_id: str, user=Depends(current_user)):
    item = _get_workflow_or_404(workflow_id)
    defn = WorkflowDefinition.model_validate(item["workflow"])
    try:
        report = _validation_report(defn)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    local_export_path = _write_workflow_definition_file(
        workflow_id=workflow_id,
        workflow=item["workflow"],
        report=report,
        created_by=item["created_by"],
    )
    _STORE.update_workflow_report(workflow_id, report)
    return {"report": report, "local_export_path": _host_visible_path(local_export_path)}


@router.get("/{workflow_id}/plan")
def workflow_plan(workflow_id: str, user=Depends(current_user)):
    item = _present_workflow_row(_get_workflow_or_404(workflow_id))
    return {
        "workflow_id": workflow_id,
        "topological_order": item["report"]["topological_order"],
        "branching_nodes": item["report"]["branching_nodes"],
        "local_export_path": item.get("local_export_path"),
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
    project_id = _workflow_project_id(workflow_id)
    external_input_state = _resolve_external_inputs(project_id, defn)
    if external_input_state["missing"]:
        raise HTTPException(
            status_code=400,
            detail={
                "message": "Workflow inputs are missing or not validated for this project",
                "project_id": project_id,
                "missing_inputs": external_input_state["missing"],
                "validated_datasets": external_input_state["validated_datasets"],
            },
        )
    start = time.time()
    run_id = f"{workflow_id}-dist-{int(start * 1000)}"

    def _execute_one(expanded_run: dict[str, Any]) -> dict[str, Any]:
        materialized = _materialize_workflow_run(
            workflow_id,
            defn,
            expanded_run,
            plan,
            run_id,
            external_input_state["bindings"],
        )
        execution = _execute_snakemake(materialized)
        payload = {
            "workflow_id": workflow_id,
            "project_id": project_id,
            "index": expanded_run["index"],
            "overrides": expanded_run["overrides"],
            "plan": plan,
            "external_inputs": {
                f"{node_id}.{port}": binding["filename"]
                for (node_id, port), binding in sorted(external_input_state["bindings"].items())
            },
            "workflow_artifacts": {
                "run_dir": materialized["run_dir"],
                "snakefile": materialized["snakefile"],
                "configfile": materialized["configfile"],
            },
        }
        return {
            "index": expanded_run["index"],
            "status": execution["status"],
            "plan": plan,
            "overrides": expanded_run["overrides"],
            "signature": _sha256_json(payload),
            "engine": execution["engine"],
            "snakemake_available": execution["snakemake_available"],
            "workflow_artifacts": payload["workflow_artifacts"],
            "external_inputs": payload["external_inputs"],
            "dry_run": execution["dry_run"],
            "run": execution["run"],
            "outputs": materialized["all_targets"],
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
        "project_id": _workflow_project_id(workflow_id),
        "created_by": user["email"],
        "submitted_runs": len(expanded),
        "completed_runs": len([row for row in results if row["status"] == "completed"]),
        "failed_runs": len([row for row in results if row["status"] != "completed"]),
        "max_workers": body.max_workers,
        "duration_ms": int((time.time() - start) * 1000),
    }
    local_record_path = _write_run_record_file(
        run_id=run_id,
        workflow_id=workflow_id,
        created_by=user["email"],
        summary=summary,
        results=results,
    )
    _STORE.save_run(
        run_id=run_id,
        workflow_id=workflow_id,
        created_by=user["email"],
        summary=summary,
        results=results,
        local_record_path=local_record_path,
    )
    return {
        "summary": _present_run_summary({**summary, "local_record_path": local_record_path}),
        "results": [_present_run_result(row) for row in results],
    }


@router.get("/runs")
def list_distributed_runs(limit: int = 50, user=Depends(current_user)):
    rows: list[dict[str, Any]] = []
    for row in _STORE.list_runs(limit=limit):
        row = _present_run_summary(row)
        if user["role"] != "admin" and row.get("created_by") != user["email"]:
            continue
        rows.append(row)
        if len(rows) >= limit:
            break
    return {"runs": rows}


@router.get("/runs/{run_id}")
def get_distributed_run(run_id: str, user=Depends(current_user)):
    item = _STORE.get_run(run_id)
    if item is None:
        raise HTTPException(status_code=404, detail="Workflow run not found")
    if user["role"] != "admin" and item.get("created_by") != user["email"]:
        raise HTTPException(status_code=403, detail="Forbidden")
    return _present_run_record(item)


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
    if _STORE.get_workflow(workflow_id) is not None:
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
    workflow_payload = defn.model_dump()
    local_export_path = _write_workflow_definition_file(
        workflow_id=defn.workflow_id,
        workflow=workflow_payload,
        report=report,
        created_by=user["email"],
    )
    created = _STORE.create_workflow(
        workflow_id=defn.workflow_id,
        workflow=workflow_payload,
        report=report,
        created_by=user["email"],
        local_export_path=local_export_path,
    )
    if not created:
        raise HTTPException(status_code=409, detail="Workflow already exists")
    return {"workflow_id": defn.workflow_id, "report": report, "local_export_path": _host_visible_path(local_export_path)}


def reset_workflow_store_for_tests() -> None:
    _STORE.reset()
    shutil.rmtree(_workflow_definition_root(), ignore_errors=True)
    shutil.rmtree(_workflow_run_root(), ignore_errors=True)
