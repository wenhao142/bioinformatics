"use client";

import { useMemo, useState } from "react";

type PluginItem = {
  plugin_id: string;
  name: string;
  version: string;
  enabled?: boolean;
  tags?: string[];
};

type NodeDraft = {
  node_id: string;
  plugin_id: string;
  version: string;
  input_types: string;
  output_types: string;
  parameters: string;
};

type EdgeDraft = {
  from_node: string;
  from_output: string;
  to_node: string;
  to_input: string;
};

type RunMonitor = {
  run_id: string;
  submitted_runs: number;
  completed_runs: number;
  duration_ms: number;
};

const defaultNodeJson = {
  input_types: "{}",
  output_types: '{"out":"report.html"}',
  parameters: "{}",
};

function parseJsonObject(raw: string): Record<string, unknown> {
  const parsed = JSON.parse(raw);
  if (!parsed || typeof parsed !== "object" || Array.isArray(parsed)) {
    throw new Error("Must be a JSON object.");
  }
  return parsed as Record<string, unknown>;
}

function downloadText(filename: string, content: string): void {
  const blob = new Blob([content], { type: "text/plain;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export default function WorkflowBuilderPage() {
  const apiBase = process.env.NEXT_PUBLIC_API_URL || "http://localhost:18000";

  const [email, setEmail] = useState("admin@example.com");
  const [password, setPassword] = useState("password");
  const [token, setToken] = useState<string | null>(null);
  const [plugins, setPlugins] = useState<PluginItem[]>([]);
  const [workflowId, setWorkflowId] = useState("wf-builder-demo");
  const [nodes, setNodes] = useState<NodeDraft[]>([]);
  const [edges, setEdges] = useState<EdgeDraft[]>([]);
  const [sweepsRaw, setSweepsRaw] = useState('{"n2.alpha":[0.1,0.2]}');
  const [status, setStatus] = useState<string>("");
  const [monitor, setMonitor] = useState<RunMonitor | null>(null);
  const [runDetails, setRunDetails] = useState<unknown>(null);
  const [edgeDraft, setEdgeDraft] = useState<EdgeDraft>({ from_node: "", from_output: "", to_node: "", to_input: "" });

  const nodePortMap = useMemo(() => {
    const table: Record<string, { outputs: string[]; inputs: string[] }> = {};
    for (const node of nodes) {
      let outputs: string[] = [];
      let inputs: string[] = [];
      try {
        outputs = Object.keys(parseJsonObject(node.output_types));
      } catch {
        outputs = [];
      }
      try {
        inputs = Object.keys(parseJsonObject(node.input_types));
      } catch {
        inputs = [];
      }
      table[node.node_id] = { outputs, inputs };
    }
    return table;
  }, [nodes]);

  const authHeaders = token
    ? {
        Authorization: `Bearer ${token}`,
        "Content-Type": "application/json",
      }
    : null;

  async function login() {
    setStatus("Signing in...");
    const res = await fetch(`${apiBase}/auth/login`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ email, password }),
    });
    if (!res.ok) {
      setStatus(`Login failed: ${res.status}`);
      return;
    }
    const body = await res.json();
    setToken(body.access_token);
    setStatus("Signed in.");
  }

  async function loadPlugins() {
    if (!authHeaders) {
      setStatus("Please sign in first.");
      return;
    }
    const res = await fetch(`${apiBase}/plugins`, { headers: authHeaders });
    if (!res.ok) {
      setStatus(`Load plugins failed: ${res.status}`);
      return;
    }
    const body = await res.json();
    setPlugins(body.plugins || []);
    setStatus(`Loaded ${body.plugins?.length || 0} plugins.`);
  }

  function addNode(plugin: PluginItem) {
    const suffix = nodes.length + 1;
    const nodeId = `n${suffix}`;
    setNodes((prev) => [
      ...prev,
      {
        node_id: nodeId,
        plugin_id: plugin.plugin_id,
        version: plugin.version,
        input_types: defaultNodeJson.input_types,
        output_types: defaultNodeJson.output_types,
        parameters: defaultNodeJson.parameters,
      },
    ]);
    setStatus(`Added node ${nodeId} (${plugin.plugin_id})`);
  }

  function onDropPlugin(event: React.DragEvent<HTMLDivElement>) {
    event.preventDefault();
    const pluginId = event.dataTransfer.getData("text/plain");
    const plugin = plugins.find((p) => p.plugin_id === pluginId);
    if (!plugin) return;
    addNode(plugin);
  }

  function updateNode(nodeId: string, field: keyof NodeDraft, value: string) {
    setNodes((prev) => prev.map((n) => (n.node_id === nodeId ? { ...n, [field]: value } : n)));
  }

  function addEdge() {
    if (!edgeDraft.from_node || !edgeDraft.to_node || !edgeDraft.from_output || !edgeDraft.to_input) {
      setStatus("Choose from/to node and ports first.");
      return;
    }
    setEdges((prev) => [...prev, edgeDraft]);
    setEdgeDraft({ from_node: "", from_output: "", to_node: "", to_input: "" });
    setStatus("Edge added.");
  }

  function buildWorkflowPayload() {
    const workflowNodes = nodes.map((n) => ({
      node_id: n.node_id,
      plugin_id: n.plugin_id,
      version: n.version || undefined,
      input_types: parseJsonObject(n.input_types),
      output_types: parseJsonObject(n.output_types),
      parameters: parseJsonObject(n.parameters),
    }));
    const parameter_sweeps = parseJsonObject(sweepsRaw);
    return {
      workflow_id: workflowId,
      nodes: workflowNodes,
      edges,
      parameter_sweeps,
    };
  }

  async function saveAndValidateWorkflow() {
    if (!authHeaders) {
      setStatus("Please sign in first.");
      return;
    }
    let payload: unknown;
    try {
      payload = buildWorkflowPayload();
    } catch (err) {
      setStatus(`JSON error: ${(err as Error).message}`);
      return;
    }
    const create = await fetch(`${apiBase}/workflows/import`, {
      method: "POST",
      headers: authHeaders,
      body: JSON.stringify(payload),
    });
    if (create.status === 409) {
      const validate = await fetch(`${apiBase}/workflows/${workflowId}/validate`, {
        method: "POST",
        headers: authHeaders,
      });
      if (!validate.ok) {
        const text = await validate.text();
        setStatus(`Validate failed: ${validate.status} ${text}`);
        return;
      }
      setStatus("Workflow already existed; validation refreshed.");
      return;
    }
    if (!create.ok) {
      const text = await create.text();
      setStatus(`Import failed: ${create.status} ${text}`);
      return;
    }
    setStatus("Workflow imported and validated.");
  }

  async function runDistributed() {
    if (!authHeaders) {
      setStatus("Please sign in first.");
      return;
    }
    const res = await fetch(`${apiBase}/workflows/${workflowId}/execute/distributed`, {
      method: "POST",
      headers: authHeaders,
      body: JSON.stringify({ max_workers: 2 }),
    });
    if (!res.ok) {
      const text = await res.text();
      setStatus(`Run failed: ${res.status} ${text}`);
      return;
    }
    const body = await res.json();
    setMonitor(body.summary);
    setRunDetails(body);
    setStatus(`Distributed run completed: ${body.summary.run_id}`);
  }

  async function refreshRun() {
    if (!authHeaders || !monitor?.run_id) return;
    const res = await fetch(`${apiBase}/workflows/runs/${monitor.run_id}`, { headers: authHeaders });
    if (!res.ok) {
      setStatus(`Refresh run failed: ${res.status}`);
      return;
    }
    const body = await res.json();
    setRunDetails(body);
    setStatus("Run details refreshed.");
  }

  async function exportReproReport() {
    if (!authHeaders) {
      setStatus("Please sign in first.");
      return;
    }
    const exportRes = await fetch(`${apiBase}/workflows/${workflowId}/export`, { headers: authHeaders });
    if (!exportRes.ok) {
      setStatus(`Export failed: ${exportRes.status}`);
      return;
    }
    const exportBody = await exportRes.json();
    const reportLines = [
      "# Workflow Reproducibility Report",
      "",
      `- workflow_id: ${workflowId}`,
      `- generated_at: ${new Date().toISOString()}`,
      `- node_count: ${exportBody.workflow.nodes?.length || 0}`,
      `- edge_count: ${exportBody.workflow.edges?.length || 0}`,
      monitor ? `- last_run_id: ${monitor.run_id}` : "- last_run_id: n/a",
      "",
      "## Workflow JSON",
      "```json",
      JSON.stringify(exportBody.workflow, null, 2),
      "```",
      "",
      "## Last Run",
      "```json",
      JSON.stringify(runDetails || {}, null, 2),
      "```",
    ];
    downloadText(`${workflowId}.reproducibility.md`, reportLines.join("\n"));
    setStatus("Reproducibility report exported.");
  }

  return (
    <main style={{ padding: "1rem", maxWidth: 1200, margin: "0 auto", fontFamily: "ui-sans-serif, system-ui" }}>
      <h1 style={{ margin: "0 0 0.5rem 0" }}>Workflow Builder</h1>
      <p style={{ marginTop: 0 }}>
        Build workflows with drag-and-drop nodes, configure parameters, run distributed execution, and export reproducibility
        reports.
      </p>

      <section style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit,minmax(260px,1fr))", gap: 12, marginBottom: 12 }}>
        <div style={{ border: "1px solid #d0d7de", borderRadius: 8, padding: 12 }}>
          <h3 style={{ marginTop: 0 }}>Session</h3>
          <label style={{ display: "block", marginBottom: 6 }}>
            Email
            <input value={email} onChange={(e) => setEmail(e.target.value)} style={{ width: "100%" }} />
          </label>
          <label style={{ display: "block", marginBottom: 6 }}>
            Password
            <input type="password" value={password} onChange={(e) => setPassword(e.target.value)} style={{ width: "100%" }} />
          </label>
          <button onClick={() => void login()}>Sign In</button>
          <button onClick={() => void loadPlugins()} style={{ marginLeft: 8 }}>
            Load Tools
          </button>
        </div>

        <div style={{ border: "1px solid #d0d7de", borderRadius: 8, padding: 12 }}>
          <h3 style={{ marginTop: 0 }}>Workflow</h3>
          <label style={{ display: "block", marginBottom: 6 }}>
            Workflow ID
            <input value={workflowId} onChange={(e) => setWorkflowId(e.target.value)} style={{ width: "100%" }} />
          </label>
          <label style={{ display: "block", marginBottom: 6 }}>
            Parameter Sweeps (JSON)
            <textarea value={sweepsRaw} onChange={(e) => setSweepsRaw(e.target.value)} rows={4} style={{ width: "100%" }} />
          </label>
          <button onClick={() => void saveAndValidateWorkflow()}>Validate & Save</button>
          <button onClick={() => void runDistributed()} style={{ marginLeft: 8 }}>
            Run Distributed
          </button>
          <button onClick={() => void exportReproReport()} style={{ marginLeft: 8 }}>
            Export Repro Report
          </button>
        </div>
      </section>

      <section style={{ display: "grid", gridTemplateColumns: "280px 1fr", gap: 12 }}>
        <aside style={{ border: "1px solid #d0d7de", borderRadius: 8, padding: 12 }}>
          <h3 style={{ marginTop: 0 }}>Tool Palette (Drag)</h3>
          {plugins.map((p) => (
            <div
              key={`${p.plugin_id}:${p.version}`}
              draggable
              onDragStart={(e) => e.dataTransfer.setData("text/plain", p.plugin_id)}
              onDoubleClick={() => addNode(p)}
              style={{
                border: "1px solid #8c959f",
                borderRadius: 6,
                padding: 8,
                marginBottom: 8,
                cursor: "grab",
                background: "#f6f8fa",
              }}
            >
              <strong>{p.plugin_id}</strong>
              <div style={{ fontSize: 12 }}>v{p.version}</div>
              <div style={{ fontSize: 12 }}>{(p.tags || []).join(", ") || "no-tags"}</div>
            </div>
          ))}
        </aside>

        <div
          onDragOver={(e) => e.preventDefault()}
          onDrop={onDropPlugin}
          style={{ border: "2px dashed #8c959f", borderRadius: 8, padding: 12, minHeight: 360 }}
        >
          <h3 style={{ marginTop: 0 }}>Canvas (Drop tool here)</h3>
          {nodes.length === 0 ? <p>No nodes yet. Drag a tool from the left panel.</p> : null}
          {nodes.map((node) => (
            <div key={node.node_id} style={{ border: "1px solid #d0d7de", borderRadius: 6, padding: 8, marginBottom: 8 }}>
              <div style={{ display: "flex", gap: 8, marginBottom: 6 }}>
                <input
                  value={node.node_id}
                  onChange={(e) => updateNode(node.node_id, "node_id", e.target.value)}
                  style={{ width: 90 }}
                />
                <input
                  value={node.plugin_id}
                  onChange={(e) => updateNode(node.node_id, "plugin_id", e.target.value)}
                  style={{ flex: 1 }}
                />
                <input
                  value={node.version}
                  onChange={(e) => updateNode(node.node_id, "version", e.target.value)}
                  style={{ width: 100 }}
                />
              </div>
              <label style={{ display: "block", marginBottom: 4 }}>
                Input types (JSON)
                <textarea
                  value={node.input_types}
                  onChange={(e) => updateNode(node.node_id, "input_types", e.target.value)}
                  rows={2}
                  style={{ width: "100%" }}
                />
              </label>
              <label style={{ display: "block", marginBottom: 4 }}>
                Output types (JSON)
                <textarea
                  value={node.output_types}
                  onChange={(e) => updateNode(node.node_id, "output_types", e.target.value)}
                  rows={2}
                  style={{ width: "100%" }}
                />
              </label>
              <label style={{ display: "block" }}>
                Parameters (JSON)
                <textarea
                  value={node.parameters}
                  onChange={(e) => updateNode(node.node_id, "parameters", e.target.value)}
                  rows={2}
                  style={{ width: "100%" }}
                />
              </label>
            </div>
          ))}
        </div>
      </section>

      <section style={{ border: "1px solid #d0d7de", borderRadius: 8, padding: 12, marginTop: 12 }}>
        <h3 style={{ marginTop: 0 }}>Connections</h3>
        <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit,minmax(180px,1fr))", gap: 8 }}>
          <select value={edgeDraft.from_node} onChange={(e) => setEdgeDraft({ ...edgeDraft, from_node: e.target.value, from_output: "" })}>
            <option value="">From node</option>
            {nodes.map((n) => (
              <option key={n.node_id} value={n.node_id}>
                {n.node_id}
              </option>
            ))}
          </select>
          <select value={edgeDraft.from_output} onChange={(e) => setEdgeDraft({ ...edgeDraft, from_output: e.target.value })}>
            <option value="">From output</option>
            {(nodePortMap[edgeDraft.from_node]?.outputs || []).map((port) => (
              <option key={port} value={port}>
                {port}
              </option>
            ))}
          </select>
          <select value={edgeDraft.to_node} onChange={(e) => setEdgeDraft({ ...edgeDraft, to_node: e.target.value, to_input: "" })}>
            <option value="">To node</option>
            {nodes.map((n) => (
              <option key={n.node_id} value={n.node_id}>
                {n.node_id}
              </option>
            ))}
          </select>
          <select value={edgeDraft.to_input} onChange={(e) => setEdgeDraft({ ...edgeDraft, to_input: e.target.value })}>
            <option value="">To input</option>
            {(nodePortMap[edgeDraft.to_node]?.inputs || []).map((port) => (
              <option key={port} value={port}>
                {port}
              </option>
            ))}
          </select>
        </div>
        <button onClick={addEdge} style={{ marginTop: 8 }}>
          Add Edge
        </button>
        <pre style={{ background: "#f6f8fa", padding: 8, marginTop: 8, overflowX: "auto" }}>{JSON.stringify(edges, null, 2)}</pre>
      </section>

      <section style={{ border: "1px solid #d0d7de", borderRadius: 8, padding: 12, marginTop: 12 }}>
        <h3 style={{ marginTop: 0 }}>Run Monitor</h3>
        <button onClick={() => void refreshRun()} disabled={!monitor?.run_id}>
          Refresh Last Run
        </button>
        {monitor ? (
          <div style={{ marginTop: 8 }}>
            <div>Run ID: {monitor.run_id}</div>
            <div>
              Completed: {monitor.completed_runs}/{monitor.submitted_runs}
            </div>
            <div>Duration: {monitor.duration_ms} ms</div>
          </div>
        ) : (
          <p>No run yet.</p>
        )}
      </section>

      <p style={{ marginTop: 12, color: "#14532d" }}>{status}</p>
    </main>
  );
}
