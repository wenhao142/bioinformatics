"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import { useRouter } from "next/navigation";

type PluginItem = {
  plugin_id: string;
  name: string;
  version: string;
  enabled?: boolean;
  tags?: string[];
};

type PaletteTool = PluginItem & {
  source: "builtin" | "registry";
  default_input_types: Record<string, string>;
  default_output_types: Record<string, string>;
  default_parameters: Record<string, unknown>;
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

type SimpleEdgeDraft = {
  from_node: string;
  to_node: string;
};

type RunMonitor = {
  run_id: string;
  submitted_runs: number;
  completed_runs: number;
  duration_ms: number;
};

type WorkflowListItem = {
  workflow_id: string;
  node_count: number;
  edge_count: number;
  created_by: string;
};

const BUILTIN_TOOLS: PaletteTool[] = [
  {
    plugin_id: "fasta_qc",
    name: "FASTA QC",
    version: "0.1.0",
    tags: ["wgs", "qc", "reference"],
    source: "builtin",
    default_input_types: {},
    default_output_types: { report: "report.html" },
    default_parameters: { mode: "fasta-basic-qc" },
  },
  {
    plugin_id: "fastqc",
    name: "FastQC",
    version: "0.1.0",
    tags: ["rna-seq", "wgs", "qc"],
    source: "builtin",
    default_input_types: { reads: "reads.fastq.gz" },
    default_output_types: { report: "report.html" },
    default_parameters: { threads: 2 },
  },
  {
    plugin_id: "cutadapt",
    name: "Cutadapt",
    version: "0.1.0",
    tags: ["rna-seq", "wgs", "trim"],
    source: "builtin",
    default_input_types: { reads: "reads.fastq.gz" },
    default_output_types: { trimmed_reads: "reads.fastq.gz" },
    default_parameters: { quality_cutoff: 20, min_length: 20 },
  },
  {
    plugin_id: "star_align",
    name: "STAR Align",
    version: "0.1.0",
    tags: ["rna-seq", "align"],
    source: "builtin",
    default_input_types: { reads: "reads.fastq.gz" },
    default_output_types: { alignment: "align.bam" },
    default_parameters: { threads: 8 },
  },
  {
    plugin_id: "samtools_sort_index",
    name: "samtools sort+index",
    version: "0.1.0",
    tags: ["rna-seq", "wgs", "align"],
    source: "builtin",
    default_input_types: { alignment: "align.bam" },
    default_output_types: { alignment_sorted: "align.bam" },
    default_parameters: { threads: 4 },
  },
  {
    plugin_id: "featurecounts",
    name: "featureCounts",
    version: "0.1.0",
    tags: ["rna-seq", "counts"],
    source: "builtin",
    default_input_types: { alignment_sorted: "align.bam" },
    default_output_types: { counts: "expression.counts.tsv" },
    default_parameters: { threads: 4, feature_type: "exon", attribute: "gene_id" },
  },
  {
    plugin_id: "deseq2_diffexp",
    name: "DESeq2",
    version: "0.1.0",
    tags: ["rna-seq", "diffexp"],
    source: "builtin",
    default_input_types: { counts: "expression.counts.tsv" },
    default_output_types: { diff_table: "expression.diff_table.tsv", report: "report.html" },
    default_parameters: { design_formula: "~ condition", alpha: 0.05 },
  },
  {
    plugin_id: "bcftools_call",
    name: "bcftools call",
    version: "0.1.0",
    tags: ["wgs", "variant"],
    source: "builtin",
    default_input_types: { alignment_sorted: "align.bam" },
    default_output_types: { variants: "variants.vcf.gz" },
    default_parameters: { ploidy: 2 },
  },
  {
    plugin_id: "variant_report",
    name: "Variant Report",
    version: "0.1.0",
    tags: ["wgs", "report"],
    source: "builtin",
    default_input_types: { variants: "variants.vcf.gz" },
    default_output_types: { report: "report.html" },
    default_parameters: {},
  },
];

const CANVAS_WIDTH = 1800;
const CANVAS_HEIGHT = 980;
const NODE_WIDTH = 220;
const NODE_HEIGHT = 76;
const NODE_MARGIN = 24;

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
  const router = useRouter();
  const apiBase = process.env.NEXT_PUBLIC_API_URL || "http://localhost:18000";

  const [token, setToken] = useState<string | null>(null);
  const [authReady, setAuthReady] = useState(false);
  const [plugins, setPlugins] = useState<PluginItem[]>([]);
  const [projectId, setProjectId] = useState("demo-project");
  const [workflowName, setWorkflowName] = useState("lego-flow");
  const [nodes, setNodes] = useState<NodeDraft[]>([]);
  const [edges, setEdges] = useState<EdgeDraft[]>([]);
  const [sweepsRaw, setSweepsRaw] = useState('{"n2.alpha":[0.1,0.2]}');
  const [status, setStatus] = useState<string>("");
  const [error, setError] = useState<string>("");
  const [monitor, setMonitor] = useState<RunMonitor | null>(null);
  const [runDetails, setRunDetails] = useState<unknown>(null);
  const [simpleEdgeDraft, setSimpleEdgeDraft] = useState<SimpleEdgeDraft>({ from_node: "", to_node: "" });
  const [savedWorkflows, setSavedWorkflows] = useState<WorkflowListItem[]>([]);
  const [selectedSavedWorkflow, setSelectedSavedWorkflow] = useState("");
  const [paletteCategory, setPaletteCategory] = useState<"all" | "wgs" | "rna-seq">("all");
  const [paletteQuery, setPaletteQuery] = useState("");
  const [canvasZoom, setCanvasZoom] = useState(1);
  const [nodePositions, setNodePositions] = useState<Record<string, { x: number; y: number }>>({});
  const [isDraggingNode, setIsDraggingNode] = useState(false);
  const [linkStartNodeId, setLinkStartNodeId] = useState<string | null>(null);
  const [linkPreviewPoint, setLinkPreviewPoint] = useState<{ x: number; y: number } | null>(null);
  const canvasRef = useRef<HTMLDivElement | null>(null);
  const dragContextRef = useRef<{
    nodeId: string;
    startClientX: number;
    startClientY: number;
    offsetX: number;
    offsetY: number;
    moved: boolean;
  } | null>(null);

  const workflowId = useMemo(() => {
    const p = projectId.trim();
    const w = workflowName.trim();
    if (!p || !w) return "";
    return `${p}__${w}`;
  }, [projectId, workflowName]);

  const draftStorageKey = useMemo(() => {
    if (!workflowId) return "";
    return `workflow-builder:draft:${workflowId}`;
  }, [workflowId]);

  const paletteTools = useMemo(() => {
    const builtinById = new Map(BUILTIN_TOOLS.map((tool) => [tool.plugin_id, tool] as const));
    const fromRegistry: PaletteTool[] = plugins
      .filter((plugin) => plugin.enabled !== false)
      .map((plugin) => {
        const fallback = builtinById.get(plugin.plugin_id);
        return {
          ...plugin,
          source: "registry",
          default_input_types: fallback?.default_input_types || {},
          default_output_types: fallback?.default_output_types || { out: "report.html" },
          default_parameters: fallback?.default_parameters || {},
        };
      });
    const merged = [...fromRegistry];
    const seen = new Set(fromRegistry.map((tool) => tool.plugin_id));
    for (const tool of BUILTIN_TOOLS) {
      if (!seen.has(tool.plugin_id)) {
        merged.push(tool);
      }
    }
    return merged;
  }, [plugins]);

  const filteredPaletteTools = useMemo(() => {
    const query = paletteQuery.trim().toLowerCase();
    return paletteTools.filter((tool) => {
      if (paletteCategory !== "all") {
        const tags = tool.tags || [];
        if (!tags.includes(paletteCategory)) return false;
      }
      if (!query) return true;
      const haystack = `${tool.plugin_id} ${tool.name} ${(tool.tags || []).join(" ")}`.toLowerCase();
      return haystack.includes(query);
    });
  }, [paletteTools, paletteCategory, paletteQuery]);

  const toolMetaById = useMemo(() => {
    const table = new Map<string, PaletteTool>();
    for (const tool of paletteTools) {
      table.set(tool.plugin_id, tool);
    }
    return table;
  }, [paletteTools]);

  const parsedNodeTypes = useMemo(() => {
    const outputs: Record<string, Record<string, string>> = {};
    const inputs: Record<string, Record<string, string>> = {};
    const errors: string[] = [];
    for (const node of nodes) {
      if (!node.node_id.trim()) {
        errors.push("Node ID cannot be empty.");
      }
      if (!node.plugin_id.trim()) {
        errors.push(`Node ${node.node_id || "<empty>"} missing plugin_id.`);
      }
      try {
        const outObj = parseJsonObject(node.output_types);
        outputs[node.node_id] = {};
        for (const [k, v] of Object.entries(outObj)) {
          if (typeof v !== "string" || !v.trim()) {
            errors.push(`Node ${node.node_id} output_types.${k} must be non-empty string canonical type.`);
          } else {
            outputs[node.node_id][k] = v.trim();
          }
        }
      } catch (err) {
        errors.push(`Node ${node.node_id || "<empty>"} output_types JSON invalid: ${(err as Error).message}`);
      }
      try {
        const inObj = parseJsonObject(node.input_types);
        inputs[node.node_id] = {};
        for (const [k, v] of Object.entries(inObj)) {
          if (typeof v !== "string" || !v.trim()) {
            errors.push(`Node ${node.node_id} input_types.${k} must be non-empty string canonical type.`);
          } else {
            inputs[node.node_id][k] = v.trim();
          }
        }
      } catch (err) {
        errors.push(`Node ${node.node_id || "<empty>"} input_types JSON invalid: ${(err as Error).message}`);
      }
      try {
        parseJsonObject(node.parameters);
      } catch (err) {
        errors.push(`Node ${node.node_id || "<empty>"} parameters JSON invalid: ${(err as Error).message}`);
      }
    }
    return { outputs, inputs, errors };
  }, [nodes]);

  const canvasNodes = useMemo(() => {
    return nodes.map((node, idx) => {
      const fallbackX = NODE_MARGIN + (idx % 5) * 280;
      const fallbackY = NODE_MARGIN + Math.floor(idx / 5) * 130;
      const pos = nodePositions[node.node_id] || { x: fallbackX, y: fallbackY };
      const meta = toolMetaById.get(node.plugin_id);
      const tags = meta?.tags || [];
      let tone = "node-generic";
      if (meta?.source === "registry") {
        tone = "node-registry";
      } else if (tags.includes("wgs") && tags.includes("rna-seq")) {
        tone = "node-mixed";
      } else if (tags.includes("wgs")) {
        tone = "node-wgs";
      } else if (tags.includes("rna-seq")) {
        tone = "node-rna";
      }
      return {
        nodeId: node.node_id,
        pluginId: node.plugin_id,
        title: meta?.name || node.plugin_id,
        tone,
        x: pos.x,
        y: pos.y,
      };
    });
  }, [nodes, nodePositions, toolMetaById]);

  const canvasNodeById = useMemo(() => {
    return new Map(canvasNodes.map((node) => [node.nodeId, node] as const));
  }, [canvasNodes]);

  const canvasEdgePaths = useMemo(() => {
    const pairCounter = new Map<string, number>();
    return edges
      .map((edge) => {
        const fromNode = canvasNodeById.get(edge.from_node);
        const toNode = canvasNodeById.get(edge.to_node);
        if (!fromNode || !toNode) return null;

        const laneKey = `${edge.from_node}->${edge.to_node}`;
        const lane = pairCounter.get(laneKey) || 0;
        pairCounter.set(laneKey, lane + 1);
        const laneShift = (lane - 1) * 10;

        const x1 = fromNode.x + NODE_WIDTH - 2;
        const y1 = fromNode.y + NODE_HEIGHT / 2 + laneShift;
        const x2 = toNode.x + 2;
        const y2 = toNode.y + NODE_HEIGHT / 2 + laneShift;
        const cp = Math.max(48, Math.abs(x2 - x1) * 0.45);
        const d = `M ${x1} ${y1} C ${x1 + cp} ${y1}, ${x2 - cp} ${y2}, ${x2} ${y2}`;

        return {
          key: `${edge.from_node}.${edge.from_output}->${edge.to_node}.${edge.to_input}`,
          d,
        };
      })
      .filter((item): item is { key: string; d: string } => Boolean(item));
  }, [edges, canvasNodeById]);

  const previewLinkPath = useMemo(() => {
    if (!linkStartNodeId || !linkPreviewPoint) return "";
    const fromNode = canvasNodeById.get(linkStartNodeId);
    if (!fromNode) return "";
    const x1 = fromNode.x + NODE_WIDTH - 2;
    const y1 = fromNode.y + NODE_HEIGHT / 2;
    const x2 = linkPreviewPoint.x;
    const y2 = linkPreviewPoint.y;
    const cp = Math.max(48, Math.abs(x2 - x1) * 0.45);
    return `M ${x1} ${y1} C ${x1 + cp} ${y1}, ${x2 - cp} ${y2}, ${x2} ${y2}`;
  }, [canvasNodeById, linkStartNodeId, linkPreviewPoint]);

  const validationErrors = useMemo(() => {
    const issues: string[] = [];
    if (!projectId.trim()) issues.push("Project ID is required.");
    if (!workflowName.trim()) issues.push("Workflow name is required.");
    if (!workflowId) issues.push("Workflow ID cannot be empty.");
    if (nodes.length === 0) issues.push("Workflow must contain at least one node.");

    const nodeIds = nodes.map((n) => n.node_id);
    const seenNodeIds = new Set<string>();
    for (const id of nodeIds) {
      if (!id.trim()) continue;
      if (seenNodeIds.has(id)) issues.push(`Duplicate node_id: ${id}`);
      seenNodeIds.add(id);
    }
    issues.push(...parsedNodeTypes.errors);

    const edgeSeen = new Set<string>();
    for (const edge of edges) {
      const sig = `${edge.from_node}.${edge.from_output}->${edge.to_node}.${edge.to_input}`;
      if (edgeSeen.has(sig)) {
        issues.push(`Duplicate edge: ${sig}`);
      }
      edgeSeen.add(sig);
      if (edge.from_node === edge.to_node) {
        issues.push(`Self-loop is not allowed: ${sig}`);
      }
      const outType = parsedNodeTypes.outputs[edge.from_node]?.[edge.from_output];
      const inType = parsedNodeTypes.inputs[edge.to_node]?.[edge.to_input];
      if (!outType) issues.push(`Edge invalid source port: ${edge.from_node}.${edge.from_output}`);
      if (!inType) issues.push(`Edge invalid target port: ${edge.to_node}.${edge.to_input}`);
      if (outType && inType && outType !== inType) {
        issues.push(`Type mismatch: ${edge.from_node}.${edge.from_output} (${outType}) -> ${edge.to_node}.${edge.to_input} (${inType})`);
      }
    }

    try {
      const sweepObj = parseJsonObject(sweepsRaw);
      for (const [key, value] of Object.entries(sweepObj)) {
        if (!key.includes(".")) {
          issues.push(`Sweep key must be node.param: ${key}`);
          continue;
        }
        const [nodeId, paramName] = key.split(".", 2);
        if (!seenNodeIds.has(nodeId)) issues.push(`Sweep references unknown node: ${nodeId}`);
        if (!paramName) issues.push(`Sweep key missing param name: ${key}`);
        if (!Array.isArray(value) || value.length === 0) {
          issues.push(`Sweep value must be non-empty array: ${key}`);
        }
      }
    } catch (err) {
      issues.push(`Parameter sweeps JSON invalid: ${(err as Error).message}`);
    }

    return issues;
  }, [projectId, workflowName, workflowId, nodes, edges, sweepsRaw, parsedNodeTypes]);

  const authHeaders = token
    ? {
        Authorization: `Bearer ${token}`,
        "Content-Type": "application/json",
      }
    : null;

  useEffect(() => {
    if (typeof window === "undefined") return;
    const stored = window.localStorage.getItem("ad_api_token");
    if (stored) {
      setToken(stored);
    }
    setAuthReady(true);
  }, [router]);

  useEffect(() => {
    if (typeof window === "undefined" || !authReady) return;
    if (!token) {
      window.localStorage.removeItem("ad_api_token");
      return;
    }
    window.localStorage.setItem("ad_api_token", token);
  }, [token, authReady]);

  function logoutSession() {
    setToken(null);
    router.replace("/login");
  }

  function clampNodePosition(x: number, y: number): { x: number; y: number } {
    const nx = Math.max(NODE_MARGIN, Math.min(CANVAS_WIDTH - NODE_WIDTH - NODE_MARGIN, x));
    const ny = Math.max(NODE_MARGIN, Math.min(CANVAS_HEIGHT - NODE_HEIGHT - NODE_MARGIN, y));
    return { x: nx, y: ny };
  }

  function clientToCanvas(clientX: number, clientY: number): { x: number; y: number } | null {
    const canvasEl = canvasRef.current;
    if (!canvasEl) return null;
    const rect = canvasEl.getBoundingClientRect();
    return {
      x: (clientX - rect.left) / canvasZoom,
      y: (clientY - rect.top) / canvasZoom,
    };
  }

  async function loadPlugins() {
    if (!authHeaders) {
      setError("Please sign in first.");
      return;
    }
    setError("");
    const res = await fetch(`${apiBase}/plugins`, { headers: authHeaders });
    if (!res.ok) {
      setError(`Load plugins failed: ${res.status}`);
      return;
    }
    const body = await res.json();
    setPlugins(body.plugins || []);
    setStatus(`Loaded ${body.plugins?.length || 0} plugins.`);
  }

  function addNode(plugin: PaletteTool, preferredPosition?: { x: number; y: number }) {
    setError("");
    const usedIds = new Set(nodes.map((n) => n.node_id));
    let suffix = nodes.length + 1;
    let nodeId = `n${suffix}`;
    while (usedIds.has(nodeId)) {
      suffix += 1;
      nodeId = `n${suffix}`;
    }
    const fallbackX = NODE_MARGIN + ((nodes.length + 1) % 5) * 280;
    const fallbackY = NODE_MARGIN + Math.floor((nodes.length + 1) / 5) * 130;
    const pos = clampNodePosition(preferredPosition?.x ?? fallbackX, preferredPosition?.y ?? fallbackY);
    setNodes((prev) => [
      ...prev,
      {
        node_id: nodeId,
        plugin_id: plugin.plugin_id,
        version: plugin.version,
        input_types: JSON.stringify(plugin.default_input_types || {}, null, 2),
        output_types: JSON.stringify(plugin.default_output_types || { out: "report.html" }, null, 2),
        parameters: JSON.stringify(plugin.default_parameters || {}, null, 2),
      },
    ]);
    setNodePositions((prev) => ({ ...prev, [nodeId]: pos }));
    setStatus(`Added node ${nodeId} (${plugin.plugin_id})`);
  }

  function onDropPlugin(event: React.DragEvent<HTMLDivElement>) {
    event.preventDefault();
    const pluginId = event.dataTransfer.getData("text/plain");
    const plugin = paletteTools.find((p) => p.plugin_id === pluginId);
    if (!plugin) return;
    const point = clientToCanvas(event.clientX, event.clientY);
    addNode(plugin, point ? { x: point.x - NODE_WIDTH / 2, y: point.y - NODE_HEIGHT / 2 } : undefined);
  }

  function removeNode(nodeId: string) {
    setNodes((prev) => prev.filter((node) => node.node_id !== nodeId));
    setEdges((prev) => prev.filter((edge) => edge.from_node !== nodeId && edge.to_node !== nodeId));
    setNodePositions((prev) => {
      const next = { ...prev };
      delete next[nodeId];
      return next;
    });
    if (linkStartNodeId === nodeId) {
      setLinkStartNodeId(null);
      setLinkPreviewPoint(null);
    }
    setStatus(`Removed node ${nodeId}.`);
  }

  function paletteColorClass(tool: PaletteTool): string {
    if (tool.source === "registry") return "palette-registry";
    const tags = tool.tags || [];
    const hasWgs = tags.includes("wgs");
    const hasRna = tags.includes("rna-seq");
    if (hasWgs && hasRna) return "palette-mixed";
    if (hasWgs) return "palette-wgs";
    if (hasRna) return "palette-rna";
    return "palette-generic";
  }

  function resolveCompatiblePorts(fromNodeId: string, toNodeId: string): { from_output: string; to_input: string; type_id: string } | null {
    const fromOutputs = Object.entries(parsedNodeTypes.outputs[fromNodeId] || {});
    const toInputs = Object.entries(parsedNodeTypes.inputs[toNodeId] || {});
    for (const [fromOutput, outType] of fromOutputs) {
      for (const [toInput, inType] of toInputs) {
        if (outType === inType) {
          return { from_output: fromOutput, to_input: toInput, type_id: outType };
        }
      }
    }
    return null;
  }

  function edgeSignature(edge: EdgeDraft): string {
    return `${edge.from_node}.${edge.from_output}->${edge.to_node}.${edge.to_input}`;
  }

  function connectNodes(fromNodeId: string, toNodeId: string): boolean {
    if (!fromNodeId || !toNodeId) {
      setError("Choose from/to node first.");
      return false;
    }
    if (fromNodeId === toNodeId) {
      setError("Self-loop is not allowed.");
      return false;
    }
    const compatible = resolveCompatiblePorts(fromNodeId, toNodeId);
    if (!compatible) {
      setError(`No compatible data type between ${fromNodeId} and ${toNodeId}.`);
      return false;
    }
    const edgeDraft: EdgeDraft = {
      from_node: fromNodeId,
      from_output: compatible.from_output,
      to_node: toNodeId,
      to_input: compatible.to_input,
    };
    const sig = edgeSignature(edgeDraft);
    if (edges.some((e) => edgeSignature(e) === sig)) {
      setError(`Duplicate edge: ${sig}`);
      return false;
    }
    setEdges((prev) => [...prev, edgeDraft]);
    setStatus(`Edge added (${compatible.type_id}).`);
    return true;
  }

  function removeEdgeBySignature(sig: string) {
    setEdges((prev) => {
      const idx = prev.findIndex((edge) => edgeSignature(edge) === sig);
      if (idx < 0) return prev;
      const next = prev.slice();
      next.splice(idx, 1);
      return next;
    });
    setStatus(`Removed edge ${sig}.`);
  }

  function addEdgeAuto() {
    setError("");
    const { from_node, to_node } = simpleEdgeDraft;
    if (connectNodes(from_node, to_node)) {
      setSimpleEdgeDraft({ from_node: "", to_node: "" });
    }
  }

  function onNodeClick(nodeId: string) {
    setError("");
    if (!linkStartNodeId) {
      setLinkStartNodeId(nodeId);
      const fromNode = canvasNodeById.get(nodeId);
      if (fromNode) {
        setLinkPreviewPoint({ x: fromNode.x + NODE_WIDTH + 20, y: fromNode.y + NODE_HEIGHT / 2 });
      }
      setStatus(`Link mode: ${nodeId} -> ...`);
      return;
    }
    if (linkStartNodeId === nodeId) {
      setLinkStartNodeId(null);
      setLinkPreviewPoint(null);
      setStatus("Link mode canceled.");
      return;
    }
    const ok = connectNodes(linkStartNodeId, nodeId);
    if (ok) {
      setLinkStartNodeId(null);
      setLinkPreviewPoint(null);
    }
  }

  function onCanvasMouseMove(event: React.MouseEvent<HTMLDivElement>) {
    if (!linkStartNodeId || isDraggingNode) return;
    const point = clientToCanvas(event.clientX, event.clientY);
    if (!point) return;
    setLinkPreviewPoint(point);
  }

  function onNodeMouseDown(event: React.MouseEvent<HTMLDivElement>, nodeId: string) {
    if (event.button !== 0) return;
    if ((event.target as HTMLElement).closest("button")) return;
    const node = canvasNodeById.get(nodeId);
    if (!node) return;
    const point = clientToCanvas(event.clientX, event.clientY);
    if (!point) return;
    event.preventDefault();

    dragContextRef.current = {
      nodeId,
      startClientX: event.clientX,
      startClientY: event.clientY,
      offsetX: point.x - node.x,
      offsetY: point.y - node.y,
      moved: false,
    };

    const onMove = (moveEvent: MouseEvent) => {
      const ctx = dragContextRef.current;
      if (!ctx) return;
      const dx = moveEvent.clientX - ctx.startClientX;
      const dy = moveEvent.clientY - ctx.startClientY;
      if (!ctx.moved && Math.hypot(dx, dy) < 5) return;
      if (!ctx.moved) {
        ctx.moved = true;
        setIsDraggingNode(true);
      }
      const canvasEl = canvasRef.current;
      if (!canvasEl) return;
      const rect = canvasEl.getBoundingClientRect();
      const pointX = (moveEvent.clientX - rect.left) / canvasZoom;
      const pointY = (moveEvent.clientY - rect.top) / canvasZoom;
      const nextX = Math.max(NODE_MARGIN, Math.min(CANVAS_WIDTH - NODE_WIDTH - NODE_MARGIN, pointX - ctx.offsetX));
      const nextY = Math.max(NODE_MARGIN, Math.min(CANVAS_HEIGHT - NODE_HEIGHT - NODE_MARGIN, pointY - ctx.offsetY));
      setNodePositions((prev) => ({ ...prev, [ctx.nodeId]: { x: nextX, y: nextY } }));
    };

    const onUp = () => {
      const ctx = dragContextRef.current;
      dragContextRef.current = null;
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
      if (!ctx) return;
      if (ctx.moved) {
        setIsDraggingNode(false);
        return;
      }
      onNodeClick(nodeId);
    };

    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
  }

  function autoConnectSequence() {
    setError("");
    if (nodes.length < 2) {
      setError("Need at least 2 nodes to auto-connect sequence.");
      return;
    }
    const nextEdges = [...edges];
    const failures: string[] = [];
    let added = 0;
    for (let i = 0; i < nodes.length - 1; i += 1) {
      const fromNodeId = nodes[i].node_id;
      const toNodeId = nodes[i + 1].node_id;
      const compatible = resolveCompatiblePorts(fromNodeId, toNodeId);
      if (!compatible) {
        failures.push(`${fromNodeId} -> ${toNodeId}`);
        continue;
      }
      const candidate: EdgeDraft = {
        from_node: fromNodeId,
        from_output: compatible.from_output,
        to_node: toNodeId,
        to_input: compatible.to_input,
      };
      const sig = `${candidate.from_node}.${candidate.from_output}->${candidate.to_node}.${candidate.to_input}`;
      if (nextEdges.some((e) => `${e.from_node}.${e.from_output}->${e.to_node}.${e.to_input}` === sig)) {
        continue;
      }
      nextEdges.push(candidate);
      added += 1;
    }
    setEdges(nextEdges);
    if (failures.length > 0) {
      setError(`Some pairs have no compatible IO: ${failures.join(", ")}`);
    }
    setStatus(`Auto-connect finished. Added ${added} edge(s).`);
  }

  function buildWorkflowPayload() {
    if (validationErrors.length > 0) {
      throw new Error(`Please fix ${validationErrors.length} validation issue(s) before save/run.`);
    }
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
      setError("Please sign in first.");
      return;
    }
    setError("");
    let payload: unknown;
    try {
      payload = buildWorkflowPayload();
    } catch (err) {
      setError((err as Error).message);
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
        setError(`Validate failed: ${validate.status} ${text}`);
        return;
      }
      setStatus("Workflow already existed; validation refreshed and reusable.");
      return;
    }
    if (!create.ok) {
      const text = await create.text();
      setError(`Import failed: ${create.status} ${text}`);
      return;
    }
    setStatus("Workflow imported and validated.");
  }

  async function runDistributed() {
    if (!authHeaders) {
      setError("Please sign in first.");
      return;
    }
    setError("");
    if (validationErrors.length > 0) {
      setError(`Cannot run with validation issues (${validationErrors.length}).`);
      return;
    }
    const res = await fetch(`${apiBase}/workflows/${workflowId}/execute/distributed`, {
      method: "POST",
      headers: authHeaders,
      body: JSON.stringify({ max_workers: 2 }),
    });
    if (!res.ok) {
      const text = await res.text();
      setError(`Run failed: ${res.status} ${text}`);
      return;
    }
    const body = await res.json();
    setMonitor(body.summary);
    setRunDetails(body);
    setStatus(`Distributed run completed: ${body.summary.run_id}`);
  }

  async function refreshRun() {
    if (!authHeaders || !monitor?.run_id) return;
    setError("");
    const res = await fetch(`${apiBase}/workflows/runs/${monitor.run_id}`, { headers: authHeaders });
    if (!res.ok) {
      setError(`Refresh run failed: ${res.status}`);
      return;
    }
    const body = await res.json();
    setRunDetails(body);
    setStatus("Run details refreshed.");
  }

  async function exportReproReport() {
    if (!authHeaders) {
      setError("Please sign in first.");
      return;
    }
    setError("");
    const exportRes = await fetch(`${apiBase}/workflows/${workflowId}/export`, { headers: authHeaders });
    if (!exportRes.ok) {
      setError(`Export failed: ${exportRes.status}`);
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

  function saveDraftToBrowser() {
    if (!draftStorageKey) {
      setError("Project ID and workflow name are required.");
      return;
    }
    const payload = {
      projectId,
      workflowName,
      nodes,
      edges,
      nodePositions,
      sweepsRaw,
      savedAt: new Date().toISOString(),
    };
    localStorage.setItem(draftStorageKey, JSON.stringify(payload));
    setError("");
    setStatus(`Draft saved locally (${workflowId}).`);
  }

  function loadDraftFromBrowser() {
    if (!draftStorageKey) {
      setError("Project ID and workflow name are required.");
      return;
    }
    const raw = localStorage.getItem(draftStorageKey);
    if (!raw) {
      setError("No local draft found for this project/workflow.");
      return;
    }
    try {
      const payload = JSON.parse(raw) as {
        nodes: NodeDraft[];
        edges: EdgeDraft[];
        nodePositions?: Record<string, { x: number; y: number }>;
        sweepsRaw: string;
      };
      setNodes(payload.nodes || []);
      setEdges(payload.edges || []);
      setNodePositions(payload.nodePositions || {});
      setSweepsRaw(payload.sweepsRaw || "{}");
      setError("");
      setStatus(`Draft loaded (${workflowId}).`);
    } catch (err) {
      setError(`Draft parse failed: ${(err as Error).message}`);
    }
  }

  function clearCanvas() {
    setNodes([]);
    setEdges([]);
    setNodePositions({});
    setSweepsRaw("{}");
    setMonitor(null);
    setRunDetails(null);
    setLinkStartNodeId(null);
    setLinkPreviewPoint(null);
    setError("");
    setStatus("Canvas cleared.");
  }

  async function refreshProjectWorkflows() {
    if (!authHeaders) {
      setError("Please sign in first.");
      return;
    }
    setError("");
    const res = await fetch(`${apiBase}/workflows`, { headers: authHeaders });
    if (!res.ok) {
      setError(`Load workflows failed: ${res.status}`);
      return;
    }
    const body = (await res.json()) as { workflows?: WorkflowListItem[] };
    const prefix = `${projectId.trim()}__`;
    const rows = (body.workflows || []).filter((w) => w.workflow_id.startsWith(prefix));
    setSavedWorkflows(rows);
    if (rows.length > 0) {
      setSelectedSavedWorkflow(rows[0].workflow_id);
    } else {
      setSelectedSavedWorkflow("");
    }
    setStatus(`Loaded ${rows.length} saved workflow(s) for project ${projectId}.`);
  }

  async function loadSelectedWorkflowFromServer() {
    if (!authHeaders) {
      setError("Please sign in first.");
      return;
    }
    if (!selectedSavedWorkflow) {
      setError("Please select a saved workflow.");
      return;
    }
    setError("");
    const res = await fetch(`${apiBase}/workflows/${encodeURIComponent(selectedSavedWorkflow)}/export`, { headers: authHeaders });
    if (!res.ok) {
      const text = await res.text();
      setError(`Load workflow failed: ${res.status} ${text}`);
      return;
    }
    const body = await res.json();
    const wf = body.workflow || {};
    const loadedNodes = Array.isArray(wf.nodes) ? wf.nodes : [];
    const loadedEdges = Array.isArray(wf.edges) ? wf.edges : [];
    const loadedSweeps = wf.parameter_sweeps || {};

    const name = String(selectedSavedWorkflow).split("__").slice(1).join("__");
    setWorkflowName(name || workflowName);
    setNodes(
      loadedNodes.map((n: any) => ({
        node_id: String(n.node_id || ""),
        plugin_id: String(n.plugin_id || ""),
        version: String(n.version || ""),
        input_types: JSON.stringify(n.input_types || {}, null, 2),
        output_types: JSON.stringify(n.output_types || {}, null, 2),
        parameters: JSON.stringify(n.parameters || {}, null, 2),
      })),
    );
    setEdges(
      loadedEdges.map((e: any) => ({
        from_node: String(e.from_node || ""),
        from_output: String(e.from_output || ""),
        to_node: String(e.to_node || ""),
        to_input: String(e.to_input || ""),
      })),
    );
    setNodePositions({});
    setLinkStartNodeId(null);
    setLinkPreviewPoint(null);
    setSweepsRaw(JSON.stringify(loadedSweeps, null, 2));
    setStatus(`Loaded workflow ${selectedSavedWorkflow}.`);
  }

  if (!authReady) {
    return (
      <main className="workflow-page">
        <section className="surface-section">
          <h2>Checking session...</h2>
          <p>Preparing workflow builder.</p>
        </section>
      </main>
    );
  }

  return (
    <main className="analysis-page workflow-page">
      <div className="ambient-shape shape-a" />
      <div className="ambient-shape shape-b" />

      <header className="hero-panel reveal reveal-1">
        <div>
          <p className="hero-kicker">AD Multi-Omics Locus Evidence Platform</p>
          <h1 className="hero-title">Workflow Builder</h1>
          <p className="hero-note">
            Compose methods like Lego blocks, validate wiring before execution, and export reproducibility artifacts.
          </p>
        </div>
        <div className={`health-pill ${token ? "health-pill-ok" : "health-pill-warn"}`}>
          {token ? `Signed in - ${workflowId || "unscoped"}` : "Guest mode - sign in to run/save"}
        </div>
      </header>

      <section className="surface-section split-shell">
        <div className="mini-card">
          <h2>Session</h2>
          <p className={token ? "ok-text" : "err-text"} style={{ margin: 0 }}>
            {token ? "Authenticated session is active." : "No active token. Tool palette is still available in guest mode."}
          </p>
          <div className="btn-row">
            {token ? (
              <>
                <button onClick={logoutSession}>Logout</button>
                <button onClick={() => router.push("/login")}>Switch Account</button>
              </>
            ) : (
              <button onClick={() => router.push("/login")}>Sign In</button>
            )}
            <button onClick={() => void loadPlugins()}>Load Tools</button>
            <button onClick={() => void refreshProjectWorkflows()}>Load Saved Workflows</button>
          </div>
        </div>

        <div className="mini-card">
          <h2>Workflow</h2>
          <label>
            Project ID
            <input value={projectId} onChange={(e) => setProjectId(e.target.value)} />
          </label>
          <label>
            Workflow Name
            <input value={workflowName} onChange={(e) => setWorkflowName(e.target.value)} />
          </label>
          <label>
            Workflow ID (scoped)
            <input value={workflowId} readOnly className="readonly" />
          </label>
          <label>
            Parameter Sweeps (JSON)
            <textarea value={sweepsRaw} onChange={(e) => setSweepsRaw(e.target.value)} rows={4} />
          </label>

          <div className="btn-row">
            <button onClick={() => void saveAndValidateWorkflow()}>Validate & Save</button>
            <button onClick={() => void runDistributed()}>Run Distributed</button>
            <button onClick={() => void exportReproReport()}>Export Repro Report</button>
          </div>
          <div className="btn-row">
            <button onClick={saveDraftToBrowser}>Save Local Draft</button>
            <button onClick={loadDraftFromBrowser}>Load Local Draft</button>
            <button onClick={clearCanvas}>Clear Canvas</button>
          </div>
          <div className="btn-row">
            <select value={selectedSavedWorkflow} onChange={(e) => setSelectedSavedWorkflow(e.target.value)} className="wide-select">
              <option value="">Saved workflows in this project</option>
              {savedWorkflows.map((w) => (
                <option key={w.workflow_id} value={w.workflow_id}>
                  {w.workflow_id} ({w.node_count} nodes, {w.edge_count} edges)
                </option>
              ))}
            </select>
            <button onClick={() => void loadSelectedWorkflowFromServer()}>Load Selected</button>
          </div>
        </div>
      </section>

      <section className="surface-section">
        <h2>Validation</h2>
        {validationErrors.length === 0 ? (
          <p className="ok-text">No validation errors. Workflow is ready to save and run.</p>
        ) : (
          <div>
            <p className="err-text">Found {validationErrors.length} issue(s):</p>
            <ul>
              {validationErrors.map((issue, i) => (
                <li key={`issue-${i}`}>{issue}</li>
              ))}
            </ul>
          </div>
        )}
      </section>

      <section className="surface-section tool-shell">
        <aside className="mini-card">
          <h2>Tool Palette</h2>
          <p style={{ margin: 0, color: "#4d626e", fontSize: 12 }}>
            Built-in RNA-seq/WGS tools are always available. Registry tools are merged after clicking Load Tools.
          </p>
          <label>
            Tool Group
            <select value={paletteCategory} onChange={(e) => setPaletteCategory(e.target.value as "all" | "wgs" | "rna-seq")}>
              <option value="all">All Tools</option>
              <option value="wgs">WGS Tools</option>
              <option value="rna-seq">RNA-seq Tools</option>
            </select>
          </label>
          <label>
            Search Tool
            <input value={paletteQuery} onChange={(e) => setPaletteQuery(e.target.value)} placeholder="Search by id/name/tag" />
          </label>
          <p style={{ margin: 0, color: "#4d626e", fontSize: 12 }}>Showing {filteredPaletteTools.length} tool(s).</p>
          {filteredPaletteTools.map((p) => (
            <div
              key={`${p.plugin_id}:${p.version}`}
              draggable
              onDragStart={(e) => e.dataTransfer.setData("text/plain", p.plugin_id)}
              onDoubleClick={() => addNode(p)}
              className={`palette-item ${paletteColorClass(p)}`}
            >
              <strong>{p.plugin_id}</strong>
              <div>{p.name}</div>
              <div>v{p.version}</div>
              <div>{(p.tags || []).join(", ") || "no-tags"}</div>
              <div>{p.source === "registry" ? "source: API registry" : "source: built-in template"}</div>
            </div>
          ))}
          {filteredPaletteTools.length === 0 ? <p style={{ margin: 0 }}>No tools matched this filter.</p> : null}
        </aside>

        <div onDragOver={(e) => e.preventDefault()} onDrop={onDropPlugin} className="canvas-box">
          <div className="canvas-head">
            <div>
              <h2>Canvas</h2>
              <p style={{ margin: 0, fontSize: 12 }}>
                Drag blocks to any position. One block can branch to multiple downstream blocks.
              </p>
            </div>
            <div className="canvas-controls">
              <button type="button" onClick={() => setCanvasZoom((z) => Math.max(0.6, Math.round((z - 0.1) * 100) / 100))}>
                -
              </button>
              <span>{Math.round(canvasZoom * 100)}%</span>
              <button type="button" onClick={() => setCanvasZoom((z) => Math.min(1.8, Math.round((z + 0.1) * 100) / 100))}>
                +
              </button>
              <button
                type="button"
                onClick={() => {
                  setLinkStartNodeId(null);
                  setLinkPreviewPoint(null);
                }}
                disabled={!linkStartNodeId}
                className="link-cancel-btn"
              >
                Cancel Link
              </button>
            </div>
          </div>
          {nodes.length === 0 ? <p>No nodes yet. Drag a tool from the left panel.</p> : null}
          {nodes.length > 0 ? (
            <div className="mindmap-shell">
              <div className="mindmap-viewport">
                <div
                  className="mindmap-zoom-layer"
                  style={{
                    width: `${Math.ceil(CANVAS_WIDTH * canvasZoom)}px`,
                    height: `${Math.ceil(CANVAS_HEIGHT * canvasZoom)}px`,
                  }}
                >
                  <div
                    ref={canvasRef}
                    className={`canvas-stage ${isDraggingNode ? "dragging" : ""}`}
                    style={{
                      width: `${CANVAS_WIDTH}px`,
                      height: `${CANVAS_HEIGHT}px`,
                      transform: `scale(${canvasZoom})`,
                    }}
                    onMouseMove={onCanvasMouseMove}
                    onMouseLeave={() => setLinkPreviewPoint(null)}
                  >
                    <svg className="canvas-links" width={CANVAS_WIDTH} height={CANVAS_HEIGHT} viewBox={`0 0 ${CANVAS_WIDTH} ${CANVAS_HEIGHT}`}>
                      {canvasEdgePaths.map((edge) => (
                        <path key={edge.key} d={edge.d} className="mind-link" onDoubleClick={() => removeEdgeBySignature(edge.key)} />
                      ))}
                      {previewLinkPath ? <path d={previewLinkPath} className="mind-link preview-link" /> : null}
                    </svg>
                    {canvasNodes.map((node) => (
                      <div
                        key={node.nodeId}
                        className={`canvas-node ${node.tone} ${linkStartNodeId === node.nodeId ? "selected-source" : ""}`}
                        style={{
                          left: `${node.x}px`,
                          top: `${node.y}px`,
                        }}
                        onMouseDown={(e) => onNodeMouseDown(e, node.nodeId)}
                      >
                        <div className="node-title-row">
                          <strong>{node.nodeId}</strong>
                          <span>{node.title}</span>
                        </div>
                        <p className="node-plugin">{node.pluginId}</p>
                        <div className="node-actions">
                          <button type="button" onClick={() => removeNode(node.nodeId)} className="small-btn danger-btn">
                            Remove
                          </button>
                        </div>
                      </div>
                    ))}
                  </div>
                </div>
              </div>
              <p className="canvas-tip">Linking: click source block, then click target block. Double-click a line to remove that connection.</p>
            </div>
          ) : null}
        </div>
      </section>

      <section className="surface-section">
        <h2>Connections</h2>
        <div className="edge-grid">
          <select value={simpleEdgeDraft.from_node} onChange={(e) => setSimpleEdgeDraft({ ...simpleEdgeDraft, from_node: e.target.value })}>
            <option value="">From node</option>
            {nodes.map((n) => (
              <option key={n.node_id} value={n.node_id}>
                {n.node_id}
              </option>
            ))}
          </select>
          <select value={simpleEdgeDraft.to_node} onChange={(e) => setSimpleEdgeDraft({ ...simpleEdgeDraft, to_node: e.target.value })}>
            <option value="">To node</option>
            {nodes.map((n) => (
              <option key={n.node_id} value={n.node_id}>
                {n.node_id}
              </option>
            ))}
          </select>
        </div>
        <div className="btn-row">
          <button onClick={addEdgeAuto}>Auto Connect Pair</button>
          <button onClick={autoConnectSequence} disabled={nodes.length < 2}>
            Auto Connect Sequence
          </button>
          <button
            onClick={() => {
              setEdges((prev) => prev.slice(0, Math.max(0, prev.length - 1)));
            }}
            disabled={edges.length === 0}
          >
            Remove Last Edge
          </button>
        </div>
        {edges.length === 0 ? (
          <p style={{ margin: 0 }}>No connections yet.</p>
        ) : (
          <ul style={{ margin: 0 }}>
            {edges.map((edge, idx) => {
              const typeId = parsedNodeTypes.outputs[edge.from_node]?.[edge.from_output] || "unknown";
              return (
                <li key={`${edge.from_node}-${edge.from_output}-${edge.to_node}-${edge.to_input}-${idx}`}>
                  {edge.from_node}
                  {" -> "}
                  {edge.to_node} ({typeId})
                </li>
              );
            })}
          </ul>
        )}
      </section>

      <section className="surface-section">
        <h2>Run Monitor</h2>
        <div className="btn-row">
          <button onClick={() => void refreshRun()} disabled={!monitor?.run_id}>
            Refresh Last Run
          </button>
        </div>
        {monitor ? (
          <div>
            <div>Run ID: {monitor.run_id}</div>
            <div>
              Completed: {monitor.completed_runs}/{monitor.submitted_runs}
            </div>
            <div>Duration: {monitor.duration_ms} ms</div>
          </div>
        ) : (
          <p>No run yet.</p>
        )}
        <pre className="code-box">{JSON.stringify(runDetails || {}, null, 2)}</pre>
      </section>

      {error ? <p className="err-text">{error}</p> : null}
      {status ? <p className="ok-text">{status}</p> : null}

      <style jsx>{`
        .workflow-page {
          --bg: #edf2f4;
          --surface: #ffffff;
          --line: #ccd6dc;
          --ink: #10222c;
          --muted: #4d626e;
          --teal: #0f7d76;
          --orange: #c8611f;
          --deep: #1f3441;
          min-height: 100vh;
          padding: 24px clamp(12px, 2vw, 28px) 40px;
          display: grid;
          grid-template-columns: repeat(12, minmax(0, 1fr));
          gap: 14px;
          background: repeating-linear-gradient(
              90deg,
              rgba(255, 255, 255, 0.26),
              rgba(255, 255, 255, 0.26) 1px,
              transparent 1px,
              transparent 26px
            ),
            repeating-linear-gradient(
              0deg,
              rgba(255, 255, 255, 0.26),
              rgba(255, 255, 255, 0.26) 1px,
              transparent 1px,
              transparent 26px
            ),
            var(--bg);
          color: var(--ink);
          font-family: "IBM Plex Sans", "Noto Sans TC", "Segoe UI", sans-serif;
          position: relative;
          overflow: hidden;
        }
        .ambient-shape {
          position: fixed;
          width: 280px;
          height: 280px;
          border-radius: 999px;
          pointer-events: none;
          z-index: 0;
          opacity: 0.18;
          filter: blur(8px);
        }
        .shape-a {
          top: -90px;
          left: -70px;
          background: #0f7d76;
        }
        .shape-b {
          top: -65px;
          right: -90px;
          background: #c8611f;
        }
        .hero-panel {
          z-index: 1;
          grid-column: 1 / -1;
          border: 1px solid var(--line);
          border-radius: 18px;
          padding: 18px;
          background: var(--surface);
          display: flex;
          justify-content: space-between;
          gap: 14px;
          align-items: flex-start;
        }
        .hero-kicker {
          margin: 0 0 4px;
          letter-spacing: 0.14em;
          text-transform: uppercase;
          font-size: 11px;
          font-weight: 700;
          color: var(--teal);
        }
        .hero-title {
          margin: 0;
          font-family: "Space Grotesk", "IBM Plex Sans", sans-serif;
          font-size: clamp(30px, 4.6vw, 48px);
          line-height: 1.03;
          letter-spacing: 0.01em;
        }
        .hero-note {
          margin: 6px 0 0;
          max-width: 720px;
          color: var(--muted);
          font-size: 13px;
        }
        .health-pill {
          border: 1px solid var(--line);
          border-radius: 999px;
          padding: 8px 14px;
          font-size: 12px;
          font-weight: 700;
          white-space: nowrap;
          background: #f7fafb;
        }
        .health-pill-ok {
          border-color: #78b38f;
          color: #135a2e;
          background: #e9f8ee;
        }
        .health-pill-warn {
          border-color: #d08f7a;
          color: #972315;
          background: #fff0ea;
        }
        .surface-section {
          z-index: 1;
          grid-column: 1 / -1;
          border: 1px solid var(--line);
          border-radius: 16px;
          background: var(--surface);
          padding: 14px;
          display: grid;
          gap: 10px;
          animation: rise 420ms ease both;
        }
        .split-shell {
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 12px;
        }
        .mini-card {
          border: 1px solid var(--line);
          border-radius: 12px;
          padding: 12px;
          display: grid;
          gap: 8px;
          background: #fbfdfe;
        }
        .tool-shell {
          grid-template-columns: 300px minmax(0, 1fr);
          gap: 12px;
        }
        .canvas-box {
          border: 1px solid #232a33;
          border-radius: 16px;
          padding: 12px;
          min-height: 420px;
          background: radial-gradient(circle at 12% 10%, #3a424f 0, #2e343f 32%, #292d36 100%);
          color: #dce5f0;
          display: grid;
          gap: 10px;
        }
        .canvas-box p {
          color: #d2dde9;
        }
        .canvas-head {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 8px;
        }
        .canvas-head h2 {
          color: #eef4ff;
        }
        .canvas-controls {
          display: inline-flex;
          align-items: center;
          gap: 6px;
          border: 1px solid #596677;
          border-radius: 999px;
          padding: 4px 6px;
          background: rgba(17, 24, 32, 0.56);
          color: #d6dfeb;
          font-size: 12px;
          font-weight: 700;
        }
        .canvas-controls button {
          min-height: 30px;
          width: 30px;
          padding: 0;
          border-radius: 999px;
          border: 1px solid #6b7c90;
          background: #1f2732;
          color: #ecf2fb;
          text-transform: none;
          font-size: 18px;
          line-height: 1;
        }
        .canvas-controls button:hover {
          border-color: #95aac0;
          background: #263242;
        }
        .canvas-controls .link-cancel-btn {
          width: auto;
          min-height: 30px;
          padding: 0 10px;
          border-radius: 999px;
          font-size: 11px;
          letter-spacing: 0.02em;
        }
        .mindmap-shell {
          display: grid;
          gap: 8px;
        }
        .mindmap-viewport {
          border: 1px solid #3e4856;
          border-radius: 12px;
          overflow: auto;
          min-height: 336px;
          background: rgba(10, 15, 22, 0.26);
        }
        .mindmap-zoom-layer {
          position: relative;
        }
        .canvas-stage {
          position: relative;
          transform-origin: 0 0;
          border: 1px dashed rgba(255, 255, 255, 0.18);
          border-radius: 12px;
          background: linear-gradient(rgba(255, 255, 255, 0.02) 1px, transparent 1px),
            linear-gradient(90deg, rgba(255, 255, 255, 0.02) 1px, transparent 1px), rgba(12, 16, 24, 0.34);
          background-size: 28px 28px;
        }
        .canvas-stage.dragging {
          cursor: grabbing;
        }
        .canvas-links {
          position: absolute;
          inset: 0;
          pointer-events: auto;
        }
        .mind-link {
          fill: none;
          stroke: #94adc8;
          stroke-width: 2;
          stroke-linecap: round;
          opacity: 0.92;
          pointer-events: stroke;
          cursor: pointer;
        }
        .preview-link {
          stroke: #f5cf7f;
          stroke-dasharray: 5 5;
          pointer-events: none;
        }
        .canvas-node {
          position: absolute;
          width: ${NODE_WIDTH}px;
          min-height: ${NODE_HEIGHT}px;
          border-radius: 12px;
          border: 1px solid #8ea3bb;
          padding: 8px 9px;
          color: #edf4fe;
          user-select: none;
          cursor: grab;
          backdrop-filter: blur(2px);
        }
        .canvas-node.selected-source {
          border-color: #ffde9c;
          box-shadow: 0 0 0 2px rgba(255, 222, 156, 0.26);
        }
        .canvas-node:active {
          cursor: grabbing;
        }
        .canvas-node.node-wgs {
          background: #554020;
          border-color: #dbb776;
        }
        .canvas-node.node-rna {
          background: #174c4a;
          border-color: #7fd4c8;
        }
        .canvas-node.node-mixed {
          background: #2a415c;
          border-color: #90b3d8;
        }
        .canvas-node.node-registry {
          background: #49325f;
          border-color: #c4a1ef;
        }
        .canvas-node.node-generic {
          background: #344558;
          border-color: #93a7be;
        }
        .node-title-row {
          display: flex;
          align-items: baseline;
          gap: 6px;
          white-space: nowrap;
          overflow: hidden;
          text-overflow: ellipsis;
        }
        .node-title-row strong {
          font-size: 12px;
          color: #f5f8ff;
        }
        .node-title-row span {
          font-size: 11px;
          color: #d7e1ee;
          overflow: hidden;
          text-overflow: ellipsis;
        }
        .node-plugin {
          margin: 4px 0 8px;
          font-size: 10px;
          color: #c8d4e3 !important;
          letter-spacing: 0.02em;
          white-space: nowrap;
          overflow: hidden;
          text-overflow: ellipsis;
        }
        .node-actions {
          display: flex;
          gap: 6px;
        }
        .small-btn {
          min-height: 24px;
          border-radius: 999px;
          padding: 0 8px;
          font-size: 10px;
          text-transform: none;
          letter-spacing: 0.02em;
          background: rgba(17, 27, 39, 0.9) !important;
          border: 1px solid rgba(191, 207, 224, 0.58);
        }
        .small-btn:hover {
          border-color: #d4e3f2;
        }
        .danger-btn {
          background: rgba(95, 18, 24, 0.95) !important;
          border-color: rgba(241, 159, 159, 0.7);
        }
        .canvas-tip {
          margin: 0;
          color: #d2dde9;
          font-size: 12px;
        }
        .palette-item {
          border: 1px solid #8c959f;
          border-radius: 6px;
          padding: 8px;
          margin-bottom: 8px;
          cursor: grab;
          background: #f6f8fa;
          font-size: 12px;
          color: var(--muted);
        }
        .palette-wgs {
          border-color: #9a6a13;
          background: #fff7e8;
        }
        .palette-rna {
          border-color: #0f7d76;
          background: #e9f8f6;
        }
        .palette-mixed {
          border-color: #1f3441;
          background: #ecf2f7;
        }
        .palette-registry {
          border-color: #7a3fc1;
          background: #f3ecff;
        }
        .palette-generic {
          border-color: #8c959f;
          background: #f6f8fa;
        }
        .palette-item strong {
          font-size: 13px;
          color: var(--ink);
        }
        h2 {
          margin: 0;
          font-family: "Space Grotesk", "IBM Plex Sans", sans-serif;
          font-size: 21px;
          letter-spacing: 0.01em;
        }
        p,
        li {
          color: var(--muted);
        }
        label {
          display: grid;
          gap: 6px;
          color: var(--muted);
          font-size: 12px;
          font-weight: 700;
          letter-spacing: 0.02em;
        }
        input,
        select,
        textarea {
          min-height: 42px;
          border-radius: 10px;
          border: 1px solid #c8d2d8;
          background: #ffffff;
          color: var(--ink);
          padding: 8px 10px;
          font-size: 14px;
          font-family: "IBM Plex Sans", "Noto Sans TC", "Segoe UI", sans-serif;
        }
        textarea {
          min-height: 68px;
          resize: vertical;
        }
        input:focus,
        select:focus,
        textarea:focus {
          outline: 2px solid #0f7d76;
          outline-offset: 1px;
          border-color: #0f7d76;
        }
        .readonly {
          background: #f6f8fa;
        }
        button {
          min-height: 42px;
          border-radius: 10px;
          border: 1px solid transparent;
          padding: 8px 12px;
          cursor: pointer;
          font-size: 13px;
          font-weight: 700;
          letter-spacing: 0.03em;
          text-transform: uppercase;
          background: var(--teal);
          color: #ffffff;
          transition: border-color 150ms ease, background-color 150ms ease, color 150ms ease;
        }
        button:nth-of-type(even) {
          background: var(--deep);
        }
        button:hover {
          border-color: #0c645f;
        }
        button:disabled {
          opacity: 0.55;
          cursor: default;
        }
        .btn-row {
          display: flex;
          flex-wrap: wrap;
          gap: 8px;
        }
        .wide-select {
          min-width: 280px;
          flex: 1;
        }
        .edge-grid {
          display: grid;
          grid-template-columns: repeat(2, minmax(180px, 1fr));
          gap: 8px;
        }
        .code-box {
          background: #f6f8fa;
          border: 1px solid #d0d7de;
          border-radius: 10px;
          padding: 10px;
          margin: 0;
          overflow-x: auto;
          color: #243746;
          font-size: 12px;
        }
        .ok-text {
          margin: 0;
          color: #14532d;
          font-weight: 700;
        }
        .err-text {
          margin: 0;
          color: #b42318;
          font-weight: 700;
        }
        *:hover {
          box-shadow: none !important;
        }
        .reveal {
          animation: rise 380ms ease both;
        }
        .reveal-1 {
          animation-delay: 20ms;
        }
        @keyframes rise {
          from {
            opacity: 0;
            transform: translateY(8px);
          }
          to {
            opacity: 1;
            transform: translateY(0);
          }
        }
        @media (max-width: 980px) {
          .hero-panel {
            flex-direction: column;
          }
          .split-shell,
          .tool-shell {
            grid-template-columns: 1fr;
          }
          .edge-grid {
            grid-template-columns: repeat(2, minmax(0, 1fr));
          }
          .wide-select {
            min-width: 0;
            width: 100%;
          }
        }
      `}</style>
    </main>
  );
}
