"use client";

import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useRouter } from "next/navigation";
import { Activity, Database, FlaskConical, Play, Workflow } from "lucide-react";

import { AppShell } from "@/components/app-shell";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Select } from "@/components/ui/select";
import { resolveApiBase } from "@/lib/api-base";
import { clearStoredToken, readStoredToken, TOKEN_STORAGE_KEY } from "@/lib/session";

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

type RunMonitor = {
  run_id: string;
  submitted_runs: number;
  completed_runs: number;
  failed_runs?: number;
  duration_ms: number;
  local_record_path?: string;
};

type WorkflowListItem = {
  workflow_id: string;
  node_count: number;
  edge_count: number;
  created_by: string;
};

type DatasetItem = {
  id: number;
  filename: string;
  project_id: string | null;
  input_role: string;
  canonical_type: string | null;
  validation_status: string;
  validation_detail: string;
};

type RunDetailResult = {
  index?: number;
  status?: string;
  engine?: string;
  workflow_artifacts?: {
    run_dir?: string;
    snakefile?: string;
    configfile?: string;
  };
  outputs?: string[];
  dry_run?: { status?: string; returncode?: number };
  run?: { status?: string; returncode?: number };
};

type RunDetailsResponse = {
  workflow_id?: string;
  created_by?: string;
  local_record_path?: string;
  summary?: RunMonitor & { local_record_path?: string };
  results?: RunDetailResult[];
};

type TemplateNodeSeed = {
  node_id: string;
  plugin_id: string;
  x: number;
  y: number;
};

const BUILTIN_TOOLS: PaletteTool[] = [
  {
    plugin_id: "fastqc",
    name: "FastQC",
    version: "0.1.0",
    tags: ["wgs", "qc", "locked"],
    source: "builtin",
    default_input_types: { reads: "reads.fastq.gz" },
    default_output_types: { fastqc_html: "report.html" },
    default_parameters: { threads: 2, stage: "raw" },
  },
  {
    plugin_id: "cutadapt",
    name: "Cutadapt",
    version: "0.1.0",
    tags: ["wgs", "trim", "locked"],
    source: "builtin",
    default_input_types: { reads: "reads.fastq.gz" },
    default_output_types: { trimmed_reads: "reads.fastq.gz" },
    default_parameters: { quality_cutoff: 20, min_length: 20, threads: 4 },
  },
  {
    plugin_id: "bwa_mem",
    name: "BWA-MEM",
    version: "0.1.0",
    tags: ["wgs", "align", "locked"],
    source: "builtin",
    default_input_types: { reads: "reads.fastq.gz", reference: "reference.genome.fasta" },
    default_output_types: { alignment: "align.bam" },
    default_parameters: { threads: 8, preset: "illumina-pe" },
  },
  {
    plugin_id: "samtools_sort_index",
    name: "samtools sort+index",
    version: "0.1.0",
    tags: ["wgs", "align", "locked"],
    source: "builtin",
    default_input_types: { alignment: "align.bam" },
    default_output_types: { alignment_sorted: "align.bam" },
    default_parameters: { threads: 4, create_index: true },
  },
  {
    plugin_id: "gatk_haplotypecaller",
    name: "GATK HaplotypeCaller",
    version: "0.1.0",
    tags: ["wgs", "variant", "locked"],
    source: "builtin",
    default_input_types: { alignment_sorted: "align.bam", reference: "reference.genome.fasta" },
    default_output_types: { gvcf: "variants.vcf.gz" },
    default_parameters: { emit_mode: "gvcf", intervals: "whole-genome" },
  },
  {
    plugin_id: "gatk_genotypegvcfs",
    name: "GATK GenotypeGVCFs",
    version: "0.1.0",
    tags: ["wgs", "variant", "locked"],
    source: "builtin",
    default_input_types: { gvcf: "variants.vcf.gz", reference: "reference.genome.fasta" },
    default_output_types: { variants: "variants.vcf.gz" },
    default_parameters: { call_conf: 30, output_mode: "sites-only" },
  },
  {
    plugin_id: "variant_report",
    name: "Variant Report",
    version: "0.1.0",
    tags: ["wgs", "report", "locked"],
    source: "builtin",
    default_input_types: { variants: "variants.vcf.gz" },
    default_output_types: { report: "report.html" },
    default_parameters: { include_qc_summary: true },
  },
];

const WGS_TEMPLATE_NODES: TemplateNodeSeed[] = [
  { node_id: "fastqc_raw", plugin_id: "fastqc", x: 50, y: 80 },
  { node_id: "cutadapt", plugin_id: "cutadapt", x: 340, y: 80 },
  { node_id: "bwa_mem", plugin_id: "bwa_mem", x: 630, y: 80 },
  { node_id: "samtools_sort", plugin_id: "samtools_sort_index", x: 920, y: 80 },
  { node_id: "gatk_hc", plugin_id: "gatk_haplotypecaller", x: 1210, y: 80 },
  { node_id: "gatk_genotype", plugin_id: "gatk_genotypegvcfs", x: 1210, y: 250 },
  { node_id: "variant_report", plugin_id: "variant_report", x: 1500, y: 250 },
];

const WGS_TEMPLATE_EDGES: EdgeDraft[] = [
  { from_node: "cutadapt", from_output: "trimmed_reads", to_node: "bwa_mem", to_input: "reads" },
  { from_node: "bwa_mem", from_output: "alignment", to_node: "samtools_sort", to_input: "alignment" },
  { from_node: "samtools_sort", from_output: "alignment_sorted", to_node: "gatk_hc", to_input: "alignment_sorted" },
  { from_node: "gatk_hc", from_output: "gvcf", to_node: "gatk_genotype", to_input: "gvcf" },
  { from_node: "gatk_genotype", from_output: "variants", to_node: "variant_report", to_input: "variants" },
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
  const apiBase = useMemo(() => resolveApiBase(), []);

  const [token, setToken] = useState<string | null>(null);
  const [authReady, setAuthReady] = useState(false);
  const [sessionChecked, setSessionChecked] = useState(false);
  const [plugins] = useState<PluginItem[]>([]);
  const [projectId, setProjectId] = useState("demo-project");
  const [workflowName, setWorkflowName] = useState("lego-flow");
  const [nodes, setNodes] = useState<NodeDraft[]>([]);
  const [edges, setEdges] = useState<EdgeDraft[]>([]);
  const [sweepsRaw, setSweepsRaw] = useState("{}");
  const [status, setStatus] = useState<string>("");
  const [error, setError] = useState<string>("");
  const [monitor, setMonitor] = useState<RunMonitor | null>(null);
  const [runDetails, setRunDetails] = useState<unknown>(null);
  const [rawDataFiles, setRawDataFiles] = useState<File[]>([]);
  const [uploadLoading, setUploadLoading] = useState(false);
  const [uploadError, setUploadError] = useState<string>("");
  const [uploadMessage, setUploadMessage] = useState<string>("");
  const [savedWorkflows, setSavedWorkflows] = useState<WorkflowListItem[]>([]);
  const [projectDatasets, setProjectDatasets] = useState<DatasetItem[]>([]);
  const [selectedSavedWorkflow, setSelectedSavedWorkflow] = useState("");
  const [canvasZoom, setCanvasZoom] = useState(1);
  const [showRawResult, setShowRawResult] = useState(false);
  const [nodePositions, setNodePositions] = useState<Record<string, { x: number; y: number }>>({});
  const [isDraggingNode, setIsDraggingNode] = useState(false);
  const canvasRef = useRef<HTMLDivElement | null>(null);
  const rawUploadInputRef = useRef<HTMLInputElement | null>(null);
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

  const filteredPaletteTools = useMemo(() => paletteTools.filter((tool) => (tool.tags || []).includes("wgs")), [paletteTools]);

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

  const validationErrors = useMemo(() => {
    const issues: string[] = [];
    const allowedToolIds = new Set(BUILTIN_TOOLS.map((tool) => tool.plugin_id));
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
    for (const node of nodes) {
      if (!allowedToolIds.has(node.plugin_id)) {
        issues.push(`Node ${node.node_id} uses unsupported tool for locked WGS mode: ${node.plugin_id}`);
      }
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

    if (sweepsRaw.trim()) {
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
    }

    return issues;
  }, [projectId, workflowName, workflowId, nodes, edges, sweepsRaw, parsedNodeTypes]);

  const handleUnauthorized = useCallback((message = "Session expired. Please sign in again.") => {
    clearStoredToken();
    setToken(null);
    setUploadError(message);
    setError(message);
    router.replace("/login");
  }, [router]);

  const authHeaders = useMemo(() => (token
    ? {
        Authorization: `Bearer ${token}`,
        "Content-Type": "application/json",
      }
    : null), [token]);
  const uploadHeaders = useMemo(() => (token
    ? {
        Authorization: `Bearer ${token}`,
      }
    : null), [token]);

  const primaryRunResult = useMemo(() => {
    const body = runDetails as RunDetailsResponse | null;
    if (!body?.results || body.results.length === 0) return null;
    return body.results[0];
  }, [runDetails]);

  const runRecordPath = useMemo(() => {
    const body = runDetails as RunDetailsResponse | null;
    return body?.local_record_path || body?.summary?.local_record_path || monitor?.local_record_path || "";
  }, [monitor, runDetails]);

  function fileKey(file: File): string {
    return `${file.name}:${file.size}:${file.lastModified}`;
  }

  function appendRawFiles(fileList: FileList | null): void {
    const nextFiles = Array.from(fileList || []);
    if (nextFiles.length === 0) return;
    setRawDataFiles((prev) => {
      const merged = [...prev];
      const seen = new Set(prev.map(fileKey));
      for (const file of nextFiles) {
        const key = fileKey(file);
        if (seen.has(key)) continue;
        merged.push(file);
        seen.add(key);
      }
      return merged;
    });
    if (rawUploadInputRef.current) {
      rawUploadInputRef.current.value = "";
    }
    setUploadError("");
    setUploadMessage("");
  }

  function removeRawFile(targetKey: string): void {
    setRawDataFiles((prev) => prev.filter((file) => fileKey(file) !== targetKey));
  }

  function clearRawFiles(): void {
    setRawDataFiles([]);
    if (rawUploadInputRef.current) {
      rawUploadInputRef.current.value = "";
    }
  }

  useEffect(() => {
    if (typeof window === "undefined") return;
    const stored = readStoredToken();
    if (stored) {
      setToken(stored);
    }
    setAuthReady(true);
  }, [router]);

  useEffect(() => {
    if (typeof window === "undefined" || !authReady) return;
    if (!token) {
      clearStoredToken();
      setSessionChecked(true);
      return;
    }
    window.localStorage.setItem(TOKEN_STORAGE_KEY, token);
  }, [token, authReady]);

  useEffect(() => {
    if (!authReady || !token) return;
    let cancelled = false;
    setSessionChecked(false);
    void fetch(`${apiBase}/auth/me`, {
      headers: { Authorization: `Bearer ${token}` },
    }).then(async (response) => {
      if (cancelled) {
        return;
      }
      if (response.ok) {
        setSessionChecked(true);
        return;
      }
      if (response.status === 401) {
        clearStoredToken();
        setToken(null);
        setError("Session expired. Please sign in again.");
        setSessionChecked(true);
        router.replace("/login");
        return;
      }
      setSessionChecked(true);
    }).catch(() => {
      if (!cancelled) {
        setSessionChecked(true);
      }
    });
    return () => {
      cancelled = true;
    };
  }, [apiBase, authReady, router, token]);

  useEffect(() => {
    if (!sessionChecked || !authHeaders) return;
    let cancelled = false;
    void fetch(`${apiBase}/datasets`, { headers: authHeaders })
      .then((response) => response.json().then((body) => ({ ok: response.ok, status: response.status, body })))
      .then(({ ok, status, body }) => {
        if (cancelled) return;
        if (!ok) {
          if (status === 401) {
            handleUnauthorized("Invalid token");
          }
          return;
        }
        const rows = Array.isArray(body?.datasets) ? body.datasets : [];
        setProjectDatasets(rows.filter((row: DatasetItem) => row.project_id === projectId));
      })
      .catch(() => {
        if (!cancelled) {
          setProjectDatasets([]);
        }
      });
    return () => {
      cancelled = true;
    };
  }, [apiBase, authHeaders, handleUnauthorized, projectId, sessionChecked]);

  function logoutSession() {
    setToken(null);
    clearStoredToken();
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

  function loadWgsTemplate() {
    const toolTable = new Map(BUILTIN_TOOLS.map((tool) => [tool.plugin_id, tool] as const));
    const seededNodes: NodeDraft[] = [];
    const seededPositions: Record<string, { x: number; y: number }> = {};

    for (const seed of WGS_TEMPLATE_NODES) {
      const plugin = toolTable.get(seed.plugin_id);
      if (!plugin) {
        setError(`Missing built-in tool definition: ${seed.plugin_id}`);
        return;
      }
      seededNodes.push({
        node_id: seed.node_id,
        plugin_id: plugin.plugin_id,
        version: plugin.version,
        input_types: JSON.stringify(plugin.default_input_types || {}, null, 2),
        output_types: JSON.stringify(plugin.default_output_types || {}, null, 2),
        parameters: JSON.stringify(plugin.default_parameters || {}, null, 2),
      });
      seededPositions[seed.node_id] = clampNodePosition(seed.x, seed.y);
    }

    setNodes(seededNodes);
    setEdges(WGS_TEMPLATE_EDGES);
    setNodePositions(seededPositions);
    setSweepsRaw("{}");
    setError("");
    setStatus("Loaded WGS fixed-node template.");
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
    const fallbackFromOutput = Object.keys(parsedNodeTypes.outputs[fromNodeId] || {})[0];
    const fallbackToInput = Object.keys(parsedNodeTypes.inputs[toNodeId] || {})[0];
    let edgeDraft: EdgeDraft;
    let warnIncompatible = false;
    if (compatible) {
      edgeDraft = {
        from_node: fromNodeId,
        from_output: compatible.from_output,
        to_node: toNodeId,
        to_input: compatible.to_input,
      };
    } else {
      if (!fallbackFromOutput || !fallbackToInput) {
        setError(`Cannot connect ${fromNodeId} -> ${toNodeId}: missing output/input ports.`);
        return false;
      }
      edgeDraft = {
        from_node: fromNodeId,
        from_output: fallbackFromOutput,
        to_node: toNodeId,
        to_input: fallbackToInput,
      };
      warnIncompatible = true;
    }
    const sig = edgeSignature(edgeDraft);
    if (edges.some((e) => edgeSignature(e) === sig)) {
      setError(`Duplicate edge: ${sig}`);
      return false;
    }
    setEdges((prev) => [...prev, edgeDraft]);
    if (warnIncompatible) {
      setStatus(`Edge added as draft (${fromNodeId} -> ${toNodeId}). Type mismatch will be shown in Validation.`);
    } else {
      setStatus(`Edge added (${compatible.type_id}).`);
    }
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
        setStatus(`Moved ${nodeId}.`);
        return;
      }
      setStatus(`Selected ${nodeId}. Drag to reposition or remove it.`);
    };

    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
  }

  async function uploadRawFiles() {
    if (!sessionChecked) {
      setUploadError("Checking session. Try again in a moment.");
      return;
    }
    if (!uploadHeaders) {
      setUploadError("Sign in first.");
      return;
    }
    if (rawDataFiles.length === 0) {
      setUploadError("Choose one or more raw files.");
      return;
    }
    setUploadLoading(true);
    setUploadError("");
    setUploadMessage("");
    try {
      const uploadedNames: string[] = [];
      for (const file of rawDataFiles) {
        const form = new FormData();
        form.append("file", file);
        const response = await fetch(`${apiBase}/datasets/upload?project_id=${encodeURIComponent(projectId)}`, {
          method: "POST",
          headers: uploadHeaders,
          body: form,
        });
        const payload = await response.json().catch(() => ({}));
        if (!response.ok) {
          if (response.status === 401) {
            handleUnauthorized(typeof payload?.detail === "string" ? payload.detail : "Invalid token");
            return;
          }
          throw new Error(typeof payload?.detail === "string" ? payload.detail : `Upload failed (${response.status})`);
        }
        uploadedNames.push(payload?.dataset?.filename ?? file.name);
      }
      setUploadMessage(`Uploaded ${uploadedNames.length} raw file(s): ${uploadedNames.join(", ")}`);
      const refresh = await fetch(`${apiBase}/datasets`, { headers: uploadHeaders });
      if (refresh.ok) {
        const body = await refresh.json().catch(() => ({}));
        const rows = Array.isArray(body?.datasets) ? body.datasets : [];
        setProjectDatasets(rows.filter((row: DatasetItem) => row.project_id === projectId));
      }
      clearRawFiles();
    } catch (err: any) {
      setUploadError(err?.message || "Upload failed");
    } finally {
      setUploadLoading(false);
    }
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
    if (!sessionChecked) {
      setError("Checking session. Try again in a moment.");
      return;
    }
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
        if (validate.status === 401) {
          handleUnauthorized("Invalid token");
          return;
        }
        const text = await validate.text();
        setError(`Validate failed: ${validate.status} ${text}`);
        return;
      }
      const body = await validate.json().catch(() => ({} as { local_export_path?: string }));
      setStatus(
        body.local_export_path
          ? `Workflow refreshed in PostgreSQL and ${body.local_export_path}.`
          : `Workflow refreshed in PostgreSQL as ${workflowId}.`,
      );
      return;
    }
    if (!create.ok) {
      if (create.status === 401) {
        handleUnauthorized("Invalid token");
        return;
      }
      const text = await create.text();
      setError(`Import failed: ${create.status} ${text}`);
      return;
    }
    const body = await create.json().catch(() => ({} as { local_export_path?: string }));
    setStatus(
      body.local_export_path
        ? `Workflow saved to PostgreSQL and ${body.local_export_path}.`
        : `Workflow saved to PostgreSQL as ${workflowId}.`,
    );
  }

  async function runDistributed() {
    if (!sessionChecked) {
      setError("Checking session. Try again in a moment.");
      return;
    }
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
      if (res.status === 401) {
        handleUnauthorized("Invalid token");
        return;
      }
      const payload = await res.json().catch(() => null);
      const detail = payload?.detail;
      if (detail && typeof detail === "object" && Array.isArray(detail.missing_inputs)) {
        const missing = detail.missing_inputs
          .map((item: { node_id?: string; port?: string; type_id?: string }) => `${item.node_id}.${item.port} (${item.type_id})`)
          .join(", ");
        setError(`Run blocked. Missing validated inputs: ${missing}`);
        return;
      }
      const text = typeof detail === "string" ? detail : await res.text();
      setError(`Run failed: ${res.status} ${text}`);
      return;
    }
    const body = await res.json();
    setMonitor(body.summary);
    setRunDetails(body);
    setShowRawResult(false);
    setStatus(
      body.summary?.local_record_path
        ? `Run completed. Record saved at ${body.summary.local_record_path}.`
        : `Distributed run completed: ${body.summary.run_id}`,
    );
  }

  async function refreshRun() {
    if (!sessionChecked) return;
    if (!authHeaders || !monitor?.run_id) return;
    setError("");
    const res = await fetch(`${apiBase}/workflows/runs/${monitor.run_id}`, { headers: authHeaders });
    if (!res.ok) {
      if (res.status === 401) {
        handleUnauthorized("Invalid token");
        return;
      }
      setError(`Refresh run failed: ${res.status}`);
      return;
    }
    const body = await res.json();
    setRunDetails(body);
    setShowRawResult(false);
    setStatus("Run details refreshed.");
  }

  async function exportReproReport() {
    if (!sessionChecked) {
      setError("Checking session. Try again in a moment.");
      return;
    }
    if (!authHeaders) {
      setError("Please sign in first.");
      return;
    }
    setError("");
    const exportRes = await fetch(`${apiBase}/workflows/${workflowId}/export`, { headers: authHeaders });
    if (!exportRes.ok) {
      if (exportRes.status === 401) {
        handleUnauthorized("Invalid token");
        return;
      }
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
    setError("");
    setStatus("Canvas cleared.");
  }

  async function refreshProjectWorkflows() {
    if (!sessionChecked) {
      setError("Checking session. Try again in a moment.");
      return;
    }
    if (!authHeaders) {
      setError("Please sign in first.");
      return;
    }
    setError("");
    const res = await fetch(`${apiBase}/workflows`, { headers: authHeaders });
    if (!res.ok) {
      if (res.status === 401) {
        handleUnauthorized("Invalid token");
        return;
      }
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
    if (!sessionChecked) {
      setError("Checking session. Try again in a moment.");
      return;
    }
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
      if (res.status === 401) {
        handleUnauthorized("Invalid token");
        return;
      }
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
    setSweepsRaw(JSON.stringify(loadedSweeps, null, 2));
    setStatus(`Loaded workflow ${selectedSavedWorkflow}.`);
  }

  if (!authReady) {
    return (
      <div className="flex min-h-screen items-center justify-center bg-background px-4">
        <Card className="w-full max-w-md">
          <CardHeader>
            <CardTitle>Checking session</CardTitle>
            <CardDescription>Preparing workflow builder.</CardDescription>
          </CardHeader>
        </Card>
      </div>
    );
  }

  return (
    <AppShell
      title="Workflow Builder"
      subtitle="Upload raw files, load the fixed WGS template, then save and run."
      badge={!authReady || (token && !sessionChecked) ? 'checking session' : token ? workflowId || 'signed in' : 'sign in'}
      navItems={[
        { label: 'Console', href: '/' },
        { label: 'Workflow Builder', href: '/workflow-builder', active: true },
        { label: 'Locus Explorer', href: `/locus/${encodeURIComponent('chr1:1-1000')}` },
        token ? { label: 'Logout', onClick: logoutSession, variant: 'ghost' } : { label: 'Sign In', href: '/login' },
      ]}
    >
      <section className="grid gap-4 md:grid-cols-3">
        <Card>
          <CardHeader>
            <CardDescription>Project</CardDescription>
            <CardTitle>{projectId || 'unset'}</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-muted-foreground">{workflowId || 'Name this workflow to create a scoped id.'}</p>
          </CardContent>
        </Card>
        <Card>
          <CardHeader>
            <CardDescription>Canvas</CardDescription>
            <CardTitle>{nodes.length} nodes / {edges.length} edges</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-muted-foreground">
              {nodes.length === 0
                ? 'Load template to start.'
                : validationErrors.length === 0
                  ? 'Ready to save and run.'
                  : `${validationErrors.length} validation issue(s).`}
            </p>
          </CardContent>
        </Card>
        <Card>
          <CardHeader>
            <CardDescription>Execution</CardDescription>
            <CardTitle>{monitor ? monitor.run_id : 'No run yet'}</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-muted-foreground">
              {monitor ? `${monitor.completed_runs}/${monitor.submitted_runs} completed` : 'Save and run from the right panel.'}
            </p>
          </CardContent>
        </Card>
      </section>

      <section className="grid gap-6 xl:grid-cols-[1.05fr_0.95fr]">
        <Card>
          <CardHeader>
            <CardTitle className="flex items-center gap-2"><FlaskConical className="h-5 w-5 text-primary" />Upload Raw</CardTitle>
            <CardDescription>Add raw files first. You can select multiple at once or append across folders.</CardDescription>
          </CardHeader>
          <CardContent className="grid gap-4">
            <Input
              ref={rawUploadInputRef}
              type="file"
              multiple
              accept=".fasta,.fa,.fna,.gtf,.gff,.gff3,.fastq,.fq,.fastq.gz,.fq.gz,.bam,.sam,.vcf.gz,.tsv"
              onChange={(e) => appendRawFiles(e.target.files)}
            />
            <p className="text-sm text-muted-foreground">
              {rawDataFiles.length > 0 ? `${rawDataFiles.length} file(s) selected` : 'No raw files selected yet.'}
            </p>
            {rawDataFiles.length > 0 ? (
              <div className="selected-files">
                {rawDataFiles.map((file) => (
                  <div key={fileKey(file)} className="selected-file">
                    <span>{file.name}</span>
                    <Button type="button" size="sm" variant="ghost" onClick={() => removeRawFile(fileKey(file))}>
                      Remove
                    </Button>
                  </div>
                ))}
              </div>
            ) : null}
            <div className="flex flex-wrap gap-2">
              <Button type="button" onClick={() => void uploadRawFiles()} disabled={!sessionChecked || uploadLoading || rawDataFiles.length === 0}>
                {uploadLoading ? 'Uploading' : 'Upload Raw'}
              </Button>
              <Button type="button" variant="outline" onClick={clearRawFiles} disabled={!sessionChecked || uploadLoading || rawDataFiles.length === 0}>
                Clear
              </Button>
            </div>
            {uploadError ? <p className="text-sm text-destructive">{uploadError}</p> : null}
            {uploadMessage ? <p className="text-sm text-emerald-700">{uploadMessage}</p> : null}
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>Flow</CardTitle>
            <CardDescription>Use the shortest path: upload, load template, save, run, review result.</CardDescription>
          </CardHeader>
          <CardContent className="grid gap-3 md:grid-cols-2">
            <div className="rounded-lg border border-border bg-muted/40 p-4">
              <p className="mb-1 text-sm font-medium text-foreground">1. Upload</p>
              <p className="text-sm text-muted-foreground">Add FASTQ, reference, or other raw files for this project.</p>
            </div>
            <div className="rounded-lg border border-border bg-muted/40 p-4">
              <p className="mb-1 text-sm font-medium text-foreground">2. Load template</p>
              <p className="text-sm text-muted-foreground">Start from the fixed WGS block set instead of wiring from scratch.</p>
            </div>
            <div className="rounded-lg border border-border bg-muted/40 p-4">
              <p className="mb-1 text-sm font-medium text-foreground">3. Save</p>
              <p className="text-sm text-muted-foreground">Store the workflow under its `workflow_id` so it can be loaded again.</p>
            </div>
            <div className="rounded-lg border border-border bg-muted/40 p-4">
              <p className="mb-1 text-sm font-medium text-foreground">4. Run</p>
              <p className="text-sm text-muted-foreground">Execute the generated workflow and inspect the saved run directory and outputs.</p>
            </div>
          </CardContent>
        </Card>
      </section>

      <section className="grid gap-6">
        <Card>
          <CardHeader>
            <CardTitle>Project Inputs</CardTitle>
            <CardDescription>Only validated files can be bound into the workflow.</CardDescription>
          </CardHeader>
          <CardContent className="grid gap-3">
            {projectDatasets.length === 0 ? (
              <p className="text-sm text-muted-foreground">No validated project inputs yet for this project.</p>
            ) : (
              projectDatasets.map((dataset) => (
                <div key={dataset.id} className="dataset-row">
                  <div>
                    <p className="font-medium text-foreground">{dataset.filename}</p>
                    <p className="text-sm text-muted-foreground">{dataset.input_role} | {dataset.canonical_type || "unknown type"}</p>
                  </div>
                  <div className="dataset-meta">
                    <span className={`dataset-badge ${dataset.validation_status === "validated" ? "dataset-badge-ok" : "dataset-badge-warn"}`}>
                      {dataset.validation_status}
                    </span>
                    <span className="text-sm text-muted-foreground">{dataset.validation_detail}</span>
                  </div>
                </div>
              ))
            )}
          </CardContent>
        </Card>
      </section>

      <section className="grid gap-6 xl:grid-cols-[260px_minmax(0,1fr)_340px]">
        <Card className="self-start">
          <CardHeader>
            <CardTitle className="flex items-center gap-2"><Workflow className="h-5 w-5 text-primary" />Blocks</CardTitle>
            <CardDescription>Fixed WGS nodes only. Double-click a block to add it.</CardDescription>
          </CardHeader>
          <CardContent className="grid gap-3">
            <div className="mini-note">FastQC / Cutadapt / BWA / samtools / GATK / Report</div>
            {filteredPaletteTools.map((p) => (
              <div
                key={`${p.plugin_id}:${p.version}`}
                draggable
                onDragStart={(e) => e.dataTransfer.setData('text/plain', p.plugin_id)}
                onDoubleClick={() => addNode(p)}
                className={`palette-item ${paletteColorClass(p)}`}
              >
                <strong>{p.name}</strong>
                <div>{p.plugin_id}</div>
              </div>
            ))}
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-start justify-between gap-4 space-y-0">
            <div className="space-y-1.5">
              <CardTitle>Canvas</CardTitle>
              <CardDescription>Drag nodes to arrange them. Use template or auto-connect for links.</CardDescription>
            </div>
            <div className="canvas-controls">
              <Button type="button" variant="outline" size="sm" onClick={() => setCanvasZoom((z) => Math.max(0.6, Math.round((z - 0.1) * 100) / 100))}>-</Button>
              <Badge variant="outline">{Math.round(canvasZoom * 100)}%</Badge>
              <Button type="button" variant="outline" size="sm" onClick={() => setCanvasZoom((z) => Math.min(1.8, Math.round((z + 0.1) * 100) / 100))}>+</Button>
            </div>
          </CardHeader>
          <CardContent>
            <div onDragOver={(e) => e.preventDefault()} onDrop={onDropPlugin} className="canvas-box">
              {nodes.length === 0 ? <p className="text-sm text-muted-foreground">No blocks yet.</p> : null}
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
                        className={`canvas-stage ${isDraggingNode ? 'dragging' : ''}`}
                        style={{
                          width: `${CANVAS_WIDTH}px`,
                          height: `${CANVAS_HEIGHT}px`,
                          transform: `scale(${canvasZoom})`,
                        }}
                      >
                        <svg className="canvas-links" width={CANVAS_WIDTH} height={CANVAS_HEIGHT} viewBox={`0 0 ${CANVAS_WIDTH} ${CANVAS_HEIGHT}`}>
                          {canvasEdgePaths.map((edge) => (
                            <path key={edge.key} d={edge.d} className="mind-link" onDoubleClick={() => removeEdgeBySignature(edge.key)} />
                          ))}
                        </svg>
                        {canvasNodes.map((node) => (
                          <div
                            key={node.nodeId}
                            className={`canvas-node ${node.tone}`}
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
                  <p className="canvas-tip">Drag to arrange. Double-click a line to remove it.</p>
                </div>
              ) : null}
            </div>
          </CardContent>
        </Card>

        <div className="grid gap-6">
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2"><Database className="h-5 w-5 text-primary" />Workflow</CardTitle>
              <CardDescription>Set identifiers, then save and run.</CardDescription>
            </CardHeader>
            <CardContent className="grid gap-4">
              <div className="flex flex-wrap gap-2">
                <Button onClick={loadWgsTemplate} variant="secondary">Load Template</Button>
                <Button onClick={autoConnectSequence} variant="outline">Auto-connect</Button>
                <Button onClick={() => void refreshProjectWorkflows()} variant="outline" disabled={!sessionChecked}>Refresh</Button>
                <Button onClick={clearCanvas} variant="ghost">Clear Canvas</Button>
              </div>
              <div className="space-y-2">
                <Label htmlFor="builder-project-id">Project ID</Label>
                <Input id="builder-project-id" value={projectId} onChange={(e) => setProjectId(e.target.value)} />
              </div>
              <div className="space-y-2">
                <Label htmlFor="builder-workflow-name">Workflow Name</Label>
                <Input id="builder-workflow-name" value={workflowName} onChange={(e) => setWorkflowName(e.target.value)} />
              </div>
              <div className="space-y-2">
                <Label htmlFor="builder-workflow-id">Workflow ID</Label>
                <Input id="builder-workflow-id" value={workflowId} readOnly className="readonly" />
              </div>
              <div className="flex flex-wrap gap-2">
                <Button onClick={() => void saveAndValidateWorkflow()} disabled={!sessionChecked}>Save</Button>
                <Button onClick={() => void runDistributed()} disabled={!sessionChecked}><Play className="h-4 w-4" />Run</Button>
              </div>
              <p className="text-sm text-muted-foreground">Save persists the workflow to PostgreSQL and `results/workflows/definitions/`.</p>
              <div className="grid gap-2">
                <Label htmlFor="saved-workflows">Saved workflows</Label>
                <Select id="saved-workflows" value={selectedSavedWorkflow} onChange={(e) => setSelectedSavedWorkflow(e.target.value)}>
                  <option value="">Saved workflows</option>
                  {savedWorkflows.map((w) => (
                    <option key={w.workflow_id} value={w.workflow_id}>
                      {w.workflow_id} ({w.node_count} nodes, {w.edge_count} edges)
                    </option>
                  ))}
                </Select>
                <div className="flex flex-wrap gap-2">
                  <Button onClick={() => void loadSelectedWorkflowFromServer()} variant="outline" disabled={!sessionChecked}>Load</Button>
                  <Button onClick={() => void exportReproReport()} variant="outline" disabled={!sessionChecked}>Export Report</Button>
                </div>
              </div>
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle>Validation</CardTitle>
              <CardDescription>Only fixed WGS-compatible links are allowed.</CardDescription>
            </CardHeader>
            <CardContent>
              {validationErrors.length === 0 ? (
                <p className="text-sm text-emerald-700">Ready.</p>
              ) : (
                <div className="grid gap-2">
                  <p className="text-sm text-destructive">{validationErrors.length} issue(s).</p>
                  <ul className="tight-list text-sm text-muted-foreground">
                    {validationErrors.map((issue, i) => (
                      <li key={`issue-${i}`}>{issue}</li>
                    ))}
                  </ul>
                </div>
              )}
            </CardContent>
          </Card>

          <Card>
            <CardHeader className="flex flex-row items-start justify-between gap-4 space-y-0">
              <div className="space-y-1.5">
                <CardTitle className="flex items-center gap-2"><Activity className="h-5 w-5 text-primary" />Result</CardTitle>
                <CardDescription>Run summary first, raw JSON only when needed.</CardDescription>
              </div>
              <Button onClick={() => void refreshRun()} variant="outline" size="sm" disabled={!sessionChecked || !monitor?.run_id}>
                Refresh
              </Button>
            </CardHeader>
            <CardContent className="grid gap-4">
              {monitor ? (
                <div className="grid gap-2 text-sm text-muted-foreground">
                  <div>Run ID: <span className="text-foreground">{monitor.run_id}</span></div>
                  <div>Done: <span className="text-foreground">{monitor.completed_runs}/{monitor.submitted_runs}</span></div>
                  <div>Failed: <span className="text-foreground">{monitor.failed_runs || 0}</span></div>
                  <div>Duration: <span className="text-foreground">{monitor.duration_ms} ms</span></div>
                  {runRecordPath ? <div>Run Record: <span className="path-text">{runRecordPath}</span></div> : null}
                </div>
              ) : (
                <p className="text-sm text-muted-foreground">No run yet.</p>
              )}
              {primaryRunResult ? (
                <div className="result-summary">
                  <div><strong>Status:</strong> {primaryRunResult.status || 'unknown'}</div>
                  <div><strong>Engine:</strong> {primaryRunResult.engine || 'unknown'}</div>
                  <div><strong>Run Dir:</strong> <span className="path-text">{primaryRunResult.workflow_artifacts?.run_dir || 'n/a'}</span></div>
                  <div><strong>Snakefile:</strong> <span className="path-text">{primaryRunResult.workflow_artifacts?.snakefile || 'n/a'}</span></div>
                  <div><strong>Config:</strong> <span className="path-text">{primaryRunResult.workflow_artifacts?.configfile || 'n/a'}</span></div>
                  <div><strong>Outputs:</strong> {primaryRunResult.outputs?.length || 0}</div>
                  {primaryRunResult.outputs && primaryRunResult.outputs.length > 0 ? (
                    <ul className="tight-list">
                      {primaryRunResult.outputs.slice(0, 5).map((output) => (
                        <li key={output} className="path-text">{output}</li>
                      ))}
                    </ul>
                  ) : null}
                </div>
              ) : null}
              <div className="flex flex-wrap gap-2">
                <Button type="button" variant="outline" onClick={() => setShowRawResult((open) => !open)} disabled={!runDetails}>
                  {showRawResult ? 'Hide Raw JSON' : 'Show Raw JSON'}
                </Button>
              </div>
              {showRawResult ? <pre className="code-box">{JSON.stringify(runDetails || {}, null, 2)}</pre> : null}
            </CardContent>
          </Card>
        </div>
      </section>

      {error ? <p className="text-sm text-destructive">{error}</p> : null}
      {status ? <p className="text-sm text-emerald-700">{status}</p> : null}

      <style jsx>{`
        .mini-note {
          border: 1px solid #d8e2e7;
          border-radius: 10px;
          padding: 9px 10px;
          background: #f6fbfb;
          color: #23404e;
          font-size: 12px;
        }
        .selected-files {
          display: grid;
          gap: 8px;
        }
        .selected-file {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 8px;
          border: 1px solid #d8e2e7;
          border-radius: 10px;
          padding: 8px 10px;
          background: #f7fafb;
          font-size: 12px;
          color: #10222c;
        }
        .selected-file span {
          overflow: hidden;
          text-overflow: ellipsis;
          white-space: nowrap;
        }
        .dataset-row {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 12px;
          border: 1px solid #d8e2e7;
          border-radius: 12px;
          background: #fbfdfc;
          padding: 12px;
        }
        .dataset-meta {
          display: grid;
          justify-items: end;
          gap: 4px;
        }
        .dataset-badge {
          display: inline-flex;
          align-items: center;
          border-radius: 999px;
          padding: 4px 8px;
          font-size: 11px;
          font-weight: 600;
        }
        .dataset-badge-ok {
          background: #e9f8ee;
          color: #135a2e;
        }
        .dataset-badge-warn {
          background: #fff0ea;
          color: #972315;
        }
        .readonly {
          background: #f7fafb;
        }
        .path-text {
          display: inline-block;
          max-width: 100%;
          word-break: break-all;
          color: #23404e;
        }
        .tight-list {
          margin: 0;
          padding-left: 18px;
          display: grid;
          gap: 4px;
        }
        .result-summary {
          display: grid;
          gap: 6px;
          border: 1px solid #d8e2e7;
          border-radius: 10px;
          padding: 10px;
          background: #f7fafb;
          font-size: 13px;
        }
        .palette-item {
          border: 1px solid #d8e2e7;
          border-radius: 12px;
          padding: 12px;
          display: grid;
          gap: 4px;
          background: #fff;
          cursor: grab;
          transition: transform 0.18s ease, box-shadow 0.18s ease;
        }
        .palette-item:hover {
          transform: translateY(-1px);
          box-shadow: 0 10px 20px rgba(20, 45, 58, 0.08);
        }
        .palette-wgs {
          background: linear-gradient(180deg, #ffffff 0%, #eef8f5 100%);
        }
        .palette-registry {
          background: linear-gradient(180deg, #ffffff 0%, #f7f1e8 100%);
        }
        .palette-generic,
        .palette-rna,
        .palette-mixed {
          background: linear-gradient(180deg, #ffffff 0%, #f5f7fb 100%);
        }
        .canvas-box {
          border: 1px solid #ccd6dc;
          border-radius: 14px;
          background: linear-gradient(180deg, #ffffff 0%, #f7faf8 100%);
          padding: 12px;
          min-height: 520px;
        }
        .canvas-controls {
          display: flex;
          flex-wrap: wrap;
          align-items: center;
          gap: 8px;
        }
        .mindmap-shell {
          display: grid;
          gap: 8px;
        }
        .mindmap-viewport {
          overflow: auto;
          border-radius: 12px;
          border: 1px solid #d8e2e7;
          background:
            radial-gradient(circle at top, rgba(208, 234, 228, 0.55), transparent 24%),
            linear-gradient(180deg, #fbfdfc 0%, #f3f6f5 100%);
        }
        .mindmap-zoom-layer {
          position: relative;
        }
        .canvas-stage {
          position: relative;
          transform-origin: top left;
        }
        .canvas-stage.dragging {
          cursor: grabbing;
        }
        .canvas-links {
          position: absolute;
          inset: 0;
          overflow: visible;
          pointer-events: none;
        }
        .mind-link {
          fill: none;
          stroke: #0f7d76;
          stroke-width: 3;
          stroke-linecap: round;
          pointer-events: stroke;
          cursor: pointer;
        }
        .canvas-node {
          position: absolute;
          width: ${NODE_WIDTH}px;
          min-height: ${NODE_HEIGHT}px;
          border-radius: 14px;
          border: 1px solid #ccd6dc;
          padding: 12px;
          background: #fff;
          box-shadow: 0 10px 24px rgba(20, 45, 58, 0.08);
          cursor: grab;
          user-select: none;
        }
        .node-wgs {
          background: linear-gradient(180deg, #ffffff 0%, #eef8f5 100%);
        }
        .node-registry {
          background: linear-gradient(180deg, #ffffff 0%, #f7f1e8 100%);
        }
        .node-generic,
        .node-mixed,
        .node-rna {
          background: linear-gradient(180deg, #ffffff 0%, #f5f7fb 100%);
        }
        .node-title-row {
          display: flex;
          flex-direction: column;
          gap: 2px;
        }
        .node-title-row strong {
          font-size: 13px;
          color: #10222c;
        }
        .node-title-row span {
          font-size: 12px;
          color: #4d626e;
        }
        .node-plugin {
          margin: 8px 0 0;
          font-size: 12px;
          color: #4d626e;
        }
        .node-actions {
          margin-top: 10px;
          display: flex;
          justify-content: flex-end;
        }
        .small-btn {
          border: 1px solid #ccd6dc;
          border-radius: 999px;
          padding: 4px 8px;
          background: #fff;
          font-size: 11px;
        }
        .danger-btn {
          color: #9d2d1d;
        }
        .canvas-tip {
          margin: 0;
          font-size: 12px;
          color: #4d626e;
        }
        .code-box {
          margin: 0;
          max-height: 280px;
          overflow: auto;
          border-radius: 10px;
          border: 1px solid #d8e2e7;
          background: #10222c;
          color: #ecf3f6;
          padding: 12px;
          font-size: 12px;
        }
      `}</style>
    </AppShell>
  );
}
