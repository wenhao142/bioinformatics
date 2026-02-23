"use client";

import { useCallback, useEffect, useMemo, useState } from "react";

type Health = { status: string };

type GeneScore = { gene: string; score: number; rank: number };

type EvidenceRankResponse = {
  run: { run_id: string; kind: string; project_id: string | null; selected_method?: string };
  result: { ranked_genes: Array<{ gene: string; score: number; rank: number }> };
};

type CausalRunResponse = {
  run: { run_id: string; kind: string; project_id: string | null; selected_method?: string };
  result: { gene_scores: GeneScore[] };
};

type PluginManifest = {
  plugin_id: string;
  name: string;
  enabled: boolean;
  image: string;
  tags: string[];
};

type PluginRun = {
  run_id: string;
  created_at: number;
  engine: string;
  counts: { ranked_genes: number; ranked_loci: number };
};

type PluginCompareResult = {
  overlap_genes: string[];
  only_left_genes: string[];
  only_right_genes: string[];
  score_deltas: Array<{ gene: string; left_score: number; right_score: number; delta: number }>;
};

type ResearchSummaryResult = {
  mode: string;
  summary: string;
  citations: Array<{ id: string; source_type: string; label: string }>;
  citation_ids: string[];
  warnings: string[];
  llm_mode_requested: "auto" | "offline";
  llm_model_used: string | null;
};

type MethodRunResult = {
  kind: "evidence" | "causal" | "plugin";
  runId: string;
  selectedMethod: string;
  topGenes: GeneScore[];
  reportUrl: string | null;
};

const API_PORT_FALLBACK = 18000;

function resolveApiBase(): string {
  const envUrl = process.env.NEXT_PUBLIC_API_URL;
  if (envUrl) return envUrl;
  if (typeof window !== "undefined") {
    const { protocol, hostname } = window.location;
    return `${protocol}//${hostname}:${API_PORT_FALLBACK}`;
  }
  return `http://localhost:${API_PORT_FALLBACK}`;
}

function parseRegion(region: string): { chr: string; start: number; end: number } {
  const match = region.match(/^([^:]+):(\d+)-(\d+)$/);
  if (!match) {
    throw new Error("Region format must be chr:start-end");
  }
  const [, chr, start, end] = match;
  return { chr, start: Number(start), end: Number(end) };
}

function normalizeGenes(raw: string): string[] {
  return raw
    .split(",")
    .map((item) => item.trim())
    .filter((item) => item.length > 0);
}

export default function Home() {
  const [apiStatus, setApiStatus] = useState<string>("checking...");
  const [token, setToken] = useState<string>("");
  const [authEmail, setAuthEmail] = useState<string>("admin@example.com");
  const [authPassword, setAuthPassword] = useState<string>("password");
  const [authLoading, setAuthLoading] = useState(false);
  const [authError, setAuthError] = useState<string | null>(null);
  const [newUserEmail, setNewUserEmail] = useState<string>("");
  const [newUserPassword, setNewUserPassword] = useState<string>("");
  const [registerLoading, setRegisterLoading] = useState(false);
  const [registerError, setRegisterError] = useState<string | null>(null);
  const [registerMessage, setRegisterMessage] = useState<string | null>(null);
  const [region, setRegion] = useState<string>("chr1:1-1000");
  const [projectId, setProjectId] = useState<string>("demo-project");

  const [plugins, setPlugins] = useState<PluginManifest[]>([]);
  const [pluginsLoading, setPluginsLoading] = useState(false);
  const [pluginsError, setPluginsError] = useState<string | null>(null);
  const [selectedPluginId, setSelectedPluginId] = useState<string>("");
  const [pluginRuns, setPluginRuns] = useState<PluginRun[]>([]);
  const [pluginRunsLoading, setPluginRunsLoading] = useState(false);
  const [pluginRunsError, setPluginRunsError] = useState<string | null>(null);
  const [leftRunId, setLeftRunId] = useState<string>("");
  const [rightRunId, setRightRunId] = useState<string>("");
  const [compareLoading, setCompareLoading] = useState(false);
  const [compareError, setCompareError] = useState<string | null>(null);
  const [compareResult, setCompareResult] = useState<PluginCompareResult | null>(null);

  const [variantFile, setVariantFile] = useState<File | null>(null);
  const [exprFile, setExprFile] = useState<File | null>(null);
  const [uploadLoading, setUploadLoading] = useState(false);
  const [uploadError, setUploadError] = useState<string | null>(null);
  const [uploadMessage, setUploadMessage] = useState<string | null>(null);

  const [analysisMethod, setAnalysisMethod] = useState<string>("causal-score");
  const [analysisLoading, setAnalysisLoading] = useState(false);
  const [analysisError, setAnalysisError] = useState<string | null>(null);
  const [methodRunResult, setMethodRunResult] = useState<MethodRunResult | null>(null);
  const [reportError, setReportError] = useState<string | null>(null);

  const [summaryDisease, setSummaryDisease] = useState("Alzheimer disease");
  const [summaryGenesInput, setSummaryGenesInput] = useState("APP, APOE");
  const [summaryLlmMode, setSummaryLlmMode] = useState<"auto" | "offline">("offline");
  const [summaryLlmModel, setSummaryLlmModel] = useState("gpt-4o-mini");
  const [summaryLoading, setSummaryLoading] = useState(false);
  const [summaryError, setSummaryError] = useState<string | null>(null);
  const [summaryResult, setSummaryResult] = useState<ResearchSummaryResult | null>(null);

  const apiBase = useMemo(() => resolveApiBase(), []);

  const authHeaders = useMemo(() => {
    if (!token) return null;
    return { Authorization: `Bearer ${token}` };
  }, [token]);

  useEffect(() => {
    fetch(`${apiBase}/health`)
      .then((r) => r.json())
      .then((data: Health) => setApiStatus(data.status))
      .catch(() => setApiStatus("offline"));
  }, [apiBase]);

  useEffect(() => {
    if (typeof window === "undefined") return;
    const stored = window.localStorage.getItem("ad_api_token");
    if (stored) setToken(stored);
  }, []);

  useEffect(() => {
    if (typeof window === "undefined") return;
    if (!token) {
      window.localStorage.removeItem("ad_api_token");
      return;
    }
    window.localStorage.setItem("ad_api_token", token);
  }, [token]);

  const loginWithPassword = async () => {
    setAuthLoading(true);
    setAuthError(null);
    try {
      const response = await fetch(`${apiBase}/auth/login`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ email: authEmail, password: authPassword }),
      });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Login failed (${response.status})`);
      }
      const accessToken = typeof payload?.access_token === "string" ? payload.access_token : "";
      if (!accessToken) {
        throw new Error("Login response missing token");
      }
      setToken(accessToken);
    } catch (error: any) {
      setAuthError(error?.message || "Login failed");
    } finally {
      setAuthLoading(false);
    }
  };

  const logoutSession = () => {
    setToken("");
    setAuthError(null);
  };

  const registerUser = async () => {
    if (!authHeaders) {
      setRegisterError("Please sign in as admin first.");
      return;
    }
    setRegisterLoading(true);
    setRegisterError(null);
    setRegisterMessage(null);
    try {
      const email = newUserEmail.trim();
      if (!email) {
        throw new Error("New user email is required");
      }
      if (newUserPassword.length < 8) {
        throw new Error("Password must be at least 8 characters");
      }
      const response = await fetch(`${apiBase}/auth/register`, {
        method: "POST",
        headers: { ...authHeaders, "Content-Type": "application/json" },
        body: JSON.stringify({ email, password: newUserPassword }),
      });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Register failed (${response.status})`);
      }
      setRegisterMessage(`Created account: ${payload.email} (${payload.role})`);
      setNewUserEmail("");
      setNewUserPassword("");
    } catch (error: any) {
      setRegisterError(error?.message || "Register failed");
    } finally {
      setRegisterLoading(false);
    }
  };

  const loadPlugins = useCallback(async () => {
    if (!authHeaders) {
      setPlugins([]);
      setSelectedPluginId("");
      return;
    }
    setPluginsLoading(true);
    setPluginsError(null);
    try {
      const response = await fetch(`${apiBase}/plugins`, { headers: authHeaders });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Request failed (${response.status})`);
      }
      const rows = Array.isArray(payload?.plugins) ? (payload.plugins as PluginManifest[]) : [];
      setPlugins(rows);
      if (rows.length > 0 && !rows.find((row) => row.plugin_id === selectedPluginId)) {
        setSelectedPluginId(rows[0].plugin_id);
      }
    } catch (error: any) {
      setPluginsError(error?.message || "Failed to load plugins");
    } finally {
      setPluginsLoading(false);
    }
  }, [apiBase, authHeaders, selectedPluginId]);

  const loadPluginRuns = useCallback(
    async (pluginId: string) => {
      if (!authHeaders || !pluginId) {
        setPluginRuns([]);
        setLeftRunId("");
        setRightRunId("");
        return;
      }
      setPluginRunsLoading(true);
      setPluginRunsError(null);
      setCompareError(null);
      setCompareResult(null);
      try {
        const response = await fetch(`${apiBase}/plugins/${encodeURIComponent(pluginId)}/runs`, {
          headers: authHeaders,
        });
        const payload = await response.json().catch(() => ({}));
        if (!response.ok) {
          throw new Error(typeof payload?.detail === "string" ? payload.detail : `Request failed (${response.status})`);
        }
        const rows = Array.isArray(payload?.runs) ? (payload.runs as PluginRun[]) : [];
        setPluginRuns(rows);
        setLeftRunId(rows[0]?.run_id ?? "");
        setRightRunId(rows[1]?.run_id ?? rows[0]?.run_id ?? "");
      } catch (error: any) {
        setPluginRunsError(error?.message || "Failed to load plugin runs");
      } finally {
        setPluginRunsLoading(false);
      }
    },
    [apiBase, authHeaders],
  );

  useEffect(() => {
    void loadPlugins();
  }, [loadPlugins]);

  useEffect(() => {
    if (!selectedPluginId) {
      setPluginRuns([]);
      return;
    }
    void loadPluginRuns(selectedPluginId);
  }, [loadPluginRuns, selectedPluginId]);

  const analysisMethodOptions = useMemo(() => {
    const pluginOptions = plugins
      .filter((plugin) => plugin.enabled)
      .map((plugin) => ({
        value: `plugin:${plugin.plugin_id}`,
        label: `Plugin: ${plugin.name} (${plugin.plugin_id})`,
      }));
    return [
      { value: "evidence-rank", label: "Baseline Evidence Rank" },
      { value: "causal-score", label: "Causal Score" },
      ...pluginOptions,
    ];
  }, [plugins]);

  useEffect(() => {
    const exists = analysisMethodOptions.some((item) => item.value === analysisMethod);
    if (!exists && analysisMethodOptions.length > 0) {
      setAnalysisMethod(analysisMethodOptions[0].value);
    }
  }, [analysisMethod, analysisMethodOptions]);

  const togglePlugin = async (plugin: PluginManifest) => {
    if (!authHeaders) return;
    setPluginsError(null);
    try {
      const response = await fetch(`${apiBase}/plugins/${encodeURIComponent(plugin.plugin_id)}/enabled`, {
        method: "PATCH",
        headers: { ...authHeaders, "Content-Type": "application/json" },
        body: JSON.stringify({ enabled: !plugin.enabled }),
      });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Request failed (${response.status})`);
      }
      setPlugins((prev) =>
        prev.map((row) => (row.plugin_id === plugin.plugin_id ? { ...row, enabled: !plugin.enabled } : row)),
      );
    } catch (error: any) {
      setPluginsError(error?.message || "Failed to toggle plugin");
    }
  };

  const registerDefaultPlugin = async () => {
    if (!authHeaders) return;
    setPluginsError(null);
    try {
      const payload = {
        plugin_id: "baseline-runner",
        name: "Baseline Evidence Runner",
        version: "1.0.0",
        image: "builtin/baseline-evidence:1.0.0",
        description: "Built-in baseline evidence ranking plugin for quick start.",
        input_schema: {
          type: "object",
          properties: {
            project_id: { type: "string" },
            chr: { type: "string" },
            start: { type: "integer" },
            end: { type: "integer" },
            top_n: { type: "integer" },
          },
          required: ["project_id", "chr", "start", "end"],
          additionalProperties: true,
        },
        output_schema: {
          type: "object",
          properties: {
            ranked_genes: { type: "array" },
            ranked_loci: { type: "array" },
            region: { type: "object" },
            meta: { type: "object" },
          },
          required: ["ranked_genes", "ranked_loci"],
          additionalProperties: true,
        },
        resources: { cpu_millicores: 1000, memory_mb: 1024, gpu: false },
        enabled: true,
        tags: ["baseline", "offline"],
      };
      const response = await fetch(`${apiBase}/plugins/register`, {
        method: "POST",
        headers: { ...authHeaders, "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      const body = await response.json().catch(() => ({}));
      if (!response.ok && response.status !== 409) {
        throw new Error(typeof body?.detail === "string" ? body.detail : `Request failed (${response.status})`);
      }
      await loadPlugins();
      setPluginsError(null);
    } catch (error: any) {
      setPluginsError(error?.message || "Failed to register default plugin");
    }
  };

  const runPluginById = useCallback(
    async (pluginId: string) => {
      if (!authHeaders || !pluginId) {
        throw new Error("Plugin and token are required");
      }
      const { chr, start, end } = parseRegion(region);
      const response = await fetch(`${apiBase}/plugins/${encodeURIComponent(pluginId)}/run`, {
        method: "POST",
        headers: { ...authHeaders, "Content-Type": "application/json" },
        body: JSON.stringify({
          project_id: projectId,
          chr,
          start,
          end,
          top_n: 5,
          parameters: { selected_method: `plugin:${pluginId}` },
        }),
      });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Request failed (${response.status})`);
      }
      return payload;
    },
    [apiBase, authHeaders, projectId, region],
  );

  const runPluginFromManagement = async () => {
    if (!selectedPluginId) return;
    setPluginRunsError(null);
    try {
      await runPluginById(selectedPluginId);
      await loadPluginRuns(selectedPluginId);
    } catch (error: any) {
      setPluginRunsError(error?.message || "Failed to run plugin");
    }
  };

  const compareRuns = async () => {
    if (!authHeaders || !selectedPluginId || !leftRunId || !rightRunId) return;
    setCompareLoading(true);
    setCompareError(null);
    setCompareResult(null);
    try {
      const response = await fetch(`${apiBase}/plugins/${encodeURIComponent(selectedPluginId)}/compare`, {
        method: "POST",
        headers: { ...authHeaders, "Content-Type": "application/json" },
        body: JSON.stringify({ left_run_id: leftRunId, right_run_id: rightRunId }),
      });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Request failed (${response.status})`);
      }
      setCompareResult(payload as PluginCompareResult);
    } catch (error: any) {
      setCompareError(error?.message || "Failed to compare runs");
    } finally {
      setCompareLoading(false);
    }
  };

  const uploadFile = async (kind: "variants" | "expr") => {
    if (!authHeaders) {
      setUploadError("API token required");
      return;
    }
    const file = kind === "variants" ? variantFile : exprFile;
    if (!file) {
      setUploadError(kind === "variants" ? "Please choose a VCF file." : "Please choose an expression TSV file.");
      return;
    }
    setUploadLoading(true);
    setUploadError(null);
    setUploadMessage(null);
    try {
      const form = new FormData();
      form.append("file", file);
      const endpoint =
        kind === "variants"
          ? `${apiBase}/variants/ingest?project_id=${encodeURIComponent(projectId)}`
          : `${apiBase}/omics/expr/upload`;
      const response = await fetch(endpoint, { method: "POST", headers: authHeaders, body: form });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Upload failed (${response.status})`);
      }
      if (kind === "variants") {
        setUploadMessage(`VCF uploaded. Ingested ${payload.ingested ?? 0} variants.`);
      } else {
        setUploadMessage(`Expression table uploaded. Ingested ${payload.ingested ?? 0} rows.`);
      }
    } catch (error: any) {
      setUploadError(error?.message || "Upload failed");
    } finally {
      setUploadLoading(false);
    }
  };

  const runSelectedMethod = async () => {
    if (!authHeaders) {
      setAnalysisError("API token required");
      return;
    }
    setAnalysisLoading(true);
    setAnalysisError(null);
    setMethodRunResult(null);
    try {
      const { chr, start, end } = parseRegion(region);
      if (analysisMethod === "evidence-rank") {
        const response = await fetch(
          `${apiBase}/runs/evidence?chr=${encodeURIComponent(chr)}&start=${start}&end=${end}&top_n=5&project_id=${encodeURIComponent(projectId)}&method=${encodeURIComponent(analysisMethod)}`,
          { method: "POST", headers: authHeaders },
        );
        const payload = (await response.json().catch(() => ({}))) as Partial<EvidenceRankResponse> & {
          detail?: string;
        };
        if (!response.ok || !payload.run || !payload.result) {
          throw new Error(typeof payload.detail === "string" ? payload.detail : `Request failed (${response.status})`);
        }
        setMethodRunResult({
          kind: "evidence",
          runId: payload.run.run_id,
          selectedMethod: payload.run.selected_method || analysisMethod,
          topGenes: payload.result.ranked_genes.slice(0, 5).map((row) => ({
            gene: row.gene,
            score: row.score,
            rank: row.rank,
          })),
          reportUrl: `${apiBase}/report/${encodeURIComponent(projectId)}/${encodeURIComponent(payload.run.run_id)}`,
        });
        return;
      }
      if (analysisMethod === "causal-score") {
        const response = await fetch(
          `${apiBase}/causal/score?chr=${encodeURIComponent(chr)}&start=${start}&end=${end}&top_n=5&project_id=${encodeURIComponent(projectId)}&method=${encodeURIComponent(analysisMethod)}`,
          { method: "POST", headers: authHeaders },
        );
        const payload = (await response.json().catch(() => ({}))) as Partial<CausalRunResponse> & { detail?: string };
        if (!response.ok || !payload.run || !payload.result) {
          throw new Error(typeof payload.detail === "string" ? payload.detail : `Request failed (${response.status})`);
        }
        setMethodRunResult({
          kind: "causal",
          runId: payload.run.run_id,
          selectedMethod: payload.run.selected_method || analysisMethod,
          topGenes: payload.result.gene_scores.slice(0, 5),
          reportUrl: `${apiBase}/report/${encodeURIComponent(projectId)}/${encodeURIComponent(payload.run.run_id)}`,
        });
        return;
      }
      if (analysisMethod.startsWith("plugin:")) {
        const pluginId = analysisMethod.slice("plugin:".length);
        if (!pluginId) {
          throw new Error("Plugin ID is missing");
        }
        const payload = await runPluginById(pluginId);
        const run = payload?.run;
        const result = payload?.result;
        const topGenes = Array.isArray(result?.ranked_genes)
          ? result.ranked_genes.slice(0, 5).map((row: any, idx: number) => ({
              gene: String(row?.gene ?? "NA"),
              score: Number(row?.score ?? 0),
              rank: Number(row?.rank ?? idx + 1),
            }))
          : [];
        setMethodRunResult({
          kind: "plugin",
          runId: String(run?.run_id ?? "unknown"),
          selectedMethod: analysisMethod,
          topGenes,
          reportUrl: null,
        });
        await loadPluginRuns(pluginId);
        return;
      }
      throw new Error("Unsupported method");
    } catch (error: any) {
      setAnalysisError(error?.message || "Failed to run selected method");
    } finally {
      setAnalysisLoading(false);
    }
  };

  const openReportWithAuth = async (format: "markdown" | "html") => {
    if (!authHeaders || !methodRunResult?.reportUrl || !token) {
      setReportError("Report is unavailable for this run.");
      return;
    }
    setReportError(null);
    try {
      const baseUrl = format === "html" ? `${methodRunResult.reportUrl}?format=html` : methodRunResult.reportUrl;
      const separator = baseUrl.includes("?") ? "&" : "?";
      const targetUrl = `${baseUrl}${separator}token=${encodeURIComponent(token)}`;
      const response = await fetch(targetUrl, { headers: authHeaders });
      if (!response.ok) {
        const body = await response.json().catch(() => ({}));
        throw new Error(typeof body?.detail === "string" ? body.detail : `Report request failed (${response.status})`);
      }
      const content = await response.text();
      const mimeType = format === "html" ? "text/html;charset=utf-8" : "text/markdown;charset=utf-8";
      const blob = new Blob([content], { type: mimeType });
      const objectUrl = URL.createObjectURL(blob);
      const popup = window.open(objectUrl, "_blank", "noopener,noreferrer");
      if (!popup) {
        throw new Error("Browser blocked popup. Please allow popups and retry.");
      }
      setTimeout(() => URL.revokeObjectURL(objectUrl), 60_000);
    } catch (error: any) {
      setReportError(error?.message || "Failed to open report");
    }
  };

  const generateResearchSummary = async () => {
    if (!authHeaders) {
      setSummaryError("API token required");
      return;
    }
    setSummaryLoading(true);
    setSummaryError(null);
    setSummaryResult(null);
    try {
      const genes = normalizeGenes(summaryGenesInput);
      if (genes.length === 0) {
        throw new Error("Please provide at least one gene");
      }
      const payload: Record<string, unknown> = {
        disease: summaryDisease,
        genes,
        top_n: Math.min(5, genes.length),
        include_pubmed: true,
        llm_mode: summaryLlmMode,
      };
      if (summaryLlmMode === "auto" && summaryLlmModel.trim()) {
        payload.llm_model = summaryLlmModel.trim();
      }
      const response = await fetch(`${apiBase}/research/summary`, {
        method: "POST",
        headers: { ...authHeaders, "Content-Type": "application/json" },
        body: JSON.stringify(payload),
      });
      const body = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof body?.detail === "string" ? body.detail : `Request failed (${response.status})`);
      }
      setSummaryResult(body as ResearchSummaryResult);
    } catch (error: any) {
      setSummaryError(error?.message || "Failed to generate summary");
    } finally {
      setSummaryLoading(false);
    }
  };

  return (
    <main className="analysis-page">
      <div className="ambient-shape shape-a" />
      <div className="ambient-shape shape-b" />
      <header className="hero-panel reveal reveal-1">
        <div>
          <p className="hero-kicker">AD Multi-Omics Locus Evidence Platform</p>
          <h1 className="hero-title">Analysis Console</h1>
          <p className="hero-note">Manage uploads, run statistical methods, and control LLM summary mode in one page.</p>
        </div>
        <div className={`health-pill ${apiStatus === "ok" ? "health-pill-ok" : "health-pill-warn"}`}>API health: {apiStatus}</div>
      </header>

      <section style={{ padding: "1rem", border: "1px solid #ddd", borderRadius: 8, display: "grid", gap: "0.6rem" }}>
        <h2 style={{ margin: 0 }}>Session</h2>
        <div style={{ display: "grid", gap: "0.45rem", gridTemplateColumns: "repeat(auto-fit, minmax(180px, 1fr))" }}>
          <input value={authEmail} onChange={(event) => setAuthEmail(event.target.value)} placeholder="email" />
          <input type="password" value={authPassword} onChange={(event) => setAuthPassword(event.target.value)} placeholder="password" />
          <button onClick={() => void loginWithPassword()} disabled={authLoading}>
            {authLoading ? "Signing In..." : token ? "Re-login" : "Sign In"}
          </button>
          <button onClick={logoutSession} disabled={!token}>
            Logout
          </button>
        </div>
        {authError ? <p style={{ margin: 0, color: "#b42318" }}>{authError}</p> : null}
        <p style={{ margin: 0, fontSize: 12, color: "#4d626e" }}>
          {token ? "Logged in. Token is stored in browser localStorage." : "Not logged in."}
        </p>
        <div
          style={{
            display: "grid",
            gap: "0.45rem",
            gridTemplateColumns: "repeat(auto-fit, minmax(170px, 1fr))",
            border: "1px solid #d6dfe4",
            borderRadius: 10,
            padding: "0.55rem",
          }}
        >
          <input
            value={newUserEmail}
            onChange={(event) => setNewUserEmail(event.target.value)}
            placeholder="new user email"
          />
          <input
            type="password"
            value={newUserPassword}
            onChange={(event) => setNewUserPassword(event.target.value)}
            placeholder="new user password"
          />
          <button onClick={() => void registerUser()} disabled={registerLoading || !token}>
            {registerLoading ? "Creating..." : "Create Account"}
          </button>
        </div>
        {registerError ? <p style={{ margin: 0, color: "#b42318" }}>{registerError}</p> : null}
        {registerMessage ? <p style={{ margin: 0, color: "#14532d" }}>{registerMessage}</p> : null}
        <input value={region} onChange={(event) => setRegion(event.target.value)} placeholder="chr:start-end" />
        <input value={projectId} onChange={(event) => setProjectId(event.target.value)} placeholder="project id" />
      </section>

      <section style={{ padding: "1rem", border: "1px solid #ddd", borderRadius: 8, display: "grid", gap: "0.6rem" }}>
        <h2 style={{ margin: 0 }}>Upload Data (T8.2)</h2>
        <p style={{ margin: 0 }}>Upload your own files from UI. No CLI required.</p>
        <div style={{ display: "grid", gap: "0.35rem" }}>
          <label>
            VCF file
            <input type="file" accept=".vcf,.vcf.gz,text/vcf" onChange={(e) => setVariantFile(e.target.files?.[0] || null)} />
          </label>
          <button disabled={uploadLoading || !variantFile} onClick={() => void uploadFile("variants")}>
            {uploadLoading ? "Uploading..." : "Upload VCF"}
          </button>
          <label>
            Expression TSV
            <input type="file" accept=".tsv,text/tab-separated-values,text/plain" onChange={(e) => setExprFile(e.target.files?.[0] || null)} />
          </label>
          <button disabled={uploadLoading || !exprFile} onClick={() => void uploadFile("expr")}>
            {uploadLoading ? "Uploading..." : "Upload Expression"}
          </button>
        </div>
        {uploadError ? <p style={{ margin: 0, color: "#b42318" }}>{uploadError}</p> : null}
        {uploadMessage ? <p style={{ margin: 0, color: "#14532d" }}>{uploadMessage}</p> : null}
      </section>

      <section style={{ padding: "1rem", border: "1px solid #ddd", borderRadius: 8, display: "grid", gap: "0.6rem" }}>
        <h2 style={{ margin: 0 }}>Select Statistical Method (T8.1)</h2>
        <p style={{ margin: 0 }}>Choose one method and run on the selected region.</p>
        <label>
          Method
          <select value={analysisMethod} onChange={(event) => setAnalysisMethod(event.target.value)}>
            {analysisMethodOptions.map((option) => (
              <option key={option.value} value={option.value}>
                {option.label}
              </option>
            ))}
          </select>
        </label>
        <button disabled={analysisLoading} onClick={() => void runSelectedMethod()}>
          {analysisLoading ? "Running..." : "Run Selected Method"}
        </button>
        {analysisError ? <p style={{ margin: 0, color: "#b42318" }}>{analysisError}</p> : null}
        {methodRunResult ? (
          <div style={{ border: "1px solid #d0d7de", borderRadius: 8, padding: "0.75rem" }}>
            <p style={{ margin: "0 0 0.25rem" }}>
              <strong>Run:</strong> {methodRunResult.runId}
            </p>
            <p style={{ margin: "0 0 0.25rem" }}>
              <strong>Method:</strong> {methodRunResult.selectedMethod}
            </p>
            <p style={{ margin: "0 0 0.25rem" }}>
              <strong>Type:</strong> {methodRunResult.kind}
            </p>
            {methodRunResult.reportUrl ? (
              <p style={{ margin: "0 0 0.25rem" }}>
                <strong>Report:</strong>{" "}
                <button type="button" onClick={() => void openReportWithAuth("markdown")} style={{ minHeight: 30, padding: "4px 8px" }}>
                  markdown
                </button>{" "}
                /{" "}
                <button type="button" onClick={() => void openReportWithAuth("html")} style={{ minHeight: 30, padding: "4px 8px" }}>
                  html
                </button>
              </p>
            ) : null}
            {reportError ? <p style={{ margin: 0, color: "#b42318" }}>{reportError}</p> : null}
            <p style={{ margin: "0.35rem 0" }}>
              <strong>Top genes:</strong>
            </p>
            <ul style={{ marginTop: 0 }}>
              {methodRunResult.topGenes.map((row) => (
                <li key={`${methodRunResult.runId}-${row.gene}-${row.rank}`}>
                  #{row.rank} {row.gene} ({row.score.toFixed(3)})
                </li>
              ))}
            </ul>
          </div>
        ) : null}
      </section>

      <section style={{ padding: "1rem", border: "1px solid #ddd", borderRadius: 8, display: "grid", gap: "0.6rem" }}>
        <h2 style={{ margin: 0 }}>Choose LLM Method (T8.3)</h2>
        <p style={{ margin: 0 }}>Switch between offline template and cloud LLM summarization.</p>
        <input value={summaryDisease} onChange={(event) => setSummaryDisease(event.target.value)} placeholder="Disease" />
        <input value={summaryGenesInput} onChange={(event) => setSummaryGenesInput(event.target.value)} placeholder="Genes (comma-separated)" />
        <label>
          LLM mode
          <select value={summaryLlmMode} onChange={(event) => setSummaryLlmMode(event.target.value as "auto" | "offline")}>
            <option value="offline">Offline template</option>
            <option value="auto">Auto cloud LLM (if enabled)</option>
          </select>
        </label>
        {summaryLlmMode === "auto" ? (
          <input value={summaryLlmModel} onChange={(event) => setSummaryLlmModel(event.target.value)} placeholder="Cloud model (optional)" />
        ) : null}
        <button disabled={summaryLoading} onClick={() => void generateResearchSummary()}>
          {summaryLoading ? "Generating..." : "Generate Summary"}
        </button>
        {summaryError ? <p style={{ margin: 0, color: "#b42318" }}>{summaryError}</p> : null}
        {summaryResult ? (
          <div style={{ border: "1px solid #d0d7de", borderRadius: 8, padding: "0.75rem", display: "grid", gap: "0.35rem" }}>
            <p style={{ margin: 0 }}>
              <strong>Returned mode:</strong> {summaryResult.mode}
            </p>
            <p style={{ margin: 0 }}>
              <strong>Requested mode:</strong> {summaryResult.llm_mode_requested}
            </p>
            <p style={{ margin: 0 }}>
              <strong>Model used:</strong> {summaryResult.llm_model_used || "none"}
            </p>
            <p style={{ margin: 0 }}>
              <strong>Summary:</strong> {summaryResult.summary}
            </p>
            <p style={{ margin: "0.3rem 0 0" }}>
              <strong>Citations:</strong> {summaryResult.citation_ids.join(", ") || "none"}
            </p>
            {summaryResult.warnings.length > 0 ? (
              <ul style={{ margin: 0, color: "#9a3412" }}>
                {summaryResult.warnings.map((warning) => (
                  <li key={warning}>{warning}</li>
                ))}
              </ul>
            ) : null}
          </div>
        ) : null}
      </section>

      <section style={{ padding: "1rem", border: "1px solid #ddd", borderRadius: 8 }}>
        <h2 style={{ marginTop: 0 }}>Method Management (T6.3)</h2>
        <p style={{ marginTop: 0 }}>Enable/disable plugin methods and compare method run outputs.</p>
        {pluginsError ? <p style={{ color: "#b42318" }}>{pluginsError}</p> : null}
        {pluginsLoading ? (
          <p>Loading plugins...</p>
        ) : (
          <div style={{ display: "grid", gap: "0.5rem" }}>
            {plugins.map((plugin) => (
              <div
                key={plugin.plugin_id}
                style={{
                  border: "1px solid #d0d7de",
                  borderRadius: 8,
                  padding: "0.6rem",
                  display: "grid",
                  gap: "0.35rem",
                }}
              >
                <div style={{ display: "flex", gap: "0.5rem", alignItems: "center", flexWrap: "wrap" }}>
                  <strong>{plugin.name}</strong>
                  <code>{plugin.plugin_id}</code>
                  <span>{plugin.enabled ? "Enabled" : "Disabled"}</span>
                </div>
                <small>{plugin.image}</small>
                <div style={{ display: "flex", gap: "0.5rem", flexWrap: "wrap" }}>
                  <button onClick={() => setSelectedPluginId(plugin.plugin_id)}>Select</button>
                  <button onClick={() => void togglePlugin(plugin)}>{plugin.enabled ? "Disable" : "Enable"}</button>
                </div>
              </div>
            ))}
            {plugins.length === 0 ? (
              <div style={{ display: "grid", gap: "0.5rem" }}>
                <p style={{ margin: 0 }}>No plugins registered.</p>
                <button onClick={() => void registerDefaultPlugin()}>Register Default Baseline Plugin</button>
              </div>
            ) : null}
          </div>
        )}

        <div style={{ marginTop: "1rem", borderTop: "1px solid #eee", paddingTop: "0.75rem" }}>
          <p style={{ margin: "0 0 0.5rem" }}>
            <strong>Selected plugin:</strong> {selectedPluginId || "none"}
          </p>
          <div style={{ display: "flex", gap: "0.5rem", flexWrap: "wrap", marginBottom: "0.5rem" }}>
            <button disabled={!selectedPluginId} onClick={() => void runPluginFromManagement()}>
              Run Selected Plugin
            </button>
            <button disabled={!selectedPluginId} onClick={() => void loadPluginRuns(selectedPluginId)}>
              Refresh Runs
            </button>
          </div>
          {pluginRunsError ? <p style={{ color: "#b42318" }}>{pluginRunsError}</p> : null}
          {pluginRunsLoading ? <p>Loading runs...</p> : null}
          {pluginRuns.length > 0 ? (
            <>
              <div style={{ display: "grid", gap: "0.4rem", maxWidth: 620 }}>
                <label>
                  Left run
                  <select value={leftRunId} onChange={(event) => setLeftRunId(event.target.value)}>
                    {pluginRuns.map((run) => (
                      <option key={`left-${run.run_id}`} value={run.run_id}>
                        {run.run_id} ({run.engine})
                      </option>
                    ))}
                  </select>
                </label>
                <label>
                  Right run
                  <select value={rightRunId} onChange={(event) => setRightRunId(event.target.value)}>
                    {pluginRuns.map((run) => (
                      <option key={`right-${run.run_id}`} value={run.run_id}>
                        {run.run_id} ({run.engine})
                      </option>
                    ))}
                  </select>
                </label>
                <button disabled={!leftRunId || !rightRunId || compareLoading} onClick={() => void compareRuns()}>
                  {compareLoading ? "Comparing..." : "Compare Runs"}
                </button>
              </div>
              {compareError ? <p style={{ color: "#b42318" }}>{compareError}</p> : null}
              {compareResult ? (
                <div style={{ marginTop: "0.75rem" }}>
                  <p style={{ margin: "0.25rem 0" }}>
                    <strong>Overlap genes:</strong> {compareResult.overlap_genes.join(", ") || "none"}
                  </p>
                  <p style={{ margin: "0.25rem 0" }}>
                    <strong>Only left:</strong> {compareResult.only_left_genes.join(", ") || "none"}
                  </p>
                  <p style={{ margin: "0.25rem 0" }}>
                    <strong>Only right:</strong> {compareResult.only_right_genes.join(", ") || "none"}
                  </p>
                  <p style={{ margin: "0.5rem 0 0.25rem" }}>
                    <strong>Top score deltas:</strong>
                  </p>
                  <ul style={{ marginTop: 0 }}>
                    {compareResult.score_deltas.slice(0, 5).map((row) => (
                      <li key={`delta-${row.gene}`}>
                        {row.gene}: {row.left_score.toFixed(3)} vs {row.right_score.toFixed(3)} (d {row.delta.toFixed(3)})
                      </li>
                    ))}
                  </ul>
                </div>
              ) : null}
            </>
          ) : (
            <p>No runs for selected plugin yet.</p>
          )}
        </div>
      </section>
      <style jsx>{`
        .analysis-page {
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
        .analysis-page > section {
          z-index: 1;
          grid-column: 1 / -1;
          border: 1px solid var(--line) !important;
          border-radius: 16px !important;
          background: var(--surface) !important;
          padding: 14px !important;
          display: grid !important;
          gap: 10px !important;
          animation: rise 420ms ease both;
        }
        .analysis-page > section:nth-of-type(1) {
          animation-delay: 70ms;
        }
        .analysis-page > section:nth-of-type(2),
        .analysis-page > section:nth-of-type(3) {
          grid-column: span 6;
        }
        .analysis-page > section:nth-of-type(2) {
          animation-delay: 120ms;
        }
        .analysis-page > section:nth-of-type(3) {
          animation-delay: 170ms;
        }
        .analysis-page > section:nth-of-type(4) {
          animation-delay: 220ms;
        }
        .analysis-page > section:nth-of-type(5) {
          animation-delay: 270ms;
        }
        .analysis-page h2 {
          margin: 0 !important;
          font-family: "Space Grotesk", "IBM Plex Sans", sans-serif;
          font-size: 21px !important;
          letter-spacing: 0.01em;
        }
        .analysis-page p {
          color: var(--muted);
        }
        .analysis-page label {
          display: grid;
          gap: 6px;
          color: var(--muted);
          font-size: 12px;
          font-weight: 700;
          letter-spacing: 0.02em;
        }
        .analysis-page input,
        .analysis-page select {
          min-height: 42px;
          border-radius: 10px;
          border: 1px solid #c8d2d8;
          background: #ffffff;
          color: var(--ink);
          padding: 8px 10px;
          font-size: 14px;
          font-family: "IBM Plex Sans", "Noto Sans TC", "Segoe UI", sans-serif;
        }
        .analysis-page input:focus,
        .analysis-page select:focus {
          outline: 2px solid #0f7d76;
          outline-offset: 1px;
          border-color: #0f7d76;
        }
        .analysis-page button {
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
        .analysis-page button:disabled {
          opacity: 0.55;
          cursor: default;
        }
        .analysis-page button:nth-of-type(even) {
          background: var(--deep);
        }
        .analysis-page button:hover {
          border-color: #0c645f;
        }
        .analysis-page a {
          color: #0f6fbd;
          text-decoration: none;
          font-weight: 700;
        }
        .analysis-page a:hover {
          text-decoration: underline;
        }
        .analysis-page ul {
          margin-top: 0;
        }
        .analysis-page *:hover {
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
          .analysis-page > section:nth-of-type(2),
          .analysis-page > section:nth-of-type(3) {
            grid-column: 1 / -1;
          }
          .hero-panel {
            flex-direction: column;
          }
        }
      `}</style>
    </main>
  );
}
