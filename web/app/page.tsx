"use client";

import { useEffect, useMemo, useState } from "react";
import { useRouter } from "next/navigation";

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
  kind: "evidence" | "causal";
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
  const router = useRouter();
  const [apiStatus, setApiStatus] = useState<string>("checking...");
  const [token, setToken] = useState<string>("");
  const [authReady, setAuthReady] = useState(false);
  const [region, setRegion] = useState<string>("chr1:1-1000");
  const [projectId, setProjectId] = useState<string>("demo-project");

  const [variantFile, setVariantFile] = useState<File | null>(null);
  const [exprFile, setExprFile] = useState<File | null>(null);
  const [rawDataFile, setRawDataFile] = useState<File | null>(null);
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
    if (stored) {
      setToken(stored);
      setAuthReady(true);
      return;
    }
    setAuthReady(true);
    router.replace("/login");
  }, [router]);

  useEffect(() => {
    if (typeof window === "undefined" || !authReady) return;
    if (!token) {
      window.localStorage.removeItem("ad_api_token");
      return;
    }
    window.localStorage.setItem("ad_api_token", token);
  }, [token, authReady]);

  const logoutSession = () => {
    setToken("");
    router.replace("/login");
  };

  const analysisMethodOptions = useMemo(
    () => [
      { value: "evidence-rank", label: "Baseline Evidence Rank" },
      { value: "causal-score", label: "Causal Score" },
    ],
    [],
  );

  useEffect(() => {
    const exists = analysisMethodOptions.some((item) => item.value === analysisMethod);
    if (!exists && analysisMethodOptions.length > 0) {
      setAnalysisMethod(analysisMethodOptions[0].value);
    }
  }, [analysisMethod, analysisMethodOptions]);

  const uploadFile = async (kind: "variants" | "expr" | "raw") => {
    if (!authHeaders) {
      setUploadError("API token required");
      return;
    }
    const file = kind === "variants" ? variantFile : kind === "expr" ? exprFile : rawDataFile;
    if (!file) {
      if (kind === "variants") {
        setUploadError("Please choose a VCF file.");
      } else if (kind === "expr") {
        setUploadError("Please choose an expression TSV file.");
      } else {
        setUploadError("Please choose a raw bioinformatics file.");
      }
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
          : kind === "expr"
            ? `${apiBase}/omics/expr/upload`
            : `${apiBase}/datasets/upload?project_id=${encodeURIComponent(projectId)}`;
      const response = await fetch(endpoint, { method: "POST", headers: authHeaders, body: form });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Upload failed (${response.status})`);
      }
      if (kind === "variants") {
        setUploadMessage(`VCF uploaded. Ingested ${payload.ingested ?? 0} variants.`);
      } else if (kind === "expr") {
        setUploadMessage(`Expression table uploaded. Ingested ${payload.ingested ?? 0} rows.`);
      } else {
        const uploadedName = payload?.dataset?.filename ?? file.name;
        setUploadMessage(`Raw dataset uploaded: ${uploadedName}`);
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

  if (!authReady || !token) {
    return (
      <main className="analysis-page">
        <section style={{ gridColumn: "1 / -1", padding: "1rem", border: "1px solid #ddd", borderRadius: 8 }}>
          <h2 style={{ margin: 0 }}>Checking session...</h2>
          <p style={{ margin: "0.5rem 0 0" }}>Redirecting to login page.</p>
        </section>
      </main>
    );
  }

  return (
    <main className="analysis-page">
      <div className="ambient-shape shape-a" />
      <div className="ambient-shape shape-b" />
      <header className="hero-panel reveal reveal-1">
        <div>
          <p className="hero-kicker">AD Multi-Omics Locus Evidence Platform</p>
          <h1 className="hero-title">Analysis Console</h1>
          <p className="hero-note">Manage uploads, run statistical methods, and control LLM summary mode in one page.</p>
          <p className="hero-note" style={{ marginTop: 6 }}>
            <a href="/workflow-builder">Open Workflow Builder</a>
          </p>
        </div>
        <div className="hero-actions">
          <button onClick={logoutSession} type="button">
            Logout
          </button>
          <div className={`health-pill ${apiStatus === "ok" ? "health-pill-ok" : "health-pill-warn"}`}>API health: {apiStatus}</div>
        </div>
      </header>

      <section style={{ padding: "1rem", border: "1px solid #ddd", borderRadius: 8, display: "grid", gap: "0.6rem" }}>
        <h2 style={{ margin: 0 }}>Session</h2>
        <p style={{ margin: 0, fontSize: 12, color: "#4d626e" }}>
          Logged in. Token is stored in browser localStorage.
        </p>
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
          <label>
            Raw Bio File (FASTA/GTF/FASTQ/BAM/etc.)
            <input
              type="file"
              accept=".fasta,.fa,.fna,.gtf,.gff,.gff3,.fastq,.fq,.fastq.gz,.fq.gz,.bam,.cram,.sam,.tsv,.csv,text/plain,application/gzip,application/octet-stream"
              onChange={(e) => setRawDataFile(e.target.files?.[0] || null)}
            />
          </label>
          <button disabled={uploadLoading || !rawDataFile} onClick={() => void uploadFile("raw")}>
            {uploadLoading ? "Uploading..." : "Upload Raw Dataset"}
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
        .hero-actions {
          display: grid;
          gap: 8px;
          justify-items: end;
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
