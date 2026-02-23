"use client";

import { useEffect, useMemo, useState } from "react";

type Health = { status: string };
type CausalRunResult = {
  run: { run_id: string; kind: string; created_by: string; project_id: string | null };
  result: { gene_scores: Array<{ gene: string; score: number; rank: number }> };
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

export default function Home() {
  const [apiStatus, setApiStatus] = useState<string>("checking...");
  const [token, setToken] = useState<string>("");
  const [region, setRegion] = useState<string>("chr1:1-1000");
  const [projectId, setProjectId] = useState<string>("demo-project");
  const [causalResult, setCausalResult] = useState<CausalRunResult | null>(null);
  const [causalLoading, setCausalLoading] = useState(false);
  const [causalError, setCausalError] = useState<string | null>(null);
  const apiBase = useMemo(() => resolveApiBase(), []);

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

  const runCausal = async () => {
    const match = region.match(/^([^:]+):(\d+)-(\d+)$/);
    if (!match) {
      setCausalError("Region format must be chr:start-end");
      return;
    }
    if (!token) {
      setCausalError("API token required");
      return;
    }
    const [, chr, start, end] = match;
    setCausalLoading(true);
    setCausalError(null);
    setCausalResult(null);
    try {
      const response = await fetch(
        `${apiBase}/causal/score?chr=${encodeURIComponent(chr)}&start=${start}&end=${end}&top_n=5&project_id=${encodeURIComponent(projectId)}`,
        { method: "POST", headers: { Authorization: `Bearer ${token}` } },
      );
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Request failed (${response.status})`);
      }
      setCausalResult(payload as CausalRunResult);
    } catch (error: any) {
      setCausalError(error?.message || "Failed to run causal scoring");
    } finally {
      setCausalLoading(false);
    }
  };

  return (
    <main style={{ padding: "2rem", fontFamily: "system-ui, -apple-system" }}>
      <h1>AD Multi-Omics Locus Evidence Platform</h1>
      <p>Scaffold ready. Web + API + worker + storage via docker compose.</p>
      <div style={{ marginTop: "1rem", padding: "1rem", border: "1px solid #ddd", borderRadius: 8 }}>
        <strong>API health:</strong> {apiStatus}
      </div>
      <section style={{ marginTop: "1rem", padding: "1rem", border: "1px solid #ddd", borderRadius: 8 }}>
        <h2 style={{ marginTop: 0 }}>Causal Scoring (T3.4 MVP)</h2>
        <p style={{ marginTop: 0 }}>Run per-variant and per-gene causal scoring from uploaded VCF (+ optional omics).</p>
        <div style={{ display: "grid", gap: "0.5rem", maxWidth: 620 }}>
          <input value={token} onChange={(event) => setToken(event.target.value)} placeholder="JWT token" />
          <input value={region} onChange={(event) => setRegion(event.target.value)} placeholder="chr:start-end" />
          <input value={projectId} onChange={(event) => setProjectId(event.target.value)} placeholder="project id" />
          <button onClick={runCausal} disabled={causalLoading}>
            {causalLoading ? "Running..." : "Run Causal Scoring"}
          </button>
        </div>
        {causalError ? <p style={{ color: "#b42318" }}>{causalError}</p> : null}
        {causalResult ? (
          <div style={{ marginTop: "0.75rem" }}>
            <p style={{ margin: 0 }}>
              <strong>Run ID:</strong> {causalResult.run.run_id}
            </p>
            <p style={{ margin: "0.35rem 0" }}>
              <strong>Top Genes:</strong>
            </p>
            <ul style={{ marginTop: 0 }}>
              {causalResult.result.gene_scores.slice(0, 5).map((row) => (
                <li key={`${row.rank}-${row.gene}`}>
                  #{row.rank} {row.gene} ({row.score.toFixed(3)})
                </li>
              ))}
            </ul>
          </div>
        ) : null}
      </section>
    </main>
  );
}
