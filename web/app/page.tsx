"use client";

import { useEffect, useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import { ArrowUpRight, Dna, FlaskConical, LogOut, MapPinned, Play, Upload, Workflow } from "lucide-react";

import { AppShell } from "@/components/app-shell";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Select } from "@/components/ui/select";
import { Separator } from "@/components/ui/separator";
import { resolveApiBase } from "@/lib/api-base";
import { clearStoredToken, readStoredToken, TOKEN_STORAGE_KEY } from "@/lib/session";


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

type MethodRunResult = {
  kind: "evidence" | "causal";
  runId: string;
  selectedMethod: string;
  topGenes: GeneScore[];
  reportUrl: string | null;
};

function parseRegion(region: string): { chr: string; start: number; end: number } {
  const match = region.match(/^([^:]+):(\d+)-(\d+)$/);
  if (!match) {
    throw new Error("Region must be chr:start-end");
  }
  const [, chr, start, end] = match;
  return { chr, start: Number(start), end: Number(end) };
}

function scoreText(value: number): string {
  return Number.isFinite(value) ? value.toFixed(3) : String(value);
}

export default function Home() {
  const router = useRouter();
  const apiBase = useMemo(() => resolveApiBase(), []);
  const [apiStatus, setApiStatus] = useState<string>("checking");
  const [token, setToken] = useState<string>("");
  const [authReady, setAuthReady] = useState(false);
  const [sessionChecked, setSessionChecked] = useState(false);
  const [projectId, setProjectId] = useState<string>("demo-project");
  const [region, setRegion] = useState<string>("chr1:1-1000");

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
    const stored = readStoredToken();
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
        setToken("");
        setAuthReady(true);
        setAnalysisError("Session expired. Please sign in again.");
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

  const logoutSession = () => {
    setToken("");
    clearStoredToken();
    router.replace("/login");
  };

  const handleUnauthorized = (message = "Session expired. Please sign in again.") => {
    clearStoredToken();
    setToken("");
    setUploadError(message);
    setAnalysisError(message);
    router.replace("/login");
  };

  const openLocus = () => {
    try {
      parseRegion(region);
      router.push(`/locus/${encodeURIComponent(region)}`);
    } catch (error: any) {
      setAnalysisError(error?.message || "Invalid region");
    }
  };

  const uploadFile = async (kind: "variants" | "expr" | "raw") => {
    if (!sessionChecked) {
      setUploadError("Checking session. Try again in a moment.");
      return;
    }
    if (!authHeaders) {
      setUploadError("Login required");
      return;
    }
    const file = kind === "variants" ? variantFile : kind === "expr" ? exprFile : rawDataFile;
    if (!file) {
      setUploadError(kind === "variants" ? "Choose a VCF file." : kind === "expr" ? "Choose an expression TSV." : "Choose a raw bio file.");
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
        if (response.status === 401) {
          handleUnauthorized(typeof payload?.detail === "string" ? payload.detail : "Invalid token");
          return;
        }
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Upload failed (${response.status})`);
      }
      if (kind === "variants") {
        setUploadMessage(`VCF uploaded: ${payload.ingested ?? 0} variants`);
      } else if (kind === "expr") {
        setUploadMessage(`Expression uploaded: ${payload.ingested ?? 0} rows`);
      } else {
        setUploadMessage(`Raw file uploaded: ${payload?.dataset?.filename ?? file.name}`);
      }
    } catch (error: any) {
      setUploadError(error?.message || "Upload failed");
    } finally {
      setUploadLoading(false);
    }
  };

  const runSelectedMethod = async () => {
    if (!sessionChecked) {
      setAnalysisError("Checking session. Try again in a moment.");
      return;
    }
    if (!authHeaders) {
      setAnalysisError("Login required");
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
        const payload = (await response.json().catch(() => ({}))) as Partial<EvidenceRankResponse> & { detail?: string };
        if (!response.ok || !payload.run || !payload.result) {
          if (response.status === 401) {
            handleUnauthorized(typeof payload.detail === "string" ? payload.detail : "Invalid token");
            return;
          }
          throw new Error(typeof payload.detail === "string" ? payload.detail : `Request failed (${response.status})`);
        }
        setMethodRunResult({
          kind: "evidence",
          runId: payload.run.run_id,
          selectedMethod: payload.run.selected_method || analysisMethod,
          topGenes: payload.result.ranked_genes.slice(0, 5).map((row) => ({ gene: row.gene, score: row.score, rank: row.rank })),
          reportUrl: `${apiBase}/report/${encodeURIComponent(projectId)}/${encodeURIComponent(payload.run.run_id)}`,
        });
        return;
      }
      const response = await fetch(
        `${apiBase}/causal/score?chr=${encodeURIComponent(chr)}&start=${start}&end=${end}&top_n=5&project_id=${encodeURIComponent(projectId)}&method=${encodeURIComponent(analysisMethod)}`,
        { method: "POST", headers: authHeaders },
      );
      const payload = (await response.json().catch(() => ({}))) as Partial<CausalRunResponse> & { detail?: string };
      if (!response.ok || !payload.run || !payload.result) {
        if (response.status === 401) {
          handleUnauthorized(typeof payload.detail === "string" ? payload.detail : "Invalid token");
          return;
        }
        throw new Error(typeof payload.detail === "string" ? payload.detail : `Request failed (${response.status})`);
      }
      setMethodRunResult({
        kind: "causal",
        runId: payload.run.run_id,
        selectedMethod: payload.run.selected_method || analysisMethod,
        topGenes: payload.result.gene_scores.slice(0, 5),
        reportUrl: `${apiBase}/report/${encodeURIComponent(projectId)}/${encodeURIComponent(payload.run.run_id)}`,
      });
    } catch (error: any) {
      setAnalysisError(error?.message || "Run failed");
    } finally {
      setAnalysisLoading(false);
    }
  };

  const openReportWithAuth = async (format: "markdown" | "html") => {
    if (!sessionChecked) {
      setReportError("Checking session. Try again in a moment.");
      return;
    }
    if (!authHeaders || !methodRunResult?.reportUrl || !token) {
      setReportError("No report yet.");
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
        if (response.status === 401) {
          handleUnauthorized(typeof body?.detail === "string" ? body.detail : "Invalid token");
          return;
        }
        throw new Error(typeof body?.detail === "string" ? body.detail : `Report request failed (${response.status})`);
      }
      const content = await response.text();
      const mimeType = format === "html" ? "text/html;charset=utf-8" : "text/markdown;charset=utf-8";
      const blob = new Blob([content], { type: mimeType });
      const objectUrl = URL.createObjectURL(blob);
      const popup = window.open(objectUrl, "_blank", "noopener,noreferrer");
      if (!popup) {
        throw new Error("Popup blocked.");
      }
      setTimeout(() => URL.revokeObjectURL(objectUrl), 60_000);
    } catch (error: any) {
      setReportError(error?.message || "Failed to open report");
    }
  };

  if (!authReady || !token) {
    return (
      <div className="flex min-h-screen items-center justify-center bg-background px-4">
        <Card className="w-full max-w-md">
          <CardHeader>
            <CardTitle>Checking session</CardTitle>
            <CardDescription>Redirecting to login.</CardDescription>
          </CardHeader>
        </Card>
      </div>
    );
  }

  return (
    <AppShell
      title="Console"
      subtitle="Upload data, run scoring, then jump into the workflow builder or locus explorer."
      badge={apiStatus}
      navItems={[
        { label: "Console", href: "/", active: true },
        { label: "Workflow Builder", href: "/workflow-builder" },
        { label: "Locus Explorer", href: `/locus/${encodeURIComponent(region)}` },
        { label: "Logout", onClick: logoutSession, variant: "ghost" },
      ]}
    >
      <section className="grid gap-4 md:grid-cols-3">
        <Card>
          <CardHeader>
            <CardDescription>Project</CardDescription>
            <CardTitle>{projectId}</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-muted-foreground">Use one project id across uploads, workflows, and reports.</p>
          </CardContent>
        </Card>
        <Card>
          <CardHeader>
            <CardDescription>Region</CardDescription>
            <CardTitle>{region}</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-muted-foreground">Current locus used for scoring and locus explorer deep link.</p>
          </CardContent>
        </Card>
        <Card>
          <CardHeader>
            <CardDescription>Latest run</CardDescription>
            <CardTitle>{methodRunResult?.runId || "Not started"}</CardTitle>
          </CardHeader>
          <CardContent>
            <p className="text-sm text-muted-foreground">{methodRunResult ? `${methodRunResult.topGenes.length} top genes ready` : "Run a method to populate results."}</p>
          </CardContent>
        </Card>
      </section>

      <section className="grid gap-6 xl:grid-cols-[1.2fr_1fr]">
        <Card>
          <CardHeader>
            <CardTitle>Quick start</CardTitle>
            <CardDescription>Set project scope once, then use the same values across pages.</CardDescription>
          </CardHeader>
          <CardContent className="grid gap-4 md:grid-cols-2">
            <div className="space-y-2">
              <Label htmlFor="project-id">Project ID</Label>
              <Input id="project-id" value={projectId} onChange={(e) => setProjectId(e.target.value)} />
            </div>
            <div className="space-y-2">
              <Label htmlFor="region">Region</Label>
              <Input id="region" value={region} onChange={(e) => setRegion(e.target.value)} />
            </div>
            <div className="md:col-span-2 flex flex-wrap gap-2">
              <Button onClick={() => router.push("/workflow-builder")}><Workflow className="h-4 w-4" />Open workflow builder</Button>
              <Button variant="outline" onClick={openLocus}><MapPinned className="h-4 w-4" />Open locus explorer</Button>
              <Button variant="ghost" onClick={() => router.push("/workflow-builder")}><ArrowUpRight className="h-4 w-4" />Go to build and run</Button>
            </div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>Run scoring</CardTitle>
            <CardDescription>Use uploaded project data and the current locus window.</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="space-y-2">
              <Label htmlFor="method">Method</Label>
              <Select id="method" value={analysisMethod} onChange={(e) => setAnalysisMethod(e.target.value)}>
                <option value="causal-score">Causal score</option>
                <option value="evidence-rank">Evidence rank</option>
              </Select>
            </div>
            <div className="flex flex-wrap gap-2">
              <Button onClick={() => void runSelectedMethod()} disabled={!sessionChecked || analysisLoading}><Play className="h-4 w-4" />{analysisLoading ? "Running" : "Run"}</Button>
              <Button variant="outline" onClick={() => void openReportWithAuth("markdown")} disabled={!sessionChecked || !methodRunResult}>Markdown</Button>
              <Button variant="outline" onClick={() => void openReportWithAuth("html")} disabled={!sessionChecked || !methodRunResult}>HTML</Button>
            </div>
            {analysisError ? <p className="text-sm text-destructive">{analysisError}</p> : null}
            {reportError ? <p className="text-sm text-destructive">{reportError}</p> : null}
          </CardContent>
        </Card>
      </section>

      <section className="grid gap-6 lg:grid-cols-3">
        {[
          { key: "variants", title: "Upload VCF", icon: Dna, accept: ".vcf,.vcf.gz,.tsv,.txt", file: variantFile, setFile: setVariantFile, run: () => uploadFile("variants") },
          { key: "expr", title: "Upload Expression", icon: FlaskConical, accept: ".tsv,.csv,.txt", file: exprFile, setFile: setExprFile, run: () => uploadFile("expr") },
          { key: "raw", title: "Upload Raw", icon: Upload, accept: ".fasta,.fa,.fna,.gtf,.gff,.gff3,.fastq,.fq,.fastq.gz,.fq.gz,.bam,.cram,.sam,.tsv,.csv,.txt,.gz", file: rawDataFile, setFile: setRawDataFile, run: () => uploadFile("raw") },
        ].map((item) => (
          <Card key={item.key}>
            <CardHeader>
              <CardTitle className="flex items-center gap-2 text-lg"><item.icon className="h-5 w-5 text-primary" />{item.title}</CardTitle>
              <CardDescription>{item.key === "raw" ? "For builder and pipeline inputs." : "Stored into the current project scope."}</CardDescription>
            </CardHeader>
            <CardContent className="space-y-3">
              <Input type="file" accept={item.accept} onChange={(e) => item.setFile(e.target.files?.[0] || null)} />
              <p className="text-sm text-muted-foreground">{item.file ? item.file.name : "No file selected"}</p>
              <Button variant="outline" onClick={() => void item.run()} disabled={!sessionChecked || uploadLoading}>{uploadLoading ? "Uploading" : item.title}</Button>
            </CardContent>
          </Card>
        ))}
      </section>

      <section className="grid gap-6 xl:grid-cols-[1.15fr_0.85fr]">
        <Card>
          <CardHeader>
            <CardTitle>Latest result</CardTitle>
            <CardDescription>Keep this page minimal. Detailed run traces live in builder and locus pages.</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {uploadError ? <p className="text-sm text-destructive">{uploadError}</p> : null}
            {uploadMessage ? <p className="text-sm text-emerald-700">{uploadMessage}</p> : null}
            {methodRunResult ? (
              <div className="space-y-4">
                <div className="flex flex-wrap items-center gap-2">
                  <Badge variant="secondary">{methodRunResult.selectedMethod}</Badge>
                  <Badge variant="outline">{methodRunResult.runId}</Badge>
                </div>
                <Separator />
                <div className="space-y-3">
                  {methodRunResult.topGenes.map((gene) => (
                    <div key={`${gene.rank}:${gene.gene}`} className="flex items-center justify-between rounded-lg border border-border bg-muted/50 px-4 py-3">
                      <div>
                        <p className="font-medium">#{gene.rank} {gene.gene}</p>
                        <p className="text-sm text-muted-foreground">Gene score</p>
                      </div>
                      <p className="font-mono text-sm">{scoreText(gene.score)}</p>
                    </div>
                  ))}
                </div>
              </div>
            ) : (
              <p className="text-sm text-muted-foreground">No result yet. Upload data and run a method or move to workflow builder for raw-first testing.</p>
            )}
          </CardContent>
        </Card>

        <Card>
          <CardHeader>
            <CardTitle>What to do next</CardTitle>
            <CardDescription>Use the shortest path based on what you are testing.</CardDescription>
          </CardHeader>
          <CardContent className="space-y-3 text-sm text-muted-foreground">
            <div className="rounded-lg border border-border bg-background p-4">
              <p className="font-medium text-foreground">Raw workflow test</p>
              <p>Go to workflow builder, upload raw files, load WGS template, save, then run.</p>
            </div>
            <div className="rounded-lg border border-border bg-background p-4">
              <p className="font-medium text-foreground">Locus review</p>
              <p>Keep the region value here, then open locus explorer to inspect tracks and causal overlays.</p>
            </div>
          </CardContent>
        </Card>
      </section>
    </AppShell>
  );
}
