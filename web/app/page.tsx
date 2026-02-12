"use client";

import { useEffect, useState } from "react";

type Health = { status: string };

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

  useEffect(() => {
    const base = resolveApiBase();
    fetch(`${base}/health`)
      .then((r) => r.json())
      .then((data: Health) => setApiStatus(data.status))
      .catch(() => setApiStatus("offline"));
  }, []);

  return (
    <main style={{ padding: "2rem", fontFamily: "system-ui, -apple-system" }}>
      <h1>AD Multi-Omics Locus Evidence Platform</h1>
      <p>Scaffold ready. Web + API + worker + storage via docker compose.</p>
      <div style={{ marginTop: "1rem", padding: "1rem", border: "1px solid #ddd", borderRadius: 8 }}>
        <strong>API health:</strong> {apiStatus}
      </div>
    </main>
  );
}
