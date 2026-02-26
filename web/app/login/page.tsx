"use client";

import { useEffect, useMemo, useState } from "react";
import { useRouter } from "next/navigation";

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

export default function LoginPage() {
  const router = useRouter();
  const apiBase = useMemo(() => resolveApiBase(), []);
  const [email, setEmail] = useState("admin@example.com");
  const [password, setPassword] = useState("password");
  const [sessionToken, setSessionToken] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [createEmail, setCreateEmail] = useState("");
  const [createPassword, setCreatePassword] = useState("");
  const [createLoading, setCreateLoading] = useState(false);
  const [createError, setCreateError] = useState<string | null>(null);
  const [createMessage, setCreateMessage] = useState<string | null>(null);
  const [ready, setReady] = useState(false);

  useEffect(() => {
    if (typeof window === "undefined") return;
    const existing = window.localStorage.getItem("ad_api_token");
    if (existing) {
      setSessionToken(existing);
    }
    setReady(true);
  }, []);

  const fetchLoginToken = async (adminEmail: string, adminPassword: string): Promise<string> => {
    const response = await fetch(`${apiBase}/auth/login`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ email: adminEmail.trim(), password: adminPassword }),
    });
    const payload = await response.json().catch(() => ({}));
    if (!response.ok) {
      throw new Error(typeof payload?.detail === "string" ? payload.detail : `Login failed (${response.status})`);
    }
    const accessToken = typeof payload?.access_token === "string" ? payload.access_token : "";
    if (!accessToken) {
      throw new Error("Login response missing token");
    }
    return accessToken;
  };

  const login = async () => {
    setLoading(true);
    setError(null);
    try {
      const accessToken = await fetchLoginToken(email, password);
      window.localStorage.setItem("ad_api_token", accessToken);
      setSessionToken(accessToken);
      router.replace("/");
    } catch (err: any) {
      setError(err?.message || "Login failed");
    } finally {
      setLoading(false);
    }
  };

  const switchAccount = () => {
    if (typeof window !== "undefined") {
      window.localStorage.removeItem("ad_api_token");
    }
    setSessionToken(null);
    setCreateMessage(null);
    setCreateError(null);
    setError(null);
  };

  const createAccount = async () => {
    setCreateLoading(true);
    setCreateError(null);
    setCreateMessage(null);
    try {
      const targetEmail = createEmail.trim();
      if (!targetEmail) {
        throw new Error("New user email is required");
      }
      if (createPassword.length < 8) {
        throw new Error("New user password must be at least 8 characters");
      }

      const adminToken = sessionToken ?? (await fetchLoginToken(email, password));
      const response = await fetch(`${apiBase}/auth/register`, {
        method: "POST",
        headers: {
          Authorization: `Bearer ${adminToken}`,
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ email: targetEmail, password: createPassword }),
      });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(typeof payload?.detail === "string" ? payload.detail : `Create account failed (${response.status})`);
      }
      setCreateMessage(`Created account: ${payload.email} (${payload.role})`);
      setCreateEmail("");
      setCreatePassword("");
    } catch (err: any) {
      setCreateError(err?.message || "Create account failed");
    } finally {
      setCreateLoading(false);
    }
  };

  if (!ready) {
    return (
      <main className="login-page">
        <section className="login-shell">
          <h1>Checking session...</h1>
        </section>
      </main>
    );
  }

  return (
    <main className="login-page">
      <div className="ambient-shape shape-a" />
      <div className="ambient-shape shape-b" />

      <section className="login-shell reveal">
        <p className="hero-kicker">AD Multi-Omics Locus Evidence Platform</p>
        <h1 className="hero-title">Sign In</h1>
        <p className="hero-note">Authenticate first, then enter the Analysis Console.</p>

        {sessionToken ? (
          <div className="token-box">
            <p className="ok-text">An active session already exists on this browser.</p>
            <div className="btn-row">
              <button type="button" onClick={() => router.push("/")}>
                Enter Analysis Console
              </button>
              <button type="button" onClick={switchAccount}>
                Switch Account
              </button>
            </div>
          </div>
        ) : null}

        <label>
          Email
          <input value={email} onChange={(e) => setEmail(e.target.value)} placeholder="email" />
        </label>
        <label>
          Password
          <input
            type="password"
            value={password}
            onChange={(e) => setPassword(e.target.value)}
            placeholder="password"
            onKeyDown={(e) => {
              if (e.key === "Enter" && !loading) {
                void login();
              }
            }}
          />
        </label>
        <button onClick={() => void login()} disabled={loading}>
          {loading ? "Signing In..." : "Sign In"}
        </button>
        {error ? <p className="err-text">{error}</p> : null}

        <div className="divider" />
        <h2 className="sub-title">Create Account</h2>
        <p className="hero-note">Use admin credentials above, then create a new user.</p>
        <label>
          New User Email
          <input value={createEmail} onChange={(e) => setCreateEmail(e.target.value)} placeholder="new user email" />
        </label>
        <label>
          New User Password
          <input
            type="password"
            value={createPassword}
            onChange={(e) => setCreatePassword(e.target.value)}
            placeholder="new user password"
          />
        </label>
        <button type="button" onClick={() => void createAccount()} disabled={createLoading}>
          {createLoading ? "Creating..." : "Create Account"}
        </button>
        {createError ? <p className="err-text">{createError}</p> : null}
        {createMessage ? <p className="ok-text">{createMessage}</p> : null}
      </section>

      <style jsx>{`
        .login-page {
          --bg: #edf2f4;
          --surface: #ffffff;
          --line: #ccd6dc;
          --ink: #10222c;
          --muted: #4d626e;
          --teal: #0f7d76;
          --orange: #c8611f;
          min-height: 100vh;
          display: grid;
          place-items: center;
          padding: 16px;
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
          bottom: -100px;
          right: -90px;
          background: #c8611f;
        }
        .login-shell {
          z-index: 1;
          width: min(460px, 100%);
          border: 1px solid var(--line);
          border-radius: 18px;
          padding: 18px;
          background: var(--surface);
          display: grid;
          gap: 10px;
        }
        .hero-kicker {
          margin: 0;
          letter-spacing: 0.14em;
          text-transform: uppercase;
          font-size: 11px;
          font-weight: 700;
          color: var(--teal);
        }
        .hero-title {
          margin: 0;
          font-family: "Space Grotesk", "IBM Plex Sans", sans-serif;
          font-size: clamp(30px, 4.6vw, 42px);
          line-height: 1.03;
          letter-spacing: 0.01em;
        }
        .hero-note {
          margin: 0;
          color: var(--muted);
          font-size: 13px;
        }
        .sub-title {
          margin: 0;
          font-family: "Space Grotesk", "IBM Plex Sans", sans-serif;
          font-size: 22px;
          letter-spacing: 0.01em;
        }
        .token-box {
          border: 1px solid #d6dfe4;
          border-radius: 10px;
          padding: 10px;
          display: grid;
          gap: 8px;
          background: #f7fafb;
        }
        .btn-row {
          display: flex;
          flex-wrap: wrap;
          gap: 8px;
        }
        .divider {
          height: 1px;
          background: #d6dfe4;
          margin: 6px 0;
        }
        label {
          display: grid;
          gap: 6px;
          color: var(--muted);
          font-size: 12px;
          font-weight: 700;
          letter-spacing: 0.02em;
        }
        input {
          min-height: 42px;
          border-radius: 10px;
          border: 1px solid #c8d2d8;
          background: #ffffff;
          color: var(--ink);
          padding: 8px 10px;
          font-size: 14px;
          font-family: "IBM Plex Sans", "Noto Sans TC", "Segoe UI", sans-serif;
        }
        input:focus {
          outline: 2px solid #0f7d76;
          outline-offset: 1px;
          border-color: #0f7d76;
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
        }
        button:disabled {
          opacity: 0.55;
          cursor: default;
        }
        .err-text {
          margin: 0;
          color: #b42318;
          font-weight: 700;
        }
        .ok-text {
          margin: 0;
          color: #14532d;
          font-weight: 700;
        }
        .reveal {
          animation: rise 380ms ease both;
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
      `}</style>
    </main>
  );
}
