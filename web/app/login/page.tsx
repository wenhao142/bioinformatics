"use client";

import { useMemo, useState } from "react";
import { useRouter } from "next/navigation";
import { KeyRound, ShieldPlus, UserCircle2 } from "lucide-react";

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { resolveApiBase } from "@/lib/api-base";
import { clearStoredToken, TOKEN_STORAGE_KEY } from "@/lib/session";

const DEFAULT_LOGIN_EMAIL = process.env.NEXT_PUBLIC_DEFAULT_LOGIN_EMAIL || "";

export default function LoginPage() {
  const router = useRouter();
  const apiBase = useMemo(() => resolveApiBase(), []);
  const [email, setEmail] = useState(DEFAULT_LOGIN_EMAIL);
  const [password, setPassword] = useState("");
  const [sessionToken, setSessionToken] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [createEmail, setCreateEmail] = useState("");
  const [createPassword, setCreatePassword] = useState("");
  const [createLoading, setCreateLoading] = useState(false);
  const [createError, setCreateError] = useState<string | null>(null);
  const [createMessage, setCreateMessage] = useState<string | null>(null);

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
      if (typeof window !== "undefined") {
        window.localStorage.setItem(TOKEN_STORAGE_KEY, accessToken);
      }
      setSessionToken(accessToken);
      router.replace("/");
    } catch (err: any) {
      setError(err?.message || "Login failed");
    } finally {
      setLoading(false);
    }
  };

  const clearSession = () => {
    if (typeof window !== "undefined") {
      clearStoredToken();
    }
    setSessionToken(null);
    setError(null);
    setCreateError(null);
    setCreateMessage(null);
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
      setSessionToken(adminToken);
      setCreateMessage(`Created account: ${payload.email} (${payload.role})`);
      setCreateEmail("");
      setCreatePassword("");
    } catch (err: any) {
      setCreateError(err?.message || "Create account failed");
    } finally {
      setCreateLoading(false);
    }
  };

  return (
    <main className="min-h-screen bg-[radial-gradient(circle_at_top,_rgba(208,234,228,0.8),_transparent_24%),linear-gradient(180deg,_#fcfcf9_0%,_#f5f3ed_100%)] px-4 py-8 text-foreground sm:px-6 lg:px-8">
      <div className="mx-auto grid min-h-[calc(100vh-4rem)] max-w-6xl gap-6 lg:grid-cols-[1.1fr_0.9fr]">
        <Card className="border-border/80 bg-card/95 backdrop-blur">
          <CardHeader className="space-y-4">
            <Badge variant="secondary" className="w-fit">Intranet sign-in</Badge>
            <CardTitle className="font-serif text-4xl leading-tight">Keep authentication separate from the work.</CardTitle>
            <CardDescription className="max-w-xl text-base">
              Login only unlocks the platform. Data upload, workflow design, and locus review each stay on their own page.
            </CardDescription>
          </CardHeader>
          <CardContent className="grid gap-4 text-sm text-muted-foreground">
            <div className="rounded-xl border border-border bg-background/80 p-4">
              <p className="mb-2 flex items-center gap-2 font-medium text-foreground"><UserCircle2 className="h-4 w-4 text-primary" />After login</p>
              <p>`/` for console and analysis runs</p>
              <p>`/workflow-builder` for raw upload and Snakemake blocks</p>
              <p>`/locus/chr:start-end` for IGV and evidence review</p>
            </div>
            <div className="rounded-xl border border-border bg-background/80 p-4">
              <p className="mb-2 flex items-center gap-2 font-medium text-foreground"><KeyRound className="h-4 w-4 text-primary" />Current admin preset</p>
              <p>{DEFAULT_LOGIN_EMAIL || "No preset admin email configured"}</p>
            </div>
          </CardContent>
        </Card>

        <div className="grid gap-6">
          <Card>
            <CardHeader>
              <CardTitle>Sign in</CardTitle>
              <CardDescription>Use the admin account or your own project account.</CardDescription>
            </CardHeader>
            <CardContent className="grid gap-4">
              <div className="space-y-2">
                <Label htmlFor="email">Email</Label>
                <Input id="email" value={email} onChange={(e) => setEmail(e.target.value)} placeholder="email" />
              </div>
              <div className="space-y-2">
                <Label htmlFor="password">Password</Label>
                <Input
                  id="password"
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
              </div>
              <div className="flex flex-wrap gap-2">
                <Button onClick={() => void login()} disabled={loading}>{loading ? "Signing in" : "Sign in"}</Button>
                <Button variant="outline" onClick={clearSession}>Clear token</Button>
              </div>
              {error ? <p className="text-sm text-destructive">{error}</p> : null}
            </CardContent>
          </Card>

          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2"><ShieldPlus className="h-5 w-5 text-primary" />Create account</CardTitle>
              <CardDescription>Admin only. This stays secondary to sign-in.</CardDescription>
            </CardHeader>
            <CardContent className="grid gap-4">
              <div className="space-y-2">
                <Label htmlFor="create-email">New user email</Label>
                <Input id="create-email" value={createEmail} onChange={(e) => setCreateEmail(e.target.value)} placeholder="new user email" />
              </div>
              <div className="space-y-2">
                <Label htmlFor="create-password">New user password</Label>
                <Input id="create-password" type="password" value={createPassword} onChange={(e) => setCreatePassword(e.target.value)} placeholder="minimum 8 characters" />
              </div>
              <Button variant="outline" onClick={() => void createAccount()} disabled={createLoading}>{createLoading ? "Creating" : "Create account"}</Button>
              {createError ? <p className="text-sm text-destructive">{createError}</p> : null}
              {createMessage ? <p className="text-sm text-emerald-700">{createMessage}</p> : null}
            </CardContent>
          </Card>
        </div>
      </div>
    </main>
  );
}
