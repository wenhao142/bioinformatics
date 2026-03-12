const API_PORT_FALLBACK = 18000;

function isLoopbackHost(hostname: string): boolean {
  return hostname === "localhost" || hostname === "127.0.0.1" || hostname === "0.0.0.0";
}

export function resolveApiBase(): string {
  const envUrl = process.env.NEXT_PUBLIC_API_URL?.trim();

  if (typeof window !== "undefined") {
    const fallback = `${window.location.protocol}//${window.location.hostname}:${API_PORT_FALLBACK}`;

    if (!envUrl) {
      return fallback;
    }

    try {
      const parsed = new URL(envUrl);
      if (isLoopbackHost(parsed.hostname) && !isLoopbackHost(window.location.hostname)) {
        parsed.protocol = window.location.protocol;
        parsed.hostname = window.location.hostname;
        if (!parsed.port) {
          parsed.port = String(API_PORT_FALLBACK);
        }
      }
      return parsed.toString().replace(/\/$/, "");
    } catch {
      return fallback;
    }
  }

  return envUrl || `http://localhost:${API_PORT_FALLBACK}`;
}
