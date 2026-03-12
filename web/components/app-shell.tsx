import Link from "next/link";
import { ReactNode } from "react";

import { Badge } from "@/components/ui/badge";
import { Button, buttonVariants } from "@/components/ui/button";
import { cn } from "@/lib/utils";

type NavItem = {
  label: string;
  href?: string;
  active?: boolean;
  onClick?: () => void;
  variant?: "default" | "secondary" | "outline" | "ghost" | "destructive";
};

export function AppShell({
  title,
  subtitle,
  badge,
  navItems,
  children,
}: {
  title: string;
  subtitle?: string;
  badge?: string;
  navItems?: NavItem[];
  children: ReactNode;
}) {
  return (
    <main className="min-h-screen bg-[radial-gradient(circle_at_top_left,_rgba(208,234,228,0.75),_transparent_28%),linear-gradient(180deg,_#fcfcf9_0%,_#f5f3ed_100%)] text-foreground">
      <div className="mx-auto flex min-h-screen w-full max-w-7xl flex-col gap-8 px-4 py-6 sm:px-6 lg:px-8">
        <header className="rounded-2xl border border-border/80 bg-card/90 p-5 shadow-sm backdrop-blur">
          <div className="flex flex-col gap-4 md:flex-row md:items-center md:justify-between">
            <div className="space-y-2">
              <p className="text-xs font-semibold uppercase tracking-[0.28em] text-muted-foreground">AD Multi-Omics Platform</p>
              <div className="flex flex-wrap items-center gap-3">
                <h1 className="font-serif text-3xl tracking-tight text-foreground md:text-4xl">{title}</h1>
                {badge ? <Badge variant="secondary">{badge}</Badge> : null}
              </div>
              {subtitle ? <p className="max-w-2xl text-sm text-muted-foreground md:text-base">{subtitle}</p> : null}
            </div>
            {navItems?.length ? (
              <nav className="flex flex-wrap gap-2">
                {navItems.map((item) => {
                  const variant = item.active ? "default" : item.variant || "outline";
                  if (item.href) {
                    return (
                      <Link key={`${item.label}:${item.href}`} href={item.href} className={cn(item.active ? "pointer-events-none" : "") }>
                        <span className={buttonVariants({ variant, size: "sm" })}>{item.label}</span>
                      </Link>
                    );
                  }
                  return (
                    <Button key={item.label} variant={variant} size="sm" onClick={item.onClick}>
                      {item.label}
                    </Button>
                  );
                })}
              </nav>
            ) : null}
          </div>
        </header>
        {children}
      </div>
    </main>
  );
}
