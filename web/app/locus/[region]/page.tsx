'use client';

import { useEffect, useRef, useState } from 'react';

const parseRegion = (region: string) => {
  const [chrPart, range] = region.split(':');
  if (!range) return { chr: chrPart, start: 1, end: 100000 };
  const [startStr, endStr] = range.split('-');
  const start = parseInt(startStr, 10) || 1;
  const end = parseInt(endStr, 10) || start + 1000;
  return { chr: chrPart, start, end };
};

export default function LocusPage({ params }: { params: { region: string } }) {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let browser: any;
    const load = async () => {
      try {
        const mod = await import('igv/dist/igv.esm'); // ESM build exposes createBrowser
        const igv = (mod as any).default ?? mod;
        const decoded = decodeURIComponent(params.region);
        const { chr, start, end } = parseRegion(decoded);
        browser = await igv.createBrowser(containerRef.current as HTMLDivElement, {
          genome: 'hg38',
          locus: `${chr}:${start}-${end}`,
          tracks: [
            {
              name: 'Genes (sample)',
              type: 'annotation',
              format: 'bed',
              url: '/genes.sample.bed', // local static file to avoid S3 CORS/forbidden
            },
          ],
        });
      } catch (e: any) {
        setError(e?.message || 'Failed to load genome viewer');
      }
    };
    load();
    return () => {
      if (browser) {
        browser.dispose();
      }
    };
  }, [params.region]);

  return (
    <main className="min-h-screen p-6 space-y-4">
      <h1 className="text-2xl font-semibold">Locus Explorer</h1>
      <p className="text-sm text-gray-600">Region: {params.region}</p>
      {error ? (
        <div className="text-red-600 text-sm">{error}</div>
      ) : (
        <div
          ref={containerRef}
          className="border rounded shadow-sm"
          style={{ minHeight: 500 }}
        />
      )}
    </main>
  );
}
