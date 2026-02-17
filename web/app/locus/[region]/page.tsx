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
  const [selected, setSelected] = useState<{ track?: string; name?: string; pos?: string } | null>(
    null,
  );
  const [showGenes, setShowGenes] = useState(true);
  const [showVariants, setShowVariants] = useState(true);
  const [showScores, setShowScores] = useState(true);
  const trackUrl =
    process.env.NEXT_PUBLIC_GENE_TRACK_URL && process.env.NEXT_PUBLIC_GENE_TRACK_URL.trim().length > 0
      ? process.env.NEXT_PUBLIC_GENE_TRACK_URL
      : '/genes.sample.bed';
  const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:18000';
  const apiToken = process.env.NEXT_PUBLIC_API_TOKEN;
  const scoreTrack =
    process.env.NEXT_PUBLIC_SCORE_TRACK_URL && process.env.NEXT_PUBLIC_SCORE_TRACK_URL.trim().length > 0
      ? process.env.NEXT_PUBLIC_SCORE_TRACK_URL
      : '/scores.sample.bed';

  useEffect(() => {
    let browser: any;
    const load = async () => {
      try {
        const mod = await import('igv/dist/igv.esm'); // ESM build exposes createBrowser
        const igv = (mod as any).default ?? mod;
        const decoded = decodeURIComponent(params.region);
        const { chr, start, end } = parseRegion(decoded);
        const tracks: any[] = [];
        if (showGenes) {
          tracks.push({
            name: 'Genes',
            type: 'annotation',
            format: 'bed',
            url: trackUrl,
          });
        }
        if (showVariants) {
          tracks.push({
            name: 'Variants',
            type: 'annotation',
            format: 'bed',
            url: `${apiUrl}/variants/bed?chr=${chr}&start=${start}&end=${end}${
              apiToken ? `&token=${apiToken}` : ''
            }`,
            headers: apiToken ? { Authorization: `Bearer ${apiToken}` } : undefined,
          });
        }
        if (showScores) {
          tracks.push({
            name: 'Scores',
            type: 'annotation',
            format: 'bed',
            url: scoreTrack,
            color: '#ef4444',
          });
        }

        browser = await igv.createBrowser(containerRef.current as HTMLDivElement, {
          genome: 'hg38',
          locus: `${chr}:${start}-${end}`,
          tracks,
        });

        browser.on('trackclick', (_t: any, pop: any) => {
          if (!pop || !pop.name) return;
          setSelected({
            track: pop.track || 'track',
            name: pop.name,
            pos: pop.locus || '',
          });
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
  }, [params.region, trackUrl, apiUrl, apiToken, scoreTrack, showGenes, showVariants, showScores]);

  return (
    <main className="min-h-screen p-6 space-y-4">
      <div className="flex items-center justify-between flex-wrap gap-3">
        <div>
          <h1 className="text-2xl font-semibold">Locus Explorer</h1>
          <p className="text-sm text-gray-600">Region: {decodeURIComponent(params.region)}</p>
        </div>
        <div className="flex items-center gap-2 text-sm">
          <label className="flex items-center gap-1">
            <input type="checkbox" checked={showGenes} onChange={() => setShowGenes(!showGenes)} />
            Genes
          </label>
          <label className="flex items-center gap-1">
            <input
              type="checkbox"
              checked={showVariants}
              onChange={() => setShowVariants(!showVariants)}
            />
            Variants
          </label>
          <label className="flex items-center gap-1">
            <input type="checkbox" checked={showScores} onChange={() => setShowScores(!showScores)} />
            Scores
          </label>
        </div>
      </div>

      {error ? (
        <div className="text-red-600 text-sm">{error}</div>
      ) : (
        <div className="grid md:grid-cols-[3fr,1fr] gap-4">
          <div
            ref={containerRef}
            className="border rounded shadow-sm"
            style={{ minHeight: 520 }}
          />
          <div className="border rounded p-3 bg-gray-50 min-h-[120px]">
            <div className="text-sm font-semibold mb-2">Feature</div>
            {selected ? (
              <div className="space-y-1 text-sm">
                <div className="text-gray-700">{selected.name}</div>
                {selected.pos && <div className="text-gray-500">{selected.pos}</div>}
                {selected.track && <div className="text-gray-500">Track: {selected.track}</div>}
              </div>
            ) : (
              <div className="text-sm text-gray-500">Click a feature to see details</div>
            )}
          </div>
        </div>
      )}
    </main>
  );
}
