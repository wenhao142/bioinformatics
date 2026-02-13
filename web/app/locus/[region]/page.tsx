'use client';

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';

const parseRegion = (region: string) => {
  const [chrPart, range] = region.split(':');
  if (!range) return { chr: chrPart, start: 1, end: 100000 };
  const [startStr, endStr] = range.split('-');
  const start = parseInt(startStr, 10) || 1;
  const end = parseInt(endStr, 10) || start + 1000;
  return { chr: chrPart, start, end };
};

type TrackKey = 'genes' | 'variants';
type TrackState = Record<TrackKey, boolean>;

type EvidenceField = {
  name: string;
  value: string;
};

type EvidenceCard = {
  id: string;
  createdAt: string;
  trackName: string;
  source: string;
  sourceKind: 'source' | 'inference';
  sourceLabel: string;
  note: string;
  fields: EvidenceField[];
  zoomLocus: string | null;
};

type LoginResponse = {
  access_token: string;
  token_type: string;
  expires_in: number;
};

type LiteratureRecord = {
  gene: string;
  pmid: string;
  title: string;
  year: number | null;
  source: string;
};

type VariantWindowSnp = {
  id: string;
  chr: string;
  pos: number;
  ref: string;
  alt: string;
  score: number;
  filter: string;
};

type TooltipPosition = {
  x: number;
  y: number;
};

const TRACK_ORDER: TrackKey[] = ['genes', 'variants'];
const TRACK_NAMES: Record<TrackKey, string> = {
  genes: 'Genes',
  variants: 'Variants',
};

const INITIAL_TRACK_STATE: TrackState = {
  genes: true,
  variants: true,
};

const SNP_HALF_WINDOW_BP = 60;
const VARIANT_WINDOW_HALF_BP = 600;

type LocusRange = {
  chr: string;
  start: number;
  end: number;
};

function toInt(value: string | number | null | undefined): number | null {
  if (value === null || value === undefined) {
    return null;
  }
  const cleaned = String(value).replace(/,/g, '').trim();
  if (!cleaned) {
    return null;
  }
  const parsed = Number.parseInt(cleaned, 10);
  return Number.isFinite(parsed) ? parsed : null;
}

function normalizeChromosome(raw: string | null | undefined): string | null {
  if (!raw) {
    return null;
  }
  const trimmed = raw.trim();
  if (!trimmed) {
    return null;
  }
  if (/^chr/i.test(trimmed)) {
    return trimmed;
  }
  if (/^(?:\d+|X|Y|M|MT)$/i.test(trimmed)) {
    return `chr${trimmed.toUpperCase()}`;
  }
  return trimmed;
}

function parseLocusText(raw: string): LocusRange | null {
  const match = raw.match(/(chr[0-9A-Za-z_]+):\s*([\d,]+)(?:\s*-\s*([\d,]+))?/i);
  if (!match) {
    return null;
  }
  const chr = normalizeChromosome(match[1]);
  const start = toInt(match[2]);
  const end = toInt(match[3]) ?? start;
  if (!chr || !start || !end) {
    return null;
  }
  const s = Math.min(start, end);
  const e = Math.max(start, end);
  return { chr, start: s, end: e };
}

function pickFieldValue(fields: EvidenceField[], keys: string[]): string | null {
  const lowered = keys.map((k) => k.toLowerCase());
  for (const field of fields) {
    if (lowered.includes(field.name.toLowerCase())) {
      return field.value;
    }
  }
  return null;
}

function deriveRefAltFields(fields: EvidenceField[]): EvidenceField[] {
  const refVal = pickFieldValue(fields, ['ref', 'reference', 'reference allele']);
  const altVal = pickFieldValue(fields, ['alt', 'alternate', 'alternate allele']);
  if (refVal && altVal) {
    return [
      { name: 'REF', value: refVal },
      { name: 'ALT', value: altVal },
    ];
  }

  const pattern = /\b([ACGTN]+)\s*>\s*([ACGTN]+)\b/i;
  for (const field of fields) {
    const m = field.value.match(pattern) ?? field.name.match(pattern);
    if (m) {
      return [
        { name: 'REF', value: m[1].toUpperCase() },
        { name: 'ALT', value: m[2].toUpperCase() },
      ];
    }
  }
  return [];
}

function extractFeatureLocus(fields: EvidenceField[]): LocusRange | null {
  const chrField = pickFieldValue(fields, ['chr', 'chrom', 'chromosome', 'contig']);
  const startField = pickFieldValue(fields, ['start']);
  const endField = pickFieldValue(fields, ['end', 'stop']);
  const posField = pickFieldValue(fields, ['pos', 'position', 'bp', 'genomic position']);
  const locusField = pickFieldValue(fields, ['locus', 'location']);

  const chr = normalizeChromosome(chrField);
  const start = toInt(startField);
  const end = toInt(endField);
  const pos = toInt(posField);

  if (chr && pos) {
    return { chr, start: pos, end: pos };
  }
  if (chr && start && end) {
    const s = Math.min(start, end);
    const e = Math.max(start, end);
    return { chr, start: s, end: e };
  }
  if (chr && start) {
    return { chr, start, end: start };
  }

  if (locusField) {
    const parsed = parseLocusText(locusField);
    if (parsed) {
      return parsed;
    }
  }

  for (const field of fields) {
    const parsed = parseLocusText(field.value);
    if (parsed) {
      return parsed;
    }
  }
  return null;
}

function getBrowserLocus(browser: any): LocusRange | null {
  const frame = browser?.referenceFrameList?.[0];
  if (!frame) {
    return null;
  }
  const chr = normalizeChromosome(frame.chr);
  const start = Math.max(1, Math.floor(frame.start ?? 1));
  const end = Math.max(start, Math.floor(frame.end ?? start + 100));
  if (!chr) {
    return null;
  }
  return { chr, start, end };
}

function buildSnpZoomLocus(featureLocus: LocusRange | null, fallback: LocusRange): LocusRange {
  const target = featureLocus ?? fallback;
  const center = Math.floor((target.start + target.end) / 2);
  const zoomStart = Math.max(1, center - SNP_HALF_WINDOW_BP);
  const zoomEnd = Math.max(zoomStart + 1, center + SNP_HALF_WINDOW_BP);
  return { chr: target.chr, start: zoomStart, end: zoomEnd };
}

function normalizePopupData(dataList: any[] | undefined): EvidenceField[] {
  if (!Array.isArray(dataList)) {
    return [];
  }

  const fields: EvidenceField[] = [];
  for (const item of dataList) {
    const key = typeof item?.name === 'string' && item.name.trim() ? item.name : 'value';
    const rawValue = item?.value;
    if (rawValue === undefined || rawValue === null) {
      continue;
    }
    const value = typeof rawValue === 'string' ? rawValue : JSON.stringify(rawValue);
    fields.push({ name: key, value });
  }
  return fields;
}

function closeStaleIgvPopups() {
  if (typeof document === 'undefined') {
    return;
  }
  // Keep only the newest popup so previous feature dialogs do not stack.
  const selectors = ['.igv-menu-popup', '.igv-popover'];
  for (const selector of selectors) {
    const nodes = Array.from(document.querySelectorAll(selector));
    const staleNodes = nodes.slice(0, Math.max(0, nodes.length - 1));
    staleNodes.forEach((node) => {
      if (node instanceof HTMLElement) {
        node.style.display = 'none';
      }
    });
  }
}

function extractHoveredGeneLocus(browser: any, event: MouseEvent): LocusRange | null {
  const target = event.target;
  if (!(target instanceof Node)) {
    return null;
  }

  const trackViews = Array.isArray(browser?.trackViews) ? browser.trackViews : [];
  for (const trackView of trackViews) {
    if (trackView?.track?.name !== TRACK_NAMES.genes) {
      continue;
    }
    const viewports = Array.isArray(trackView?.viewports) ? trackView.viewports : [];
    for (const viewport of viewports) {
      const viewportElement = viewport?.$viewport?.get?.(0) ?? viewport?.viewport;
      if (!(viewportElement instanceof HTMLElement) || !viewportElement.contains(target)) {
        continue;
      }
      if (typeof viewport.createClickState !== 'function' || typeof trackView.track.popupData !== 'function') {
        continue;
      }
      try {
        const clickState = viewport.createClickState(event);
        if (!clickState) {
          continue;
        }
        const popupData = trackView.track.popupData(clickState);
        const fields = normalizePopupData(popupData);
        const locus = extractFeatureLocus(fields);
        if (locus) {
          return locus;
        }
      } catch {
        // Ignore transient hover parsing errors from IGV internals.
        continue;
      }
    }
  }
  return null;
}

function trackSourceKind(trackName: string): 'source' | 'inference' {
  return 'source';
}

export default function LocusPage({ params }: { params: { region: string } }) {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const browserRef = useRef<any>(null);
  const viewerInitRef = useRef(0);
  const collapsingMultiLocusRef = useRef(false);
  const hoverWindowKeyRef = useRef('');
  const [error, setError] = useState<string | null>(null);
  const [viewerReady, setViewerReady] = useState(false);
  const [viewerClosed, setViewerClosed] = useState(false);
  const [viewerExpanded, setViewerExpanded] = useState(false);
  const [tracks, setTracks] = useState<TrackState>(INITIAL_TRACK_STATE);
  const [evidenceCards, setEvidenceCards] = useState<EvidenceCard[]>([]);
  const [currentLocus, setCurrentLocus] = useState('');
  const [apiToken, setApiToken] = useState<string | null>(process.env.NEXT_PUBLIC_API_TOKEN || null);
  const [authEmail, setAuthEmail] = useState('viewer@example.com');
  const [authPassword, setAuthPassword] = useState('password');
  const [authLoading, setAuthLoading] = useState(false);
  const [authError, setAuthError] = useState<string | null>(null);
  const [literatureRecords, setLiteratureRecords] = useState<LiteratureRecord[]>([]);
  const [literatureLoading, setLiteratureLoading] = useState(false);
  const [literatureError, setLiteratureError] = useState<string | null>(null);
  const [variantWindowOpen, setVariantWindowOpen] = useState(false);
  const [variantWindowLoading, setVariantWindowLoading] = useState(false);
  const [variantWindowError, setVariantWindowError] = useState<string | null>(null);
  const [variantWindowLabel, setVariantWindowLabel] = useState('');
  const [variantWindowSnps, setVariantWindowSnps] = useState<VariantWindowSnp[]>([]);
  const [variantWindowPos, setVariantWindowPos] = useState<TooltipPosition>({ x: 20, y: 20 });
  const genesTrackUrl =
    process.env.NEXT_PUBLIC_GENE_TRACK_URL && process.env.NEXT_PUBLIC_GENE_TRACK_URL.trim().length > 0
      ? process.env.NEXT_PUBLIC_GENE_TRACK_URL
      : '/genes.sample.bed';
  const apiUrl = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:18000';
  const decodedRegion = useMemo(() => decodeURIComponent(params.region), [params.region]);
  const { chr, start, end } = useMemo(() => parseRegion(decodedRegion), [decodedRegion]);

  useEffect(() => {
    setCurrentLocus(decodedRegion);
  }, [decodedRegion]);

  useEffect(() => {
    if (typeof window === 'undefined') {
      return;
    }
    if (apiToken) {
      return;
    }
    const localToken = window.localStorage.getItem('ad_api_token');
    if (localToken) {
      setApiToken(localToken);
    }
  }, [apiToken]);

  const invalidateToken = useCallback((reason?: string) => {
    setApiToken(null);
    if (reason) {
      setAuthError(reason);
    }
    if (typeof window !== 'undefined') {
      window.localStorage.removeItem('ad_api_token');
    }
  }, []);

  const handleLogin = async () => {
    try {
      setAuthError(null);
      setAuthLoading(true);
      const response = await fetch(`${apiUrl}/auth/login`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ email: authEmail, password: authPassword }),
      });
      if (!response.ok) {
        const text = await response.text();
        throw new Error(text || 'Login failed');
      }
      const payload: LoginResponse = await response.json();
      setApiToken(payload.access_token);
      if (typeof window !== 'undefined') {
        window.localStorage.setItem('ad_api_token', payload.access_token);
      }
    } catch (e: any) {
      setAuthError(e?.message || 'Login failed');
    } finally {
      setAuthLoading(false);
    }
  };

  const clearToken = () => {
    setAuthError(null);
    invalidateToken();
  };

  const hideVariantWindow = useCallback(() => {
    setVariantWindowOpen(false);
    setVariantWindowLoading(false);
    setVariantWindowError(null);
    hoverWindowKeyRef.current = '';
  }, []);

  const loadVariantWindow = useCallback(
    async (target: LocusRange) => {
      if (!apiToken) {
        setAuthError('Variants window needs login token.');
        return;
      }
      const center = Math.floor((target.start + target.end) / 2);
      const startBp = Math.max(1, center - VARIANT_WINDOW_HALF_BP);
      const endBp = Math.max(startBp + 1, center + VARIANT_WINDOW_HALF_BP);

      setVariantWindowOpen(true);
      setVariantWindowLoading(true);
      setVariantWindowError(null);
      setVariantWindowLabel(`${target.chr}:${startBp}-${endBp}`);

      try {
        const response = await fetch(
          `${apiUrl}/variants?chr=${encodeURIComponent(target.chr)}&start=${startBp}&end=${endBp}`,
          {
            headers: {
              Authorization: `Bearer ${apiToken}`,
            },
          }
        );
        if (response.status === 401) {
          invalidateToken('Token expired or invalid. Please sign in again.');
          setVariantWindowError('Token expired. Please sign in again.');
          setVariantWindowSnps([]);
          return;
        }
        if (!response.ok) {
          throw new Error(`Failed to load variants (${response.status})`);
        }
        const payload = await response.json().catch(() => ({}));
        const rows = Array.isArray(payload?.variants) ? payload.variants : [];
        const snps: VariantWindowSnp[] = rows
          .map((row: any, index: number) => {
            const pos = Number(row?.pos);
            if (!Number.isFinite(pos)) {
              return null;
            }
            const score = Number(row?.qual);
            return {
              id: `${row?.chr || target.chr}-${pos}-${row?.ref || 'N'}-${row?.alt || 'N'}-${index}`,
              chr: String(row?.chr || target.chr),
              pos,
              ref: String(row?.ref || 'N'),
              alt: String(row?.alt || 'N'),
              score: Number.isFinite(score) ? score : 0,
              filter: String(row?.filter || 'NA'),
            };
          })
          .filter(Boolean) as VariantWindowSnp[];
        snps.sort((a, b) => a.pos - b.pos);
        setVariantWindowSnps(snps);
        if (snps.length === 0) {
          setVariantWindowError('No SNPs found in this sliding window.');
        }
      } catch (e: any) {
        setVariantWindowSnps([]);
        setVariantWindowError(e?.message || 'Failed to load SNP window');
      } finally {
        setVariantWindowLoading(false);
      }
    },
    [apiToken, apiUrl, invalidateToken]
  );

  useEffect(() => {
    if (!viewerExpanded) {
      return;
    }
    const onKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') {
        setViewerExpanded(false);
      }
    };
    window.addEventListener('keydown', onKeyDown);
    return () => window.removeEventListener('keydown', onKeyDown);
  }, [viewerExpanded]);

  useEffect(() => {
    if (!apiToken) {
      return;
    }
    let cancelled = false;

    const verifyToken = async () => {
      try {
        const probeEnd = Math.min(end, start + 1000);
        const response = await fetch(
          `${apiUrl}/variants/bed?chr=${encodeURIComponent(chr)}&start=${start}&end=${probeEnd}&token=${encodeURIComponent(apiToken)}`,
          {
            headers: {
              Authorization: `Bearer ${apiToken}`,
            },
          }
        );
        if (cancelled) {
          return;
        }
        if (response.status === 401) {
          invalidateToken('Token expired or invalid. Please sign in again.');
          setError(null);
          return;
        }
        if (!response.ok) {
          setError(`Variants authorization check failed (${response.status}).`);
          return;
        }
        setAuthError(null);
      } catch (e: any) {
        if (!cancelled) {
          setError(e?.message || 'Failed to verify API token');
        }
      }
    };

    void verifyToken();

    return () => {
      cancelled = true;
    };
  }, [apiToken, apiUrl, chr, end, invalidateToken, start]);

  const selectedGene = useMemo(() => {
    const card = evidenceCards[0];
    if (!card) {
      return null;
    }
    const gene =
      pickFieldValue(card.fields, ['gene', 'gene_name', 'symbol']) ??
      (card.trackName === TRACK_NAMES.genes ? pickFieldValue(card.fields, ['name']) : null);
    if (!gene) {
      return null;
    }
    const trimmed = gene.trim();
    if (!trimmed || !/^[A-Za-z0-9._-]+$/.test(trimmed)) {
      return null;
    }
    return trimmed;
  }, [evidenceCards]);

  const variantScoreMax = useMemo(() => {
    let max = 0;
    for (const row of variantWindowSnps) {
      if (row.score > max) {
        max = row.score;
      }
    }
    return max > 0 ? max : 1;
  }, [variantWindowSnps]);

  useEffect(() => {
    if (!selectedGene) {
      setLiteratureRecords([]);
      setLiteratureError(null);
      setLiteratureLoading(false);
      return;
    }
    if (!apiToken) {
      setLiteratureRecords([]);
      setLiteratureError('Sign in to load PubMed metadata.');
      setLiteratureLoading(false);
      return;
    }

    let cancelled = false;
    const loadLiterature = async () => {
      try {
        setLiteratureLoading(true);
        setLiteratureError(null);
        const response = await fetch(
          `${apiUrl}/literature/pubmed/fetch?genes=${encodeURIComponent(selectedGene)}&max_per_gene=3`,
          {
            headers: {
              Authorization: `Bearer ${apiToken}`,
            },
          }
        );
        if (cancelled) {
          return;
        }
        if (response.status === 401) {
          invalidateToken('Token expired or invalid. Please sign in again.');
          setLiteratureRecords([]);
          return;
        }
        const payload = await response.json().catch(() => ({}));
        if (!response.ok) {
          if (response.status === 503) {
            setLiteratureRecords([]);
            setLiteratureError('PubMed is disabled in offline mode.');
            return;
          }
          const detail = typeof payload?.detail === 'string' ? payload.detail : null;
          throw new Error(detail || `Failed to load PubMed metadata (${response.status})`);
        }
        const records = Array.isArray(payload?.records)
          ? (payload.records.filter((row: any) => typeof row?.pmid === 'string') as LiteratureRecord[])
          : [];
        setLiteratureRecords(records);
        if (records.length === 0) {
          setLiteratureError(`No PubMed records found for ${selectedGene}.`);
        }
      } catch (e: any) {
        if (!cancelled) {
          setLiteratureRecords([]);
          setLiteratureError(e?.message || 'Failed to load PubMed metadata');
        }
      } finally {
        if (!cancelled) {
          setLiteratureLoading(false);
        }
      }
    };

    void loadLiterature();

    return () => {
      cancelled = true;
    };
  }, [apiToken, apiUrl, invalidateToken, selectedGene]);

  useEffect(() => {
    let cancelled = false;
    const initId = ++viewerInitRef.current;
    const containerEl = containerRef.current;
    let mountEl: HTMLDivElement | null = null;

    if (viewerClosed) {
      if (browserRef.current) {
        browserRef.current.dispose();
        browserRef.current = null;
      }
      if (containerEl) {
        containerEl.innerHTML = '';
      }
      setViewerReady(false);
      return;
    }

    const load = async () => {
      try {
        if (!containerEl) {
          return;
        }
        containerEl.innerHTML = '';
        mountEl = document.createElement('div');
        mountEl.className = 'igv-mount';
        mountEl.style.minHeight = '100%';
        mountEl.style.height = '100%';
        containerEl.appendChild(mountEl);

        setError(null);
        setViewerReady(false);
        setEvidenceCards([]);

        if (browserRef.current) {
          browserRef.current.dispose();
          browserRef.current = null;
        }

        const mod = await import('igv/dist/igv.esm'); // ESM build exposes createBrowser
        const igv = (mod as any).default ?? mod;
        const browser = await igv.createBrowser(mountEl as HTMLDivElement, {
          genome: 'hg38',
          locus: `${chr}:${start}-${end}`,
          tracks: [],
        });
        if (cancelled || initId !== viewerInitRef.current) {
          browser.dispose();
          if (mountEl?.parentElement === containerEl) {
            mountEl.remove();
          }
          return;
        }

        const locusChangeHandler = (frames: any[]) => {
          const frame = Array.isArray(frames) && frames.length > 0 ? frames[0] : null;
          const activeChr = normalizeChromosome(frame?.chr);
          const activeStart = toInt(frame?.start);
          const activeEnd = toInt(frame?.end);
          if (!activeChr || !activeStart || !activeEnd) {
            return;
          }
          const s = Math.max(1, Math.floor(activeStart));
          const e = Math.max(s, Math.floor(activeEnd));
          setCurrentLocus(`${activeChr}:${s}-${e}`);

          // Force single-locus mode to avoid duplicated IGV locus bars in UI.
          if (
            Array.isArray(frames) &&
            frames.length > 1 &&
            browserRef.current === browser &&
            !collapsingMultiLocusRef.current
          ) {
            collapsingMultiLocusRef.current = true;
            window.setTimeout(() => {
              void (async () => {
                try {
                  if (browserRef.current !== browser) {
                    return;
                  }
                  if (browser.referenceFrameList && browser.referenceFrameList.length > 1) {
                    await browser.search(`${activeChr}:${s}-${e}`);
                  }
                } catch {
                  // Ignore collapse failure and keep viewer alive.
                } finally {
                  collapsingMultiLocusRef.current = false;
                }
              })();
            }, 0);
          }
        };

        const trackClickHandler = (track: any, dataList: any[]) => {
          window.setTimeout(() => closeStaleIgvPopups(), 0);
          const trackName = typeof track?.name === 'string' ? track.name : 'Unknown';
          const source =
            trackName === TRACK_NAMES.genes
              ? genesTrackUrl
              : trackName === TRACK_NAMES.variants
                ? `${apiUrl}/variants/bed`
                : genesTrackUrl;

          const sourceKind = trackSourceKind(trackName);
          const popupFields = normalizePopupData(dataList);
          const fallbackLocus = getBrowserLocus(browser) ?? { chr, start, end };
          const clickedLocus = extractFeatureLocus(popupFields);
          const zoomTarget = buildSnpZoomLocus(clickedLocus, fallbackLocus);
          const zoomLocus = `${zoomTarget.chr}:${zoomTarget.start}-${zoomTarget.end}`;
          setCurrentLocus(zoomLocus);

          void (async () => {
            try {
              if (typeof browser.goto === 'function') {
                await browser.goto(zoomTarget.chr, zoomTarget.start, zoomTarget.end);
              } else if (typeof browser.search === 'function') {
                await browser.search(zoomLocus);
              }
            } catch {
              // Keep evidence capture even when zoom fails on a specific click.
            }
          })();

          const derivedRefAlt = deriveRefAltFields(popupFields);
          const enrichedFields =
            derivedRefAlt.length > 0 ? [...derivedRefAlt, ...popupFields] : popupFields;
          const note =
            sourceKind === 'source'
              ? 'Hover for quick tooltip; click pins details and zooms to SNP-level window.'
              : 'Inference signal. Hover for tooltip; click pins details and zooms to SNP-level window.';
          const card: EvidenceCard = {
            id: `${Date.now()}-${Math.random().toString(36).slice(2, 8)}`,
            createdAt: new Date().toISOString(),
            trackName,
            source,
            sourceKind,
            sourceLabel: sourceKind === 'source' ? 'Source Evidence' : 'Inference',
            note,
            fields: enrichedFields,
            zoomLocus,
          };
          setEvidenceCards([card]);
          return true;
        };

        browser.on('locuschange', locusChangeHandler);
        browser.on('trackclick', trackClickHandler);
        browserRef.current = browser;
        if (!cancelled) {
          setViewerReady(true);
        }
      } catch (e: any) {
        if (!cancelled) {
          setError(e?.message || 'Failed to load genome viewer');
        }
      }
    };

    load();

    return () => {
      cancelled = true;
      if (browserRef.current) {
        browserRef.current.dispose();
        browserRef.current = null;
      }
      collapsingMultiLocusRef.current = false;
      if (mountEl && mountEl.parentElement === containerEl) {
        mountEl.remove();
      }
      if (containerEl) {
        containerEl.innerHTML = '';
      }
      setViewerReady(false);
    };
  }, [apiToken, apiUrl, chr, end, genesTrackUrl, loadVariantWindow, start, viewerClosed]);

  useEffect(() => {
    if (!viewerReady || viewerClosed || !browserRef.current || !containerRef.current) {
      return;
    }
    const browser = browserRef.current;
    const containerEl = containerRef.current;
    let cancelled = false;
    let lastHoverAt = 0;

    const placeVariantWindow = (event: MouseEvent) => {
      const tooltipWidth = Math.min(560, Math.max(340, Math.floor(window.innerWidth * 0.42)));
      const tooltipHeight = 290;
      const gap = 14;
      const nextX = Math.min(
        Math.max(10, event.clientX + gap),
        Math.max(10, window.innerWidth - tooltipWidth - 10)
      );
      const nextY = Math.min(
        Math.max(10, event.clientY + gap),
        Math.max(10, window.innerHeight - tooltipHeight - 10)
      );
      setVariantWindowPos({ x: nextX, y: nextY });
    };

    const onMouseMove = (event: MouseEvent) => {
      if (cancelled || browserRef.current !== browser) {
        return;
      }
      const now = Date.now();
      if (now - lastHoverAt < 60) {
        return;
      }
      lastHoverAt = now;

      const locus = extractHoveredGeneLocus(browser, event);
      if (!locus) {
        hideVariantWindow();
        return;
      }
      placeVariantWindow(event);
      const center = Math.floor((locus.start + locus.end) / 2);
      const key = `${locus.chr}:${center}`;
      if (hoverWindowKeyRef.current === key) {
        setVariantWindowOpen(true);
        return;
      }
      hoverWindowKeyRef.current = key;
      void loadVariantWindow(locus);
    };

    const onMouseLeave = () => {
      hideVariantWindow();
    };

    containerEl.addEventListener('mousemove', onMouseMove);
    containerEl.addEventListener('mouseleave', onMouseLeave);
    return () => {
      cancelled = true;
      containerEl.removeEventListener('mousemove', onMouseMove);
      containerEl.removeEventListener('mouseleave', onMouseLeave);
    };
  }, [hideVariantWindow, loadVariantWindow, viewerClosed, viewerReady]);

  useEffect(() => {
    if (!viewerReady || !browserRef.current) {
      return;
    }

    let cancelled = false;

    const syncTracks = async () => {
      const browser = browserRef.current;
      if (!browser) {
        return;
      }

      const tokenPart = apiToken ? `&token=${encodeURIComponent(apiToken)}` : '';
      const configs: Record<TrackKey, any> = {
        genes: {
          name: TRACK_NAMES.genes,
          type: 'annotation',
          format: 'bed',
          url: genesTrackUrl,
        },
        variants: {
          name: TRACK_NAMES.variants,
          type: 'annotation',
          format: 'bed',
          url: `${apiUrl}/variants/bed?chr=${chr}&start=${start}&end=${end}${tokenPart}`,
          headers: apiToken ? { Authorization: `Bearer ${apiToken}` } : undefined,
        },
      };

      for (const key of TRACK_ORDER) {
        if (cancelled) {
          return;
        }
        if (key === 'variants' && !apiToken) {
          const existingNoToken = browser.findTracks('name', TRACK_NAMES[key])[0];
          if (existingNoToken) {
            browser.removeTrack(existingNoToken);
          }
          continue;
        }
        const existing = browser.findTracks('name', TRACK_NAMES[key])[0];
        if (tracks[key] && !existing) {
          await browser.loadTrack(configs[key]);
        } else if (!tracks[key] && existing) {
          browser.removeTrack(existing);
        }
      }
    };

    syncTracks().catch((e: any) => {
      if (!cancelled) {
        const message = e?.message || 'Failed to update tracks';
        if (/401|unauthorized/i.test(String(message))) {
          invalidateToken('Token expired or invalid. Please sign in again.');
          setError(null);
          return;
        }
        setError(message);
      }
    });

    return () => {
      cancelled = true;
    };
  }, [apiToken, apiUrl, chr, end, genesTrackUrl, invalidateToken, start, tracks, viewerReady]);

  return (
    <main className="locus-page">
      <h1 className="page-title">Locus Explorer</h1>
      <p className="region-line">Region: {decodedRegion}</p>
      <section className="auth-bar">
        {apiToken ? (
          <>
            <p className="auth-ok">API token loaded. Variants track is enabled.</p>
            <button type="button" className="auth-clear" onClick={clearToken}>
              Clear Token
            </button>
          </>
        ) : (
          <>
            <p className="auth-warn">Variants track needs login token.</p>
            <input
              className="auth-input"
              value={authEmail}
              onChange={(event) => setAuthEmail(event.target.value)}
              placeholder="email"
            />
            <input
              className="auth-input"
              type="password"
              value={authPassword}
              onChange={(event) => setAuthPassword(event.target.value)}
              placeholder="password"
            />
            <button type="button" className="auth-login" onClick={handleLogin} disabled={authLoading}>
              {authLoading ? 'Signing In...' : 'Sign In'}
            </button>
            {authError ? <span className="auth-error">{authError}</span> : null}
          </>
        )}
      </section>
      <section className="toggle-bar">
        {TRACK_ORDER.map((key) => (
          <label key={key} className="toggle-item">
            <input
              type="checkbox"
              checked={tracks[key]}
              disabled={key === 'variants' && !apiToken}
              onChange={(event) => {
                const checked = event.target.checked;
                setTracks((prev) => ({ ...prev, [key]: checked }));
              }}
            />
            <span>{TRACK_NAMES[key]}</span>
          </label>
        ))}
      </section>
      <section className="viewer-toolbar">
        <div className="viewer-toolbar-actions">
          <button
            type="button"
            className="viewer-toggle"
            onClick={() => {
              setError(null);
              setViewerClosed((prev) => !prev);
              if (!viewerClosed) {
                setViewerExpanded(false);
              }
            }}
          >
            {viewerClosed ? 'Open IGV' : 'Close IGV'}
          </button>
          {!viewerClosed ? (
            <button
              type="button"
              className="viewer-expand"
              onClick={() => setViewerExpanded((prev) => !prev)}
            >
              {viewerExpanded ? 'Restore Size' : 'Expand Viewer'}
            </button>
          ) : null}
        </div>
      </section>
      {error ? (
        <div className="error">{error}</div>
      ) : (
        <section className={`viewer-grid ${viewerExpanded ? 'viewer-grid-expanded' : ''}`}>
          {viewerExpanded ? (
            <button type="button" className="viewer-overlay-close" onClick={() => setViewerExpanded(false)}>
              Exit Fullscreen
            </button>
          ) : null}
          {viewerClosed ? (
            <div className="viewer viewer-closed">
              <p>IGV viewer is closed.</p>
              <button type="button" className="viewer-reopen" onClick={() => setViewerClosed(false)}>
                Reopen IGV
              </button>
            </div>
          ) : (
            <div ref={containerRef} className="viewer" />
          )}
          <aside className="evidence">
            <div className="evidence-head">
              <div>
                <h2>Feature Detail</h2>
                <p className="evidence-focus">Focus: {currentLocus || decodedRegion}</p>
              </div>
              {evidenceCards.length > 0 ? (
                <button type="button" className="evidence-clear" onClick={() => setEvidenceCards([])}>
                  Clear Selection
                </button>
              ) : null}
            </div>
            {evidenceCards.length === 0 ? (
              <p className="placeholder">
                Hover on a gene/variant in IGV to see quick info, then click to pin REF/ALT and locus details here.
              </p>
            ) : (
              <div className="evidence-list">
                {evidenceCards.map((card) => (
                  <article key={card.id} className="evidence-card">
                    <div className="evidence-meta">
                      <div className={`evidence-pill ${card.sourceKind === 'source' ? 'pill-source' : 'pill-infer'}`}>
                        {card.sourceLabel}
                      </div>
                      <p>
                        <strong>Track:</strong> {card.trackName}
                      </p>
                      <p>
                        <strong>Type:</strong> {card.sourceKind}
                      </p>
                      <p>
                        <strong>Source:</strong> {card.source}
                      </p>
                      <p>
                        <strong>Captured:</strong> {new Date(card.createdAt).toLocaleString()}
                      </p>
                      {card.zoomLocus ? (
                        <p>
                          <strong>Zoom:</strong> {card.zoomLocus}
                        </p>
                      ) : null}
                      <p>{card.note}</p>
                    </div>
                    <div className="evidence-fields">
                      {card.fields.length === 0 ? (
                        <p className="placeholder">No popup fields returned for this click.</p>
                      ) : (
                        card.fields.map((field, index) => (
                          <div key={`${card.id}-${field.name}-${index}`} className="field-row">
                            <span className="field-name">{field.name}</span>
                            <span className="field-value">{field.value}</span>
                          </div>
                        ))
                      )}
                    </div>
                  </article>
                ))}
                <section className="literature-card">
                  <h3>PubMed</h3>
                  <p className="literature-subtitle">
                    Gene: <strong>{selectedGene ?? 'N/A'}</strong>
                  </p>
                  {literatureLoading ? (
                    <p className="placeholder">Loading PubMed metadata...</p>
                  ) : literatureError ? (
                    <p className="literature-error">{literatureError}</p>
                  ) : literatureRecords.length === 0 ? (
                    <p className="placeholder">No linked PubMed records.</p>
                  ) : (
                    <div className="literature-list">
                      {literatureRecords.map((record) => (
                        <article key={`${record.gene}-${record.pmid}`} className="literature-item">
                          <a
                            href={`https://pubmed.ncbi.nlm.nih.gov/${record.pmid}/`}
                            target="_blank"
                            rel="noreferrer"
                          >
                            PMID {record.pmid}
                          </a>
                          <p>{record.title}</p>
                          <span>{record.year ?? 'Year unknown'}</span>
                        </article>
                      ))}
                    </div>
                  )}
                </section>
              </div>
            )}
          </aside>
        </section>
      )}
      {variantWindowOpen ? (
        <section
          className="variant-window"
          style={{ left: `${variantWindowPos.x}px`, top: `${variantWindowPos.y}px` }}
          aria-live="polite"
        >
          <div className="variant-window-head">
            <div>
              <h3>SNP Score</h3>
              <p>{variantWindowLabel}</p>
            </div>
          </div>
          {variantWindowLoading ? (
            <p className="placeholder">Loading SNP scores...</p>
          ) : variantWindowError ? (
            <p className="variant-window-error">{variantWindowError}</p>
          ) : (
            <div className="variant-window-scroll">
              <div className="snp-chart-strip">
                {variantWindowSnps.map((snp) => (
                  <div key={`bar-${snp.id}`} className="snp-bar-item" title={`${snp.chr}:${snp.pos}`}>
                    <span>{snp.pos}</span>
                    <div
                      className="snp-bar"
                      style={{ height: `${Math.max(10, Math.round((snp.score / variantScoreMax) * 140))}px` }}
                    />
                    <span>{snp.score.toFixed(1)}</span>
                  </div>
                ))}
              </div>
            </div>
          )}
        </section>
      ) : null}
      <style jsx>{`
        .locus-page {
          --paper: #eef4f4;
          --paper-soft: #f8fcfc;
          --ink: #142328;
          --ink-muted: #3d565d;
          --panel: #ffffff;
          --line: #9fb8bf;
          --accent: #cc4a2d;
          --accent-soft: #ef8c4b;
          --deep: #0f3f4d;
          min-height: 100vh;
          padding: 28px;
          font-family: 'IBM Plex Sans', 'Segoe UI', sans-serif;
          font-size: 14pt;
          line-height: 1.35;
          color: var(--ink);
          background: #e6f0f2;
          animation: pageIn 360ms ease-out;
          overflow-x: hidden;
        }
        .locus-page :is(p, span, label, input, button, a, h1, h2, h3) {
          font-size: 14pt !important;
        }
        .page-title {
          margin: 0 0 4px;
          font-family: 'IBM Plex Serif', Georgia, serif;
          font-size: clamp(40px, 4.8vw, 58px);
          line-height: 1.15;
          letter-spacing: 0.01em;
          color: #163641;
          text-shadow: 0 1px 0 #fff;
        }
        .region-line {
          margin: 0 0 18px;
          display: inline-block;
          padding: 12px 18px;
          border: 1px solid var(--line);
          border-radius: 999px;
          font-weight: 700;
          color: var(--ink-muted);
          background: rgba(250, 255, 255, 0.9);
          backdrop-filter: blur(2px);
        }
        .toggle-bar {
          display: flex;
          gap: 10px;
          flex-wrap: wrap;
          margin-bottom: 16px;
          padding: 14px;
          border: 1px solid var(--line);
          border-radius: 14px;
          background: #f8efdc;
          box-shadow:
            0 8px 24px rgba(44, 63, 52, 0.08),
            inset 0 1px 0 rgba(255, 255, 255, 0.7);
        }
        .viewer-toolbar {
          display: flex;
          justify-content: flex-end;
          margin-bottom: 12px;
        }
        .viewer-toolbar-actions {
          display: flex;
          gap: 8px;
          align-items: center;
        }
        .viewer-toggle,
        .viewer-reopen,
        .viewer-expand {
          min-height: 56px;
          border-radius: 10px;
          border: 1px solid #2f4a3a;
          background: #2f4a3a;
          color: #fff6e8;
          font-weight: 700;
          letter-spacing: 0.05em;
          text-transform: uppercase;
          padding: 8px 16px;
          cursor: pointer;
        }
        .viewer-reopen {
          background: #b84f28;
          border-color: #b84f28;
        }
        .viewer-expand {
          background: #0f3f4d;
          border-color: #0f3f4d;
        }
        .auth-bar {
          display: flex;
          align-items: center;
          gap: 8px;
          flex-wrap: wrap;
          margin-bottom: 12px;
          padding: 10px;
          border: 1px solid #ccb891;
          border-radius: 12px;
          background: #fff8e8;
        }
        .auth-ok,
        .auth-warn {
          margin: 0;
          font-size: 12px;
          font-weight: 700;
          letter-spacing: 0.04em;
          text-transform: uppercase;
        }
        .auth-ok {
          color: #234031;
        }
        .auth-warn {
          color: #8a3f20;
        }
        .auth-input {
          min-height: 56px;
          min-width: 150px;
          border: 1px solid #c9b284;
          border-radius: 8px;
          padding: 8px 12px;
          background: #fffcf5;
          color: #213329;
        }
        .auth-input:focus {
          outline: 2px solid #e98e4d;
          outline-offset: 1px;
        }
        .auth-login,
        .auth-clear {
          min-height: 56px;
          border-radius: 8px;
          border: 1px solid #b84f28;
          background: #b84f28;
          color: #fff6e8;
          font-weight: 700;
          letter-spacing: 0.04em;
          text-transform: uppercase;
          padding: 8px 14px;
          cursor: pointer;
        }
        .auth-login:disabled {
          opacity: 0.7;
          cursor: default;
        }
        .auth-clear {
          background: #234031;
          border-color: #234031;
        }
        .auth-error {
          color: #9a2b23;
          font-size: 12px;
          font-weight: 700;
        }
        .toggle-item {
          display: flex;
          align-items: center;
          gap: 9px;
          font-weight: 700;
          text-transform: uppercase;
          letter-spacing: 0.06em;
          color: #2f4035;
          border: 1px solid #c6b48d;
          border-radius: 999px;
          padding: 10px 16px;
          background: #fff9ec;
          transition: transform 180ms ease, border-color 180ms ease, box-shadow 180ms ease;
        }
        .toggle-item input:disabled {
          opacity: 0.5;
        }
        .toggle-item:hover {
          transform: translateY(-1px);
          border-color: #bb6b3e;
          box-shadow: 0 4px 14px rgba(184, 79, 40, 0.18);
        }
        .toggle-item input {
          width: 24px;
          height: 24px;
          accent-color: #b84f28;
        }
        .error {
          color: #7f1d1d;
          font-size: 14px;
          font-weight: 700;
          padding: 12px 14px;
          border: 2px solid #d57a51;
          border-radius: 10px;
          background: #fceadf;
        }
        .viewer-grid {
          display: grid;
          grid-template-columns: minmax(0, 2fr) minmax(280px, 1fr);
          gap: 18px;
          align-items: start;
          overflow-x: hidden;
        }
        .viewer-grid-expanded {
          position: fixed;
          inset: 10px;
          z-index: 120;
          background: rgba(238, 244, 244, 0.98);
          border: 2px solid #9fb8bf;
          border-radius: 16px;
          padding: 12px;
          grid-template-columns: 1fr;
          box-shadow: 0 24px 64px rgba(11, 31, 37, 0.35);
          backdrop-filter: blur(3px);
          overflow: auto;
        }
        .viewer-overlay-close {
          position: sticky;
          top: 0;
          margin-left: auto;
          z-index: 2;
          min-height: 44px;
          border-radius: 8px;
          border: 1px solid #0f3f4d;
          background: #0f3f4d;
          color: #ffffff;
          padding: 6px 12px;
          cursor: pointer;
        }
        .viewer-grid-expanded .evidence {
          display: none;
        }
        .viewer-grid-expanded .viewer {
          min-height: calc(100vh - 56px);
        }
        .viewer {
          min-height: 560px;
          border: 2px solid #cbb78c;
          border-radius: 18px;
          background: var(--panel);
          box-shadow:
            0 16px 34px rgba(26, 45, 34, 0.14),
            inset 0 0 0 1px rgba(255, 255, 255, 0.7);
          overflow: hidden;
        }
        .viewer-closed {
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          gap: 10px;
          color: #425249;
          font-size: 14px;
        }
        .viewer-closed p {
          margin: 0;
          font-weight: 700;
          letter-spacing: 0.02em;
        }
        .evidence {
          border: 2px solid #cbb78c;
          border-radius: 18px;
          background: var(--panel);
          padding: 16px;
          box-shadow:
            0 16px 34px rgba(26, 45, 34, 0.14),
            inset 0 0 0 1px rgba(255, 255, 255, 0.7);
          min-height: 560px;
          display: flex;
          flex-direction: column;
          overflow: hidden;
        }
        .evidence-head {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 10px;
          margin-bottom: 10px;
        }
        .evidence-clear {
          min-height: 52px;
          border-radius: 8px;
          border: 1px solid #9a3f21;
          background: #9a3f21;
          color: #fff6e8;
          padding: 8px 14px;
          font-weight: 700;
          letter-spacing: 0.04em;
          text-transform: uppercase;
          cursor: pointer;
        }
        .evidence h2 {
          margin: 0;
          font-family: 'Palatino Linotype', 'Book Antiqua', Georgia, serif;
          font-size: 22px;
          color: var(--deep);
        }
        .evidence-focus {
          margin: 4px 0 0;
          font-size: 12px;
          color: #5f6f64;
          font-weight: 700;
          letter-spacing: 0.03em;
        }
        .evidence-list {
          display: flex;
          flex-direction: column;
          gap: 10px;
          flex: 1;
          min-height: 0;
          overflow-y: auto;
          overflow-x: hidden;
          padding-right: 6px;
        }
        .evidence-list::-webkit-scrollbar {
          width: 10px;
        }
        .evidence-list::-webkit-scrollbar-track {
          background: #f3e7cd;
          border-radius: 999px;
        }
        .evidence-list::-webkit-scrollbar-thumb {
          background: #b56b3d;
          border-radius: 999px;
          border: 2px solid #f3e7cd;
        }
        .evidence-list::-webkit-scrollbar-thumb:hover {
          background: #934f27;
        }
        .evidence-card {
          border: 1px solid #d8c59d;
          border-radius: 12px;
          background: #fffdf7;
          padding: 10px;
          box-shadow: 0 3px 12px rgba(32, 47, 38, 0.08);
        }
        .evidence-meta {
          margin-bottom: 14px;
          padding: 12px;
          background: #f8f3e7;
          border-radius: 10px;
          border: 1px solid #dbcaa6;
          font-size: 13px;
        }
        .evidence-pill {
          display: inline-flex;
          align-items: center;
          border-radius: 999px;
          padding: 3px 10px;
          font-size: 11px;
          font-weight: 800;
          letter-spacing: 0.06em;
          text-transform: uppercase;
          margin-bottom: 8px;
          border: 1px solid transparent;
        }
        .pill-source {
          background: #e8f5ec;
          color: #1f5a35;
          border-color: #84c69b;
        }
        .pill-infer {
          background: #fff0e2;
          color: #8a3f20;
          border-color: #e5a676;
        }
        .evidence-meta p {
          margin: 0 0 6px;
          overflow-wrap: anywhere;
        }
        .evidence-meta p:last-child {
          margin-bottom: 0;
        }
        .evidence-fields {
          display: flex;
          flex-direction: column;
          gap: 10px;
        }
        .field-row {
          border: 1px solid #d8c59d;
          border-radius: 10px;
          padding: 10px;
          background: #fffcf4;
          box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.85);
        }
        .field-name {
          display: block;
          font-size: 11px;
          color: #6c7a70;
          margin-bottom: 5px;
          text-transform: uppercase;
          letter-spacing: 0.07em;
          font-weight: 700;
        }
        .field-value {
          display: block;
          font-size: 13px;
          color: #18261f;
          overflow-wrap: anywhere;
        }
        .placeholder {
          margin: 0;
          color: #5f6f64;
          font-size: 13px;
          line-height: 1.5;
        }
        .literature-card {
          border: 1px solid #cdb792;
          border-radius: 12px;
          background: #fef7e7;
          padding: 12px;
        }
        .literature-card h3 {
          margin: 0 0 6px;
          font-size: 15px;
          font-weight: 800;
          letter-spacing: 0.04em;
          text-transform: uppercase;
          color: #2c4638;
        }
        .literature-subtitle {
          margin: 0 0 10px;
          font-size: 12px;
          color: #4a5c52;
        }
        .literature-list {
          display: flex;
          flex-direction: column;
          gap: 10px;
        }
        .literature-item {
          border: 1px solid #d8c59d;
          border-radius: 10px;
          padding: 10px;
          background: #fffdf8;
        }
        .literature-item a {
          color: #9a3f21;
          font-size: 12px;
          font-weight: 700;
          text-decoration: none;
        }
        .literature-item a:hover {
          text-decoration: underline;
        }
        .literature-item p {
          margin: 6px 0;
          font-size: 13px;
          color: #18261f;
          line-height: 1.4;
        }
        .literature-item span {
          font-size: 12px;
          color: #5f6f64;
          font-weight: 700;
        }
        .literature-error {
          margin: 0;
          color: #9a2b23;
          font-size: 13px;
          font-weight: 700;
        }
        .variant-window {
          position: fixed;
          width: min(560px, calc(100vw - 20px));
          max-height: 290px;
          border: 2px solid #8cb1ba;
          border-radius: 16px;
          background: #f7fcfd;
          box-shadow: 0 20px 56px rgba(12, 34, 41, 0.32);
          z-index: 2147483647;
          display: flex;
          flex-direction: column;
          isolation: isolate;
          pointer-events: none;
        }
        .variant-window-head {
          display: flex;
          align-items: center;
          justify-content: flex-start;
          gap: 10px;
          padding: 10px 12px;
          border-bottom: 1px solid #b7d0d6;
          background: #edf7f9;
        }
        .variant-window-head h3 {
          margin: 0;
          font-size: 12pt !important;
          font-weight: 800;
          color: #103946;
          text-transform: uppercase;
          letter-spacing: 0.04em;
        }
        .variant-window-head p {
          margin: 2px 0 0;
          font-size: 10pt !important;
          color: #3a5961;
        }
        .variant-window-error {
          margin: 0;
          padding: 10px 12px;
          font-size: 11pt !important;
          color: #8f251f;
          font-weight: 700;
        }
        .variant-window-scroll {
          overflow-x: auto;
          overflow-y: hidden;
          padding: 8px 10px 12px;
          display: flex;
          flex-direction: row;
          align-items: flex-end;
        }
        .snp-chart-strip {
          display: flex;
          align-items: flex-end;
          gap: 8px;
          min-width: max-content;
          padding: 4px 6px;
        }
        .snp-bar-item {
          display: flex;
          flex-direction: column;
          align-items: center;
          gap: 4px;
          min-width: 62px;
          font-size: 9pt !important;
        }
        .snp-bar {
          width: 36px;
          border-radius: 7px 7px 3px 3px;
          background: #cc4a2d;
          box-shadow: inset 0 0 0 1px rgba(255, 255, 255, 0.35);
        }
        @keyframes pageIn {
          from {
            opacity: 0.4;
            transform: translateY(5px);
          }
          to {
            opacity: 1;
            transform: translateY(0);
          }
        }
        @media (max-width: 980px) {
          .locus-page {
            padding: 18px;
          }
          .auth-input {
            min-width: 130px;
          }
          .viewer-grid {
            grid-template-columns: 1fr;
          }
          .viewer-grid-expanded {
            inset: 0;
            border-radius: 0;
            padding: 8px;
          }
          .viewer,
          .evidence {
            min-height: 420px;
          }
          .variant-window {
            left: 8px;
            width: auto;
            max-height: 240px;
          }
        }
      `}</style>
      <style jsx global>{`
        /* IGV visual language: slate + orange, non-purple */
        .igv-root-div,
        .igv-container,
        .igv-column-container {
          font-family: 'IBM Plex Sans', 'Segoe UI', sans-serif !important;
        }
        .igv-navbar {
          min-height: 48px !important;
          height: 48px !important;
          background: #1f4335 !important;
          border-color: #1b392c !important;
          color: #fff6e4 !important;
        }
        .igv-navbar .igv-current-genome {
          color: #fff0d8 !important;
          font-size: 14pt !important;
          font-weight: 700 !important;
          letter-spacing: 0.03em !important;
        }
        .igv-navbar .igv-navbar-left-container .igv-navbar-genomic-location .igv-locus-size-group .igv-search-container {
          height: 30px !important;
        }
        .igv-navbar
          .igv-navbar-left-container
          .igv-navbar-genomic-location
          .igv-locus-size-group
          .igv-search-container
          input.igv-search-input {
          height: 44px !important;
          border-color: #d58d45 !important;
          border-width: 2px !important;
          font-size: 14pt !important;
          color: #1f2a24 !important;
          background-color: #fff9ec !important;
        }
        .igv-navbar .igv-navbar-button {
          min-height: 44px !important;
          line-height: 42px !important;
          border-radius: 8px !important;
          border-width: 2px !important;
          border-color: #d58d45 !important;
          background-color: #fff2d6 !important;
          color: #2f3d34 !important;
          font-size: 14pt !important;
          font-weight: 700 !important;
          letter-spacing: 0.03em !important;
          padding-left: 10px !important;
          padding-right: 10px !important;
        }
        .igv-navbar .igv-navbar-button.igv-navbar-button-clicked {
          background-color: #b84f28 !important;
          border-color: #b84f28 !important;
          color: #fff9ef !important;
        }
        .igv-navbar .igv-navbar-right-container .igv-zoom-widget div:first-child,
        .igv-navbar .igv-navbar-right-container .igv-zoom-widget div:last-child {
          width: 44px !important;
          height: 44px !important;
          color: #ffd9a3 !important;
        }
        .igv-navbar .igv-navbar-right-container .igv-zoom-widget svg {
          width: 30px !important;
          height: 30px !important;
        }
        .igv-gear-menu-column > div > div {
          width: 24px !important;
          height: 24px !important;
          color: #234031 !important;
        }
        .igv-gear-menu-column > div > div > svg {
          width: 24px !important;
          height: 24px !important;
        }
        .igv-track-drag-column > .igv-track-drag-handle {
          background-color: #b84f28 !important;
          border-radius: 8px !important;
        }
        .igv-track-label {
          font-size: 14pt !important;
          font-weight: 700 !important;
          border-width: 2px !important;
          border-color: #2d493a !important;
          background-color: #fff5df !important;
          color: #1f2a24 !important;
          padding: 2px 8px !important;
          border-radius: 8px !important;
        }
        .igv-menu-popup {
          border-width: 2px !important;
          border-color: #2d493a !important;
          box-shadow: 0 12px 26px rgba(20, 33, 27, 0.2) !important;
        }
        .igv-menu-popup-header {
          background-color: #e9d7b8 !important;
        }
        .igv-axis-column,
        .igv-scrollbar-column,
        .igv-track-drag-column,
        .igv-gear-menu-column {
          background-color: #fff9ec !important;
        }
        .igv-scrollbar-column {
          width: 16px !important;
        }
        .igv-scrollbar-column > div {
          width: 16px !important;
        }
        .igv-scrollbar-column > div > div {
          width: 12px !important;
          left: 1px !important;
          border-color: #89958f !important;
          background-color: rgba(184, 79, 40, 0.28) !important;
        }
        .igv-scrollbar-column > div > div:hover {
          background-color: rgba(184, 79, 40, 0.46) !important;
        }
      `}</style>
    </main>
  );
}
