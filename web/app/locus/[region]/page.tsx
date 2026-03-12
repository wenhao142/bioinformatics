'use client';

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { resolveApiBase } from '@/lib/api-base';

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

type CausalRunRecord = {
  run_id: string;
  project_id: string | null;
  params: {
    chr: string;
    start: number;
    end: number;
    top_n: number;
    ld_window_bp: number;
  };
  counts?: {
    variant_scores?: number;
    gene_scores?: number;
  };
};

type CausalVariantScore = {
  chr: string;
  pos: number;
  ref: string;
  alt: string;
  gene: string | null;
  score: number;
  ld_proxy_signal: number;
  expr_signal: number;
  annotation_weight: number;
};

type CausalGeneScore = {
  gene: string;
  score: number;
  rank: number;
  variant_count: number;
  top_variant_score: number;
  mean_top3_variant_score: number;
  expr_signal: number;
};

type CausalResult = {
  region: {
    chr: string;
    start: number;
    end: number;
  };
  lead_variant?: {
    chr: string;
    pos: number;
    qual: number | null;
  };
  variant_scores: CausalVariantScore[];
  gene_scores: CausalGeneScore[];
  meta?: {
    ld_window_bp?: number;
    variants_considered?: number;
    genes_scored?: number;
  };
};

type WorkflowRunSummary = {
  run_id: string;
  workflow_id: string;
  created_by?: string;
  submitted_runs: number;
  completed_runs: number;
  failed_runs?: number;
  duration_ms: number;
};

type WorkflowRunDetail = {
  summary: WorkflowRunSummary;
  results: Array<{
    index: number;
    status: string;
    engine: string;
    workflow_artifacts?: {
      run_dir: string;
      snakefile: string;
      configfile: string;
    };
    outputs?: string[];
  }>;
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
const VARIANT_WINDOW_HALF_BP = 1000;
const VARIANT_WINDOW_STEP_BP = 25;

type LocusRange = {
  chr: string;
  start: number;
  end: number;
};

type HoveredGeneContext = {
  featureLocus: LocusRange;
  hoverPos: number;
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

function extractHoveredGeneContext(browser: any, event: MouseEvent): HoveredGeneContext | null {
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
        const featureLocus = extractFeatureLocus(fields);
        const referenceFrame = clickState?.referenceFrame ?? viewport?.referenceFrame ?? browser?.referenceFrameList?.[0];
        const frameChr = normalizeChromosome(referenceFrame?.chr);

        let hoverPos =
          toInt(clickState?.genomicLocation) ?? toInt(clickState?.genomicStart) ?? toInt(clickState?.bp);
        if (
          (!hoverPos || hoverPos < 1) &&
          Number.isFinite(Number(referenceFrame?.start)) &&
          Number.isFinite(Number(referenceFrame?.bpPerPixel))
        ) {
          const rect = viewportElement.getBoundingClientRect();
          const offsetX = Math.max(0, Math.min(rect.width, event.clientX - rect.left));
          const guessedPos = Math.round(Number(referenceFrame.start) + offsetX * Number(referenceFrame.bpPerPixel));
          if (Number.isFinite(guessedPos) && guessedPos >= 1) {
            hoverPos = guessedPos;
          }
        }

        const chr = frameChr ?? featureLocus?.chr ?? null;
        if (!chr) {
          continue;
        }
        const locus =
          featureLocus !== null
            ? { chr, start: featureLocus.start, end: featureLocus.end }
            : hoverPos
              ? { chr, start: hoverPos, end: hoverPos }
              : null;
        if (!locus) {
          continue;
        }
        const center = hoverPos ?? Math.floor((locus.start + locus.end) / 2);
        return { featureLocus: locus, hoverPos: center };
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
  const hoverRequestSeqRef = useRef(0);
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
  const [causalRuns, setCausalRuns] = useState<CausalRunRecord[]>([]);
  const [causalRunsLoading, setCausalRunsLoading] = useState(false);
  const [causalRunsError, setCausalRunsError] = useState<string | null>(null);
  const [causalRunsRefreshKey, setCausalRunsRefreshKey] = useState(0);
  const [selectedCausalRunId, setSelectedCausalRunId] = useState('');
  const [causalResult, setCausalResult] = useState<CausalResult | null>(null);
  const [causalResultLoading, setCausalResultLoading] = useState(false);
  const [causalResultError, setCausalResultError] = useState<string | null>(null);
  const [workflowRuns, setWorkflowRuns] = useState<WorkflowRunSummary[]>([]);
  const [workflowRunsLoading, setWorkflowRunsLoading] = useState(false);
  const [workflowRunsError, setWorkflowRunsError] = useState<string | null>(null);
  const [workflowRunsRefreshKey, setWorkflowRunsRefreshKey] = useState(0);
  const [selectedWorkflowRunId, setSelectedWorkflowRunId] = useState('');
  const [workflowRunDetail, setWorkflowRunDetail] = useState<WorkflowRunDetail | null>(null);
  const [workflowRunDetailLoading, setWorkflowRunDetailLoading] = useState(false);
  const [workflowRunDetailError, setWorkflowRunDetailError] = useState<string | null>(null);
  const [variantWindowOpen, setVariantWindowOpen] = useState(false);
  const [variantWindowLoading, setVariantWindowLoading] = useState(false);
  const [variantWindowError, setVariantWindowError] = useState<string | null>(null);
  const [variantWindowLabel, setVariantWindowLabel] = useState('');
  const [variantWindowSnps, setVariantWindowSnps] = useState<VariantWindowSnp[]>([]);
  const genesTrackUrl =
    process.env.NEXT_PUBLIC_GENE_TRACK_URL && process.env.NEXT_PUBLIC_GENE_TRACK_URL.trim().length > 0
      ? process.env.NEXT_PUBLIC_GENE_TRACK_URL
      : '/genes.sample.bed';
  const apiUrl = useMemo(() => resolveApiBase(), []);
  const decodedRegion = useMemo(() => decodeURIComponent(params.region), [params.region]);
  const { chr, start, end } = useMemo(() => parseRegion(decodedRegion), [decodedRegion]);
  const [navLocusInput, setNavLocusInput] = useState(decodedRegion);

  useEffect(() => {
    setCurrentLocus(decodedRegion);
    setNavLocusInput(decodedRegion);
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

  useEffect(() => {
    if (!currentLocus) {
      return;
    }
    setNavLocusInput(currentLocus);
  }, [currentLocus]);

  const goToLocus = useCallback(
    async (rawLocus?: string) => {
      const browser = browserRef.current;
      if (!browser) {
        return;
      }
      const locus = (rawLocus ?? navLocusInput).trim();
      if (!locus) {
        return;
      }
      try {
        if (typeof browser.search === 'function') {
          await browser.search(locus);
        }
        setError(null);
      } catch (e: any) {
        setError(e?.message || 'Failed to navigate locus');
      }
    },
    [navLocusInput]
  );

  const zoomBy = useCallback(async (factor: number) => {
    const browser = browserRef.current;
    if (!browser) {
      return;
    }
    const locus = getBrowserLocus(browser);
    if (!locus) {
      return;
    }
    const span = Math.max(2, locus.end - locus.start + 1);
    const nextSpan = Math.max(2, Math.floor(span * factor));
    const center = Math.floor((locus.start + locus.end) / 2);
    const nextStart = Math.max(1, center - Math.floor(nextSpan / 2));
    const nextEnd = Math.max(nextStart + 1, nextStart + nextSpan - 1);
    try {
      if (typeof browser.goto === 'function') {
        await browser.goto(locus.chr, nextStart, nextEnd);
      } else if (typeof browser.search === 'function') {
        await browser.search(`${locus.chr}:${nextStart}-${nextEnd}`);
      }
      setError(null);
    } catch (e: any) {
      setError(e?.message || 'Failed to zoom');
    }
  }, []);

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
    hoverRequestSeqRef.current += 1;
    setVariantWindowOpen(false);
    setVariantWindowLoading(false);
    setVariantWindowError(null);
    hoverWindowKeyRef.current = '';
  }, []);

  const loadVariantWindow = useCallback(
    async (target: LocusRange, centerOverride?: number) => {
      if (!apiToken) {
        setAuthError('Variants window needs login token.');
        return;
      }
      const requestSeq = hoverRequestSeqRef.current + 1;
      hoverRequestSeqRef.current = requestSeq;
      const center = centerOverride ?? Math.floor((target.start + target.end) / 2);
      const startBp = Math.max(1, center - VARIANT_WINDOW_HALF_BP);
      const endBp = Math.max(startBp + 1, center + VARIANT_WINDOW_HALF_BP);

      setVariantWindowOpen(true);
      setVariantWindowLoading(true);
      setVariantWindowError(null);
      setVariantWindowLabel(`${target.chr}:${startBp}-${endBp} (hover ${center} bp)`);

      try {
        const response = await fetch(
          `${apiUrl}/variants?chr=${encodeURIComponent(target.chr)}&start=${startBp}&end=${endBp}`,
          {
            headers: {
              Authorization: `Bearer ${apiToken}`,
            },
          }
        );
        if (requestSeq !== hoverRequestSeqRef.current) {
          return;
        }
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
        if (requestSeq !== hoverRequestSeqRef.current) {
          return;
        }
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
        if (requestSeq !== hoverRequestSeqRef.current) {
          return;
        }
        setVariantWindowSnps([]);
        setVariantWindowError(e?.message || 'Failed to load SNP window');
      } finally {
        if (requestSeq === hoverRequestSeqRef.current) {
          setVariantWindowLoading(false);
        }
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

  const activeTrackCount = useMemo(() => TRACK_ORDER.filter((key) => tracks[key]).length, [tracks]);

  const selectedCard = useMemo(() => evidenceCards[0] ?? null, [evidenceCards]);

  const selectedFeatureLocus = useMemo(() => {
    if (!selectedCard) {
      return null;
    }
    return extractFeatureLocus(selectedCard.fields);
  }, [selectedCard]);

  const selectedRefAlt = useMemo(() => {
    if (!selectedCard) {
      return { ref: null as string | null, alt: null as string | null };
    }
    return {
      ref: pickFieldValue(selectedCard.fields, ['ref', 'reference', 'reference allele']),
      alt: pickFieldValue(selectedCard.fields, ['alt', 'alternate', 'alternate allele']),
    };
  }, [selectedCard]);

  const variantSignalSummary = useMemo(() => {
    if (variantWindowSnps.length === 0) {
      return {
        strongest: null as VariantWindowSnp | null,
        meanScore: null as number | null,
      };
    }
    let strongest = variantWindowSnps[0];
    let total = 0;
    for (const row of variantWindowSnps) {
      total += row.score;
      if (row.score > strongest.score) {
        strongest = row;
      }
    }
    return {
      strongest,
      meanScore: total / variantWindowSnps.length,
    };
  }, [variantWindowSnps]);

  const selectedCausalRun = useMemo(
    () => causalRuns.find((run) => run.run_id === selectedCausalRunId) ?? null,
    [causalRuns, selectedCausalRunId]
  );

  const selectedRunVariant = useMemo(() => {
    if (!causalResult) {
      return null;
    }
    const targetLocus = selectedFeatureLocus ?? getBrowserLocus(browserRef.current);
    const rows = causalResult.variant_scores || [];
    if (rows.length === 0) {
      return null;
    }
    if (!targetLocus) {
      return rows[0];
    }
    let best: CausalVariantScore | null = null;
    let bestDistance = Number.POSITIVE_INFINITY;
    const center = Math.floor((targetLocus.start + targetLocus.end) / 2);
    for (const row of rows) {
      if (normalizeChromosome(row.chr) !== normalizeChromosome(targetLocus.chr)) {
        continue;
      }
      const distance = Math.abs(row.pos - center);
      if (distance < bestDistance) {
        bestDistance = distance;
        best = row;
      }
    }
    return best ?? rows[0];
  }, [causalResult, selectedFeatureLocus, currentLocus]);

  const selectedRunGene = useMemo(() => {
    if (!causalResult) {
      return null;
    }
    const rows = causalResult.gene_scores || [];
    if (rows.length === 0) {
      return null;
    }
    if (selectedGene) {
      const matched = rows.find((row) => row.gene.toUpperCase() === selectedGene.toUpperCase());
      if (matched) {
        return matched;
      }
    }
    if (selectedRunVariant?.gene) {
      const matched = rows.find((row) => row.gene.toUpperCase() === selectedRunVariant.gene!.toUpperCase());
      if (matched) {
        return matched;
      }
    }
    return rows[0];
  }, [causalResult, selectedGene, selectedRunVariant]);

  const visibleRange = useMemo(() => parseRegion(currentLocus || decodedRegion), [currentLocus, decodedRegion]);

  const visibleCausalVariants = useMemo(() => {
    if (!causalResult) {
      return [];
    }
    const span = Math.max(1, visibleRange.end - visibleRange.start);
    return (causalResult.variant_scores || [])
      .filter((row) => normalizeChromosome(row.chr) === normalizeChromosome(visibleRange.chr))
      .filter((row) => row.pos >= visibleRange.start && row.pos <= visibleRange.end)
      .map((row) => ({
        ...row,
        leftPct: ((row.pos - visibleRange.start) / span) * 100,
      }))
      .sort((a, b) => a.pos - b.pos);
  }, [causalResult, visibleRange]);

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
    if (!apiToken) {
      setCausalRuns([]);
      setSelectedCausalRunId('');
      setCausalResult(null);
      setCausalRunsError(null);
      setCausalResultError(null);
      setCausalRunsLoading(false);
      setCausalResultLoading(false);
      return;
    }

    let cancelled = false;
    const loadRuns = async () => {
      try {
        setCausalRunsLoading(true);
        setCausalRunsError(null);
        const response = await fetch(`${apiUrl}/causal/runs?limit=50`, {
          headers: {
            Authorization: `Bearer ${apiToken}`,
          },
        });
        if (cancelled) {
          return;
        }
        if (response.status === 401) {
          invalidateToken('Token expired or invalid. Please sign in again.');
          return;
        }
        if (!response.ok) {
          throw new Error(`Failed to load causal runs (${response.status})`);
        }
        const payload = await response.json().catch(() => ({}));
        const rows = Array.isArray(payload?.runs) ? (payload.runs as CausalRunRecord[]) : [];
        const filtered = rows.filter((run) => normalizeChromosome(run.params?.chr) === normalizeChromosome(chr));
        setCausalRuns(filtered);
        setSelectedCausalRunId((prev) => {
          if (prev && filtered.some((run) => run.run_id === prev)) {
            return prev;
          }
          return filtered[0]?.run_id || '';
        });
      } catch (e: any) {
        if (!cancelled) {
          setCausalRuns([]);
          setSelectedCausalRunId('');
          setCausalRunsError(e?.message || 'Failed to load causal runs');
        }
      } finally {
        if (!cancelled) {
          setCausalRunsLoading(false);
        }
      }
    };

    void loadRuns();
    return () => {
      cancelled = true;
    };
  }, [apiToken, apiUrl, chr, invalidateToken, causalRunsRefreshKey]);

  useEffect(() => {
    if (!apiToken || !selectedCausalRunId) {
      setCausalResult(null);
      setCausalResultError(null);
      setCausalResultLoading(false);
      return;
    }

    let cancelled = false;
    const loadResult = async () => {
      try {
        setCausalResultLoading(true);
        setCausalResultError(null);
        const response = await fetch(`${apiUrl}/causal/runs/${encodeURIComponent(selectedCausalRunId)}/result`, {
          headers: {
            Authorization: `Bearer ${apiToken}`,
          },
        });
        if (cancelled) {
          return;
        }
        if (response.status === 401) {
          invalidateToken('Token expired or invalid. Please sign in again.');
          return;
        }
        if (!response.ok) {
          throw new Error(`Failed to load causal result (${response.status})`);
        }
        const payload = await response.json().catch(() => ({}));
        setCausalResult((payload?.result as CausalResult) || null);
      } catch (e: any) {
        if (!cancelled) {
          setCausalResult(null);
          setCausalResultError(e?.message || 'Failed to load causal result');
        }
      } finally {
        if (!cancelled) {
          setCausalResultLoading(false);
        }
      }
    };

    void loadResult();
    return () => {
      cancelled = true;
    };
  }, [apiToken, apiUrl, invalidateToken, selectedCausalRunId]);

  useEffect(() => {
    if (!apiToken) {
      setWorkflowRuns([]);
      setSelectedWorkflowRunId('');
      setWorkflowRunDetail(null);
      setWorkflowRunsError(null);
      setWorkflowRunDetailError(null);
      setWorkflowRunsLoading(false);
      setWorkflowRunDetailLoading(false);
      return;
    }

    let cancelled = false;
    const loadWorkflowRuns = async () => {
      try {
        setWorkflowRunsLoading(true);
        setWorkflowRunsError(null);
        const response = await fetch(`${apiUrl}/workflows/runs`, {
          headers: {
            Authorization: `Bearer ${apiToken}`,
          },
        });
        if (cancelled) {
          return;
        }
        if (response.status === 401) {
          invalidateToken('Token expired or invalid. Please sign in again.');
          return;
        }
        if (!response.ok) {
          throw new Error(`Failed to load workflow runs (${response.status})`);
        }
        const payload = await response.json().catch(() => ({}));
        const rows = Array.isArray(payload?.runs) ? (payload.runs as WorkflowRunSummary[]) : [];
        setWorkflowRuns(rows);
        setSelectedWorkflowRunId((prev) => {
          if (prev && rows.some((run) => run.run_id === prev)) {
            return prev;
          }
          return rows[0]?.run_id || '';
        });
      } catch (e: any) {
        if (!cancelled) {
          setWorkflowRuns([]);
          setSelectedWorkflowRunId('');
          setWorkflowRunsError(e?.message || 'Failed to load workflow runs');
        }
      } finally {
        if (!cancelled) {
          setWorkflowRunsLoading(false);
        }
      }
    };

    void loadWorkflowRuns();
    return () => {
      cancelled = true;
    };
  }, [apiToken, apiUrl, invalidateToken, workflowRunsRefreshKey]);

  useEffect(() => {
    if (!apiToken || !selectedWorkflowRunId) {
      setWorkflowRunDetail(null);
      setWorkflowRunDetailError(null);
      setWorkflowRunDetailLoading(false);
      return;
    }

    let cancelled = false;
    const loadWorkflowRunDetail = async () => {
      try {
        setWorkflowRunDetailLoading(true);
        setWorkflowRunDetailError(null);
        const response = await fetch(`${apiUrl}/workflows/runs/${encodeURIComponent(selectedWorkflowRunId)}`, {
          headers: {
            Authorization: `Bearer ${apiToken}`,
          },
        });
        if (cancelled) {
          return;
        }
        if (response.status === 401) {
          invalidateToken('Token expired or invalid. Please sign in again.');
          return;
        }
        if (!response.ok) {
          throw new Error(`Failed to load workflow run (${response.status})`);
        }
        const payload = await response.json().catch(() => ({}));
        setWorkflowRunDetail((payload as WorkflowRunDetail) || null);
      } catch (e: any) {
        if (!cancelled) {
          setWorkflowRunDetail(null);
          setWorkflowRunDetailError(e?.message || 'Failed to load workflow run');
        }
      } finally {
        if (!cancelled) {
          setWorkflowRunDetailLoading(false);
        }
      }
    };

    void loadWorkflowRunDetail();
    return () => {
      cancelled = true;
    };
  }, [apiToken, apiUrl, invalidateToken, selectedWorkflowRunId]);

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
          showNavigation: false,
          showIdeogram: false,
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
    let lastWheelAt = 0;

    const onMouseMove = (event: MouseEvent) => {
      if (cancelled || browserRef.current !== browser) {
        return;
      }
      const now = Date.now();
      if (now - lastHoverAt < 60) {
        return;
      }
      lastHoverAt = now;

      const hovered = extractHoveredGeneContext(browser, event);
      if (!hovered) {
        hideVariantWindow();
        return;
      }
      const center = Math.max(1, Math.round(hovered.hoverPos));
      const snappedCenter = Math.max(1, Math.round(center / VARIANT_WINDOW_STEP_BP) * VARIANT_WINDOW_STEP_BP);
      const key = `${hovered.featureLocus.chr}:${snappedCenter}`;
      if (hoverWindowKeyRef.current === key) {
        setVariantWindowOpen(true);
        return;
      }
      hoverWindowKeyRef.current = key;
      void loadVariantWindow(hovered.featureLocus, snappedCenter);
    };

    const onMouseLeave = () => {
      hideVariantWindow();
    };

    const onWheel = (event: WheelEvent) => {
      if (cancelled || browserRef.current !== browser) {
        return;
      }
      event.preventDefault();
      const now = Date.now();
      if (now - lastWheelAt < 90) {
        return;
      }
      lastWheelAt = now;
      const factor = event.deltaY > 0 ? 1.8 : 0.55;
      void zoomBy(factor);
    };

    containerEl.addEventListener('mousemove', onMouseMove);
    containerEl.addEventListener('mouseleave', onMouseLeave);
    containerEl.addEventListener('wheel', onWheel, { passive: false });
    return () => {
      cancelled = true;
      containerEl.removeEventListener('mousemove', onMouseMove);
      containerEl.removeEventListener('mouseleave', onMouseLeave);
      containerEl.removeEventListener('wheel', onWheel);
    };
  }, [hideVariantWindow, loadVariantWindow, viewerClosed, viewerReady, zoomBy]);

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
      <header className="hero-shell">
        <div className="hero-copy">
          <p className="hero-kicker">Genome Workspace</p>
          <h1 className="page-title">Locus Analysis Dock</h1>
          <p className="hero-note">
            Custom shell around IGV with synchronized evidence and statistics panels for downstream scoring.
          </p>
        </div>
        <div className="hero-status">
          <div className="status-chip">
            <span>Focus</span>
            <strong>{currentLocus || decodedRegion}</strong>
          </div>
          <div className="status-chip">
            <span>Tracks</span>
            <strong>
              {activeTrackCount}/{TRACK_ORDER.length}
            </strong>
          </div>
          <div className={`status-chip ${apiToken ? 'status-ok' : 'status-warn'}`}>
            <span>Variants</span>
            <strong>{apiToken ? 'Authorized' : 'Token Required'}</strong>
          </div>
        </div>
      </header>

      <section className="command-grid">
        <div className="command-card auth-card">
          <div className="section-head">
            <h2>Access</h2>
            <p>Variants and linked PubMed metadata require an API token.</p>
          </div>
          {apiToken ? (
            <div className="auth-inline">
              <p className="auth-ok">API token loaded. Variants track is enabled.</p>
              <button type="button" className="auth-clear" onClick={clearToken}>
                Clear Token
              </button>
            </div>
          ) : (
            <div className="auth-form">
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
            </div>
          )}
        </div>

        <div className="command-card">
          <div className="section-head">
            <h2>Navigate</h2>
            <p>Drive the IGV frame from a controlled command bar instead of the default toolbar.</p>
          </div>
          <div className="viewer-nav">
            <input
              className="viewer-locus-input"
              value={navLocusInput}
              onChange={(event) => setNavLocusInput(event.target.value)}
              onKeyDown={(event) => {
                if (event.key === 'Enter') {
                  void goToLocus();
                }
              }}
              placeholder="chr:start-end"
              disabled={viewerClosed || !viewerReady}
            />
            <button
              type="button"
              className="viewer-nav-btn"
              onClick={() => void goToLocus()}
              disabled={viewerClosed || !viewerReady}
            >
              Go
            </button>
            <button
              type="button"
              className="viewer-nav-btn"
              onClick={() => void zoomBy(2)}
              disabled={viewerClosed || !viewerReady}
            >
              Zoom Out
            </button>
            <button
              type="button"
              className="viewer-nav-btn"
              onClick={() => void zoomBy(0.5)}
              disabled={viewerClosed || !viewerReady}
            >
              Zoom In
            </button>
            <button
              type="button"
              className="viewer-nav-btn"
              onClick={() => void goToLocus(decodedRegion)}
              disabled={viewerClosed || !viewerReady}
            >
              Reset
            </button>
          </div>
        </div>

        <div className="command-card">
          <div className="section-head">
            <h2>Causal Runs</h2>
            <p>Select a stored causal scoring run to sync gene and variant statistics with the locus view.</p>
          </div>
          {apiToken ? (
            <div className="run-selector-shell">
              <div className="run-selector-row">
                <select
                  className="run-select"
                  value={selectedCausalRunId}
                  onChange={(event) => setSelectedCausalRunId(event.target.value)}
                >
                  <option value="">No causal run selected</option>
                  {causalRuns.map((run) => (
                    <option key={run.run_id} value={run.run_id}>
                      {run.run_id} {run.project_id ? `- ${run.project_id}` : ''} ({run.params.start}-{run.params.end})
                    </option>
                  ))}
                </select>
                <button
                  type="button"
                  className="viewer-nav-btn"
                  onClick={() => setCausalRunsRefreshKey((value) => value + 1)}
                >
                  Refresh
                </button>
              </div>
              {causalRunsLoading ? <p className="placeholder">Loading causal runs...</p> : null}
              {causalRunsError ? <p className="literature-error">{causalRunsError}</p> : null}
              {selectedCausalRun ? (
                <div className="run-summary-chip">
                  <span>Region</span>
                  <strong>
                    {selectedCausalRun.params.chr}:{selectedCausalRun.params.start}-{selectedCausalRun.params.end}
                  </strong>
                </div>
              ) : (
                <p className="placeholder">No causal run selected for this chromosome.</p>
              )}
            </div>
          ) : (
            <p className="placeholder">Sign in to load saved causal scoring runs.</p>
          )}
        </div>

        <div className="command-card">
          <div className="section-head">
            <h2>Workflow Runs</h2>
            <p>Attach a generated Snakemake workflow run so the result dock can show execution context and artifacts.</p>
          </div>
          {apiToken ? (
            <div className="run-selector-shell">
              <div className="run-selector-row">
                <select
                  className="run-select"
                  value={selectedWorkflowRunId}
                  onChange={(event) => setSelectedWorkflowRunId(event.target.value)}
                >
                  <option value="">No workflow run selected</option>
                  {workflowRuns.map((run) => (
                    <option key={run.run_id} value={run.run_id}>
                      {run.run_id} - {run.workflow_id}
                    </option>
                  ))}
                </select>
                <button
                  type="button"
                  className="viewer-nav-btn"
                  onClick={() => setWorkflowRunsRefreshKey((value) => value + 1)}
                >
                  Refresh
                </button>
              </div>
              {workflowRunsLoading ? <p className="placeholder">Loading workflow runs...</p> : null}
              {workflowRunsError ? <p className="literature-error">{workflowRunsError}</p> : null}
              {workflowRunDetail ? (
                <div className="run-summary-chip">
                  <span>Workflow</span>
                  <strong>
                    {workflowRunDetail.summary.workflow_id} · {workflowRunDetail.summary.completed_runs}/
                    {workflowRunDetail.summary.submitted_runs} completed
                  </strong>
                </div>
              ) : (
                <p className="placeholder">No workflow run attached.</p>
              )}
            </div>
          ) : (
            <p className="placeholder">Sign in to load generated workflow runs.</p>
          )}
        </div>
      </section>

      {error ? (
        <div className="error">{error}</div>
      ) : (
        <section className={`workspace-shell ${viewerExpanded ? 'viewer-grid-expanded' : ''}`}>
          {viewerExpanded ? (
            <button type="button" className="viewer-overlay-close" onClick={() => setViewerExpanded(false)}>
              Exit Fullscreen
            </button>
          ) : null}
          <div className="viewer-column">
            <section className="track-ribbon">
              <div className="section-head">
                <h2>Tracks</h2>
                <p>Toggle genome layers and keep the shell synchronized with feature selection.</p>
              </div>
              <div className="toggle-bar">
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
              </div>
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

            <section className="overlay-rail-card">
              <div className="section-head">
                <h2>Causal Rail</h2>
                <p>
                  Overlaying {visibleCausalVariants.length} scored variants in {visibleRange.chr}:{visibleRange.start}-{visibleRange.end}
                </p>
              </div>
              {causalResultLoading ? (
                <p className="placeholder">Loading causal overlay...</p>
              ) : !causalResult ? (
                <p className="placeholder">Select a causal run to project score markers onto the active locus window.</p>
              ) : visibleCausalVariants.length === 0 ? (
                <p className="placeholder">No causal variants fall inside the current IGV window.</p>
              ) : (
                <div className="overlay-rail">
                  {visibleCausalVariants.map((row) => (
                    <button
                      key={`overlay-${row.chr}-${row.pos}-${row.alt}`}
                      type="button"
                      className={`overlay-marker ${selectedRunVariant?.pos === row.pos ? 'overlay-marker-active' : ''}`}
                      style={{
                        left: `${Math.max(0, Math.min(100, row.leftPct))}%`,
                        height: `${Math.max(18, Math.round(row.score * 72))}px`,
                      }}
                      title={`${row.chr}:${row.pos} score ${row.score.toFixed(3)}`}
                      onClick={() => void goToLocus(`${row.chr}:${Math.max(1, row.pos - 120)}-${row.pos + 120}`)}
                    >
                      <span>{row.pos}</span>
                    </button>
                  ))}
                </div>
              )}
            </section>

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

            <section className="stats-shell">
              <div className="metric-grid">
                <article className="metric-card">
                  <span>Selected Gene</span>
                  <strong>{selectedGene ?? 'N/A'}</strong>
                </article>
                <article className="metric-card">
                  <span>Feature Locus</span>
                  <strong>
                    {selectedFeatureLocus
                      ? `${selectedFeatureLocus.chr}:${selectedFeatureLocus.start}-${selectedFeatureLocus.end}`
                      : currentLocus || decodedRegion}
                  </strong>
                </article>
                <article className="metric-card">
                  <span>Hover Window SNPs</span>
                  <strong>{variantWindowSnps.length}</strong>
                </article>
                <article className="metric-card">
                  <span>Top Score</span>
                  <strong>{variantSignalSummary.strongest ? variantSignalSummary.strongest.score.toFixed(1) : 'N/A'}</strong>
                </article>
                <article className="metric-card">
                  <span>Mean Score</span>
                  <strong>{variantSignalSummary.meanScore !== null ? variantSignalSummary.meanScore.toFixed(1) : 'N/A'}</strong>
                </article>
                <article className="metric-card">
                  <span>Causal Gene Score</span>
                  <strong>{selectedRunGene ? selectedRunGene.score.toFixed(3) : 'N/A'}</strong>
                </article>
                <article className="metric-card">
                  <span>Causal Variant Score</span>
                  <strong>{selectedRunVariant ? selectedRunVariant.score.toFixed(3) : 'N/A'}</strong>
                </article>
              </div>
              <div className="signal-board">
                <div className="signal-board-head">
                  <div>
                    <h2>Signal Strip</h2>
                    <p>{variantWindowOpen ? variantWindowLabel : 'Move across a gene track to stream local variant statistics.'}</p>
                  </div>
                  {variantSignalSummary.strongest ? (
                    <div className="signal-highlight">
                      Peak {variantSignalSummary.strongest.chr}:{variantSignalSummary.strongest.pos}
                    </div>
                  ) : null}
                </div>
                {variantWindowLoading ? (
                  <p className="placeholder">Loading SNP scores...</p>
                ) : variantWindowError ? (
                  <p className="variant-window-error">{variantWindowError}</p>
                ) : variantWindowSnps.length === 0 ? (
                  <p className="placeholder">
                    No synchronized variant statistics yet. Hover a gene track or sign in to activate the variants stream.
                  </p>
                ) : (
                  <div className="variant-window-scroll integrated-scroll">
                    <div className="snp-chart-strip">
                      {variantWindowSnps.map((snp) => (
                        <div key={`bar-${snp.id}`} className="snp-bar-item" title={`${snp.chr}:${snp.pos}`}>
                          <span>{snp.pos}</span>
                          <div
                            className={`snp-bar ${variantSignalSummary.strongest?.id === snp.id ? 'snp-bar-peak' : ''}`}
                            style={{ height: `${Math.max(10, Math.round((snp.score / variantScoreMax) * 160))}px` }}
                          />
                          <span>{snp.score.toFixed(1)}</span>
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            </section>
          </div>

          <aside className="analysis-dock">
            <section className="dock-card">
              <div className="evidence-head">
                <div>
                  <h2>Snapshot</h2>
                  <p className="evidence-focus">Focus: {currentLocus || decodedRegion}</p>
                </div>
                {evidenceCards.length > 0 ? (
                  <button type="button" className="evidence-clear" onClick={() => setEvidenceCards([])}>
                    Clear Selection
                  </button>
                ) : null}
              </div>
              <div className="snapshot-grid">
                <div className="snapshot-cell">
                  <span>Track</span>
                  <strong>{selectedCard?.trackName ?? 'N/A'}</strong>
                </div>
                <div className="snapshot-cell">
                  <span>Source</span>
                  <strong>{selectedCard?.sourceLabel ?? 'Pending'}</strong>
                </div>
                <div className="snapshot-cell">
                  <span>REF / ALT</span>
                  <strong>{selectedRefAlt.ref && selectedRefAlt.alt ? `${selectedRefAlt.ref} -> ${selectedRefAlt.alt}` : 'N/A'}</strong>
                </div>
                <div className="snapshot-cell">
                  <span>Literature</span>
                  <strong>{literatureRecords.length}</strong>
                </div>
              </div>
            </section>

            <section className="dock-card dock-scroll">
              <div className="section-head">
                <h2>Feature Detail</h2>
                <p>Click a track feature to pin genomic fields, provenance, and zoom target.</p>
              </div>
              {evidenceCards.length === 0 ? (
                <p className="placeholder">
                  Hover on a gene or variant in IGV to preview signal, then click to pin details into the dock.
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
                </div>
              )}
            </section>

            <section className="dock-card">
              <div className="section-head">
                <h2>Causal Overlay</h2>
                <p>Run-linked scores synchronized with the active locus and selected feature.</p>
              </div>
              {causalResultLoading ? (
                <p className="placeholder">Loading causal run result...</p>
              ) : causalResultError ? (
                <p className="literature-error">{causalResultError}</p>
              ) : !causalResult ? (
                <p className="placeholder">No causal run loaded. Use the run selector above to attach a scoring result.</p>
              ) : (
                <div className="overlay-stack">
                  <div className="overlay-grid">
                    <div className="snapshot-cell">
                      <span>Run</span>
                      <strong>{selectedCausalRunId}</strong>
                    </div>
                    <div className="snapshot-cell">
                      <span>Region</span>
                      <strong>
                        {causalResult.region.chr}:{causalResult.region.start}-{causalResult.region.end}
                      </strong>
                    </div>
                    <div className="snapshot-cell">
                      <span>Gene Rank</span>
                      <strong>{selectedRunGene ? `#${selectedRunGene.rank}` : 'N/A'}</strong>
                    </div>
                    <div className="snapshot-cell">
                      <span>Lead Variant</span>
                      <strong>
                        {causalResult.lead_variant
                          ? `${causalResult.lead_variant.chr}:${causalResult.lead_variant.pos}`
                          : 'N/A'}
                      </strong>
                    </div>
                  </div>
                  <div className="overlay-grid overlay-grid-single">
                    <div className="snapshot-cell">
                      <span>Focused Gene</span>
                      <strong>
                        {selectedRunGene
                          ? `${selectedRunGene.gene} (${selectedRunGene.score.toFixed(3)})`
                          : 'No matching gene score'}
                      </strong>
                    </div>
                    <div className="snapshot-cell">
                      <span>Focused Variant</span>
                      <strong>
                        {selectedRunVariant
                          ? `${selectedRunVariant.chr}:${selectedRunVariant.pos} (${selectedRunVariant.score.toFixed(3)})`
                          : 'No matching variant score'}
                      </strong>
                    </div>
                  </div>
                  <div className="overlay-rank-list">
                    {(causalResult.gene_scores || []).slice(0, 5).map((row) => (
                      <div key={`gene-score-${row.gene}`} className="overlay-rank-item">
                        <span>#{row.rank}</span>
                        <strong>{row.gene}</strong>
                        <em>{row.score.toFixed(3)}</em>
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </section>

            <section className="dock-card">
              <div className="section-head">
                <h2>Workflow Context</h2>
                <p>Execution summary and generated artifacts from the attached Snakemake workflow run.</p>
              </div>
              {workflowRunDetailLoading ? (
                <p className="placeholder">Loading workflow run detail...</p>
              ) : workflowRunDetailError ? (
                <p className="literature-error">{workflowRunDetailError}</p>
              ) : !workflowRunDetail ? (
                <p className="placeholder">No workflow run selected.</p>
              ) : (
                <div className="overlay-stack">
                  <div className="overlay-grid">
                    <div className="snapshot-cell">
                      <span>Workflow ID</span>
                      <strong>{workflowRunDetail.summary.workflow_id}</strong>
                    </div>
                    <div className="snapshot-cell">
                      <span>Run ID</span>
                      <strong>{workflowRunDetail.summary.run_id}</strong>
                    </div>
                    <div className="snapshot-cell">
                      <span>Completed</span>
                      <strong>
                        {workflowRunDetail.summary.completed_runs}/{workflowRunDetail.summary.submitted_runs}
                      </strong>
                    </div>
                    <div className="snapshot-cell">
                      <span>Duration</span>
                      <strong>{workflowRunDetail.summary.duration_ms} ms</strong>
                    </div>
                  </div>
                  <div className="overlay-rank-list">
                    {workflowRunDetail.results.slice(0, 3).map((row) => (
                      <div key={`wf-result-${row.index}`} className="overlay-rank-item">
                        <span>Task {row.index}</span>
                        <strong>{row.engine}</strong>
                        <em>{row.status}</em>
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </section>

            <section className="dock-card">
              <div className="section-head">
                <h2>PubMed</h2>
                <p>
                  Gene: <strong>{selectedGene ?? 'N/A'}</strong>
                </p>
              </div>
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
                      <a href={`https://pubmed.ncbi.nlm.nih.gov/${record.pmid}/`} target="_blank" rel="noreferrer">
                        PMID {record.pmid}
                      </a>
                      <p>{record.title}</p>
                      <span>{record.year ?? 'Year unknown'}</span>
                    </article>
                  ))}
                </div>
              )}
            </section>
          </aside>
        </section>
      )}
      <style jsx>{`
        .locus-page {
          --paper: #edf3ef;
          --paper-soft: #ffffff;
          --ink: #0f172a;
          --ink-muted: #475569;
          --panel: #ffffff;
          --line: #d2dae4;
          --accent: #0f766e;
          --accent-soft: #7dc7c1;
          --deep: #1e293b;
          min-height: 100vh;
          padding: 26px;
          font-family: 'IBM Plex Sans', 'Segoe UI', sans-serif;
          line-height: 1.35;
          color: var(--ink);
          background:
            radial-gradient(circle at top left, rgba(15, 118, 110, 0.12), transparent 24%),
            radial-gradient(circle at top right, rgba(194, 120, 39, 0.12), transparent 20%),
            #edf3ef;
          animation: pageIn 360ms ease-out;
          overflow-x: hidden;
        }
        .hero-shell,
        .command-grid,
        .workspace-shell {
          position: relative;
          z-index: 1;
        }
        .hero-shell {
          display: grid;
          grid-template-columns: minmax(0, 1.6fr) minmax(280px, 0.8fr);
          gap: 16px;
          margin-bottom: 16px;
          padding: 20px;
          border: 1px solid var(--line);
          border-radius: 20px;
          background: linear-gradient(135deg, rgba(255, 255, 255, 0.96), rgba(242, 247, 245, 0.98));
        }
        .hero-kicker {
          margin: 0 0 8px;
          font-size: 11px;
          letter-spacing: 0.2em;
          text-transform: uppercase;
          font-weight: 800;
          color: var(--accent);
        }
        .page-title {
          margin: 0;
          font-family: 'Space Grotesk', 'IBM Plex Sans', 'Segoe UI', sans-serif;
          font-size: clamp(30px, 4vw, 52px);
          line-height: 1.02;
          letter-spacing: -0.03em;
        }
        .hero-note {
          margin: 10px 0 0;
          max-width: 720px;
          color: var(--ink-muted);
          font-size: 14px;
        }
        .hero-status {
          display: grid;
          gap: 10px;
          align-content: start;
        }
        .status-chip {
          display: grid;
          gap: 2px;
          padding: 12px 14px;
          border-radius: 14px;
          border: 1px solid #d7e2de;
          background: rgba(255, 255, 255, 0.88);
        }
        .status-chip span {
          font-size: 11px;
          letter-spacing: 0.12em;
          text-transform: uppercase;
          color: #5a6c74;
          font-weight: 700;
        }
        .status-chip strong {
          font-size: 14px;
          color: #12212a;
        }
        .status-ok {
          border-color: #87bca9;
          background: #effaf5;
        }
        .status-warn {
          border-color: #dfb294;
          background: #fff4eb;
        }
        .command-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
          gap: 14px;
          margin-bottom: 16px;
        }
        .command-card,
        .track-ribbon,
        .stats-shell,
        .dock-card {
          padding: 16px;
          border: 1px solid var(--line);
          border-radius: 18px;
          background: rgba(255, 255, 255, 0.95);
        }
        .section-head h2 {
          margin: 0;
          font-size: 22px;
          color: var(--deep);
        }
        .section-head p {
          margin: 4px 0 0;
          font-size: 13px;
          color: var(--ink-muted);
        }
        .auth-card {
          display: grid;
          gap: 12px;
        }
        .run-selector-shell {
          display: grid;
          gap: 10px;
        }
        .run-selector-row {
          display: flex;
          gap: 8px;
          align-items: center;
        }
        .run-select {
          min-height: 46px;
          flex: 1;
          border-radius: 10px;
          border: 1px solid #cbd5e1;
          background: #ffffff;
          color: #0f172a;
          padding: 8px 12px;
        }
        .run-summary-chip {
          display: grid;
          gap: 2px;
          padding: 10px 12px;
          border-radius: 12px;
          border: 1px solid #d8e4de;
          background: #f7fbf9;
        }
        .run-summary-chip span {
          font-size: 11px;
          letter-spacing: 0.12em;
          text-transform: uppercase;
          color: #62747d;
          font-weight: 700;
        }
        .run-summary-chip strong {
          font-size: 13px;
          color: #10222c;
          overflow-wrap: anywhere;
        }
        .auth-inline,
        .auth-form {
          display: flex;
          flex-wrap: wrap;
          gap: 8px;
        }
        .viewer-nav {
          display: flex;
          align-items: center;
          gap: 8px;
          flex-wrap: wrap;
        }
        .viewer-locus-input {
          min-height: 46px;
          min-width: min(420px, 70vw);
          border-radius: 10px;
          border: 1px solid #cbd5e1;
          background: #ffffff;
          color: #0f172a;
          padding: 8px 12px;
        }
        .viewer-nav-btn {
          min-height: 46px;
          border-radius: 10px;
          border: 1px solid #cbd5e1;
          background: #ffffff;
          color: #0f172a;
          padding: 8px 12px;
          font-weight: 700;
          letter-spacing: 0.02em;
          text-transform: uppercase;
          cursor: pointer;
        }
        .viewer-nav-btn:disabled {
          opacity: 0.55;
          cursor: default;
        }
        .viewer-toolbar-actions {
          display: flex;
          gap: 8px;
          align-items: center;
          flex-wrap: wrap;
        }
        .viewer-toggle,
        .viewer-reopen,
        .viewer-expand {
          min-height: 56px;
          border-radius: 10px;
          border: 1px solid #c7d2e1;
          background: #ffffff;
          color: #0f172a;
          font-weight: 700;
          letter-spacing: 0.02em;
          text-transform: uppercase;
          padding: 8px 16px;
          cursor: pointer;
        }
        .viewer-reopen {
          background: #0f6cbd;
          border-color: #0f6cbd;
          color: #ffffff;
        }
        .viewer-expand {
          background: #1e293b;
          border-color: #1e293b;
          color: #ffffff;
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
          color: #14532d;
        }
        .auth-warn {
          color: #9a3412;
        }
        .auth-input {
          min-height: 56px;
          min-width: 150px;
          border: 1px solid #cbd5e1;
          border-radius: 8px;
          padding: 8px 12px;
          background: #ffffff;
          color: #0f172a;
        }
        .auth-input:focus {
          outline: 2px solid #0f6cbd;
          outline-offset: 1px;
        }
        .auth-login,
        .auth-clear {
          min-height: 56px;
          border-radius: 8px;
          border: 1px solid #0f6cbd;
          background: #0f6cbd;
          color: #ffffff;
          font-weight: 700;
          letter-spacing: 0.02em;
          text-transform: uppercase;
          padding: 8px 14px;
          cursor: pointer;
        }
        .auth-login:disabled {
          opacity: 0.7;
          cursor: default;
        }
        .auth-clear {
          background: #1e293b;
          border-color: #1e293b;
        }
        .auth-error {
          color: #9a2b23;
          font-size: 12px;
          font-weight: 700;
        }
        .workspace-shell {
          display: grid;
          grid-template-columns: minmax(0, 1.45fr) minmax(320px, 0.75fr);
          gap: 18px;
          align-items: start;
        }
        .viewer-column {
          display: grid;
          gap: 14px;
        }
        .overlay-rail-card {
          padding: 16px;
          border: 1px solid var(--line);
          border-radius: 18px;
          background: rgba(255, 255, 255, 0.95);
          display: grid;
          gap: 12px;
        }
        .overlay-rail {
          position: relative;
          height: 118px;
          border-radius: 16px;
          border: 1px solid #dbe5e1;
          background:
            linear-gradient(180deg, rgba(15, 118, 110, 0.06), rgba(15, 118, 110, 0)),
            linear-gradient(90deg, rgba(15, 118, 110, 0.04) 1px, transparent 1px);
          background-size: 100% 100%, 36px 100%;
          overflow: hidden;
        }
        .overlay-marker {
          position: absolute;
          bottom: 0;
          width: 14px;
          margin-left: -7px;
          border: none;
          border-radius: 999px 999px 4px 4px;
          background: linear-gradient(180deg, #d4833f, #0f766e);
          cursor: pointer;
          display: flex;
          align-items: flex-start;
          justify-content: center;
          padding: 4px 0 0;
          transition: transform 140ms ease, opacity 140ms ease;
        }
        .overlay-marker:hover {
          transform: translateY(-2px);
        }
        .overlay-marker span {
          writing-mode: vertical-rl;
          transform: rotate(180deg);
          font-size: 9px;
          line-height: 1;
          color: rgba(255, 255, 255, 0.9);
          font-weight: 700;
          letter-spacing: 0.04em;
        }
        .overlay-marker-active {
          background: linear-gradient(180deg, #ff9b4b, #0b5d57);
          box-shadow: inset 0 0 0 2px rgba(255, 255, 255, 0.4);
        }
        .track-ribbon {
          display: grid;
          gap: 12px;
        }
        .toggle-item {
          display: flex;
          align-items: center;
          gap: 9px;
          font-weight: 700;
          text-transform: uppercase;
          letter-spacing: 0.06em;
          color: #1e293b;
          border: 1px solid #d2dae4;
          border-radius: 999px;
          padding: 10px 16px;
          background: #ffffff;
          transition: border-color 160ms ease, background-color 160ms ease;
        }
        .toggle-item input:disabled {
          opacity: 0.5;
        }
        .toggle-item:hover {
          border-color: #0f6cbd;
          background: #f8fbff;
        }
        .toggle-item input {
          width: 24px;
          height: 24px;
          accent-color: #0f6cbd;
        }
        .error {
          color: #7f1d1d;
          font-size: 14px;
          font-weight: 700;
          padding: 12px 14px;
          border: 1px solid #fda4af;
          border-radius: 10px;
          background: #fff1f2;
        }
        .viewer-grid-expanded {
          position: fixed;
          inset: 10px;
          z-index: 120;
          background: rgba(237, 243, 239, 0.99);
          border: 1px solid #d2dae4;
          border-radius: 16px;
          padding: 12px;
          grid-template-columns: 1fr;
          overflow: auto;
        }
        .viewer-overlay-close {
          position: sticky;
          top: 0;
          margin-left: auto;
          z-index: 2;
          min-height: 44px;
          border-radius: 8px;
          border: 1px solid #1e293b;
          background: #1e293b;
          color: #ffffff;
          padding: 6px 12px;
          cursor: pointer;
        }
        .viewer-grid-expanded .analysis-dock {
          display: none;
        }
        .viewer-grid-expanded .viewer {
          min-height: calc(100vh - 56px);
        }
        .viewer {
          min-height: 620px;
          border: 1px solid #d2dae4;
          border-radius: 18px;
          background: var(--panel);
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
        .stats-shell {
          display: grid;
          gap: 14px;
        }
        .metric-grid {
          display: grid;
          grid-template-columns: repeat(5, minmax(0, 1fr));
          gap: 10px;
        }
        .metric-card {
          display: grid;
          gap: 6px;
          padding: 14px;
          border-radius: 14px;
          background: linear-gradient(180deg, #ffffff, #f6fbf9);
          border: 1px solid #d9e4df;
        }
        .metric-card span {
          font-size: 11px;
          letter-spacing: 0.12em;
          text-transform: uppercase;
          color: #5f7078;
          font-weight: 700;
        }
        .metric-card strong {
          font-size: 16px;
          color: #10222c;
          overflow-wrap: anywhere;
        }
        .signal-board {
          border: 1px solid #d9e4df;
          border-radius: 18px;
          background: linear-gradient(180deg, #fbfffd, #f3f9f6);
          padding: 16px;
        }
        .signal-board-head {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 10px;
          margin-bottom: 12px;
        }
        .signal-board-head h2 {
          margin: 0;
          font-size: 22px;
          color: var(--deep);
        }
        .signal-board-head p {
          margin: 4px 0 0;
          font-size: 13px;
          color: #64748b;
        }
        .signal-highlight {
          border-radius: 999px;
          padding: 7px 12px;
          border: 1px solid #c8ddd7;
          background: #ffffff;
          color: #12564e;
          font-size: 12px;
          font-weight: 700;
        }
        .analysis-dock {
          display: grid;
          gap: 14px;
          min-height: 0;
        }
        .dock-card {
          display: grid;
          gap: 12px;
        }
        .dock-scroll {
          min-height: 420px;
        }
        .evidence-head {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 10px;
        }
        .evidence-clear {
          min-height: 52px;
          border-radius: 8px;
          border: 1px solid #1e293b;
          background: #1e293b;
          color: #ffffff;
          padding: 8px 14px;
          font-weight: 700;
          letter-spacing: 0.04em;
          text-transform: uppercase;
          cursor: pointer;
        }
        .evidence-head h2 {
          margin: 0;
          font-size: 22px;
          color: var(--deep);
        }
        .evidence-focus {
          margin: 4px 0 0;
          font-size: 12px;
          color: #64748b;
          font-weight: 700;
          letter-spacing: 0.03em;
        }
        .snapshot-grid {
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 10px;
        }
        .overlay-stack {
          display: grid;
          gap: 10px;
        }
        .overlay-grid {
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 10px;
        }
        .overlay-grid-single {
          grid-template-columns: 1fr;
        }
        .overlay-rank-list {
          display: grid;
          gap: 8px;
        }
        .overlay-rank-item {
          display: grid;
          grid-template-columns: auto 1fr auto;
          gap: 10px;
          align-items: center;
          padding: 10px 12px;
          border-radius: 12px;
          border: 1px solid #dbe5e1;
          background: #f8fcfa;
        }
        .overlay-rank-item span {
          font-size: 12px;
          color: #5f7179;
          font-weight: 700;
        }
        .overlay-rank-item strong {
          font-size: 14px;
          color: #10222c;
        }
        .overlay-rank-item em {
          font-style: normal;
          font-size: 13px;
          color: #0f766e;
          font-weight: 700;
        }
        .snapshot-cell {
          display: grid;
          gap: 5px;
          padding: 12px;
          border-radius: 12px;
          border: 1px solid #dbe5e1;
          background: #f8fcfa;
        }
        .snapshot-cell span {
          font-size: 11px;
          letter-spacing: 0.12em;
          text-transform: uppercase;
          color: #6a7a81;
          font-weight: 700;
        }
        .snapshot-cell strong {
          font-size: 14px;
          color: #10222c;
          overflow-wrap: anywhere;
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
          background: #eef2f7;
          border-radius: 999px;
        }
        .evidence-list::-webkit-scrollbar-thumb {
          background: #94a3b8;
          border-radius: 999px;
          border: 2px solid #eef2f7;
        }
        .evidence-list::-webkit-scrollbar-thumb:hover {
          background: #64748b;
        }
        .evidence-card {
          border: 1px solid #dbe2ea;
          border-radius: 12px;
          background: #ffffff;
          padding: 10px;
        }
        .evidence-meta {
          margin-bottom: 14px;
          padding: 12px;
          background: #f8fafc;
          border-radius: 10px;
          border: 1px solid #dbe2ea;
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
          background: #e0f2fe;
          color: #0c4a6e;
          border-color: #93c5fd;
        }
        .pill-infer {
          background: #eef2ff;
          color: #3730a3;
          border-color: #c7d2fe;
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
          border: 1px solid #dbe2ea;
          border-radius: 10px;
          padding: 10px;
          background: #ffffff;
        }
        .field-name {
          display: block;
          font-size: 11px;
          color: #64748b;
          margin-bottom: 5px;
          text-transform: uppercase;
          letter-spacing: 0.07em;
          font-weight: 700;
        }
        .field-value {
          display: block;
          font-size: 13px;
          color: #0f172a;
          overflow-wrap: anywhere;
        }
        .placeholder {
          margin: 0;
          color: #64748b;
          font-size: 13px;
          line-height: 1.5;
        }
        .literature-card {
          border: 1px solid #dbe2ea;
          border-radius: 12px;
          background: #ffffff;
          padding: 12px;
        }
        .literature-card h3 {
          margin: 0 0 6px;
          font-size: 15px;
          font-weight: 800;
          letter-spacing: 0.04em;
          text-transform: uppercase;
          color: #1e293b;
        }
        .literature-subtitle {
          margin: 0 0 10px;
          font-size: 12px;
          color: #64748b;
        }
        .literature-list {
          display: flex;
          flex-direction: column;
          gap: 10px;
        }
        .literature-item {
          border: 1px solid #dbe2ea;
          border-radius: 10px;
          padding: 10px;
          background: #ffffff;
        }
        .literature-item a {
          color: #0f6cbd;
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
          color: #0f172a;
          line-height: 1.4;
        }
        .literature-item span {
          font-size: 12px;
          color: #64748b;
          font-weight: 700;
        }
        .literature-error {
          margin: 0;
          color: #9a2b23;
          font-size: 13px;
          font-weight: 700;
        }
        .variant-window-error {
          margin: 0;
          padding: 10px 0;
          font-size: 13px;
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
        .integrated-scroll {
          padding: 6px 0 2px;
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
          background: linear-gradient(180deg, #c46d2e, #0f766e);
        }
        .snp-bar-peak {
          background: linear-gradient(180deg, #f18b3a, #0b5d57);
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
          .hero-shell,
          .command-grid,
          .workspace-shell {
            grid-template-columns: 1fr;
          }
          .metric-grid {
            grid-template-columns: repeat(2, minmax(0, 1fr));
          }
          .snapshot-grid {
            grid-template-columns: 1fr;
          }
          .overlay-grid {
            grid-template-columns: 1fr;
          }
          .run-selector-row {
            flex-direction: column;
            align-items: stretch;
          }
          .auth-input {
            min-width: 130px;
          }
          .viewer-locus-input {
            min-width: 100%;
          }
          .viewer-grid-expanded {
            inset: 0;
            border-radius: 0;
            padding: 8px;
          }
          .viewer,
          .dock-scroll {
            min-height: 420px;
          }
        }
      `}</style>
      <style jsx global>{`
        /* IGV minimal modern theme (flat, no hover shadow) */
        .igv-root-div,
        .igv-container,
        .igv-column-container {
          font-family: 'IBM Plex Sans', 'Segoe UI', sans-serif !important;
        }
        .igv-navbar {
          display: none !important;
        }
        .igv-navbar .igv-current-genome {
          color: #e2e8f0 !important;
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
          border-color: #cbd5e1 !important;
          border-width: 1px !important;
          font-size: 14pt !important;
          color: #0f172a !important;
          background-color: #ffffff !important;
        }
        .igv-navbar .igv-navbar-button {
          min-height: 44px !important;
          line-height: 42px !important;
          border-radius: 8px !important;
          border-width: 1px !important;
          border-color: #cbd5e1 !important;
          background-color: #ffffff !important;
          color: #0f172a !important;
          font-size: 14pt !important;
          font-weight: 700 !important;
          letter-spacing: 0.03em !important;
          padding-left: 10px !important;
          padding-right: 10px !important;
        }
        .igv-navbar .igv-navbar-button.igv-navbar-button-clicked {
          background-color: #0f6cbd !important;
          border-color: #0f6cbd !important;
          color: #ffffff !important;
        }
        .igv-navbar .igv-navbar-right-container .igv-zoom-widget div:first-child,
        .igv-navbar .igv-navbar-right-container .igv-zoom-widget div:last-child {
          width: 44px !important;
          height: 44px !important;
          color: #cbd5e1 !important;
        }
        .igv-navbar .igv-navbar-right-container .igv-zoom-widget svg {
          width: 30px !important;
          height: 30px !important;
        }
        .igv-gear-menu-column > div > div {
          width: 24px !important;
          height: 24px !important;
          color: #475569 !important;
        }
        .igv-gear-menu-column > div > div > svg {
          width: 24px !important;
          height: 24px !important;
        }
        .igv-track-drag-column > .igv-track-drag-handle {
          background-color: #0f6cbd !important;
          border-radius: 6px !important;
        }
        .igv-track-label {
          font-size: 14pt !important;
          font-weight: 700 !important;
          border-width: 1px !important;
          border-color: #cbd5e1 !important;
          background-color: #ffffff !important;
          color: #0f172a !important;
          padding: 2px 8px !important;
          border-radius: 8px !important;
        }
        .igv-menu-popup {
          border-width: 1px !important;
          border-color: #cbd5e1 !important;
          box-shadow: none !important;
        }
        .igv-menu-popup-header {
          background-color: #f1f5f9 !important;
        }
        .igv-axis-column,
        .igv-scrollbar-column,
        .igv-track-drag-column,
        .igv-gear-menu-column {
          background-color: #f8fafc !important;
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
          border-color: #94a3b8 !important;
          background-color: rgba(148, 163, 184, 0.4) !important;
        }
        .igv-scrollbar-column > div > div:hover {
          background-color: rgba(148, 163, 184, 0.4) !important;
        }
        .igv-root-div *:hover {
          box-shadow: none !important;
        }
      `}</style>
    </main>
  );
}
