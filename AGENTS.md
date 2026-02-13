# AGENTS.md — AD Multi-Omics Locus Evidence Platform (Intranet-first SaaS)

## Mission
Build an intranet-deployable (offline-capable) SaaS-style web app for Alzheimer’s disease (AD) multi-omics evidence integration:
- Upload genomics/transcriptomics/proteomics data
- Visualize loci on chromosomes (embedded genome viewer)
- Rank candidate loci/genes with explainable baseline scoring
- Generate evidence-backed research directions (offline fallback; online adds PubMed + LLM)
- Support plugin-based methods (new AI/statistical models) via Dockerized runners
- Produce auditable, reproducible, exportable reports for biotech clients

## Non-goals (MVP)
- Single-cell pipelines, raw FASTQ/BAM/CRAM processing
- Full fine-mapping/coloc/MR (unless explicitly added as a plugin later)
- Fully automated “no-human-review” deployment to production

## Operating modes
- Offline mode: MUST work without external network (no PubMed, no cloud LLM)
- Online mode: Optional enhancements (PubMed metadata retrieval, cloud LLM summaries)

## Hard constraints
- Evidence integrity: Never fabricate citations. If an output is an inference, label it clearly.
- Reproducibility: Every analysis run stores params, tool versions, input hashes.
- Security: Project-based RBAC + audit log for upload/run/export events.
- One task per PR. PR must pass CI (lint + tests) before merge.

## Repository layout (expected)
- `web/`      Next.js UI (Locus Explorer, Dashboard, Upload)
- `api/`      FastAPI services (auth, RBAC, upload, runs, evidence)
- `workers/`  background tasks (parsing, scoring, report, plugins)
- `infra/`    docker/helm/compose, env templates, db migrations

## Local development
### Required commands (update to match repo)
- Start: `docker compose up --build`
- Web lint/test: `cd web && npm run lint && npm test`
- API lint/test: `cd api && ruff check . && pytest -q`
- Workers test: `cd workers && pytest -q`

### Definition of Done (DoD)
- Feature implemented with tests
- Works in offline mode (if relevant)
- Documentation updated (README or relevant docs)
- PR includes verification steps + screenshots for UI changes

---

# Roadmap (MVP)

## EPIC 0 — Scaffolding & CI
- [ ] T0.1 Monorepo scaffold + docker compose (web+api+postgres+minio+worker)
  - Acceptance:
    - `docker compose up` brings up services and UI loads
- [ ] T0.2 CI gate (lint + unit tests + container build)
  - Acceptance:
    - PR cannot merge unless CI green

## EPIC 1 — Auth / RBAC / Audit
- [ ] T1.1 Auth (local accounts) + stub for OIDC
  - Acceptance:
    - Login works; API protected by tokens
- [ ] T1.2 Project RBAC (admin/analyst/viewer)
  - Acceptance:
    - Viewer cannot upload or run analyses
- [ ] T1.3 Audit log
  - Acceptance:
    - Upload/run/export recorded and queryable

## EPIC 2 — Data ingest
- [ ] T2.1 Upload API + object storage (S3 interface; MinIO in intranet)
  - Acceptance:
    - Upload creates dataset record with hash + storage URI
- [ ] T2.2 VCF parser (MVP fields) + indexing
  - Acceptance:
    - Variants list can be queried by chr; basic stats visible
- [ ] T2.3 Transcriptomics/Proteomics diff-table parser
  - Acceptance:
    - Summary stats shown; standardized schema stored

## EPIC 3 — Baseline scoring (explainable)
- [ ] T3.1 Variant→Gene mapping (MVP: nearest gene / annotation-based)
  - Acceptance:
    - Each top locus links to candidate gene(s)
- [ ] T3.2 Evidence join (omics + genomics) + rank aggregation
  - Acceptance:
    - Top genes/loci list generated with feature breakdown
- [ ] T3.3 Run reproducibility metadata
  - Acceptance:
    - Run stores params, tool versions, input hashes; rerun stable

## EPIC 4 — Locus Explorer (genome viewer)
- [x] T4.1 Embed genome viewer (choose igv.js first)
  - Acceptance:
    - `/locus/chr:start-end` shows tracks UI
- [x] T4.2 Tracks loader (genes, variants, scores)
  - Acceptance:
    - Tracks toggle; clicking feature updates evidence panel
- [x] T4.3 Evidence panel cards
  - Acceptance:
    - Cards show sources or clearly labeled inferences
- [ ] T4.x Host gene annotation track
  - Acceptance:
    - Genes track served from our own MinIO/S3 (or signed URL), not public igv.org; CORS configured and loads in `/locus/*`

## EPIC 5 — Literature + Research direction
- [x] T5.1 PubMed provider (online mode only)
  - Acceptance:
    - Top genes retrieve PubMed metadata; store pmid + title + year
- [x] T5.2 Offline fallback research-direction generator (template/rules)
  - Acceptance:
    - Offline mode still produces hypotheses/experiments section
- [ ] T5.3 Optional cloud LLM summarizer (must cite evidence bundle)
  - Acceptance:
    - Output contains citations list; no fabricated references

## EPIC 6 — Method Registry (plugins)
- [ ] T6.1 Plugin manifest spec (input/output schema, version, resources)
  - Acceptance:
    - Register method via JSON; validate schema
- [ ] T6.2 Plugin runner (Docker image execution)
  - Acceptance:
    - Baseline plugin runs and writes ranked loci/genes back to DB
- [ ] T6.3 UI method management
  - Acceptance:
    - Enable/disable methods; compare runs

## EPIC 7 — Deliverables
- [ ] T7.1 Report export (Markdown/HTML; PDF later)
  - Acceptance:
    - Export includes top loci/genes, evidence tables, research direction, run metadata
- [ ] T7.2 Demo dataset + scripted demo
  - Acceptance:
    - Fresh install can demo end-to-end in ≤10 minutes
