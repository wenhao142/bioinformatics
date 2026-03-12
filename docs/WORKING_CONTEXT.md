# Working Context

## Current State

- Product direction is a focused internal bioinformatics workflow console, not a general workflow SaaS.
- Backend remains `FastAPI + Snakemake`.
- `workflow-builder` is a constrained WGS builder with fixed tool blocks.
- `workflow-builder` now follows the main user flow:
  - upload raw
  - load template
  - save
  - run
  - result
- The builder upload area now keeps only `Raw` upload.
- Raw upload now supports adding files across multiple picker actions without replacing previous selections.
- The builder now shows the selected raw-file list before upload and allows removing individual files.
- Builder raw upload no longer sends `Content-Type: application/json` for multipart requests; this fixed the previous `422` upload failure.
- Empty builder canvas now shows a neutral start state instead of reporting hidden parameter-sweep validation errors.
- Builder `Save` now explains that workflows are stored in the API workflow store under `workflow_id`.
- Builder `Result` now shows a readable run summary first and only shows raw JSON on demand.
- Added `docs/flowcharts/` with local-only and future connected execution flowcharts derived from `AGENTS.md`.
- Added `docs/CLEANUP_CANDIDATES.md` to track redundant or transitional files safely before deletion.
- Removed generated `web/tsconfig.tsbuildinfo`.

## Newly Added Durable Rule

- This project should not keep saved workflows or run results only in transient API memory.
- Target persistence model:
  - local filesystem for workflow/result artifacts
  - PostgreSQL for workflow definitions, run metadata, and queryable history
- Future cloud execution should dispatch from locally persisted state, not replace local persistence.
- The main console page was simplified to only keep currently used inputs and latest results.
- Login page was simplified to focus on authentication only.
- Frontend now includes a `shadcn/ui` base layer:
  - `web/components.json`
  - `web/components/ui/*`
  - `web/lib/utils.ts`
  - Tailwind v4 via `web/app/globals.css` and `web/postcss.config.mjs`
- The current UI refresh follows local `impeccable-style-universal` guidance:
  - quieter surfaces
  - shorter copy
  - clearer hierarchy
  - desktop/mobile adaptation first
- `AGENTS.md` now explicitly treats `impeccable-style-universal/` as the primary design-guidance source for future web UI work.
- A Tailwind global-style regression was fixed:
  - `web/app/globals.css` no longer uses the broken `@apply` pattern that caused the page to render like mostly unstyled text
  - `web` was rebuilt after the fix
- The frontend styling pipeline was repaired:
  - Tailwind v4 setup was removed because utilities were not being generated in this `Next.js 14.1` environment
  - The project now uses a stable Tailwind v3 setup with `tailwind.config.ts` and `postcss.config.js`
  - Verified generated CSS now contains actual utility classes such as `min-h-screen`, `rounded-xl`, and `bg-card`
- Frontend API base resolution is now centralized in `web/lib/api-base.ts`.
- `Console`, `Login`, `Workflow Builder`, and `Locus` now all use the same API resolver.
- The resolver now rewrites `localhost` API targets to the current browser hostname when the page is opened from a non-loopback host.
- Frontend session handling is now centralized in `web/lib/session.ts`.
- `Console` and `Workflow Builder` now block upload/run actions until `/auth/me` finishes validating the stored token.
- If the browser still has a stale token, the UI clears it and redirects back to `/login` instead of letting upload/run fail immediately with `Invalid token`.
- `AGENTS.md` now includes a stricter reporting rule:
  - no success claim before local verification passes
  - verification result must be reported before the change summary
- `AGENTS.md` now explicitly states that:
  - upload success does not mean workflow compatibility
  - invalid or unbound workflow inputs must fail before execution
  - the system must not return successful placeholder results for invalid inputs
- `workflow-builder` now uses the same `AppShell + shadcn/ui` page structure as the console instead of the older custom topbar/hero shell.
- The builder layout is now organized as:
  - summary cards
  - raw upload
  - short flow guidance
  - blocks / canvas / workflow panel / validation / result
- The builder no longer relies on the old click-to-enter-link-mode interaction for basic use.
- The workflow panel now exposes `Load Template`, `Auto-connect`, `Save`, `Run`, and `Saved workflows` in a clearer order.
- Workflow persistence is now wired through `_STORE` instead of the old direct `_WORKFLOWS/_WORKFLOW_RUNS` endpoint path.
- Workflow import now writes a local definition JSON file under `results/workflows/definitions/`.
- Workflow execution now writes a local run record JSON file under `results/workflows/runs/<run_id>/run-record.json`.
- `infra/docker-compose.yml` now bind-mounts `../results` into the API container at `/app/results`, so generated workflow artifacts are visible on the host machine.
- Workflow API responses now convert container-only `/app/results/...` paths into host-visible `results/...` paths.
- `workflow-builder` save/run copy now describes persistence as `PostgreSQL + results/...` instead of the old in-memory/API-store wording.
- Dataset upload now validates supported workflow input types before saving:
  - FASTQ
  - FASTA
  - GTF/GFF
  - BAM
  - VCF.GZ
  - TSV
- Unsupported raw files now fail with `400 Unsupported raw file type for workflow input` instead of being accepted as placeholder data.
- Workflow execution now resolves validated project datasets into unbound node inputs by canonical type.
- Workflow execution now fails fast with `400` when required validated project inputs are missing.
- `workflow-builder` now shows `Project Inputs` with:
  - filename
  - input role
  - canonical type
  - validation status
  - validation detail
- `workflow_runs` and `workflow_definitions` in PostgreSQL now persist `project_id` as a real column, not only inside JSON payloads.
- `api/requirements.txt` now includes `httpx` and `pytest` so the documented FastAPI test path can actually run in a clean environment.
- Verified persistence run:
  - workflow id: `demo-project__verify-host-persist`
  - run id: `demo-project__verify-host-persist-dist-1773323067792`
  - host-visible files:
    - `results/workflows/definitions/demo-project__verify-host-persist.json`
    - `results/workflows/runs/demo-project__verify-host-persist-dist-1773323067792/run-record.json`
    - `results/workflows/runs/demo-project__verify-host-persist-dist-1773323067792/task-000/Snakefile`
    - `results/workflows/runs/demo-project__verify-host-persist-dist-1773323067792/task-000/config.json`
- Verified PostgreSQL rows exist in:
  - `workflow_definitions`
  - `workflow_runs`
- Verified input-validation run:
  - project id: `verify-input-binding`
  - workflow id: `verify-input-binding__wgs-min`
  - invalid upload rejected with `400`
  - run blocked with missing `reference.genome.fasta`
  - run succeeded after uploading `demo_reference.fasta`
- Verified `project_id` column persistence run:
  - project id: `verify-project-column`
  - workflow id: `verify-project-column__wgs-min`
  - `workflow_runs.project_id` stored as `verify-project-column`
- IGV/locus page was already redesigned into a results workspace and connected to workflow/causal run data.
- Confirmed-unused files removed:
  - `api/api/rbac 2.py`
  - `login.json`
  - `infra/login.json`
  - `tmp_workflow_builder.html`
  - `tmp_workflow_chunk.js`
  - `web/public/scores.sample.bed`
  - `web/tsconfig.tsbuildinfo`

## Current UI Intent

- Only expose controls that are useful right now.
- Avoid showing speculative or disconnected features.
- Keep wording short and obvious.
- Keep workflow actions inside the builder instead of sending the user back to the console.

## Main Files Recently Changed

- `AGENTS.md`
- `docs/CLEANUP_CANDIDATES.md`
- `web/app/page.tsx`
- `web/app/login/page.tsx`
- `web/app/workflow-builder/page.tsx`
- `web/app/locus/[region]/page.tsx`
- `web/app/globals.css`
- `web/app/layout.tsx`
- `web/components/app-shell.tsx`
- `web/components/ui/button.tsx`
- `web/components/ui/card.tsx`
- `web/components/ui/input.tsx`
- `web/components/ui/label.tsx`
- `web/components/ui/badge.tsx`
- `web/components/ui/select.tsx`
- `web/components/ui/separator.tsx`
- `web/components/ui/textarea.tsx`
- `web/components.json`
- `web/tailwind.config.ts`
- `web/postcss.config.js`
- `web/lib/api-base.ts`
- `web/lib/session.ts`
- `web/lib/utils.ts`
- `web/Dockerfile`
- `api/api/workflows.py`
- `api/api/datasets.py`
- `api/api/canonical_types.py`
- `api/tests/test_datasets.py`
- `api/tests/test_workflows.py`
- `api/requirements.txt`
- `infra/docker-compose.yml`
- `infra/raw_demo.fastq`
- `infra/sample_data/demo_reads_R1.fastq`
- `infra/sample_data/demo_reads_R2.fastq`
- `infra/sample_data/README.md`

## Current Test Entry

### Login

- URL: `http://localhost:13000/login`
- Admin account:
  - email: `wenhow14@gmail.com`
  - password: `2929q36s82ds`

### Workflow Builder

- URL: `http://localhost:13000/workflow-builder`
- Current expected test flow:
  1. Set `Project ID` to `demo-project`
  2. Use `Upload Raw`
  3. Select one or more files. You can add them in multiple rounds:
     - `infra/raw_demo.fastq`
     - `infra/sample_data/demo_reads_R1.fastq`
     - `infra/sample_data/demo_reads_R2.fastq`
  4. Click `Load Template`
  5. Click `Save`
  6. Click `Run`
  7. Check `Result`

### Console

- URL: `http://localhost:13000/`
- Use when testing direct upload or simple method runs outside the builder.

### Locus

- Example URL: `http://localhost:13000/locus/chr1%3A1-1000`

## Known Issues

- `web` build still shows one existing eslint warning in `web/app/locus/[region]/page.tsx` about a `useMemo` dependency.
- The compose file still warns that the `version` field is obsolete.
- The WGS builder is intentionally constrained and not yet a generic pipeline editor.
- API container and `/health` are currently healthy; the recent API-access fix was on the frontend URL-resolution side rather than the FastAPI service itself.
- `workflow-builder` and `console` still need a full end-to-end browser retest after the newest session-gating change.
- `workflow-builder` page shell and layout have been rebuilt and verified by `lint`, `next build`, Docker rebuild, and HTTP `200`, but still need manual browser UX review for drag behavior.
- Direct runtime persistence verification passed through the Docker API service.
- Repo tests now pass in a clean temporary venv with `api/requirements.txt`, but the local checked-in `api/.venv` is still not the verified path.
- `workflow-builder` frontend passed `lint` and `next build` after the persistence-copy update.
- Some files appear redundant or transitional; see `docs/CLEANUP_CANDIDATES.md` before deleting anything ambiguous.

## Next Recommended Step

1. Connect uploaded raw/reference files more explicitly to workflow execution inputs in the Snakemake materialization layer.
2. Surface saved workflow definition paths in the builder UI, not just run-record paths.
3. Fix the local Python test environment so `api/.venv/bin/pytest` can run without dependency/SSL friction.
