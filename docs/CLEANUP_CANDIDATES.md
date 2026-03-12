# Cleanup Candidates

This file tracks files that look redundant, transitional, generated, or demo-specific.

The goal is safe cleanup. Only clearly generated files should be deleted without further review.

## Already Cleaned

- `web/tsconfig.tsbuildinfo`
  - generated TypeScript build cache
  - safe to delete
- `api/api/rbac 2.py`
  - duplicate file of `api/api/rbac.py`
  - deleted after confirming there were no references
- `login.json`
- `infra/login.json`
  - stale auth-output artifacts
  - deleted
- `tmp_workflow_builder.html`
- `tmp_workflow_chunk.js`
  - temporary debug artifacts
  - deleted
- `web/public/scores.sample.bed`
  - unused public fallback asset
  - deleted

## Review Before Deleting

### Demo output artifacts

- `infra/demo-output/report-run-1.html`
- `infra/demo-output/report-run-1.md`
  - likely generated demo reports
  - currently referenced by `README.md`
  - keep if they are part of demo collateral
  - remove if they are no longer needed as sample output

### Sample/public assets

- `infra/assets/genes.sample.bed`
- `web/public/genes.sample.bed`
  - `infra/assets/genes.sample.bed` is currently used by `infra/docker-compose.yml` for MinIO seeding
  - `web/public/genes.sample.bed` is currently used by the locus page fallback path
  - some duplication is intentional for MinIO seeding vs direct web serving
  - verify actual runtime path usage before cleanup

## Suggested Next Cleanup Pass

1. Decide whether demo output files should stay as demo fixtures or move to a dedicated demo-fixtures folder.
2. Consolidate duplicated BED assets only after confirming which runtime path is authoritative.
