# Local-Only Flow

This flow reflects the local execution direction required by `AGENTS.md`.

Key rules:

- Everything runs locally
- Workflow definitions persist to PostgreSQL
- Workflow/result artifacts persist to local disk
- Every major stage has validation before continuing

## Mermaid

```mermaid
flowchart TD
    A[User opens Workflow Builder] --> A1[Validate login and project context]
    A1 -->|pass| B[Upload raw files]
    A1 -->|fail| AERR[Show auth or project error]

    B --> B1[Validate files selected]
    B1 -->|pass| B2[Upload raw files to local dataset system]
    B1 -->|fail| BERR[Show file selection error]

    B2 --> B3[Validate upload success and dataset records]
    B3 -->|pass| C[Load WGS template]
    B3 -->|fail| BERR2[Show upload error]

    C --> C1[Validate fixed blocks and fixed ports]
    C1 -->|pass| D[User edits workflow on canvas]
    C1 -->|fail| CERR[Show template validation error]

    D --> D1[Validate DAG, node ids, I/O compatibility, required fields]
    D1 -->|pass| E[Save workflow]
    D1 -->|fail| DERR[Show workflow validation issues]

    E --> E1[Persist workflow JSON to PostgreSQL]
    E1 --> E2[Export workflow JSON to local disk]
    E2 --> E3[Validate DB row and local workflow file]
    E3 -->|pass| F[Run workflow]
    E3 -->|fail| EERR[Show save failure]

    F --> F1[Load saved workflow definition]
    F1 --> F2[Resolve uploaded datasets and local references]
    F2 --> F3[Materialize local run directory]
    F3 --> F4[Generate Snakefile and config]
    F4 --> F5[Validate generated artifacts exist]
    F5 -->|pass| G[Snakemake dry-run]
    F5 -->|fail| FERR[Show materialization error]

    G --> G1[Validate dry-run success]
    G1 -->|pass| H[Run Snakemake locally]
    G1 -->|fail| GERR[Show dry-run error]

    H --> H1[Capture logs, exit code, outputs, hashes]
    H1 --> H2[Persist run metadata to PostgreSQL]
    H2 --> H3[Persist artifacts and logs on local disk]
    H3 --> H4[Validate result paths and DB metadata]
    H4 -->|pass| I[Show result summary in UI]
    H4 -->|fail| HERR[Show persistence error]

    I --> I1[Show local run directory, Snakefile, config, outputs]
    I1 --> I2[Allow raw JSON only on demand]
    I2 --> I3[Optional open IGV/locus workspace]
```

## Step List

1. User opens builder
2. Validate login and project context
3. Upload raw files
4. Validate upload
5. Load fixed WGS template
6. Validate builder state
7. Save workflow
8. Write workflow to PostgreSQL
9. Write workflow export JSON to local disk
10. Validate persistence
11. Run workflow
12. Resolve local inputs
13. Materialize run directory
14. Generate `Snakefile` and config
15. Validate generated files
16. Run `snakemake -n`
17. Validate dry-run
18. Run Snakemake locally
19. Capture outputs/logs/hashes
20. Write run metadata to PostgreSQL
21. Write artifacts to local disk
22. Validate stored results
23. Show readable result summary in UI

## Storage Targets

- PostgreSQL:
  - workflow definition
  - workflow metadata
  - run metadata
  - provenance metadata
  - artifact paths

- Local disk:
  - workflow export JSON
  - generated `Snakefile`
  - generated config
  - logs
  - result files
  - reproducibility metadata

## UI Requirements

- `Save` must tell the user the workflow was written to PostgreSQL and exported locally
- `Run` result must show local run directory and output paths
- Validation failures must stop the flow before execution
