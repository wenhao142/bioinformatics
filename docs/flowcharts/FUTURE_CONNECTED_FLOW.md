# Future Connected Flow

This flow extends the local-first model from `AGENTS.md`.

Key rules:

- Local persistence happens first
- Cloud execution is optional
- Cloud execution never replaces local storage or PostgreSQL indexing
- Every stage still has validation

## Mermaid

```mermaid
flowchart TD
    A[User opens Workflow Builder] --> A1[Validate login, project, and execution target]
    A1 -->|pass| B[Upload raw files locally]
    A1 -->|fail| AERR[Show auth or target error]

    B --> B1[Validate upload success]
    B1 -->|pass| C[Build workflow locally]
    B1 -->|fail| BERR[Show upload error]

    C --> C1[Validate DAG and fixed ports]
    C1 -->|pass| D[Save workflow locally first]
    C1 -->|fail| CERR[Show workflow error]

    D --> D1[Persist workflow JSON to PostgreSQL]
    D1 --> D2[Persist workflow export JSON to local disk]
    D2 --> D3[Validate local persistence]
    D3 -->|pass| E[Choose execution target]
    D3 -->|fail| DERR[Show save failure]

    E -->|Local| L1[Run local flow]
    E -->|Cloud| F[Prepare cloud submission]

    L1 --> L2[Run local dry-run]
    L2 --> L3[Validate dry-run]
    L3 -->|pass| L4[Run local Snakemake]
    L3 -->|fail| LERR[Show local dry-run error]
    L4 --> L5[Persist local run metadata and artifacts]
    L5 --> L6[Validate local persistence]
    L6 -->|pass| Z[Show result summary in UI]
    L6 -->|fail| LPERR[Show local persistence error]

    F --> F1[Resolve persisted workflow and dataset references]
    F1 --> F2[Create cloud run payload from persisted state]
    F2 --> F3[Validate cloud payload]
    F3 -->|pass| G[Submit to cloud runner]
    F3 -->|fail| FERR[Show submission validation error]

    G --> G1[Persist cloud submission record to PostgreSQL]
    G1 --> G2[Validate submission record]
    G2 -->|pass| H[Cloud runner materializes workflow]
    G2 -->|fail| GERR[Show submission persistence error]

    H --> H1[Cloud dry-run]
    H1 --> H2[Validate cloud dry-run]
    H2 -->|pass| I[Cloud execution]
    H2 -->|fail| HERR[Show cloud dry-run failure]

    I --> I1[Cloud runner captures logs, outputs, hashes]
    I1 --> I2[Cloud returns status and artifact references]
    I2 --> I3[Local system stores cloud run metadata in PostgreSQL]
    I3 --> I4[Optional local sync of selected outputs]
    I4 --> I5[Validate returned status and stored references]
    I5 -->|pass| Z[Show result summary in UI]
    I5 -->|fail| IERR[Show result sync error]

    Z --> Z1[Show execution target: local or cloud]
    Z1 --> Z2[Show local artifact paths or cloud artifact references]
    Z2 --> Z3[Optional open IGV/locus workspace]
```

## Step List

1. Validate local session and project
2. Upload raw data locally
3. Validate raw upload
4. Build workflow locally
5. Validate builder graph
6. Save workflow to PostgreSQL
7. Export workflow JSON to local disk
8. Validate local persistence
9. Choose execution target

### Local target

10. Dry-run locally
11. Validate dry-run
12. Execute locally
13. Store local run metadata
14. Store local artifacts
15. Validate local result persistence

### Cloud target

10. Resolve saved workflow and input references
11. Build cloud submission payload
12. Validate payload
13. Submit job
14. Store cloud submission record in PostgreSQL
15. Validate submission record
16. Cloud dry-run
17. Validate cloud dry-run
18. Cloud execution
19. Receive status, logs, and artifact references
20. Store cloud run metadata locally in PostgreSQL
21. Optionally sync selected outputs back to local disk
22. Validate stored cloud results

## Cloud-Specific Requirements

- Workflow and dataset references must come from persisted local records
- Cloud run must be traceable back to:
  - workflow_id
  - run_id
  - project_id
  - local dataset records
- UI must clearly label:
  - local run
  - cloud run
- Cloud result view must still be queryable from PostgreSQL

## UX Requirements

- User must know before pressing `Run` whether it will execute locally or in the cloud
- Result screen must show:
  - execution target
  - persisted run id
  - where artifacts are stored
  - whether files are local, remote, or synced back
