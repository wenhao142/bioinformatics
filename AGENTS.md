# AGENTS.md

# Snakemake Visual Workflow Console

This document defines the development specification for an intranet-first bioinformatics workflow console centered on visualizing and running Snakemake pipelines.

The AI agent must follow this document strictly.

---

# 1. System Vision

Build a modular, extensible, reproducible bioinformatics workflow console where:

- The primary user experience is a web-based drag-and-drop interface for Snakemake workflows
- The first product goal is not a general workflow SaaS, but a focused internal console for biotech teams
- The first workflow is a single end-to-end pipeline rendered as blocks with fixed ports
- Execution is driven by Snakemake CLI, with generated `Snakefile` and config artifacts
- All runs are reproducible
- New tools can be added through structured tool definitions without rewriting the builder core
- IGV is used as a result exploration surface, not as the orchestration engine

---

# 2. Non-Negotiable Architectural Rules

AI MUST NOT:

- Treat IGV as part of the pipeline scheduler
- Expose arbitrary user-defined ports in the MVP builder
- Let users wire incompatible inputs/outputs manually in ways the platform cannot validate
- Modify canonical data types without version bump
- Bypass workflow validation or Snakemake dry-run before execution
- Execute tools outside the controlled runner path
- Break backward compatibility

AI MUST:

- Use schema-driven validation for workflow JSON and tool definitions
- Keep the first builder constrained to fixed tool cards and fixed I/O ports
- Generate deterministic `Snakefile` and config outputs from the visual DAG
- Record full provenance metadata
- Maintain deterministic execution
- Version all breaking changes

---

# 3. Core Components

## 3.1 Method Registry

Responsible for:

- Storing tool definitions used by the visual builder
- Validating tool definition schema
- Version control of tool specs
- Querying available blocks for the UI palette

## 3.2 Plugin Manifest Specification

Each tool definition must define:

- id
- version (semantic versioning)
- execution backend (`snakemake-rule`, later `container-wrapper` if needed)
- inputs (typed)
- outputs (typed)
- params_schema
- snakemake rule template or command template
- UI metadata for labels and grouping

Definition must pass platform JSON schema validation.

---

# 4. Canonical Data Types

Platform-defined types:

- reads.fastq.gz
- align.bam (+ index)
- variants.vcf.gz (+ index)
- expression.counts.tsv
- expression.diff_table.tsv
- report.html

New types require:

- Version increment
- Migration documentation

---

# 5. Workflow Engine

Workflows are Directed Acyclic Graphs (DAG).

Engine responsibilities:

- Validate DAG acyclic property
- Validate input/output compatibility
- Resolve execution order (topological sort)
- Support branching
- Generate `Snakefile` and config from validated graph
- Run `snakemake -n` before real execution
- Support workflow export/import JSON

---

# 6. Execution Layer

Execution runner must:

1. Materialize workflow JSON into `Snakefile`, config, and run directory
2. Resolve input datasets and references
3. Inject parameters into Snakemake config
4. Execute `snakemake` in controlled mode
5. Capture logs and exit code
6. Collect outputs
7. Compute output hashes

Execution must be isolated and reproducible.

---

# 7. Provenance & Reproducibility

Each run must record:

- Input dataset SHA256
- Tool IDs + versions used in the workflow
- Snakemake version
- Generated workflow artifact hashes
- Parameter JSON
- Output hashes
- Execution timestamp

Runs must be re-executable.

## 7.1 Persistence Rules

Workflow definitions and run results must not exist only in transient API memory.

Required persistence targets:

- Local filesystem on the deployed machine
- PostgreSQL for queryable metadata and run indexing

Minimum persistence behavior:

- `Save workflow` must persist the workflow definition to PostgreSQL
- Workflow export/import JSON must also be available on local disk
- `Run workflow` must persist run metadata to PostgreSQL
- Generated runtime artifacts must be stored on local disk
- Result summaries shown in UI must be reconstructable from PostgreSQL + local artifacts

Transient in-memory storage may be used only as a short-lived cache, never as the system of record.

## 7.2 Local Artifact Requirements

Completed or failed runs must write a local run directory containing at least:

- generated `Snakefile`
- generated config
- stdout/stderr or tool logs
- output file paths
- reproducibility metadata

The UI must make it clear where the run was saved locally.

## 7.3 Database Requirements

PostgreSQL must store at least:

- workflow_id
- project_id
- created_by
- workflow JSON
- run_id
- run status
- local artifact paths
- timestamps
- summary metrics
- provenance metadata

The database layer is the canonical query source for saved workflows and prior runs.

## 7.4 Cloud Execution Direction

The system should support a future flow where saved workflows and selected inputs can be dispatched from the local console to a cloud execution environment.

Design rules:

- Local save and PostgreSQL persistence happen first
- Cloud execution is an additional execution target, not a replacement for local persistence
- Cloud submission must reference persisted workflow definitions and persisted input records
- Cloud runs must return status and result references back to the local system
- Local UI must be able to distinguish local runs vs cloud runs

---

# 8. Adapter Pattern

If new tools use incompatible formats:

- Implement an adapter tool
- Convert to canonical type
- Do not alter core data model

---

# 9. Product Scope

System must support:

- WGS pipeline visualization and execution first
- RNA-seq as the next pipeline after WGS stabilizes
- Result exploration through an IGV-based analysis workspace
- Future statistical overlays inside the locus viewer

System does not need, in the MVP:

- A fully generic plugin marketplace
- A public SaaS-style multi-tenant product shell
- Nextflow parity before Snakemake execution is production-ready

---

# 9.1 Frontend Design Source

For web UI redesign or visual refinement work, AI must treat the local `impeccable-style-universal/` folder as the primary design-guidance source.

Scope:

- Applies to `web/` UI structure, copy, layout, spacing, hierarchy, responsiveness, and interaction polish
- Does not override backend architecture, bioinformatics workflow rules, provenance rules, or persistence requirements

Required behavior:

- Read the relevant guidance in `impeccable-style-universal/` before major UI redesign work
- Prefer quieter, clearer, more restrained interfaces over loud or generic dashboard styling
- Keep copy short and task-oriented
- Adapt layouts for both desktop and mobile instead of designing desktop-only shells
- Use the guidance as a refinement system, not as permission to add decorative UI that does not support the workflow

Current preferred guidance themes:

- quieter
- adapt
- polish

If the UI direction from `impeccable-style-universal/` conflicts with explicit product workflow rules in this file, product workflow rules win.

---

# 10. Autonomous Development Protocol

AI must follow this loop:

1. Read AGENTS.md
2. Locate first unchecked task
3. Break task into subtasks
4. Implement
5. Validate with tests
6. If tests pass ??check box
7. Log changes
8. Repeat

AI must not skip validation.

Strict reporting rule:

- Do not tell the user a change is complete until the relevant verification has passed locally.
- If verification is still running, incomplete, timed out, or failed, say that directly instead of claiming success.
- Report in this order only:
  1. verification result
  2. what changed
- Never reverse that order.

Minimum verification expectations:

- Frontend change:
  - type check
  - lint
  - build or runtime check when rendered behavior changed
- Backend change:
  - syntax/lint or targeted test
  - targeted API check when behavior changed
- Workflow or execution change:
  - targeted run, dry-run, or artifact verification

If a change is only partial, AI must say it is partial and continue until the verified state is clear.

Additional required rule:

- After any meaningful code, UI, workflow, demo-data, or test-flow change, AI must update `docs/WORKING_CONTEXT.md`.
- `docs/WORKING_CONTEXT.md` must describe current product state, changed files, how to test, known issues, and the next recommended step.
- `AGENTS.md` stores durable rules. `docs/WORKING_CONTEXT.md` stores current project state.

---

# 11. Development Modes

## Safe Mode
- Small changes only
- No schema modification

## Refactor Mode
- Internal improvements
- No external behavior change

## Expansion Mode
- Add plugins
- Add workflow templates
- Extend ecosystem

Mode must be explicitly declared before major updates.

---

# 12. MVP Visual Snakemake Builder

Goal: Implement the first end-to-end visual Snakemake builder that supports one real workflow with 5-7 common tools and can run it from the UI.

Constraints:
- The MVP UI is a fixed-node builder, not a generic plugin canvas.
- Each block has fixed input/output ports defined by the platform.
- The initial pipeline target is WGS variant calling.
- Keep the first version small; extend later by adding more tool specs and templates.

## 12.1 Builder Model

The builder is a thin web console over Snakemake, not a replacement workflow language.

MVP rules:
- Users drag predefined blocks from the palette onto a canvas.
- Users may edit parameters, but not invent new ports or arbitrary execution semantics.
- Connections are only allowed when canonical input and output types match.
- The saved workflow JSON is transformed into `Snakefile` + config by the backend.
- The first run action must support `dry-run` and `real run`.

This is described as "fixed I/O ports, no generic plugin UI" because:
- it reduces invalid states,
- makes Snakemake generation deterministic,
- keeps the first production feature understandable by non-programmers,
- and still remains extensible by adding new tool definitions later.

## 12.1.a Product UX Rules

UI rules for the current product direction:

- Show only fields and controls that are actively used in the current workflow.
- Hide unfinished, unconnected, or speculative features from the main UI.
- Keep copy short, direct, and understandable by non-engineering users.
- The primary builder flow is: `upload raw -> build -> save -> run -> result`.
- Do not mix login/session/account-switching controls into the main workflow workspace.
- In the builder, prefer raw data upload over exposing intermediate artifacts that would confuse the workflow mental model.
- Save actions must explain clearly where workflow definitions are stored.
- Run result screens must show human-readable summaries first, raw JSON only on demand.
- Run result screens must show where outputs were saved locally.

## 12.1.b Input Validation Rules

Workflow input handling must be strict.

Required rules:

- Upload success does not mean workflow compatibility.
- The system must not treat any arbitrary uploaded file as a valid workflow input just because it was stored successfully.
- Raw upload must validate:
  - allowed file category for the selected workflow mode
  - basic content or structure when practical
- The UI must distinguish clearly between:
  - file uploaded
  - file validated for workflow use
  - file bound to a workflow input
- Workflow execution must fail before execution if required inputs are:
  - missing
  - unbound
  - incompatible
  - structurally invalid for the selected workflow
- The system must not return a successful placeholder workflow result for invalid or unbound inputs.

## 12.2 Canonical Types used in MVP

**Core types (already defined):**
- `reads.fastq.gz`
- `align.bam` (+ index)
- `expression.counts.tsv`
- `expression.diff_table.tsv`
- `report.html`

**MVP additions (allowed as aliases without changing canonical registry):**
- `reads.fastq_pair.gz` (pair of `reads.fastq.gz`)
- `qc.fastqc.zip` (FastQC output bundle)
- `qc.multiqc.html` (MultiQC report)
- `reference.genome.fasta` (FASTA + index bundle)
- `annotation.gtf` (gene annotation)

If these are promoted to canonical types later, do a version bump and document migration.

---

## 12.3 Minimal WGS DAG (v0) - primary MVP

Pipeline:
1. FastQC
2. Cutadapt
3. BWA-MEM
4. samtools sort+index
5. GATK HaplotypeCaller
6. GATK GenotypeGVCFs or a simplified variant call step
7. Optional report/QC aggregation

Tool blocks for MVP:
- `fastqc`
- `cutadapt`
- `bwa_mem`
- `samtools_sort_index`
- `gatk_haplotypecaller`
- `gatk_genotypegvcfs`
- `variant_report` or `multiqc` (optional)

Fixed-port examples:
- `fastqc`
  - inputs: `reads`
  - outputs: `fastqc_zip`, `fastqc_html`
- `bwa_mem`
  - inputs: `reads`, `reference`
  - outputs: `alignment`
- `samtools_sort_index`
  - inputs: `alignment`
  - outputs: `alignment_sorted`, `alignment_index`
- `gatk_haplotypecaller`
  - inputs: `alignment_sorted`, `reference`
  - outputs: `gvcf`
- `gatk_genotypegvcfs`
  - inputs: `gvcf`, `reference`
  - outputs: `variants`

Acceptance (WGS v0):
- DAG validates and can only connect compatible ports.
- Backend generates deterministic `Snakefile` and config files.
- `snakemake -n` works from the generated workflow.
- Real execution can be launched from API/UI.
- Each node stores params, versions, input hashes, output hashes, logs, and exit status.
- Final artifacts visible in UI:
  - FastQC HTML
  - Sorted/indexed BAM
  - VCF / gVCF outputs
  - Basic run report

## 12.4 Expansion path after WGS

After the WGS path is stable:
- add RNA-seq as a second fixed workflow template,
- add more tools by extending the tool definition catalog,
- keep the builder constrained until at least two real pipelines execute end-to-end reliably.

## 12.5 IGV Result Workspace

IGV is part of the result viewer, not the scheduler.

Result workspace requirements:
- Hide most raw IGV chrome and wrap `igv.js` in a custom UI shell.
- Place locus navigation, run selector, sample selector, and track toggles in the app shell.
- Add a right-side evidence panel for variant, gene, and provenance details.
- Add synchronized statistical panels below or beside the genome viewer.
- Support future overlays such as p-values, posterior probabilities, effect sizes, gene scores, and expression/protein summaries.

Acceptance (IGV redesign v0):
- Viewer looks like an integrated product surface, not a default IGV embed.
- Track selection and feature click update custom evidence panels.
- Statistical summaries can sync to the currently viewed locus.
- IGV remains optional for pipeline execution; failure in viewer code must not block workflow runs.

---

# 13. AI Task Checklist

## Core Infrastructure

- [x] Method Registry implemented
- [x] Manifest JSON schema validator
- [x] Canonical data type system
- [x] Plugin runner completed
- [x] Provenance tracking
- [x] Run comparison engine

---

## Workflow Engine

- [x] DAG validation
- [x] I/O compatibility checker
- [x] Topological execution resolver
- [x] Branching support
- [x] Parameter sweep mode
- [x] Workflow import/export

---

## Plugin Ecosystem

- [x] Adapter plugin framework
- [x] Semantic versioning enforcement
- [x] Deprecated plugin policy
- [x] Plugin isolation validation
- [x] Registry search API

---

## Bioinformatics Templates

- [x] WGS template
- [x] RNA-seq template
- [x] GWAS template
- [x] eQTL template
- [x] Multi-omics template

---

## Advanced Integration

- [ ] Snakemake execution from generated workflow files
- [x] Nextflow wrapper
- [x] Distributed execution support
- [x] Cloud storage integration

---

## UI Layer

- [ ] Fixed-node Snakemake workflow builder
- [ ] Drag-and-drop composition with fixed ports only
- [ ] Parameter configuration panel for tool-specific params
- [ ] Run monitoring dashboard with Snakemake job state/logs
- [x] Reproducibility report export
- [x] Web raw bio file upload (FASTA/GTF/FASTQ/BAM via dataset API)
- [ ] IGV redesign with evidence dock and statistics panels
- [ ] Statistical overlays synchronized with viewed locus

---

# 13. Completion Policy

A checkbox may only be marked complete if:

- Implementation exists
- Unit tests pass
- Integration tests pass
- No architectural rules violated
- Every execution cycle after reading AGENTS.md must include verification testing before marking completion.

## 13.1 Latest Verification Log (2026-02-26)

- API: `ruff check .` + `pytest -q` => passed (`64 passed`)
- Workers: `pytest -q` => passed (`1 passed`)
- Workflow utils: `pytest -q workflow/tests/test_platform_tools.py` => passed (`4 passed`)
- Web: `npm run lint` + `npm run build` => passed
- E2E smoke: `infra/demo_smoke.ps1` => passed (uploads raw FASTQ/BAM + sample FASTA/GTF/paired FASTQ/sample sheet, runs evidence/causal/report)
- UI upload extension: analysis console now supports browser upload for FASTA/GTF/FASTQ/BAM (and related raw files) through `/datasets/upload?project_id=...`

## 13.2 Execution Log Update (2026-02-26)

- Web verification re-run: `cd web && npm run lint && npm run build` => passed
- API verification re-run: `cd api && pytest -q tests/test_datasets.py::test_upload_generated_bio_sample_files` => passed (`1 passed`)

---

End of AGENTS.md
