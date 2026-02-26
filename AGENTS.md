# AGENTS.md

# Autonomous Bioinformatics Workflow Engine

This document defines the full autonomous development specification for an AI-driven, plugin-based bioinformatics workflow platform.

The AI agent must follow this document strictly.

---

# 1. System Vision

Build a modular, extensible, reproducible bioinformatics workflow engine where:

- Each bioinformatics tool is a plugin (building block)
- Workflows are DAG-based compositions of plugins
- Execution is containerized (Docker-first)
- All runs are reproducible
- New tools can be added without modifying core logic
- AI can autonomously extend the system

---

# 2. Non-Negotiable Architectural Rules

AI MUST NOT:

- Hardcode tool names in core logic
- Modify canonical data types without version bump
- Bypass manifest validation
- Execute tools outside container runner
- Break backward compatibility

AI MUST:

- Use schema-driven validation
- Record full provenance metadata
- Maintain deterministic execution
- Version all breaking changes

---

# 3. Core Components

## 3.1 Method Registry

Responsible for:

- Registering plugins
- Validating plugin manifest JSON schema
- Version control of plugins
- Querying available tools

## 3.2 Plugin Manifest Specification

Each plugin must define:

- id
- version (semantic versioning)
- container_image (pinned tag)
- inputs (typed)
- outputs (typed)
- params_schema
- resource requirements
- execution command template

Manifest must pass platform JSON schema validation.

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
- Support parameter sweep mode
- Support workflow export/import JSON

---

# 6. Execution Layer

Plugin Runner must:

1. Pull container image
2. Mount inputs
3. Inject parameters
4. Execute command
5. Capture logs
6. Collect outputs
7. Compute output hashes

Execution must be isolated and reproducible.

---

# 7. Provenance & Reproducibility

Each run must record:

- Input dataset SHA256
- Plugin ID + version
- Container image digest
- Parameter JSON
- Output hashes
- Execution timestamp

Runs must be re-executable.

---

# 8. Adapter Plugin Pattern

If new tools use incompatible formats:

- Implement adapter plugin
- Convert to canonical type
- Do not alter core data model

---

# 9. Advanced Workflow Support

System must support:

- WGS pipelines
- RNA-seq pipelines
- GWAS pipelines
- eQTL pipelines
- Multi-omics integration
- Snakemake wrapper plugin
- Nextflow wrapper plugin

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

# 12. MVP Workflow Composer (Minimal DAG)

Goal: Implement the **first end-to-end workflow composer** that supports **only 5?? common tools** and can run a complete **RNA-seq** (preferred) or **WGS** pipeline as a DAG.

Constraints:
- Tools are executed ONLY via plugin runner (containerized).
- Workflows are validated by canonical IO types.
- Keep the first version small; extend later via registry.

## 12.1 Canonical Types used in MVP

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

## 12.2 Minimal RNA-seq DAG (v0) ??6 tools

**Pipeline:**
1) FastQC ??2) Cutadapt ??3) STAR (or HISAT2) ??4) samtools sort+index ??5) featureCounts ??6) DESeq2 ??(optional) MultiQC

### Tools (Plugins) ??Detailed Specs

#### P1. `fastqc`
- Purpose: raw read QC
- Container: `biocontainers/fastqc` (pinned tag)
- Inputs:
  - `reads`: `reads.fastq.gz` OR `reads.fastq_pair.gz`
- Outputs:
  - `fastqc_zip`: `qc.fastqc.zip`
  - `fastqc_html`: `report.html` (single-sample FastQC HTML)
- Params:
  - `threads` (int, default 2)
- Command template:
  - `fastqc -t {threads} -o {out_dir} {reads...}`

#### P2. `cutadapt`
- Purpose: adapter/quality trimming
- Container: `quay.io/biocontainers/cutadapt` (pinned tag)
- Inputs:
  - `reads`: `reads.fastq.gz` OR `reads.fastq_pair.gz`
- Outputs:
  - `trimmed_reads`: `reads.fastq.gz` OR `reads.fastq_pair.gz`
  - `trim_log`: `report.html` (or `text/plain` stored as artifact; rendered in UI)
- Params (MVP subset):
  - `adapter_fwd` (string, optional)
  - `adapter_rev` (string, optional)
  - `quality_cutoff` (int, default 20)
  - `min_length` (int, default 20)
  - `threads` (int, default 4)
- Command template (paired example):
  - `cutadapt -j {threads} -q {quality_cutoff} -m {min_length} -a {adapter_fwd} -A {adapter_rev} -o {out1} -p {out2} {in1} {in2}`

#### P3. `star_align` (preferred) OR `hisat2_align` (alternative)
- Purpose: RNA-seq alignment
- Container:
  - STAR: `quay.io/biocontainers/star` (pinned)
  - HISAT2: `quay.io/biocontainers/hisat2` (pinned)
- Inputs:
  - `reads`: `reads.fastq.gz` OR `reads.fastq_pair.gz`
  - `reference`: `reference.genome.fasta` (includes required indices)
  - `annotation` (optional for STAR quant modes): `annotation.gtf`
- Outputs:
  - `alignment_unsorted`: `align.bam`
  - `align_log`: `report.html` (or text artifact)
- Params:
  - `threads` (int, default 8)
  - `read_group` (string, optional)
- Command template (STAR BAM unsorted):
  - `STAR --runThreadN {threads} --genomeDir {ref_index} --readFilesIn {reads...} --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix {out_prefix}`

#### P4. `samtools_sort_index`
- Purpose: sort + index BAM
- Container: `quay.io/biocontainers/samtools` (pinned)
- Inputs:
  - `alignment`: `align.bam`
- Outputs:
  - `alignment_sorted`: `align.bam` (+ index)
- Params:
  - `threads` (int, default 4)
  - `memory` (string, default `1G`)
- Command template:
  - `samtools sort -@ {threads} -m {memory} -o {sorted_bam} {bam} && samtools index {sorted_bam}`

#### P5. `featurecounts`
- Purpose: gene-level read counting
- Container: `quay.io/biocontainers/subread` (pinned)
- Inputs:
  - `alignment_sorted`: `align.bam` (+ index)
  - `annotation`: `annotation.gtf`
- Outputs:
  - `counts`: `expression.counts.tsv`
  - `counts_summary`: `report.html` (or text artifact)
- Params:
  - `threads` (int, default 4)
  - `feature_type` (string, default `exon`)
  - `attribute` (string, default `gene_id`)
  - `strand` (int, default 0)
- Command template:
  - `featureCounts -T {threads} -t {feature_type} -g {attribute} -s {strand} -a {gtf} -o {counts_tsv} {bam}`

#### P6. `deseq2_diffexp`
- Purpose: differential expression analysis
- Container: `bioconductor/bioconductor_docker` (pinned) OR custom image with R + DESeq2
- Inputs:
  - `counts`: `expression.counts.tsv`
  - `sample_sheet`: (CSV/TSV; stored as dataset artifact)
- Outputs:
  - `diff_table`: `expression.diff_table.tsv`
  - `qc_plots`: `report.html` (MA plot / PCA / dispersion)
- Params:
  - `design_formula` (string, default `~ condition`)
  - `contrast` (string, required, e.g., `condition,treat,ctrl`)
  - `alpha` (float, default 0.05)
- Command template:
  - `Rscript /app/run_deseq2.R --counts {counts} --samples {samples} --design "{design_formula}" --contrast "{contrast}" --alpha {alpha} --out {out_dir}`

#### P7 (optional). `multiqc`
- Purpose: aggregate QC across steps
- Container: `quay.io/biocontainers/multiqc` (pinned)
- Inputs:
  - `qc_inputs`: list of artifacts from FastQC/Cutadapt/STAR/featureCounts
- Outputs:
  - `multiqc_report`: `qc.multiqc.html`
- Params:
  - `title` (string, default `RNA-seq QC`)
- Command template:
  - `multiqc {in_dir} -o {out_dir} --title "{title}"`

### Minimal RNA-seq Workflow DAG (example)

Nodes:
- `n1_fastqc_raw` (fastqc)
- `n2_cutadapt` (cutadapt)
- `n3_fastqc_trimmed` (fastqc)
- `n4_align` (star_align OR hisat2_align)
- `n5_sort_index` (samtools_sort_index)
- `n6_featurecounts` (featurecounts)
- `n7_deseq2` (deseq2_diffexp)
- `n8_multiqc` (multiqc, optional)

Edges:
- raw reads ??fastqc_raw
- raw reads ??cutadapt ??fastqc_trimmed
- cutadapt ??align ??sort_index ??featurecounts ??deseq2
- (optional) qc artifacts ??multiqc

Acceptance (RNA-seq v0):
- DAG validates (acyclic + IO types compatible)
- Runner executes all nodes in topological order
- Each node stores: params + tool versions + input hashes + output hashes + logs
- Final artifacts visible in UI:
  - FastQC HTML
  - Alignment BAM (download link)
  - Counts TSV
  - DE diff table TSV
  - Report HTML

---

## 12.3 Minimal WGS DAG (v0) ??6 tools (alternative)

**Pipeline:**
1) FastQC ??2) Cutadapt ??3) BWA-MEM ??4) samtools sort+index ??5) bcftools mpileup+call ??6) bcftools filter ??(optional) MultiQC

Tools (plugins): `fastqc`, `cutadapt`, `bwa_mem`, `samtools_sort_index`, `bcftools_call`, `bcftools_filter`.

Acceptance (WGS v0):
- Produces `variants.vcf.gz(+tbi)` and basic QC report.

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

- [x] Snakemake wrapper
- [x] Nextflow wrapper
- [x] Distributed execution support
- [x] Cloud storage integration

---

## UI Layer

- [x] Visual workflow builder
- [x] Drag-and-drop composition
- [x] Parameter configuration panel
- [x] Run monitoring dashboard
- [x] Reproducibility report export
- [x] Web raw bio file upload (FASTA/GTF/FASTQ/BAM via dataset API)

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


