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
- Deliver implementations compatible with both Windows and macOS

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
6. If tests pass → check box
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

# 12. AI Task Checklist

## Core Infrastructure

- [x] Method Registry implemented
- [x] Manifest JSON schema validator
- [x] Windows + macOS compatibility
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
- [x] Docker raw-data end-to-end execution
  - Acceptance:
    - User can start Docker stack and submit raw bioinformatics input data (e.g., FASTQ/BAM/VCF) without manual file surgery.
    - Pipeline runs inside containers and produces queryable outputs in platform storage.
    - Verified on both Windows and macOS.

---

## UI Layer

- [x] Visual workflow builder
- [x] Drag-and-drop composition
- [x] Parameter configuration panel
- [x] Run monitoring dashboard
- [x] Reproducibility report export
- [x] Lego-style tool selection flow
  - Acceptance:
    - User can compose/select tools like building blocks in UI and launch a run.
    - UI prevents invalid connections/inputs and shows clear validation errors.
    - Composed workflow is saved and reusable per project.

---

# 13. Completion Policy

A checkbox may only be marked complete if:

- Implementation exists
- Unit tests pass
- Integration tests pass
- Windows and macOS compatibility is verified
- No architectural rules violated

---

End of AGENTS.md
