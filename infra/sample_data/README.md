# Demo Bioinformatics Sample Data

This folder contains tiny synthetic files for local smoke/integration tests.

- `demo_reference.fasta`: toy reference genome
- `demo_annotation.gtf`: toy gene annotation
- `demo_reads_R1.fastq`: paired-end read 1
- `demo_reads_R2.fastq`: paired-end read 2
- `demo_samplesheet.tsv`: minimal sample metadata
- `demo_counts.tsv`: minimal count matrix
- `demo_diff.tsv`: minimal differential expression table

Related raw upload fixture outside this folder:

- `../raw_demo.fastq`: single-end FASTQ for the "upload raw data first" flow

Suggested from-scratch test order:

1. Upload `infra/raw_demo.fastq`
2. Upload `infra/sample_data/demo_reference.fasta`
3. Upload `infra/sample_data/demo_annotation.gtf`
4. Upload `infra/sample_data/demo_reads_R1.fastq`
5. Upload `infra/sample_data/demo_reads_R2.fastq`
6. Upload `infra/sample_data/demo_samplesheet.tsv`
7. Upload `infra/example.vcf` and `infra/expr.tsv` if you also want evidence/causal/report testing

All files are intentionally small and non-biological, for test flow validation only.
