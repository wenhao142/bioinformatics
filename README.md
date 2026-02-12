# AD sex-bias fine-mapping & menopause transcriptome scaffold

> 整理：人類 ESR 相關 raw/result 已移至 `projects/human_esr/`；小鼠微膠流程維持在 `projects/mouse_microglia/`.

## Quick map of what’s ready
- **Streamed WGS (no local alignment):** NA12878 30X CRAM is streamed to ESR1/ESR2 windows → `bam/NA12878.ESR1/2.bam` + VCFs in `results/vcf/`.
- **PLINK conversion (single-sample):** `tools/plink2/plink2` (ARM) converts those VCFs to bed/bim/fam in `results/plink/NA12878.ESR1/2.*`.
- **Expression pipeline (menopause ↔ female AD blood):** new Snakemake file `workflow/menopause_ad.smk` + R scripts to fetch GEO, run limma DE, build signatures, and score AD samples.
- **Existing WGS pipeline scaffold:** `workflow/Snakefile` keeps the FASTQ→GATK→PLINK→GCTA/FINEMAP placeholders if you later add aligned BAMs/gVCFs.
- **References:** Only chr6/chr14 reduced FASTA kept: `ref/GRCh38.chr6_14.fa` (+ index/bwa index). Full genome removed to save space.

## Local paths you’ll interact with
- IGV local session (no remote CRAM): `igv_esr_local.xml`
- BAM/VCF for ESR1/2: `bam/NA12878.ESR1.bam`, `bam/NA12878.ESR2.bam`, `results/vcf/NA12878.ESR1/2.vcf.gz`
- PLINK outputs: `results/plink/NA12878.ESR1.*`, `results/plink/NA12878.ESR2.*`
- Transcriptome scripts: `scripts/fetch_geo.R`, `scripts/de_limma.R`, `scripts/build_signature.R`, `scripts/score_signature.R`
- Transcriptome workflow: `workflow/menopause_ad.smk`, config at `config/menopause_ad.yaml`

## How to view ESR1/ESR2 in IGV
1) Genome = `hg38`.  
2) File → Open Session → `igv_esr_local.xml`.  
3) Go to:
   - ESR1 `chr6:146,809,107-157,129,619`
   - ESR2 `chr14:59,226,707-69,338,613`
Zoom in (<50 kb) to see reads; the BAMs only contain these windows.

## Transcriptome (menopause ↔ female AD) pipeline
Prereqs: R with packages `data.table`, `GEOquery`, `limma` (install via `BiocManager::install(c("GEOquery","limma")); install.packages("data.table")`).

Config (`config/menopause_ad.yaml`):
- `gse_ids`: GSE3492, E-GEOD-2208, GSE97760 (female AD blood).
- `designs`: columns to use as group labels; edit if GEO field names differ.
- `menopause_sets`: which studies build the menopause signature (defaults: GSE3492, E-GEOD-2208).

Run:
```bash
snakemake -s workflow/menopause_ad.smk --cores 4
```
Outputs:
- GEO downloads under `data/geo/<GSE>/` (`expression.tsv.gz`, `pheno.tsv`).
- DE tables + volcano: `results/de/<GSE>.de.tsv/.pdf`
- Menopause signature gene lists: `results/signatures/<GSE>.up.txt` / `.down.txt`
- Signature scores applied to AD dataset: `results/ssgsea/ad_signature_scores.tsv`

Notes:
- Group columns in GEO pheno often look like `characteristics_ch1`… adjust in config if mismatch.
- DE uses limma, contrast = case − control.
- Signature score is rank-based (up mean rank – down mean rank).

## WGS/finemap scaffold (unchanged but pared to chr6/14)
- Main Snakefile: `workflow/Snakefile`
- Config: `config/config.yaml`; sample manifest: `samples.tsv`
- Steps: prefetch SRA → FASTQ → FastQC → BWA-MEM2 → markdup → GATK HaplotypeCaller (gVCF) → PLINK region → GCTA/FINEMAP placeholders.
- To rerun alignment you’d need full FASTQs; currently removed to save space. Use `scripts/stream_call.sh` + CRAM URL to avoid local alignment.

## Container option
- Dockerfile at `docker/Dockerfile` builds an ARM64 image with samtools/bcftools/plink/plink2 for reproducible runs (note: host Docker must support linux/arm64).
```bash
docker build -t ad-finemap -f docker/Dockerfile .
docker run --rm -v $PWD:/work -w /work ad-finemap bash scripts/stream_call.sh chr6:146809107-157129619 NA12878.ESR1
```

## Space & cleanup
- Full genome removed; repo footprint ~2.3 GB. Largest items now: `ref/GRCh38.chr6_14.fa*` (~1.8 GB total) and ESR1/ESR2 BAMs (~120 MB each).
- If tight on space, you can delete `results/vcf` or `bam/NA12878.*` and regenerate via `scripts/stream_call.sh`.

## Run log (2026-02-02)
- Built 1000G EUR LD panels for ESR1/ESR2 (503 samples) → `projects/human_esr/ref/1kg_eur_ld/1kg.EUR.ESR1/2.*`.
- Filtered AD GWAS GCST90027158 (GRCh38) to ESR windows → `results/gcta/GCST90027158.EUR.ESR1/2.cojo_input.txt`.
- COJO with p=1e-4:
  - ESR1: 11 signals → `1kgEUR.ESR1.p1e-4.jma.cojo`, `1kgEUR.ESR1.p1e-4.cma.cojo`, plot `plots/1kgEUR_ESR1_p1e4.manhattan.png`.
  - ESR2: 6 signals → `1kgEUR.ESR2.p1e-4.jma.cojo`, `1kgEUR.ESR2.p1e-4.cma.cojo`, plot `plots/1kgEUR_ESR2_p1e4.manhattan.png`.
- COJO with p=5e-6: no signals for ESR1/ESR2 (logs + badsnp/freq files only).
- Checked p=1e-4 lead SNPs against single-sample NA12878 VCFs (`results/vcf/NA12878.ESR1/2.vcf.gz`): none present (likely ref/ref; NA12878 is healthy control).
## Quick commands reference
```bash
# Stream CRAM region → BAM+VCF (already done for ESR1/2)
./scripts/stream_call.sh chr6:146809107-157129619 NA12878.ESR1
./scripts/stream_call.sh chr14:59226707-69338613 NA12878.ESR2

# PLINK conversion (single-sample VCF)
tools/plink2/plink2 --vcf results/vcf/NA12878.ESR1.vcf.gz \
  --vcf-allow-no-nonvar --max-alleles 2 --set-missing-var-ids @:# \
  --double-id --make-bed --out results/plink/NA12878.ESR1

# Transcriptome pipeline
snakemake -s workflow/menopause_ad.smk --cores 4
```
# bioinformatics
