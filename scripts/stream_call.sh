#!/usr/bin/env bash
set -euo pipefail

# Stream a CRAM from a URL, extract a region, sort/index BAM, and call variants.
# Defaults target ESR1/ESR2 when region/label omitted.
#
# Examples (run inside container; repo mounted at /work):
#   bash scripts/stream_call.sh chr6:146809107-157129619 NA12878.ESR1
#   bash scripts/stream_call.sh chr14:59226707-69338613 NA12878.ESR2
#
# Environment overrides:
#   CRAM_URL   - source CRAM (default: NA12878 30X at EBI)
#   REF_FASTA  - reference FASTA path inside container (default: ref/GRCh38.chr6_14.fa)

REGION=${1:-chr6:146809107-157129619}
LABEL=${2:-NA12878.ESR1}

CRAM_URL=${CRAM_URL:-https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram}
REF_FASTA=${REF_FASTA:-ref/GRCh38.chr6_14.fa}

out_bam="bam/${LABEL}.bam"
out_vcf="results/vcf/${LABEL}.vcf.gz"

mkdir -p bam results/vcf

echo "[stream_call] Region: ${REGION}"
echo "[stream_call] Label:  ${LABEL}"
echo "[stream_call] CRAM:   ${CRAM_URL}"
echo "[stream_call] REF:    ${REF_FASTA}"

# Stream region to BAM
samtools view -T "${REF_FASTA}" -b "${CRAM_URL}" "${REGION}" | \
  samtools sort -o "${out_bam}"
samtools index "${out_bam}"

# Variant calling
bcftools mpileup -f "${REF_FASTA}" -r "${REGION}" "${out_bam}" | \
  bcftools call -mv -Oz -o "${out_vcf}"
bcftools index -t "${out_vcf}"

echo "[stream_call] Done."
echo "BAM : ${out_bam}"
echo "VCF : ${out_vcf}"
