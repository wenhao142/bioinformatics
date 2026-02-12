#!/usr/bin/env Rscript

# Build up/down gene signatures from a DE table (limma topTable output).
# Usage: Rscript build_signature.R de.tsv logfc_cutoff fdr_cutoff out_prefix

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: build_signature.R de.tsv logfc_cutoff fdr_cutoff out_prefix")
}
de_path <- args[[1]]
logfc_cut <- as.numeric(args[[2]])
fdr_cut <- as.numeric(args[[3]])
out_prefix <- args[[4]]

library(data.table)
de <- fread(de_path)
if (!"gene" %in% colnames(de)) stop("DE table must contain 'gene' column")

up <- de[logFC >= logfc_cut & adj.P.Val <= fdr_cut, gene]
down <- de[logFC <= -logfc_cut & adj.P.Val <= fdr_cut, gene]

writeLines(up, paste0(out_prefix, ".up.txt"))
writeLines(down, paste0(out_prefix, ".down.txt"))

message("[signature] up genes: ", length(up), " | down genes: ", length(down))
