#!/usr/bin/env Rscript

# Download a GEO SeriesMatrix and save expression + phenotype tables.
# Usage: Rscript fetch_geo.R GSE_ID out_dir

suppressPackageStartupMessages({
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Please install GEOquery: BiocManager::install('GEOquery')")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: fetch_geo.R GSE_ID out_dir")
}
gse_id <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("[fetch_geo] downloading ", gse_id)
gset_list <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)
gset <- gset_list[[1]]

expr <- Biobase::exprs(gset)
pheno <- Biobase::pData(gset)
fdata <- Biobase::fData(gset)

# Add sample column for clarity
pheno$sample <- rownames(pheno)

# Write tables
expr_path <- file.path(out_dir, "expression.tsv.gz")
pheno_path <- file.path(out_dir, "pheno.tsv")
feature_path <- file.path(out_dir, "features.tsv")

message("[fetch_geo] writing expression to ", expr_path)
write.table(expr, file = gzfile(expr_path), sep = "\t", quote = FALSE, col.names = NA)

message("[fetch_geo] writing pheno to ", pheno_path)
write.table(pheno, file = pheno_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("[fetch_geo] writing features to ", feature_path)
write.table(fdata, file = feature_path, sep = "\t", quote = FALSE, col.names = NA)

message("[fetch_geo] done")
