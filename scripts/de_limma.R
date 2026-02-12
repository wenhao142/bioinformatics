#!/usr/bin/env Rscript

# Differential expression with limma for two-group comparison.
# Usage:
#   Rscript de_limma.R expr.tsv.gz pheno.tsv group_col case_label ctrl_label out_prefix
#
# expr: genes x samples table (rownames = gene, colnames = samples)
# pheno: must contain column "sample" matching expr colnames, plus group_col

suppressPackageStartupMessages({
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Please install limma: BiocManager::install('limma')")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: de_limma.R expr.tsv.gz pheno.tsv group_col case_label ctrl_label out_prefix")
}
expr_path <- args[[1]]
pheno_path <- args[[2]]
group_col <- args[[3]]
case_label <- args[[4]]
ctrl_label <- args[[5]]
out_prefix <- args[[6]]

library(data.table)
library(limma)

expr <- as.matrix(fread(expr_path), rownames = 1)
pheno <- fread(pheno_path)

if (!"sample" %in% colnames(pheno)) stop("pheno needs a 'sample' column")

# keep samples present in both tables
common <- intersect(colnames(expr), pheno$sample)
expr <- expr[, common, drop = FALSE]
pheno <- pheno[match(common, pheno$sample)]

if (!(group_col %in% colnames(pheno))) stop("group_col not in pheno")
pheno[[group_col]] <- factor(pheno[[group_col]], levels = c(ctrl_label, case_label))

design <- model.matrix(~ 0 + pheno[[group_col]])
colnames(design) <- c("CTRL", "CASE")
contrast <- makeContrasts(CASEvsCTRL = CASE - CTRL, levels = design)

fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, coef = "CASEvsCTRL", number = Inf, sort.by = "P")
tt$gene <- rownames(tt)

out_table <- paste0(out_prefix, ".de.tsv")
fwrite(tt, file = out_table, sep = "\t")

# simple volcano
pdf(paste0(out_prefix, ".volcano.pdf"), width = 6, height = 5)
with(tt, plot(logFC, -log10(adj.P.Val), pch = 16, cex = 0.5,
              xlab = "logFC (case vs ctrl)", ylab = "-log10 FDR"))
abline(h = -log10(0.05), col = "red", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)
dev.off()

message("[de_limma] wrote: ", out_table)
