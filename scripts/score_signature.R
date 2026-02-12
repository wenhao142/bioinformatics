#!/usr/bin/env Rscript

# Simple single-sample signature score (rank-based) for up/down gene sets.
# Usage: Rscript score_signature.R expr.tsv.gz up.txt down.txt out_tsv

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: score_signature.R expr.tsv.gz up.txt down.txt out_tsv")
}
expr_path <- args[[1]]
up_path <- args[[2]]
down_path <- args[[3]]
out_path <- args[[4]]

library(data.table)
expr <- as.matrix(fread(expr_path), rownames = 1)
up <- unique(readLines(up_path))
down <- unique(readLines(down_path))

all_genes <- rownames(expr)
up <- intersect(up, all_genes)
down <- intersect(down, all_genes)

rank_mat <- apply(expr, 2, rank, ties.method = "average")
n <- nrow(expr)

score_one <- function(col_ranks) {
  up_score <- if (length(up)) mean(col_ranks[up]) / n else 0
  down_score <- if (length(down)) mean(col_ranks[down]) / n else 0
  up_score - down_score
}

scores <- apply(rank_mat, 2, score_one)
out <- data.table(sample = names(scores), signature_score = scores)
fwrite(out, file = out_path, sep = "\t")

message("[score_signature] wrote: ", out_path)
