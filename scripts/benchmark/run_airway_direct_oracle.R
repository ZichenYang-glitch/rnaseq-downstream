#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: run_airway_direct_oracle.R R_LIBRARY COUNTS_TSV",
      "METADATA_TSV OUTPUT_DIRECTORY"
    ),
    call. = FALSE
  )
}

r_library <- normalizePath(args[[1L]], mustWork = TRUE)
counts_path <- normalizePath(args[[2L]], mustWork = TRUE)
metadata_path <- normalizePath(args[[3L]], mustWork = TRUE)
output_directory <- normalizePath(args[[4L]], mustWork = FALSE)
.libPaths(c(r_library, .libPaths()))

suppressPackageStartupMessages(library(edgeR))
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

counts_frame <- read.delim(
  counts_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)
metadata <- read.delim(
  metadata_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (!identical(colnames(counts_frame)[-1L], metadata$sample_id)) {
  stop("Oracle input samples do not align exactly.", call. = FALSE)
}
if (anyDuplicated(counts_frame$gene_id) || anyDuplicated(metadata$sample_id)) {
  stop("Oracle input identifiers must be unique.", call. = FALSE)
}

count_matrix <- as.matrix(counts_frame[-1L])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- counts_frame$gene_id
metadata$cell <- factor(metadata$cell, levels = sort(unique(metadata$cell)))
metadata$dex <- factor(metadata$dex, levels = c("untrt", "trt"))
if (anyNA(metadata$cell) || anyNA(metadata$dex)) {
  stop("Oracle metadata contains an undeclared factor level.", call. = FALSE)
}

# This is deliberately a direct edgeR implementation, independent from the
# toolkit adapter. Keep the route literal so an adapter regression cannot alter
# both sides of the parity comparison.
design <- model.matrix(~ cell + dex, data = metadata)
y <- DGEList(counts = count_matrix)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)
fit <- glmQLFit(
  y,
  design,
  abundance.trend = TRUE,
  robust = TRUE,
  winsor.tail.p = c(0.05, 0.1),
  legacy = FALSE,
  top.proportion = NULL,
  keep.unit.mat = FALSE
)
test <- glmQLFTest(fit, coef = "dextrt")
table <- topTags(test, n = Inf, sort.by = "none")$table

results <- data.frame(
  gene_id = rownames(table),
  logFC = table$logFC,
  logCPM = table$logCPM,
  F = table$F,
  PValue = table$PValue,
  FDR = table$FDR,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
coefficients <- data.frame(
  gene_id = rownames(fit$coefficients),
  fit$coefficients,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
old_options <- options(digits = 17)
on.exit(options(old_options), add = TRUE)
write.table(
  results,
  file.path(output_directory, "results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  coefficients,
  file.path(output_directory, "coefficients.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

diagnostics <- list(
  implementation = "independent-direct-edgeR",
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  edgeR_version = as.character(packageVersion("edgeR")),
  design = "~ cell + dex",
  design_columns = colnames(design),
  tested_coefficient = "dextrt",
  input_gene_count = nrow(count_matrix),
  tested_gene_count = nrow(table),
  filtered_gene_count = nrow(count_matrix) - nrow(table)
)
jsonlite::write_json(
  diagnostics,
  file.path(output_directory, "diagnostics.json"),
  auto_unbox = TRUE,
  pretty = FALSE,
  digits = NA,
  na = "null"
)
