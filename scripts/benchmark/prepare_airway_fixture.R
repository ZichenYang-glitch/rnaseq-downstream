#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: prepare_airway_fixture.R R_LIBRARY OUTPUT_DIRECTORY", call. = FALSE)
}

r_library <- normalizePath(args[[1L]], mustWork = TRUE)
output_directory <- normalizePath(args[[2L]], mustWork = FALSE)
.libPaths(c(r_library, .libPaths()))

suppressPackageStartupMessages(library(airway))
suppressPackageStartupMessages(library(SummarizedExperiment))

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
data("airway", package = "airway")

count_matrix <- assay(airway)
metadata <- as.data.frame(colData(airway), stringsAsFactors = FALSE)

if (is.null(rownames(count_matrix)) || anyDuplicated(rownames(count_matrix))) {
  stop("The locked airway assay does not have unique gene identifiers.", call. = FALSE)
}
if (!identical(colnames(count_matrix), rownames(metadata))) {
  stop("The locked airway assay and metadata samples do not align.", call. = FALSE)
}
if (!all(c("cell", "dex") %in% colnames(metadata))) {
  stop("The locked airway metadata is missing cell or dex.", call. = FALSE)
}

counts_output <- data.frame(
  gene_id = rownames(count_matrix),
  count_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata_output <- data.frame(
  sample_id = rownames(metadata),
  cell = as.character(metadata$cell),
  dex = as.character(metadata$dex),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.table(
  counts_output,
  file.path(output_directory, "counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  metadata_output,
  file.path(output_directory, "metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

fixture <- list(
  fixture_id = "airway-1.32.0-cell-plus-dex",
  assay = "airway::airway",
  airway_version = as.character(packageVersion("airway")),
  gene_count = nrow(count_matrix),
  sample_count = ncol(count_matrix),
  design = "~ cell + dex",
  factor_levels = list(
    cell = sort(unique(as.character(metadata$cell))),
    dex = c("untrt", "trt")
  ),
  coefficient = "dextrt"
)
jsonlite::write_json(
  fixture,
  file.path(output_directory, "fixture.json"),
  auto_unbox = TRUE,
  pretty = FALSE,
  digits = NA,
  na = "null"
)
