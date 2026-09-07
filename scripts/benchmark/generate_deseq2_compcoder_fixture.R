#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: generate_deseq2_compcoder_fixture.R",
      "R_LIBRARY MODE REPLICATE OUTPUT_DIRECTORY"
    ),
    call. = FALSE
  )
}

r_library <- normalizePath(args[[1L]], mustWork = TRUE)
mode <- args[[2L]]
replicate_id <- suppressWarnings(as.integer(args[[3L]]))
output_directory <- normalizePath(args[[4L]], mustWork = FALSE)
.libPaths(c(r_library, .libPaths()))

if (!(mode %in% c("exploratory", "certification"))) {
  stop("Mode must be exactly 'exploratory' or 'certification'.", call. = FALSE)
}
if (is.na(replicate_id) || replicate_id < 1L || replicate_id > 20L) {
  stop("Replicate must be an integer from 1 through 20.", call. = FALSE)
}

suppressPackageStartupMessages(library(compcodeR))
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

seed_base <- if (identical(mode, "exploratory")) 61000L else 62000L
seed <- seed_base + replicate_id
set.seed(seed)
simulation <- suppressMessages(suppressWarnings(generateSyntheticData(
  dataset = sprintf("p1_deseq2_nb_balanced_%s", mode),
  n.vars = 5000,
  samples.per.cond = 6,
  n.diffexp = 500,
  repl.id = replicate_id,
  fraction.upregulated = 0.5,
  between.group.diffdisp = FALSE,
  filter.threshold.total = 0,
  random.outlier.high.prob = 0,
  random.outlier.low.prob = 0,
  single.outlier.high.prob = 0,
  single.outlier.low.prob = 0,
  effect.size = 2
)))

count_matrix <- simulation@count.matrix
sample_annotations <- simulation@sample.annotations
truth_annotations <- simulation@variable.annotations
if (anyDuplicated(rownames(count_matrix))) {
  stop("compcodeR generated duplicate gene identifiers.", call. = FALSE)
}
if (any(count_matrix < 0) || any(count_matrix != round(count_matrix))) {
  stop("compcodeR generated values outside the integer count domain.", call. = FALSE)
}

counts_output <- data.frame(
  gene_id = rownames(count_matrix),
  count_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata_output <- data.frame(
  sample_id = colnames(count_matrix),
  condition = ifelse(sample_annotations$condition == 1, "control", "treated"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
truth_output <- data.frame(
  gene_id = rownames(count_matrix),
  differential_expression = as.integer(truth_annotations$differential.expression),
  true_log2_fold_change = truth_annotations$truelog2foldchanges,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
old_options <- options(digits = 17)
on.exit(options(old_options), add = TRUE)

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
write.table(
  truth_output,
  file.path(output_directory, "truth.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

fixture <- list(
  fixture_id = sprintf(
    "compcodeR-1.48.0-p1-deseq2-nb-balanced-%s-replicate-%02d",
    mode,
    replicate_id
  ),
  mode = mode,
  compcodeR_version = as.character(packageVersion("compcodeR")),
  seed = seed,
  replicate = replicate_id,
  n_genes = nrow(count_matrix),
  samples_per_condition = 6L,
  n_differential = sum(truth_annotations$differential.expression == 1),
  n_upregulated = sum(truth_annotations$upregulation == 1),
  n_downregulated = sum(truth_annotations$downregulation == 1),
  effect_size = 2,
  fraction_upregulated = 0.5,
  filter_threshold_total = 0,
  outlier_probabilities = list(
    random_high = 0,
    random_low = 0,
    single_high = 0,
    single_low = 0
  ),
  design = "~ condition",
  factor_levels = list(condition = c("control", "treated")),
  coefficient = "conditiontreated",
  input_route = "featurecounts_integer_combined_matrix_self_attested_simulation"
)
jsonlite::write_json(
  fixture,
  file.path(output_directory, "fixture.json"),
  auto_unbox = TRUE,
  pretty = FALSE,
  digits = NA,
  na = "null"
)
