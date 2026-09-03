#!/usr/bin/env Rscript

# Generate the frozen C2 negative-binomial fixtures and outcome-independent
# gene-set memberships used by the pathway regression gate.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: generate_pathway_compcoder_fixture.R R_LIBRARY",
      "SCENARIO REPLICATE OUTPUT_DIRECTORY"
    ),
    call. = FALSE
  )
}

r_library <- normalizePath(args[[1L]], mustWork = TRUE)
scenario <- args[[2L]]
replicate_id <- suppressWarnings(as.integer(args[[3L]]))
output_directory <- normalizePath(args[[4L]], mustWork = FALSE)
.libPaths(c(r_library, .libPaths()))

if (!scenario %in% c("mixed", "complete_null")) {
  stop("Scenario must be 'mixed' or 'complete_null'.", call. = FALSE)
}
maximum_replicate <- if (identical(scenario, "mixed")) 20L else 40L
if (is.na(replicate_id) || replicate_id < 1L || replicate_id > maximum_replicate) {
  stop(
    sprintf("Replicate must be an integer from 1 through %d.", maximum_replicate),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(compcodeR))
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
simulation_seed <- if (identical(scenario, "mixed")) {
  5000L + replicate_id
} else {
  15000L + replicate_id
}
membership_seed <- if (identical(scenario, "mixed")) {
  25000L + replicate_id
} else {
  45000L + replicate_id
}
n_differential <- if (identical(scenario, "mixed")) 500L else 0L

set.seed(simulation_seed)
simulation <- suppressMessages(suppressWarnings(generateSyntheticData(
  dataset = paste0("c2_pathway_", scenario),
  n.vars = 5000,
  samples.per.cond = 6,
  n.diffexp = n_differential,
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
gene_ids <- rownames(count_matrix)
if (is.null(gene_ids) || anyDuplicated(gene_ids)) {
  stop("compcodeR generated invalid gene identifiers.", call. = FALSE)
}

symbols <- sprintf("SYM%05d", seq_along(gene_ids))
names(symbols) <- gene_ids
differential <- as.integer(truth_annotations$differential.expression)
true_lfc <- as.numeric(truth_annotations$truelog2foldchanges)
if (length(differential) != length(gene_ids) || length(true_lfc) != length(gene_ids)) {
  stop("compcodeR truth annotations do not align with the count matrix.", call. = FALSE)
}

# Membership uses a distinct RNG stream and only frozen truth classes. It never
# inspects observed counts, fitted statistics, p-values, or gene rank.
set.seed(membership_seed)
null_symbols <- unname(symbols[differential == 0L])
up_symbols <- unname(symbols[differential == 1L & true_lfc > 0])
down_symbols <- unname(symbols[differential == 1L & true_lfc < 0])

set_records <- list()
add_set <- function(identifier, members, truth_class, true_direction) {
  if (length(members) == 0L || anyNA(members) || anyDuplicated(members)) {
    stop(sprintf("Generated invalid membership for %s.", identifier), call. = FALSE)
  }
  set_records[[length(set_records) + 1L]] <<- list(
    gene_set_id = identifier,
    description = paste0("frozen_", scenario, "_", identifier),
    truth_class = truth_class,
    true_direction = true_direction,
    nominal_size = length(members),
    members = as.character(members)
  )
}

for (size in c(20L, 40L, 80L, 160L)) {
  for (set_index in seq_len(25L)) {
    add_set(
      sprintf("NULL_%03d_%02d", size, set_index),
      sample(null_symbols, size = size, replace = FALSE),
      "null",
      "None"
    )
  }
}

if (identical(scenario, "mixed")) {
  if (length(up_symbols) != 250L || length(down_symbols) != 250L) {
    stop("The mixed scenario does not contain exactly 250 Up and 250 Down genes.", call. = FALSE)
  }
  selected_up <- sample(up_symbols, size = 200L, replace = FALSE)
  selected_down <- sample(down_symbols, size = 200L, replace = FALSE)
  for (set_index in seq_len(5L)) {
    range <- ((set_index - 1L) * 40L + 1L):(set_index * 40L)
    add_set(
      sprintf("POSITIVE_UP_%02d", set_index),
      selected_up[range],
      "positive",
      "Up"
    )
  }
  for (set_index in seq_len(5L)) {
    range <- ((set_index - 1L) * 40L + 1L):(set_index * 40L)
    add_set(
      sprintf("POSITIVE_DOWN_%02d", set_index),
      selected_down[range],
      "positive",
      "Down"
    )
  }
}

expected_set_count <- if (identical(scenario, "mixed")) 110L else 100L
if (length(set_records) != expected_set_count) {
  stop("Generated gene-set count differs from the frozen scenario.", call. = FALSE)
}
set_record_ids <- vapply(set_records, `[[`, character(1L), "gene_set_id")
set_records <- set_records[order(set_record_ids, method = "radix")]

counts_output <- data.frame(
  gene_id = gene_ids,
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
  gene_id = gene_ids,
  symbol = unname(symbols),
  differential_expression = differential,
  true_log2_fold_change = true_lfc,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotation_output <- data.frame(
  gene_id = gene_ids,
  symbol = unname(symbols),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
membership_output <- do.call(
  rbind,
  lapply(set_records, function(record) {
    data.frame(
      gene_set_id = record$gene_set_id,
      truth_class = record$truth_class,
      true_direction = record$true_direction,
      nominal_size = record$nominal_size,
      symbol = record$members,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)
gmt_lines <- vapply(
  set_records,
  function(record) {
    paste(c(record$gene_set_id, record$description, record$members), collapse = "\t")
  },
  character(1L)
)

old_options <- options(digits = 17)
on.exit(options(old_options), add = TRUE)
write.table(
  counts_output,
  file.path(output_directory, "counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
write.table(
  metadata_output,
  file.path(output_directory, "metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
write.table(
  truth_output,
  file.path(output_directory, "truth.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
write.table(
  annotation_output,
  file.path(output_directory, "annotation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
write.table(
  membership_output,
  file.path(output_directory, "membership.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
writeLines(gmt_lines, file.path(output_directory, "sets.gmt"), sep = "\n")

fixture <- list(
  fixture_id = sprintf(
    "compcodeR-1.48.0-c2-pathway-%s-replicate-%02d",
    scenario,
    replicate_id
  ),
  scenario = scenario,
  compcodeR_version = as.character(packageVersion("compcodeR")),
  simulation_seed = simulation_seed,
  membership_seed = membership_seed,
  replicate = replicate_id,
  rng_kind = list(
    generator = "Mersenne-Twister",
    normal = "Inversion",
    sample = "Rejection"
  ),
  n_genes = nrow(count_matrix),
  samples_per_condition = 6L,
  n_differential = sum(differential == 1L),
  n_upregulated = sum(differential == 1L & true_lfc > 0),
  n_downregulated = sum(differential == 1L & true_lfc < 0),
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
  null_set_count = 100L,
  null_set_sizes = list(`20` = 25L, `40` = 25L, `80` = 25L, `160` = 25L),
  positive_set_count = if (identical(scenario, "mixed")) 10L else 0L,
  positive_set_size = if (identical(scenario, "mixed")) 40L else 0L,
  minimum_tested_genes = 10L,
  construction_policy = paste(
    "membership_rng_is_separate_from_simulation_rng;",
    "null_sets_sample_truth_null_genes_without_replacement_within_set;",
    "positive_sets_are_disjoint_within_direction_and_sample_truth_DE_genes;",
    "no_observed_count_statistic_pvalue_or_gene_rank_is_used",
    sep = ""
  )
)
jsonlite::write_json(
  fixture,
  file.path(output_directory, "fixture.json"),
  auto_unbox = TRUE,
  pretty = FALSE,
  digits = NA,
  na = "null"
)
