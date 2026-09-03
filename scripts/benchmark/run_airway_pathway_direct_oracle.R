#!/usr/bin/env Rscript

# Independent C2 same-engine oracle. This script deliberately does not source
# or call any toolkit R helper; the complete edgeR/limma route is written here.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    paste(
      "Usage: run_airway_pathway_direct_oracle.R R_LIBRARY COUNTS_TSV",
      "METADATA_TSV GMT ANNOTATION MINIMUM_TESTED_GENES OUTPUT_DIRECTORY"
    ),
    call. = FALSE
  )
}

r_library <- normalizePath(args[[1L]], mustWork = TRUE)
counts_path <- normalizePath(args[[2L]], mustWork = TRUE)
metadata_path <- normalizePath(args[[3L]], mustWork = TRUE)
gmt_path <- normalizePath(args[[4L]], mustWork = TRUE)
annotation_path <- normalizePath(args[[5L]], mustWork = TRUE)
minimum_tested_genes <- suppressWarnings(as.integer(args[[6L]]))
output_directory <- normalizePath(args[[7L]], mustWork = FALSE)
.libPaths(c(r_library, .libPaths()))

if (is.na(minimum_tested_genes) || minimum_tested_genes < 1L) {
  stop("MINIMUM_TESTED_GENES must be a positive integer.", call. = FALSE)
}

suppressPackageStartupMessages(library(edgeR))
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

counts_frame <- read.delim(
  counts_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
metadata <- read.delim(
  metadata_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
annotation <- read.delim(
  annotation_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
if (!identical(colnames(counts_frame)[-1L], metadata$sample_id)) {
  stop("Oracle input samples do not align exactly.", call. = FALSE)
}
if (!identical(colnames(annotation), c("gene_id", "symbol"))) {
  stop("Oracle annotation does not have the frozen two-column schema.", call. = FALSE)
}
if (anyDuplicated(counts_frame$gene_id) || anyDuplicated(metadata$sample_id) ||
    anyDuplicated(annotation$gene_id) || anyDuplicated(annotation$symbol)) {
  stop("Oracle identifiers must be unique for this one-to-one fixture.", call. = FALSE)
}

gmt_lines <- readLines(gmt_path, warn = FALSE, encoding = "UTF-8")
if (length(gmt_lines) == 0L || any(!nzchar(gmt_lines))) {
  stop("Oracle GMT must contain non-empty records.", call. = FALSE)
}
gmt_fields <- strsplit(gmt_lines, "\t", fixed = TRUE)
if (any(lengths(gmt_fields) < 3L)) {
  stop("Oracle GMT records must contain an identifier, description, and members.", call. = FALSE)
}
gene_set_ids <- vapply(gmt_fields, `[[`, character(1L), 1L)
gene_set_descriptions <- vapply(gmt_fields, `[[`, character(1L), 2L)
if (anyDuplicated(gene_set_ids) || !identical(gene_set_ids, sort(gene_set_ids))) {
  stop("Oracle GMT identifiers must be unique and lexically ordered.", call. = FALSE)
}
gene_set_symbols <- lapply(gmt_fields, function(fields) fields[-c(1L, 2L)])
if (any(vapply(gene_set_symbols, anyDuplicated, integer(1L)) > 0L)) {
  stop("Oracle GMT records may not contain duplicate members.", call. = FALSE)
}
names(gene_set_descriptions) <- gene_set_ids
names(gene_set_symbols) <- gene_set_ids

count_matrix <- as.matrix(counts_frame[-1L])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- counts_frame$gene_id
metadata$cell <- factor(metadata$cell, levels = sort(unique(metadata$cell)))
metadata$dex <- factor(metadata$dex, levels = c("untrt", "trt"))
if (anyNA(metadata$cell) || anyNA(metadata$dex)) {
  stop("Oracle metadata contains an undeclared factor level.", call. = FALSE)
}

# Literal independent edgeR v4 QL path, kept aligned to the locked production
# defaults so the pathway comparison can detect a backend regression.
design <- model.matrix(~ cell + dex, data = metadata)
y <- DGEList(counts = count_matrix)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y, method = "TMM")
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
contrast <- rep(0, ncol(design))
names(contrast) <- colnames(design)
contrast[["dextrt"]] <- 1

symbol_to_gene <- setNames(annotation$gene_id, annotation$symbol)
all_gene_ids <- rownames(count_matrix)
tested_gene_ids <- rownames(fit$coefficients)
mapping <- lapply(gene_set_ids, function(identifier) {
  symbols <- gene_set_symbols[[identifier]]
  mapped <- unname(symbol_to_gene[symbols])
  mapped <- mapped[!is.na(mapped)]
  mapped <- sort(unique(mapped))
  tested <- sort(intersect(mapped, tested_gene_ids))
  list(
    symbols = symbols,
    mapped = mapped,
    tested = tested,
    filtered = sort(intersect(mapped, setdiff(all_gene_ids, tested_gene_ids)))
  )
})
names(mapping) <- gene_set_ids
eligible_ids <- gene_set_ids[vapply(
  mapping,
  function(record) length(record$tested) >= minimum_tested_genes,
  logical(1L)
)]
indices <- lapply(
  eligible_ids,
  function(identifier) match(mapping[[identifier]]$tested, tested_gene_ids)
)
names(indices) <- eligible_ids
if (length(indices) == 0L || any(vapply(indices, anyNA, logical(1L)))) {
  stop("Oracle fixture has no valid self-contained gene sets.", call. = FALSE)
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(1729L)
roast <- limma::mroast(
  fit,
  index = indices,
  design = design,
  contrast = contrast,
  geneid = rownames(fit$counts),
  set.statistic = "mean",
  gene.weights = NULL,
  nrot = 9999L,
  adjust.method = "BH",
  midp = FALSE,
  sort = "none"
)
fry_result <- limma::fry(
  fit,
  index = indices,
  design = design,
  contrast = contrast,
  geneid = rownames(fit$counts),
  sort = "none"
)
camera_ids <- eligible_ids[vapply(
  mapping[eligible_ids],
  function(record) {
    length(record$tested) >= 2L && length(record$tested) < length(tested_gene_ids)
  },
  logical(1L)
)]
camera_indices <- indices[camera_ids]
camera_result <- limma::camera(
  fit,
  index = camera_indices,
  design = design,
  contrast = contrast,
  weights = NULL,
  use.ranks = FALSE,
  allow.neg.cor = FALSE,
  inter.gene.cor = NA_real_,
  sort = FALSE
)

rows <- list()
add_row <- function(
    gene_set_id, method_id, hypothesis, direction, p_value, fdr,
    method_ngenes, proportion_down = NA_real_, proportion_up = NA_real_,
    correlation_raw = NA_real_, correlation_effective = NA_real_,
    vif_used = NA_real_) {
  rows[[length(rows) + 1L]] <<- data.frame(
    contrast_id = "trt_vs_untrt",
    gene_set_id = gene_set_id,
    method_id = method_id,
    hypothesis = hypothesis,
    direction = direction,
    proportion_down = proportion_down,
    proportion_up = proportion_up,
    p_value = p_value,
    fdr = fdr,
    method_ngenes = as.integer(method_ngenes),
    correlation_estimate_raw = correlation_raw,
    correlation_effective = correlation_effective,
    vif_used = vif_used,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
for (identifier in eligible_ids) {
  add_row(
    identifier, "limma_mroast", "directional", roast[identifier, "Direction"],
    roast[identifier, "PValue"], roast[identifier, "FDR"],
    roast[identifier, "NGenes"], roast[identifier, "PropDown"],
    roast[identifier, "PropUp"]
  )
  add_row(
    identifier, "limma_mroast", "mixed", "Mixed", roast[identifier, "PValue.Mixed"],
    roast[identifier, "FDR.Mixed"], roast[identifier, "NGenes"],
    roast[identifier, "PropDown"], roast[identifier, "PropUp"]
  )
  add_row(
    identifier, "limma_fry", "directional", fry_result[identifier, "Direction"],
    fry_result[identifier, "PValue"], fry_result[identifier, "FDR"],
    fry_result[identifier, "NGenes"]
  )
  add_row(
    identifier, "limma_fry", "mixed", "Mixed", fry_result[identifier, "PValue.Mixed"],
    fry_result[identifier, "FDR.Mixed"], fry_result[identifier, "NGenes"]
  )
  if (identifier %in% camera_ids) {
    raw_correlation <- camera_result[identifier, "Correlation"]
    effective_correlation <- max(raw_correlation, 0)
    n_genes <- camera_result[identifier, "NGenes"]
    add_row(
      identifier, "limma_camera", "directional",
      camera_result[identifier, "Direction"],
      camera_result[identifier, "PValue"],
      camera_result[identifier, "FDR"], n_genes,
      correlation_raw = raw_correlation,
      correlation_effective = effective_correlation,
      vif_used = 1 + (n_genes - 1) * effective_correlation
    )
  }
}
oracle <- do.call(rbind, rows)
oracle <- oracle[order(
  oracle$contrast_id,
  oracle$gene_set_id,
  oracle$method_id,
  oracle$hypothesis
), , drop = FALSE]

old_options <- options(digits = 17)
on.exit(options(old_options), add = TRUE)
write.table(
  oracle,
  file.path(output_directory, "pathway_oracle.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = ""
)

diagnostics <- list(
  implementation = "independent-direct-edgeR-limma",
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  edgeR_version = as.character(packageVersion("edgeR")),
  limma_version = as.character(packageVersion("limma")),
  design = "~ cell + dex",
  design_columns = colnames(design),
  contrast_id = "trt_vs_untrt",
  contrast = as.list(contrast),
  input_gene_count = nrow(count_matrix),
  tested_universe_gene_count = length(tested_gene_ids),
  gene_set_count = length(gene_set_ids),
  eligible_gene_set_count = length(eligible_ids),
  mroast = list(
    seed = 1729L,
    nrot = 9999L,
    rng_kind = list(
      generator = "Mersenne-Twister",
      normal = "Inversion",
      sample = "Rejection"
    ),
    set_statistic = "mean",
    midp = FALSE,
    adjust_method = "BH"
  ),
  camera = list(
    inter_gene_correlation = "estimated_set_specific",
    use_ranks = FALSE,
    allow_negative_correlation = FALSE
  )
)
jsonlite::write_json(
  diagnostics,
  file.path(output_directory, "diagnostics.json"),
  auto_unbox = TRUE, pretty = FALSE, digits = NA, na = "null"
)
