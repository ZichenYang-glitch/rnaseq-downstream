#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop(
        paste(
            "Usage: run_deseq2_airway_direct_oracle.R R_LIBRARY COUNTS_TSV",
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

suppressPackageStartupMessages(library(DESeq2))

counts_frame <- utils::read.delim(
    counts_path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
)
metadata <- utils::read.delim(
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
metadata$cell <- factor(
    metadata$cell,
    levels = c("N052611", "N061011", "N080611", "N61311")
)
metadata$dex <- factor(metadata$dex, levels = c("untrt", "trt"))
if (anyNA(metadata$cell) || anyNA(metadata$dex)) {
    stop("Oracle metadata contains an undeclared factor level.", call. = FALSE)
}
rownames(metadata) <- metadata$sample_id

# This is deliberately a direct DESeq2 implementation. It does not import or
# call toolkit code. Both modes consume the same serialized airway matrix as
# the public toolkit route, but this fixture is not evidence of featureCounts
# provenance.
full_design <- stats::model.matrix(~ cell + dex, data = metadata)
reduced_design <- stats::model.matrix(~ cell, data = metadata)
contrast <- rep(0, ncol(full_design))
names(contrast) <- colnames(full_design)
contrast[["dextrt"]] <- 1

write_tsv <- function(value, path) {
    old_options <- options(digits = 17)
    on.exit(options(old_options), add = TRUE)
    utils::write.table(
        value,
        path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = "",
        eol = "\n",
        fileEncoding = "UTF-8"
    )
}

run_mode <- function(mode) {
    dataset <- DESeq2::DESeqDataSetFromMatrix(
        countData = count_matrix,
        colData = S4Vectors::DataFrame(metadata),
        design = full_design
    )
    if (identical(mode, "wald")) {
        dataset <- DESeq2::DESeq(
            dataset,
            test = "Wald",
            fitType = "parametric",
            sfType = "ratio",
            betaPrior = FALSE,
            quiet = TRUE,
            minReplicatesForReplace = 7,
            useT = FALSE,
            minmu = 0.5,
            parallel = FALSE
        )
        converged <- as.logical(S4Vectors::mcols(dataset)$betaConv)
    } else {
        dataset <- DESeq2::DESeq(
            dataset,
            test = "LRT",
            fitType = "parametric",
            sfType = "ratio",
            betaPrior = FALSE,
            full = full_design,
            reduced = reduced_design,
            quiet = TRUE,
            minReplicatesForReplace = 7,
            useT = FALSE,
            minmu = 0.5,
            parallel = FALSE
        )
        full_converged <- as.logical(S4Vectors::mcols(dataset)$fullBetaConv)
        reduced_converged <- as.logical(
            S4Vectors::mcols(dataset)$reducedBetaConv
        )
        converged <- full_converged & reduced_converged
        converged[is.na(full_converged) | is.na(reduced_converged)] <- NA
    }

    result <- DESeq2::results(
        dataset,
        contrast = contrast,
        lfcThreshold = 0,
        altHypothesis = "greaterAbs",
        independentFiltering = TRUE,
        alpha = 0.1,
        pAdjustMethod = "BH",
        parallel = FALSE
    )
    result_frame <- as.data.frame(result)
    tested <- !is.na(result_frame$pvalue) & !is.na(result_frame$padj)
    tested <- tested & is.finite(result_frame$pvalue)
    tested <- tested & is.finite(result_frame$padj)
    tested <- tested & !is.na(converged) & converged
    tested_results <- data.frame(
        gene_id = rownames(result_frame)[tested],
        logFC = result_frame$log2FoldChange[tested],
        statistic = result_frame$stat[tested],
        PValue = result_frame$pvalue[tested],
        FDR = result_frame$padj[tested],
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    result_names <- DESeq2::resultsNames(dataset)
    if (length(result_names) != ncol(full_design)) {
        stop("Oracle coefficient mapping is not one-to-one.", call. = FALSE)
    }
    coefficient_matrix <- as.matrix(
        S4Vectors::mcols(dataset)[, result_names, drop = FALSE]
    )
    colnames(coefficient_matrix) <- colnames(full_design)
    coefficient_rows <- vector("list", length(coefficient_matrix))
    cursor <- 1L
    for (gene_index in seq_len(nrow(coefficient_matrix))) {
        if (!is.na(converged[[gene_index]]) && converged[[gene_index]]) {
            for (column_index in seq_len(ncol(coefficient_matrix))) {
                estimate <- coefficient_matrix[gene_index, column_index]
                if (is.finite(estimate)) {
                    coefficient_rows[[cursor]] <- data.frame(
                        gene_id = rownames(coefficient_matrix)[[gene_index]],
                        coefficient = colnames(coefficient_matrix)[[column_index]],
                        estimate = estimate,
                        stringsAsFactors = FALSE
                    )
                    cursor <- cursor + 1L
                }
            }
        }
    }
    coefficient_rows <- coefficient_rows[seq_len(cursor - 1L)]
    coefficients <- do.call(rbind, coefficient_rows)

    mode_directory <- file.path(output_directory, mode)
    dir.create(mode_directory, recursive = TRUE, showWarnings = FALSE)
    write_tsv(tested_results, file.path(mode_directory, "results.tsv"))
    write_tsv(coefficients, file.path(mode_directory, "coefficients.tsv"))

    diagnostics <- list(
        implementation = "independent-direct-DESeq2",
        mode = mode,
        r_version = paste(R.version$major, R.version$minor, sep = "."),
        bioconductor_version = as.character(packageVersion("BiocVersion")),
        deseq2_version = as.character(packageVersion("DESeq2")),
        apeglm_version = as.character(packageVersion("apeglm")),
        tximport_version = as.character(packageVersion("tximport")),
        full_design = "~ cell + dex",
        reduced_design = if (identical(mode, "lrt")) "~ cell" else NULL,
        design_columns = colnames(full_design),
        reporting_contrast = as.list(contrast),
        inference = if (identical(mode, "lrt")) {
            "full_vs_reduced_omnibus"
        } else {
            "wald_dextrt_equals_zero"
        },
        reporting_effect = "dextrt",
        input_gene_count = nrow(count_matrix),
        tested_gene_count = nrow(tested_results),
        coefficient_value_count = nrow(coefficients)
    )
    jsonlite::write_json(
        diagnostics,
        file.path(mode_directory, "diagnostics.json"),
        auto_unbox = TRUE,
        pretty = FALSE,
        digits = NA,
        na = "null"
    )
}

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
run_mode("wald")
run_mode("lrt")
