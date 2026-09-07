#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop(
        paste(
            "Usage: run_deseq2_compcoder_mechanism_diagnostic.R",
            "R_LIBRARY REPLICATE FIXTURE_DIRECTORY OUTPUT_DIRECTORY"
        ),
        call. = FALSE
    )
}

r_library <- normalizePath(args[[1L]], mustWork = TRUE)
replicate_id <- suppressWarnings(as.integer(args[[2L]]))
fixture_directory <- normalizePath(args[[3L]], mustWork = TRUE)
output_directory <- normalizePath(args[[4L]], mustWork = FALSE)
.libPaths(c(r_library, .libPaths()))

if (is.na(replicate_id) || replicate_id < 1L || replicate_id > 20L) {
    stop("Replicate must be an integer from 1 through 20.", call. = FALSE)
}

suppressPackageStartupMessages({
    library(compcodeR)
    library(DESeq2)
    library(edgeR)
})

seed <- 61000L + replicate_id
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

read_fixture <- function() {
    count_frame <- utils::read.delim(
        file.path(fixture_directory, "counts.tsv"),
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    metadata <- utils::read.delim(
        file.path(fixture_directory, "metadata.tsv"),
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    truth <- utils::read.delim(
        file.path(fixture_directory, "truth.tsv"),
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    if (!identical(colnames(count_frame)[[1L]], "gene_id") ||
        !identical(colnames(truth), c(
            "gene_id", "differential_expression", "true_log2_fold_change"
        )) ||
        !identical(colnames(metadata), c("sample_id", "condition"))) {
        stop("The exploratory fixture has an unexpected schema.", call. = FALSE)
    }
    gene_ids <- as.character(count_frame$gene_id)
    counts <- as.matrix(count_frame[, -1L, drop = FALSE])
    storage.mode(counts) <- "integer"
    rownames(counts) <- gene_ids
    if (!identical(colnames(counts), as.character(metadata$sample_id)) ||
        !identical(gene_ids, as.character(truth$gene_id))) {
        stop("The exploratory fixture is not identically aligned.", call. = FALSE)
    }
    metadata$condition <- factor(
        metadata$condition,
        levels = c("control", "treated")
    )
    rownames(metadata) <- metadata$sample_id
    list(counts = counts, metadata = metadata, truth = truth)
}

regenerate_truth <- function(fixture) {
    set.seed(seed)
    simulation <- suppressMessages(suppressWarnings(generateSyntheticData(
        dataset = "p1_deseq2_nb_balanced_exploratory",
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
    annotations <- simulation@variable.annotations
    regenerated_counts <- simulation@count.matrix
    if (!identical(dimnames(regenerated_counts), dimnames(fixture$counts)) ||
        !isTRUE(all.equal(
            unname(regenerated_counts),
            unname(fixture$counts),
            tolerance = 0,
            check.attributes = FALSE
        ))) {
        stop(
            "Regenerated exploratory counts differ from the archived fixture.",
            call. = FALSE
        )
    }
    if (!identical(
        as.integer(annotations$differential.expression),
        as.integer(fixture$truth$differential_expression)
    ) || !isTRUE(all.equal(
        as.numeric(annotations$truelog2foldchanges),
        as.numeric(fixture$truth$true_log2_fold_change),
        tolerance = 1e-12,
        check.attributes = FALSE
    ))) {
        stop(
            "Regenerated exploratory truth differs from the archived fixture.",
            call. = FALSE
        )
    }
    if (!isTRUE(all.equal(
        annotations$truedispersions.S1,
        annotations$truedispersions.S2,
        tolerance = 0,
        check.attributes = FALSE
    ))) {
        stop("The frozen simulation unexpectedly has group-specific dispersion.", call. = FALSE)
    }
    as.numeric(annotations$truedispersions.S1)
}

capture_notices <- function(expression) {
    notices <- list(warnings = character(), messages = character())
    value <- withCallingHandlers(
        expression,
        warning = function(condition) {
            notices$warnings <<- c(notices$warnings, conditionMessage(condition))
            invokeRestart("muffleWarning")
        },
        message = function(condition) {
            notices$messages <<- c(notices$messages, conditionMessage(condition))
            invokeRestart("muffleMessage")
        }
    )
    list(value = value, notices = notices)
}

run_deseq <- function(fixture, mode, full_design, reduced_design) {
    dataset <- DESeq2::DESeqDataSetFromMatrix(
        countData = fixture$counts,
        colData = S4Vectors::DataFrame(fixture$metadata),
        design = full_design
    )
    arguments <- list(
        object = dataset,
        test = if (identical(mode, "wald")) "Wald" else "LRT",
        fitType = "parametric",
        sfType = "ratio",
        betaPrior = FALSE,
        quiet = TRUE,
        minReplicatesForReplace = 7,
        useT = FALSE,
        minmu = 0.5,
        parallel = FALSE
    )
    if (identical(mode, "lrt")) {
        arguments$full <- full_design
        arguments$reduced <- reduced_design
    }
    captured <- capture_notices(do.call(DESeq2::DESeq, arguments))
    list(dataset = captured$value, notices = captured$notices)
}

extract_results <- function(dataset, contrast, independent_filtering) {
    as.data.frame(DESeq2::results(
        dataset,
        contrast = contrast,
        lfcThreshold = 0,
        altHypothesis = "greaterAbs",
        independentFiltering = independent_filtering,
        alpha = 0.1,
        pAdjustMethod = "BH",
        parallel = FALSE
    ))
}

fixture <- read_fixture()
true_dispersion <- regenerate_truth(fixture)
full_design <- stats::model.matrix(~ condition, data = fixture$metadata)
reduced_design <- stats::model.matrix(~ 1, data = fixture$metadata)
contrast <- c(0, 1)

wald <- run_deseq(fixture, "wald", full_design, reduced_design)
wald_if_on <- extract_results(wald$dataset, contrast, TRUE)
wald_if_off <- extract_results(wald$dataset, contrast, FALSE)
lrt <- run_deseq(fixture, "lrt", full_design, reduced_design)
lrt_if_on <- extract_results(lrt$dataset, contrast, TRUE)

y_all <- edgeR::DGEList(
    counts = fixture$counts,
    genes = data.frame(
        gene_id = rownames(fixture$counts),
        stringsAsFactors = FALSE
    )
)
keep <- edgeR::filterByExpr(y_all, design = full_design)
if (!any(keep)) {
    stop("edgeR::filterByExpr removed every gene.", call. = FALSE)
}
y <- y_all[keep, , keep.lib.sizes = FALSE]
y <- edgeR::normLibSizes(y, method = "TMM")
edge_fit <- edgeR::glmQLFit(
    y,
    design = full_design,
    abundance.trend = TRUE,
    robust = TRUE,
    winsor.tail.p = c(0.05, 0.1),
    legacy = FALSE,
    top.proportion = NULL,
    keep.unit.mat = FALSE
)
edge_test <- edgeR::glmQLFTest(edge_fit, contrast = contrast)
edge_pvalue <- rep(NA_real_, nrow(fixture$counts))
edge_fdr <- rep(NA_real_, nrow(fixture$counts))
names(edge_pvalue) <- rownames(fixture$counts)
names(edge_fdr) <- rownames(fixture$counts)
edge_ids <- rownames(edge_test$table)
edge_pvalue[edge_ids] <- edge_test$table$PValue
edge_fdr[edge_ids] <- stats::p.adjust(edge_test$table$PValue, method = "BH")

wald_columns <- S4Vectors::mcols(wald$dataset)
fit_function <- DESeq2::dispersionFunction(wald$dataset)
fit_type <- attr(fit_function, "fitType")
fit_coefficients <- attr(fit_function, "coefficients")
lrt_fit_function <- DESeq2::dispersionFunction(lrt$dataset)
lrt_fit_type <- attr(lrt_fit_function, "fitType")
lrt_fit_coefficients <- attr(lrt_fit_function, "coefficients")
if (!fit_type %in% c("parametric", "local", "mean") ||
    !lrt_fit_type %in% c("parametric", "local", "mean")) {
    stop("DESeq2 returned an unrecognized dispersion fit type.", call. = FALSE)
}
assert_fit_coefficients <- function(values, resolved_type, mode) {
    if (identical(resolved_type, "parametric") &&
        (is.null(values) || length(values) != 2L ||
            any(!is.finite(values)) || any(values <= 0))) {
        stop(
            sprintf("The %s parametric fit coefficients are invalid.", mode),
            call. = FALSE
        )
    }
}
assert_fit_coefficients(fit_coefficients, fit_type, "Wald")
assert_fit_coefficients(lrt_fit_coefficients, lrt_fit_type, "LRT")
max_dispersion_bound <- max(10, ncol(wald$dataset))

gene_table <- data.frame(
    gene_id = rownames(fixture$counts),
    truth_class = ifelse(
        fixture$truth$differential_expression == 1L, "DE", "null"
    ),
    true_dispersion = true_dispersion,
    base_mean = as.numeric(wald_if_on$baseMean),
    dispersion_gene_estimate = as.numeric(wald_columns$dispGeneEst),
    dispersion_fitted_trend = as.numeric(wald_columns$dispFit),
    dispersion_final = as.numeric(wald_columns$dispersion),
    dispersion_outlier = as.integer(wald_columns$dispOutlier %in% TRUE),
    true_dispersion_above_deseq2_max_disp = as.integer(
        true_dispersion > max_dispersion_bound
    ),
    wald_pvalue = as.numeric(wald_if_on$pvalue),
    wald_fdr_if_on = as.numeric(wald_if_on$padj),
    wald_fdr_if_off = as.numeric(wald_if_off$padj),
    lrt_pvalue = as.numeric(lrt_if_on$pvalue),
    lrt_fdr_if_on = as.numeric(lrt_if_on$padj),
    edger_tested = as.integer(keep),
    edger_pvalue = unname(edge_pvalue),
    edger_fdr_native = unname(edge_fdr),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

old_options <- options(digits = 17)
on.exit(options(old_options), add = TRUE)
utils::write.table(
    gene_table,
    file.path(output_directory, "gene-diagnostics.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "",
    eol = "\n",
    fileEncoding = "UTF-8"
)

fit_summary <- list(
    replicate = replicate_id,
    seed = seed,
    fixture_regeneration_match = TRUE,
    requested_fit_type = "parametric",
    resolved_fit_type = fit_type,
    lrt_resolved_fit_type = lrt_fit_type,
    automatic_fallback = !identical(fit_type, "parametric") ||
        !identical(lrt_fit_type, "parametric"),
    fit_coefficients = if (is.null(fit_coefficients)) {
        NULL
    } else {
        as.list(as.numeric(fit_coefficients)) |>
            stats::setNames(names(fit_coefficients))
    },
    lrt_fit_coefficients = if (is.null(lrt_fit_coefficients)) {
        NULL
    } else {
        as.list(as.numeric(lrt_fit_coefficients)) |>
            stats::setNames(names(lrt_fit_coefficients))
    },
    max_dispersion_bound = max_dispersion_bound,
    dispersion_outlier_count = sum(wald_columns$dispOutlier %in% TRUE),
    true_dispersion_above_max_count = sum(
        true_dispersion > max_dispersion_bound
    ),
    wald_notices = wald$notices,
    lrt_notices = lrt$notices,
    n_genes = nrow(gene_table),
    n_truth_null = sum(gene_table$truth_class == "null"),
    n_truth_de = sum(gene_table$truth_class == "DE"),
    n_edger_tested = sum(keep),
    runtime = list(
        R = as.character(getRversion()),
        DESeq2 = as.character(utils::packageVersion("DESeq2")),
        edgeR = as.character(utils::packageVersion("edgeR")),
        compcodeR = as.character(utils::packageVersion("compcodeR"))
    )
)
jsonlite::write_json(
    fit_summary,
    file.path(output_directory, "fit-summary.json"),
    auto_unbox = TRUE,
    pretty = FALSE,
    digits = NA,
    na = "null"
)
