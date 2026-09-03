#!/usr/bin/env Rscript

# Locked edgeR v4 quasi-likelihood kernel.  Python owns input evidence,
# orchestration, publication, and error mapping; this process owns every
# statistical operation.

EXPECTED_RUNTIME <- list(
    R = "4.6.1",
    Bioconductor = "3.23",
    edgeR = "4.10.0",
    tximport = "1.40.0",
    limma = "3.68.0"
)
QR_TOLERANCE <- 1e-7
ACTIVE_SCHEMA_VERSION <- "1.0"
STATUS_VOCABULARY <- c(
    "filtered", "not_tested", "not_estimable", "failed", "tested"
)

backend_abort <- function(code, message, details = list(), exit_status = 4L) {
    condition <- structure(
        list(
            message = as.character(message),
            call = NULL,
            code = as.character(code),
            details = details,
            exit_status = as.integer(exit_status)
        ),
        class = c("rnaseq_backend_error", "error", "condition")
    )
    stop(condition)
}

emit_document <- function(document) {
    payload <- jsonlite::toJSON(
        document,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = NA,
        pretty = FALSE
    )
    cat(payload, "\n", sep = "")
    flush.console()
}

runtime_identity <- function() {
    bioc_version <- as.character(utils::packageVersion("BiocVersion"))
    list(
        R = paste(R.version$major, R.version$minor, sep = "."),
        Bioconductor = paste(strsplit(bioc_version, "\\.")[[1L]][1:2], collapse = "."),
        BiocVersion_package = bioc_version,
        edgeR = as.character(utils::packageVersion("edgeR")),
        tximport = as.character(utils::packageVersion("tximport")),
        limma = as.character(utils::packageVersion("limma"))
    )
}

assert_runtime <- function() {
    observed <- runtime_identity()
    mismatches <- list()
    for (name in names(EXPECTED_RUNTIME)) {
        if (!identical(observed[[name]], EXPECTED_RUNTIME[[name]])) {
            mismatches[[name]] <- list(
                expected = EXPECTED_RUNTIME[[name]],
                observed = observed[[name]]
            )
        }
    }
    if (length(mismatches) > 0L) {
        backend_abort(
            "BACKEND_FAILED",
            "The R/Bioconductor runtime does not match the locked P0 identity.",
            list(reason = "runtime_identity_mismatch", mismatches = mismatches)
        )
    }
    observed
}

read_request <- function(path) {
    if (!file.exists(path) || file.info(path)$isdir) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized backend request is not a regular file.",
            list(reason = "backend_request_missing", path = path)
        )
    }
    tryCatch(
        jsonlite::fromJSON(path, simplifyVector = FALSE),
        error = function(error) backend_abort(
            "BACKEND_FAILED",
            "The normalized backend request is not valid JSON.",
            list(reason = "backend_request_invalid", cause = class(error)[1L])
        )
    )
}

assert_request_envelope <- function(request) {
    required <- c(
        "schema_version", "kind", "execution_scope", "analysis_request",
        "input_evidence", "input", "design", "contrasts"
    )
    allowed <- c(
        required, "validated_input_bundle", "display_export", "gene_sets"
    )
    observed <- names(request)
    if (!is.list(request) || is.null(observed) ||
        length(setdiff(required, observed)) > 0L ||
        length(setdiff(observed, allowed)) > 0L ||
        anyDuplicated(observed) > 0L) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized backend request envelope is incompatible.",
            list(
                reason = "backend_request_schema_invalid",
                missing_fields = sort(setdiff(required, observed)),
                unexpected_fields = sort(setdiff(observed, allowed))
            )
        )
    }
    expected_schema <- if ("gene_sets" %in% observed) "1.1" else "1.0"
    if (!identical(request$schema_version, expected_schema) ||
        !identical(request$kind, "edger_ql_backend_request") ||
        !is.character(request$execution_scope) ||
        length(request$execution_scope) != 1L ||
        !nzchar(request$execution_scope)) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized backend request identity is incompatible.",
            list(reason = "backend_request_identity_invalid")
        )
    }
    ACTIVE_SCHEMA_VERSION <<- expected_schema
    invisible(request)
}

as_character_vector <- function(value, field) {
    result <- unlist(value, recursive = FALSE, use.names = FALSE)
    if (!is.character(result) || length(result) == 0L || anyNA(result)) {
        backend_abort(
            "BACKEND_FAILED",
            paste0("The normalized field '", field, "' is invalid."),
            list(reason = "normalized_request_invalid", field = field)
        )
    }
    result
}

strict_numeric <- function(values, field) {
    values <- as.character(values)
    decimal_pattern <- "^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
    if (any(!grepl(decimal_pattern, values, perl = TRUE))) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "A continuous design variable contains a non-numeric value.",
            list(reason = "continuous_value_invalid", field = field),
            exit_status = 3L
        )
    }
    result <- suppressWarnings(as.numeric(values))
    if (any(!is.finite(result))) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "A continuous design variable contains a non-finite value.",
            list(reason = "continuous_value_nonfinite", field = field),
            exit_status = 3L
        )
    }
    result
}

build_design <- function(request) {
    input <- request$input
    specification <- request$design
    sample_order <- as_character_vector(input$sample_order, "input.sample_order")
    terms <- as_character_vector(specification$terms, "design.terms")
    metadata <- input$metadata_values
    if (!is.list(metadata)) {
        backend_abort(
            "BACKEND_FAILED",
            "Normalized metadata values are absent.",
            list(reason = "normalized_metadata_missing")
        )
    }
    sample_values <- as_character_vector(
        metadata[[input$metadata$sample_id_column]],
        paste0("metadata.", input$metadata$sample_id_column)
    )
    if (!identical(sample_values, sample_order)) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "Normalized metadata is not aligned to the assay sample order.",
            list(
                reason = "metadata_sample_order_mismatch",
                assay_sample_order = sample_order,
                metadata_sample_order = sample_values
            ),
            exit_status = 3L
        )
    }

    design_data <- data.frame(row.names = sample_order, check.names = FALSE)
    factor_levels <- list()
    variable_types <- list()
    for (term in terms) {
        variable <- specification$variables[[term]]
        values <- as_character_vector(metadata[[term]], paste0("metadata.", term))
        if (length(values) != length(sample_order)) {
            backend_abort(
                "DESIGN_RANK_DEFICIENT",
                "A design variable does not contain one value per sample.",
                list(reason = "design_variable_length", term = term),
                exit_status = 3L
            )
        }
        variable_types[[term]] <- variable$type
        if (identical(variable$type, "factor")) {
            levels <- as_character_vector(
                variable$levels, paste0("design.variables.", term, ".levels")
            )
            encoded <- factor(values, levels = levels, ordered = FALSE)
            if (anyNA(encoded)) {
                backend_abort(
                    "DESIGN_RANK_DEFICIENT",
                    "A factor contains a value outside its declared levels.",
                    list(
                        reason = "factor_level_mismatch",
                        term = term,
                        declared_levels = levels,
                        observed_values = sort(unique(values))
                    ),
                    exit_status = 3L
                )
            }
            design_data[[term]] <- encoded
            factor_levels[[term]] <- levels
        } else if (identical(variable$type, "continuous")) {
            design_data[[term]] <- strict_numeric(values, term)
        } else {
            backend_abort(
                "BACKEND_FAILED",
                "A normalized design variable has an unknown type.",
                list(reason = "normalized_request_invalid", term = term)
            )
        }
    }

    formula <- stats::reformulate(
        termlabels = terms,
        response = NULL,
        intercept = isTRUE(specification$intercept)
    )
    design <- stats::model.matrix(formula, data = design_data)
    rownames(design) <- sample_order
    storage.mode(design) <- "double"
    duplicated_columns <- unique(
        colnames(design)[duplicated(colnames(design))]
    )
    if (length(duplicated_columns) > 0L) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "The statistical design produces ambiguous duplicate coefficient names.",
            list(
                reason = "design_column_names_duplicated",
                duplicated_columns = duplicated_columns,
                design_columns = colnames(design)
            ),
            exit_status = 3L
        )
    }
    qr_design <- qr(design, tol = QR_TOLERANCE, LAPACK = FALSE)
    rank <- qr_design$rank
    column_count <- ncol(design)
    aliases <- if (rank < column_count) limma::nonEstimable(design) else NULL
    if (rank < column_count) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "The statistical design matrix is not full column rank.",
            list(
                reason = "complete_confounding_or_redundant_columns",
                sample_count = nrow(design),
                column_count = column_count,
                rank = rank,
                design_columns = colnames(design),
                non_estimable_coefficients = aliases,
                qr_tolerance = QR_TOLERANCE
            ),
            exit_status = 3L
        )
    }
    residual_df <- nrow(design) - rank
    if (residual_df <= 0L) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "The statistical design has no residual degrees of freedom.",
            list(
                reason = "residual_df_nonpositive",
                sample_count = nrow(design),
                column_count = column_count,
                rank = rank,
                residual_df = residual_df,
                design_columns = colnames(design),
                qr_tolerance = QR_TOLERANCE
            ),
            exit_status = 3L
        )
    }

    list(
        matrix = design,
        terms = terms,
        intercept = isTRUE(specification$intercept),
        variable_types = variable_types,
        factor_levels = factor_levels,
        rank = rank,
        residual_df = residual_df,
        qr_tolerance = QR_TOLERANCE
    )
}

build_contrasts <- function(request, design_info) {
    contrast_specs <- request$contrasts
    design <- design_info$matrix
    if (!is.list(contrast_specs) || length(contrast_specs) == 0L) {
        backend_abort(
            "CONTRAST_NOT_ESTIMABLE",
            "No contrasts were supplied to the backend.",
            list(reason = "contrast_set_empty"),
            exit_status = 3L
        )
    }
    result <- vector("list", length(contrast_specs))
    identifiers <- character(length(contrast_specs))
    svd_design <- svd(design, nu = 0L, nv = ncol(design))
    row_basis <- svd_design$v[, seq_len(design_info$rank), drop = FALSE]
    for (index in seq_along(contrast_specs)) {
        specification <- contrast_specs[[index]]
        identifier <- as.character(specification$contrast_id)
        weights <- unlist(specification$weights, recursive = FALSE, use.names = TRUE)
        if (is.null(names(weights)) || any(!is.finite(as.numeric(weights)))) {
            backend_abort(
                "CONTRAST_NOT_ESTIMABLE",
                "A contrast contains invalid coefficient weights.",
                list(reason = "contrast_weights_invalid", contrast_id = identifier),
                exit_status = 3L
            )
        }
        unknown <- setdiff(names(weights), colnames(design))
        if (length(unknown) > 0L) {
            backend_abort(
                "CONTRAST_NOT_ESTIMABLE",
                "A contrast names coefficients absent from the design matrix.",
                list(
                    reason = "contrast_coefficients_unknown",
                    contrast_id = identifier,
                    unknown_coefficients = unknown,
                    design_columns = colnames(design)
                ),
                exit_status = 3L
            )
        }
        vector <- setNames(numeric(ncol(design)), colnames(design))
        vector[names(weights)] <- as.numeric(weights)
        if (all(vector == 0)) {
            backend_abort(
                "CONTRAST_NOT_ESTIMABLE",
                "A contrast is identically zero.",
                list(reason = "contrast_zero", contrast_id = identifier),
                exit_status = 3L
            )
        }
        projected <- as.vector(row_basis %*% crossprod(row_basis, vector))
        estimability_residual <- max(abs(vector - projected))
        estimability_tolerance <- QR_TOLERANCE * max(1, max(abs(vector)))
        if (!is.finite(estimability_residual) ||
            estimability_residual > estimability_tolerance) {
            backend_abort(
                "CONTRAST_NOT_ESTIMABLE",
                "A contrast is outside the row space of the design matrix.",
                list(
                    reason = "contrast_outside_design_row_space",
                    contrast_id = identifier,
                    estimability_residual = estimability_residual,
                    estimability_tolerance = estimability_tolerance
                ),
                exit_status = 3L
            )
        }
        threshold <- as.numeric(specification$lfc_threshold)
        if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
            backend_abort(
                "CONTRAST_NOT_ESTIMABLE",
                "A contrast has an invalid log2-fold-change threshold.",
                list(reason = "contrast_threshold_invalid", contrast_id = identifier),
                exit_status = 3L
            )
        }
        result[[index]] <- list(
            contrast_id = identifier,
            vector = vector,
            weights = as.list(vector[vector != 0]),
            lfc_threshold = threshold,
            test_method = if (threshold > 0) "glmTreat" else "glmQLFTest",
            estimability_residual = estimability_residual,
            estimability_tolerance = estimability_tolerance
        )
        identifiers[[index]] <- identifier
    }
    if (anyDuplicated(identifiers)) {
        backend_abort(
            "CONTRAST_NOT_ESTIMABLE",
            "Contrast identifiers are duplicated.",
            list(reason = "contrast_ids_duplicated"),
            exit_status = 3L
        )
    }
    result
}

strip_gene_versions <- function(values, enabled) {
    if (isTRUE(enabled)) sub("\\.[0-9]+$", "", values, perl = TRUE) else values
}

strict_counts <- function(values, context) {
    raw <- as.character(values)
    if (any(!grepl("^[0-9]+$", raw, perl = TRUE))) {
        backend_abort(
            "BACKEND_FAILED",
            "A normalized integer-count input contains an invalid value.",
            list(reason = "count_value_invalid", context = context)
        )
    }
    numeric <- suppressWarnings(as.numeric(raw))
    if (any(!is.finite(numeric)) || any(numeric < 0) || any(numeric > 2^53 - 1)) {
        backend_abort(
            "BACKEND_FAILED",
            "A normalized integer-count input exceeds its exact numeric domain.",
            list(reason = "count_value_out_of_domain", context = context)
        )
    }
    numeric
}

assert_gene_ids <- function(gene_ids) {
    if (length(gene_ids) == 0L || anyNA(gene_ids) || any(gene_ids == "") ||
        anyDuplicated(gene_ids)) {
        backend_abort(
            "BACKEND_FAILED",
            "The backend input does not contain unique stable gene identifiers.",
            list(reason = "stable_gene_id_invalid")
        )
    }
}

strict_boolean_scalar <- function(value, field) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
        backend_abort(
            "BACKEND_FAILED",
            paste0("The normalized field '", field, "' must be boolean."),
            list(reason = "normalized_request_invalid", field = field)
        )
    }
    value
}

strict_nonnegative_integer_scalar <- function(value, field) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value != floor(value)) {
        backend_abort(
            "BACKEND_FAILED",
            paste0(
                "The normalized field '", field,
                "' must be a nonnegative integer."
            ),
            list(reason = "normalized_request_invalid", field = field)
        )
    }
    as.numeric(value)
}

strict_character_scalar <- function(value, field) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
        backend_abort(
            "BACKEND_FAILED",
            paste0(
                "The normalized field '", field,
                "' must be a non-empty string."
            ),
            list(reason = "normalized_request_invalid", field = field)
        )
    }
    value
}

validate_inferential_replicate_summary <- function(input, samples) {
    summary <- input$salmon$inferential_replicates
    if (!is.list(summary)) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized inferential-replicate summary is incomplete.",
            list(reason = "normalized_replicate_summary_invalid")
        )
    }
    records <- summary$per_sample
    if (!is.list(records) || length(records) != length(samples) ||
        any(!vapply(records, is.list, logical(1L)))) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized inferential-replicate summary is incomplete.",
            list(reason = "normalized_replicate_summary_invalid")
        )
    }
    record_samples <- vapply(
        seq_along(records),
        function(index) strict_character_scalar(
            records[[index]]$sample_id,
            paste0("input.salmon.inferential_replicates.per_sample[", index,
                   "].sample_id")
        ),
        character(1L)
    )
    if (!identical(record_samples, samples)) {
        backend_abort(
            "BACKEND_FAILED",
            "The inferential-replicate summary is not aligned to sample order.",
            list(
                reason = "normalized_replicate_sample_order_mismatch",
                expected = samples,
                observed = record_samples
            )
        )
    }
    present <- vapply(
        seq_along(records),
        function(index) strict_boolean_scalar(
            records[[index]]$present,
            paste0("input.salmon.inferential_replicates.per_sample[", index,
                   "].present")
        ),
        logical(1L)
    )
    counts <- vapply(
        seq_along(records),
        function(index) strict_nonnegative_integer_scalar(
            records[[index]]$count,
            paste0("input.salmon.inferential_replicates.per_sample[", index,
                   "].count")
        ),
        numeric(1L)
    )
    methods <- vapply(
        seq_along(records),
        function(index) {
            method <- records[[index]]$method
            if (is.null(method)) return(NA_character_)
            strict_character_scalar(
                method,
                paste0("input.salmon.inferential_replicates.per_sample[", index,
                       "].method")
            )
        },
        character(1L)
    )
    state <- strict_character_scalar(
        summary$state, "input.salmon.inferential_replicates.state"
    )
    consistent <- strict_boolean_scalar(
        summary$consistent_method_and_count,
        "input.salmon.inferential_replicates.consistent_method_and_count"
    )
    expected_state <- if (all(present)) {
        "all"
    } else if (any(present)) {
        "mixed"
    } else {
        "none"
    }
    if (!identical(state, expected_state) || !consistent || state == "mixed") {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized inferential-replicate summary is inconsistent.",
            list(
                reason = "normalized_replicate_summary_inconsistent",
                declared_state = state,
                observed_state = expected_state,
                declared_consistent = consistent
            )
        )
    }

    if (state == "none") {
        if (any(counts != 0L) || any(!is.na(methods)) ||
            !is.null(summary$replicate_count) || !is.null(summary$method)) {
            backend_abort(
                "BACKEND_FAILED",
                "A no-replicate summary contains positive replicate evidence.",
                list(reason = "normalized_replicate_summary_inconsistent")
            )
        }
        count <- 0L
        method <- NULL
    } else {
        if (any(counts < 1L) || length(unique(counts)) != 1L ||
            any(is.na(methods)) || length(unique(methods)) != 1L) {
            backend_abort(
                "BACKEND_FAILED",
                "Replicate method and count must agree across every sample.",
                list(reason = "normalized_replicate_summary_inconsistent")
            )
        }
        count <- counts[[1L]]
        method <- methods[[1L]]
        declared_count <- strict_nonnegative_integer_scalar(
            summary$replicate_count,
            "input.salmon.inferential_replicates.replicate_count"
        )
        declared_method <- strict_character_scalar(
            summary$method, "input.salmon.inferential_replicates.method"
        )
        if (declared_count != count || !identical(declared_method, method)) {
            backend_abort(
                "BACKEND_FAILED",
                "The aggregate and per-sample replicate summaries disagree.",
                list(reason = "normalized_replicate_summary_inconsistent")
            )
        }
    }

    if (identical(input$input_semantics, "salmon_quant_dirs_full_length") &&
        state == "all" && count < 2L) {
        backend_abort(
            "BACKEND_FAILED",
            paste(
                "Full-length Salmon inferential overdispersion requires at",
                "least two replicates per sample."
            ),
            list(
                reason = "inferential_replicate_count_below_minimum",
                observed_replicates_per_sample = count,
                minimum_replicates_per_sample = 2L
            )
        )
    }
    list(state = state, count = count, method = method, present = present)
}

validate_salmon_route <- function(input, samples) {
    replicate_info <- validate_inferential_replicate_summary(input, samples)
    route <- input$route
    if (!is.list(route)) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized Salmon route is invalid.",
            list(reason = "normalized_route_invalid")
        )
    }
    tximport_options <- route$tximport
    if (!is.list(tximport_options) ||
        !identical(tximport_options$countsFromAbundance, "no") ||
        strict_boolean_scalar(
            tximport_options$dropInfReps, "input.route.tximport.dropInfReps"
        )) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized Salmon tximport route is invalid.",
            list(reason = "normalized_route_invalid")
        )
    }

    semantics <- input$input_semantics
    if (identical(semantics, "salmon_quant_dirs_full_length")) {
        expected_divide <- identical(replicate_info$state, "all")
        edgeR_options <- route$edgeR_options
        overdispersion <- route$inferential_overdispersion
        if (!is.list(edgeR_options) || !is.list(overdispersion)) {
            backend_abort(
                "BACKEND_FAILED",
                "The normalized full-length Salmon route is incomplete.",
                list(reason = "normalized_route_invalid")
            )
        }
        divide <- strict_boolean_scalar(
            edgeR_options$divide, "input.route.edgeR_options.divide"
        )
        enabled <- strict_boolean_scalar(
            overdispersion$enabled,
            "input.route.inferential_overdispersion.enabled"
        )
        relative_adjustment <- strict_boolean_scalar(
            overdispersion$relative_abundance_adjustment,
            paste0(
                "input.route.inferential_overdispersion.",
                "relative_abundance_adjustment"
            )
        )
        if (!identical(route$edgeR_constructor, "edgeR::DGEListFromTximport") ||
            !identical(route$count_source, "txi$counts") ||
            !identical(route$transcript_length_offset, TRUE) ||
            !identical(overdispersion$source, "tximport.infReps") ||
            !identical(divide, expected_divide) ||
            !identical(enabled, divide) ||
            !identical(relative_adjustment, divide)) {
            backend_abort(
                "BACKEND_FAILED",
                "The full-length Salmon route contradicts its input evidence.",
                list(
                    reason = "normalized_route_replicate_mismatch",
                    expected_divide = expected_divide,
                    observed_divide = divide
                )
            )
        }
        return(list(divide = divide, replicates = replicate_info))
    }

    if (identical(semantics, "salmon_quant_dirs_three_prime")) {
        if (!identical(route$edgeR_constructor, "edgeR::DGEList") ||
            !identical(route$count_source, "txi$counts") ||
            !identical(route$transcript_length_offset, FALSE) ||
            !identical(route$gene_length_correction, FALSE) ||
            !identical(route$certified_path_execution_permitted, TRUE) ||
            !is.null(route$edgeR_options) ||
            !is.null(route$inferential_overdispersion)) {
            backend_abort(
                "BACKEND_FAILED",
                "The three-prime Salmon route contradicts its declared semantics.",
                list(reason = "normalized_route_invalid")
            )
        }
        return(list(divide = FALSE, replicates = replicate_info))
    }

    backend_abort(
        "BACKEND_FAILED",
        "The normalized Salmon route is unsupported.",
        list(reason = "normalized_route_invalid", input_semantics = semantics)
    )
}

validate_imported_replicates <- function(txi, route_info, samples) {
    replicate_info <- route_info$replicates
    if (identical(replicate_info$state, "none")) {
        if (!is.null(txi$infReps)) {
            backend_abort(
                "BACKEND_FAILED",
                "tximport imported replicates absent from validated evidence.",
                list(reason = "imported_replicates_evidence_mismatch")
            )
        }
        return(FALSE)
    }
    if (!is.list(txi$infReps) || length(txi$infReps) != length(samples)) {
        backend_abort(
            "BACKEND_FAILED",
            "tximport did not import one replicate matrix per sample.",
            list(reason = "inferential_replicates_missing_after_tximport")
        )
    }
    if (!is.null(names(txi$infReps)) &&
        !identical(names(txi$infReps), samples)) {
        backend_abort(
            "BACKEND_FAILED",
            "tximport inferential replicates are out of sample order.",
            list(reason = "imported_replicate_sample_order_mismatch")
        )
    }
    imported_counts <- vapply(
        txi$infReps,
        function(replicates) {
            if (length(dim(replicates)) != 2L) return(NA_integer_)
            ncol(replicates)
        },
        integer(1L)
    )
    if (anyNA(imported_counts) ||
        any(imported_counts != replicate_info$count)) {
        backend_abort(
            "BACKEND_FAILED",
            "Imported replicate dimensions differ from validated evidence.",
            list(
                reason = "imported_replicate_count_mismatch",
                expected_replicates_per_sample = replicate_info$count,
                observed_replicates_per_sample = imported_counts
            )
        )
    }
    if (isTRUE(route_info$divide) && any(imported_counts < 2L)) {
        backend_abort(
            "BACKEND_FAILED",
            "The divided-count route requires at least two replicates per sample.",
            list(
                reason = "inferential_replicate_count_below_minimum",
                observed_replicates_per_sample = imported_counts,
                minimum_replicates_per_sample = 2L
            )
        )
    }
    TRUE
}

read_combined_featurecounts <- function(input) {
    info <- input$featurecounts
    table <- utils::read.delim(
        info$matrix_path,
        header = TRUE,
        sep = "\t",
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        colClasses = "character"
    )
    samples <- as_character_vector(input$sample_order, "input.sample_order")
    displays <- unlist(info$display_columns, recursive = FALSE, use.names = FALSE)
    expected <- c(info$gene_id_column, displays, samples)
    if (!identical(colnames(table), expected)) {
        backend_abort(
            "BACKEND_FAILED",
            "The private featureCounts snapshot header changed unexpectedly.",
            list(reason = "snapshot_header_mismatch")
        )
    }
    gene_ids <- strip_gene_versions(
        table[[1L]], input$gene_id_policy$strip_version
    )
    assert_gene_ids(gene_ids)
    count_columns <- table[samples]
    counts <- vapply(
        samples,
        function(sample) strict_counts(count_columns[[sample]], sample),
        numeric(nrow(table))
    )
    if (length(samples) == 1L) counts <- matrix(counts, ncol = 1L)
    rownames(counts) <- gene_ids
    colnames(counts) <- samples
    edgeR::DGEList(
        counts = counts,
        genes = data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
    )
}

read_per_sample_featurecounts <- function(input) {
    info <- input$featurecounts
    samples <- as_character_vector(input$sample_order, "input.sample_order")
    matrices <- vector("list", length(samples))
    gene_ids <- NULL
    for (index in seq_along(samples)) {
        sample <- samples[[index]]
        record <- info$files[[index]]
        if (!identical(record$sample_id, sample)) {
            backend_abort(
                "BACKEND_FAILED",
                "Per-sample featureCounts snapshots are out of order.",
                list(reason = "snapshot_sample_order_mismatch", sample_id = sample)
            )
        }
        table <- utils::read.delim(
            record$path,
            header = TRUE,
            sep = "\t",
            quote = "\"",
            comment.char = "#",
            check.names = FALSE,
            stringsAsFactors = FALSE,
            colClasses = "character"
        )
        expected <- c(
            "Geneid", "Chr", "Start", "End", "Strand", "Length",
            record$count_column
        )
        if (!identical(colnames(table), expected)) {
            backend_abort(
                "BACKEND_FAILED",
                "A private per-sample featureCounts snapshot has an invalid header.",
                list(reason = "snapshot_header_mismatch", sample_id = sample)
            )
        }
        observed_ids <- strip_gene_versions(
            table$Geneid, input$gene_id_policy$strip_version
        )
        assert_gene_ids(observed_ids)
        if (is.null(gene_ids)) {
            gene_ids <- observed_ids
        } else if (!identical(observed_ids, gene_ids)) {
            backend_abort(
                "BACKEND_FAILED",
                "Per-sample snapshots contain different gene inventories.",
                list(reason = "snapshot_gene_order_mismatch", sample_id = sample)
            )
        }
        matrices[[index]] <- strict_counts(table[[record$count_column]], sample)
    }
    counts <- do.call(cbind, matrices)
    rownames(counts) <- gene_ids
    colnames(counts) <- samples
    edgeR::DGEList(
        counts = counts,
        genes = data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
    )
}

read_benchmark_counts <- function(input) {
    table <- utils::read.delim(
        input$benchmark_counts$matrix_path,
        header = TRUE,
        sep = "\t",
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        colClasses = "character"
    )
    samples <- as_character_vector(input$sample_order, "input.sample_order")
    if (!identical(colnames(table), c("gene_id", samples))) {
        backend_abort(
            "BACKEND_FAILED",
            "The benchmark count snapshot has an invalid header.",
            list(reason = "benchmark_count_header_mismatch")
        )
    }
    gene_ids <- table$gene_id
    assert_gene_ids(gene_ids)
    counts <- vapply(
        samples,
        function(sample) strict_counts(table[[sample]], sample),
        numeric(nrow(table))
    )
    if (length(samples) == 1L) counts <- matrix(counts, ncol = 1L)
    rownames(counts) <- gene_ids
    colnames(counts) <- samples
    edgeR::DGEList(
        counts = counts,
        genes = data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
    )
}

read_salmon <- function(input) {
    samples <- as_character_vector(input$sample_order, "input.sample_order")
    route_info <- validate_salmon_route(input, samples)
    records <- input$salmon$samples
    quant_dirs <- vapply(records, function(record) record$quant_dir, character(1L))
    record_samples <- vapply(records, function(record) record$sample_id, character(1L))
    if (!identical(record_samples, samples)) {
        backend_abort(
            "BACKEND_FAILED",
            "Salmon snapshots are not aligned to the analysis sample order.",
            list(reason = "snapshot_sample_order_mismatch")
        )
    }
    files <- stats::setNames(file.path(quant_dirs, "quant.sf"), samples)
    mapping <- utils::read.delim(
        input$salmon$tx2gene_path,
        header = TRUE,
        sep = "\t",
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE,
        colClasses = "character"
    )
    if (!identical(colnames(mapping)[1:2], c("transcript_id", "gene_id"))) {
        backend_abort(
            "BACKEND_FAILED",
            "The private tx2gene snapshot has an invalid header.",
            list(reason = "snapshot_tx2gene_header_mismatch")
        )
    }
    mapping$gene_id <- strip_gene_versions(
        mapping$gene_id, input$gene_id_policy$strip_version
    )
    tx2gene <- mapping[, c("transcript_id", "gene_id"), drop = FALSE]
    txi <- tximport::tximport(
        files,
        type = "salmon",
        tx2gene = tx2gene,
        countsFromAbundance = "no",
        dropInfReps = FALSE
    )
    if (!identical(txi$countsFromAbundance, "no")) {
        backend_abort(
            "BACKEND_FAILED",
            "tximport did not preserve countsFromAbundance='no'.",
            list(reason = "tximport_count_semantics_mismatch")
        )
    }
    inferential_replicates_imported <- validate_imported_replicates(
        txi, route_info, samples
    )
    semantics <- input$input_semantics
    if (identical(semantics, "salmon_quant_dirs_full_length")) {
        divide <- route_info$divide
        y <- edgeR::DGEListFromTximport(txi, divide = divide)
        if (!identical(isTRUE(y$divided.counts), divide)) {
            backend_abort(
                "BACKEND_FAILED",
                "DGEListFromTximport did not honor the validated divide setting.",
                list(
                    reason = "dgelist_divide_mismatch",
                    requested = divide,
                    observed = isTRUE(y$divided.counts)
                )
            )
        }
        if (divide) {
            overdispersion <- y$genes$Overdispersion
            if (is.null(overdispersion) ||
                length(overdispersion) != nrow(y$counts) ||
                any(!is.finite(overdispersion)) || any(overdispersion <= 0)) {
                backend_abort(
                    "BACKEND_FAILED",
                    paste(
                        "DGEListFromTximport did not produce finite positive",
                        "inferential overdispersion estimates."
                    ),
                    list(reason = "inferential_overdispersion_invalid")
                )
            }
        }
        if (is.null(y$offset.prior)) {
            backend_abort(
                "BACKEND_FAILED",
                "DGEListFromTximport did not construct the required length offset.",
                list(reason = "salmon_length_offset_missing")
            )
        }
        route_observed <- list(
            constructor = "edgeR::DGEListFromTximport",
            transcript_length_offset = TRUE,
            countsFromAbundance = "no",
            dropInfReps = FALSE,
            divide = divide,
            inferential_replicates_imported = inferential_replicates_imported
        )
    } else if (identical(semantics, "salmon_quant_dirs_three_prime")) {
        gene_ids <- rownames(txi$counts)
        assert_gene_ids(gene_ids)
        y <- edgeR::DGEList(
            counts = txi$counts,
            genes = data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
        )
        if (!is.null(y$offset) || !is.null(y$offset.prior)) {
            backend_abort(
                "BACKEND_FAILED",
                "The three-prime route unexpectedly contains a length offset.",
                list(reason = "three_prime_offset_present")
            )
        }
        route_observed <- list(
            constructor = "edgeR::DGEList",
            count_source = "txi$counts",
            transcript_length_offset = FALSE,
            gene_length_correction = FALSE,
            countsFromAbundance = "no",
            dropInfReps = FALSE,
            divide = FALSE,
            inferential_replicates_imported = inferential_replicates_imported
        )
    } else {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized Salmon route is unsupported.",
            list(reason = "normalized_route_invalid", input_semantics = semantics)
        )
    }
    assert_gene_ids(rownames(y$counts))
    if (is.null(y$genes$gene_id)) y$genes$gene_id <- rownames(y$counts)
    list(dge = y, route_observed = route_observed)
}

construct_dge <- function(input) {
    semantics <- input$input_semantics
    if (identical(semantics, "featurecounts_integer")) {
        if (identical(input$featurecounts$layout, "combined_matrix")) {
            y <- read_combined_featurecounts(input)
        } else if (identical(input$featurecounts$layout, "per_sample_files")) {
            y <- read_per_sample_featurecounts(input)
        } else {
            backend_abort(
                "BACKEND_FAILED",
                "The normalized featureCounts layout is unsupported.",
                list(reason = "normalized_route_invalid")
            )
        }
        return(list(
            dge = y,
            route_observed = list(
                constructor = "edgeR::DGEList",
                count_semantics = "integer",
                transcript_length_offset = FALSE
            )
        ))
    }
    if (semantics %in% c(
        "salmon_quant_dirs_full_length", "salmon_quant_dirs_three_prime"
    )) return(read_salmon(input))
    if (identical(semantics, "benchmark_kernel_integer_counts")) {
        return(list(
            dge = read_benchmark_counts(input),
            route_observed = list(
                constructor = "edgeR::DGEList",
                benchmark_scope = "backend_kernel_only",
                transcript_length_offset = FALSE
            )
        ))
    }
    backend_abort(
        "BACKEND_FAILED",
        "The normalized backend input route is unsupported.",
        list(reason = "normalized_route_invalid", input_semantics = semantics)
    )
}

long_design_table <- function(design) {
    rows <- vector("list", nrow(design) * ncol(design))
    cursor <- 1L
    for (sample_index in seq_len(nrow(design))) {
        for (column_index in seq_len(ncol(design))) {
            rows[[cursor]] <- data.frame(
                sample_id = rownames(design)[sample_index],
                coefficient = colnames(design)[column_index],
                value = design[sample_index, column_index],
                stringsAsFactors = FALSE
            )
            cursor <- cursor + 1L
        }
    }
    do.call(rbind, rows)
}

long_coefficient_table <- function(gene_ids, keep, fit, design_columns) {
    expected_gene_ids <- gene_ids[keep]
    if (!identical(rownames(fit$coefficients), expected_gene_ids)) {
        backend_abort(
            "BACKEND_FAILED",
            "edgeR fitted coefficients are not aligned to filtered gene IDs.",
            list(reason = "fit_gene_order_mismatch")
        )
    }
    fitted <- matrix(
        NA_real_,
        nrow = length(gene_ids),
        ncol = length(design_columns),
        dimnames = list(gene_ids, design_columns)
    )
    fitted[keep, ] <- fit$coefficients[, design_columns, drop = FALSE]
    if (any(!is.finite(fitted[keep, , drop = FALSE]))) {
        backend_abort(
            "BACKEND_FAILED",
            "edgeR returned a non-finite fitted coefficient.",
            list(reason = "nonfinite_fitted_coefficient")
        )
    }
    rows <- vector("list", length(gene_ids) * length(design_columns))
    cursor <- 1L
    for (gene_index in seq_along(gene_ids)) {
        for (column_index in seq_along(design_columns)) {
            rows[[cursor]] <- data.frame(
                gene_id = gene_ids[[gene_index]],
                status = if (keep[[gene_index]]) "tested" else "filtered",
                coefficient = design_columns[[column_index]],
                estimate = fitted[gene_index, column_index],
                scale = "natural_log",
                stringsAsFactors = FALSE
            )
            cursor <- cursor + 1L
        }
    }
    do.call(rbind, rows)
}

contrast_results <- function(gene_ids, keep, fit, contrast) {
    test <- if (contrast$lfc_threshold > 0) {
        edgeR::glmTreat(
            fit,
            contrast = contrast$vector,
            lfc = contrast$lfc_threshold,
            null = "interval"
        )
    } else {
        edgeR::glmQLFTest(fit, contrast = contrast$vector)
    }
    table <- test$table
    expected_gene_ids <- gene_ids[keep]
    if (!identical(rownames(table), expected_gene_ids)) {
        backend_abort(
            "BACKEND_FAILED",
            "The edgeR test table is not aligned to filtered gene IDs.",
            list(
                reason = "test_gene_order_mismatch",
                contrast_id = contrast$contrast_id
            )
        )
    }
    is_treat <- identical(contrast$test_method, "glmTreat")
    if (!is_treat && !"F" %in% colnames(table)) {
        backend_abort(
            "BACKEND_FAILED",
            "The edgeR test result lacks its quasi-likelihood F statistic.",
            list(
                reason = "unexpected_edger_result_schema",
                contrast_id = contrast$contrast_id,
                columns = colnames(table)
            )
        )
    }
    fdr <- stats::p.adjust(table$PValue, method = "BH")
    required_numeric <- c("logFC", "logCPM", "PValue")
    if (!is_treat) required_numeric <- c(required_numeric, "F")
    nonfinite_columns <- required_numeric[vapply(
        required_numeric,
        function(column) any(!is.finite(table[[column]])),
        logical(1L)
    )]
    if (any(!is.finite(fdr))) nonfinite_columns <- c(nonfinite_columns, "FDR")
    if ("unshrunk.logFC" %in% colnames(table) &&
        any(!is.finite(table$unshrunk.logFC))) {
        nonfinite_columns <- c(nonfinite_columns, "unshrunk.logFC")
    }
    if (length(nonfinite_columns) > 0L) {
        backend_abort(
            "BACKEND_FAILED",
            "edgeR returned a non-finite value for a tested gene.",
            list(
                reason = "nonfinite_test_result",
                contrast_id = contrast$contrast_id,
                columns = unique(nonfinite_columns)
            )
        )
    }
    output <- data.frame(
        gene_id = gene_ids,
        contrast_id = contrast$contrast_id,
        status = ifelse(keep, "tested", "filtered"),
        logFC = NA_real_,
        unshrunk_logFC = NA_real_,
        logCPM = NA_real_,
        statistic = NA_real_,
        statistic_type = if (is_treat) "not_reported_by_glmTreat" else "F",
        statistic_status = ifelse(
            keep,
            if (is_treat) "not_reported" else "reported",
            "not_applicable_filtered"
        ),
        PValue = NA_real_,
        FDR = NA_real_,
        test_method = contrast$test_method,
        lfc_threshold = contrast$lfc_threshold,
        stringsAsFactors = FALSE
    )
    output$logFC[keep] <- table$logFC
    if ("unshrunk.logFC" %in% colnames(table)) {
        output$unshrunk_logFC[keep] <- table$unshrunk.logFC
    }
    output$logCPM[keep] <- table$logCPM
    if (!is_treat) output$statistic[keep] <- table$F
    output$PValue[keep] <- table$PValue
    output$FDR[keep] <- fdr
    output
}

write_tsv <- function(document, path) {
    utils::write.table(
        document,
        file = path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = "",
        eol = "\n",
        fileEncoding = "UTF-8"
    )
}

file_record <- function(path, relative_path) {
    list(
        path = relative_path,
        sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE),
        size_bytes = unname(file.info(path)$size)
    )
}

write_display_logcpm <- function(request, y) {
    specification <- request$display_export
    if (is.null(specification)) return(NULL)
    if (!is.list(specification) ||
        !identical(sort(names(specification)), c("path"))) {
        backend_abort(
            "BACKEND_FAILED",
            "The private display export instruction is invalid.",
            list(reason = "display_export_invalid")
        )
    }
    path <- specification$path
    if (!is.character(path) || length(path) != 1L || !nzchar(path) ||
        file.exists(path) || !dir.exists(dirname(path))) {
        backend_abort(
            "BACKEND_FAILED",
            "The private display export target is unavailable.",
            list(reason = "display_export_target_invalid")
        )
    }
    values <- edgeR::cpm(
        y,
        normalized.lib.sizes = TRUE,
        log = TRUE,
        prior.count = 2
    )
    if (!is.matrix(values) || any(!is.finite(values)) ||
        !identical(colnames(values), colnames(y$counts)) ||
        !identical(rownames(values), rownames(y$counts))) {
        backend_abort(
            "BACKEND_FAILED",
            "The display-only logCPM matrix is invalid.",
            list(reason = "display_logcpm_invalid")
        )
    }
    table <- data.frame(
        gene_id = as.character(y$genes$gene_id),
        as.data.frame(values, check.names = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    write_tsv(table, path)
    list(
        method = "edgeR::cpm.DGEList",
        source = "post_filter_post_TMM_observed_DGEList",
        purpose = "display_only_not_for_inference",
        arguments = list(
            normalized.lib.sizes = TRUE,
            log = TRUE,
            prior.count = 2
        ),
        scale = "log2",
        gene_count = nrow(values),
        sample_count = ncol(values)
    )
}

run_analysis <- function(request_path, output_dir) {
    if (file.exists(output_dir)) {
        backend_abort(
            "BACKEND_FAILED",
            "The private backend output stage already exists.",
            list(reason = "output_stage_exists", output_dir = output_dir)
        )
    }
    runtime <- assert_runtime()
    request <- read_request(request_path)
    assert_request_envelope(request)
    output_schema_version <- if (is.null(request$gene_sets)) "1.0" else "1.1"
    design_info <- build_design(request)
    contrasts <- build_contrasts(request, design_info)
    constructed <- construct_dge(request$input)
    y_all <- constructed$dge
    design <- design_info$matrix
    if (!identical(colnames(y_all$counts), rownames(design))) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "The design rows and count columns are not identically aligned.",
            list(
                reason = "design_assay_order_mismatch",
                count_samples = colnames(y_all$counts),
                design_samples = rownames(design)
            ),
            exit_status = 3L
        )
    }

    gene_ids <- as.character(y_all$genes$gene_id)
    assert_gene_ids(gene_ids)
    keep <- edgeR::filterByExpr(y_all, design = design)
    if (!any(keep)) {
        backend_abort(
            "BACKEND_FAILED",
            "filterByExpr removed every gene; the QL model cannot be fitted.",
            list(reason = "all_genes_filtered", gene_count = length(keep))
        )
    }
    y <- y_all[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::normLibSizes(y, method = "TMM")
    fit <- edgeR::glmQLFit(
        y,
        design = design,
        abundance.trend = TRUE,
        robust = TRUE,
        winsor.tail.p = c(0.05, 0.1),
        legacy = FALSE,
        top.proportion = NULL,
        keep.unit.mat = FALSE
    )
    display_export <- write_display_logcpm(request, y)

    result_tables <- lapply(
        contrasts,
        function(contrast) contrast_results(gene_ids, keep, fit, contrast)
    )
    results <- do.call(rbind, result_tables)
    coefficients <- long_coefficient_table(
        gene_ids, keep, fit, colnames(design)
    )
    design_table <- long_design_table(design)
    contrast_provenance <- lapply(contrasts, function(contrast) list(
        contrast_id = contrast$contrast_id,
        weights = contrast$weights,
        lfc_threshold = contrast$lfc_threshold,
        test_method = contrast$test_method,
        treat_null = if (contrast$lfc_threshold > 0) "interval" else NULL,
        estimability_residual = contrast$estimability_residual,
        estimability_tolerance = contrast$estimability_tolerance
    ))
    pathway <- if (is.null(request$gene_sets)) {
        NULL
    } else {
        run_pathway_tests(
            request$gene_sets,
            request$input,
            gene_ids,
            keep,
            fit,
            design,
            contrasts
        )
    }
    analysis <- list(
        schema_version = output_schema_version,
        kind = "edger_ql_analysis",
        backend = "edgeR_QL",
        execution_scope = request$execution_scope,
        analysis_request = request$analysis_request,
        input_evidence = request$input_evidence,
        runtime_identity = runtime,
        input_semantics = request$input$input_semantics,
        route_observed = constructed$route_observed,
        pipeline = list(
            list(step = "filterByExpr", arguments = list(design = TRUE)),
            list(step = "normLibSizes", arguments = list(method = "TMM")),
            list(
                step = "glmQLFit",
                arguments = list(legacy = FALSE, robust = TRUE)
            ),
            list(
                step = "contrast_test",
                dispatch = "lfc_threshold == 0: glmQLFTest; lfc_threshold > 0: glmTreat"
            )
        ),
        design = list(
            intercept = design_info$intercept,
            terms = design_info$terms,
            variable_types = design_info$variable_types,
            factor_levels = design_info$factor_levels,
            columns = colnames(design),
            sample_count = nrow(design),
            rank = design_info$rank,
            residual_df = design_info$residual_df,
            qr_tolerance = design_info$qr_tolerance
        ),
        contrasts = contrast_provenance,
        genes = list(
            total = length(gene_ids),
            tested = sum(keep),
            filtered = sum(!keep)
        ),
        status_vocabulary = STATUS_VOCABULARY,
        result_logFC_scale = "log2",
        coefficient_scale = "natural_log",
        multiple_testing = "Benjamini-Hochberg within each contrast"
    )
    if (!is.null(pathway)) {
        analysis$pathway_analysis <- pathway$provenance
    }

    if (!dir.create(output_dir, recursive = FALSE, mode = "0700")) {
        backend_abort(
            "BACKEND_FAILED",
            "The private backend output stage could not be created.",
            list(reason = "output_stage_create_failed", output_dir = output_dir)
        )
    }
    results_path <- file.path(output_dir, "results.tsv")
    coefficients_path <- file.path(output_dir, "coefficients.tsv")
    design_path <- file.path(output_dir, "design.tsv")
    analysis_path <- file.path(output_dir, "analysis.json")
    write_tsv(results, results_path)
    write_tsv(coefficients, coefficients_path)
    write_tsv(design_table, design_path)
    if (!is.null(pathway)) {
        write_tsv(
            pathway$results,
            file.path(output_dir, "pathway_results.tsv")
        )
    }
    writeLines(
        jsonlite::toJSON(
            analysis,
            auto_unbox = TRUE,
            null = "null",
            na = "null",
            digits = NA,
            pretty = TRUE
        ),
        con = analysis_path,
        useBytes = TRUE
    )
    member_paths <- c(
        "analysis.json", "coefficients.tsv", "design.tsv", "results.tsv"
    )
    if (!is.null(pathway)) {
        member_paths <- c(member_paths, "pathway_results.tsv")
    }
    members <- lapply(
        member_paths,
        function(name) file_record(file.path(output_dir, name), name)
    )
    manifest <- list(
        schema_version = output_schema_version,
        kind = "edger_ql_backend_manifest",
        backend = "edgeR_QL",
        runtime_identity = runtime,
        execution_scope = request$execution_scope,
        input_evidence = request$input_evidence,
        members = members
    )
    manifest_path <- file.path(output_dir, "backend_manifest.json")
    writeLines(
        jsonlite::toJSON(
            manifest,
            auto_unbox = TRUE,
            null = "null",
            na = "null",
            digits = NA,
            pretty = TRUE
        ),
        con = manifest_path,
        useBytes = TRUE
    )

    response_data <- list(
        execution_scope = request$execution_scope,
        runtime_identity = runtime,
        input_semantics = request$input$input_semantics,
        route_observed = constructed$route_observed,
        design_columns = colnames(design),
        design_rank = design_info$rank,
        residual_df = design_info$residual_df,
        gene_count = length(gene_ids),
        tested_gene_count = sum(keep),
        filtered_gene_count = sum(!keep),
        contrasts = contrast_provenance,
        ql_fit_parameters = list(
            abundance.trend = TRUE,
            robust = TRUE,
            winsor.tail.p = c(0.05, 0.1),
            legacy = FALSE,
            top.proportion = NULL,
            resolved_top.proportion = fit$top.proportion,
            keep.unit.mat = FALSE
        ),
        display_export = display_export
    )
    if (!is.null(pathway)) {
        response_data$pathway_analysis <- list(
            enabled = TRUE,
            result_row_count = nrow(pathway$results),
            gene_set_count = length(unique(pathway$results$gene_set_id)),
            methods = c("limma_mroast", "limma_fry", "limma_camera")
        )
    }
    emit_document(list(
        schema_version = output_schema_version,
        status = "success",
        backend = "edgeR_QL",
        data = response_data,
        warnings = list(),
        errors = list(),
        artifacts = c(members, list(file_record(manifest_path, "backend_manifest.json")))
    ))
}

main <- function() {
    file_argument <- grep(
        "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
    )
    if (length(file_argument) != 1L) {
        backend_abort(
            "BACKEND_FAILED",
            "The locked backend script location cannot be resolved.",
            list(reason = "backend_script_path_unavailable")
        )
    }
    script_path <- normalizePath(
        sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE
    )
    source(
        file.path(dirname(script_path), "pathway_tests.R"),
        local = .GlobalEnv,
        encoding = "UTF-8"
    )
    arguments <- commandArgs(trailingOnly = TRUE)
    if (length(arguments) != 2L) {
        backend_abort(
            "BACKEND_FAILED",
            "Usage: edger_ql.R NORMALIZED_REQUEST_JSON PRIVATE_OUTPUT_STAGE",
            list(reason = "backend_arguments_invalid")
        )
    }
    run_analysis(arguments[[1L]], arguments[[2L]])
}

exit_status <- 0L
tryCatch(
    main(),
    rnaseq_backend_error = function(error) {
        emit_document(list(
            schema_version = ACTIVE_SCHEMA_VERSION,
            status = "error",
            backend = "edgeR_QL",
            data = NULL,
            warnings = list(),
            errors = list(list(
                code = error$code,
                message = error$message,
                details = error$details
            )),
            artifacts = list()
        ))
        exit_status <<- error$exit_status
    },
    error = function(error) {
        emit_document(list(
            schema_version = ACTIVE_SCHEMA_VERSION,
            status = "error",
            backend = "edgeR_QL",
            data = NULL,
            warnings = list(),
            errors = list(list(
                code = "BACKEND_FAILED",
                message = "The locked edgeR backend failed unexpectedly.",
                details = list(
                    reason = "unexpected_r_error",
                    cause_type = class(error)[1L],
                    cause_message = conditionMessage(error)
                )
            )),
            artifacts = list()
        ))
        exit_status <<- 4L
    }
)
quit(save = "no", status = exit_status, runLast = FALSE)
