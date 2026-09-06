#!/usr/bin/env Rscript

# Locked DESeq2 differential-expression backend. Python owns public contracts,
# evidence capture, orchestration, publication, and independent verification;
# this process owns every statistical operation and count conversion.

EXPECTED_RUNTIME <- list(
    R = "4.6.1",
    Bioconductor = "3.23",
    BiocVersion_package = "3.23.1",
    DESeq2 = "1.52.0",
    apeglm = "1.34.0",
    tximport = "1.40.0"
)
SCHEMA_VERSION <- "1.0"
BACKEND_NAME <- "DESeq2"
EXECUTION_SCOPE <- "validated_p1_deseq2_input"
QR_TOLERANCE <- 1e-7
INTEGER_MAX <- .Machine$integer.max
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
        Bioconductor = paste(
            strsplit(bioc_version, "\\.")[[1L]][1:2], collapse = "."
        ),
        BiocVersion_package = bioc_version,
        DESeq2 = as.character(utils::packageVersion("DESeq2")),
        apeglm = as.character(utils::packageVersion("apeglm")),
        tximport = as.character(utils::packageVersion("tximport"))
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
            "The R/Bioconductor runtime does not match the locked DESeq2 identity.",
            list(reason = "runtime_identity_mismatch", mismatches = mismatches)
        )
    }
    observed
}

read_request <- function(path) {
    if (!file.exists(path) || file.info(path)$isdir) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized DESeq2 request is not a regular file.",
            list(reason = "backend_request_missing", path = path)
        )
    }
    tryCatch(
        jsonlite::fromJSON(path, simplifyVector = FALSE),
        error = function(error) backend_abort(
            "BACKEND_FAILED",
            "The normalized DESeq2 request is not valid JSON.",
            list(reason = "backend_request_invalid", cause = class(error)[1L])
        )
    )
}

assert_exact_names <- function(value, required, allowed = required, context) {
    observed <- names(value)
    if (!is.list(value) || is.null(observed) || anyDuplicated(observed) > 0L ||
        length(setdiff(required, observed)) > 0L ||
        length(setdiff(observed, allowed)) > 0L) {
        backend_abort(
            "BACKEND_FAILED",
            paste0("The normalized ", context, " is incompatible."),
            list(
                reason = "backend_request_schema_invalid",
                context = context,
                missing_fields = sort(setdiff(required, observed)),
                unexpected_fields = sort(setdiff(observed, allowed))
            )
        )
    }
    invisible(value)
}

assert_request_envelope <- function(request) {
    required <- c(
        "schema_version", "kind", "execution_scope", "analysis_request",
        "input_evidence", "input", "design", "contrasts", "deseq2"
    )
    allowed <- c(required, "validated_input_bundle")
    assert_exact_names(request, required, allowed, "backend request envelope")
    if (!identical(request$schema_version, SCHEMA_VERSION) ||
        !identical(request$kind, "deseq2_backend_request") ||
        !identical(request$execution_scope, EXECUTION_SCOPE)) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized DESeq2 request identity is incompatible.",
            list(reason = "backend_request_identity_invalid")
        )
    }
    invisible(request)
}

as_character_vector <- function(value, field, allow_empty = FALSE) {
    result <- unlist(value, recursive = FALSE, use.names = FALSE)
    if (allow_empty && length(result) == 0L) return(character(0L))
    if (!is.character(result) || anyNA(result) ||
        (!allow_empty && length(result) == 0L) || any(!nzchar(result))) {
        backend_abort(
            "BACKEND_FAILED",
            paste0("The normalized field '", field, "' is invalid."),
            list(reason = "normalized_request_invalid", field = field)
        )
    }
    result
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
            paste0("The normalized field '", field, "' must be a nonnegative integer."),
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
            paste0("The normalized field '", field, "' must be a non-empty string."),
            list(reason = "normalized_request_invalid", field = field)
        )
    }
    value
}

strict_numeric <- function(values, field) {
    values <- as.character(values)
    pattern <- "^[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
    if (any(!grepl(pattern, values, perl = TRUE))) {
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
            "Normalized metadata is not aligned to assay sample order.",
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
    matrix <- stats::model.matrix(formula, data = design_data)
    rownames(matrix) <- sample_order
    storage.mode(matrix) <- "double"
    if (anyDuplicated(colnames(matrix))) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "The design produces duplicate coefficient names.",
            list(
                reason = "design_column_names_duplicated",
                design_columns = colnames(matrix)
            ),
            exit_status = 3L
        )
    }
    qr_matrix <- qr(matrix, tol = QR_TOLERANCE, LAPACK = FALSE)
    rank <- qr_matrix$rank
    if (rank < ncol(matrix)) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "The statistical design matrix is not full column rank.",
            list(
                reason = "complete_confounding_or_redundant_columns",
                sample_count = nrow(matrix),
                column_count = ncol(matrix),
                rank = rank,
                design_columns = colnames(matrix),
                non_estimable_coefficients = limma::nonEstimable(matrix),
                qr_tolerance = QR_TOLERANCE
            ),
            exit_status = 3L
        )
    }
    residual_df <- nrow(matrix) - rank
    if (residual_df <= 0L) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "The statistical design has no residual degrees of freedom.",
            list(
                reason = "residual_df_nonpositive",
                sample_count = nrow(matrix),
                column_count = ncol(matrix),
                rank = rank,
                residual_df = residual_df,
                design_columns = colnames(matrix),
                qr_tolerance = QR_TOLERANCE
            ),
            exit_status = 3L
        )
    }
    list(
        matrix = matrix,
        data = design_data,
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
    specifications <- request$contrasts
    design <- design_info$matrix
    if (!is.list(specifications) || length(specifications) == 0L) {
        backend_abort(
            "CONTRAST_NOT_ESTIMABLE",
            "No contrasts were supplied to the DESeq2 backend.",
            list(reason = "contrast_set_empty"),
            exit_status = 3L
        )
    }
    result <- vector("list", length(specifications))
    identifiers <- character(length(specifications))
    svd_design <- svd(design, nu = 0L, nv = ncol(design))
    row_basis <- svd_design$v[, seq_len(design_info$rank), drop = FALSE]
    for (index in seq_along(specifications)) {
        specification <- specifications[[index]]
        identifier <- strict_character_scalar(
            specification$contrast_id,
            paste0("contrasts[", index, "].contrast_id")
        )
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
        vector <- stats::setNames(numeric(ncol(design)), colnames(design))
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
        residual <- max(abs(vector - projected))
        tolerance <- QR_TOLERANCE * max(1, max(abs(vector)))
        if (!is.finite(residual) || residual > tolerance) {
            backend_abort(
                "CONTRAST_NOT_ESTIMABLE",
                "A contrast is outside the design row space.",
                list(
                    reason = "contrast_outside_design_row_space",
                    contrast_id = identifier,
                    estimability_residual = residual,
                    estimability_tolerance = tolerance
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
            estimability_residual = residual,
            estimability_tolerance = tolerance
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

build_test_specification <- function(request, design_info, contrasts) {
    specification <- request$deseq2
    if (!is.list(specification)) {
        backend_abort(
            "INVALID_REQUEST",
            "The normalized DESeq2 options are missing.",
            list(reason = "deseq2_options_missing"),
            exit_status = 2L
        )
    }
    mode <- strict_character_scalar(specification$test_mode, "deseq2.test_mode")
    shrinkage <- strict_character_scalar(
        specification$shrinkage, "deseq2.shrinkage"
    )
    if (!mode %in% c("wald", "lrt") || !shrinkage %in% c("none", "apeglm")) {
        backend_abort(
            "INVALID_REQUEST",
            "The normalized DESeq2 test or shrinkage mode is unsupported.",
            list(reason = "deseq2_options_invalid"),
            exit_status = 2L
        )
    }
    if (identical(mode, "wald")) {
        assert_exact_names(
            specification,
            c("test_mode", "shrinkage"),
            c("test_mode", "shrinkage"),
            "Wald DESeq2 options"
        )
        return(list(mode = mode, shrinkage = shrinkage, reduced = NULL))
    }

    assert_exact_names(
        specification,
        c("test_mode", "shrinkage", "reduced"),
        c("test_mode", "shrinkage", "reduced"),
        "LRT DESeq2 options"
    )
    if (length(contrasts) != 1L || contrasts[[1L]]$lfc_threshold != 0) {
        backend_abort(
            "INVALID_REQUEST",
            "DESeq2 LRT requires one reporting contrast with lfc_threshold=0.",
            list(
                reason = "lrt_reporting_contract_invalid",
                contrast_count = length(contrasts),
                lfc_thresholds = lapply(contrasts, `[[`, "lfc_threshold")
            ),
            exit_status = 2L
        )
    }
    reduced <- specification$reduced
    assert_exact_names(
        reduced,
        c("intercept", "terms"),
        c("intercept", "terms"),
        "LRT reduced design"
    )
    reduced_intercept <- strict_boolean_scalar(
        reduced$intercept, "deseq2.reduced.intercept"
    )
    reduced_terms <- as_character_vector(
        reduced$terms, "deseq2.reduced.terms", allow_empty = TRUE
    )
    if (!identical(reduced_intercept, design_info$intercept) ||
        length(setdiff(reduced_terms, design_info$terms)) > 0L ||
        anyDuplicated(reduced_terms) ||
        length(reduced_terms) >= length(design_info$terms)) {
        backend_abort(
            "INVALID_REQUEST",
            "The LRT reduced design is not a proper additive subset of the full design.",
            list(
                reason = "lrt_reduced_not_proper_subset",
                full_intercept = design_info$intercept,
                reduced_intercept = reduced_intercept,
                full_terms = design_info$terms,
                reduced_terms = reduced_terms
            ),
            exit_status = 2L
        )
    }
    formula <- stats::reformulate(
        termlabels = reduced_terms,
        response = NULL,
        intercept = reduced_intercept
    )
    matrix <- stats::model.matrix(formula, data = design_info$data)
    rownames(matrix) <- rownames(design_info$matrix)
    storage.mode(matrix) <- "double"
    qr_reduced <- qr(matrix, tol = QR_TOLERANCE, LAPACK = FALSE)
    rank <- qr_reduced$rank
    if (ncol(matrix) == 0L || rank < ncol(matrix) || rank >= design_info$rank) {
        backend_abort(
            "INVALID_REQUEST",
            "The LRT reduced design is rank deficient or not strictly smaller.",
            list(
                reason = "lrt_reduced_rank_invalid",
                reduced_rank = rank,
                reduced_columns = ncol(matrix),
                full_rank = design_info$rank
            ),
            exit_status = 2L
        )
    }
    fitted <- design_info$matrix %*% qr.coef(qr(design_info$matrix), matrix)
    nesting_residual <- max(abs(matrix - fitted))
    nesting_tolerance <- QR_TOLERANCE * max(1, max(abs(matrix)))
    if (!is.finite(nesting_residual) || nesting_residual > nesting_tolerance) {
        backend_abort(
            "INVALID_REQUEST",
            "The LRT reduced design is not nested in the full design column space.",
            list(
                reason = "lrt_reduced_not_nested",
                nesting_residual = nesting_residual,
                nesting_tolerance = nesting_tolerance
            ),
            exit_status = 2L
        )
    }
    list(
        mode = mode,
        shrinkage = shrinkage,
        reduced = list(
            intercept = reduced_intercept,
            terms = reduced_terms,
            columns = colnames(matrix),
            rank = rank,
            residual_df = nrow(matrix) - rank,
            nesting_residual = nesting_residual,
            nesting_tolerance = nesting_tolerance,
            matrix = matrix
        )
    )
}

strip_gene_versions <- function(values, enabled) {
    if (isTRUE(enabled)) sub("\\.[0-9]+$", "", values, perl = TRUE) else values
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

strict_integer_counts <- function(values, context) {
    raw <- as.character(values)
    if (any(!grepl("^[0-9]+$", raw, perl = TRUE))) {
        backend_abort(
            "COUNT_VALUES_INVALID",
            "An integer-count input contains an invalid value.",
            list(reason = "count_value_invalid", context = context),
            exit_status = 3L
        )
    }
    numeric <- suppressWarnings(as.numeric(raw))
    if (any(!is.finite(numeric)) || any(numeric < 0) ||
        any(numeric > INTEGER_MAX)) {
        backend_abort(
            "COUNT_VALUES_INVALID",
            "A count exceeds the DESeq2 signed-integer domain.",
            list(
                reason = "count_value_out_of_deseq2_integer_domain",
                context = context,
                maximum = INTEGER_MAX
            ),
            exit_status = 3L
        )
    }
    as.integer(numeric)
}

matrix_sha256 <- function(values) {
    payload <- serialize(
        list(
            dimensions = dim(values),
            dimnames = dimnames(values),
            storage_mode = storage.mode(values),
            values = as.vector(values)
        ),
        connection = NULL,
        version = 3
    )
    digest::digest(payload, algo = "sha256", serialize = FALSE)
}

rounding_audit <- function(source, rounded, mode) {
    if (!is.matrix(source) || !is.matrix(rounded) ||
        !identical(dim(source), dim(rounded)) ||
        !identical(dimnames(source), dimnames(rounded))) {
        backend_abort(
            "BACKEND_FAILED",
            "The count-rounding audit matrices are not aligned.",
            list(reason = "rounding_audit_matrix_mismatch")
        )
    }
    source_numeric <- matrix(
        as.numeric(source),
        nrow = nrow(source),
        ncol = ncol(source),
        dimnames = dimnames(source)
    )
    rounded_numeric <- matrix(
        as.numeric(rounded),
        nrow = nrow(rounded),
        ncol = ncol(rounded),
        dimnames = dimnames(rounded)
    )
    delta <- rounded_numeric - source_numeric
    absolute <- abs(delta)
    per_sample <- lapply(seq_len(ncol(source_numeric)), function(index) {
        sample_delta <- delta[, index]
        sample_absolute <- absolute[, index]
        list(
            sample_id = colnames(source_numeric)[[index]],
            cell_count = nrow(source_numeric),
            changed_cell_count = sum(sample_delta != 0),
            max_absolute_delta = if (length(sample_absolute)) {
                max(sample_absolute)
            } else {
                0
            },
            absolute_delta_sum = sum(sample_absolute),
            total_count_before = sum(source_numeric[, index]),
            total_count_after = sum(rounded_numeric[, index]),
            total_count_delta = sum(sample_delta)
        )
    })
    maximum <- if (length(absolute)) max(absolute) else 0
    if (!is.finite(maximum) || maximum > 0.5) {
        backend_abort(
            "BACKEND_FAILED",
            "The observed base::round conversion exceeded its mathematical bound.",
            list(reason = "rounding_delta_invalid", max_absolute_delta = maximum)
        )
    }
    list(
        mode = mode,
        round_function = if (identical(mode, "none_integer_input")) {
            "none"
        } else {
            "base::round"
        },
        rule = if (identical(mode, "none_integer_input")) {
            "not_applicable_exact_integer_input"
        } else {
            "R_base_round_IEC_60559_ties_to_even"
        },
        hash_serialization = "R_serialize_version_3_locked_runtime",
        source_matrix_sha256 = matrix_sha256(source),
        rounded_matrix_sha256 = matrix_sha256(rounded),
        gene_count = nrow(source_numeric),
        sample_count = ncol(source_numeric),
        cell_count = length(source_numeric),
        changed_cell_count = sum(delta != 0),
        max_absolute_delta = maximum,
        absolute_delta_sum = sum(absolute),
        per_sample = per_sample
    )
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
            "Inferential-replicate evidence is out of sample order.",
            list(reason = "normalized_replicate_sample_order_mismatch")
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
    expected_state <- if (all(present)) {
        "all"
    } else if (any(present)) {
        "mixed"
    } else {
        "none"
    }
    consistent <- strict_boolean_scalar(
        summary$consistent_method_and_count,
        "input.salmon.inferential_replicates.consistent_method_and_count"
    )
    if (!identical(state, expected_state) || !consistent || state == "mixed") {
        backend_abort(
            "BACKEND_FAILED",
            "The inferential-replicate summary is inconsistent.",
            list(reason = "normalized_replicate_summary_inconsistent")
        )
    }
    if (identical(state, "none")) {
        if (any(counts != 0) || any(!is.na(methods))) {
            backend_abort(
                "BACKEND_FAILED",
                "A no-replicate summary contains positive replicate evidence.",
                list(reason = "normalized_replicate_summary_inconsistent")
            )
        }
        return(list(state = state, count = 0, method = NULL))
    }
    if (any(counts < 1) || length(unique(counts)) != 1L || anyNA(methods) ||
        length(unique(methods)) != 1L ||
        strict_nonnegative_integer_scalar(
            summary$replicate_count,
            "input.salmon.inferential_replicates.replicate_count"
        ) != counts[[1L]] ||
        !identical(
            strict_character_scalar(
                summary$method,
                "input.salmon.inferential_replicates.method"
            ),
            methods[[1L]]
        )) {
        backend_abort(
            "BACKEND_FAILED",
            "Replicate method and count do not agree across samples.",
            list(reason = "normalized_replicate_summary_inconsistent")
        )
    }
    list(state = state, count = counts[[1L]], method = methods[[1L]])
}

validate_imported_replicates <- function(txi, evidence, samples) {
    if (identical(evidence$state, "none")) {
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
    if (!is.null(names(txi$infReps)) && !identical(names(txi$infReps), samples)) {
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
    if (anyNA(imported_counts) || any(imported_counts != evidence$count)) {
        backend_abort(
            "BACKEND_FAILED",
            "Imported replicate dimensions differ from validated evidence.",
            list(
                reason = "imported_replicate_count_mismatch",
                expected = evidence$count,
                observed = imported_counts
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
    gene_ids <- strip_gene_versions(table[[1L]], input$gene_id_policy$strip_version)
    assert_gene_ids(gene_ids)
    columns <- lapply(
        samples,
        function(sample) strict_integer_counts(table[[sample]], sample)
    )
    counts <- do.call(cbind, columns)
    rownames(counts) <- gene_ids
    colnames(counts) <- samples
    counts
}

read_per_sample_featurecounts <- function(input) {
    info <- input$featurecounts
    samples <- as_character_vector(input$sample_order, "input.sample_order")
    columns <- vector("list", length(samples))
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
                "Per-sample featureCounts snapshots have different gene inventories.",
                list(reason = "snapshot_gene_order_mismatch", sample_id = sample)
            )
        }
        columns[[index]] <- strict_integer_counts(
            table[[record$count_column]], sample
        )
    }
    counts <- do.call(cbind, columns)
    rownames(counts) <- gene_ids
    colnames(counts) <- samples
    counts
}

read_salmon <- function(input) {
    samples <- as_character_vector(input$sample_order, "input.sample_order")
    evidence <- validate_inferential_replicate_summary(input, samples)
    if (identical(input$input_semantics, "salmon_quant_dirs_full_length") &&
        evidence$count == 1) {
        backend_abort(
            "BACKEND_FAILED",
            "Full-length Salmon input with one inferential replicate is unsupported.",
            list(
                reason = "inferential_replicate_count_below_minimum",
                observed_replicates_per_sample = 1,
                minimum_replicates_per_sample = 2
            )
        )
    }
    records <- input$salmon$samples
    quant_dirs <- vapply(records, function(record) record$quant_dir, character(1L))
    record_samples <- vapply(records, function(record) record$sample_id, character(1L))
    if (!identical(record_samples, samples)) {
        backend_abort(
            "BACKEND_FAILED",
            "Salmon snapshots are not aligned to analysis sample order.",
            list(reason = "snapshot_sample_order_mismatch")
        )
    }
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
    files <- stats::setNames(file.path(quant_dirs, "quant.sf"), samples)
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
    imported <- validate_imported_replicates(txi, evidence, samples)
    assert_gene_ids(rownames(txi$counts))
    if (!is.matrix(txi$counts) || any(!is.finite(txi$counts)) ||
        any(txi$counts < 0)) {
        backend_abort(
            "COUNT_VALUES_INVALID",
            "tximport returned invalid estimated gene counts.",
            list(reason = "tximport_count_value_invalid"),
            exit_status = 3L
        )
    }
    list(txi = txi, replicate_evidence = evidence, imported = imported)
}

construct_dataset <- function(request, design_info) {
    input <- request$input
    semantics <- input$input_semantics
    col_data <- S4Vectors::DataFrame(design_info$data)
    if (identical(semantics, "featurecounts_integer")) {
        counts <- if (identical(input$featurecounts$layout, "combined_matrix")) {
            read_combined_featurecounts(input)
        } else if (identical(input$featurecounts$layout, "per_sample_files")) {
            read_per_sample_featurecounts(input)
        } else {
            backend_abort(
                "BACKEND_FAILED",
                "The normalized featureCounts layout is unsupported.",
                list(reason = "normalized_route_invalid")
            )
        }
        audit <- rounding_audit(counts, counts, "none_integer_input")
        dataset <- DESeq2::DESeqDataSetFromMatrix(
            countData = counts,
            colData = col_data,
            design = design_info$matrix
        )
        route <- list(
            constructor = "DESeq2::DESeqDataSetFromMatrix",
            count_source = "validated_featureCounts_integer_counts",
            count_semantics = "integer",
            transcript_length_offset = FALSE,
            gene_length_correction = FALSE,
            rounding_audit = audit,
            inferential_replicates_imported = FALSE,
            inferential_replicates_used_for_inference = FALSE,
            inferential_replicates_unused_reason = "not_applicable"
        )
        return(list(dataset = dataset, route_observed = route))
    }
    if (!semantics %in% c(
        "salmon_quant_dirs_full_length", "salmon_quant_dirs_three_prime"
    )) {
        backend_abort(
            "BACKEND_FAILED",
            "The normalized DESeq2 input route is unsupported.",
            list(reason = "normalized_route_invalid", input_semantics = semantics)
        )
    }
    imported <- read_salmon(input)
    txi <- imported$txi
    rounded <- round(txi$counts)
    if (any(!is.finite(rounded)) || any(rounded < 0) || any(rounded > INTEGER_MAX)) {
        backend_abort(
            "COUNT_VALUES_INVALID",
            "Rounded Salmon counts exceed the DESeq2 signed-integer domain.",
            list(
                reason = "count_value_out_of_deseq2_integer_domain",
                maximum = INTEGER_MAX
            ),
            exit_status = 3L
        )
    }
    storage.mode(rounded) <- "integer"
    mode <- if (identical(semantics, "salmon_quant_dirs_full_length")) {
        "deseq2_constructor"
    } else {
        "explicit_before_matrix_constructor"
    }
    audit <- rounding_audit(txi$counts, rounded, mode)
    if (identical(semantics, "salmon_quant_dirs_full_length")) {
        dataset <- DESeq2::DESeqDataSetFromTximport(
            txi = txi,
            colData = col_data,
            design = design_info$matrix
        )
        if (!identical(DESeq2::counts(dataset), rounded) ||
            !"avgTxLength" %in% SummarizedExperiment::assayNames(dataset)) {
            backend_abort(
                "BACKEND_FAILED",
                "DESeqDataSetFromTximport did not apply the audited constructor route.",
                list(reason = "deseq2_tximport_constructor_mismatch")
            )
        }
        constructor <- "DESeq2::DESeqDataSetFromTximport"
        offset <- TRUE
        correction <- TRUE
        rounding_note <- paste(
            "DESeq2::DESeqDataSetFromTximport internally calls base::round;",
            "the pre/post conversion is audited here."
        )
    } else {
        dataset <- DESeq2::DESeqDataSetFromMatrix(
            countData = rounded,
            colData = col_data,
            design = design_info$matrix
        )
        if ("avgTxLength" %in% SummarizedExperiment::assayNames(dataset)) {
            backend_abort(
                "BACKEND_FAILED",
                "The three-prime DESeq2 route unexpectedly contains avgTxLength.",
                list(reason = "three_prime_offset_present")
            )
        }
        constructor <- "DESeq2::DESeqDataSetFromMatrix"
        offset <- FALSE
        correction <- FALSE
        rounding_note <- paste(
            "The toolkit explicitly calls base::round before",
            "DESeqDataSetFromMatrix; no length offset is attached."
        )
    }
    route <- list(
        constructor = constructor,
        count_source = "txi$counts",
        count_semantics = "salmon_estimated_counts_rounded_for_DESeq2",
        transcript_length_offset = offset,
        gene_length_correction = correction,
        countsFromAbundance = "no",
        dropInfReps = FALSE,
        rounding_audit = audit,
        rounding_disclosure = rounding_note,
        inferential_replicates_imported = imported$imported,
        inferential_replicate_state = imported$replicate_evidence$state,
        inferential_replicates_per_sample = imported$replicate_evidence$count,
        inferential_replicate_method = imported$replicate_evidence$method,
        inferential_replicates_used_for_inference = FALSE,
        inferential_replicates_unused_reason = paste(
            "DESeq2 1.52.0 does not consume tximport infReps in this backend;",
            "they are imported only to verify declared evidence."
        )
    )
    list(dataset = dataset, route_observed = route)
}

long_design_table <- function(design) {
    rows <- vector("list", nrow(design) * ncol(design))
    cursor <- 1L
    for (sample_index in seq_len(nrow(design))) {
        for (column_index in seq_len(ncol(design))) {
            rows[[cursor]] <- data.frame(
                sample_id = rownames(design)[[sample_index]],
                coefficient = colnames(design)[[column_index]],
                value = design[sample_index, column_index],
                stringsAsFactors = FALSE
            )
            cursor <- cursor + 1L
        }
    }
    do.call(rbind, rows)
}

model_convergence <- function(dataset, test_mode) {
    metadata <- S4Vectors::mcols(dataset)
    all_zero <- as.logical(metadata$allZero)
    if (identical(test_mode, "wald")) {
        convergence <- as.logical(metadata$betaConv)
        reason <- rep("wald_model_nonconvergence", length(convergence))
    } else {
        full <- as.logical(metadata$fullBetaConv)
        reduced <- as.logical(metadata$reducedBetaConv)
        convergence <- full & reduced
        convergence[is.na(full) | is.na(reduced)] <- NA
        reason <- rep("lrt_model_nonconvergence", length(convergence))
        reason[!is.na(full) & !full & !is.na(reduced) & reduced] <- (
            "lrt_full_model_nonconvergence"
        )
        reason[!is.na(full) & full & !is.na(reduced) & !reduced] <- (
            "lrt_reduced_model_nonconvergence"
        )
        reason[!is.na(full) & !full & !is.na(reduced) & !reduced] <- (
            "lrt_full_and_reduced_model_nonconvergence"
        )
    }
    if (length(convergence) != nrow(dataset)) {
        backend_abort(
            "BACKEND_FAILED",
            "DESeq2 did not return one convergence diagnostic per gene.",
            list(reason = "model_convergence_diagnostic_missing")
        )
    }
    list(all_zero = all_zero, converged = convergence, failure_reason = reason)
}

coefficient_mapping <- function(dataset, design) {
    result_names <- DESeq2::resultsNames(dataset)
    if (length(result_names) != ncol(design) || anyDuplicated(result_names)) {
        backend_abort(
            "BACKEND_FAILED",
            "DESeq2 coefficient names cannot be mapped one-to-one to the design.",
            list(
                reason = "deseq2_coefficient_mapping_invalid",
                design_columns = colnames(design),
                results_names = result_names
            )
        )
    }
    lapply(seq_along(result_names), function(index) list(
        design_coefficient = colnames(design)[[index]],
        deseq2_result_name = result_names[[index]],
        position = index
    ))
}

long_coefficient_table <- function(dataset, design, diagnostics, mapping) {
    gene_ids <- rownames(dataset)
    result_names <- vapply(mapping, `[[`, character(1L), "deseq2_result_name")
    design_names <- vapply(mapping, `[[`, character(1L), "design_coefficient")
    values <- as.matrix(S4Vectors::mcols(dataset)[, result_names, drop = FALSE])
    if (!identical(dim(values), c(length(gene_ids), length(result_names)))) {
        backend_abort(
            "BACKEND_FAILED",
            "The fitted DESeq2 coefficient matrix has an invalid shape.",
            list(reason = "deseq2_coefficient_shape_invalid")
        )
    }
    rows <- vector("list", length(gene_ids) * length(result_names))
    cursor <- 1L
    for (gene_index in seq_along(gene_ids)) {
        if (isTRUE(diagnostics$all_zero[[gene_index]])) {
            status <- "not_tested"
            status_reason <- "all_zero"
        } else if (!isTRUE(diagnostics$converged[[gene_index]])) {
            status <- "failed"
            status_reason <- diagnostics$failure_reason[[gene_index]]
        } else {
            status <- "tested"
            status_reason <- "fitted"
        }
        for (column_index in seq_along(result_names)) {
            estimate <- values[gene_index, column_index]
            row_status <- status
            row_reason <- status_reason
            if (identical(row_status, "tested") && !is.finite(estimate)) {
                row_status <- "failed"
                row_reason <- "coefficient_unavailable"
            }
            rows[[cursor]] <- data.frame(
                gene_id = gene_ids[[gene_index]],
                status = row_status,
                status_reason = row_reason,
                coefficient = design_names[[column_index]],
                estimate = if (is.finite(estimate)) estimate else NA_real_,
                scale = "log2",
                stringsAsFactors = FALSE
            )
            cursor <- cursor + 1L
        }
    }
    do.call(rbind, rows)
}

same_optional_numeric <- function(left, right) {
    if (length(left) != length(right) || !identical(is.na(left), is.na(right))) {
        return(FALSE)
    }
    finite <- !is.na(left)
    all(left[finite] == right[finite])
}

single_shrinkage_coefficient <- function(contrast, design, mapping) {
    nonzero <- which(contrast$vector != 0)
    if (length(nonzero) != 1L || contrast$vector[[nonzero]] != 1 ||
        identical(colnames(design)[[nonzero]], "(Intercept)")) {
        backend_abort(
            "INVALID_REQUEST",
            "apeglm shrinkage requires one +1 non-intercept coefficient.",
            list(
                reason = "apeglm_contrast_not_single_coefficient",
                contrast_id = contrast$contrast_id
            ),
            exit_status = 2L
        )
    }
    mapping[[nonzero]]$deseq2_result_name
}

contrast_result <- function(
    dataset,
    design,
    contrast,
    test_specification,
    diagnostics,
    mapping,
    cooks_cutoff
) {
    shrinkage <- test_specification$shrinkage
    coefficient_name <- NULL
    if (identical(shrinkage, "apeglm")) {
        coefficient_name <- single_shrinkage_coefficient(
            contrast, design, mapping
        )
    }
    result_arguments <- list(
        object = dataset,
        lfcThreshold = contrast$lfc_threshold,
        altHypothesis = "greaterAbs",
        independentFiltering = TRUE,
        alpha = 0.1,
        pAdjustMethod = "BH",
        parallel = FALSE
    )
    if (is.null(coefficient_name)) {
        result_arguments$contrast <- as.numeric(contrast$vector)
    } else {
        result_arguments$name <- coefficient_name
    }
    result <- do.call(DESeq2::results, result_arguments)

    diagnostic_arguments <- result_arguments
    diagnostic_arguments$independentFiltering <- FALSE
    diagnostic_arguments$cooksCutoff <- FALSE
    diagnostic <- do.call(DESeq2::results, diagnostic_arguments)

    unshrunk <- as.numeric(result$log2FoldChange)
    published <- unshrunk
    shrinkage_method <- "none"
    shrinkage_converged <- rep(TRUE, length(unshrunk))
    if (identical(shrinkage, "apeglm")) {
        shrinkage_fit <- DESeq2::lfcShrink(
            dds = dataset,
            coef = coefficient_name,
            res = result,
            type = "apeglm",
            # The formal threshold belongs to results() above.  Keep apeglm
            # strictly in its effect-estimation role: a positive threshold
            # makes DESeq2 replace PValue/FDR with FSOS s-values.
            lfcThreshold = 0,
            svalue = FALSE,
            quiet = TRUE,
            parallel = FALSE,
            returnList = TRUE
        )
        if (!is.list(shrinkage_fit) || is.null(shrinkage_fit$res) ||
            !is.list(shrinkage_fit$fit) || is.null(shrinkage_fit$fit$diag) ||
            !is.matrix(shrinkage_fit$fit$diag) ||
            !"conv" %in% colnames(shrinkage_fit$fit$diag) ||
            nrow(shrinkage_fit$fit$diag) != length(unshrunk) ||
            (!is.null(rownames(shrinkage_fit$fit$diag)) &&
             !identical(rownames(shrinkage_fit$fit$diag), rownames(dataset)))) {
            backend_abort(
                "BACKEND_FAILED",
                "apeglm did not return aligned convergence diagnostics.",
                list(
                    reason = "apeglm_convergence_diagnostic_invalid",
                    contrast_id = contrast$contrast_id
                )
            )
        }
        shrunken <- shrinkage_fit$res
        convergence_codes <- as.numeric(shrinkage_fit$fit$diag[, "conv"])
        shrinkage_converged <- !is.na(convergence_codes) & convergence_codes == 0
        if (!same_optional_numeric(result$pvalue, shrunken$pvalue) ||
            !same_optional_numeric(result$padj, shrunken$padj)) {
            backend_abort(
                "BACKEND_FAILED",
                "apeglm changed p-values or adjusted p-values.",
                list(
                    reason = "shrinkage_inference_columns_changed",
                    contrast_id = contrast$contrast_id
                )
            )
        }
        published <- as.numeric(shrunken$log2FoldChange)
        shrinkage_method <- "apeglm"
    }

    metadata <- S4Vectors::metadata(result)
    filter_threshold <- as.numeric(metadata$filterThreshold)
    filter_theta <- as.numeric(metadata$filterTheta)
    if (length(filter_threshold) != 1L || !is.finite(filter_threshold) ||
        length(filter_theta) != 1L || !is.finite(filter_theta)) {
        backend_abort(
            "BACKEND_FAILED",
            "DESeq2 independent-filter metadata is unavailable.",
            list(
                reason = "independent_filter_metadata_missing",
                contrast_id = contrast$contrast_id
            )
        )
    }

    base_mean <- as.numeric(result$baseMean)
    standard_error <- as.numeric(result$lfcSE)
    statistic <- as.numeric(result$stat)
    p_value <- as.numeric(result$pvalue)
    fdr <- as.numeric(result$padj)
    diagnostic_p <- as.numeric(diagnostic$pvalue)
    status <- rep(NA_character_, length(p_value))
    status_reason <- rep(NA_character_, length(p_value))
    for (index in seq_along(status)) {
        if (isTRUE(diagnostics$all_zero[[index]])) {
            status[[index]] <- "not_tested"
            status_reason[[index]] <- "all_zero"
        } else if (!isTRUE(diagnostics$converged[[index]])) {
            status[[index]] <- "failed"
            status_reason[[index]] <- diagnostics$failure_reason[[index]]
        } else if (is.na(p_value[[index]]) && is.finite(diagnostic_p[[index]])) {
            status[[index]] <- "failed"
            status_reason[[index]] <- "cooks_outlier"
        } else if (is.finite(p_value[[index]]) && is.na(fdr[[index]])) {
            status[[index]] <- "filtered"
            status_reason[[index]] <- "independent_filtering"
        } else if (is.finite(p_value[[index]]) && is.finite(fdr[[index]])) {
            status[[index]] <- "tested"
            status_reason[[index]] <- "tested"
        } else {
            backend_abort(
                "BACKEND_FAILED",
                "DESeq2 returned an unclassified missing-value pattern.",
                list(
                    reason = "result_status_unclassified",
                    contrast_id = contrast$contrast_id,
                    gene_id = rownames(dataset)[[index]]
                )
            )
        }
        if (identical(shrinkage_method, "apeglm") &&
            status[[index]] %in% c("tested", "filtered") &&
            (!isTRUE(shrinkage_converged[[index]]) ||
             !is.finite(published[[index]]))) {
            status[[index]] <- "failed"
            status_reason[[index]] <- "apeglm_nonconvergence"
            published[[index]] <- NA_real_
        }
    }

    if (identical(test_specification$mode, "lrt")) {
        statistic_type <- "LRT"
        statistic_hypothesis <- "full_vs_reduced_omnibus"
        fdr_basis <- "omnibus_pvalue_BH"
        test_method <- "DESeq2::results_LRT"
        reporting_role <- "reported_effect_not_lrt_hypothesis"
    } else if (contrast$lfc_threshold > 0) {
        statistic_type <- "Wald"
        statistic_hypothesis <- "abs_log2_fold_change_greater_than_threshold"
        fdr_basis <- "contrast_threshold_pvalue_BH_after_independent_filtering"
        test_method <- "DESeq2::results_Wald_greaterAbs"
        reporting_role <- "tested_contrast"
    } else {
        statistic_type <- "Wald"
        statistic_hypothesis <- "contrast_equals_zero"
        fdr_basis <- "contrast_pvalue_BH_after_independent_filtering"
        test_method <- "DESeq2::results_Wald"
        reporting_role <- "tested_contrast"
    }

    table <- data.frame(
        gene_id = rownames(dataset),
        contrast_id = contrast$contrast_id,
        status = status,
        status_reason = status_reason,
        baseMean = base_mean,
        logFC = published,
        unshrunk_logFC = unshrunk,
        lfcSE = standard_error,
        statistic = statistic,
        statistic_type = statistic_type,
        statistic_hypothesis = statistic_hypothesis,
        PValue = p_value,
        FDR = fdr,
        fdr_basis = fdr_basis,
        test_method = test_method,
        lfc_threshold = contrast$lfc_threshold,
        shrinkage_method = shrinkage_method,
        stringsAsFactors = FALSE
    )
    provenance <- list(
        contrast_id = contrast$contrast_id,
        weights = contrast$weights,
        lfc_threshold = contrast$lfc_threshold,
        test_method = test_method,
        alternative_hypothesis = if (contrast$lfc_threshold > 0) {
            "greaterAbs"
        } else if (identical(test_specification$mode, "wald")) {
            "greaterAbs_at_zero_equivalent_two_sided"
        } else {
            "full_vs_reduced_omnibus"
        },
        estimability_residual = contrast$estimability_residual,
        estimability_tolerance = contrast$estimability_tolerance,
        independent_filter_threshold = filter_threshold,
        independent_filter_theta = filter_theta,
        independent_filter_alpha = 0.1,
        cooks_filter_applied = TRUE,
        cooks_cutoff = cooks_cutoff,
        coefficient_name = coefficient_name,
        shrinkage_method = shrinkage_method,
        shrinkage_nonconvergence_count = sum(
            status_reason == "apeglm_nonconvergence"
        )
    )
    reporting_effect <- list(
        contrast_id = contrast$contrast_id,
        role = reporting_role,
        weights = contrast$weights,
        coefficient_name = coefficient_name
    )
    list(
        table = table,
        provenance = provenance,
        reporting_effect = reporting_effect
    )
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

run_analysis <- function(request_path, output_dir) {
    if (file.exists(output_dir)) {
        backend_abort(
            "BACKEND_FAILED",
            "The private DESeq2 output stage already exists.",
            list(reason = "output_stage_exists", output_dir = output_dir)
        )
    }
    runtime <- assert_runtime()
    request <- read_request(request_path)
    assert_request_envelope(request)
    design_info <- build_design(request)
    contrasts <- build_contrasts(request, design_info)
    test_specification <- build_test_specification(
        request, design_info, contrasts
    )
    constructed <- construct_dataset(request, design_info)
    dataset <- constructed$dataset
    if (!identical(colnames(dataset), rownames(design_info$matrix))) {
        backend_abort(
            "DESIGN_RANK_DEFICIENT",
            "The DESeq2 dataset and design sample order differ.",
            list(reason = "design_assay_order_mismatch"),
            exit_status = 3L
        )
    }

    if (identical(test_specification$mode, "wald")) {
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
    } else {
        dataset <- DESeq2::DESeq(
            dataset,
            test = "LRT",
            fitType = "parametric",
            sfType = "ratio",
            betaPrior = FALSE,
            full = design_info$matrix,
            reduced = test_specification$reduced$matrix,
            quiet = TRUE,
            minReplicatesForReplace = 7,
            useT = FALSE,
            minmu = 0.5,
            parallel = FALSE
        )
    }
    resolved_fit_type <- attr(DESeq2::dispersionFunction(dataset), "fitType")
    if (!resolved_fit_type %in% c("parametric", "local")) {
        backend_abort(
            "BACKEND_FAILED",
            "DESeq2 resolved an unsupported dispersion fit type.",
            list(
                reason = "dispersion_fit_type_invalid",
                observed = resolved_fit_type
            )
        )
    }
    mapping <- coefficient_mapping(dataset, design_info$matrix)
    diagnostics <- model_convergence(dataset, test_specification$mode)
    cooks_cutoff <- stats::qf(
        0.99,
        ncol(design_info$matrix),
        nrow(design_info$matrix) - ncol(design_info$matrix)
    )
    result_objects <- lapply(
        contrasts,
        function(contrast) contrast_result(
            dataset,
            design_info$matrix,
            contrast,
            test_specification,
            diagnostics,
            mapping,
            cooks_cutoff
        )
    )
    results <- do.call(rbind, lapply(result_objects, `[[`, "table"))
    coefficients <- long_coefficient_table(
        dataset, design_info$matrix, diagnostics, mapping
    )
    design_table <- long_design_table(design_info$matrix)
    contrast_provenance <- lapply(result_objects, `[[`, "provenance")
    reporting_effect <- lapply(result_objects, `[[`, "reporting_effect")
    status_counts <- as.list(table(factor(results$status, levels = STATUS_VOCABULARY)))
    status_counts <- lapply(status_counts, as.integer)
    names(status_counts) <- STATUS_VOCABULARY
    reduced_public <- if (is.null(test_specification$reduced)) {
        NULL
    } else {
        test_specification$reduced[setdiff(
            names(test_specification$reduced), "matrix"
        )]
    }
    replace_count <- sum(S4Vectors::mcols(dataset)$replace %in% TRUE, na.rm = TRUE)
    defaults <- list(
        fit_type_requested = "parametric",
        fit_type_resolved = resolved_fit_type,
        size_factor_type = "ratio",
        beta_prior = FALSE,
        min_replicates_for_replace = 7,
        use_t = FALSE,
        minmu = 0.5,
        parallel = FALSE,
        independent_filtering = TRUE,
        cooks_cutoff = list(
            requested = "automatic",
            resolved_f_quantile = 0.99,
            resolved_value = cooks_cutoff,
            numerator_df = ncol(design_info$matrix),
            denominator_df = nrow(design_info$matrix) - ncol(design_info$matrix)
        ),
        alpha = 0.1,
        p_adjust_method = "BH",
        results_alt_hypothesis = "greaterAbs",
        outlier_replacement_count = replace_count
    )
    test <- list(
        mode = test_specification$mode,
        shrinkage = test_specification$shrinkage,
        reduced = reduced_public
    )
    analysis <- list(
        schema_version = SCHEMA_VERSION,
        kind = "deseq2_analysis",
        backend = BACKEND_NAME,
        execution_scope = request$execution_scope,
        analysis_request = request$analysis_request,
        input_evidence = request$input_evidence,
        runtime_identity = runtime,
        input_semantics = request$input$input_semantics,
        route_observed = constructed$route_observed,
        pipeline = list(
            list(
                step = "construct_DESeqDataSet",
                constructor = constructed$route_observed$constructor
            ),
            list(
                step = "DESeq",
                arguments = list(
                    test = if (identical(test_specification$mode, "wald")) {
                        "Wald"
                    } else {
                        "LRT"
                    },
                    fitType = "parametric",
                    sfType = "ratio",
                    betaPrior = FALSE,
                    minReplicatesForReplace = 7,
                    useT = FALSE,
                    minmu = 0.5,
                    parallel = FALSE
                )
            ),
            list(
                step = "results",
                arguments = list(
                    independentFiltering = TRUE,
                    cooksCutoff = "automatic",
                    alpha = 0.1,
                    pAdjustMethod = "BH",
                    altHypothesis = "greaterAbs",
                    parallel = FALSE
                )
            ),
            list(
                step = "lfcShrink",
                method = test_specification$shrinkage,
                arguments = if (identical(test_specification$shrinkage, "apeglm")) {
                    list(lfcThreshold = 0, svalue = FALSE)
                } else {
                    NULL
                },
                inferential_columns_changed = FALSE
            )
        ),
        defaults = defaults,
        design = list(
            intercept = design_info$intercept,
            terms = design_info$terms,
            variable_types = design_info$variable_types,
            factor_levels = design_info$factor_levels,
            columns = colnames(design_info$matrix),
            coefficient_mapping = mapping,
            sample_count = nrow(design_info$matrix),
            rank = design_info$rank,
            residual_df = design_info$residual_df,
            qr_tolerance = design_info$qr_tolerance
        ),
        test = test,
        contrasts = contrast_provenance,
        genes = list(
            total = nrow(dataset),
            result_rows = nrow(results),
            status_counts = status_counts
        ),
        status_vocabulary = STATUS_VOCABULARY,
        result_logFC_scale = "log2",
        coefficient_scale = "log2",
        multiple_testing = paste(
            "Benjamini-Hochberg within each reporting contrast after",
            "DESeq2 independent filtering"
        ),
        reporting_effect = reporting_effect
    )

    if (!dir.create(output_dir, recursive = FALSE, mode = "0700")) {
        backend_abort(
            "BACKEND_FAILED",
            "The private DESeq2 output stage could not be created.",
            list(reason = "output_stage_create_failed", output_dir = output_dir)
        )
    }
    paths <- list(
        analysis = file.path(output_dir, "analysis.json"),
        coefficients = file.path(output_dir, "coefficients.tsv"),
        design = file.path(output_dir, "design.tsv"),
        results = file.path(output_dir, "results.tsv")
    )
    write_tsv(results, paths$results)
    write_tsv(coefficients, paths$coefficients)
    write_tsv(design_table, paths$design)
    writeLines(
        jsonlite::toJSON(
            analysis,
            auto_unbox = TRUE,
            null = "null",
            na = "null",
            digits = NA,
            pretty = TRUE
        ),
        con = paths$analysis,
        useBytes = TRUE
    )
    member_names <- c(
        "analysis.json", "coefficients.tsv", "design.tsv", "results.tsv"
    )
    members <- lapply(
        member_names,
        function(name) file_record(file.path(output_dir, name), name)
    )
    manifest <- list(
        schema_version = SCHEMA_VERSION,
        kind = "deseq2_backend_manifest",
        backend = BACKEND_NAME,
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
        design_columns = colnames(design_info$matrix),
        design_rank = design_info$rank,
        residual_df = design_info$residual_df,
        gene_count = nrow(dataset),
        result_status_counts = status_counts,
        test = test,
        defaults = defaults,
        contrasts = contrast_provenance,
        reporting_effect = reporting_effect
    )
    emit_document(list(
        schema_version = SCHEMA_VERSION,
        status = "success",
        backend = BACKEND_NAME,
        data = response_data,
        warnings = list(),
        errors = list(),
        artifacts = c(
            members,
            list(file_record(manifest_path, "backend_manifest.json"))
        )
    ))
}

main <- function() {
    arguments <- commandArgs(trailingOnly = TRUE)
    if (length(arguments) != 2L) {
        backend_abort(
            "BACKEND_FAILED",
            "Usage: deseq2.R NORMALIZED_REQUEST_JSON PRIVATE_OUTPUT_STAGE",
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
            schema_version = SCHEMA_VERSION,
            status = "error",
            backend = BACKEND_NAME,
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
            schema_version = SCHEMA_VERSION,
            status = "error",
            backend = BACKEND_NAME,
            data = NULL,
            warnings = list(),
            errors = list(list(
                code = "BACKEND_FAILED",
                message = "The locked DESeq2 backend failed unexpectedly.",
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
