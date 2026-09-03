# Frozen limma gene-set testing layer for the edgeR QL backend.
#
# This file is sourced by edger_ql.R.  It deliberately contains no top-level
# execution: the caller supplies the already-filtered, TMM-normalized DGEGLM
# fit, its exact design, and the validated contrast vectors.

PATHWAY_RESULT_COLUMNS <- c(
    "contrast_id", "gene_set_id", "gene_set_description", "method_id",
    "test_class", "hypothesis", "inference_role", "status",
    "status_reason", "direction", "proportion_down", "proportion_up",
    "p_value", "fdr", "fdr_family_id",
    "gmt_member_count_raw", "gmt_symbol_count_unique",
    "mapped_symbol_count_unique", "ambiguous_symbol_count_unique",
    "unmapped_symbol_count_unique", "mapping_rate",
    "mapped_gene_id_count_unique", "tested_gene_count",
    "filtered_gene_count", "tested_universe_gene_count", "method_ngenes",
    "correlation_status", "correlation_estimate_raw",
    "correlation_effective", "vif_used", "rotation_status"
)

PATHWAY_MROAST_NROT <- 9999L
PATHWAY_BASE_SEED <- 1729L
PATHWAY_RNG <- list(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
)

pathway_assert_exact_keys <- function(value, required, field) {
    observed <- names(value)
    if (!is.list(value) || is.null(observed) || anyDuplicated(observed) > 0L ||
        !identical(sort(observed), sort(required))) {
        backend_abort(
            "BACKEND_FAILED",
            paste0("The normalized pathway field '", field, "' is incompatible."),
            list(
                reason = "normalized_pathway_request_invalid",
                field = field,
                missing_fields = sort(setdiff(required, observed)),
                unexpected_fields = sort(setdiff(observed, required))
            )
        )
    }
    invisible(value)
}

pathway_has_control <- function(value) {
    any(grepl("[[:cntrl:]]", value, perl = TRUE))
}

pathway_string <- function(value, field) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value) || pathway_has_control(value)) {
        backend_abort(
            "BACKEND_FAILED",
            paste0("The normalized pathway field '", field, "' is invalid."),
            list(
                reason = "normalized_pathway_request_invalid",
                field = field
            )
        )
    }
    value
}

pathway_positive_integer <- function(value, field) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0 || value != floor(value) || value > .Machine$integer.max) {
        backend_abort(
            "BACKEND_FAILED",
            paste0("The normalized pathway field '", field, "' is invalid."),
            list(
                reason = "normalized_pathway_request_invalid",
                field = field
            )
        )
    }
    as.integer(value)
}

pathway_sha256 <- function(value, field) {
    value <- pathway_string(value, field)
    if (!grepl("^[0-9a-f]{64}$", value, perl = TRUE)) {
        backend_abort(
            "BACKEND_FAILED",
            paste0("The normalized pathway digest '", field, "' is invalid."),
            list(
                reason = "normalized_pathway_request_invalid",
                field = field
            )
        )
    }
    value
}

pathway_read_source <- function(source, role) {
    path <- pathway_string(source$path, paste0("gene_sets.", role, ".path"))
    if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path, perl = TRUE) ||
        nzchar(Sys.readlink(path)) || !file.exists(path) ||
        isTRUE(file.info(path)$isdir)) {
        backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source is not a regular local file.",
            list(
                reason = "pathway_source_invalid",
                role = role
            )
        )
    }
    expected_digest <- pathway_sha256(
        source$sha256, paste0("gene_sets.", role, ".sha256")
    )
    expected_size <- pathway_positive_integer(
        source$size_bytes, paste0("gene_sets.", role, ".size_bytes")
    )
    before <- file.info(path)
    observed_size <- unname(before$size)
    if (!identical(as.numeric(observed_size), as.numeric(expected_size))) {
        backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source does not match its evidence.",
            list(
                reason = "pathway_source_identity_mismatch",
                role = role,
                expected_sha256 = expected_digest,
                observed_sha256 = NULL,
                expected_size_bytes = expected_size,
                observed_size_bytes = observed_size
            )
        )
    }
    raw <- readBin(path, what = "raw", n = expected_size)
    if (length(raw) != expected_size) {
        backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source could not be captured completely.",
            list(reason = "pathway_source_short_read", role = role)
        )
    }
    after <- file.info(path)
    observed_digest <- digest::digest(
        raw, algo = "sha256", serialize = FALSE
    )
    source_stable <- identical(
        c(as.numeric(before$size), as.numeric(before$mtime)),
        c(as.numeric(after$size), as.numeric(after$mtime))
    )
    if (!source_stable || !identical(observed_digest, expected_digest)) {
        backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source does not match its evidence.",
            list(
                reason = "pathway_source_identity_mismatch",
                role = role,
                expected_sha256 = expected_digest,
                observed_sha256 = observed_digest,
                expected_size_bytes = expected_size,
                observed_size_bytes = unname(after$size),
                source_identity_stable = source_stable
            )
        )
    }
    text <- tryCatch(
        rawToChar(raw),
        error = function(error) backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source is not valid UTF-8 text.",
            list(
                reason = "pathway_source_utf8_invalid",
                role = role,
                cause_type = class(error)[1L]
            )
        )
    )
    if (!validUTF8(text)) {
        backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source is not valid UTF-8 text.",
            list(reason = "pathway_source_utf8_invalid", role = role)
        )
    }
    Encoding(text) <- "UTF-8"
    text
}

pathway_lines <- function(text, role) {
    if (!nzchar(text)) {
        backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source is empty.",
            list(reason = "pathway_source_empty", role = role)
        )
    }
    lines <- strsplit(paste0(text, "\n.__RNASEQ_EOF__"), "\n", fixed = TRUE)[[1L]]
    lines <- lines[-length(lines)]
    if (length(lines) > 0L && !nzchar(lines[[length(lines)]])) {
        lines <- lines[-length(lines)]
    }
    if (length(lines) == 0L || any(!nzchar(lines))) {
        backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source contains an empty row.",
            list(reason = "pathway_source_empty_row", role = role)
        )
    }
    lines
}

pathway_split_tsv_line <- function(line) {
    fields <- strsplit(
        paste0(line, "\t.__RNASEQ_FIELD_END__"), "\t", fixed = TRUE
    )[[1L]]
    fields[-length(fields)]
}

pathway_assert_fields <- function(fields, role, row_number) {
    if (length(fields) == 0L || any(!nzchar(fields)) ||
        pathway_has_control(fields)) {
        backend_abort(
            "BACKEND_FAILED",
            "A private frozen pathway source contains an invalid field.",
            list(
                reason = "pathway_source_field_invalid",
                role = role,
                row_number = row_number
            )
        )
    }
    invisible(fields)
}

pathway_parse_gmt <- function(text) {
    lines <- pathway_lines(text, "gmt")
    records <- vector("list", length(lines))
    identifiers <- character(length(lines))
    for (row_number in seq_along(lines)) {
        fields <- pathway_split_tsv_line(lines[[row_number]])
        if (length(fields) < 3L) {
            backend_abort(
                "BACKEND_FAILED",
                "A GMT row must contain an ID, description, and at least one symbol.",
                list(
                    reason = "gmt_row_width_invalid",
                    row_number = row_number,
                    observed_field_count = length(fields)
                )
            )
        }
        pathway_assert_fields(fields, "gmt", row_number)
        members <- fields[3:length(fields)]
        if (anyDuplicated(members) > 0L) {
            backend_abort(
                "BACKEND_FAILED",
                "A GMT gene set contains a duplicate symbol.",
                list(
                    reason = "gmt_duplicate_symbol",
                    row_number = row_number,
                    gene_set_id = fields[[1L]],
                    duplicate_symbols = sort(unique(members[duplicated(members)]))
                )
            )
        }
        identifiers[[row_number]] <- fields[[1L]]
        records[[row_number]] <- list(
            gene_set_id = fields[[1L]],
            description = fields[[2L]],
            symbols = members,
            member_count_raw = length(members),
            symbol_count_unique = length(unique(members))
        )
    }
    duplicates <- sort(unique(identifiers[duplicated(identifiers)]))
    if (length(duplicates) > 0L) {
        backend_abort(
            "BACKEND_FAILED",
            "GMT gene-set identifiers must be unique.",
            list(
                reason = "gmt_duplicate_gene_set_id",
                duplicate_gene_set_ids = duplicates
            )
        )
    }
    records[order(identifiers, method = "radix")]
}

pathway_parse_annotation <- function(text, strip_version) {
    lines <- pathway_lines(text, "annotation")
    header <- pathway_split_tsv_line(lines[[1L]])
    if (!identical(header, c("gene_id", "symbol"))) {
        backend_abort(
            "BACKEND_FAILED",
            "The frozen annotation header must be exactly gene_id and symbol.",
            list(
                reason = "annotation_header_invalid",
                observed_header = header
            )
        )
    }
    if (length(lines) < 2L) {
        backend_abort(
            "BACKEND_FAILED",
            "The frozen annotation contains no mapping rows.",
            list(reason = "annotation_empty")
        )
    }
    gene_ids <- character(length(lines) - 1L)
    symbols <- character(length(lines) - 1L)
    for (row_number in 2:length(lines)) {
        fields <- pathway_split_tsv_line(lines[[row_number]])
        if (length(fields) != 2L) {
            backend_abort(
                "BACKEND_FAILED",
                "A frozen annotation row must contain exactly two fields.",
                list(
                    reason = "annotation_row_width_invalid",
                    row_number = row_number,
                    observed_field_count = length(fields)
                )
            )
        }
        pathway_assert_fields(fields, "annotation", row_number)
        position <- row_number - 1L
        gene_ids[[position]] <- strip_gene_versions(fields[[1L]], strip_version)
        symbols[[position]] <- fields[[2L]]
    }
    if (any(!nzchar(gene_ids)) || pathway_has_control(gene_ids)) {
        backend_abort(
            "BACKEND_FAILED",
            "A frozen annotation gene ID is invalid after version handling.",
            list(reason = "annotation_gene_id_invalid_after_version_policy")
        )
    }
    duplicate_gene_ids <- sort(unique(gene_ids[duplicated(gene_ids)]))
    if (length(duplicate_gene_ids) > 0L) {
        backend_abort(
            "BACKEND_FAILED",
            "Frozen annotation gene IDs must be unique after version handling.",
            list(
                reason = "annotation_gene_id_duplicated_after_version_policy",
                duplicate_gene_ids = duplicate_gene_ids
            )
        )
    }
    by_symbol <- split(gene_ids, symbols)
    ambiguous <- names(by_symbol)[lengths(by_symbol) > 1L]
    unique_groups <- by_symbol[lengths(by_symbol) == 1L]
    unique_map <- if (length(unique_groups) == 0L) {
        setNames(character(), character())
    } else {
        setNames(vapply(unique_groups, `[[`, character(1L), 1L), names(unique_groups))
    }
    list(
        unique_map = unique_map,
        ambiguous_symbols = sort(ambiguous, method = "radix"),
        row_count = length(gene_ids),
        unique_symbol_count = length(by_symbol),
        uniquely_mapped_symbol_count = length(unique_map),
        ambiguous_symbol_count = length(ambiguous)
    )
}

pathway_validate_specification <- function(specification) {
    pathway_assert_exact_keys(
        specification,
        c("gmt", "annotation", "minimum_tested_genes"),
        "gene_sets"
    )
    pathway_assert_exact_keys(
        specification$gmt,
        c(
            "path", "declared_path", "sha256", "size_bytes", "collection",
            "version", "identifier_type"
        ),
        "gene_sets.gmt"
    )
    pathway_assert_exact_keys(
        specification$annotation,
        c(
            "path", "declared_path", "sha256", "size_bytes", "name",
            "version", "gene_id_column", "symbol_column"
        ),
        "gene_sets.annotation"
    )
    gmt <- specification$gmt
    annotation <- specification$annotation
    pathway_string(gmt$declared_path, "gene_sets.gmt.declared_path")
    pathway_string(gmt$collection, "gene_sets.gmt.collection")
    pathway_string(gmt$version, "gene_sets.gmt.version")
    if (!identical(
        pathway_string(gmt$identifier_type, "gene_sets.gmt.identifier_type"),
        "symbol"
    )) {
        backend_abort(
            "BACKEND_FAILED",
            "The GMT identifier type must be symbol.",
            list(reason = "pathway_identifier_type_invalid")
        )
    }
    pathway_string(annotation$declared_path, "gene_sets.annotation.declared_path")
    pathway_string(annotation$name, "gene_sets.annotation.name")
    pathway_string(annotation$version, "gene_sets.annotation.version")
    if (!identical(
        pathway_string(
            annotation$gene_id_column,
            "gene_sets.annotation.gene_id_column"
        ),
        "gene_id"
    ) || !identical(
        pathway_string(
            annotation$symbol_column,
            "gene_sets.annotation.symbol_column"
        ),
        "symbol"
    )) {
        backend_abort(
            "BACKEND_FAILED",
            "The frozen annotation columns must be gene_id and symbol.",
            list(reason = "pathway_annotation_columns_invalid")
        )
    }
    list(
        gmt_text = pathway_read_source(gmt, "gmt"),
        annotation_text = pathway_read_source(annotation, "annotation"),
        minimum_tested_genes = pathway_positive_integer(
            specification$minimum_tested_genes,
            "gene_sets.minimum_tested_genes"
        )
    )
}

pathway_map_sets <- function(records, annotation, all_gene_ids, tested_gene_ids) {
    all_gene_set <- unique(all_gene_ids)
    tested_gene_set <- unique(tested_gene_ids)
    lapply(records, function(record) {
        symbols <- record$symbols
        ambiguous <- symbols[symbols %in% annotation$ambiguous_symbols]
        uniquely_mapped <- symbols[symbols %in% names(annotation$unique_map)]
        unmapped <- setdiff(symbols, c(ambiguous, uniquely_mapped))
        mapped_gene_ids <- sort(
            unique(unname(annotation$unique_map[uniquely_mapped])),
            method = "radix"
        )
        tested <- mapped_gene_ids[mapped_gene_ids %in% tested_gene_set]
        filtered <- mapped_gene_ids[
            mapped_gene_ids %in% all_gene_set & !mapped_gene_ids %in% tested_gene_set
        ]
        c(record, list(
            mapped_symbol_count_unique = length(unique(uniquely_mapped)),
            ambiguous_symbol_count_unique = length(unique(ambiguous)),
            unmapped_symbol_count_unique = length(unique(unmapped)),
            mapping_rate = length(unique(uniquely_mapped)) /
                record$symbol_count_unique,
            mapped_gene_ids = mapped_gene_ids,
            tested_gene_ids = sort(unique(tested), method = "radix"),
            filtered_gene_ids = sort(unique(filtered), method = "radix")
        ))
    })
}

pathway_index_list <- function(mapped_sets, tested_gene_ids, eligible) {
    selected <- mapped_sets[vapply(mapped_sets, eligible, logical(1L))]
    if (length(selected) == 0L) return(list())
    indices <- lapply(
        selected,
        function(record) match(record$tested_gene_ids, tested_gene_ids)
    )
    names(indices) <- vapply(selected, `[[`, character(1L), "gene_set_id")
    if (any(vapply(indices, anyNA, logical(1L)))) {
        backend_abort(
            "BACKEND_FAILED",
            "A mapped pathway index is outside the tested gene universe.",
            list(reason = "pathway_index_invalid")
        )
    }
    indices
}

pathway_ordered_set_hash <- function(index, tested_gene_ids) {
    document <- lapply(names(index), function(identifier) list(
        gene_set_id = identifier,
        gene_ids = as.list(unname(tested_gene_ids[index[[identifier]]]))
    ))
    canonical <- jsonlite::toJSON(
        document,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = NA,
        pretty = FALSE
    )
    digest::digest(charToRaw(canonical), algo = "sha256", serialize = FALSE)
}

pathway_assert_method_table <- function(table, index, columns, method_id) {
    expected_ids <- names(index)
    if (!is.data.frame(table) || !identical(rownames(table), expected_ids) ||
        !all(columns %in% colnames(table))) {
        backend_abort(
            "BACKEND_FAILED",
            "A limma pathway method returned an incompatible result table.",
            list(
                reason = "pathway_method_result_schema_invalid",
                method_id = method_id,
                expected_gene_set_ids = expected_ids,
                observed_gene_set_ids = rownames(table),
                observed_columns = colnames(table)
            )
        )
    }
    numeric_columns <- intersect(
        c(
            "NGenes", "PropDown", "PropUp", "PValue", "FDR",
            "PValue.Mixed", "FDR.Mixed", "Correlation"
        ),
        columns
    )
    bad <- numeric_columns[vapply(
        numeric_columns,
        function(column) any(!is.finite(table[[column]])),
        logical(1L)
    )]
    if (length(bad) > 0L || any(!table$Direction %in% c("Up", "Down"))) {
        backend_abort(
            "BACKEND_FAILED",
            "A limma pathway method returned a non-finite or invalid result.",
            list(
                reason = "pathway_method_result_invalid",
                method_id = method_id,
                invalid_columns = bad
            )
        )
    }
    if (any(table$NGenes != lengths(index))) {
        backend_abort(
            "BACKEND_FAILED",
            "A limma pathway method tested an unexpected number of genes.",
            list(reason = "pathway_method_gene_count_mismatch", method_id = method_id)
        )
    }
    invisible(table)
}

pathway_verify_bh <- function(table, p_column, fdr_column, method_id, hypothesis) {
    expected <- stats::p.adjust(table[[p_column]], method = "BH")
    observed <- table[[fdr_column]]
    difference <- if (length(expected) == 0L) 0 else max(abs(expected - observed))
    if (!is.finite(difference) || difference > 1e-12) {
        backend_abort(
            "BACKEND_FAILED",
            "A limma pathway FDR column disagrees with an independent BH calculation.",
            list(
                reason = "pathway_fdr_mismatch",
                method_id = method_id,
                hypothesis = hypothesis,
                maximum_absolute_difference = difference
            )
        )
    }
    expected
}

pathway_empty_row <- function(
    contrast_id, record, method_id, test_class, hypothesis, inference_role,
    tested_universe_count
) {
    data.frame(
        contrast_id = contrast_id,
        gene_set_id = record$gene_set_id,
        gene_set_description = record$description,
        method_id = method_id,
        test_class = test_class,
        hypothesis = hypothesis,
        inference_role = inference_role,
        status = "not_tested",
        status_reason = "",
        direction = "",
        proportion_down = NA_real_,
        proportion_up = NA_real_,
        p_value = NA_real_,
        fdr = NA_real_,
        fdr_family_id = paste(contrast_id, method_id, hypothesis, sep = "|"),
        gmt_member_count_raw = record$member_count_raw,
        gmt_symbol_count_unique = record$symbol_count_unique,
        mapped_symbol_count_unique = record$mapped_symbol_count_unique,
        ambiguous_symbol_count_unique = record$ambiguous_symbol_count_unique,
        unmapped_symbol_count_unique = record$unmapped_symbol_count_unique,
        mapping_rate = record$mapping_rate,
        mapped_gene_id_count_unique = length(record$mapped_gene_ids),
        tested_gene_count = length(record$tested_gene_ids),
        filtered_gene_count = length(record$filtered_gene_ids),
        tested_universe_gene_count = tested_universe_count,
        method_ngenes = NA_real_,
        correlation_status = if (identical(method_id, "limma_camera")) {
            "not_estimated_not_tested"
        } else {
            "not_applicable"
        },
        correlation_estimate_raw = NA_real_,
        correlation_effective = NA_real_,
        vif_used = NA_real_,
        rotation_status = if (identical(method_id, "limma_mroast")) {
            "not_performed_not_tested"
        } else if (identical(method_id, "limma_fry")) {
            "not_applicable_analytic_approximation"
        } else {
            "not_applicable"
        },
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

pathway_not_tested_reason <- function(record, method_id, minimum, universe_count) {
    count <- length(record$tested_gene_ids)
    if (count < minimum) return("tested_gene_count_below_minimum")
    if (identical(method_id, "limma_camera") && count < 2L) {
        return("camera_requires_at_least_two_tested_genes")
    }
    if (identical(method_id, "limma_camera") && count == universe_count) {
        return("competitive_test_requires_background_genes")
    }
    NULL
}

pathway_method_rows <- function(
    contrast_id, mapped_sets, method_id, test_class, hypotheses,
    inference_role, tested_universe_count, minimum, table = NULL
) {
    rows <- list()
    cursor <- 1L
    for (record in mapped_sets) {
        reason <- pathway_not_tested_reason(
            record, method_id, minimum, tested_universe_count
        )
        for (hypothesis in hypotheses) {
            row <- pathway_empty_row(
                contrast_id, record, method_id, test_class, hypothesis,
                inference_role, tested_universe_count
            )
            if (!is.null(reason)) {
                row$status_reason <- reason
            } else {
                result <- table[record$gene_set_id, , drop = FALSE]
                row$status <- "tested"
                if (identical(hypothesis, "directional")) {
                    row$direction <- result$Direction[[1L]]
                } else {
                    row$direction <- "Mixed"
                }
                row$method_ngenes <- result$NGenes[[1L]]
                if (identical(hypothesis, "directional")) {
                    row$p_value <- result$PValue[[1L]]
                    row$fdr <- result$FDR[[1L]]
                } else {
                    row$p_value <- result$PValue.Mixed[[1L]]
                    row$fdr <- result$FDR.Mixed[[1L]]
                }
                if (identical(method_id, "limma_mroast")) {
                    row$proportion_down <- result$PropDown[[1L]]
                    row$proportion_up <- result$PropUp[[1L]]
                    row$rotation_status <- "performed_fixed_seed_9999_rotations"
                }
                if (identical(method_id, "limma_camera")) {
                    raw <- result$Correlation[[1L]]
                    effective <- max(0, raw)
                    row$correlation_status <- "estimated_set_specific"
                    row$correlation_estimate_raw <- raw
                    row$correlation_effective <- effective
                    row$vif_used <- max(1, 1 + (result$NGenes[[1L]] - 1) * raw)
                }
            }
            rows[[cursor]] <- row
            cursor <- cursor + 1L
        }
    }
    do.call(rbind, rows)
}

pathway_run_methods <- function(fit, design, contrast, mapped_sets, minimum, seed) {
    tested_gene_ids <- rownames(fit$counts)
    universe_count <- length(tested_gene_ids)
    self_index <- pathway_index_list(
        mapped_sets,
        tested_gene_ids,
        function(record) length(record$tested_gene_ids) >= minimum
    )
    camera_index <- pathway_index_list(
        mapped_sets,
        tested_gene_ids,
        function(record) {
            count <- length(record$tested_gene_ids)
            count >= minimum && count >= 2L && count < universe_count
        }
    )

    mroast_table <- NULL
    fry_table <- NULL
    camera_table <- NULL
    if (length(self_index) > 0L) {
        RNGkind(
            kind = PATHWAY_RNG$kind,
            normal.kind = PATHWAY_RNG$normal.kind,
            sample.kind = PATHWAY_RNG$sample.kind
        )
        set.seed(seed)
        mroast_table <- limma::mroast(
            fit,
            index = self_index,
            design = design,
            contrast = contrast$vector,
            geneid = rownames(fit$counts),
            set.statistic = "mean",
            gene.weights = NULL,
            nrot = PATHWAY_MROAST_NROT,
            adjust.method = "BH",
            midp = FALSE,
            sort = "none"
        )
        pathway_assert_method_table(
            mroast_table,
            self_index,
            c(
                "NGenes", "PropDown", "PropUp", "Direction", "PValue",
                "FDR", "PValue.Mixed", "FDR.Mixed"
            ),
            "limma_mroast"
        )
        mroast_table$FDR <- pathway_verify_bh(
            mroast_table, "PValue", "FDR", "limma_mroast", "directional"
        )
        mroast_table$FDR.Mixed <- pathway_verify_bh(
            mroast_table, "PValue.Mixed", "FDR.Mixed", "limma_mroast", "mixed"
        )

        fry_table <- limma::fry(
            fit,
            index = self_index,
            design = design,
            contrast = contrast$vector,
            geneid = rownames(fit$counts),
            sort = "none"
        )
        pathway_assert_method_table(
            fry_table,
            self_index,
            c(
                "NGenes", "Direction", "PValue", "FDR", "PValue.Mixed",
                "FDR.Mixed"
            ),
            "limma_fry"
        )
        fry_table$FDR <- pathway_verify_bh(
            fry_table, "PValue", "FDR", "limma_fry", "directional"
        )
        fry_table$FDR.Mixed <- pathway_verify_bh(
            fry_table, "PValue.Mixed", "FDR.Mixed", "limma_fry", "mixed"
        )
    }
    if (length(camera_index) > 0L) {
        camera_table <- limma::camera(
            fit,
            index = camera_index,
            design = design,
            contrast = contrast$vector,
            weights = NULL,
            use.ranks = FALSE,
            allow.neg.cor = FALSE,
            inter.gene.cor = NA_real_,
            sort = FALSE
        )
        pathway_assert_method_table(
            camera_table,
            camera_index,
            c("NGenes", "Correlation", "Direction", "PValue", "FDR"),
            "limma_camera"
        )
        camera_table$FDR <- pathway_verify_bh(
            camera_table, "PValue", "FDR", "limma_camera", "directional"
        )
    }

    rows <- rbind(
        pathway_method_rows(
            contrast$contrast_id, mapped_sets, "limma_mroast",
            "self_contained", c("directional", "mixed"), "corroborative",
            universe_count, minimum, mroast_table
        ),
        pathway_method_rows(
            contrast$contrast_id, mapped_sets, "limma_fry",
            "self_contained", c("directional", "mixed"), "primary",
            universe_count, minimum, fry_table
        ),
        pathway_method_rows(
            contrast$contrast_id, mapped_sets, "limma_camera", "competitive",
            "directional", "supplementary", universe_count, minimum,
            camera_table
        )
    )
    rows <- rows[
        order(
            rows$contrast_id, rows$gene_set_id, rows$method_id,
            rows$hypothesis, method = "radix"
        ),
        PATHWAY_RESULT_COLUMNS,
        drop = FALSE
    ]
    rownames(rows) <- NULL
    list(
        rows = rows,
        provenance = list(
            contrast_id = contrast$contrast_id,
            gene_level_lfc_threshold = contrast$lfc_threshold,
            pathway_statistical_null = "zero_effect",
            gene_level_lfc_threshold_applied_to_pathways = FALSE,
            ordered_set_lists = list(
                self_contained = list(
                    gene_set_ids = as.list(names(self_index)),
                    sha256 = pathway_ordered_set_hash(self_index, tested_gene_ids)
                ),
                competitive = list(
                    gene_set_ids = as.list(names(camera_index)),
                    sha256 = pathway_ordered_set_hash(camera_index, tested_gene_ids)
                )
            ),
            rotation = list(
                method_id = "limma_mroast",
                seed = seed,
                seed_policy = "base_1729_plus_zero_based_declared_contrast_index",
                rng = PATHWAY_RNG,
                nrot = PATHWAY_MROAST_NROT,
                reset_before_each_contrast = TRUE
            )
        )
    )
}

run_pathway_tests <- function(
    specification, input, all_gene_ids, keep, fit, design, contrasts
) {
    validated <- pathway_validate_specification(specification)
    parsed_gmt <- pathway_parse_gmt(validated$gmt_text)
    annotation <- pathway_parse_annotation(
        validated$annotation_text,
        strict_boolean_scalar(
            input$gene_id_policy$strip_version,
            "input.gene_id_policy.strip_version"
        )
    )
    tested_gene_ids <- rownames(fit$counts)
    if (!identical(tested_gene_ids, all_gene_ids[keep])) {
        backend_abort(
            "BACKEND_FAILED",
            "The pathway tested universe is not aligned to the QL fit.",
            list(reason = "pathway_tested_universe_mismatch")
        )
    }
    mapped_sets <- pathway_map_sets(
        parsed_gmt, annotation, all_gene_ids, tested_gene_ids
    )
    executions <- lapply(seq_along(contrasts), function(index) {
        run <- pathway_run_methods(
            fit, design, contrasts[[index]], mapped_sets,
            validated$minimum_tested_genes,
            PATHWAY_BASE_SEED + index - 1L
        )
        run
    })
    rows <- do.call(rbind, lapply(executions, `[[`, "rows"))
    rows <- rows[
        order(
            rows$contrast_id, rows$gene_set_id, rows$method_id,
            rows$hypothesis, method = "radix"
        ),
        PATHWAY_RESULT_COLUMNS,
        drop = FALSE
    ]
    rownames(rows) <- NULL
    if (!identical(colnames(rows), PATHWAY_RESULT_COLUMNS) ||
        any(!rows$status %in% c("tested", "not_tested")) ||
        any(rows$status == "tested" & rows$method_ngenes != rows$tested_gene_count) ||
        any(!is.finite(rows$mapping_rate)) ||
        any(rows$mapping_rate < 0 | rows$mapping_rate > 1) ||
        any(
            rows$status == "tested" & rows$method_id == "limma_mroast" &
            (!is.finite(rows$proportion_down) |
                !is.finite(rows$proportion_up) |
                rows$proportion_down < 0 | rows$proportion_down > 1 |
                rows$proportion_up < 0 | rows$proportion_up > 1)
        ) ||
        any(
            rows$status == "tested" &
            (!is.finite(rows$p_value) | !is.finite(rows$fdr) |
                rows$p_value < 0 | rows$p_value > 1 |
                rows$fdr < 0 | rows$fdr > 1)
        )) {
        backend_abort(
            "BACKEND_FAILED",
            "The pathway result table violates its internal output contract.",
            list(reason = "pathway_result_contract_invalid")
        )
    }

    gmt <- specification$gmt
    annotation_source <- specification$annotation
    provenance <- list(
        gene_sets = list(
            gmt = list(
                collection = gmt$collection,
                version = gmt$version,
                identifier_type = gmt$identifier_type,
                sha256 = gmt$sha256,
                size_bytes = gmt$size_bytes,
                gene_set_count = length(parsed_gmt)
            ),
            annotation = list(
                name = annotation_source$name,
                version = annotation_source$version,
                gene_id_column = annotation_source$gene_id_column,
                symbol_column = annotation_source$symbol_column,
                sha256 = annotation_source$sha256,
                size_bytes = annotation_source$size_bytes,
                row_count = annotation$row_count
            ),
            minimum_tested_genes = validated$minimum_tested_genes,
            sets = lapply(mapped_sets, function(record) list(
                gene_set_id = record$gene_set_id,
                gene_set_description = record$description,
                gmt_member_count_raw = record$member_count_raw,
                gmt_symbol_count_unique = record$symbol_count_unique,
                mapped_symbol_count_unique = record$mapped_symbol_count_unique,
                ambiguous_symbol_count_unique = record$ambiguous_symbol_count_unique,
                unmapped_symbol_count_unique = record$unmapped_symbol_count_unique,
                mapping_rate = record$mapping_rate,
                mapped_gene_id_count_unique = length(record$mapped_gene_ids),
                tested_gene_count = length(record$tested_gene_ids),
                filtered_gene_count = length(record$filtered_gene_ids),
                absent_from_assay_gene_count = length(record$mapped_gene_ids) -
                    length(record$tested_gene_ids) -
                    length(record$filtered_gene_ids)
            ))
        ),
        mapping_policy = list(
            source_identifier = "symbol",
            target_identifier = "stable_gene_id",
            annotation_gene_id_version_stripping = isTRUE(
                input$gene_id_policy$strip_version
            ),
            one_to_many_symbols = "ambiguous_excluded",
            duplicate_gmt_members = "hard_fail",
            mapping_rate = paste(
                "uniquely_mapped_unique_symbols divided by",
                "unique_GMT_symbols"
            ),
            tested_membership = "intersection_with_filterByExpr_tested_universe",
            not_tested_policy = "tested_gene_count_below_minimum"
        ),
        tested_universe_gene_count = length(tested_gene_ids),
        methods = list(
            limma_mroast = list(
                generic = "limma::mroast",
                dispatch = "edgeR::mroast.DGEGLM",
                test_class = "self_contained",
                inference_role = "corroborative",
                statistical_null = "zero_effect",
                parameters = list(
                    set.statistic = "mean",
                    gene.weights = NULL,
                    nrot = PATHWAY_MROAST_NROT,
                    adjust.method = "BH",
                    midp = FALSE,
                    sort = "none"
                )
            ),
            limma_fry = list(
                generic = "limma::fry",
                dispatch = "edgeR::fry.DGEGLM",
                test_class = "self_contained",
                inference_role = "primary",
                statistical_null = "zero_effect",
                parameters = list(sort = "none")
            ),
            limma_camera = list(
                generic = "limma::camera",
                dispatch = "edgeR::camera.DGEGLM",
                test_class = "competitive",
                inference_role = "supplementary",
                statistical_null = "no_enrichment_relative_to_complement",
                parameters = list(
                    weights = NULL,
                    use.ranks = FALSE,
                    allow.neg.cor = FALSE,
                    inter.gene.cor = "estimated_per_gene_set",
                    sort = FALSE
                )
            )
        ),
        multiple_testing = list(
            method = "Benjamini-Hochberg",
            scope = paste(
                "separately within contrast x method x hypothesis across",
                "tested gene sets only"
            ),
            family_id_format = "contrast_id|method_id|hypothesis",
            python_independent_recalculation_required = TRUE
        ),
        contrasts = lapply(executions, `[[`, "provenance")
    )
    list(results = rows, provenance = provenance)
}
