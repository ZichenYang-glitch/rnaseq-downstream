"""Data-integrity regression tests for the retained legacy loaders."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest

pd = pytest.importorskip(
    "pandas",
    reason="Legacy data-integrity tests require the optional pandas runtime.",
)
pytest.importorskip(
    "yaml",
    reason="Importing the legacy data module requires the optional PyYAML runtime.",
)

from modules.data import (  # noqa: E402
    load_annotation_table,
    load_counts,
    load_expression_matrix,
    load_metadata,
    merge_gene_name_annotations,
)
from rnaseq_downstream.errors import (  # noqa: E402
    CountValuesInvalidError,
    GeneIdentifierError,
    InputReadError,
    InputValidationError,
    SampleSetMismatchError,
)


def _write_tsv(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.unit
@pytest.mark.parametrize(
    ("loader", "kwargs", "operation"),
    [
        (load_expression_matrix, {}, "load_expression_matrix"),
        (load_annotation_table, {}, "load_annotation_table"),
    ],
)
def test_direct_table_read_failures_remain_typed_io_errors(
    tmp_path: Path,
    loader: Callable[..., object],
    kwargs: dict[str, object],
    operation: str,
) -> None:
    missing = tmp_path / "missing.tsv"

    with pytest.raises(InputReadError) as captured:
        loader(missing, **kwargs)

    assert captured.value.details["path"] == str(missing)
    assert captured.value.details["operation"] == operation


@pytest.mark.unit
@pytest.mark.parametrize(
    "header",
    [
        "gene_name\ts1\nTP53\t10\n",
        "identifier\ts1\nENSG000001\t10\n",
    ],
)
def test_expression_matrix_requires_explicit_gene_id(
    tmp_path: Path,
    header: str,
) -> None:
    path = _write_tsv(tmp_path / "matrix.tsv", header)

    with pytest.raises(GeneIdentifierError) as captured:
        load_expression_matrix(path)

    assert captured.value.details == {
        "path": str(path),
        "column": "gene_id",
        "reason": "column_missing",
        "available_columns": header.splitlines()[0].split("\t"),
    }


@pytest.mark.unit
def test_gene_name_is_display_only_and_never_becomes_the_index(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "matrix.tsv",
        "gene_id\tgene_name\ts1\nENSG000002\tSAME\t2\nENSG000001\tSAME\t1\n",
    )

    matrix = load_expression_matrix(path, metadata_samples=["s1"])

    assert matrix.index.tolist() == ["ENSG000002", "ENSG000001"]
    assert matrix.columns.tolist() == ["s1"]
    assert matrix["s1"].tolist() == [2, 1]


@pytest.mark.unit
def test_embedded_gene_names_are_preserved_as_one_to_one_display_metadata(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "matrix.tsv",
        (
            "gene_id\tgene_name\ts1\n"
            "ENSG000002\tSAME\t2\n"
            "ENSG000001\tNA\t1\n"
            "ENSG000003\t\t3\n"
        ),
    )

    matrix, annotation = load_expression_matrix(
        path,
        metadata_samples=["s1"],
        return_gene_annotations=True,
    )

    assert matrix.index.tolist() == ["ENSG000002", "ENSG000001", "ENSG000003"]
    assert annotation.index.tolist() == matrix.index.tolist()
    assert annotation.loc["ENSG000002", "gene_name"] == "SAME"
    assert annotation.loc["ENSG000001", "gene_name"] == "NA"
    assert pd.isna(annotation.loc["ENSG000003", "gene_name"])


@pytest.mark.unit
def test_conflicting_gene_name_display_mappings_hard_fail() -> None:
    first = pd.DataFrame(
        {"gene_name": ["TP53", "SHARED"]},
        index=pd.Index(["ENSG1", "ENSG2"], name="gene_id"),
    )
    second = pd.DataFrame(
        {"gene_name": ["P53", "SHARED"]},
        index=pd.Index(["ENSG1", "ENSG3"], name="gene_id"),
    )

    with pytest.raises(GeneIdentifierError) as captured:
        merge_gene_name_annotations(first, second)

    assert captured.value.details["reason"] == "ambiguous_gene_name_mapping"
    assert captured.value.details["conflicting_gene_ids"] == ["ENSG1"]


@pytest.mark.unit
@pytest.mark.parametrize("loader", [load_expression_matrix, load_counts])
def test_duplicate_sample_headers_fail_before_pandas_can_mangle_them(
    tmp_path: Path,
    loader: Callable[..., object],
) -> None:
    path = _write_tsv(
        tmp_path / "duplicate-samples.tsv",
        "gene_id\ts1\ts1\nENSG000001\t1\t2\n",
    )

    with pytest.raises(SampleSetMismatchError) as captured:
        loader(path, metadata_samples=["s1"])

    assert captured.value.details == {
        "path": str(path),
        "reason": "duplicate_matrix_samples",
        "duplicate_samples": ["s1"],
        "occurrences": {"s1": 2},
    }


@pytest.mark.unit
def test_duplicate_gene_id_headers_fail_before_pandas_can_mangle_them(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "duplicate-gene-id.tsv",
        "gene_id\tgene_id\ts1\nENSG000001\tENSG000002\t1\n",
    )

    with pytest.raises(GeneIdentifierError) as captured:
        load_expression_matrix(path, metadata_samples=["s1"])

    assert captured.value.details["reason"] == "duplicate_gene_id_column"
    assert captured.value.details["occurrences"] == 2


@pytest.mark.unit
def test_empty_expression_header_fails_before_pandas_can_mangle_it(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "empty-header.tsv",
        "gene_id\t\nENSG000001\t1\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_expression_matrix(path)

    assert captured.value.details["reason"] == "empty_expression_column_name"
    assert captured.value.details["column_positions"] == [2]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("filename", "contents", "loader", "kwargs"),
    [
        (
            "expression.tsv",
            "gene_id\ts1\nENSG000001\t1\textra\n",
            load_expression_matrix,
            {},
        ),
        (
            "annotation.tsv",
            "gene_id\tgene_name\nENSG000001\tTP53\textra\n",
            load_annotation_table,
            {},
        ),
        (
            "metadata.tsv",
            "sample_id\tcondition\ns1\tcontrol\textra\n",
            load_metadata,
            {"design_col": "condition"},
        ),
        (
            "expression.tsv",
            "gene_id\ts1\nENSG000001\n",
            load_expression_matrix,
            {},
        ),
        (
            "annotation.tsv",
            "gene_id\tgene_name\nENSG000001\n",
            load_annotation_table,
            {},
        ),
        (
            "metadata.tsv",
            "sample_id\tcondition\ns1\n",
            load_metadata,
            {"design_col": "condition"},
        ),
    ],
)
def test_mismatched_fields_cannot_be_silently_reinterpreted_or_dropped(
    tmp_path: Path,
    filename: str,
    contents: str,
    loader: Callable[..., object],
    kwargs: dict[str, object],
) -> None:
    path = _write_tsv(tmp_path / filename, contents)

    with pytest.raises(InputValidationError) as captured:
        loader(path, **kwargs)

    details = captured.value.details
    assert details["path"] == str(path)
    assert details["reason"] == "row_field_count_mismatch"
    assert details["expected_field_count"] == 2
    assert details["invalid_row_count"] == 1
    assert details["examples"][0]["row_number"] == 2
    assert details["examples"][0]["observed_field_count"] in {1, 3}


@pytest.mark.unit
def test_row_width_validation_respects_quoted_tsv_fields(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "expression.tsv",
        'gene_id\tgene_name\ts1\nENSG000001\t"TP\t53"\t1\n',
    )

    matrix, annotation = load_expression_matrix(
        path,
        metadata_samples=["s1"],
        return_gene_annotations=True,
    )

    assert matrix.loc["ENSG000001", "s1"] == 1
    assert annotation.loc["ENSG000001", "gene_name"] == "TP\t53"


@pytest.mark.unit
def test_gene_ids_are_read_as_identifiers_without_numeric_coercion(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "matrix.tsv",
        "gene_id\ts1\n00123\t1\n",
    )

    matrix = load_expression_matrix(path, strip_gene_version=False)

    assert matrix.index.tolist() == ["00123"]


@pytest.mark.unit
def test_expression_matrix_rejects_empty_gene_or_sample_axes(tmp_path: Path) -> None:
    no_genes = _write_tsv(tmp_path / "no-genes.tsv", "gene_id\ts1\n")
    no_samples = _write_tsv(tmp_path / "no-samples.tsv", "gene_id\nENSG000001\n")

    with pytest.raises(GeneIdentifierError) as gene_error:
        load_expression_matrix(no_genes)
    with pytest.raises(InputValidationError) as sample_error:
        load_expression_matrix(no_samples)

    assert gene_error.value.details["reason"] == "no_gene_rows"
    assert sample_error.value.details["reason"] == "no_sample_columns"


@pytest.mark.unit
@pytest.mark.parametrize(
    ("rows", "strip_gene_version", "collisions"),
    [
        (
            "ENSG000001\t1\nENSG000001\t2\n",
            False,
            {"ENSG000001": ["ENSG000001"]},
        ),
        (
            "ENSG000001.1\t1\nENSG000001.2\t2\n",
            True,
            {"ENSG000001": ["ENSG000001.1", "ENSG000001.2"]},
        ),
    ],
)
def test_duplicate_or_version_colliding_gene_ids_hard_fail(
    tmp_path: Path,
    rows: str,
    strip_gene_version: bool,
    collisions: dict[str, list[str]],
) -> None:
    path = _write_tsv(tmp_path / "matrix.tsv", f"gene_id\ts1\n{rows}")

    with pytest.raises(GeneIdentifierError) as captured:
        load_expression_matrix(path, strip_gene_version=strip_gene_version)

    assert captured.value.details["reason"] == "duplicate_or_normalization_collision"
    assert captured.value.details["collisions"] == collisions


@pytest.mark.unit
@pytest.mark.parametrize(
    ("gene_id", "reason", "strip_gene_version"),
    [
        ("", "missing_or_empty", True),
        (" ENSG000001", "surrounding_whitespace", True),
        (".1", "empty_after_normalization", True),
    ],
)
def test_empty_or_whitespace_padded_gene_ids_hard_fail(
    tmp_path: Path,
    gene_id: str,
    reason: str,
    strip_gene_version: bool,
) -> None:
    path = _write_tsv(tmp_path / "matrix.tsv", f"gene_id\ts1\n{gene_id}\t1\n")

    with pytest.raises(GeneIdentifierError) as captured:
        load_expression_matrix(path, strip_gene_version=strip_gene_version)

    assert captured.value.details["reason"] == reason


@pytest.mark.unit
def test_sample_sets_must_match_exactly_and_report_both_directions(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    path = _write_tsv(
        tmp_path / "matrix.tsv",
        "gene_id\ts1\ts3\nENSG000001\t1\t3\n",
    )

    with pytest.raises(SampleSetMismatchError) as captured:
        load_expression_matrix(path, metadata_samples=["s2", "s1"])

    assert captured.value.details == {
        "path": str(path),
        "missing_from_matrix": ["s2"],
        "extra_in_matrix": ["s3"],
        "metadata_sample_count": 2,
        "matrix_sample_count": 2,
    }
    assert capsys.readouterr() == ("", "")


@pytest.mark.unit
def test_count_sample_mismatch_is_not_reclassified_as_an_io_error(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "matrix.tsv",
        "gene_id\ts1\nENSG000001\t1\n",
    )

    with pytest.raises(SampleSetMismatchError) as captured:
        load_counts(path, metadata_samples=["s2"], min_counts=0)

    assert not isinstance(captured.value, InputReadError)
    assert captured.value.details["missing_from_matrix"] == ["s2"]
    assert captured.value.details["extra_in_matrix"] == ["s1"]


@pytest.mark.unit
def test_metadata_order_deterministically_controls_expression_columns(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "matrix.tsv",
        "gene_id\ts1\ts2\nENSG000001\t1\t2\n",
    )

    matrix = load_expression_matrix(path, metadata_samples=["s2", "s1"])

    assert matrix.columns.tolist() == ["s2", "s1"]
    assert matrix.loc["ENSG000001"].tolist() == [2, 1]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("value", "reason"),
    [
        ("", "missing_expression_value"),
        ("not-numeric", "non_numeric_expression_value"),
        ("inf", "non_finite_expression_value"),
    ],
)
def test_generic_expression_values_are_numeric_finite_and_never_imputed(
    tmp_path: Path,
    value: str,
    reason: str,
) -> None:
    path = _write_tsv(
        tmp_path / "expression.tsv",
        f"gene_id\ts1\nENSG000001\t{value}\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_expression_matrix(path, metadata_samples=["s1"])

    assert captured.value.details["reason"] == reason
    assert captured.value.details["invalid_value_count"] == 1


@pytest.mark.unit
def test_generic_expression_preserves_fractional_values(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "expression.tsv",
        "gene_id\ts1\nENSG000001\t1.25\n",
    )

    matrix = load_expression_matrix(path, metadata_samples=["s1"])

    assert matrix.loc["ENSG000001", "s1"] == pytest.approx(1.25)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("value", "reason"),
    [
        ("", "missing"),
        ("not-a-count", "non_numeric"),
        ("inf", "non_finite"),
        ("1.25", "fractional"),
        ("-1", "negative"),
    ],
)
def test_counts_reject_invalid_values_before_integer_conversion(
    tmp_path: Path,
    value: str,
    reason: str,
) -> None:
    path = _write_tsv(
        tmp_path / f"{reason}.tsv",
        f"gene_id\ts1\nENSG000001\t{value}\n",
    )

    with pytest.raises(CountValuesInvalidError) as captured:
        load_counts(path, metadata_samples=["s1"], min_counts=0)

    assert captured.value.details["reason"] == reason
    assert captured.value.details["invalid_value_count"] == 1
    assert captured.value.details["examples"][0]["gene_id"] == "ENSG000001"
    assert captured.value.details["examples"][0]["sample_id"] == "s1"


@pytest.mark.unit
@pytest.mark.parametrize(
    ("value", "reason"),
    [
        ("True", "non_numeric"),
        ("False", "non_numeric"),
        ("1.0", "non_integer_lexeme"),
        ("1e0", "non_integer_lexeme"),
        ("1.0000000000000001", "fractional"),
        ("9007199254740993.0", "non_integer_lexeme"),
        (str(2**63), "integer_out_of_range"),
    ],
)
def test_counts_are_parsed_from_exact_integer_lexemes_without_loss(
    tmp_path: Path,
    value: str,
    reason: str,
) -> None:
    path = _write_tsv(
        tmp_path / "counts.tsv",
        f"gene_id\ts1\nENSG000001\t{value}\n",
    )

    with pytest.raises(CountValuesInvalidError) as captured:
        load_counts(path, metadata_samples=["s1"], min_counts=0)

    assert captured.value.details["reason"] == reason
    assert captured.value.details["examples"][0]["value"] == repr(value)


@pytest.mark.unit
def test_counts_preserve_exact_large_integers_and_leading_zero_lexemes(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "counts.tsv",
        f"gene_id\ts1\nsmall\t0007\nlarge\t{2**53 + 1}\n",
    )

    counts = load_counts(path, metadata_samples=["s1"], min_counts=0)

    assert counts.loc["s1", "small"] == 7
    assert counts.loc["s1", "large"] == 2**53 + 1


@pytest.mark.unit
def test_extremely_long_count_lexeme_returns_a_structured_range_error(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "counts.tsv",
        f"gene_id\ts1\nENSG000001\t{'9' * 5000}\n",
    )

    with pytest.raises(CountValuesInvalidError) as captured:
        load_counts(path, metadata_samples=["s1"], min_counts=0)

    assert captured.value.details["reason"] == "integer_out_of_range"
    assert captured.value.details["invalid_value_count"] == 1


@pytest.mark.unit
def test_count_filter_is_an_inclusive_total_count_threshold(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "counts.tsv",
        "gene_id\ts2\ts1\nlow\t2\t2\nboundary\t4\t6\nkeep\t7\t8\n",
    )

    counts = load_counts(path, metadata_samples=["s1", "s2"], min_counts=10)

    assert counts.index.tolist() == ["s1", "s2"]
    assert counts.columns.tolist() == ["boundary", "keep"]
    assert counts.loc["s1"].tolist() == [6, 8]
    assert all(dtype == "int64" for dtype in counts.dtypes.astype(str))


@pytest.mark.unit
@pytest.mark.parametrize("min_counts", [-1, 1.5, True])
def test_count_filter_threshold_must_be_a_nonnegative_integer(
    tmp_path: Path,
    min_counts: object,
) -> None:
    path = _write_tsv(
        tmp_path / "counts.tsv",
        "gene_id\ts1\nENSG000001\t1\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_counts(path, metadata_samples=["s1"], min_counts=min_counts)

    assert captured.value.details["reason"] == "invalid_min_counts"


@pytest.mark.unit
def test_filtered_counts_return_only_matching_display_annotations(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "counts.tsv",
        ("gene_id\tgene_name\ts1\n" "filtered\tDROP\t1\n" "retained\tKEEP\t10\n"),
    )

    counts, annotation = load_counts(
        path,
        metadata_samples=["s1"],
        min_counts=10,
        return_gene_annotations=True,
    )

    assert counts.columns.tolist() == ["retained"]
    assert annotation.index.tolist() == ["retained"]
    assert annotation.loc["retained", "gene_name"] == "KEEP"


@pytest.mark.unit
def test_annotation_symbols_are_optional_and_missing_values_are_preserved(
    tmp_path: Path,
) -> None:
    with_symbols = _write_tsv(
        tmp_path / "annotation-with-symbols.tsv",
        "stable_id\tsymbol\nENSG000001\tTP53\nENSG000002\t\n",
    )
    without_symbols = _write_tsv(
        tmp_path / "annotation-without-symbols.tsv",
        "stable_id\nENSG000001\nENSG000002\n",
    )

    annotated = load_annotation_table(
        with_symbols,
        gene_id_col="stable_id",
        gene_name_col="symbol",
    )
    unannotated = load_annotation_table(
        without_symbols,
        gene_id_col="stable_id",
        gene_name_col="symbol",
    )

    assert annotated.index.tolist() == ["ENSG000001", "ENSG000002"]
    assert annotated.loc["ENSG000001", "gene_name"] == "TP53"
    assert pd.isna(annotated.loc["ENSG000002", "gene_name"])
    assert unannotated.index.tolist() == ["ENSG000001", "ENSG000002"]
    assert unannotated["gene_name"].isna().all()


@pytest.mark.unit
@pytest.mark.parametrize(
    ("suffix", "delimiter"),
    [(".tsv", "\t"), (".txt", "\t"), (".csv", ",")],
)
def test_duplicate_annotation_gene_id_headers_fail_before_pandas_mangling(
    tmp_path: Path,
    suffix: str,
    delimiter: str,
) -> None:
    path = _write_tsv(
        tmp_path / f"annotation{suffix}",
        delimiter.join(["gene_id", "gene_id", "gene_name"])
        + "\n"
        + delimiter.join(["ENSG000001", "ENSG000002", "TP53"])
        + "\n",
    )

    with pytest.raises(GeneIdentifierError) as captured:
        load_annotation_table(path)

    assert captured.value.details == {
        "path": str(path),
        "column": "gene_id",
        "reason": "duplicate_gene_id_column",
        "occurrences": 2,
    }


@pytest.mark.unit
def test_duplicate_annotation_display_headers_fail_before_pandas_mangling(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "annotation.tsv",
        "gene_id\tgene_name\tgene_name\nENSG000001\tTP53\tP53\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_annotation_table(path)

    assert captured.value.details == {
        "path": str(path),
        "reason": "duplicate_annotation_columns",
        "duplicate_columns": ["gene_name"],
        "occurrences": {"gene_name": 2},
    }


@pytest.mark.unit
def test_empty_annotation_header_fails_before_pandas_can_mangle_it(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "annotation.tsv",
        "gene_id\t\nENSG000001\tTP53\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_annotation_table(path)

    assert captured.value.details["reason"] == "empty_annotation_column_name"
    assert captured.value.details["column_positions"] == [2]


@pytest.mark.unit
def test_annotation_duplicate_primary_ids_are_not_silently_dropped(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "annotation.tsv",
        "gene_id\tgene_name\nENSG000001\tTP53\nENSG000001\tP53\n",
    )

    with pytest.raises(GeneIdentifierError) as captured:
        load_annotation_table(path)

    assert captured.value.details["collisions"] == {"ENSG000001": ["ENSG000001"]}


@pytest.mark.unit
def test_metadata_semantic_errors_are_not_mislabeled_as_io_failures(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "metadata.tsv",
        "sample_id\tbatch\ns1\tb1\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_metadata(path, design_col="condition")

    assert not isinstance(captured.value, InputReadError)
    assert captured.value.details["reason"] == "design_column_missing"


@pytest.mark.unit
def test_valid_metadata_preserves_lexical_sample_ids_and_explicit_types(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "metadata.tsv",
        "sample_id\tcondition\tbatch\tage\n001\tcontrol\t1\t20\n002\ttreated\t2\t30\n",
    )

    metadata = load_metadata(path, design_col="condition")

    assert metadata.index.tolist() == ["001", "002"]
    assert metadata["batch"].tolist() == ["1", "2"]
    assert metadata["age"].tolist() == ["20", "30"]


@pytest.mark.unit
def test_metadata_leading_zero_sample_id_never_matches_a_different_header(
    tmp_path: Path,
) -> None:
    metadata_path = _write_tsv(
        tmp_path / "metadata.tsv",
        "sample_id\tcondition\n001\tcontrol\n",
    )
    matrix_path = _write_tsv(
        tmp_path / "matrix.tsv",
        "gene_id\t1\nENSG000001\t1\n",
    )
    metadata = load_metadata(metadata_path, design_col="condition")

    with pytest.raises(SampleSetMismatchError):
        load_expression_matrix(matrix_path, metadata_samples=metadata.index)


@pytest.mark.unit
def test_metadata_blank_design_value_is_a_semantic_failure(tmp_path: Path) -> None:
    path = _write_tsv(
        tmp_path / "metadata.tsv",
        "sample_id\tcondition\ns1\t\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_metadata(path, design_col="condition")

    assert captured.value.details["reason"] == "missing_design_values"
    assert captured.value.details["samples"] == ["s1"]


@pytest.mark.unit
def test_tabular_metadata_does_not_shift_a_missing_middle_field(
    tmp_path: Path,
) -> None:
    path = _write_tsv(
        tmp_path / "metadata.tsv",
        "sample_id\tcondition\tbatch\ns1\t\tb1\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_metadata(path, design_col="condition")

    assert captured.value.details["reason"] == "missing_design_values"
    assert captured.value.details["samples"] == ["s1"]


@pytest.mark.unit
@pytest.mark.parametrize("delimiter", ["\t", " "])
def test_duplicate_metadata_headers_fail_before_pandas_mangling(
    tmp_path: Path,
    delimiter: str,
) -> None:
    path = _write_tsv(
        tmp_path / "metadata.tsv",
        delimiter.join(["sample_id", "condition", "condition"])
        + "\n"
        + delimiter.join(["s1", "control", "treated"])
        + "\n",
    )

    with pytest.raises(InputValidationError) as captured:
        load_metadata(path, design_col="condition")

    assert captured.value.details == {
        "path": str(path),
        "reason": "duplicate_metadata_columns",
        "duplicate_columns": ["condition"],
        "occurrences": {"condition": 2},
    }


@pytest.mark.unit
def test_metadata_rejects_a_header_without_sample_rows(tmp_path: Path) -> None:
    path = _write_tsv(tmp_path / "metadata.tsv", "sample_id\tcondition\n")

    with pytest.raises(InputValidationError) as captured:
        load_metadata(path, design_col="condition")

    assert captured.value.details["reason"] == "no_sample_rows"
