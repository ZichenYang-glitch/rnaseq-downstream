"""NumPy-only tests for the visualization QC mathematics."""

from __future__ import annotations

from pathlib import Path

import pytest

np = pytest.importorskip(
    "numpy",
    reason="Visualization QC math tests require the separately reported NumPy lane.",
)

from rnaseq_downstream.errors import (  # noqa: E402
    CovariateConfoundedError,
    ErrorCode,
    QCValidationError,
    ToolkitError,
)
from rnaseq_downstream.qc_math import (  # noqa: E402
    FactorSpec,
    build_joint_design,
    centered_unscaled_pca,
    metadata_reorder_indices,
    remove_nuisance_effects,
    resolve_adjustment_biology_factors,
    select_top_variable_features,
    validate_factor_names,
)

pytestmark = pytest.mark.unit
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _assert_reason(
    exc_info: pytest.ExceptionInfo[ToolkitError],
    reason: str,
    code: ErrorCode = ErrorCode.QC_VALIDATION_FAILED,
) -> None:
    assert exc_info.value.code is code
    assert exc_info.value.details["reason"] == reason


def test_top_variable_selection_uses_positive_sample_variance() -> None:
    values = np.asarray(
        [
            [0.0, 3.0, 10.0, 8.0],
            [1.0, 3.0, 8.0, 4.0],
            [2.0, 3.0, 6.0, 0.0],
        ]
    )

    selected = select_top_variable_features(
        values, ["low", "constant", "medium", "high"], top_n=2
    )

    assert selected.names == ("high", "medium")
    np.testing.assert_array_equal(selected.indices, [3, 2])
    np.testing.assert_allclose(selected.variances, [16.0, 4.0])


def test_top_variable_selection_breaks_variance_ties_by_feature_id() -> None:
    values = np.asarray([[0.0, 3.0], [1.0, 2.0], [2.0, 1.0]])

    selected = select_top_variable_features(values, ["gene_b", "gene_a"], top_n=2)

    assert selected.names == ("gene_a", "gene_b")
    np.testing.assert_array_equal(selected.indices, [1, 0])


def test_top_variable_selection_rejects_duplicate_feature_ids() -> None:
    with pytest.raises(QCValidationError) as exc_info:
        select_top_variable_features([[0.0, 1.0], [2.0, 3.0]], ["g", "g"], 2)

    _assert_reason(exc_info, "duplicate_feature_ids")


def test_centered_pca_is_unscaled_and_has_stable_signs() -> None:
    values = np.asarray(
        [
            [-100.0, -1.0, 5.0],
            [-50.0, 1.0, 5.0],
            [50.0, -1.0, 5.0],
            [100.0, 1.0, 5.0],
        ]
    )

    result = centered_unscaled_pca(values, n_components=2)

    np.testing.assert_allclose(result.coordinates.mean(axis=0), 0.0, atol=1e-12)
    assert abs(result.components[0, 0]) > 0.999
    assert result.explained_variance_ratio[0] > 0.999
    np.testing.assert_allclose(
        result.coordinates @ result.components,
        values - result.feature_means,
        atol=1e-12,
    )
    for component in result.components:
        pivot = int(np.argmax(np.abs(component)))
        assert component[pivot] >= 0.0


def test_centered_pca_does_not_mutate_input() -> None:
    values = np.asarray([[0.0, 2.0, 5.0], [1.0, 0.0, 5.0], [3.0, 1.0, 5.0]])
    original = values.copy()

    centered_unscaled_pca(values)

    np.testing.assert_array_equal(values, original)


def test_centered_pca_rejects_second_component_for_two_samples() -> None:
    values = np.asarray(
        [
            [0.0, 1.0, 2.0],
            [3.0, 4.0, 5.0],
        ]
    )

    with pytest.raises(QCValidationError) as exc_info:
        centered_unscaled_pca(values, n_components=2)

    _assert_reason(exc_info, "insufficient_pca_dimensions")
    assert exc_info.value.details["shape"] == [2, 3]
    assert exc_info.value.details["maximum_components"] == 1


def test_centered_pca_accepts_the_only_direction_for_two_samples() -> None:
    result = centered_unscaled_pca(
        [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]],
        n_components=1,
    )

    assert result.coordinates.shape == (2, 1)
    assert result.components.shape == (1, 3)
    assert result.explained_variance_ratio.shape == (1,)


def test_centered_pca_preserves_feature_directions_when_samples_are_plentiful() -> None:
    values = np.asarray(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [2.0, 1.0],
            [1.0, 3.0],
        ]
    )

    result = centered_unscaled_pca(values, n_components=2)

    assert result.coordinates.shape == (5, 2)
    assert result.components.shape == (2, 2)
    assert result.explained_variance_ratio.shape == (2,)


def test_metadata_alignment_is_lossless_and_explicit() -> None:
    reorder = metadata_reorder_indices(
        ["sample_2", "sample_1", "sample_3"],
        ["sample_1", "sample_2", "sample_3"],
    )

    np.testing.assert_array_equal(reorder, [1, 0, 2])


@pytest.mark.parametrize(
    ("expression_ids", "metadata_ids", "reason"),
    [
        (["s1", "s1"], ["s1", "s2"], "duplicate_expression_samples"),
        (["s1", "s2"], ["s1", "s1"], "duplicate_metadata_samples"),
        (["s1", "s2"], ["s1", "s3"], "sample_mismatch"),
    ],
)
def test_metadata_alignment_rejects_ambiguous_or_missing_samples(
    expression_ids: list[str], metadata_ids: list[str], reason: str
) -> None:
    with pytest.raises(QCValidationError) as exc_info:
        metadata_reorder_indices(expression_ids, metadata_ids)

    _assert_reason(exc_info, reason)


@pytest.mark.parametrize(
    ("available", "biology", "nuisance", "reason"),
    [
        ({"group", "batch"}, ["group"], ["missing"], "unknown_factors"),
        ({"group", "batch"}, ["group"], ["batch", "batch"], "duplicate_factors"),
        ({"group", "batch"}, ["group"], ["group"], "overlapping_factors"),
    ],
)
def test_factor_name_validation_hard_fails(
    available: set[str], biology: list[str], nuisance: list[str], reason: str
) -> None:
    with pytest.raises(QCValidationError) as exc_info:
        validate_factor_names(available, biology, nuisance)

    _assert_reason(exc_info, reason)


def test_adjustment_biology_resolution_always_protects_primary_factor() -> None:
    protected = resolve_adjustment_biology_factors(
        "~ group + sex + batch",
        "group",
        ["sex"],
        ["batch"],
    )

    assert protected == ("group", "sex")


def test_adjustment_biology_resolution_requires_complete_additive_coverage() -> None:
    with pytest.raises(QCValidationError) as exc_info:
        resolve_adjustment_biology_factors(
            "~ group + sex + batch",
            "group",
            ["group"],
            ["batch"],
        )

    _assert_reason(exc_info, "incomplete_biology_protection")
    assert exc_info.value.details["missing_biology_factors"] == ["sex"]


@pytest.mark.parametrize(
    "design",
    [
        "~ group * sex",
        "~ group:sex",
        "~ group + splines(age)",
        "~ group - 1",
        "~ 0 + group",
        "~ group + 0",
        "~ 1 + group",
        "~ group + 1",
    ],
)
def test_adjustment_biology_resolution_rejects_complex_formulae(design: str) -> None:
    with pytest.raises(QCValidationError) as exc_info:
        resolve_adjustment_biology_factors(
            design,
            "group",
            ["group", "sex", "age"],
            ["batch"],
        )

    _assert_reason(exc_info, "unsupported_adjustment_design")


def test_primary_design_factor_cannot_be_declared_nuisance() -> None:
    with pytest.raises(QCValidationError) as exc_info:
        resolve_adjustment_biology_factors(
            "~ group",
            "group",
            ["group"],
            ["group"],
        )

    _assert_reason(exc_info, "primary_design_factor_is_nuisance")


def test_joint_adjustment_removes_only_categorical_nuisance() -> None:
    groups = np.asarray(["A", "A", "A", "B", "B", "B"] * 2)
    batches = np.asarray(["X"] * 6 + ["Y"] * 6)
    design = build_joint_design(
        [FactorSpec("group", groups, "categorical")],
        [FactorSpec("batch", batches, "categorical")],
    )
    group_effect = (groups == "B").astype(float)
    batch_effect = (batches == "Y").astype(float)
    expression = np.column_stack(
        [
            10.0 + 5.0 * group_effect + 7.0 * batch_effect,
            3.0 - 2.0 * group_effect + 4.0 * batch_effect,
        ]
    )

    adjusted = remove_nuisance_effects(expression, design)

    for group in ("A", "B"):
        x_mean = adjusted[(groups == group) & (batches == "X")].mean(axis=0)
        y_mean = adjusted[(groups == group) & (batches == "Y")].mean(axis=0)
        np.testing.assert_allclose(x_mean, y_mean, atol=1e-12)
    observed_biology = adjusted[groups == "B"].mean(axis=0) - adjusted[
        groups == "A"
    ].mean(axis=0)
    np.testing.assert_allclose(observed_biology, [5.0, -2.0], atol=1e-12)
    assert design.biology_columns == ("intercept", "group[B]")
    assert design.nuisance_columns == ("batch[sum:X]",)
    assert design.residual_df == 9


def test_joint_adjustment_encodes_centered_numeric_nuisance() -> None:
    groups = np.asarray(["A", "B"] * 4)
    age = np.asarray([20.0, 21.0, 25.0, 26.0, 30.0, 31.0, 35.0, 36.0])
    design = build_joint_design(
        [FactorSpec("group", groups, "categorical")],
        [FactorSpec("age", age, "numeric")],
    )
    expression = np.column_stack(
        [2.0 + 3.0 * (groups == "B") + 0.5 * (age - age.mean())]
    )

    adjusted = remove_nuisance_effects(expression, design)

    np.testing.assert_allclose(adjusted[groups == "A"], 2.0, atol=1e-12)
    np.testing.assert_allclose(adjusted[groups == "B"], 5.0, atol=1e-12)
    np.testing.assert_allclose(design.nuisance.mean(axis=0), 0.0, atol=1e-12)


def test_joint_adjustment_handles_unbalanced_but_estimable_factors() -> None:
    groups = np.asarray(["A", "A", "A", "A", "B", "B", "B", "B", "B"])
    batches = np.asarray(["X", "X", "X", "Y", "X", "Y", "Y", "Y", "Y"])
    group_effect = (groups == "B").astype(float)
    batch_effect = (batches == "Y").astype(float)
    expression = np.column_stack(
        [
            8.0 + 4.0 * group_effect + 6.0 * batch_effect,
            5.0 - 3.0 * group_effect + 2.0 * batch_effect,
        ]
    )
    design = build_joint_design(
        [FactorSpec("group", groups, "categorical")],
        [FactorSpec("batch", batches, "categorical")],
    )

    adjusted = remove_nuisance_effects(expression, design)

    # Locked limma 3.68 removeBatchEffect oracle: categorical batch uses
    # contr.sum, which centers over levels rather than observed sample counts.
    expected = np.column_stack(
        [
            11.0 + 4.0 * group_effect,
            6.0 - 3.0 * group_effect,
        ]
    )
    np.testing.assert_allclose(adjusted, expected, atol=1e-12)
    for group in ("A", "B"):
        x_mean = adjusted[(groups == group) & (batches == "X")].mean(axis=0)
        y_mean = adjusted[(groups == group) & (batches == "Y")].mean(axis=0)
        np.testing.assert_allclose(x_mean, y_mean, atol=1e-12)
    observed_biology = adjusted[groups == "B"].mean(axis=0) - adjusted[
        groups == "A"
    ].mean(axis=0)
    np.testing.assert_allclose(observed_biology, [4.0, -3.0], atol=1e-12)
    np.testing.assert_array_equal(
        design.nuisance[:, 0],
        np.where(batches == "X", 1.0, -1.0),
    )


def test_joint_adjustment_preserves_multiple_biology_factors() -> None:
    groups = np.asarray(["A", "A", "B", "B"] * 3)
    sexes = np.asarray(["F", "M", "F", "M"] * 3)
    batches = np.asarray(["X"] * 4 + ["Y"] * 4 + ["Z"] * 4)
    group_effect = (groups == "B").astype(float)
    sex_effect = (sexes == "M").astype(float)
    batch_levels = {"X": -2.0, "Y": 1.0, "Z": 4.0}
    batch_effect = np.asarray([batch_levels[value] for value in batches])
    expression = np.column_stack(
        [7.0 + 5.0 * group_effect - 3.0 * sex_effect + batch_effect]
    )
    design = build_joint_design(
        [
            FactorSpec("group", groups, "categorical"),
            FactorSpec("sex", sexes, "categorical"),
        ],
        [FactorSpec("batch", batches, "categorical")],
    )

    adjusted = remove_nuisance_effects(expression, design)

    for group in ("A", "B"):
        for sex in ("F", "M"):
            cell = adjusted[(groups == group) & (sexes == sex)]
            np.testing.assert_allclose(cell - cell[0], 0.0, atol=1e-12)
    assert design.biology_columns == ("intercept", "group[B]", "sex[M]")
    assert design.nuisance_columns == ("batch[sum:X]", "batch[sum:Y]")
    np.testing.assert_allclose(
        adjusted.mean(axis=0), expression.mean(axis=0), atol=1e-12
    )


def test_joint_design_rejects_complete_confounding() -> None:
    groups = np.asarray(["A", "A", "B", "B", "A", "B"])

    with pytest.raises(CovariateConfoundedError) as exc_info:
        build_joint_design(
            [FactorSpec("group", groups, "categorical")],
            [FactorSpec("batch", groups.copy(), "categorical")],
        )

    _assert_reason(
        exc_info,
        "complete_confounding_or_rank_deficiency",
        ErrorCode.COVARIATE_CONFOUNDED,
    )


def test_joint_design_rejects_rank_deficient_biology() -> None:
    groups = np.asarray(["A", "A", "B", "B", "A", "B"])

    with pytest.raises(CovariateConfoundedError) as exc_info:
        build_joint_design(
            [
                FactorSpec("group", groups, "categorical"),
                FactorSpec("duplicate_group", groups.copy(), "categorical"),
            ],
            [],
        )

    _assert_reason(
        exc_info,
        "biology_rank_deficient",
        ErrorCode.COVARIATE_CONFOUNDED,
    )


@pytest.mark.parametrize(
    "spec",
    [
        FactorSpec("age", [2.0, 2.0, 2.0, 2.0], "numeric"),
        FactorSpec("batch", ["X", "X", "X", "X"], "categorical"),
    ],
)
def test_joint_design_rejects_zero_variance_factors(spec: FactorSpec) -> None:
    with pytest.raises(QCValidationError) as exc_info:
        build_joint_design(
            [FactorSpec("group", ["A", "B", "A", "B"], "categorical")],
            [spec],
        )

    _assert_reason(exc_info, "zero_variance_factor")


def test_joint_design_rejects_nonfinite_numeric_factor() -> None:
    with pytest.raises(QCValidationError) as exc_info:
        build_joint_design(
            [FactorSpec("group", ["A", "B", "A", "B"], "categorical")],
            [FactorSpec("age", [1.0, np.nan, 3.0, 4.0], "numeric")],
        )

    _assert_reason(exc_info, "nonfinite_factor")


@pytest.mark.parametrize(
    "missing_value",
    [None, float("nan"), np.datetime64("NaT"), np.timedelta64("NaT"), np.ma.masked],
)
def test_joint_design_rejects_na_like_categorical_values_before_encoding(
    missing_value: object,
) -> None:
    with pytest.raises(QCValidationError) as exc_info:
        build_joint_design(
            [FactorSpec("group", ["A", missing_value, "B", "B"], "categorical")],
            [],
        )

    _assert_reason(exc_info, "missing_factor_value")


@pytest.mark.parametrize("sentinel_name", ["NA", "NaT"])
def test_joint_design_rejects_pandas_na_instead_of_encoding_ghost_level(
    sentinel_name: str,
) -> None:
    pandas = pytest.importorskip("pandas")
    missing_value = getattr(pandas, sentinel_name)

    with pytest.raises(QCValidationError) as exc_info:
        build_joint_design(
            [FactorSpec("group", ["A", missing_value, "B", "B"], "categorical")],
            [],
        )

    _assert_reason(exc_info, "missing_factor_value")


def test_joint_design_rejects_pandas_na_in_numeric_factor_before_cast() -> None:
    pandas = pytest.importorskip("pandas")

    with pytest.raises(QCValidationError) as exc_info:
        build_joint_design(
            [FactorSpec("group", ["A", "A", "B", "B"], "categorical")],
            [FactorSpec("age", [1.0, pandas.NA, 3.0, 4.0], "numeric")],
        )

    _assert_reason(exc_info, "nonfinite_factor")


def test_joint_design_requires_positive_residual_df() -> None:
    with pytest.raises(CovariateConfoundedError) as exc_info:
        build_joint_design(
            [FactorSpec("group", ["A", "B"], "categorical")],
            [],
        )

    _assert_reason(
        exc_info,
        "nonpositive_residual_df",
        ErrorCode.COVARIATE_CONFOUNDED,
    )


def test_legacy_qc_module_no_longer_imports_scaled_sklearn_pca() -> None:
    source = (PROJECT_ROOT / "legacy" / "modules" / "deseq.py").read_text(
        encoding="utf-8"
    )

    assert "StandardScaler" not in source
    assert "sklearn.decomposition" not in source
    assert "centered_unscaled_pca" in source
    assert "unsupported_qc_transform" in source
    assert "QC_Adjusted_Transformed_Counts.tsv" in source


def test_qc_pca_top_n_is_explicitly_wired_through_both_legacy_callers() -> None:
    legacy_root = PROJECT_ROOT / "legacy"
    config_source = (legacy_root / "config.py").read_text(encoding="utf-8")
    yaml_source = (legacy_root / "workflow_config.yaml").read_text(encoding="utf-8")
    main_source = (legacy_root / "main.py").read_text(encoding="utf-8")
    qc_script_source = (legacy_root / "scripts" / "run_qc.py").read_text(
        encoding="utf-8"
    )
    snakemake_source = (legacy_root / "workflow" / "rules" / "rnaseq.smk").read_text(
        encoding="utf-8"
    )

    assert "'QC_PCA_TOP_N': 500" in config_source
    assert "'QC_BIOLOGY_FACTORS': None" in config_source
    assert "_CONFIG['QC_BIOLOGY_FACTORS'] = [_CONFIG['DESIGN_FACTOR']]" in config_source
    assert "QC_PCA_TOP_N: 500" in yaml_source
    assert "QC_BIOLOGY_FACTORS: null" in yaml_source
    assert "pca_top_n=cfg.QC_PCA_TOP_N" in main_source
    assert "pca_top_n=cfg.QC_PCA_TOP_N" in qc_script_source
    assert "biology_factors=cfg.QC_BIOLOGY_FACTORS" in main_source
    assert "biology_factors=cfg.QC_BIOLOGY_FACTORS" in qc_script_source
    assert "QC_PCA_Method.json" in snakemake_source
    assert "QC_PCA_Selected_Genes.tsv" in snakemake_source
