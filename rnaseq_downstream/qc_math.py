"""Dependency-light mathematics for experimental QC visualizations.

This module intentionally depends only on NumPy and the toolkit error model.  It
does not implement differential expression and its outputs must never be used as
inputs to a statistical test.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import re
from typing import Iterable, NoReturn, Sequence

import numpy as np

from rnaseq_downstream.errors import CovariateConfoundedError, QCValidationError

PCA_METHOD_ID = "top-positive-variance-centered-unscaled-pca-v1"
_SIMPLE_FORMULA_FACTOR = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]*$")


@dataclass(frozen=True)
class TopFeatureSelection:
    """A deterministic positive-variance feature selection."""

    indices: np.ndarray
    names: tuple[str, ...]
    variances: np.ndarray


@dataclass(frozen=True)
class PCAResult:
    """Centered, unscaled principal-component coordinates and evidence."""

    coordinates: np.ndarray
    components: np.ndarray
    explained_variance_ratio: np.ndarray
    feature_means: np.ndarray


@dataclass(frozen=True)
class FactorSpec:
    """One explicitly typed model factor in sample order."""

    name: str
    values: Sequence[object] | np.ndarray
    kind: str


@dataclass(frozen=True)
class JointDesign:
    """Joint biology and nuisance design matrices."""

    biology: np.ndarray
    nuisance: np.ndarray
    combined: np.ndarray
    biology_columns: tuple[str, ...]
    nuisance_columns: tuple[str, ...]
    rank: int
    residual_df: int


def _validation_error(message: str, *, reason: str, **details: object) -> NoReturn:
    raise QCValidationError(
        message,
        details={"reason": reason, **details},
    )


def _confounding_error(message: str, *, reason: str, **details: object) -> NoReturn:
    raise CovariateConfoundedError(
        message,
        details={"reason": reason, **details},
    )


def _finite_matrix(values: object, *, name: str) -> np.ndarray:
    try:
        matrix = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        _validation_error(
            f"{name} must be a numeric matrix.",
            reason="non_numeric_matrix",
            matrix=name,
            cause_type=type(exc).__name__,
        )
    if matrix.ndim != 2:
        _validation_error(
            f"{name} must be two-dimensional.",
            reason="invalid_matrix_dimensions",
            matrix=name,
            dimensions=int(matrix.ndim),
        )
    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        _validation_error(
            f"{name} must not be empty.",
            reason="empty_matrix",
            matrix=name,
            shape=list(matrix.shape),
        )
    if not np.isfinite(matrix).all():
        _validation_error(
            f"{name} contains non-finite values.",
            reason="nonfinite_values",
            matrix=name,
        )
    return matrix


def select_top_variable_features(
    values: object,
    feature_names: Sequence[object],
    top_n: int,
) -> TopFeatureSelection:
    """Select top positive-variance columns with deterministic tie breaking.

    Variances are sample variances (``ddof=1``).  Ties are ordered by stable
    feature identifier, not by incidental input column order.
    """

    matrix = _finite_matrix(values, name="transformed_expression")
    if matrix.shape[0] < 2:
        _validation_error(
            "PCA feature selection requires at least two samples.",
            reason="insufficient_samples",
            sample_count=int(matrix.shape[0]),
        )
    if (
        isinstance(top_n, bool)
        or not isinstance(top_n, (int, np.integer))
        or top_n <= 0
    ):
        _validation_error(
            "QC PCA top_n must be a positive integer.",
            reason="invalid_top_n",
            top_n=repr(top_n),
        )

    names = tuple(str(name) for name in feature_names)
    if len(names) != matrix.shape[1]:
        _validation_error(
            "Feature identifiers do not match the expression columns.",
            reason="feature_count_mismatch",
            feature_count=len(names),
            matrix_columns=int(matrix.shape[1]),
        )
    if any(not name for name in names):
        _validation_error(
            "Feature identifiers must be non-empty.",
            reason="missing_feature_id",
        )
    duplicate_names = sorted(
        name for name, count in Counter(names).items() if count > 1
    )
    if duplicate_names:
        _validation_error(
            "Feature identifiers must be unique before PCA selection.",
            reason="duplicate_feature_ids",
            duplicate_feature_ids=duplicate_names,
        )

    variances = np.var(matrix, axis=0, ddof=1)
    if not np.isfinite(variances).all():
        _validation_error(
            "Feature variance calculation produced non-finite values.",
            reason="nonfinite_feature_variance",
        )
    positive = [index for index, variance in enumerate(variances) if variance > 0.0]
    positive.sort(key=lambda index: (-float(variances[index]), names[index]))
    selected = positive[: min(int(top_n), len(positive))]
    if not selected:
        _validation_error(
            "No positive-variance genes are available for PCA.",
            reason="no_positive_variance_features",
        )

    indices = np.asarray(selected, dtype=int)
    return TopFeatureSelection(
        indices=indices,
        names=tuple(names[index] for index in selected),
        variances=variances[indices].copy(),
    )


def centered_unscaled_pca(values: object, n_components: int = 2) -> PCAResult:
    """Run PCA by SVD after centering features without variance scaling."""

    matrix = _finite_matrix(values, name="pca_expression")
    if (
        isinstance(n_components, bool)
        or not isinstance(n_components, (int, np.integer))
        or n_components <= 0
    ):
        _validation_error(
            "n_components must be a positive integer.",
            reason="invalid_component_count",
            n_components=repr(n_components),
        )
    # Centering removes at most one sample-space direction. Feature space is
    # not reduced merely because there are fewer features than samples.
    maximum = min(matrix.shape[0] - 1, matrix.shape[1])
    if n_components > maximum:
        _validation_error(
            "The requested PCA dimensions exceed the centered-matrix limit.",
            reason="insufficient_pca_dimensions",
            n_components=int(n_components),
            shape=list(matrix.shape),
            maximum_components=int(maximum),
        )

    means = matrix.mean(axis=0)
    centered = matrix - means
    try:
        left, singular_values, right = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError as exc:
        _validation_error(
            "PCA singular-value decomposition did not converge.",
            reason="pca_svd_failed",
            cause_type=type(exc).__name__,
        )
    total = float(np.square(singular_values).sum())
    if total <= 0.0:
        _validation_error(
            "PCA requires at least one positive-variance feature.",
            reason="zero_total_variance",
        )

    components = right[:n_components].copy()
    coordinates = left[:, :n_components] * singular_values[:n_components]

    # SVD signs are arbitrary.  Make the loading with greatest absolute value
    # positive; np.argmax makes ties deterministic by choosing the first feature.
    for component_index in range(int(n_components)):
        pivot = int(np.argmax(np.abs(components[component_index])))
        if components[component_index, pivot] < 0.0:
            components[component_index] *= -1.0
            coordinates[:, component_index] *= -1.0

    explained = np.square(singular_values[:n_components]) / total
    return PCAResult(
        coordinates=coordinates,
        components=components,
        explained_variance_ratio=explained,
        feature_means=means,
    )


def metadata_reorder_indices(
    expression_sample_ids: Sequence[object],
    metadata_sample_ids: Sequence[object],
) -> np.ndarray:
    """Return the lossless metadata reorder needed for expression sample order."""

    expression_ids = tuple(expression_sample_ids)
    metadata_ids = tuple(metadata_sample_ids)
    if not expression_ids or not metadata_ids:
        _validation_error(
            "Expression and metadata sample identifiers must not be empty.",
            reason="empty_sample_ids",
        )
    try:
        expression_unique = set(expression_ids)
        metadata_unique = set(metadata_ids)
    except TypeError:
        _validation_error(
            "Sample identifiers must be hashable scalar values.",
            reason="invalid_sample_ids",
        )
    if len(expression_unique) != len(expression_ids):
        _validation_error(
            "Expression sample identifiers contain duplicates.",
            reason="duplicate_expression_samples",
        )
    if len(metadata_unique) != len(metadata_ids):
        _validation_error(
            "Metadata sample identifiers contain duplicates.",
            reason="duplicate_metadata_samples",
        )
    missing_metadata = expression_unique - metadata_unique
    extra_metadata = metadata_unique - expression_unique
    if missing_metadata or extra_metadata:
        _validation_error(
            "Expression and metadata sample identifiers must match exactly.",
            reason="sample_mismatch",
            missing_metadata=sorted(str(value) for value in missing_metadata),
            extra_metadata=sorted(str(value) for value in extra_metadata),
        )

    positions = {sample_id: index for index, sample_id in enumerate(metadata_ids)}
    return np.asarray([positions[sample_id] for sample_id in expression_ids], dtype=int)


def validate_factor_names(
    available_factors: Iterable[str],
    biology_factors: Sequence[str],
    nuisance_factors: Sequence[str],
) -> None:
    """Validate factor identity before any categorical or numeric encoding."""

    available = set(available_factors)
    biology = tuple(biology_factors)
    nuisance = tuple(nuisance_factors)
    if not biology:
        _validation_error(
            "At least one biology factor is required for adjusted PCA.",
            reason="missing_biology_factor",
        )
    for role, names in (("biology", biology), ("nuisance", nuisance)):
        duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
        if duplicates:
            _validation_error(
                f"{role.capitalize()} factors contain duplicates.",
                reason="duplicate_factors",
                role=role,
                duplicate_factors=duplicates,
            )
        unknown = sorted(set(names) - available)
        if unknown:
            _validation_error(
                f"Unknown {role} factors were requested.",
                reason="unknown_factors",
                role=role,
                unknown_factors=unknown,
            )
    overlap = sorted(set(biology) & set(nuisance))
    if overlap:
        _validation_error(
            "Biology and nuisance factors must not overlap.",
            reason="overlapping_factors",
            overlapping_factors=overlap,
        )


def resolve_adjustment_biology_factors(
    design: object,
    design_factor: object,
    biology_factors: Sequence[object] | object | None,
    nuisance_factors: Sequence[object] | object | None,
) -> tuple[str, ...]:
    """Resolve a safe additive biology model for nuisance-adjusted QC.

    The checkpoint-A adjustment encoder supports only simple additive factor
    names. The primary design factor is always protected, while every other
    non-nuisance term in the DE design must be named explicitly by the caller.
    """

    primary = str(design_factor)
    if not primary or not _SIMPLE_FORMULA_FACTOR.fullmatch(primary):
        _validation_error(
            "The primary design factor must be a simple non-empty factor name.",
            reason="invalid_design_factor",
            design_factor=primary,
        )

    def normalize_names(values: Sequence[object] | object | None, *, role: str):
        if values is None:
            return []
        raw = [values] if isinstance(values, str) else list(values)
        names = [str(value) for value in raw]
        invalid = [name for name in names if not _SIMPLE_FORMULA_FACTOR.fullmatch(name)]
        if invalid:
            _validation_error(
                f"{role.capitalize()} factor names must use simple identifiers.",
                reason="invalid_factor_name",
                role=role,
                invalid_factors=invalid,
            )
        duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
        if duplicates:
            _validation_error(
                f"{role.capitalize()} factors contain duplicates.",
                reason="duplicate_factors",
                role=role,
                duplicate_factors=duplicates,
            )
        return names

    explicit_biology = normalize_names(biology_factors, role="biology")
    nuisance = normalize_names(nuisance_factors, role="nuisance")
    if primary in nuisance:
        _validation_error(
            "The primary design factor cannot be removed as a nuisance effect.",
            reason="primary_design_factor_is_nuisance",
            design_factor=primary,
        )

    formula = f"~ {primary}" if design is None else str(design).strip()
    if not formula.startswith("~") or formula.count("~") != 1:
        _validation_error(
            "Adjusted PCA requires a simple additive design formula.",
            reason="unsupported_adjustment_design",
            design=formula,
            supported_example=f"~ {primary} + covariate",
        )
    right_hand_side = formula[1:].strip()
    raw_terms = [term.strip() for term in right_hand_side.split("+")]
    if not right_hand_side or any(not term for term in raw_terms):
        _validation_error(
            "Adjusted PCA requires a non-empty simple additive design formula.",
            reason="unsupported_adjustment_design",
            design=formula,
        )

    design_terms: list[str] = []
    for term in raw_terms:
        if term in {"0", "1"}:
            _validation_error(
                "Adjusted PCA does not support explicit intercept-control syntax.",
                reason="unsupported_adjustment_design",
                design=formula,
                unsupported_term=term,
                supported_syntax="simple additive factor names with an implicit intercept",
            )
        if not _SIMPLE_FORMULA_FACTOR.fullmatch(term):
            _validation_error(
                "Adjusted PCA does not support interactions or complex formula terms.",
                reason="unsupported_adjustment_design",
                design=formula,
                unsupported_term=term,
                supported_syntax="simple additive factor names separated by '+'",
            )
        if term not in design_terms:
            design_terms.append(term)

    if primary not in design_terms:
        _validation_error(
            "The adjusted-PCA design must include the primary design factor.",
            reason="primary_design_factor_missing",
            design=formula,
            design_factor=primary,
        )

    protected = [primary]
    for factor in explicit_biology:
        if factor not in protected:
            protected.append(factor)
    overlap = sorted(set(protected) & set(nuisance))
    if overlap:
        _validation_error(
            "Biology and nuisance factors must not overlap.",
            reason="overlapping_factors",
            overlapping_factors=overlap,
        )

    required = [factor for factor in design_terms if factor not in nuisance]
    missing = [factor for factor in required if factor not in protected]
    if missing:
        _validation_error(
            "Adjusted PCA requires explicit protection of every biology design term.",
            reason="incomplete_biology_protection",
            design=formula,
            design_factor=primary,
            missing_biology_factors=missing,
            protected_biology_factors=protected,
            nuisance_factors=nuisance,
        )
    return tuple(protected)


def _factor_matrix(
    spec: FactorSpec,
    sample_count: int,
    *,
    categorical_coding: str,
) -> tuple[np.ndarray, tuple[str, ...]]:
    name = str(spec.name)
    if not name:
        _validation_error(
            "Factor names must be non-empty.", reason="missing_factor_name"
        )
    # Use an object view for the boundary check so nullable scalar sentinels
    # such as pandas.NA are not stringified into an apparently valid level by
    # NumPy's common-dtype inference.
    values = np.asarray(spec.values, dtype=object)
    if values.ndim != 1 or len(values) != sample_count:
        _validation_error(
            f"Factor '{name}' does not align with the samples.",
            reason="factor_length_mismatch",
            factor=name,
            factor_length=int(values.size),
            sample_count=int(sample_count),
        )

    for value in values.tolist():
        if _is_na_like(value):
            reason = (
                "nonfinite_factor" if spec.kind == "numeric" else "missing_factor_value"
            )
            _validation_error(
                f"Factor '{name}' contains missing values.",
                reason=reason,
                factor=name,
            )

    if spec.kind == "numeric":
        try:
            numeric = values.astype(float)
        except (TypeError, ValueError) as exc:
            _validation_error(
                f"Numeric factor '{name}' contains non-numeric values.",
                reason="invalid_numeric_factor",
                factor=name,
                cause_type=type(exc).__name__,
            )
        if not np.isfinite(numeric).all():
            _validation_error(
                f"Numeric factor '{name}' contains missing or non-finite values.",
                reason="nonfinite_factor",
                factor=name,
            )
        if np.all(numeric == numeric[0]):
            _validation_error(
                f"Numeric factor '{name}' has zero variance.",
                reason="zero_variance_factor",
                factor=name,
            )
        centered = numeric - numeric.mean()
        return centered[:, np.newaxis], (name,)

    if spec.kind != "categorical":
        _validation_error(
            f"Factor '{name}' has an unsupported kind.",
            reason="invalid_factor_kind",
            factor=name,
            kind=str(spec.kind),
        )

    labels: list[str] = []
    for value in values.tolist():
        if isinstance(value, str) and not value.strip():
            _validation_error(
                f"Categorical factor '{name}' contains blank values.",
                reason="missing_factor_value",
                factor=name,
            )
        try:
            is_nonfinite_number = (
                np.isscalar(value)
                and np.asarray(value).dtype.kind in "fc"
                and not np.isfinite(value)
            )
            if bool(is_nonfinite_number):
                _validation_error(
                    f"Categorical factor '{name}' contains missing values.",
                    reason="missing_factor_value",
                    factor=name,
                )
        except (TypeError, ValueError):
            _validation_error(
                f"Categorical factor '{name}' contains invalid values.",
                reason="invalid_factor_value",
                factor=name,
            )
        labels.append(str(value))

    levels = sorted(set(labels))
    if len(levels) < 2:
        _validation_error(
            f"Categorical factor '{name}' has zero variance.",
            reason="zero_variance_factor",
            factor=name,
        )
    label_array = np.asarray(labels)
    if categorical_coding == "treatment":
        encoded = np.column_stack(
            [label_array == level for level in levels[1:]]
        ).astype(float)
        columns = tuple(f"{name}[{level}]" for level in levels[1:])
    elif categorical_coding == "sum":
        reference = levels[-1]
        encoded = np.column_stack(
            [
                np.where(
                    label_array == level,
                    1.0,
                    np.where(label_array == reference, -1.0, 0.0),
                )
                for level in levels[:-1]
            ]
        )
        columns = tuple(f"{name}[sum:{level}]" for level in levels[:-1])
    else:
        _validation_error(
            "Categorical factor coding must be treatment or sum.",
            reason="invalid_categorical_coding",
            factor=name,
            categorical_coding=categorical_coding,
        )
    return encoded, columns


def _is_na_like(value: object) -> bool:
    """Return true for scalar missing-value sentinels without importing pandas."""

    if value is None or np.ma.is_masked(value):
        return True

    value_type = type(value)
    if value_type.__module__.startswith("pandas.") and value_type.__name__ in {
        "NAType",
        "NaTType",
    }:
        return True

    try:
        scalar = np.asarray(value)
    except (TypeError, ValueError):
        return False
    if scalar.ndim != 0:
        return False
    if scalar.dtype.kind in "fc":
        return bool(np.isnan(scalar))
    if scalar.dtype.kind in "mM":
        return bool(np.isnat(scalar))
    if scalar.dtype.kind != "O":
        return False

    # Decimal NaN and similar scalar sentinels are non-reflexive. Keep this
    # deliberately scalar-only so array-valued objects are rejected by the
    # ordinary factor-shape/value checks rather than treated as missing.
    try:
        unequal = value != value
        if isinstance(unequal, (bool, np.bool_)):
            return bool(unequal)
    except (TypeError, ValueError):
        return False
    return False


def _rank(matrix: np.ndarray) -> int:
    try:
        singular_values = np.linalg.svd(matrix, compute_uv=False)
    except np.linalg.LinAlgError as exc:
        _validation_error(
            "Joint-design rank calculation did not converge.",
            reason="rank_calculation_failed",
            cause_type=type(exc).__name__,
        )
    if singular_values.size == 0:
        return 0
    tolerance = (
        max(matrix.shape)
        * np.finfo(singular_values.dtype).eps
        * float(singular_values[0])
    )
    return int(np.count_nonzero(singular_values > tolerance))


def build_joint_design(
    biology_factors: Sequence[FactorSpec],
    nuisance_factors: Sequence[FactorSpec],
) -> JointDesign:
    """Encode and validate ``[intercept, biology, nuisance]`` jointly."""

    if not biology_factors:
        _validation_error(
            "At least one biology factor is required for adjusted PCA.",
            reason="missing_biology_factor",
        )
    biology_names = tuple(str(spec.name) for spec in biology_factors)
    nuisance_names = tuple(str(spec.name) for spec in nuisance_factors)
    validate_factor_names(
        set(biology_names) | set(nuisance_names), biology_names, nuisance_names
    )

    first_values = np.asarray(biology_factors[0].values)
    if first_values.ndim != 1 or len(first_values) == 0:
        _validation_error(
            "Model factors must contain at least one sample.",
            reason="empty_factor",
            factor=biology_names[0],
        )
    sample_count = len(first_values)
    biology_parts = [np.ones((sample_count, 1), dtype=float)]
    biology_columns = ["intercept"]
    for spec in biology_factors:
        encoded, columns = _factor_matrix(
            spec,
            sample_count,
            categorical_coding="treatment",
        )
        biology_parts.append(encoded)
        biology_columns.extend(columns)

    nuisance_parts: list[np.ndarray] = []
    nuisance_columns: list[str] = []
    for spec in nuisance_factors:
        encoded, columns = _factor_matrix(
            spec,
            sample_count,
            categorical_coding="sum",
        )
        nuisance_parts.append(encoded)
        nuisance_columns.extend(columns)

    biology = np.column_stack(biology_parts)
    nuisance = (
        np.column_stack(nuisance_parts)
        if nuisance_parts
        else np.empty((sample_count, 0), dtype=float)
    )
    combined = np.column_stack([biology, nuisance])
    biology_rank = _rank(biology)
    combined_rank = _rank(combined)
    if biology_rank < biology.shape[1]:
        _confounding_error(
            "The biology design is rank deficient.",
            reason="biology_rank_deficient",
            rank=biology_rank,
            columns=int(biology.shape[1]),
        )
    if combined_rank < combined.shape[1]:
        _confounding_error(
            "Biology and nuisance effects are rank deficient or completely confounded.",
            reason="complete_confounding_or_rank_deficiency",
            rank=combined_rank,
            columns=int(combined.shape[1]),
            biology_factors=list(biology_names),
            nuisance_factors=list(nuisance_names),
        )
    residual_df = sample_count - combined_rank
    if residual_df <= 0:
        _confounding_error(
            "The joint adjustment model has no residual degrees of freedom.",
            reason="nonpositive_residual_df",
            sample_count=int(sample_count),
            rank=combined_rank,
            residual_df=int(residual_df),
        )

    return JointDesign(
        biology=biology,
        nuisance=nuisance,
        combined=combined,
        biology_columns=tuple(biology_columns),
        nuisance_columns=tuple(nuisance_columns),
        rank=combined_rank,
        residual_df=residual_df,
    )


def remove_nuisance_effects(values: object, design: JointDesign) -> np.ndarray:
    """Fit the joint model and subtract only the fitted nuisance term."""

    expression = _finite_matrix(values, name="transformed_expression")
    if expression.shape[0] != design.combined.shape[0]:
        _validation_error(
            "The joint design does not align with the expression samples.",
            reason="design_sample_mismatch",
            expression_samples=int(expression.shape[0]),
            design_samples=int(design.combined.shape[0]),
        )
    if design.nuisance.shape[1] == 0:
        return expression.copy()

    try:
        coefficients, _, fitted_rank, _ = np.linalg.lstsq(
            design.combined, expression, rcond=None
        )
    except np.linalg.LinAlgError as exc:
        _validation_error(
            "The joint adjustment model could not be fitted.",
            reason="joint_model_fit_failed",
            cause_type=type(exc).__name__,
        )
    if int(fitted_rank) != design.rank:
        _validation_error(
            "The joint adjustment model lost rank during fitting.",
            reason="fit_rank_mismatch",
            expected_rank=int(design.rank),
            fitted_rank=int(fitted_rank),
        )
    nuisance_start = design.biology.shape[1]
    nuisance_coefficients = coefficients[nuisance_start:, :]
    return expression - design.nuisance @ nuisance_coefficients
