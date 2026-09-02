"""Executable integration checks around the retained visualization QC wrapper."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import types
from types import SimpleNamespace

import pytest

np = pytest.importorskip(
    "numpy",
    reason="Legacy QC integration tests require the separately reported NumPy lane.",
)
pd = pytest.importorskip(
    "pandas",
    reason="Legacy QC integration tests require the optional pandas runtime.",
)

from rnaseq_downstream.errors import (  # noqa: E402
    CovariateConfoundedError,
    QCValidationError,
)

pytestmark = pytest.mark.unit
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _stub_module(name: str, *, package: bool = False) -> types.ModuleType:
    module = types.ModuleType(name)
    if package:
        module.__path__ = []
    return module


@pytest.fixture
def legacy_deseq(monkeypatch: pytest.MonkeyPatch):
    """Import legacy/modules/deseq.py with plotting/DE packages replaced."""

    pydeseq2 = _stub_module("pydeseq2", package=True)
    pydeseq2_dds = _stub_module("pydeseq2.dds")
    pydeseq2_ds = _stub_module("pydeseq2.ds")
    pydeseq2_dds.DeseqDataSet = object
    pydeseq2_ds.DeseqStats = object

    sklearn = _stub_module("sklearn", package=True)
    sklearn_metrics = _stub_module("sklearn.metrics")
    sklearn_metrics.pairwise_distances = lambda *args, **kwargs: None
    sklearn_metrics.silhouette_score = lambda *args, **kwargs: None

    matplotlib = _stub_module("matplotlib", package=True)
    pyplot = _stub_module("matplotlib.pyplot")
    seaborn = _stub_module("seaborn")

    for name, module in {
        "pydeseq2": pydeseq2,
        "pydeseq2.dds": pydeseq2_dds,
        "pydeseq2.ds": pydeseq2_ds,
        "sklearn": sklearn,
        "sklearn.metrics": sklearn_metrics,
        "matplotlib": matplotlib,
        "matplotlib.pyplot": pyplot,
        "seaborn": seaborn,
    }.items():
        monkeypatch.setitem(sys.modules, name, module)

    spec = importlib.util.spec_from_file_location(
        "checkpoint_a_legacy_deseq",
        PROJECT_ROOT / "legacy" / "modules" / "deseq.py",
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_numeric_coded_batch_is_categorical_unless_explicitly_continuous(
    legacy_deseq,
) -> None:
    groups = np.asarray(["A", "B"] * 6)
    batches = np.asarray([1, 1, 2, 2, 3, 3] * 2)
    batch_effect = {1: -2.0, 2: 4.0, 3: -1.0}
    effects = np.asarray([batch_effect[int(value)] for value in batches])
    transformed = pd.DataFrame(
        {"gene": 10.0 + 3.0 * (groups == "B") + effects},
        index=[f"sample_{index}" for index in range(len(groups))],
    )
    metadata = pd.DataFrame(
        {"group": groups, "batch": batches},
        index=transformed.index,
    )

    categorical = legacy_deseq._residualize_covariates(
        transformed,
        metadata,
        ["batch"],
        ["group"],
        [],
    )
    continuous = legacy_deseq._residualize_covariates(
        transformed,
        metadata,
        ["batch"],
        ["group"],
        ["batch"],
    )

    categorical_means = [
        categorical.loc[metadata["batch"].eq(level), "gene"].mean()
        for level in (1, 2, 3)
    ]
    continuous_means = [
        continuous.loc[metadata["batch"].eq(level), "gene"].mean()
        for level in (1, 2, 3)
    ]
    np.testing.assert_allclose(categorical_means, categorical_means[0], atol=1e-12)
    assert np.ptp(continuous_means) > 1.0


def test_pca_evidence_records_the_fitted_joint_design(
    legacy_deseq,
    tmp_path: Path,
) -> None:
    design = legacy_deseq.build_joint_design(
        [
            legacy_deseq.FactorSpec(
                "group", ["A", "A", "B", "B", "A", "B"], "categorical"
            )
        ],
        [
            legacy_deseq.FactorSpec(
                "batch", ["X", "Y", "X", "Y", "X", "Y"], "categorical"
            )
        ],
    )
    selection = SimpleNamespace(
        names=("gene_1", "gene_2"),
        variances=np.asarray([4.0, 2.0]),
    )

    legacy_deseq._write_pca_evidence(
        tmp_path,
        selection,
        transform="vst",
        requested_top_n=2,
        adjust_factors=["batch"],
        biology_factors=["group"],
        continuous_factors=[],
        adjustment_design=design,
    )

    evidence = json.loads((tmp_path / "QC_PCA_Method.json").read_text(encoding="utf-8"))
    assert evidence["biology_design_columns"] == ["intercept", "group[B]"]
    assert evidence["nuisance_design_columns"] == ["batch[sum:X]"]
    assert evidence["joint_design_rank"] == 3
    assert evidence["joint_design_residual_df"] == 3


@pytest.mark.parametrize("design", ["~ group * sex", "~ group:sex"])
def test_complex_adjustment_design_fails_before_creating_artifacts(
    legacy_deseq,
    tmp_path: Path,
    design: str,
) -> None:
    out_dir = tmp_path / "qc"
    counts = pd.DataFrame(
        {"gene_1": [1, 2, 3, 4], "gene_2": [4, 3, 2, 1]},
        index=["s1", "s2", "s3", "s4"],
    )
    metadata = pd.DataFrame(
        {
            "group": ["A", "A", "B", "B"],
            "sex": ["F", "M", "F", "M"],
            "batch": ["X", "Y", "X", "Y"],
        },
        index=counts.index,
    )

    with pytest.raises(QCValidationError) as captured:
        legacy_deseq.run_qc(
            counts,
            metadata,
            "group",
            out_dir,
            design=design,
            biology_factors=["group", "sex"],
            adjust_factors=["batch"],
        )

    assert captured.value.details["reason"] == "unsupported_adjustment_design"
    assert not out_dir.exists()


def test_complex_design_does_not_block_unadjusted_ordinary_pca(
    legacy_deseq,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    class TransformReached(RuntimeError):
        pass

    def reached_transform(*args, **kwargs):
        raise TransformReached

    monkeypatch.setattr(legacy_deseq, "_get_transformed_counts", reached_transform)
    counts = pd.DataFrame(
        {"gene_1": [1, 2], "gene_2": [2, 1]},
        index=["s1", "s2"],
    )
    metadata = pd.DataFrame(
        {"group": ["A", "B"], "sex": ["F", "M"]},
        index=counts.index,
    )

    with pytest.raises(TransformReached):
        legacy_deseq.run_qc(
            counts,
            metadata,
            "group",
            tmp_path / "qc",
            design="~ group * sex",
            biology_factors=["group", "sex"],
            adjust_factors=[],
        )


def test_stale_adjustment_only_biology_does_not_block_unadjusted_data(
    legacy_deseq,
) -> None:
    transformed = pd.DataFrame(
        {"gene": [1.0, 2.0]},
        index=["s1", "s2"],
    )
    metadata = pd.DataFrame(
        {"group": ["A", "B"]},
        index=transformed.index,
    )

    observed = legacy_deseq._residualize_covariates(
        transformed,
        metadata,
        [],
        ["group", "stale_adjustment_column"],
        [],
    )

    pd.testing.assert_frame_equal(observed, transformed)
    assert observed is not transformed


def test_complete_confounding_fails_before_any_qc_file_is_written(
    legacy_deseq,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    counts = pd.DataFrame(
        {"gene_1": [1, 2, 3, 4], "gene_2": [4, 3, 2, 1]},
        index=["s1", "s2", "s3", "s4"],
    )
    metadata = pd.DataFrame(
        {
            "group": ["A", "A", "B", "B"],
            "batch": ["A", "A", "B", "B"],
        },
        index=counts.index,
    )
    monkeypatch.setattr(
        legacy_deseq,
        "_get_transformed_counts",
        lambda counts, *args, **kwargs: counts.astype(float),
    )
    out_dir = tmp_path / "qc"

    with pytest.raises(CovariateConfoundedError):
        legacy_deseq.run_qc(
            counts,
            metadata,
            "group",
            out_dir,
            design="~ group",
            biology_factors=["group"],
            adjust_factors=["batch"],
            pca_top_n=2,
        )

    assert out_dir.exists()
    assert list(out_dir.iterdir()) == []


@pytest.mark.parametrize(
    ("state", "reason"),
    [
        ("legacy", "legacy_pca_output_conflict"),
        ("wrong_method", "pca_method_conflict"),
        ("missing_selection", "pca_selection_evidence_missing"),
    ],
)
def test_corrected_pca_overwrite_guard_fails_closed(
    legacy_deseq,
    tmp_path: Path,
    state: str,
    reason: str,
) -> None:
    (tmp_path / "PCA_Coordinates.tsv").write_text("legacy\n", encoding="utf-8")
    if state != "legacy":
        method_id = (
            legacy_deseq.PCA_METHOD_ID
            if state == "missing_selection"
            else "legacy-method"
        )
        (tmp_path / "QC_PCA_Method.json").write_text(
            json.dumps({"method_id": method_id}),
            encoding="utf-8",
        )

    with pytest.raises(QCValidationError) as captured:
        legacy_deseq._guard_corrected_pca_outputs(tmp_path)

    assert captured.value.details["reason"] == reason
