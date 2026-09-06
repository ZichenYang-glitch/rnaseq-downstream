"""Source-fidelity regressions for the optional C1 display sidecar."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
import io
import json
import math
from pathlib import Path
import re
from typing import Any

import pytest

np = pytest.importorskip(
    "numpy",
    reason="Display bundle tests require the separately locked NumPy runtime.",
)

from rnaseq_downstream.display_bundle import (  # noqa: E402
    DISPLAY_RENDERER_ID,
    FDR_FLOOR,
    _scatter_svg,
    build_display_bundle,
    verify_display_bundle,
)
from rnaseq_downstream.errors import (  # noqa: E402
    BackendFailedError,
    InputIntegrityError,
    QCValidationError,
)
from rnaseq_downstream.qc_math import (  # noqa: E402
    centered_unscaled_pca,
    select_top_variable_features,
)


pytestmark = pytest.mark.unit

_CORE_FILES = (
    "analysis.json",
    "backend_manifest.json",
    "coefficients.tsv",
    "design.tsv",
    "results.tsv",
)
_SAMPLES = ("sample_1", "sample_2", "sample_3", "sample_4")
_TESTED_GENES = ("gene_1", "gene_2", "gene_3", "gene_4")
_CONTRAST_ID = "treated_vs_control_treat"
_PLAN_ID = "1" * 64
_REQUEST_SHA256 = "2" * 64
_RESULT_HEADER = (
    "gene_id",
    "contrast_id",
    "status",
    "logFC",
    "unshrunk_logFC",
    "logCPM",
    "statistic",
    "statistic_type",
    "statistic_status",
    "PValue",
    "FDR",
    "test_method",
    "lfc_threshold",
)


@dataclass(frozen=True)
class _SyntheticBundle:
    core_dir: Path
    display_dir: Path
    logcpm_path: Path
    analysis: dict[str, Any]
    backend_manifest: dict[str, Any]
    backend_document: dict[str, Any]
    backend_data: dict[str, Any]
    configuration: dict[str, Any]
    core_captured: dict[str, tuple[bytes, str, int]]
    build_artifacts: list[dict[str, Any]]


def _canonical_json(document: object) -> bytes:
    return (
        json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")


def _write_json(path: Path, document: object) -> None:
    path.write_bytes(_canonical_json(document))


def _capture_core(core_dir: Path) -> dict[str, tuple[bytes, str, int]]:
    captured: dict[str, tuple[bytes, str, int]] = {}
    for name in _CORE_FILES:
        payload = (core_dir / name).read_bytes()
        captured[name] = (payload, sha256(payload).hexdigest(), len(payload))
    return captured


def _core_artifacts(
    captured: dict[str, tuple[bytes, str, int]],
) -> list[dict[str, Any]]:
    return [
        {
            "kind": "generated_analysis_artifact",
            "relative_path": name,
            "sha256": captured[name][1],
            "size_bytes": captured[name][2],
        }
        for name in _CORE_FILES
    ]


def _results_payload() -> bytes:
    tested = (
        # The first point has FDR=0 and |logFC| below the glmTreat threshold.
        # It must remain a tested, highlighted point rather than being filtered
        # post hoc by effect size.
        ("gene_1", "0.25", "5.0", "0", "0"),
        ("gene_2", "2.0", "6.0", "0.001", "0.01"),
        ("gene_3", "-1.5", "7.0", "0.1", "0.2"),
        ("gene_4", "0.5", "8.0", "0.9", "1"),
    )
    rows: list[list[str]] = [list(_RESULT_HEADER)]
    for gene_id, logfc, logcpm, p_value, fdr in tested:
        rows.append(
            [
                gene_id,
                _CONTRAST_ID,
                "tested",
                logfc,
                logfc,
                logcpm,
                "",
                "not_reported_by_glmTreat",
                "not_reported",
                p_value,
                fdr,
                "glmTreat",
                "1",
            ]
        )
    rows.append(
        [
            "gene_filtered",
            _CONTRAST_ID,
            "filtered",
            "",
            "",
            "",
            "",
            "not_reported_by_glmTreat",
            "not_applicable_filtered",
            "",
            "",
            "glmTreat",
            "1",
        ]
    )
    buffer = io.StringIO(newline="")
    writer = csv.writer(
        buffer,
        delimiter="\t",
        lineterminator="\n",
        quoting=csv.QUOTE_NONE,
        escapechar=None,
    )
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _make_sources(
    root: Path,
) -> tuple[
    Path,
    Path,
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
    dict[str, tuple[bytes, str, int]],
]:
    core_dir = root / "core"
    core_dir.mkdir(parents=True)
    logcpm_path = root / "private-display-logcpm.tsv"

    contrast = {
        "contrast_id": _CONTRAST_ID,
        "weights": {"conditiontreated": 1.0},
        "lfc_threshold": 1.0,
    }
    analysis = {
        "analysis_request": {"sha256": _REQUEST_SHA256},
        "contrasts": [dict(contrast)],
        "design": {"sample_count": len(_SAMPLES)},
        "genes": {"total": 5, "tested": 4, "filtered": 1},
    }
    backend_manifest = {"input_evidence": {"plan_id": _PLAN_ID}}
    backend_document = {
        "execution_scope": "validated_p0_input",
        "analysis_request": {"sha256": _REQUEST_SHA256},
        "input_evidence": {"plan_id": _PLAN_ID},
        "contrasts": [dict(contrast)],
    }
    logcpm_method = {
        "method": "edgeR::cpm.DGEList",
        "source": "post_filter_post_TMM_observed_DGEList",
        "purpose": "display_only_not_for_inference",
        "arguments": {
            "normalized.lib.sizes": True,
            "log": True,
            "prior.count": 2,
        },
        "scale": "log2",
        "gene_count": len(_TESTED_GENES),
        "sample_count": len(_SAMPLES),
    }
    backend_data = {
        "design_columns": ["(Intercept)", "conditiontreated"],
        "gene_count": 5,
        "tested_gene_count": len(_TESTED_GENES),
        "display_export": logcpm_method,
        "ql_fit_parameters": {
            "abundance.trend": True,
            "robust": True,
            "winsor.tail.p": [0.05, 0.1],
            "legacy": False,
            "top.proportion": None,
            "resolved_top.proportion": 0.4,
            "keep.unit.mat": False,
        },
    }
    configuration = {
        "fdr_threshold": 0.05,
        "pca_top_n": 4,
        "pca_components": [1, 2],
    }

    _write_json(core_dir / "analysis.json", analysis)
    _write_json(core_dir / "backend_manifest.json", backend_manifest)
    (core_dir / "coefficients.tsv").write_text(
        "gene_id\tstatus\tcoefficient\testimate\tscale\n"
        "gene_1\ttested\tconditiontreated\t0.1\tnatural_log\n",
        encoding="utf-8",
    )
    (core_dir / "design.tsv").write_text(
        "sample_id\tcoefficient\tvalue\n"
        + "".join(
            f"{sample}\tconditiontreated\t{index % 2}\n"
            for index, sample in enumerate(_SAMPLES)
        ),
        encoding="utf-8",
    )
    (core_dir / "results.tsv").write_bytes(_results_payload())
    logcpm_path.write_text(
        "gene_id\t" + "\t".join(_SAMPLES) + "\n"
        "gene_1\t1\t2\t3\t4\n"
        "gene_2\t0\t3\t0\t3\n"
        "gene_3\t10\t10\t11\t9\n"
        "gene_4\t5\t8\t2\t7\n",
        encoding="utf-8",
    )
    return (
        core_dir,
        logcpm_path,
        analysis,
        backend_manifest,
        backend_document,
        backend_data,
        configuration,
        _capture_core(core_dir),
    )


def _build(
    root: Path,
    *,
    display_name: str = "display",
    analysis_request_schema_version: str = "1.1",
) -> _SyntheticBundle:
    (
        core_dir,
        logcpm_path,
        analysis,
        backend_manifest,
        backend_document,
        backend_data,
        configuration,
        core_captured,
    ) = _make_sources(root)
    display_dir = root / display_name
    artifacts = build_display_bundle(
        display_dir=display_dir,
        logcpm_path=logcpm_path,
        core_dir=core_dir,
        core_artifacts=_core_artifacts(core_captured),
        backend_document=backend_document,
        backend_data=backend_data,
        configuration=configuration,
        analysis_request_schema_version=analysis_request_schema_version,
    )
    return _SyntheticBundle(
        core_dir=core_dir,
        display_dir=display_dir,
        logcpm_path=logcpm_path,
        analysis=analysis,
        backend_manifest=backend_manifest,
        backend_document=backend_document,
        backend_data=backend_data,
        configuration=configuration,
        core_captured=core_captured,
        build_artifacts=artifacts,
    )


def _display_bytes(display_dir: Path) -> dict[str, bytes]:
    return {
        path.name: path.read_bytes()
        for path in sorted(display_dir.iterdir(), key=lambda item: item.name)
    }


def _manifest(display_dir: Path) -> dict[str, Any]:
    value = json.loads(
        (display_dir / "display_manifest.json").read_text(encoding="utf-8")
    )
    assert isinstance(value, dict)
    return value


def _verify(
    bundle: _SyntheticBundle,
    *,
    sample_ids: tuple[str, ...] = _SAMPLES,
    core_captured: dict[str, tuple[bytes, str, int]] | None = None,
):
    return verify_display_bundle(
        display_dir=bundle.display_dir,
        core_dir=bundle.core_dir,
        core_captured=(
            bundle.core_captured if core_captured is None else core_captured
        ),
        analysis=bundle.analysis,
        backend_manifest=bundle.backend_manifest,
        expected_sample_ids=sample_ids,
    )


def test_build_display_bundle_reproduces_source_values_and_metadata(
    tmp_path: Path,
) -> None:
    bundle = _build(tmp_path)
    manifest = _manifest(bundle.display_dir)

    logcpm = np.asarray(
        [
            [1.0, 2.0, 3.0, 4.0],
            [0.0, 3.0, 0.0, 3.0],
            [10.0, 10.0, 11.0, 9.0],
            [5.0, 8.0, 2.0, 7.0],
        ]
    )
    selection = select_top_variable_features(
        logcpm.T, _TESTED_GENES, bundle.configuration["pca_top_n"]
    )
    expected_pca = centered_unscaled_pca(logcpm.T[:, selection.indices], n_components=2)
    with (bundle.display_dir / "pca_coordinates.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        coordinate_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sample_id"] for row in coordinate_rows] == list(_SAMPLES)
    np.testing.assert_allclose(
        [[float(row["PC1"]), float(row["PC2"])] for row in coordinate_rows],
        expected_pca.coordinates,
        rtol=0.0,
        atol=1e-15,
    )
    with (bundle.display_dir / "pca_selected_genes.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        selected_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["gene_id"] for row in selected_rows] == list(selection.names)

    volcano = (bundle.display_dir / f"volcano--{_CONTRAST_ID}.svg").read_text(
        encoding="utf-8"
    )
    ma_plot = (bundle.display_dir / f"ma--{_CONTRAST_ID}.svg").read_text(
        encoding="utf-8"
    )
    volcano_circles = re.findall(r"<circle\b[^>]*/>", volcano)
    ma_circles = re.findall(r"<circle\b[^>]*/>", ma_plot)
    assert len(volcano_circles) == len(_TESTED_GENES)
    assert len(ma_circles) == len(_TESTED_GENES)
    # Independently hand-calculated from the fixed canvas margins and padded
    # fixture ranges; these do not call the renderer or verifier math.
    assert 'cx="510.000" cy="92.679"' in volcano_circles[0]
    assert 'cx="136.786" cy="315.000"' in ma_circles[0]
    assert 'fill="#c2415d"' in volcano_circles[0]
    assert 'fill="#c2415d"' in ma_circles[0]
    assert sum('fill="#c2415d"' in point for point in volcano_circles) == 2
    assert sum('fill="#c2415d"' in point for point in ma_circles) == 2
    plotted_coordinates = re.findall(r'\b(?:cx|cy)="([^"]+)"', volcano)
    assert plotted_coordinates
    assert all(math.isfinite(float(value)) for value in plotted_coordinates)

    assert manifest["renderer"]["renderer_id"] == DISPLAY_RENDERER_ID
    source_members = [
        {
            "path": name,
            "sha256": bundle.core_captured[name][1],
            "size_bytes": bundle.core_captured[name][2],
        }
        for name in _CORE_FILES
    ]
    assert manifest["source_bundle"]["members"] == source_members
    expected_bundle_id = sha256(
        json.dumps(
            source_members,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()
    assert manifest["source_bundle"] == {
        "source_bundle_id": expected_bundle_id,
        "backend": "edgeR_QL",
        "execution_scope": "validated_p0_input",
        "plan_id": _PLAN_ID,
        "analysis_request_sha256": _REQUEST_SHA256,
        "analysis_request_schema_version": "1.1",
        "members": source_members,
    }
    assert manifest["methods"]["logcpm"]["arguments"] == {
        "normalized.lib.sizes": True,
        "log": True,
        "prior.count": 2,
    }
    assert manifest["methods"]["logcpm"]["source"] == (
        "post_filter_post_TMM_observed_DGEList"
    )
    assert manifest["methods"]["volcano"] == {
        "x": "results.tsv:logFC",
        "y": "-log10(results.tsv:FDR)",
        "fdr_floor_for_rendering": FDR_FLOOR,
        "statistical_recalculation": False,
    }
    plots = {plot["path_id"]: plot for plot in manifest["plots"]}
    volcano_metadata = plots[f"edger_ql_p0_v1.display.volcano.{_CONTRAST_ID}"]
    assert volcano_metadata["gene_count"] == len(_TESTED_GENES)
    assert volcano_metadata["excluded_gene_count"] == 1
    assert volcano_metadata["highlight_thresholds"] == {"FDR": 0.05}
    assert volcano_metadata["source_artifacts"] == [
        {
            "path": "../results.tsv",
            "sha256": bundle.core_captured["results.tsv"][1],
            "size_bytes": bundle.core_captured["results.tsv"][2],
        }
    ]
    assert volcano_metadata["threshold_lines"] == [
        {"axis": "y", "kind": "FDR", "value": 0.05},
        {
            "axis": "x",
            "kind": "absolute_log2_fold_change_test_threshold",
            "values": [-1.0, 1.0],
        },
    ]
    pca_metadata = plots["edger_ql_p0_v1.display.pca"]
    assert pca_metadata["gene_count"] == len(_TESTED_GENES)
    assert pca_metadata["excluded_gene_count"] == 1
    assert pca_metadata["sample_count"] == len(_SAMPLES)
    assert (
        pca_metadata["source_artifacts"][0]["sha256"]
        == sha256(bundle.logcpm_path.read_bytes()).hexdigest()
    )


def test_display_build_is_byte_deterministic(tmp_path: Path) -> None:
    first = _build(tmp_path / "first")
    second = _build(
        tmp_path / "second", analysis_request_schema_version="1.1"
    )

    assert first.build_artifacts == second.build_artifacts
    assert _display_bytes(first.display_dir) == _display_bytes(second.display_dir)


def test_display_bundle_binds_v12_edger_request_schema(tmp_path: Path) -> None:
    bundle = _build(tmp_path, analysis_request_schema_version="1.2")

    assert (
        _manifest(bundle.display_dir)["source_bundle"][
            "analysis_request_schema_version"
        ]
        == "1.2"
    )
    assert _verify(bundle)["summary"]["status"] == "verified_complete"


def test_display_bundle_rejects_unsupported_public_request_schema(
    tmp_path: Path,
) -> None:
    with pytest.raises(
        BackendFailedError, match="request schema identity is incompatible"
    ):
        _build(tmp_path, analysis_request_schema_version="1.0")


def test_verify_display_bundle_accepts_reproducible_bundle(tmp_path: Path) -> None:
    bundle = _build(tmp_path)

    verified = _verify(bundle)

    assert verified["summary"] == {
        "schema_version": "1.0",
        "status": "verified_complete",
        "source_bundle_id": _manifest(bundle.display_dir)["source_bundle"][
            "source_bundle_id"
        ],
        "plot_count": 3,
        "plot_types": {"pca": 1, "volcano": 1, "ma": 1},
        "tested_gene_count": len(_TESTED_GENES),
        "pca_gene_count": len(_TESTED_GENES),
        "sample_count": len(_SAMPLES),
        "configuration": bundle.configuration,
    }
    assert len(verified["artifacts"]) == 7


@pytest.mark.parametrize(
    "mutation",
    ["boolean_as_integer", "integer_as_float", "member_size_as_float"],
)
def test_verify_display_bundle_rejects_json_type_confusion(
    tmp_path: Path,
    mutation: str,
) -> None:
    bundle = _build(tmp_path)
    manifest_path = bundle.display_dir / "display_manifest.json"
    manifest = _manifest(bundle.display_dir)
    if mutation == "boolean_as_integer":
        manifest["methods"]["ql_fit_parameters"]["robust"] = 1
    elif mutation == "integer_as_float":
        manifest["renderer"]["width"] = 960.0
    elif mutation == "member_size_as_float":
        manifest["members"][0]["size_bytes"] = float(
            manifest["members"][0]["size_bytes"]
        )
    else:  # pragma: no cover - keeps additions to the parameter list explicit.
        raise AssertionError(f"unknown mutation: {mutation}")
    manifest_path.chmod(0o600)
    _write_json(manifest_path, manifest)

    with pytest.raises(InputIntegrityError):
        _verify(bundle)


def test_svg_renderer_rejects_unrepresentable_finite_range() -> None:
    with pytest.raises(QCValidationError) as exc_info:
        _scatter_svg(
            title="extreme",
            x_label="x",
            y_label="y",
            x=np.asarray([-1e308, 1e308]),
            y=np.asarray([0.0, 1.0]),
            highlighted=np.asarray([False, False]),
        )

    assert exc_info.value.details["reason"] == "unrepresentable_display_range"


@pytest.mark.parametrize(
    "mutation",
    ["tamper", "extra", "symlink", "missing", "source_hash", "sample_order"],
)
def test_verify_display_bundle_rejects_mutation(
    tmp_path: Path,
    mutation: str,
) -> None:
    bundle = _build(tmp_path)
    sample_ids = _SAMPLES
    core_captured = bundle.core_captured

    if mutation == "tamper":
        target = bundle.display_dir / f"volcano--{_CONTRAST_ID}.svg"
        target.chmod(0o600)
        target.write_bytes(target.read_bytes() + b"<!-- tampered -->\n")
    elif mutation == "extra":
        (bundle.display_dir / "unexpected.txt").write_text(
            "unexpected\n", encoding="utf-8"
        )
    elif mutation == "symlink":
        target = bundle.display_dir / "pca.svg"
        target.unlink()
        target.symlink_to("logcpm.tsv")
    elif mutation == "missing":
        (bundle.display_dir / f"ma--{_CONTRAST_ID}.svg").unlink()
    elif mutation == "source_hash":
        source = bundle.core_dir / "coefficients.tsv"
        source.write_bytes(source.read_bytes() + b"gene_2\ttested\tx\t0\tnatural_log\n")
        core_captured = _capture_core(bundle.core_dir)
    elif mutation == "sample_order":
        sample_ids = tuple(reversed(_SAMPLES))
    else:  # pragma: no cover - keeps additions to the parameter list explicit.
        raise AssertionError(f"unknown mutation: {mutation}")

    with pytest.raises(InputIntegrityError):
        _verify(
            bundle,
            sample_ids=sample_ids,
            core_captured=core_captured,
        )
