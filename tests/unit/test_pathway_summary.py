"""Fail-closed verification for schema 1.1 pathway evidence bundles."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
from typing import Any

import pytest

from rnaseq_downstream import display_bundle
from rnaseq_downstream.errors import InputIntegrityError
from rnaseq_downstream.run_summary import summarize_run

RUNTIME = {
    "R": "4.6.1",
    "Bioconductor": "3.23",
    "BiocVersion_package": "3.23.1",
    "edgeR": "4.10.0",
    "tximport": "1.40.0",
    "limma": "3.68.0",
}
PLAN_ID = "a" * 64
PATHWAY_HEADER = [
    "contrast_id",
    "gene_set_id",
    "gene_set_description",
    "method_id",
    "test_class",
    "hypothesis",
    "inference_role",
    "status",
    "status_reason",
    "direction",
    "proportion_down",
    "proportion_up",
    "p_value",
    "fdr",
    "fdr_family_id",
    "gmt_member_count_raw",
    "gmt_symbol_count_unique",
    "mapped_symbol_count_unique",
    "ambiguous_symbol_count_unique",
    "unmapped_symbol_count_unique",
    "mapping_rate",
    "mapped_gene_id_count_unique",
    "tested_gene_count",
    "filtered_gene_count",
    "tested_universe_gene_count",
    "method_ngenes",
    "correlation_status",
    "correlation_estimate_raw",
    "correlation_effective",
    "vif_used",
    "rotation_status",
]


def _input_evidence() -> dict[str, object]:
    return {
        "kind": "validated_input_bundle",
        "plan_id": PLAN_ID,
        "bundle_path": "/captured/input-bundle",
        "validation_bundle_artifacts": [
            {
                "role": "bundle_manifest",
                "path": "/captured/input-bundle/bundle_manifest.json",
                "sha256": "c" * 64,
                "size_bytes": 1,
            }
        ],
        "r_input_snapshots": [
            {
                "role": "metadata",
                "source_path": "/captured/metadata.tsv",
                "sha256": "d" * 64,
                "size_bytes": 1,
            },
            {
                "role": "pathways.gmt",
                "source_path": "/captured/sets.gmt",
                "sha256": "e" * 64,
                "size_bytes": 123,
            },
            {
                "role": "pathways.annotation",
                "source_path": "/captured/annotation.tsv",
                "sha256": "f" * 64,
                "size_bytes": 456,
            },
        ],
        "digest_coupling": "private_copy_and_hash_same_source_stream",
    }


def _write_json(path: Path, document: object) -> None:
    path.write_text(
        json.dumps(document, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    stream = io.StringIO(newline="")
    writer = csv.writer(
        stream,
        delimiter="\t",
        lineterminator="\n",
        quoting=csv.QUOTE_NONE,
        escapechar=None,
    )
    writer.writerow(header)
    writer.writerows(rows)
    path.write_text(stream.getvalue(), encoding="utf-8")


def _member(path: Path) -> dict[str, object]:
    payload = path.read_bytes()
    return {
        "path": path.name,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


def _refresh_manifest(run_dir: Path) -> None:
    manifest = {
        "schema_version": "1.1",
        "kind": "edger_ql_backend_manifest",
        "backend": "edgeR_QL",
        "runtime_identity": RUNTIME,
        "execution_scope": "validated_p0_input",
        "input_evidence": _input_evidence(),
        "members": [
            _member(run_dir / name)
            for name in (
                "analysis.json",
                "coefficients.tsv",
                "design.tsv",
                "pathway_results.tsv",
                "results.tsv",
            )
        ],
    }
    _write_json(run_dir / "backend_manifest.json", manifest)


def _set_record(
    gene_set_id: str,
    description: str,
    *,
    raw: int,
    mapped_symbols: int,
    unmapped_symbols: int,
    mapped_genes: int,
    tested: int,
    filtered: int,
) -> dict[str, Any]:
    return {
        "gene_set_id": gene_set_id,
        "gene_set_description": description,
        "gmt_member_count_raw": raw,
        "gmt_symbol_count_unique": raw,
        "mapped_symbol_count_unique": mapped_symbols,
        "ambiguous_symbol_count_unique": 0,
        "unmapped_symbol_count_unique": unmapped_symbols,
        "mapping_rate": mapped_symbols / raw,
        "mapped_gene_id_count_unique": mapped_genes,
        "tested_gene_count": tested,
        "filtered_gene_count": filtered,
        "absent_from_assay_gene_count": mapped_genes - tested - filtered,
    }


def _pathway_analysis(sets: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "gene_sets": {
            "gmt": {
                "collection": "fixture-collection",
                "version": "1.0",
                "identifier_type": "symbol",
                "sha256": "e" * 64,
                "size_bytes": 123,
                "gene_set_count": len(sets),
            },
            "annotation": {
                "name": "fixture-annotation",
                "version": "1.0",
                "gene_id_column": "gene_id",
                "symbol_column": "symbol",
                "sha256": "f" * 64,
                "size_bytes": 456,
                "row_count": 10,
            },
            "minimum_tested_genes": 2,
            "sets": sets,
        },
        "mapping_policy": {
            "source_identifier": "symbol",
            "target_identifier": "stable_gene_id",
            "annotation_gene_id_version_stripping": False,
            "one_to_many_symbols": "ambiguous_excluded",
            "duplicate_gmt_members": "hard_fail",
            "mapping_rate": (
                "uniquely_mapped_unique_symbols divided by unique_GMT_symbols"
            ),
            "tested_membership": "intersection_with_filterByExpr_tested_universe",
            "not_tested_policy": "tested_gene_count_below_minimum",
        },
        "tested_universe_gene_count": 4,
        "methods": {
            "limma_mroast": {
                "generic": "limma::mroast",
                "dispatch": "edgeR::mroast.DGEGLM",
                "test_class": "self_contained",
                "inference_role": "corroborative",
                "statistical_null": "zero_effect",
                "parameters": {
                    "set.statistic": "mean",
                    "gene.weights": None,
                    "nrot": 9999,
                    "adjust.method": "BH",
                    "midp": False,
                    "sort": "none",
                },
            },
            "limma_fry": {
                "generic": "limma::fry",
                "dispatch": "edgeR::fry.DGEGLM",
                "test_class": "self_contained",
                "inference_role": "primary",
                "statistical_null": "zero_effect",
                "parameters": {"sort": "none"},
            },
            "limma_camera": {
                "generic": "limma::camera",
                "dispatch": "edgeR::camera.DGEGLM",
                "test_class": "competitive",
                "inference_role": "supplementary",
                "statistical_null": "no_enrichment_relative_to_complement",
                "parameters": {
                    "weights": None,
                    "use.ranks": False,
                    "allow.neg.cor": False,
                    "inter.gene.cor": "estimated_per_gene_set",
                    "sort": False,
                },
            },
        },
        "multiple_testing": {
            "method": "Benjamini-Hochberg",
            "scope": (
                "separately within contrast x method x hypothesis across tested "
                "gene sets only"
            ),
            "family_id_format": "contrast_id|method_id|hypothesis",
            "python_independent_recalculation_required": True,
        },
        "contrasts": [
            {
                "contrast_id": "treated_vs_control",
                "gene_level_lfc_threshold": 0,
                "pathway_statistical_null": "zero_effect",
                "gene_level_lfc_threshold_applied_to_pathways": False,
                "ordered_set_lists": {
                    "self_contained": {
                        "gene_set_ids": ["set_a", "set_b"],
                        "sha256": "1" * 64,
                    },
                    "competitive": {
                        "gene_set_ids": ["set_a", "set_b"],
                        "sha256": "2" * 64,
                    },
                },
                "rotation": {
                    "method_id": "limma_mroast",
                    "seed": 1729,
                    "seed_policy": (
                        "base_1729_plus_zero_based_declared_contrast_index"
                    ),
                    "rng": {
                        "kind": "Mersenne-Twister",
                        "normal.kind": "Inversion",
                        "sample.kind": "Rejection",
                    },
                    "nrot": 9999,
                    "reset_before_each_contrast": True,
                },
            }
        ],
    }


def _analysis(sets: list[dict[str, Any]]) -> dict[str, Any]:
    evidence = _input_evidence()
    return {
        "schema_version": "1.1",
        "kind": "edger_ql_analysis",
        "backend": "edgeR_QL",
        "execution_scope": "validated_p0_input",
        "analysis_request": {"path": "/captured/request.json", "sha256": "b" * 64},
        "input_evidence": evidence,
        "runtime_identity": RUNTIME,
        "input_semantics": "featurecounts_integer",
        "route_observed": {
            "constructor": "edgeR::DGEList",
            "count_semantics": "integer",
            "transcript_length_offset": False,
        },
        "pipeline": [
            {"step": "filterByExpr", "arguments": {"design": True}},
            {"step": "normLibSizes", "arguments": {"method": "TMM"}},
            {"step": "glmQLFit", "arguments": {"legacy": False, "robust": True}},
            {
                "step": "contrast_test",
                "dispatch": (
                    "lfc_threshold == 0: glmQLFTest; lfc_threshold > 0: glmTreat"
                ),
            },
        ],
        "design": {
            "intercept": True,
            "terms": "condition",
            "variable_types": {"condition": "factor"},
            "factor_levels": {"condition": ["control", "treated"]},
            "columns": ["(Intercept)", "conditiontreated"],
            "sample_count": 4,
            "rank": 2,
            "residual_df": 2,
            "qr_tolerance": 1e-7,
        },
        "contrasts": [
            {
                "contrast_id": "treated_vs_control",
                "weights": {"conditiontreated": 1},
                "lfc_threshold": 0,
                "test_method": "glmQLFTest",
                "treat_null": None,
                "estimability_residual": 0,
                "estimability_tolerance": 1e-7,
            }
        ],
        "genes": {"total": 5, "tested": 4, "filtered": 1},
        "status_vocabulary": [
            "filtered",
            "not_tested",
            "not_estimable",
            "failed",
            "tested",
        ],
        "result_logFC_scale": "log2",
        "coefficient_scale": "natural_log",
        "multiple_testing": "Benjamini-Hochberg within each contrast",
        "pathway_analysis": _pathway_analysis(sets),
    }


def _pathway_row(
    record: dict[str, Any],
    method_id: str,
    hypothesis: str,
    *,
    p_value: float | None,
    fdr: float | None,
    direction: str = "Up",
    correlation: float | None = None,
) -> list[str]:
    method_properties = {
        "limma_mroast": ("self_contained", "corroborative"),
        "limma_fry": ("self_contained", "primary"),
        "limma_camera": ("competitive", "supplementary"),
    }
    test_class, inference_role = method_properties[method_id]
    tested = p_value is not None
    if tested:
        status = "tested"
        reason = ""
    else:
        status = "not_tested"
        reason = "tested_gene_count_below_minimum"
        direction = ""
    is_mroast = method_id == "limma_mroast"
    is_camera = method_id == "limma_camera"
    if tested and hypothesis == "mixed":
        direction = "Mixed"
    if is_camera and tested:
        assert correlation is not None
        correlation_status = "estimated_set_specific"
        correlation_raw = str(correlation)
        correlation_effective = str(max(0.0, correlation))
        vif = str(max(1.0, 1 + (record["tested_gene_count"] - 1) * correlation))
    elif is_camera:
        correlation_status = "not_estimated_not_tested"
        correlation_raw = correlation_effective = vif = ""
    else:
        correlation_status = "not_applicable"
        correlation_raw = correlation_effective = vif = ""
    rotation = {
        "limma_mroast": (
            "performed_fixed_seed_9999_rotations"
            if tested
            else "not_performed_not_tested"
        ),
        "limma_fry": "not_applicable_analytic_approximation",
        "limma_camera": "not_applicable",
    }[method_id]
    return [
        "treated_vs_control",
        record["gene_set_id"],
        record["gene_set_description"],
        method_id,
        test_class,
        hypothesis,
        inference_role,
        status,
        reason,
        direction,
        "0.25" if tested and is_mroast else "",
        "0.75" if tested and is_mroast else "",
        "" if p_value is None else str(p_value),
        "" if fdr is None else str(fdr),
        f"treated_vs_control|{method_id}|{hypothesis}",
        str(record["gmt_member_count_raw"]),
        str(record["gmt_symbol_count_unique"]),
        str(record["mapped_symbol_count_unique"]),
        str(record["ambiguous_symbol_count_unique"]),
        str(record["unmapped_symbol_count_unique"]),
        str(record["mapping_rate"]),
        str(record["mapped_gene_id_count_unique"]),
        str(record["tested_gene_count"]),
        str(record["filtered_gene_count"]),
        "4",
        str(record["tested_gene_count"]) if tested else "",
        correlation_status,
        correlation_raw,
        correlation_effective,
        vif,
        rotation,
    ]


def _pathway_rows(sets: list[dict[str, Any]]) -> list[list[str]]:
    p_values = {
        "limma_mroast": {
            "directional": ((0.01, 0.02), (0.2, 0.2)),
            "mixed": ((0.04, 0.08), (0.5, 0.5)),
        },
        "limma_fry": {
            "directional": ((0.03, 0.06), (0.4, 0.4)),
            "mixed": ((0.001, 0.002), (0.8, 0.8)),
        },
        "limma_camera": {
            "directional": ((0.02, 0.04), (0.3, 0.3)),
        },
    }
    rows: list[list[str]] = []
    for method_id, hypotheses in (
        ("limma_camera", ("directional",)),
        ("limma_fry", ("directional", "mixed")),
        ("limma_mroast", ("directional", "mixed")),
    ):
        for record in sets:
            for hypothesis in hypotheses:
                if record["gene_set_id"] == "set_low":
                    p_value = fdr = None
                else:
                    set_index = 0 if record["gene_set_id"] == "set_a" else 1
                    p_value, fdr = p_values[method_id][hypothesis][set_index]
                correlation = None
                if method_id == "limma_camera" and p_value is not None:
                    correlation = -0.2 if record["gene_set_id"] == "set_a" else 0.1
                rows.append(
                    _pathway_row(
                        record,
                        method_id,
                        hypothesis,
                        p_value=p_value,
                        fdr=fdr,
                        direction="Up" if record["gene_set_id"] == "set_a" else "Down",
                        correlation=correlation,
                    )
                )
    return rows


def _bundle(tmp_path: Path) -> Path:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "display").mkdir()
    sets = [
        _set_record(
            "set_a",
            "Set A",
            raw=3,
            mapped_symbols=3,
            unmapped_symbols=0,
            mapped_genes=3,
            tested=2,
            filtered=1,
        ),
        _set_record(
            "set_b",
            "Set B",
            raw=4,
            mapped_symbols=3,
            unmapped_symbols=1,
            mapped_genes=3,
            tested=3,
            filtered=0,
        ),
        _set_record(
            "set_low",
            "Set low",
            raw=2,
            mapped_symbols=1,
            unmapped_symbols=1,
            mapped_genes=1,
            tested=1,
            filtered=0,
        ),
    ]
    _write_json(run_dir / "analysis.json", _analysis(sets))
    _write_tsv(
        run_dir / "design.tsv",
        ["sample_id", "coefficient", "value"],
        [
            [sample, coefficient, value]
            for sample, values in (
                ("s1", ("1", "0")),
                ("s2", ("1", "0")),
                ("s3", ("1", "1")),
                ("s4", ("1", "1")),
            )
            for coefficient, value in zip(
                ("(Intercept)", "conditiontreated"), values, strict=True
            )
        ],
    )
    coefficient_rows = []
    for gene_index in range(1, 6):
        status = "tested" if gene_index <= 4 else "filtered"
        for coefficient, estimate in (
            ("(Intercept)", "2.3"),
            ("conditiontreated", "1.2"),
        ):
            coefficient_rows.append(
                [
                    f"gene{gene_index}",
                    status,
                    coefficient,
                    estimate if status == "tested" else "",
                    "natural_log",
                ]
            )
    _write_tsv(
        run_dir / "coefficients.tsv",
        ["gene_id", "status", "coefficient", "estimate", "scale"],
        coefficient_rows,
    )
    result_header = [
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
    ]
    result_rows = []
    for gene_index in range(1, 6):
        if gene_index <= 4:
            result_rows.append(
                [
                    f"gene{gene_index}",
                    "treated_vs_control",
                    "tested",
                    "1.2",
                    "",
                    "5.1",
                    "12.4",
                    "F",
                    "reported",
                    str(gene_index / 100),
                    "0.04",
                    "glmQLFTest",
                    "0",
                ]
            )
        else:
            result_rows.append(
                [
                    "gene5",
                    "treated_vs_control",
                    "filtered",
                    "",
                    "",
                    "",
                    "",
                    "F",
                    "not_applicable_filtered",
                    "",
                    "",
                    "glmQLFTest",
                    "0",
                ]
            )
    _write_tsv(run_dir / "results.tsv", result_header, result_rows)
    _write_tsv(run_dir / "pathway_results.tsv", PATHWAY_HEADER, _pathway_rows(sets))
    _refresh_manifest(run_dir)
    return run_dir


@pytest.fixture(autouse=True)
def _stub_display_verification(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        display_bundle,
        "verify_display_bundle",
        lambda **kwargs: {
            "summary": {"status": "verified_complete"},
            "artifacts": [],
        },
    )


def _mutate_pathway_cell(
    run_dir: Path,
    *,
    gene_set_id: str,
    method_id: str,
    hypothesis: str,
    field: str,
    value: str,
) -> None:
    path = run_dir / "pathway_results.tsv"
    rows = list(
        csv.reader(io.StringIO(path.read_text(encoding="utf-8")), delimiter="\t")
    )
    header = rows[0]
    field_index = header.index(field)
    for row in rows[1:]:
        record = dict(zip(header, row, strict=True))
        if (
            record["gene_set_id"] == gene_set_id
            and record["method_id"] == method_id
            and record["hypothesis"] == hypothesis
        ):
            row[field_index] = value
            break
    else:  # pragma: no cover - fixture programming error
        raise AssertionError("requested pathway fixture row was not found")
    _write_tsv(path, header, rows[1:])
    _refresh_manifest(run_dir)


@pytest.mark.unit
def test_pathway_summary_separates_self_contained_and_competitive_tests(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)

    summary = summarize_run(run_dir)

    assert summary["schema_version"] == "1.1"
    assert summary["pathways"]["gene_set_count"] == 3
    assert summary["pathways"]["pathway_result_row_count"] == 15
    assert summary["pathways"]["self_contained"]["limma_fry"]["directional"] == {
        "status_counts": {"tested": 2, "not_tested": 1},
        "significance_counts": {
            "fdr_le_0_05": 0,
            "fdr_gt_0_05": 2,
            "not_tested": 1,
        },
    }
    assert summary["pathways"]["competitive"]["limma_camera"]["directional"][
        "status_counts"
    ] == {"tested": 2, "not_tested": 1}
    assert "limma_camera" not in summary["pathways"]["self_contained"]
    assert {item["role"] for item in summary["artifacts"]} >= {"pathway_results"}


@pytest.mark.unit
def test_pathway_summary_rejects_fdr_family_pooling(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    _mutate_pathway_cell(
        run_dir,
        gene_set_id="set_a",
        method_id="limma_mroast",
        hypothesis="directional",
        field="fdr",
        value="0.03",
    )

    with pytest.raises(InputIntegrityError, match="Benjamini-Hochberg"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_summary_rejects_method_test_class_mislabel(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    _mutate_pathway_cell(
        run_dir,
        gene_set_id="set_a",
        method_id="limma_camera",
        hypothesis="directional",
        field="test_class",
        value="self_contained",
    )

    with pytest.raises(InputIntegrityError, match="mislabels"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_summary_rejects_mapping_arithmetic_tampering(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    _mutate_pathway_cell(
        run_dir,
        gene_set_id="set_a",
        method_id="limma_fry",
        hypothesis="directional",
        field="mapping_rate",
        value="0.5",
    )

    with pytest.raises(InputIntegrityError, match="mapping"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_summary_cross_checks_frozen_source_snapshot_identity(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["pathway_analysis"]["gene_sets"]["gmt"]["sha256"] = "9" * 64
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="captured input snapshot"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_summary_rejects_missing_grid_row(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "pathway_results.tsv"
    rows = list(
        csv.reader(io.StringIO(path.read_text(encoding="utf-8")), delimiter="\t")
    )
    _write_tsv(path, rows[0], rows[2:])
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="complete declared"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_schema_1_1_requires_display(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    (run_dir / "display").rmdir()

    with pytest.raises(InputIntegrityError, match="display directory"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_schema_identity_must_match_manifest(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["schema_version"] = "1.0"
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="analysis identity"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_summary_rejects_camera_correlation_vif_tampering(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    _mutate_pathway_cell(
        run_dir,
        gene_set_id="set_b",
        method_id="limma_camera",
        hypothesis="directional",
        field="vif_used",
        value="1.1",
    )

    with pytest.raises(InputIntegrityError, match="variance-inflation"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_summary_rejects_not_tested_inferential_outcome(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    _mutate_pathway_cell(
        run_dir,
        gene_set_id="set_low",
        method_id="limma_fry",
        hypothesis="mixed",
        field="p_value",
        value="0.5",
    )

    with pytest.raises(InputIntegrityError, match="blank inferential outcomes"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_pathway_summary_rejects_inconsistent_mroast_proportions(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    _mutate_pathway_cell(
        run_dir,
        gene_set_id="set_a",
        method_id="limma_mroast",
        hypothesis="mixed",
        field="proportion_up",
        value="0.5",
    )

    with pytest.raises(InputIntegrityError, match="rotation proportions"):
        summarize_run(run_dir)
