from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

import pytest

from rnaseq_downstream.edger_backend import run_edger_ql
from rnaseq_downstream.run_summary import summarize_run
from tests.integration.test_edger_backend import (
    _analysis_request,
    _contrasts,
    _design,
    _display_configuration,
    _validated_featurecounts_bundle,
)


_PATHWAY_HEADER = [
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


def _locked_runtime() -> tuple[str, str]:
    rscript = os.environ.get("RNASEQ_P0_RSCRIPT")
    r_library = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not rscript or not r_library:
        pytest.skip("set RNASEQ_P0_RSCRIPT and RNASEQ_P0_R_LIBRARY")
    return rscript, r_library


def _write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(root: Path, *, duplicate_gmt_member: bool = False) -> Path:
    samples = [f"s{index}" for index in range(1, 9)]
    count_lines = ["\t".join(["gene_id", *samples])]
    for gene_index in range(1, 121):
        values: list[int] = []
        for sample_index in range(8):
            value = 35 + gene_index % 13 + (gene_index * (sample_index + 3)) % 9
            if sample_index >= 4 and gene_index <= 24:
                value *= 4
            values.append(value)
        count_lines.append(
            "\t".join([f"g{gene_index}", *(str(value) for value in values)])
        )
    counts = _write(root / "counts.tsv", "\n".join(count_lines) + "\n")

    annotation_lines = ["gene_id\tsymbol"]
    annotation_lines.extend(f"g{index}\tS{index}" for index in range(1, 121))
    annotation_lines.extend(["outside_a\tAMBIG", "outside_b\tAMBIG"])
    annotation = _write(
        root / "annotation.tsv", "\n".join(annotation_lines) + "\n"
    )

    set_z = [*(f"S{index}" for index in range(1, 11)), "AMBIG", "MISSING"]
    if duplicate_gmt_member:
        set_z.append("S1")
    gmt_lines = [
        "\t".join(["SET_Z", "mapping audit set", *set_z]),
        "\t".join(
            ["SET_UNIVERSE", "all tested genes", *(f"S{i}" for i in range(1, 121))]
        ),
        "SET_SMALL\tbelow minimum\tS1",
        "\t".join(
            ["SET_A", "null-like set", *(f"S{i}" for i in range(31, 41))]
        ),
    ]
    gmt = _write(root / "sets.gmt", "\n".join(gmt_lines) + "\n")

    request = {
        "schema_version": "1.1",
        "kind": "edger_ql_backend_request",
        "execution_scope": "pathway_r_integration_test",
        "analysis_request": {
            "kind": "pathway_r_integration_fixture",
            "sha256": "a" * 64,
        },
        "input_evidence": {
            "kind": "pathway_r_integration_fixture",
            "sha256": "b" * 64,
        },
        "input": {
            "input_semantics": "benchmark_kernel_integer_counts",
            "sample_order": samples,
            "metadata": {
                "path": str(root / "metadata.tsv"),
                "sample_id_column": "sample_id",
            },
            "metadata_values": {
                "sample_id": samples,
                "condition": ["control"] * 4 + ["treated"] * 4,
            },
            "gene_id_policy": {
                "internal_key": "gene_id",
                "symbol_role": "display_only",
                "strip_version": False,
            },
            "benchmark_counts": {"matrix_path": str(counts)},
        },
        "design": {
            "intercept": True,
            "terms": ["condition"],
            "variables": {
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                }
            },
        },
        "contrasts": [
            {
                "contrast_id": "treated_vs_control",
                "weights": {"conditiontreated": 1},
                "lfc_threshold": 0.5,
            }
        ],
        "gene_sets": {
            "gmt": {
                "path": str(gmt),
                "declared_path": "sets.gmt",
                "sha256": _sha256(gmt),
                "size_bytes": gmt.stat().st_size,
                "collection": "integration_fixture",
                "version": "1",
                "identifier_type": "symbol",
            },
            "annotation": {
                "path": str(annotation),
                "declared_path": "annotation.tsv",
                "sha256": _sha256(annotation),
                "size_bytes": annotation.stat().st_size,
                "name": "integration_fixture_annotation",
                "version": "1",
                "gene_id_column": "gene_id",
                "symbol_column": "symbol",
            },
            "minimum_tested_genes": 5,
        },
    }
    request_path = root / "request.json"
    request_path.write_text(
        json.dumps(request, allow_nan=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return request_path


def _invoke(request: Path, output: Path) -> subprocess.CompletedProcess[str]:
    rscript, r_library = _locked_runtime()
    environment = os.environ.copy()
    environment["R_LIBS_USER"] = r_library
    script = (
        Path(__file__).resolve().parents[2]
        / "rnaseq_downstream"
        / "r_scripts"
        / "edger_ql.R"
    )
    return subprocess.run(
        [rscript, "--vanilla", str(script), str(request), str(output)],
        check=False,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="strict",
        env=environment,
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == _PATHWAY_HEADER
        return list(reader)


def _replace_source(request_path: Path, role: str, content: str) -> None:
    document = json.loads(request_path.read_text(encoding="utf-8"))
    source = Path(document["gene_sets"][role]["path"])
    source.write_text(content, encoding="utf-8")
    document["gene_sets"][role]["sha256"] = _sha256(source)
    document["gene_sets"][role]["size_bytes"] = source.stat().st_size
    request_path.write_text(
        json.dumps(document, allow_nan=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_pathway_backend_uses_one_fit_and_emits_audited_long_results(
    tmp_path: Path,
) -> None:
    request = _fixture(tmp_path / "valid")
    output = tmp_path / "valid" / "output"

    completed = _invoke(request, output)

    assert completed.returncode == 0, completed.stderr
    response = json.loads(completed.stdout)
    assert response["schema_version"] == "1.1"
    assert response["status"] == "success"
    assert {path.name for path in output.iterdir()} == {
        "analysis.json",
        "backend_manifest.json",
        "coefficients.tsv",
        "design.tsv",
        "pathway_results.tsv",
        "results.tsv",
    }
    manifest = json.loads((output / "backend_manifest.json").read_text())
    assert manifest["schema_version"] == "1.1"
    assert [member["path"] for member in manifest["members"]] == [
        "analysis.json",
        "coefficients.tsv",
        "design.tsv",
        "results.tsv",
        "pathway_results.tsv",
    ]

    rows = _read_tsv(output / "pathway_results.tsv")
    assert len(rows) == 4 * 5
    assert [
        (row["contrast_id"], row["gene_set_id"], row["method_id"], row["hypothesis"])
        for row in rows
    ] == sorted(
        (
            row["contrast_id"],
            row["gene_set_id"],
            row["method_id"],
            row["hypothesis"],
        )
        for row in rows
    )

    mapping_rows = [row for row in rows if row["gene_set_id"] == "SET_Z"]
    assert len(mapping_rows) == 5
    for row in mapping_rows:
        assert row["gmt_member_count_raw"] == "12"
        assert row["gmt_symbol_count_unique"] == "12"
        assert row["mapped_symbol_count_unique"] == "10"
        assert row["ambiguous_symbol_count_unique"] == "1"
        assert row["unmapped_symbol_count_unique"] == "1"
        assert float(row["mapping_rate"]) == pytest.approx(10 / 12)
        assert row["mapped_gene_id_count_unique"] == "10"
        assert row["tested_gene_count"] == "10"
        assert row["filtered_gene_count"] == "0"

    small_rows = [row for row in rows if row["gene_set_id"] == "SET_SMALL"]
    assert all(row["status"] == "not_tested" for row in small_rows)
    assert {
        row["status_reason"] for row in small_rows
    } == {"tested_gene_count_below_minimum"}
    assert all(row["p_value"] == row["fdr"] == "" for row in small_rows)

    universe_camera = next(
        row
        for row in rows
        if row["gene_set_id"] == "SET_UNIVERSE"
        and row["method_id"] == "limma_camera"
    )
    assert universe_camera["status"] == "not_tested"
    assert (
        universe_camera["status_reason"]
        == "competitive_test_requires_background_genes"
    )
    tested = [row for row in rows if row["status"] == "tested"]
    assert tested
    assert all(
        0 <= float(row["p_value"]) <= 1 and 0 <= float(row["fdr"]) <= 1
        for row in tested
    )
    assert all(
        row["direction"] == "Mixed"
        for row in tested
        if row["hypothesis"] == "mixed"
    )
    mroast = [row for row in tested if row["method_id"] == "limma_mroast"]
    assert all(
        0 <= float(row["proportion_down"]) <= 1
        and 0 <= float(row["proportion_up"]) <= 1
        and row["rotation_status"] == "performed_fixed_seed_9999_rotations"
        for row in mroast
    )
    for gene_set_id in {row["gene_set_id"] for row in mroast}:
        pair = [row for row in mroast if row["gene_set_id"] == gene_set_id]
        assert len(pair) == 2
        assert {row["proportion_down"] for row in pair} == {
            pair[0]["proportion_down"]
        }
        assert {row["proportion_up"] for row in pair} == {pair[0]["proportion_up"]}
    assert all(
        row["proportion_down"] == row["proportion_up"] == ""
        for row in tested
        if row["method_id"] != "limma_mroast"
    )
    camera = [row for row in tested if row["method_id"] == "limma_camera"]
    assert all(row["correlation_status"] == "estimated_set_specific" for row in camera)
    for row in camera:
        raw = float(row["correlation_estimate_raw"])
        effective = float(row["correlation_effective"])
        vif = float(row["vif_used"])
        size = int(row["method_ngenes"])
        assert effective == pytest.approx(max(0, raw))
        assert vif == pytest.approx(max(1, 1 + (size - 1) * raw))
        assert all(math.isfinite(value) for value in (raw, effective, vif))

    analysis = json.loads((output / "analysis.json").read_text())
    assert analysis["schema_version"] == "1.1"
    pathways = analysis["pathway_analysis"]
    assert pathways["tested_universe_gene_count"] == 120
    assert [
        item["gene_set_id"] for item in pathways["gene_sets"]["sets"]
    ] == ["SET_A", "SET_SMALL", "SET_UNIVERSE", "SET_Z"]
    assert pathways["methods"]["limma_mroast"]["generic"] == "limma::mroast"
    assert pathways["methods"]["limma_mroast"]["dispatch"] == (
        "edgeR::mroast.DGEGLM"
    )
    assert pathways["methods"]["limma_fry"]["dispatch"] == "edgeR::fry.DGEGLM"
    assert pathways["methods"]["limma_camera"]["dispatch"] == (
        "edgeR::camera.DGEGLM"
    )
    contrast = pathways["contrasts"][0]
    assert contrast["gene_level_lfc_threshold"] == 0.5
    assert contrast["pathway_statistical_null"] == "zero_effect"
    assert contrast["gene_level_lfc_threshold_applied_to_pathways"] is False
    assert contrast["rotation"] == {
        "method_id": "limma_mroast",
        "seed": 1729,
        "seed_policy": "base_1729_plus_zero_based_declared_contrast_index",
        "rng": {
            "kind": "Mersenne-Twister",
            "normal.kind": "Inversion",
            "sample.kind": "Rejection",
        },
        "nrot": 9999,
        "reset_before_each_contrast": True,
    }

    repeated_output = tmp_path / "valid" / "repeated-output"
    repeated = _invoke(request, repeated_output)
    assert repeated.returncode == 0, repeated.stderr
    assert (repeated_output / "pathway_results.tsv").read_bytes() == (
        output / "pathway_results.tsv"
    ).read_bytes()


def test_pathway_backend_hard_fails_duplicate_gmt_member_without_outputs(
    tmp_path: Path,
) -> None:
    request = _fixture(tmp_path / "duplicate", duplicate_gmt_member=True)
    output = tmp_path / "duplicate" / "must-not-exist"

    completed = _invoke(request, output)

    assert completed.returncode == 4
    response = json.loads(completed.stdout)
    assert response["schema_version"] == "1.1"
    assert response["status"] == "error"
    assert response["errors"][0]["details"]["reason"] == "gmt_duplicate_symbol"
    assert not output.exists()


def test_pathway_backend_rejects_annotation_ids_duplicated_after_version_policy(
    tmp_path: Path,
) -> None:
    request = _fixture(tmp_path / "duplicate-annotation-id")
    _replace_source(
        request,
        "annotation",
        "gene_id\tsymbol\ng1.1\tA\ng1.2\tB\n",
    )
    document = json.loads(request.read_text(encoding="utf-8"))
    document["input"]["gene_id_policy"]["strip_version"] = True
    request.write_text(
        json.dumps(document, allow_nan=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "duplicate-annotation-id" / "must-not-exist"

    completed = _invoke(request, output)

    assert completed.returncode == 4
    response = json.loads(completed.stdout)
    assert response["schema_version"] == "1.1"
    assert response["errors"][0]["details"]["reason"] == (
        "annotation_gene_id_duplicated_after_version_policy"
    )
    assert not output.exists()


def test_pathway_backend_requires_schema_11_exactly_when_gene_sets_are_present(
    tmp_path: Path,
) -> None:
    request = _fixture(tmp_path / "schema-mismatch")
    document = json.loads(request.read_text(encoding="utf-8"))
    document["schema_version"] = "1.0"
    request.write_text(
        json.dumps(document, allow_nan=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "schema-mismatch" / "must-not-exist"

    completed = _invoke(request, output)

    assert completed.returncode == 4
    response = json.loads(completed.stdout)
    assert response["schema_version"] == "1.0"
    assert response["errors"][0]["details"]["reason"] == (
        "backend_request_identity_invalid"
    )
    assert not output.exists()


def test_public_pathway_run_is_recaptured_verified_and_published_atomically(
    tmp_path: Path,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "public"
    bundle = _validated_featurecounts_bundle(root, layout="combined_matrix")
    core_output = root / "core-only"
    core_result = run_edger_ql(
        _analysis_request(root, bundle),
        core_output,
        rscript=rscript,
        r_library=r_library,
    )
    assert core_result["schema_version"] == "1.0"
    assert "pathway_analysis" not in core_result["data"]
    assert "pathway_analysis" not in core_result["analysis"]
    annotation = _write(
        root / "annotation.tsv",
        "gene_id\tsymbol\n"
        + "".join(f"gene_{index}\tS{index}\n" for index in range(1, 41)),
    )
    gmt = _write(
        root / "sets.gmt",
        "SET_DE\tfixture DE genes\t"
        + "\t".join(f"S{index}" for index in range(1, 9))
        + "\nSET_NULL\tfixture null genes\t"
        + "\t".join(f"S{index}" for index in range(21, 31))
        + "\n",
    )
    request = {
        "schema_version": "1.1",
        "validated_input_bundle": str(bundle),
        "design": _design(),
        "contrasts": _contrasts(),
        "display": _display_configuration(),
        "gene_sets": {
            "gmt": {
                "path": str(gmt),
                "sha256": _sha256(gmt),
                "collection": "integration_fixture",
                "version": "1",
                "identifier_type": "symbol",
            },
            "annotation": {
                "path": str(annotation),
                "sha256": _sha256(annotation),
                "name": "integration_fixture_annotation",
                "version": "1",
                "gene_id_column": "gene_id",
                "symbol_column": "symbol",
            },
            "minimum_tested_genes": 5,
        },
    }
    request_path = root / "analysis-request-pathways.json"
    request_path.write_text(
        json.dumps(request, allow_nan=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    output = root / "published"

    result = run_edger_ql(
        request_path,
        output,
        rscript=rscript,
        r_library=r_library,
    )

    assert result["schema_version"] == "1.1"
    assert result["status"] == "success"
    assert (output / "pathway_results.tsv").is_file()
    assert (output / "display" / "display_manifest.json").is_file()
    for name in ("results.tsv", "coefficients.tsv", "design.tsv"):
        assert (output / name).read_bytes() == (core_output / name).read_bytes()
    summary = summarize_run(output)
    assert summary["status"] == "verified_complete"
    assert summary["pathways"]["gene_set_count"] == 2
    assert summary["pathways"]["pathway_result_row_count"] == 10
