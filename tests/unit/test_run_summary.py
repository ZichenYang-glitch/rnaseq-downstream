"""Unit tests for fail-closed result-bundle summarization."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import pytest

from rnaseq_downstream.errors import InputIntegrityError, InvalidRequestError
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


def _input_evidence() -> dict[str, object]:
    return {
        "kind": "validated_input_bundle",
        "plan_id": PLAN_ID,
        "bundle_path": "/captured/input-bundle",
        "validation_bundle_artifacts": [
            {
                "kind": "consumed_validation_evidence",
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
                "private_relative_path": "metadata.tsv",
                "coupling": "copied_and_hashed_from_one_source_stream",
            }
        ],
        "digest_coupling": "private_copy_and_hash_same_source_stream",
    }


def _tsv(header: list[str], rows: list[list[str]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _write_json(path: Path, document: object) -> None:
    path.write_text(
        json.dumps(document, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _member(path: Path) -> dict[str, object]:
    content = path.read_bytes()
    return {
        "path": path.name,
        "sha256": hashlib.sha256(content).hexdigest(),
        "size_bytes": len(content),
    }


def _refresh_manifest(run_dir: Path, *, scope: str = "validated_p0_input") -> None:
    evidence = _input_evidence()
    manifest = {
        "schema_version": "1.0",
        "kind": "edger_ql_backend_manifest",
        "backend": "edgeR_QL",
        "runtime_identity": RUNTIME,
        "execution_scope": scope,
        "input_evidence": evidence,
        "members": [
            _member(run_dir / name)
            for name in (
                "analysis.json",
                "coefficients.tsv",
                "design.tsv",
                "results.tsv",
            )
        ],
    }
    _write_json(run_dir / "backend_manifest.json", manifest)


def _bundle(tmp_path: Path) -> Path:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    evidence = _input_evidence()
    contrasts = [
        {
            "contrast_id": "treated_vs_control",
            "weights": {"conditiontreated": 1},
            "lfc_threshold": 0,
            "test_method": "glmQLFTest",
            "treat_null": None,
            "estimability_residual": 0,
            "estimability_tolerance": 1e-7,
        },
        {
            "contrast_id": "treated_vs_control_treat",
            "weights": {"conditiontreated": 1},
            "lfc_threshold": 0.5,
            "test_method": "glmTreat",
            "treat_null": "interval",
            "estimability_residual": 0,
            "estimability_tolerance": 1e-7,
        },
    ]
    analysis = {
        "schema_version": "1.0",
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
            {
                "step": "glmQLFit",
                "arguments": {"legacy": False, "robust": True},
            },
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
        "contrasts": contrasts,
        "genes": {"total": 2, "tested": 1, "filtered": 1},
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
    }
    _write_json(run_dir / "analysis.json", analysis)
    (run_dir / "design.tsv").write_bytes(
        _tsv(
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
    )
    (run_dir / "coefficients.tsv").write_bytes(
        _tsv(
            ["gene_id", "status", "coefficient", "estimate", "scale"],
            [
                ["gene1", "tested", "(Intercept)", "2.3", "natural_log"],
                ["gene1", "tested", "conditiontreated", "1.2", "natural_log"],
                ["gene2", "filtered", "(Intercept)", "", "natural_log"],
                ["gene2", "filtered", "conditiontreated", "", "natural_log"],
            ],
        )
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
    (run_dir / "results.tsv").write_bytes(
        _tsv(
            result_header,
            [
                [
                    "gene1",
                    "treated_vs_control",
                    "tested",
                    "1.2",
                    "",
                    "5.1",
                    "12.4",
                    "F",
                    "reported",
                    "0.001",
                    "0.001",
                    "glmQLFTest",
                    "0",
                ],
                [
                    "gene2",
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
                ],
                [
                    "gene1",
                    "treated_vs_control_treat",
                    "tested",
                    "1.1",
                    "1.2",
                    "5.1",
                    "",
                    "not_reported_by_glmTreat",
                    "not_reported",
                    "0.02",
                    "0.02",
                    "glmTreat",
                    "0.5",
                ],
                [
                    "gene2",
                    "treated_vs_control_treat",
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
                    "0.5",
                ],
            ],
        )
    )
    _refresh_manifest(run_dir)
    return run_dir


@pytest.mark.unit
def test_summarize_verifies_complete_bundle_and_counts(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)

    summary = summarize_run(run_dir)

    assert summary["status"] == "verified_complete"
    assert summary["execution_scope"] == "validated_p0_input"
    assert summary["plan_id"] == PLAN_ID
    assert summary["gene_count"] == 2
    assert summary["result_row_count"] == 4
    assert summary["contrasts"][0]["status_counts"]["tested"] == 1
    assert summary["contrasts"][0]["status_counts"]["filtered"] == 1
    assert summary["contrasts"][0]["significance_counts"] == {
        "fdr_le_0_05": 1,
        "fdr_gt_0_05": 0,
        "not_tested": 1,
    }
    assert len(summary["artifacts"]) == 5


@pytest.mark.unit
def test_summarize_detects_member_tampering(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    with (run_dir / "results.tsv").open("ab") as handle:
        handle.write(b"\n")

    with pytest.raises(InputIntegrityError, match="does not match"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_rejects_extra_bundle_entry(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    (run_dir / "unmanifested.txt").write_text("unsafe\n", encoding="utf-8")

    with pytest.raises(InputIntegrityError) as caught:
        summarize_run(run_dir)

    assert caught.value.details["unexpected_files"] == ["unmanifested.txt"]


@pytest.mark.unit
def test_summarize_rejects_glmtreat_fabricated_statistic(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "results.tsv"
    content = path.read_text(encoding="utf-8")
    content = content.replace(
        "1.1\t1.2\t5.1\t\tnot_reported_by_glmTreat",
        "1.1\t1.2\t5.1\t9.9\tnot_reported_by_glmTreat",
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="unreported statistic"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_rejects_duplicate_gene_contrast_row(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "results.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join([*lines, lines[1]]) + "\n", encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="duplicate gene/contrast"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_recomputes_within_contrast_bh_fdr(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "results.tsv"
    content = path.read_text(encoding="utf-8")
    content = content.replace(
        "0.001\t0.001\tglmQLFTest",
        "0.001\t0.049\tglmQLFTest",
        1,
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="Benjamini-Hochberg"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_public_summarize_rejects_benchmark_kernel_output(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    analysis_path = run_dir / "analysis.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    analysis["execution_scope"] = "backend_kernel_only"
    analysis["input_evidence"] = {
        "kind": "benchmark_kernel_inputs",
        "benchmark_id": "b" * 64,
    }
    _write_json(analysis_path, analysis)
    _refresh_manifest(run_dir, scope="backend_kernel_only")

    with pytest.raises(InvalidRequestError, match="benchmark-kernel"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_rejects_duplicate_json_keys_with_valid_manifest_hash(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    content = path.read_text(encoding="utf-8")
    content = content.replace(
        '"kind": "edger_ql_analysis",',
        '"kind": "edger_ql_analysis",\n  "kind": "edger_ql_analysis",',
        1,
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="strict UTF-8 JSON"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_rejects_failed_row_in_verified_complete_bundle(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "results.tsv"
    content = path.read_text(encoding="utf-8")
    content = content.replace(
        "gene2\ttreated_vs_control\tfiltered\t\t\t\t\tF\tnot_applicable_filtered",
        "gene2\ttreated_vs_control\tfailed\t\t\t\t\tF\tnot_applicable_failed",
        1,
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="only filtered or tested"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_rejects_duplicate_design_columns(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["design"]["columns"] = ["(Intercept)", "(Intercept)"]
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="design columns are duplicated"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_rejects_boolean_gene_counts(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["genes"] = {"total": 2, "tested": True, "filtered": True}
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="gene-count fields are invalid"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_rejects_contrast_outside_recorded_estimability_tolerance(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["contrasts"][0]["estimability_residual"] = 1e-3
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="estimability tolerance"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_summarize_rejects_full_length_divide_without_imported_replicates(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["input_semantics"] = "salmon_quant_dirs_full_length"
    analysis["route_observed"] = {
        "constructor": "edgeR::DGEListFromTximport",
        "transcript_length_offset": True,
        "countsFromAbundance": "no",
        "dropInfReps": False,
        "divide": True,
        "inferential_replicates_imported": False,
    }
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="divide policy"):
        summarize_run(run_dir)
