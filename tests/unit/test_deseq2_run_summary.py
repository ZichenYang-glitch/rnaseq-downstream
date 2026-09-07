"""Unit tests for independent DESeq2 bundle verification."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import pytest

from rnaseq_downstream.errors import InputIntegrityError
from rnaseq_downstream.run_summary import summarize_run


RUNTIME = {
    "R": "4.6.1",
    "Bioconductor": "3.23",
    "BiocVersion_package": "3.23.1",
    "DESeq2": "1.52.0",
    "apeglm": "1.34.0",
    "tximport": "1.40.0",
}
PLAN_ID = "a" * 64
RESULT_HEADER = [
    "gene_id",
    "contrast_id",
    "status",
    "status_reason",
    "baseMean",
    "logFC",
    "unshrunk_logFC",
    "lfcSE",
    "statistic",
    "statistic_type",
    "statistic_hypothesis",
    "PValue",
    "FDR",
    "fdr_basis",
    "test_method",
    "lfc_threshold",
    "shrinkage_method",
]


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


def _rounding_audit(*, mode: str = "none_integer_input") -> dict[str, object]:
    rule = (
        "not_applicable_exact_integer_input"
        if mode == "none_integer_input"
        else "R_base_round_IEC_60559_ties_to_even"
    )
    return {
        "mode": mode,
        "round_function": "none" if mode == "none_integer_input" else "base::round",
        "rule": rule,
        "hash_serialization": "R_serialize_version_3_locked_runtime",
        "source_matrix_sha256": "e" * 64,
        "rounded_matrix_sha256": "e" * 64,
        "gene_count": 3,
        "sample_count": 4,
        "cell_count": 12,
        "changed_cell_count": 0,
        "max_absolute_delta": 0.0,
        "absolute_delta_sum": 0.0,
        "per_sample": [
            {
                "sample_id": sample,
                "cell_count": 3,
                "changed_cell_count": 0,
                "max_absolute_delta": 0.0,
                "absolute_delta_sum": 0.0,
                "total_count_before": 100.0,
                "total_count_after": 100.0,
                "total_count_delta": 0.0,
            }
            for sample in ("s1", "s2", "s3", "s4")
        ],
    }


def _three_prime_route(
    *,
    replicate_count: int = 0,
    rounding_audit: dict[str, object] | None = None,
) -> dict[str, object]:
    return {
        "constructor": "DESeq2::DESeqDataSetFromMatrix",
        "count_source": "txi$counts",
        "count_semantics": "salmon_estimated_counts_rounded_for_DESeq2",
        "transcript_length_offset": False,
        "gene_length_correction": False,
        "countsFromAbundance": "no",
        "dropInfReps": False,
        "inferential_replicates_imported": replicate_count > 0,
        "inferential_replicate_state": "all" if replicate_count > 0 else "none",
        "inferential_replicates_per_sample": replicate_count,
        "inferential_replicate_method": "Gibbs" if replicate_count > 0 else None,
        "inferential_replicates_used_for_inference": False,
        "inferential_replicates_unused_reason": (
            "DESeq2 1.52.0 does not consume tximport infReps in this backend; "
            "they are imported only to verify declared evidence."
        ),
        "rounding_disclosure": (
            "The toolkit explicitly calls base::round before "
            "DESeqDataSetFromMatrix; no length offset is attached."
        ),
        "rounding_audit": rounding_audit
        or _rounding_audit(mode="explicit_before_matrix_constructor"),
    }


def _refresh_manifest(run_dir: Path) -> None:
    manifest = {
        "schema_version": "1.0",
        "kind": "deseq2_backend_manifest",
        "backend": "DESeq2",
        "runtime_identity": RUNTIME,
        "execution_scope": "validated_p1_deseq2_input",
        "input_evidence": _input_evidence(),
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


def _bundle(tmp_path: Path, *, mode: str = "wald") -> Path:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    contrast = {
        "contrast_id": "treated_vs_control",
        "weights": {"conditiontreated": 1.0},
        "lfc_threshold": 0.0,
        "test_method": (
            "DESeq2::results_LRT" if mode == "lrt" else "DESeq2::results_Wald"
        ),
        "alternative_hypothesis": (
            "full_vs_reduced_omnibus"
            if mode == "lrt"
            else "greaterAbs_at_zero_equivalent_two_sided"
        ),
        "independent_filter_threshold": 0.5,
        "independent_filter_theta": 0.2,
        "independent_filter_alpha": 0.1,
        "cooks_filter_applied": True,
        "cooks_cutoff": 99.0,
        "coefficient_name": None,
        "shrinkage_method": "none",
        "shrinkage_nonconvergence_count": 0,
        "estimability_residual": 0.0,
        "estimability_tolerance": 1e-7,
    }
    test: dict[str, object] = {
        "mode": mode,
        "shrinkage": "none",
        "reduced": None,
    }
    reporting_role = "tested_contrast"
    if mode == "lrt":
        test = {
            "mode": "lrt",
            "shrinkage": "none",
            "reduced": {
                "intercept": True,
                "terms": [],
                "columns": ["(Intercept)"],
                "rank": 1,
                "residual_df": 3,
                "nesting_residual": 0.0,
                "nesting_tolerance": 1e-7,
            },
        }
        reporting_role = "reported_effect_not_lrt_hypothesis"
    analysis = {
        "schema_version": "1.0",
        "kind": "deseq2_analysis",
        "backend": "DESeq2",
        "execution_scope": "validated_p1_deseq2_input",
        "analysis_request": {
            "path": "/captured/request.json",
            "sha256": "b" * 64,
            "schema_version": "1.2",
            "backend": "deseq2",
        },
        "input_evidence": _input_evidence(),
        "runtime_identity": RUNTIME,
        "input_semantics": "featurecounts_integer",
        "route_observed": {
            "constructor": "DESeq2::DESeqDataSetFromMatrix",
            "count_source": "validated_featureCounts_integer_counts",
            "count_semantics": "integer",
            "transcript_length_offset": False,
            "gene_length_correction": False,
            "rounding_audit": _rounding_audit(),
            "inferential_replicates_imported": False,
            "inferential_replicates_used_for_inference": False,
            "inferential_replicates_unused_reason": "not_applicable",
        },
        "pipeline": [
            {
                "step": "construct_DESeqDataSet",
                "constructor": "DESeq2::DESeqDataSetFromMatrix",
            },
            {
                "step": "DESeq",
                "arguments": {
                    "test": "LRT" if mode == "lrt" else "Wald",
                    "fitType": "parametric",
                    "sfType": "ratio",
                    "betaPrior": False,
                    "minReplicatesForReplace": 7,
                    "useT": False,
                    "minmu": 0.5,
                    "parallel": False,
                },
            },
            {
                "step": "results",
                "arguments": {
                    "independentFiltering": True,
                    "cooksCutoff": "automatic",
                    "alpha": 0.1,
                    "pAdjustMethod": "BH",
                    "altHypothesis": "greaterAbs",
                    "parallel": False,
                },
            },
            {
                "step": "lfcShrink",
                "method": "none",
                "arguments": None,
                "inferential_columns_changed": False,
            },
        ],
        "defaults": {
            "fit_type_requested": "parametric",
            "fit_type_resolved": "parametric",
            "size_factor_type": "ratio",
            "beta_prior": False,
            "min_replicates_for_replace": 7,
            "independent_filtering": True,
            "cooks_cutoff": {
                "requested": "automatic",
                "resolved_f_quantile": 0.99,
                "resolved_value": 99.0,
                "numerator_df": 2,
                "denominator_df": 2,
            },
            "alpha": 0.1,
            "p_adjust_method": "BH",
            "results_alt_hypothesis": "greaterAbs",
            "use_t": False,
            "minmu": 0.5,
            "parallel": False,
            "outlier_replacement_count": 0,
        },
        "design": {
            "intercept": True,
            "terms": ["condition"],
            "variable_types": {"condition": "factor"},
            "factor_levels": {"condition": ["control", "treated"]},
            "columns": ["(Intercept)", "conditiontreated"],
            "sample_count": 4,
            "rank": 2,
            "coefficient_mapping": [
                {
                    "design_coefficient": "(Intercept)",
                    "deseq2_result_name": "Intercept",
                    "position": 1,
                },
                {
                    "design_coefficient": "conditiontreated",
                    "deseq2_result_name": "condition_treated_vs_control",
                    "position": 2,
                },
            ],
            "residual_df": 2,
            "qr_tolerance": 1e-7,
        },
        "test": test,
        "contrasts": [contrast],
        "reporting_effect": [
            {
                "contrast_id": "treated_vs_control",
                "role": reporting_role,
                "weights": {"conditiontreated": 1.0},
                "coefficient_name": None,
            }
        ],
        "genes": {
            "total": 3,
            "result_rows": 3,
            "status_counts": {
                "failed": 1,
                "filtered": 1,
                "not_estimable": 0,
                "not_tested": 0,
                "tested": 1,
            },
        },
        "status_vocabulary": [
            "filtered",
            "not_tested",
            "not_estimable",
            "failed",
            "tested",
        ],
        "result_logFC_scale": "log2",
        "coefficient_scale": "log2",
        "multiple_testing": (
            "Benjamini-Hochberg within each reporting contrast after "
            "DESeq2 independent filtering"
        ),
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
    failed_reason = (
        "lrt_full_model_nonconvergence"
        if mode == "lrt"
        else "wald_model_nonconvergence"
    )
    (run_dir / "coefficients.tsv").write_bytes(
        _tsv(
            [
                "gene_id",
                "status",
                "status_reason",
                "coefficient",
                "estimate",
                "scale",
            ],
            [
                ["gene1", "tested", "fitted", "(Intercept)", "2", "log2"],
                ["gene1", "tested", "fitted", "conditiontreated", "1", "log2"],
                ["gene2", "failed", failed_reason, "(Intercept)", "1", "log2"],
                ["gene2", "failed", failed_reason, "conditiontreated", "0.5", "log2"],
                ["gene3", "tested", "fitted", "(Intercept)", "0.5", "log2"],
                ["gene3", "tested", "fitted", "conditiontreated", "0.25", "log2"],
            ],
        )
    )
    statistic_type = "LRT" if mode == "lrt" else "Wald"
    statistic_hypothesis = (
        "full_vs_reduced_omnibus" if mode == "lrt" else "contrast_equals_zero"
    )
    fdr_basis = (
        "omnibus_pvalue_BH"
        if mode == "lrt"
        else "contrast_pvalue_BH_after_independent_filtering"
    )
    method = "DESeq2::results_LRT" if mode == "lrt" else "DESeq2::results_Wald"
    (run_dir / "results.tsv").write_bytes(
        _tsv(
            RESULT_HEADER,
            [
                [
                    "gene1",
                    "treated_vs_control",
                    "tested",
                    "tested",
                    "10",
                    "1",
                    "1",
                    "0.2",
                    "5",
                    statistic_type,
                    statistic_hypothesis,
                    "0.01",
                    "0.02",
                    fdr_basis,
                    method,
                    "0",
                    "none",
                ],
                [
                    "gene2",
                    "treated_vs_control",
                    "failed",
                    failed_reason,
                    "8",
                    "0.5",
                    "0.5",
                    "0.2",
                    "2.5",
                    statistic_type,
                    statistic_hypothesis,
                    "0.02",
                    "0.02",
                    fdr_basis,
                    method,
                    "0",
                    "none",
                ],
                [
                    "gene3",
                    "treated_vs_control",
                    "filtered",
                    "independent_filtering",
                    "0.1",
                    "0.25",
                    "0.25",
                    "0.3",
                    "0.83",
                    statistic_type,
                    statistic_hypothesis,
                    "0.001",
                    "",
                    fdr_basis,
                    method,
                    "0",
                    "none",
                ],
            ],
        )
    )
    _refresh_manifest(run_dir)
    return run_dir


def _declare_apeglm_nonconvergence(run_dir: Path) -> None:
    analysis_path = run_dir / "analysis.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    analysis["test"]["shrinkage"] = "apeglm"
    analysis["contrasts"][0].update(
        {
            "coefficient_name": "condition_treated_vs_control",
            "shrinkage_method": "apeglm",
            "shrinkage_nonconvergence_count": 1,
        }
    )
    analysis["reporting_effect"][0]["coefficient_name"] = "condition_treated_vs_control"
    analysis["pipeline"][-1] = {
        "step": "lfcShrink",
        "method": "apeglm",
        "arguments": {"lfcThreshold": 0, "svalue": False},
        "inferential_columns_changed": False,
    }
    analysis["genes"]["status_counts"].update({"failed": 2, "tested": 0})
    _write_json(analysis_path, analysis)

    results_path = run_dir / "results.tsv"
    content = results_path.read_text(encoding="utf-8")
    content = content.replace(
        "gene1\ttreated_vs_control\ttested\ttested\t10\t1\t1\t",
        "gene1\ttreated_vs_control\tfailed\tapeglm_nonconvergence\t10\t\t1\t",
    ).replace("\tnone\n", "\tapeglm\n")
    results_path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)


@pytest.mark.unit
def test_summarize_verifies_deseq2_wald_bundle(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)

    summary = summarize_run(run_dir)

    assert summary["status"] == "verified_complete"
    assert summary["backend"] == "DESeq2"
    assert summary["runtime_identity"] == RUNTIME
    assert summary["test"]["mode"] == "wald"
    assert summary["gene_count"] == 3
    assert summary["contrasts"][0]["status_counts"]["failed"] == 1
    assert len(summary["artifacts"]) == 5


@pytest.mark.unit
def test_summarize_verifies_explicit_apeglm_nonconvergence(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    _declare_apeglm_nonconvergence(run_dir)

    summary = summarize_run(run_dir)

    assert summary["contrasts"][0]["status_counts"]["failed"] == 2
    assert summary["test"]["shrinkage"] == "apeglm"


@pytest.mark.unit
def test_deseq2_rejects_apeglm_nonconvergence_count_mismatch(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    _declare_apeglm_nonconvergence(run_dir)
    analysis_path = run_dir / "analysis.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    analysis["contrasts"][0]["shrinkage_nonconvergence_count"] = 0
    _write_json(analysis_path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="convergence count"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_bh_includes_failed_rows_with_finite_pvalues(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "results.tsv"
    content = path.read_text(encoding="utf-8")
    content = content.replace(
        "0.01\t0.02\tcontrast_pvalue_BH_after_independent_filtering",
        "0.01\t0.01\tcontrast_pvalue_BH_after_independent_filtering",
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="Benjamini-Hochberg"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_accepts_cooks_outlier_without_inferential_pvalues(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    coefficients_path = run_dir / "coefficients.tsv"
    coefficient_content = coefficients_path.read_text(encoding="utf-8").replace(
        "gene2\tfailed\twald_model_nonconvergence",
        "gene2\ttested\tfitted",
    )
    coefficients_path.write_text(coefficient_content, encoding="utf-8")

    results_path = run_dir / "results.tsv"
    result_content = results_path.read_text(encoding="utf-8")
    result_content = result_content.replace(
        "\t0.01\t0.02\tcontrast_pvalue_BH_after_independent_filtering",
        "\t0.01\t0.01\tcontrast_pvalue_BH_after_independent_filtering",
        1,
    ).replace(
        "gene2\ttreated_vs_control\tfailed\twald_model_nonconvergence\t"
        "8\t0.5\t0.5\t0.2\t2.5\tWald\tcontrast_equals_zero\t0.02\t0.02\t",
        "gene2\ttreated_vs_control\tfailed\tcooks_outlier\t"
        "8\t0.5\t0.5\t0.2\t2.5\tWald\tcontrast_equals_zero\t\t\t",
    )
    results_path.write_text(result_content, encoding="utf-8")
    _refresh_manifest(run_dir)

    summary = summarize_run(run_dir)

    assert summary["contrasts"][0]["status_counts"]["failed"] == 1


@pytest.mark.unit
def test_deseq2_malformed_shrinkage_is_structured_integrity_error(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["test"]["shrinkage"] = {"malicious": True}
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="shrinkage"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_filtered_raw_pvalue_is_excluded_from_bh(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "results.tsv"
    content = path.read_text(encoding="utf-8")
    content = content.replace(
        "0.001\t\tcontrast_pvalue_BH_after_independent_filtering",
        "0.001\t0.003\tcontrast_pvalue_BH_after_independent_filtering",
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="filtered|FDR"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_filter_threshold_equality_is_included(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "results.tsv"
    content = path.read_text(encoding="utf-8").replace(
        "gene3\ttreated_vs_control\tfiltered\tindependent_filtering\t0.1",
        "gene3\ttreated_vs_control\tfiltered\tindependent_filtering\t0.5",
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="filtered|FDR"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_accepts_clean_contrast_all_zero_exception(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    analysis_path = run_dir / "analysis.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    weights = {"(Intercept)": -1.0, "conditiontreated": 1.0}
    analysis["contrasts"][0]["weights"] = weights
    analysis["reporting_effect"][0]["weights"] = weights
    _write_json(analysis_path, analysis)

    coefficients_path = run_dir / "coefficients.tsv"
    coefficient_content = coefficients_path.read_text(encoding="utf-8")
    coefficient_content = coefficient_content.replace(
        "gene2\tfailed\twald_model_nonconvergence\t(Intercept)\t1\tlog2",
        "gene2\tfailed\twald_model_nonconvergence\t(Intercept)\t0\tlog2",
    ).replace(
        "gene3\ttested\tfitted\t(Intercept)\t0.5\tlog2",
        "gene3\ttested\tfitted\t(Intercept)\t0\tlog2",
    )
    coefficients_path.write_text(coefficient_content, encoding="utf-8")

    results_path = run_dir / "results.tsv"
    content = results_path.read_text(encoding="utf-8")
    content = content.replace(
        "gene1\ttreated_vs_control\ttested\ttested\t10\t1\t1\t0.2\t5\t",
        "gene1\ttreated_vs_control\ttested\ttested\t10\t0\t0\t0.2\t0\t",
    )
    content = content.replace("\t0.01\t0.02\t", "\t1\t1\t", 1)
    content = content.replace("\t0.02\t0.02\t", "\t0.02\t0.04\t", 1)
    results_path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    assert summarize_run(run_dir)["status"] == "verified_complete"


@pytest.mark.unit
def test_summarize_verifies_deseq2_lrt_omnibus_semantics(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path, mode="lrt")

    summary = summarize_run(run_dir)

    assert summary["test"]["mode"] == "lrt"
    assert summary["contrasts"][0]["test_method"] == "DESeq2::results_LRT"


@pytest.mark.unit
def test_deseq2_lrt_rejects_wald_statistic_label(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path, mode="lrt")
    path = run_dir / "results.tsv"
    content = path.read_text(encoding="utf-8").replace(
        "\tLRT\tfull_vs_reduced_omnibus\t",
        "\tWald\tcontrast_equals_zero\t",
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="declared test"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_rejects_runtime_identity_drift(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "backend_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["runtime_identity"]["DESeq2"] = "1.52.1"
    _write_json(path, manifest)

    with pytest.raises(InputIntegrityError, match="manifest identity"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_rejects_reporting_effect_not_reproduced_by_coefficients(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "results.tsv"
    content = path.read_text(encoding="utf-8").replace(
        "\t10\t1\t1\t0.2", "\t10\t1.1\t1.1\t0.2"
    )
    path.write_text(content, encoding="utf-8")
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="reproduce from coefficients"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_rejects_three_prime_rounding_delta_over_half(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["input_semantics"] = "salmon_quant_dirs_three_prime"
    audit = _rounding_audit(mode="explicit_before_matrix_constructor")
    audit["max_absolute_delta"] = 0.6
    audit["per_sample"][0]["max_absolute_delta"] = 0.6
    analysis["route_observed"] = _three_prime_route(rounding_audit=audit)
    analysis["pipeline"][0]["constructor"] = "DESeq2::DESeqDataSetFromMatrix"
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="exceeds one half"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_rejects_three_prime_single_inferential_replicate(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["input_semantics"] = "salmon_quant_dirs_three_prime"
    analysis["route_observed"] = _three_prime_route(replicate_count=1)
    analysis["pipeline"][0]["constructor"] = "DESeq2::DESeqDataSetFromMatrix"
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="zero or at least two"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_rejects_full_length_single_inferential_replicate(
    tmp_path: Path,
) -> None:
    run_dir = _bundle(tmp_path)
    path = run_dir / "analysis.json"
    analysis = json.loads(path.read_text(encoding="utf-8"))
    analysis["input_semantics"] = "salmon_quant_dirs_full_length"
    analysis["route_observed"] = {
        "constructor": "DESeq2::DESeqDataSetFromTximport",
        "count_source": "txi$counts",
        "count_semantics": "salmon_estimated_counts_rounded_for_DESeq2",
        "transcript_length_offset": True,
        "gene_length_correction": True,
        "countsFromAbundance": "no",
        "dropInfReps": False,
        "inferential_replicates_imported": True,
        "inferential_replicate_state": "all",
        "inferential_replicates_per_sample": 1,
        "inferential_replicate_method": "Gibbs",
        "inferential_replicates_used_for_inference": False,
        "inferential_replicates_unused_reason": (
            "DESeq2 1.52.0 does not consume tximport infReps in this backend; "
            "they are imported only to verify declared evidence."
        ),
        "rounding_disclosure": (
            "DESeq2::DESeqDataSetFromTximport internally calls base::round; "
            "the pre/post conversion is audited here."
        ),
        "rounding_audit": _rounding_audit(mode="deseq2_constructor"),
    }
    analysis["pipeline"][0]["constructor"] = "DESeq2::DESeqDataSetFromTximport"
    _write_json(path, analysis)
    _refresh_manifest(run_dir)

    with pytest.raises(InputIntegrityError, match="zero or at least two"):
        summarize_run(run_dir)


@pytest.mark.unit
def test_deseq2_rejects_edgeR_display_sidecar(tmp_path: Path) -> None:
    run_dir = _bundle(tmp_path)
    (run_dir / "display").mkdir()

    with pytest.raises(InputIntegrityError, match="cannot contain an edgeR display"):
        summarize_run(run_dir)
