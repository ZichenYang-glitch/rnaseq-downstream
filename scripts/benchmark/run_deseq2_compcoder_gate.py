#!/usr/bin/env python3
"""Run disclosed exploratory or held-out DESeq2 compcodeR simulations."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from pathlib import Path
import os
import platform
import statistics
import subprocess
import sys
import tempfile
from typing import Any, Iterator, Mapping, Sequence

from common import (
    BenchmarkError,
    PROJECT_ROOT,
    SCHEMA_VERSION,
    assert_report_shape,
    file_evidence,
    finite_float,
    read_json,
    read_tsv,
    rscript_from_runtime,
    run_rscript,
    sha256_file,
    write_json,
)

EXPLORATORY_BENCHMARK_ID = "compcoder-deseq2-nb-exploratory-v1"
CERTIFICATION_BENCHMARK_ID = "compcoder-deseq2-nb-fdr-tpr-v1"
REPLICATES = tuple(range(1, 21))
SEED_BASES = {"exploratory": 61000, "certification": 62000}
WORKERS = 2
NOMINAL_FDR = 0.05

# These limits are intentionally declared in code before the exploratory seed
# grid is inspected.  They are coarse regression bounds around a nominal 0.05
# BH target and the established conservative behavior of DESeq2 in this fixed
# 6-vs-6 compcodeR scenario.  The held-out certification grid is never used to
# select or revise them.
MAX_MEAN_REPLICATE_FDP = 0.065
MAX_REPLICATE_FDP = 0.12
MIN_MEAN_TPR = 0.45
MIN_REPLICATE_TPR = 0.35

CONTRAST_ID = "treated_vs_control"
EXPECTED_RUNTIME = {
    "R": "4.6.1",
    "Bioconductor": "3.23",
    "BiocVersion_package": "3.23.1",
    "DESeq2": "1.52.0",
    "apeglm": "1.34.0",
    "tximport": "1.40.0",
}
EXPECTED_RESULT_SEMANTICS = {
    "statistic_type": "Wald",
    "statistic_hypothesis": "contrast_equals_zero",
    "fdr_basis": "contrast_pvalue_BH_after_independent_filtering",
    "test_method": "DESeq2::results_Wald",
    "lfc_threshold": "0",
    "shrinkage_method": "none",
}
EXPECTED_STATUSES = {"filtered", "not_tested", "failed", "tested"}
FIXTURE_MEMBERS = ("counts.tsv", "metadata.tsv", "truth.tsv", "fixture.json")
REQUEST_MEMBERS = (
    "reference.fa",
    "counts.manifest.json",
    "input-request.json",
    "analysis-request.json",
)
NUMERIC_ARTIFACTS = ("results.tsv", "coefficients.tsv", "design.tsv")


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode", choices=("exploratory", "certification"), required=True
    )
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument("--rscript")
    parser.add_argument("--workspace", type=Path)
    return parser.parse_args()


def _workspace(path: Path | None, *, mode: str) -> Iterator[Path]:
    if path is None:
        return tempfile.TemporaryDirectory(  # type: ignore[return-value]
            prefix=f"rnaseq-deseq2-{mode}-"
        )
    path.mkdir(parents=True, exist_ok=False)
    return nullcontext(path)  # type: ignore[return-value]


def _benchmark_id(mode: str) -> str:
    if mode == "exploratory":
        return EXPLORATORY_BENCHMARK_ID
    if mode == "certification":
        return CERTIFICATION_BENCHMARK_ID
    raise BenchmarkError(f"Unsupported simulation mode: {mode!r}")


def _thresholds(*, mode: str, exploratory_report: Path | None) -> dict[str, object]:
    values: dict[str, object] = {
        "nominal_bh_fdr": NOMINAL_FDR,
        "maximum_mean_replicate_fdp": MAX_MEAN_REPLICATE_FDP,
        "maximum_replicate_fdp": MAX_REPLICATE_FDP,
        "minimum_mean_tpr": MIN_MEAN_TPR,
        "minimum_replicate_tpr": MIN_REPLICATE_TPR,
        "threshold_basis": (
            "nominal_0.05_BH_plus_finite_simulation_tolerance_and_"
            "conservative_DESeq2_power_regression_floor"
        ),
        "held_out_policy": (
            "certification_seeds_are_disjoint_and_thresholds_are_not_revised_"
            "after_certification_results_are_observed"
        ),
        "applied_to_this_report": mode == "certification",
    }
    if mode == "exploratory":
        values["selection_context"] = (
            "candidate_limits_recorded_before_held_out_certification"
        )
    else:
        if exploratory_report is None:
            raise BenchmarkError("Certification requires an exploratory report.")
        values.update(
            {
                "selection_context": (
                    "frozen_after_disclosed_exploratory_run_before_held_out_run"
                ),
                "exploratory_report_sha256": sha256_file(exploratory_report),
            }
        )
    return values


def _strict_json_line(stdout: str, *, stage: str) -> dict[str, Any]:
    if not stdout.endswith("\n") or stdout.count("\n") != 1:
        raise BenchmarkError(
            f"Public CLI {stage} did not emit exactly one newline-terminated JSON line."
        )

    def reject_constant(value: str) -> None:
        raise ValueError(f"non-finite constant {value}")

    def unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        document: dict[str, Any] = {}
        for key, value in pairs:
            if key in document:
                raise ValueError(f"duplicate key {key!r}")
            document[key] = value
        return document

    try:
        import json

        document = json.loads(
            stdout,
            parse_constant=reject_constant,
            object_pairs_hook=unique_pairs,
        )
    except (ValueError, TypeError) as error:
        raise BenchmarkError(
            f"Public CLI {stage} emitted invalid JSON: {error}"
        ) from error
    if not isinstance(document, dict):
        raise BenchmarkError(f"Public CLI {stage} response is not an object.")
    return document


def _run_cli(
    arguments: Sequence[str],
    *,
    cwd: Path,
    stage: str,
    timeout_seconds: int = 300,
) -> dict[str, Any]:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [sys.executable, "-m", "rnaseq_downstream", *arguments],
        cwd=cwd,
        env=environment,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="strict",
        timeout=timeout_seconds,
    )
    response = _strict_json_line(completed.stdout, stage=stage)
    if completed.returncode != 0 or response.get("status") != "success":
        raise BenchmarkError(
            f"Public CLI {stage} failed (exit {completed.returncode}): {response!r}; "
            f"stderr={completed.stderr[-4000:]!r}"
        )
    if response.get("errors") != []:
        raise BenchmarkError(f"Public CLI {stage} success response contains errors.")
    return response


def _sample_ids(metadata_path: Path) -> list[str]:
    metadata = read_tsv(metadata_path, key="sample_id")
    samples = list(metadata)
    conditions = [row.get("condition") for row in metadata.values()]
    if conditions != ["control"] * 6 + ["treated"] * 6:
        raise BenchmarkError(
            "Simulation metadata does not have the frozen 6-vs-6 order."
        )
    return samples


def _write_public_requests(job_root: Path, fixture: Path) -> tuple[Path, Path]:
    samples = _sample_ids(fixture / "metadata.tsv")
    reference = fixture / "reference.fa"
    reference.write_bytes(b">compcodeR_synthetic_reference\nACGT\n")
    manifest = fixture / "counts.manifest.json"
    write_json(
        manifest,
        {
            "schema_version": "1.0",
            "artifact_type": "featurecounts_integer_matrix",
            "artifact": {
                "path": "counts.tsv",
                "sha256": sha256_file(fixture / "counts.tsv"),
            },
            "gene_id_column": "gene_id",
            "display_columns": [],
            "sample_columns": samples,
            "producer": {"name": "featureCounts", "version": "2.0.6"},
            "reference": {
                "name": "compcodeR_synthetic_reference",
                "version": "1",
                "source": "reference.fa",
                "sha256": sha256_file(reference),
            },
        },
    )
    input_request = job_root / "input-request.json"
    write_json(
        input_request,
        {
            "schema_version": "1.0",
            "input_semantics": "featurecounts_integer",
            "metadata": {
                "path": "fixture/metadata.tsv",
                "sample_id_column": "sample_id",
            },
            "producer": {"name": "featureCounts", "version": "2.0.6"},
            "reference": {
                "name": "compcodeR_synthetic_reference",
                "version": "1",
                "source": "fixture/reference.fa",
            },
            "gene_id": {"strip_version": False},
            "samples": [{"sample_id": sample} for sample in samples],
            "featurecounts": {
                "layout": "combined_matrix",
                "matrix": "fixture/counts.tsv",
                "manifest": "fixture/counts.manifest.json",
            },
        },
    )
    analysis_request = job_root / "analysis-request.json"
    write_json(
        analysis_request,
        {
            "schema_version": "1.2",
            "backend": "deseq2",
            "validated_input_bundle": "validated-input",
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
                    "contrast_id": CONTRAST_ID,
                    "weights": {"conditiontreated": 1},
                    "lfc_threshold": 0,
                }
            ],
            "deseq2": {"test_mode": "wald", "shrinkage": "none"},
        },
    )
    return input_request, analysis_request


def _assert_frozen_fixture(mode: str, replicate: int, fixture: Path) -> dict[str, Any]:
    metadata = read_json(fixture / "fixture.json")
    expected = {
        "fixture_id": (
            f"compcodeR-1.48.0-p1-deseq2-nb-balanced-{mode}-replicate-{replicate:02d}"
        ),
        "mode": mode,
        "compcodeR_version": "1.48.0",
        "seed": SEED_BASES[mode] + replicate,
        "replicate": replicate,
        "n_genes": 5000,
        "samples_per_condition": 6,
        "n_differential": 500,
        "n_upregulated": 250,
        "n_downregulated": 250,
        "effect_size": 2,
        "fraction_upregulated": 0.5,
        "filter_threshold_total": 0,
        "outlier_probabilities": {
            "random_high": 0,
            "random_low": 0,
            "single_high": 0,
            "single_low": 0,
        },
        "design": "~ condition",
        "factor_levels": {"condition": ["control", "treated"]},
        "coefficient": "conditiontreated",
        "input_route": (
            "featurecounts_integer_combined_matrix_self_attested_simulation"
        ),
    }
    if metadata != expected:
        raise BenchmarkError(
            f"{mode} replicate {replicate} differs from the frozen fixture metadata."
        )
    counts = read_tsv(fixture / "counts.tsv")
    truth = read_tsv(fixture / "truth.tsv")
    sample_metadata = read_tsv(fixture / "metadata.tsv", key="sample_id")
    differential = [
        row for row in truth.values() if row.get("differential_expression") == "1"
    ]
    directions = [
        finite_float(
            row["true_log2_fold_change"],
            field="true_log2_fold_change",
            gene_id=f"{mode}-replicate-{replicate}",
        )
        for row in differential
    ]
    observed = {
        "count_genes": len(counts),
        "truth_genes": len(truth),
        "samples": len(sample_metadata),
        "true_differential": len(differential),
        "true_upregulated": sum(value > 0 for value in directions),
        "true_downregulated": sum(value < 0 for value in directions),
    }
    required = {
        "count_genes": 5000,
        "truth_genes": 5000,
        "samples": 12,
        "true_differential": 500,
        "true_upregulated": 250,
        "true_downregulated": 250,
    }
    if observed != required:
        raise BenchmarkError(
            f"{mode} replicate {replicate} dimensions/truth differ: {observed!r}"
        )
    return metadata


def _score_replicate(
    mode: str,
    replicate: int,
    fixture: Path,
    output: Path,
    *,
    cli_stages: Mapping[str, str],
) -> dict[str, object]:
    truth = read_tsv(fixture / "truth.tsv")
    results = read_tsv(output / "results.tsv")
    if set(truth) != set(results):
        raise BenchmarkError(
            f"{mode} replicate {replicate} result gene set differs from truth."
        )

    true_positives = 0
    false_positives = 0
    true_differential = 0
    status_counts: Counter[str] = Counter()
    reason_counts: Counter[str] = Counter()
    unexpected_statuses: set[str] = set()
    for gene_id, truth_row in truth.items():
        truth_value = truth_row.get("differential_expression")
        if truth_value not in {"0", "1"}:
            raise BenchmarkError(
                f"{mode} replicate {replicate} has invalid truth for {gene_id}."
            )
        is_differential = truth_value == "1"
        true_differential += int(is_differential)
        result = results[gene_id]
        if result.get("contrast_id") != CONTRAST_ID:
            raise BenchmarkError(
                f"{mode} replicate {replicate} has an unexpected contrast."
            )
        for field, expected in EXPECTED_RESULT_SEMANTICS.items():
            observed = result.get(field)
            if field == "lfc_threshold":
                try:
                    matches = float(observed or "nan") == 0.0
                except ValueError:
                    matches = False
            else:
                matches = observed == expected
            if not matches:
                raise BenchmarkError(
                    f"{mode} replicate {replicate} has wrong {field} for {gene_id}: "
                    f"{observed!r}"
                )
        status = result.get("status", "")
        reason = result.get("status_reason", "")
        status_counts[status] += 1
        reason_counts[reason] += 1
        unexpected_statuses.update({status} - EXPECTED_STATUSES)
        called = False
        if status == "tested":
            fdr = finite_float(result["FDR"], field="FDR", gene_id=gene_id)
            if not 0 <= fdr <= 1:
                raise BenchmarkError(
                    f"{mode} replicate {replicate} FDR is outside [0, 1] for {gene_id}."
                )
            called = fdr <= NOMINAL_FDR
        elif status in {"filtered", "not_tested", "failed"}:
            if result.get("FDR") not in {"", None}:
                raise BenchmarkError(
                    f"{mode} replicate {replicate} non-tested row has FDR for {gene_id}."
                )
        if called and is_differential:
            true_positives += 1
        elif called:
            false_positives += 1

    if unexpected_statuses:
        raise BenchmarkError(
            f"{mode} replicate {replicate} has unsupported statuses: "
            f"{sorted(unexpected_statuses)!r}"
        )
    discoveries = true_positives + false_positives
    return {
        "replicate": replicate,
        "seed": SEED_BASES[mode] + replicate,
        "cli_stages": dict(cli_stages),
        "status_counts": dict(sorted(status_counts.items())),
        "status_reason_counts": dict(sorted(reason_counts.items())),
        "true_differential": true_differential,
        "discoveries": discoveries,
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_discovery_proportion": (
            false_positives / discoveries if discoveries else 0.0
        ),
        "tpr": true_positives / true_differential,
    }


def _run_public_chain(
    job_root: Path,
    *,
    rscript: Path,
    r_library: Path,
) -> tuple[Path, dict[str, str]]:
    input_request = job_root / "input-request.json"
    analysis_request = job_root / "analysis-request.json"
    validated = job_root / "validated-input"
    output = job_root / "toolkit"
    inspect = _run_cli(
        ["inspect", "--request", input_request.name],
        cwd=job_root,
        stage="inspect",
    )
    validate = _run_cli(
        [
            "validate",
            "--request",
            input_request.name,
            "--output-dir",
            validated.name,
        ],
        cwd=job_root,
        stage="validate",
    )
    run = _run_cli(
        [
            "run",
            "--request",
            analysis_request.name,
            "--output-dir",
            output.name,
            "--rscript",
            str(rscript),
            "--r-library",
            str(r_library),
        ],
        cwd=job_root,
        stage="run",
    )
    summarize = _run_cli(
        ["summarize", "--run-dir", output.name],
        cwd=job_root,
        stage="summarize",
    )
    run_data = run.get("data")
    summary_data = summarize.get("data")
    if not isinstance(run_data, Mapping) or run_data.get("backend") != "DESeq2":
        raise BenchmarkError("Public CLI run did not identify the DESeq2 backend.")
    scope = run_data.get("scope")
    if not isinstance(scope, Mapping) or any(
        scope.get(field) != "passed"
        for field in ("input_semantics", "design", "backend")
    ):
        raise BenchmarkError(
            "Public CLI run did not report all validation stages passed."
        )
    if (
        not isinstance(summary_data, Mapping)
        or summary_data.get("status") != "verified_complete"
        or summary_data.get("backend") != "DESeq2"
    ):
        raise BenchmarkError("Public CLI summarize did not verify the DESeq2 bundle.")
    if inspect.get("status") != "success" or validate.get("status") != "success":
        raise BenchmarkError(
            "Public CLI inspect/validate did not complete successfully."
        )
    return output, {
        "inspect": "success",
        "validate": "success",
        "run": "success",
        "summarize": "verified_complete",
    }


def _run_replicate(
    replicate: int,
    *,
    mode: str,
    workspace: Path,
    generator: Path,
    rscript: Path,
    r_library: Path,
) -> dict[str, Any]:
    job_root = workspace / f"{mode}-replicate-{replicate:02d}"
    fixture = job_root / "fixture"
    run_rscript(
        rscript,
        generator,
        [mode, replicate, fixture],
        r_library=r_library,
    )
    fixture_metadata = _assert_frozen_fixture(mode, replicate, fixture)
    _write_public_requests(job_root, fixture)
    output, cli_stages = _run_public_chain(
        job_root, rscript=rscript, r_library=r_library
    )
    analysis = read_json(output / "analysis.json")
    runtime_identity = analysis.get("runtime_identity")
    if not isinstance(runtime_identity, dict):
        raise BenchmarkError(f"{mode} replicate {replicate} lacks runtime identity.")
    if (
        analysis.get("backend") != "DESeq2"
        or analysis.get("execution_scope") != "validated_p1_deseq2_input"
        or analysis.get("input_semantics") != "featurecounts_integer"
    ):
        raise BenchmarkError(
            f"{mode} replicate {replicate} did not traverse the public DESeq2 input route."
        )
    test = analysis.get("test")
    if not isinstance(test, Mapping) or test.get("mode") != "wald":
        raise BenchmarkError(f"{mode} replicate {replicate} did not use Wald testing.")
    score = _score_replicate(
        mode,
        replicate,
        fixture,
        output,
        cli_stages=cli_stages,
    )
    return {
        "replicate": replicate,
        "fixture_metadata": fixture_metadata,
        "runtime_identity": runtime_identity,
        "score": score,
        "input_files": [
            file_evidence(fixture / name, name=f"replicate-{replicate:02d}/{name}")
            for name in FIXTURE_MEMBERS
        ]
        + [
            file_evidence(
                fixture / name
                if name in {"reference.fa", "counts.manifest.json"}
                else job_root / name,
                name=f"replicate-{replicate:02d}/{name}",
            )
            for name in REQUEST_MEMBERS
        ],
        "artifacts": [
            file_evidence(output / name, name=f"replicate-{replicate:02d}/{name}")
            for name in NUMERIC_ARTIFACTS
        ],
    }


def _aggregate(
    mode: str,
    replicate_metrics: Sequence[Mapping[str, object]],
    *,
    runtime_identities: Sequence[Mapping[str, object]],
    fixture_determinism: Mapping[str, bool],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    fdps = [float(item["false_discovery_proportion"]) for item in replicate_metrics]
    tprs = [float(item["tpr"]) for item in replicate_metrics]
    true_positives = sum(int(item["true_positives"]) for item in replicate_metrics)
    false_positives = sum(int(item["false_positives"]) for item in replicate_metrics)
    discoveries = true_positives + false_positives
    mean_fdp = statistics.fmean(fdps)
    maximum_fdp = max(fdps)
    mean_tpr = statistics.fmean(tprs)
    minimum_tpr = min(tprs)
    all_cli_stages = all(
        item.get("cli_stages")
        == {
            "inspect": "success",
            "validate": "success",
            "run": "success",
            "summarize": "verified_complete",
        }
        for item in replicate_metrics
    )
    assertions: list[dict[str, object]] = [
        {
            "name": "benchmark_execution",
            "passed": len(replicate_metrics) == len(REPLICATES),
            "observed": len(replicate_metrics),
            "requirement": "all 20 frozen simulation replicates complete",
        },
        {
            "name": "public_cli_chain",
            "passed": all_cli_stages,
            "observed": {
                "replicates": len(replicate_metrics),
                "stages": ["inspect", "validate", "run", "summarize"],
            },
            "requirement": (
                "every replicate traverses the public inspect-to-summarize chain"
            ),
        },
        {
            "name": "locked_backend_identity",
            "passed": (
                len(runtime_identities) == len(REPLICATES)
                and all(dict(item) == EXPECTED_RUNTIME for item in runtime_identities)
            ),
            "observed": {
                "expected": EXPECTED_RUNTIME,
                "replicate_identities_equal": all(
                    dict(item) == EXPECTED_RUNTIME for item in runtime_identities
                ),
            },
            "requirement": (
                "all replicates use R 4.6.1, Bioconductor 3.23, and DESeq2 1.52.0"
            ),
        },
        {
            "name": "fixture_byte_determinism",
            "passed": bool(fixture_determinism) and all(fixture_determinism.values()),
            "observed": dict(fixture_determinism),
            "requirement": ("a representative fixture regenerates byte-identically"),
        },
        {
            "name": "exact_fixture_grid",
            "passed": (
                len(replicate_metrics) == 20
                and [int(item["seed"]) for item in replicate_metrics]
                == [SEED_BASES[mode] + replicate for replicate in REPLICATES]
                and all(
                    int(item["true_differential"]) == 500 for item in replicate_metrics
                )
            ),
            "observed": {
                "mode": mode,
                "seed_start": SEED_BASES[mode] + 1,
                "seed_end": SEED_BASES[mode] + 20,
                "replicates": len(replicate_metrics),
                "genes_per_replicate": 5000,
                "true_differential_per_replicate": 500,
                "samples_per_condition": 6,
            },
            "requirement": (
                "exactly 20 disjoint-mode 5000-gene/500-DE balanced 6-vs-6 fixtures"
            ),
        },
    ]
    if mode == "certification":
        assertions.extend(
            [
                {
                    "name": "mean_replicate_fdp_gate",
                    "passed": mean_fdp <= MAX_MEAN_REPLICATE_FDP,
                    "observed": mean_fdp,
                    "requirement": "mean replicate FDP <= 0.065",
                },
                {
                    "name": "worst_replicate_fdp_gate",
                    "passed": maximum_fdp <= MAX_REPLICATE_FDP,
                    "observed": maximum_fdp,
                    "requirement": "maximum replicate FDP <= 0.12",
                },
                {
                    "name": "mean_tpr_gate",
                    "passed": mean_tpr >= MIN_MEAN_TPR,
                    "observed": mean_tpr,
                    "requirement": "mean TPR >= 0.45",
                },
                {
                    "name": "worst_replicate_tpr_gate",
                    "passed": minimum_tpr >= MIN_REPLICATE_TPR,
                    "observed": minimum_tpr,
                    "requirement": "minimum replicate TPR >= 0.35",
                },
            ]
        )
    metrics: dict[str, object] = {
        "scope": "public_cli_inspect_validate_run_summarize",
        "evidence_role": (
            "threshold_selection_only_not_certification"
            if mode == "exploratory"
            else "held_out_certification_gate"
        ),
        "mode": mode,
        "replicate_count": len(replicate_metrics),
        "seed_grid": [int(item["seed"]) for item in replicate_metrics],
        "total_discoveries": discoveries,
        "total_true_positives": true_positives,
        "total_false_positives": false_positives,
        "mean_replicate_fdp": mean_fdp,
        "median_replicate_fdp": statistics.median(fdps),
        "minimum_replicate_fdp": min(fdps),
        "maximum_replicate_fdp": maximum_fdp,
        "discovery_weighted_pooled_fdp_descriptive_only": (
            false_positives / discoveries if discoveries else 0.0
        ),
        "mean_tpr": mean_tpr,
        "median_tpr": statistics.median(tprs),
        "minimum_replicate_tpr": minimum_tpr,
        "maximum_replicate_tpr": max(tprs),
        "replicates": list(replicate_metrics),
        "limitations": [
            "finite gross-regression study, not proof of universal FDR control",
            "compcodeR genes are simulated independently conditional on library factors",
            "combined-matrix featureCounts manifest is self-attested simulation routing evidence, not authenticated producer provenance",
            "power applies only to this fixed 6-vs-6, 500-DE, effect-size-2 scenario",
        ],
    }
    return assertions, metrics


def _validate_exploratory_report(path: Path) -> dict[str, Any]:
    report = read_json(path)
    assert_report_shape(report, benchmark_id=EXPLORATORY_BENCHMARK_ID)
    if report.get("status") != "pass":
        raise BenchmarkError("The disclosed exploratory report is not complete.")
    metrics = report.get("metrics")
    thresholds = report.get("thresholds")
    if not isinstance(metrics, Mapping) or not isinstance(thresholds, Mapping):
        raise BenchmarkError("The disclosed exploratory report has invalid metrics.")
    if (
        metrics.get("evidence_role") != "threshold_selection_only_not_certification"
        or metrics.get("mode") != "exploratory"
        or metrics.get("seed_grid") != [61000 + item for item in REPLICATES]
        or thresholds.get("applied_to_this_report") is not False
    ):
        raise BenchmarkError(
            "The disclosed exploratory report has the wrong role/grid."
        )
    numeric = {
        "nominal_bh_fdr": NOMINAL_FDR,
        "maximum_mean_replicate_fdp": MAX_MEAN_REPLICATE_FDP,
        "maximum_replicate_fdp": MAX_REPLICATE_FDP,
        "minimum_mean_tpr": MIN_MEAN_TPR,
        "minimum_replicate_tpr": MIN_REPLICATE_TPR,
    }
    if any(thresholds.get(key) != value for key, value in numeric.items()):
        raise BenchmarkError("The disclosed candidate thresholds differ from code.")
    forbidden = {
        "mean_replicate_fdp_gate",
        "worst_replicate_fdp_gate",
        "mean_tpr_gate",
        "worst_replicate_tpr_gate",
    }
    assertion_names = {
        item.get("name")
        for item in report.get("assertions", [])
        if isinstance(item, Mapping)
    }
    if forbidden & assertion_names:
        raise BenchmarkError(
            "The exploratory report improperly applies certification gates."
        )
    return report


def _implementation_inventory(
    *, mode: str, exploratory_report: Path | None
) -> list[dict[str, object]]:
    runner = Path(__file__).resolve()
    paths = [
        runner,
        runner.with_name("common.py"),
        runner.with_name("benchmark-report-v1.schema.json"),
        runner.with_name("generate_deseq2_compcoder_fixture.R"),
        PROJECT_ROOT / "rnaseq_downstream/cli.py",
        PROJECT_ROOT / "rnaseq_downstream/input_semantics.py",
        PROJECT_ROOT / "rnaseq_downstream/validation_bundle.py",
        PROJECT_ROOT / "rnaseq_downstream/analysis_contract_v12.py",
        PROJECT_ROOT / "rnaseq_downstream/deseq2_backend.py",
        PROJECT_ROOT / "rnaseq_downstream/run_summary.py",
        PROJECT_ROOT / "rnaseq_downstream/r_scripts/deseq2.R",
        PROJECT_ROOT / "conda-lock.yml",
        PROJECT_ROOT / "renv.lock",
        PROJECT_ROOT / "environment.p0.yml",
        PROJECT_ROOT / "environment/r-sources.lock",
        PROJECT_ROOT / "environment/verify.R",
    ]
    if mode == "certification":
        if exploratory_report is None:
            raise BenchmarkError("Certification requires an exploratory report.")
        paths.extend(
            [
                exploratory_report,
                PROJECT_ROOT / "tests/simulation/deseq2-compcoder-method.md",
            ]
        )
    names: list[str] = []
    evidence: list[dict[str, object]] = []
    for path in paths:
        name = (
            f"environment/{path.name}"
            if path.parent.name == "environment"
            else f"tests/simulation/{path.name}"
            if path.parent.name == "simulation"
            else path.name
        )
        if name in names:
            raise BenchmarkError(f"Duplicate implementation evidence name: {name}")
        names.append(name)
        evidence.append(file_evidence(path, name=name))
    return evidence


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    mode = arguments.mode
    benchmark_id = _benchmark_id(mode)
    rscript = rscript_from_runtime(arguments.rscript)
    r_library = arguments.r_library.resolve(strict=True)
    generator = (
        Path(__file__).resolve().with_name("generate_deseq2_compcoder_fixture.R")
    )
    exploratory_report = (
        PROJECT_ROOT / "tests/simulation/deseq2-compcoder-exploratory-report.json"
        if mode == "certification"
        else None
    )
    if mode == "certification":
        _validate_exploratory_report(exploratory_report)  # type: ignore[arg-type]
    implementation = _implementation_inventory(
        mode=mode, exploratory_report=exploratory_report
    )
    package_inputs = [
        file_evidence(
            r_library / "compcodeR/DESCRIPTION", name="compcodeR/DESCRIPTION"
        ),
        file_evidence(r_library / "DESeq2/DESCRIPTION", name="DESeq2/DESCRIPTION"),
    ]
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "benchmark_id": benchmark_id,
        "status": "fail",
        "runtime": {
            "python": platform.python_version(),
            "python_executable": Path(sys.executable).name,
            "rscript_executable": rscript.name,
            "rscript_sha256": sha256_file(rscript),
            "r_library_identity": "exact_package_versions_recorded_below",
            "replicate_workers": WORKERS,
        },
        "implementation": implementation,
        "inputs": package_inputs,
        "artifacts": [],
        "thresholds": _thresholds(mode=mode, exploratory_report=exploratory_report),
        "metrics": {
            "scope": "public_cli_inspect_validate_run_summarize",
            "mode": mode,
        },
        "assertions": [
            {
                "name": "benchmark_execution",
                "passed": False,
                "observed": "not_completed",
                "requirement": "all 20 frozen simulation replicates complete",
            }
        ],
    }
    try:
        if platform.python_version() != "3.12.12":
            raise BenchmarkError("The simulation gate requires locked Python 3.12.12.")
        if Path(sys.executable).resolve().parent != rscript.parent:
            raise BenchmarkError(
                "Python and Rscript are not from the same locked prefix."
            )
        with _workspace(arguments.workspace, mode=mode) as workspace_value:
            workspace = Path(workspace_value)
            worker = lambda replicate: _run_replicate(  # noqa: E731
                replicate,
                mode=mode,
                workspace=workspace,
                generator=generator,
                rscript=rscript,
                r_library=r_library,
            )
            with ThreadPoolExecutor(max_workers=WORKERS) as executor:
                jobs = list(executor.map(worker, REPLICATES))

            deterministic_fixture = workspace / "determinism-check"
            run_rscript(
                rscript,
                generator,
                [mode, 1, deterministic_fixture],
                r_library=r_library,
            )
            representative = workspace / f"{mode}-replicate-01/fixture"
            fixture_determinism = {
                name: sha256_file(representative / name)
                == sha256_file(deterministic_fixture / name)
                for name in FIXTURE_MEMBERS
            }
            replicate_metrics = [job["score"] for job in jobs]
            runtime_identities = [job["runtime_identity"] for job in jobs]
            assertions, metrics = _aggregate(
                mode,
                replicate_metrics,
                runtime_identities=runtime_identities,
                fixture_determinism=fixture_determinism,
            )
            compcoder_versions = {
                str(job["fixture_metadata"].get("compcodeR_version")) for job in jobs
            }
            assertions.insert(
                2,
                {
                    "name": "compcoder_version_consistency",
                    "passed": compcoder_versions == {"1.48.0"},
                    "observed": sorted(compcoder_versions),
                    "requirement": "every fixture is generated by compcodeR 1.48.0",
                },
            )
            report["inputs"] = [
                *package_inputs,
                *(record for job in jobs for record in job["input_files"]),
            ]
            report["artifacts"] = [
                record for job in jobs for record in job["artifacts"]
            ]
            report["runtime"].update(
                {
                    "compcodeR": (
                        next(iter(compcoder_versions))
                        if len(compcoder_versions) == 1
                        else sorted(compcoder_versions)
                    ),
                    "toolkit_backend": "DESeq2",
                    "backend_runtime_identity": EXPECTED_RUNTIME,
                }
            )
            report["metrics"] = metrics
            report["assertions"] = assertions
            report["status"] = (
                "pass" if all(bool(item["passed"]) for item in assertions) else "fail"
            )
    except Exception as error:
        report["assertions"] = [
            {
                "name": "benchmark_execution",
                "passed": False,
                "observed": {
                    "error_type": type(error).__name__,
                    "message": str(error),
                },
                "requirement": "all 20 frozen simulation replicates complete",
            }
        ]
    assert_report_shape(report, benchmark_id=benchmark_id)
    write_json(arguments.report, report)
    return report


def main() -> int:
    report = run(_arguments())
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
