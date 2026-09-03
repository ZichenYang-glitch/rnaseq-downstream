#!/usr/bin/env python3
"""Run and report the locked compcodeR negative-binomial FDR/TPR gate."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
from pathlib import Path
import platform
import statistics
import tempfile
from typing import Any, Iterator

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

BENCHMARK_ID = "compcoder-edger-ql-nb-fdr-tpr-v1"
REPLICATES = tuple(range(1, 11))
NOMINAL_FDR = 0.05
MAX_MEAN_REPLICATE_FDP = 0.065
MAX_REPLICATE_FDP = 0.10
MIN_MEAN_TPR = 0.50
MIN_REPLICATE_TPR = 0.45


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument("--rscript")
    parser.add_argument("--workspace", type=Path)
    return parser.parse_args()


def _workspace(path: Path | None) -> Iterator[Path]:
    if path is None:
        return tempfile.TemporaryDirectory(prefix="rnaseq-compcoder-gate-")  # type: ignore[return-value]
    path.mkdir(parents=True, exist_ok=False)
    return nullcontext(path)  # type: ignore[return-value]


def _private_kernel(
    counts_path: Path,
    metadata_path: Path,
    output_dir: Path,
    *,
    rscript: Path,
    r_library: Path,
) -> Any:
    try:
        from rnaseq_downstream.edger_backend import _run_edger_ql_benchmark_kernel
    except (ImportError, AttributeError) as exc:
        raise BenchmarkError(
            "The private edgeR benchmark-kernel adapter is unavailable."
        ) from exc
    return _run_edger_ql_benchmark_kernel(
        counts_path,
        metadata_path,
        output_dir,
        design={
            "intercept": True,
            "terms": ["condition"],
            "variables": {
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                }
            },
        },
        contrasts=[
            {
                "contrast_id": "treated_vs_control",
                "weights": {"conditiontreated": 1.0},
                "lfc_threshold": 0.0,
            }
        ],
        rscript=str(rscript),
        r_library=r_library,
    )


def _score_replicate(
    replicate: int,
    fixture: Path,
    output: Path,
) -> dict[str, object]:
    truth = read_tsv(fixture / "truth.tsv")
    results = read_tsv(output / "results.tsv")
    if set(truth) != set(results):
        raise BenchmarkError(
            f"Replicate {replicate} backend result gene set differs from truth."
        )

    true_positives = 0
    false_positives = 0
    true_differential = 0
    tested = 0
    filtered = 0
    status_counts: dict[str, int] = {}
    for gene_id, truth_row in truth.items():
        truth_value = truth_row["differential_expression"]
        if truth_value not in {"0", "1"}:
            raise BenchmarkError(
                f"Replicate {replicate} has invalid truth for gene {gene_id}."
            )
        is_differential = truth_value == "1"
        true_differential += int(is_differential)
        result = results[gene_id]
        if result.get("contrast_id") != "treated_vs_control":
            raise BenchmarkError(
                f"Replicate {replicate} has an unexpected contrast identifier."
            )
        status = result.get("status", "")
        status_counts[status] = status_counts.get(status, 0) + 1
        called = False
        if status == "tested":
            tested += 1
            fdr = finite_float(result["FDR"], field="FDR", gene_id=gene_id)
            if not 0 <= fdr <= 1:
                raise BenchmarkError(
                    f"Replicate {replicate} FDR is outside [0, 1] for {gene_id}."
                )
            called = fdr <= NOMINAL_FDR
        elif status == "filtered":
            filtered += 1
        else:
            raise BenchmarkError(
                f"Replicate {replicate} has unsupported gene status {status!r}."
            )
        if called and is_differential:
            true_positives += 1
        elif called:
            false_positives += 1

    discoveries = true_positives + false_positives
    false_discovery_proportion = false_positives / discoveries if discoveries else 0.0
    tpr = true_positives / true_differential
    return {
        "replicate": replicate,
        "seed": 5000 + replicate,
        "tested": tested,
        "filtered": filtered,
        "status_counts": status_counts,
        "true_differential": true_differential,
        "discoveries": discoveries,
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_discovery_proportion": false_discovery_proportion,
        "tpr": tpr,
    }


def _assert_frozen_fixture(replicate: int, fixture: Path) -> dict[str, Any]:
    """Fail before backend execution unless one replicate is exactly in scope."""

    metadata = read_json(fixture / "fixture.json")
    expected = {
        "fixture_id": (f"compcodeR-1.48.0-p0-nb-balanced-replicate-{replicate:02d}"),
        "compcodeR_version": "1.48.0",
        "seed": 5000 + replicate,
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
    }
    if metadata != expected:
        raise BenchmarkError(
            f"Replicate {replicate} fixture metadata differs from the frozen scenario."
        )

    counts = read_tsv(fixture / "counts.tsv")
    sample_metadata = read_tsv(fixture / "metadata.tsv", key="sample_id")
    truth = read_tsv(fixture / "truth.tsv")
    conditions = [row["condition"] for row in sample_metadata.values()]
    differential_rows = [
        row for row in truth.values() if row["differential_expression"] == "1"
    ]
    directions = [
        finite_float(
            row["true_log2_fold_change"],
            field="true_log2_fold_change",
            gene_id=f"replicate-{replicate}",
        )
        for row in differential_rows
    ]
    observed = {
        "count_genes": len(counts),
        "truth_genes": len(truth),
        "samples": len(sample_metadata),
        "control_samples": conditions.count("control"),
        "treated_samples": conditions.count("treated"),
        "true_differential": len(differential_rows),
        "true_upregulated": sum(value > 0 for value in directions),
        "true_downregulated": sum(value < 0 for value in directions),
    }
    required = {
        "count_genes": 5000,
        "truth_genes": 5000,
        "samples": 12,
        "control_samples": 6,
        "treated_samples": 6,
        "true_differential": 500,
        "true_upregulated": 250,
        "true_downregulated": 250,
    }
    if observed != required:
        raise BenchmarkError(
            f"Replicate {replicate} fixture dimensions/truth differ: {observed!r}"
        )
    return metadata


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    rscript = rscript_from_runtime(arguments.rscript)
    r_library = arguments.r_library.resolve(strict=True)
    generator = PROJECT_ROOT / "scripts/benchmark/generate_compcoder_fixture.R"
    runner_script = Path(__file__).resolve()
    common_script = runner_script.with_name("common.py")
    schema_path = runner_script.with_name("benchmark-report-v1.schema.json")
    adapter_script = PROJECT_ROOT / "rnaseq_downstream/edger_backend.py"
    contract_script = PROJECT_ROOT / "rnaseq_downstream/analysis_contract.py"
    toolkit_script = PROJECT_ROOT / "rnaseq_downstream/r_scripts/edger_ql.R"
    pathway_script = PROJECT_ROOT / "rnaseq_downstream/r_scripts/pathway_tests.R"
    compcoder_description = r_library / "compcodeR/DESCRIPTION"
    environment_evidence = [
        PROJECT_ROOT / "conda-lock.yml",
        PROJECT_ROOT / "renv.lock",
        PROJECT_ROOT / "environment.p0.yml",
        PROJECT_ROOT / "environment/r-sources.lock",
        PROJECT_ROOT / "environment/verify.R",
    ]
    implementation = [
        file_evidence(runner_script),
        file_evidence(common_script),
        file_evidence(schema_path),
        file_evidence(generator),
        file_evidence(adapter_script),
        file_evidence(contract_script),
        file_evidence(toolkit_script),
        file_evidence(pathway_script),
        *(file_evidence(path) for path in environment_evidence),
    ]
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "benchmark_id": BENCHMARK_ID,
        "status": "fail",
        "runtime": {
            "python": platform.python_version(),
            "rscript_executable": rscript.name,
            "rscript_sha256": sha256_file(rscript),
            "r_library_identity": "exact_package_versions_recorded_below",
        },
        "implementation": implementation,
        "inputs": [file_evidence(compcoder_description, name="compcodeR/DESCRIPTION")],
        "artifacts": [],
        "thresholds": {
            "nominal_bh_fdr": NOMINAL_FDR,
            "maximum_mean_replicate_fdp": MAX_MEAN_REPLICATE_FDP,
            "maximum_replicate_fdp": MAX_REPLICATE_FDP,
            "minimum_mean_tpr": MIN_MEAN_TPR,
            "minimum_replicate_tpr": MIN_REPLICATE_TPR,
            "selection_context": "frozen_after_disclosed_exploratory_run",
        },
        "metrics": {"scope": "backend_kernel"},
        "assertions": [
            {
                "name": "benchmark_execution",
                "passed": False,
                "observed": "not_completed",
                "requirement": "all ten simulation replicates complete",
            }
        ],
    }
    try:
        with _workspace(arguments.workspace) as workspace_value:
            workspace = Path(workspace_value)
            replicate_metrics: list[dict[str, object]] = []
            inputs: list[dict[str, object]] = []
            artifacts: list[dict[str, object]] = []
            toolkit_backend: object = None
            compcoder_versions: set[str] = set()
            runtime_identities: list[dict[str, Any]] = []
            for replicate in REPLICATES:
                fixture = workspace / f"replicate-{replicate:02d}" / "fixture"
                output = workspace / f"replicate-{replicate:02d}" / "toolkit"
                run_rscript(
                    rscript,
                    generator,
                    [replicate, fixture],
                    r_library=r_library,
                )
                fixture_metadata = _assert_frozen_fixture(replicate, fixture)
                _private_kernel(
                    fixture / "counts.tsv",
                    fixture / "metadata.tsv",
                    output,
                    rscript=rscript,
                    r_library=r_library,
                )
                label = f"replicate-{replicate:02d}"
                for name in ("counts.tsv", "metadata.tsv", "truth.tsv", "fixture.json"):
                    inputs.append(file_evidence(fixture / name, name=f"{label}/{name}"))
                analysis = read_json(output / "analysis.json")
                compcoder_versions.add(str(fixture_metadata["compcodeR_version"]))
                toolkit_backend = analysis.get("backend")
                runtime_identity = analysis.get("runtime_identity")
                if not isinstance(runtime_identity, dict):
                    raise BenchmarkError(
                        f"Replicate {replicate} lacks backend runtime identity."
                    )
                runtime_identities.append(runtime_identity)
                if (
                    analysis.get("execution_scope") != "backend_kernel_only"
                    or analysis.get("input_semantics")
                    != "benchmark_kernel_integer_counts"
                ):
                    raise BenchmarkError(
                        f"Replicate {replicate} did not use the private benchmark kernel."
                    )
                for name in ("results.tsv", "coefficients.tsv", "design.tsv"):
                    artifacts.append(
                        file_evidence(output / name, name=f"{label}/{name}")
                    )
                replicate_metrics.append(_score_replicate(replicate, fixture, output))

            total_true_positives = sum(
                int(item["true_positives"]) for item in replicate_metrics
            )
            total_false_positives = sum(
                int(item["false_positives"]) for item in replicate_metrics
            )
            total_discoveries = total_true_positives + total_false_positives
            pooled_fdp = (
                total_false_positives / total_discoveries if total_discoveries else 0.0
            )
            replicate_fdps = [
                float(item["false_discovery_proportion"]) for item in replicate_metrics
            ]
            replicate_tprs = [float(item["tpr"]) for item in replicate_metrics]
            mean_replicate_fdp = statistics.fmean(replicate_fdps)
            mean_tpr = statistics.fmean(replicate_tprs)
            maximum_fdp = max(replicate_fdps)
            minimum_tpr = min(replicate_tprs)
            expected_runtime = {
                "R": "4.6.1",
                "Bioconductor": "3.23",
                "BiocVersion_package": "3.23.1",
                "edgeR": "4.10.0",
                "tximport": "1.40.0",
                "limma": "3.68.0",
            }
            assertions: list[dict[str, object]] = [
                {
                    "name": "benchmark_execution",
                    "passed": True,
                    "observed": len(replicate_metrics),
                    "requirement": "all ten simulation replicates complete",
                },
                {
                    "name": "fixture_identity_and_dimensions",
                    "passed": len(replicate_metrics) == 10,
                    "observed": {
                        "validated_replicates": len(replicate_metrics),
                        "genes_per_replicate": 5000,
                        "true_differential_per_replicate": 500,
                        "samples_per_replicate": 12,
                    },
                    "requirement": "ten exact 5000-gene/500-DE balanced frozen fixtures",
                },
                {
                    "name": "compcoder_version_consistency",
                    "passed": compcoder_versions == {"1.48.0"},
                    "observed": sorted(compcoder_versions),
                    "requirement": "every replicate is generated by compcodeR 1.48.0",
                },
                {
                    "name": "locked_backend_identity",
                    "passed": (
                        toolkit_backend == "edgeR_QL"
                        and all(item == expected_runtime for item in runtime_identities)
                    ),
                    "observed": {
                        "backend": toolkit_backend,
                        "runtime_identities": runtime_identities,
                    },
                    "requirement": "all replicates use the exact locked edgeR backend",
                },
                {
                    "name": "mean_replicate_fdp_gate",
                    "passed": mean_replicate_fdp <= MAX_MEAN_REPLICATE_FDP,
                    "observed": mean_replicate_fdp,
                    "requirement": "mean replicate FDP (Monte Carlo FDR estimate) <= 0.065",
                },
                {
                    "name": "worst_replicate_fdp_gate",
                    "passed": maximum_fdp <= MAX_REPLICATE_FDP,
                    "observed": maximum_fdp,
                    "requirement": "maximum replicate FDP <= 0.10",
                },
                {
                    "name": "mean_tpr_gate",
                    "passed": mean_tpr >= MIN_MEAN_TPR,
                    "observed": mean_tpr,
                    "requirement": "mean TPR >= 0.50",
                },
                {
                    "name": "worst_replicate_tpr_gate",
                    "passed": minimum_tpr >= MIN_REPLICATE_TPR,
                    "observed": minimum_tpr,
                    "requirement": "minimum replicate TPR >= 0.45",
                },
            ]
            report["inputs"] = inputs
            report["artifacts"] = artifacts
            report["runtime"].update(
                {
                    "compcodeR": (
                        next(iter(compcoder_versions))
                        if len(compcoder_versions) == 1
                        else sorted(compcoder_versions)
                    ),
                    "toolkit_backend": toolkit_backend,
                }
            )
            report["metrics"] = {
                "scope": "backend_kernel",
                "replicate_count": len(replicate_metrics),
                "total_discoveries": total_discoveries,
                "total_true_positives": total_true_positives,
                "total_false_positives": total_false_positives,
                "mean_replicate_fdp": mean_replicate_fdp,
                "discovery_weighted_pooled_fdp": pooled_fdp,
                "maximum_replicate_fdp": maximum_fdp,
                "mean_tpr": mean_tpr,
                "minimum_replicate_tpr": minimum_tpr,
                "replicates": replicate_metrics,
            }
            report["assertions"] = assertions
            report["status"] = (
                "pass" if all(item["passed"] for item in assertions) else "fail"
            )
    except Exception as exc:
        report["assertions"] = [
            {
                "name": "benchmark_execution",
                "passed": False,
                "observed": {"error_type": type(exc).__name__, "message": str(exc)},
                "requirement": "all ten simulation replicates complete",
            }
        ]
    assert_report_shape(report, benchmark_id=BENCHMARK_ID)
    write_json(arguments.report, report)
    return report


def main() -> int:
    report = run(_arguments())
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
