#!/usr/bin/env python3
"""Run and report the locked airway same-engine edgeR oracle gate."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
import math
from pathlib import Path
import platform
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
    read_tsv_records,
    rscript_from_runtime,
    run_rscript,
    sha256_file,
    write_json,
)

BENCHMARK_ID = "airway-edger-ql-same-engine-v1"
RELATIVE_TOLERANCE = 1e-6
ABSOLUTE_TOLERANCE = 1e-10


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument("--rscript")
    parser.add_argument(
        "--workspace",
        type=Path,
        help="Retain intermediate fixture, oracle, and toolkit outputs here.",
    )
    return parser.parse_args()


def _workspace(path: Path | None) -> Iterator[Path]:
    if path is None:
        return tempfile.TemporaryDirectory(prefix="rnaseq-airway-oracle-")  # type: ignore[return-value]
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
    """Invoke the non-public raw-count adapter used only by certification tests."""

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
            "terms": ["cell", "dex"],
            "variables": {
                "cell": {
                    "type": "factor",
                    "levels": ["N052611", "N061011", "N080611", "N61311"],
                },
                "dex": {"type": "factor", "levels": ["untrt", "trt"]},
            },
        },
        contrasts=[
            {
                "contrast_id": "trt_vs_untrt",
                "weights": {"dextrt": 1.0},
                "lfc_threshold": 0.0,
            }
        ],
        rscript=str(rscript),
        r_library=r_library,
    )


def _comparison(
    oracle_directory: Path,
    toolkit_directory: Path,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    oracle_results = read_tsv(oracle_directory / "results.tsv")
    oracle_coefficients = read_tsv(oracle_directory / "coefficients.tsv")
    toolkit_rows = read_tsv(toolkit_directory / "results.tsv")
    coefficient_header, coefficient_records = read_tsv_records(
        toolkit_directory / "coefficients.tsv"
    )

    toolkit_tested: dict[str, dict[str, str]] = {}
    unexpected_contrasts: set[str] = set()
    for gene_id, row in toolkit_rows.items():
        if row.get("contrast_id") != "trt_vs_untrt":
            unexpected_contrasts.add(row.get("contrast_id", ""))
        if row.get("status") == "tested":
            toolkit_tested[gene_id] = row

    gene_sets_equal = (
        set(oracle_results) == set(toolkit_tested) and not unexpected_contrasts
    )
    required_long_columns = {"gene_id", "status", "coefficient", "estimate", "scale"}
    if set(coefficient_header) != required_long_columns:
        raise BenchmarkError(
            "Toolkit coefficients.tsv does not have the required strict long schema."
        )
    toolkit_coefficients: dict[tuple[str, str], str] = {}
    toolkit_coefficient_columns: list[str] = []
    toolkit_coefficient_genes: set[str] = set()
    observed_coefficient_keys: set[tuple[str, str]] = set()
    for row in coefficient_records:
        gene_id = row["gene_id"]
        coefficient = row["coefficient"]
        key = (gene_id, coefficient)
        if not gene_id or not coefficient or key in observed_coefficient_keys:
            raise BenchmarkError(
                "Toolkit coefficients.tsv has a blank or duplicate gene/coefficient key."
            )
        observed_coefficient_keys.add(key)
        if row["status"] == "tested":
            if row["scale"] != "natural_log":
                raise BenchmarkError(
                    "A tested toolkit coefficient is not on the declared natural-log scale."
                )
            toolkit_coefficients[key] = row["estimate"]
            toolkit_coefficient_genes.add(gene_id)
            if coefficient not in toolkit_coefficient_columns:
                toolkit_coefficient_columns.append(coefficient)
        elif row["status"] != "filtered":
            raise BenchmarkError(
                f"Toolkit coefficient row has unsupported status {row['status']!r}."
            )

    oracle_coefficient_columns = list(next(iter(oracle_coefficients.values())))
    oracle_coefficient_columns.remove("gene_id")
    coefficient_columns_equal = (
        oracle_coefficient_columns == toolkit_coefficient_columns
    )
    coefficient_gene_sets_equal = set(oracle_coefficients) == toolkit_coefficient_genes

    coefficient_violations = 0
    coefficient_max_absolute = 0.0
    coefficient_max_relative = 0.0
    if coefficient_columns_equal and coefficient_gene_sets_equal:
        for gene_id in oracle_coefficients:
            for column in oracle_coefficient_columns:
                expected = finite_float(
                    oracle_coefficients[gene_id][column],
                    field=column,
                    gene_id=gene_id,
                )
                observed = finite_float(
                    toolkit_coefficients[(gene_id, column)],
                    field=column,
                    gene_id=gene_id,
                )
                difference = abs(observed - expected)
                relative = difference / max(abs(expected), 1e-300)
                coefficient_max_absolute = max(coefficient_max_absolute, difference)
                coefficient_max_relative = max(coefficient_max_relative, relative)
                if not math.isclose(
                    observed,
                    expected,
                    rel_tol=RELATIVE_TOLERANCE,
                    abs_tol=ABSOLUTE_TOLERANCE,
                ):
                    coefficient_violations += 1

    logfc_violations = 0
    logfc_max_absolute = 0.0
    logfc_max_relative = 0.0
    if gene_sets_equal:
        for gene_id, oracle_row in oracle_results.items():
            expected = finite_float(oracle_row["logFC"], field="logFC", gene_id=gene_id)
            observed = finite_float(
                toolkit_tested[gene_id]["logFC"], field="logFC", gene_id=gene_id
            )
            difference = abs(observed - expected)
            relative = difference / max(abs(expected), 1e-300)
            logfc_max_absolute = max(logfc_max_absolute, difference)
            logfc_max_relative = max(logfc_max_relative, relative)
            if not math.isclose(
                observed,
                expected,
                rel_tol=RELATIVE_TOLERANCE,
                abs_tol=ABSOLUTE_TOLERANCE,
            ):
                logfc_violations += 1

    assertions: list[dict[str, object]] = [
        {
            "name": "tested_gene_set_parity",
            "passed": gene_sets_equal,
            "observed": {
                "oracle_tested": len(oracle_results),
                "toolkit_tested": len(toolkit_tested),
                "unexpected_contrasts": sorted(unexpected_contrasts),
            },
            "requirement": "identical tested gene identifiers and one declared contrast",
        },
        {
            "name": "coefficient_shape_parity",
            "passed": coefficient_columns_equal and coefficient_gene_sets_equal,
            "observed": {
                "oracle_columns": oracle_coefficient_columns,
                "toolkit_columns": toolkit_coefficient_columns,
                "oracle_genes": len(oracle_coefficients),
                "toolkit_genes": len(toolkit_coefficient_genes),
            },
            "requirement": "identical fitted coefficient columns and gene identifiers",
        },
        {
            "name": "coefficient_numeric_parity",
            "passed": coefficient_violations == 0,
            "observed": {
                "violations": coefficient_violations,
                "max_absolute_difference": coefficient_max_absolute,
                "max_relative_difference": coefficient_max_relative,
            },
            "requirement": "all fitted coefficients satisfy rtol=1e-6 and atol=1e-10",
        },
        {
            "name": "logfc_numeric_parity",
            "passed": logfc_violations == 0,
            "observed": {
                "violations": logfc_violations,
                "max_absolute_difference": logfc_max_absolute,
                "max_relative_difference": logfc_max_relative,
            },
            "requirement": "all tested logFC values satisfy rtol=1e-6 and atol=1e-10",
        },
    ]
    metrics: dict[str, object] = {
        "scope": "backend_kernel",
        "input_gene_count": len(toolkit_rows),
        "tested_gene_count": len(toolkit_tested),
        "filtered_gene_count": len(toolkit_rows) - len(toolkit_tested),
        "coefficient_value_count": len(toolkit_coefficients),
        "coefficient_max_absolute_difference": coefficient_max_absolute,
        "coefficient_max_relative_difference": coefficient_max_relative,
        "logfc_max_absolute_difference": logfc_max_absolute,
        "logfc_max_relative_difference": logfc_max_relative,
    }
    return assertions, metrics


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    rscript = rscript_from_runtime(arguments.rscript)
    r_library = arguments.r_library.resolve(strict=True)
    prepare_script = PROJECT_ROOT / "scripts/benchmark/prepare_airway_fixture.R"
    oracle_script = PROJECT_ROOT / "scripts/benchmark/run_airway_direct_oracle.R"
    runner_script = Path(__file__).resolve()
    common_script = runner_script.with_name("common.py")
    schema_path = runner_script.with_name("benchmark-report-v1.schema.json")
    adapter_script = PROJECT_ROOT / "rnaseq_downstream/edger_backend.py"
    contract_script = PROJECT_ROOT / "rnaseq_downstream/analysis_contract.py"
    toolkit_script = PROJECT_ROOT / "rnaseq_downstream/r_scripts/edger_ql.R"
    airway_description = r_library / "airway/DESCRIPTION"
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
        file_evidence(prepare_script),
        file_evidence(oracle_script),
        file_evidence(adapter_script),
        file_evidence(contract_script),
        file_evidence(toolkit_script),
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
        "inputs": [file_evidence(airway_description, name="airway/DESCRIPTION")],
        "artifacts": [],
        "thresholds": {
            "relative_tolerance": RELATIVE_TOLERANCE,
            "absolute_tolerance": ABSOLUTE_TOLERANCE,
        },
        "metrics": {"scope": "backend_kernel"},
        "assertions": [
            {
                "name": "benchmark_execution",
                "passed": False,
                "observed": "not_completed",
                "requirement": "fixture, direct oracle, and toolkit route all complete",
            }
        ],
    }
    try:
        with _workspace(arguments.workspace) as workspace_value:
            workspace = Path(workspace_value)
            fixture = workspace / "fixture"
            oracle = workspace / "direct-oracle"
            toolkit = workspace / "toolkit"
            run_rscript(
                rscript,
                prepare_script,
                [fixture],
                r_library=r_library,
            )
            run_rscript(
                rscript,
                oracle_script,
                [fixture / "counts.tsv", fixture / "metadata.tsv", oracle],
                r_library=r_library,
            )
            _private_kernel(
                fixture / "counts.tsv",
                fixture / "metadata.tsv",
                toolkit,
                rscript=rscript,
                r_library=r_library,
            )
            report["inputs"] = [
                file_evidence(fixture / "counts.tsv"),
                file_evidence(fixture / "metadata.tsv"),
                file_evidence(fixture / "fixture.json"),
            ]
            direct_diagnostics = read_json(oracle / "diagnostics.json")
            toolkit_analysis = read_json(toolkit / "analysis.json")
            fixture_metadata = read_json(fixture / "fixture.json")
            expected_fixture = {
                "fixture_id": "airway-1.32.0-cell-plus-dex",
                "assay": "airway::airway",
                "airway_version": "1.32.0",
                "gene_count": 63677,
                "sample_count": 8,
                "design": "~ cell + dex",
                "factor_levels": {
                    "cell": ["N052611", "N061011", "N080611", "N61311"],
                    "dex": ["untrt", "trt"],
                },
                "coefficient": "dextrt",
            }
            expected_runtime = {
                "R": "4.6.1",
                "Bioconductor": "3.23",
                "BiocVersion_package": "3.23.1",
                "edgeR": "4.10.0",
                "tximport": "1.40.0",
                "limma": "3.68.0",
            }
            report["runtime"].update(
                {
                    "r": direct_diagnostics["r_version"],
                    "edgeR": direct_diagnostics["edgeR_version"],
                    "airway": fixture_metadata["airway_version"],
                    "toolkit_backend": toolkit_analysis.get("backend"),
                    "toolkit_runtime_identity": toolkit_analysis.get(
                        "runtime_identity"
                    ),
                }
            )
            assertions, metrics = _comparison(oracle, toolkit)
            report["assertions"] = [
                {
                    "name": "benchmark_execution",
                    "passed": True,
                    "observed": "completed",
                    "requirement": "fixture, direct oracle, and toolkit route all complete",
                },
                {
                    "name": "locked_fixture_identity",
                    "passed": fixture_metadata == expected_fixture,
                    "observed": fixture_metadata,
                    "requirement": "exact airway 1.32.0 integer assay and frozen design identity",
                },
                {
                    "name": "locked_backend_identity",
                    "passed": (
                        direct_diagnostics.get("r_version") == "4.6.1"
                        and direct_diagnostics.get("edgeR_version") == "4.10.0"
                        and toolkit_analysis.get("backend") == "edgeR_QL"
                        and toolkit_analysis.get("execution_scope")
                        == "backend_kernel_only"
                        and toolkit_analysis.get("runtime_identity") == expected_runtime
                    ),
                    "observed": {
                        "r": direct_diagnostics.get("r_version"),
                        "edgeR": direct_diagnostics.get("edgeR_version"),
                        "backend": toolkit_analysis.get("backend"),
                        "execution_scope": toolkit_analysis.get("execution_scope"),
                        "runtime_identity": toolkit_analysis.get("runtime_identity"),
                    },
                    "requirement": "R 4.6.1, edgeR 4.10.0, and private backend-kernel scope",
                },
                *assertions,
            ]
            report["artifacts"] = [
                file_evidence(oracle / "results.tsv", name="direct-oracle/results.tsv"),
                file_evidence(
                    oracle / "coefficients.tsv",
                    name="direct-oracle/coefficients.tsv",
                ),
                file_evidence(toolkit / "results.tsv", name="toolkit/results.tsv"),
                file_evidence(
                    toolkit / "coefficients.tsv", name="toolkit/coefficients.tsv"
                ),
                file_evidence(toolkit / "design.tsv", name="toolkit/design.tsv"),
            ]
            report["metrics"] = metrics
            report["status"] = (
                "pass"
                if all(item["passed"] for item in report["assertions"])
                else "fail"
            )
    except Exception as exc:
        report["assertions"] = [
            {
                "name": "benchmark_execution",
                "passed": False,
                "observed": {"error_type": type(exc).__name__, "message": str(exc)},
                "requirement": "fixture, direct oracle, and toolkit route all complete",
            }
        ]
    assert_report_shape(report, benchmark_id=BENCHMARK_ID)
    write_json(arguments.report, report)
    return report


def main() -> int:
    arguments = _arguments()
    report = run(arguments)
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
