#!/usr/bin/env python3
"""Run the locked DESeq2 airway Wald/LRT same-engine oracle gate."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
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

BENCHMARK_ID = "airway-deseq2-wald-lrt-same-engine-v1"
RELATIVE_TOLERANCE = 1e-6
ABSOLUTE_TOLERANCE = 1e-10
CONTRAST_ID = "trt_vs_untrt"
DESIGN = {
    "intercept": True,
    "terms": ["cell", "dex"],
    "variables": {
        "cell": {
            "type": "factor",
            "levels": ["N052611", "N061011", "N080611", "N61311"],
        },
        "dex": {"type": "factor", "levels": ["untrt", "trt"]},
    },
}
REDUCED_DESIGN = {
    "intercept": True,
    "terms": ["cell"],
}


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument("--rscript")
    parser.add_argument(
        "--workspace",
        type=Path,
        help="Retain fixture, direct-oracle, and public-toolkit outputs here.",
    )
    return parser.parse_args()


def _workspace(path: Path | None) -> Iterator[Path]:
    if path is None:
        return tempfile.TemporaryDirectory(  # type: ignore[return-value]
            prefix="rnaseq-deseq2-airway-oracle-"
        )
    path.mkdir(parents=True, exist_ok=False)
    return nullcontext(path)  # type: ignore[return-value]


def _strict_cli_document(stdout: str, *, command: str) -> dict[str, Any]:
    lines = stdout.splitlines()
    if len(lines) != 1 or not lines[0]:
        raise BenchmarkError(
            f"Public {command} violated the one-document stdout contract"
        )

    def reject_constant(value: str) -> None:
        raise ValueError(f"non-finite JSON constant {value}")

    def unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        document = json.loads(
            lines[0],
            parse_constant=reject_constant,
            object_pairs_hook=unique_pairs,
        )
    except (json.JSONDecodeError, ValueError) as exc:
        raise BenchmarkError(f"Public {command} returned invalid JSON: {exc}") from exc
    if not isinstance(document, dict):
        raise BenchmarkError(f"Public {command} response is not an object")
    return document


def _run_cli(
    command: str,
    arguments: list[object],
    *,
    timeout_seconds: int = 1200,
) -> dict[str, Any]:
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "rnaseq_downstream",
            command,
            *(str(value) for value in arguments),
        ],
        cwd=PROJECT_ROOT,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    document = _strict_cli_document(completed.stdout, command=command)
    if completed.returncode != 0 or document.get("status") != "success":
        raise BenchmarkError(
            f"Public {command} failed with exit {completed.returncode}: {document}"
        )
    if document.get("command") != command:
        raise BenchmarkError(f"Public {command} returned the wrong command identity")
    return document


def _write_public_input_fixture(fixture: Path) -> Path:
    """Wrap the serialized airway assay in the documented self-attested route."""

    counts_path = fixture / "counts.tsv"
    metadata_path = fixture / "metadata.tsv"
    header, _ = read_tsv_records(counts_path)
    if not header or header[0] != "gene_id":
        raise BenchmarkError("The airway fixture count header is invalid")
    samples = header[1:]
    metadata_header, metadata_rows = read_tsv_records(metadata_path)
    if metadata_header != ["sample_id", "cell", "dex"]:
        raise BenchmarkError("The airway fixture metadata header is invalid")
    if [row["sample_id"] for row in metadata_rows] != samples:
        raise BenchmarkError("The airway fixture sample order is invalid")

    reference = fixture / "benchmark-reference.fa"
    reference.write_text(
        ">airway_benchmark_route_fixture\nACGT\n",
        encoding="utf-8",
    )
    producer = {"name": "featureCounts", "version": "2.0.6"}
    reference_identity = {
        "name": "airway-benchmark-route-fixture",
        "version": "1",
        "source": reference.name,
    }
    manifest = {
        "schema_version": "1.0",
        "artifact_type": "featurecounts_integer_matrix",
        "artifact": {
            "path": counts_path.name,
            "sha256": sha256_file(counts_path),
        },
        "gene_id_column": "gene_id",
        "display_columns": [],
        "sample_columns": samples,
        "producer": producer,
        "reference": {
            **reference_identity,
            "sha256": sha256_file(reference),
        },
    }
    write_json(fixture / "counts.manifest.json", manifest)
    request = {
        "schema_version": "1.0",
        "input_semantics": "featurecounts_integer",
        "metadata": {
            "path": metadata_path.name,
            "sample_id_column": "sample_id",
        },
        "producer": producer,
        "reference": reference_identity,
        "gene_id": {"strip_version": False},
        "samples": [{"sample_id": sample} for sample in samples],
        "featurecounts": {
            "layout": "combined_matrix",
            "matrix": counts_path.name,
            "manifest": "counts.manifest.json",
        },
    }
    request_path = fixture / "input.request.json"
    write_json(request_path, request)
    return request_path


def _write_analysis_request(workspace: Path, *, mode: str) -> Path:
    deseq2: dict[str, object] = {"test_mode": mode, "shrinkage": "none"}
    if mode == "lrt":
        deseq2["reduced"] = REDUCED_DESIGN
    request = {
        "schema_version": "1.2",
        "backend": "deseq2",
        "validated_input_bundle": "validated-input",
        "design": DESIGN,
        "contrasts": [
            {
                "contrast_id": CONTRAST_ID,
                "weights": {"dextrt": 1},
                "lfc_threshold": 0,
            }
        ],
        "deseq2": deseq2,
    }
    request_path = workspace / f"{mode}.analysis-request.json"
    write_json(request_path, request)
    return request_path


def _run_public_toolkit_route(
    workspace: Path,
    *,
    rscript: Path,
    r_library: Path,
) -> dict[str, dict[str, Any]]:
    fixture = workspace / "fixture"
    input_request = _write_public_input_fixture(fixture)
    validated = workspace / "validated-input"
    receipts: dict[str, dict[str, Any]] = {}
    receipts["inspect"] = _run_cli("inspect", ["--request", input_request])
    receipts["validate"] = _run_cli(
        "validate",
        ["--request", input_request, "--output-dir", validated],
    )
    for mode in ("wald", "lrt"):
        request = _write_analysis_request(workspace, mode=mode)
        output = workspace / f"toolkit-{mode}"
        receipts[f"run_{mode}"] = _run_cli(
            "run",
            [
                "--request",
                request,
                "--output-dir",
                output,
                "--rscript",
                rscript,
                "--r-library",
                r_library,
            ],
        )
        receipts[f"summarize_{mode}"] = _run_cli(
            "summarize", ["--run-dir", output]
        )
    return receipts


def _numeric_parity(
    oracle: dict[str, dict[str, str]],
    toolkit: dict[str, dict[str, str]],
    *,
    oracle_field: str,
    toolkit_field: str,
) -> dict[str, int | float]:
    violations = 0
    maximum_absolute = 0.0
    maximum_relative = 0.0
    value_count = 0
    if set(oracle) == set(toolkit):
        for gene_id, expected_row in oracle.items():
            expected = finite_float(
                expected_row[oracle_field], field=oracle_field, gene_id=gene_id
            )
            observed = finite_float(
                toolkit[gene_id][toolkit_field],
                field=toolkit_field,
                gene_id=gene_id,
            )
            difference = abs(observed - expected)
            relative = difference / max(abs(expected), 1e-300)
            maximum_absolute = max(maximum_absolute, difference)
            maximum_relative = max(maximum_relative, relative)
            value_count += 1
            if not math.isclose(
                observed,
                expected,
                rel_tol=RELATIVE_TOLERANCE,
                abs_tol=ABSOLUTE_TOLERANCE,
            ):
                violations += 1
    return {
        "value_count": value_count,
        "violations": violations,
        "max_absolute_difference": maximum_absolute,
        "max_relative_difference": maximum_relative,
    }


def _coefficient_parity(
    oracle_directory: Path,
    toolkit_directory: Path,
) -> tuple[bool, dict[str, int | float]]:
    oracle_header, oracle_rows = read_tsv_records(
        oracle_directory / "coefficients.tsv"
    )
    toolkit_header, toolkit_rows = read_tsv_records(
        toolkit_directory / "coefficients.tsv"
    )
    if oracle_header != ["gene_id", "coefficient", "estimate"]:
        raise BenchmarkError("Direct coefficient output has an invalid schema")
    if toolkit_header != [
        "gene_id",
        "status",
        "status_reason",
        "coefficient",
        "estimate",
        "scale",
    ]:
        raise BenchmarkError("Toolkit coefficient output has an invalid schema")

    direct: dict[tuple[str, str], str] = {}
    for row in oracle_rows:
        key = (row["gene_id"], row["coefficient"])
        if not all(key) or key in direct:
            raise BenchmarkError("Direct coefficients contain a duplicate/blank key")
        direct[key] = row["estimate"]
    routed: dict[tuple[str, str], str] = {}
    for row in toolkit_rows:
        if row["status"] != "tested":
            continue
        if row["scale"] != "log2":
            raise BenchmarkError("Toolkit DESeq2 coefficient scale is not log2")
        key = (row["gene_id"], row["coefficient"])
        if not all(key) or key in routed:
            raise BenchmarkError("Toolkit coefficients contain a duplicate/blank key")
        routed[key] = row["estimate"]

    same_keys = set(direct) == set(routed)
    violations = 0
    maximum_absolute = 0.0
    maximum_relative = 0.0
    if same_keys:
        for key, expected_value in direct.items():
            gene_id, coefficient = key
            expected = finite_float(
                expected_value,
                field=f"coefficient:{coefficient}",
                gene_id=gene_id,
            )
            observed = finite_float(
                routed[key],
                field=f"coefficient:{coefficient}",
                gene_id=gene_id,
            )
            difference = abs(observed - expected)
            relative = difference / max(abs(expected), 1e-300)
            maximum_absolute = max(maximum_absolute, difference)
            maximum_relative = max(maximum_relative, relative)
            if not math.isclose(
                observed,
                expected,
                rel_tol=RELATIVE_TOLERANCE,
                abs_tol=ABSOLUTE_TOLERANCE,
            ):
                violations += 1
    metrics: dict[str, int | float] = {
        "oracle_value_count": len(direct),
        "toolkit_value_count": len(routed),
        "violations": violations,
        "max_absolute_difference": maximum_absolute,
        "max_relative_difference": maximum_relative,
    }
    return same_keys and violations == 0, metrics


def _mode_comparison(
    oracle_directory: Path,
    toolkit_directory: Path,
    *,
    mode: str,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    oracle = read_tsv(oracle_directory / "results.tsv")
    toolkit_rows = read_tsv(toolkit_directory / "results.tsv")
    toolkit_tested: dict[str, dict[str, str]] = {}
    unexpected_contrasts: set[str] = set()
    for gene_id, row in toolkit_rows.items():
        if row.get("contrast_id") != CONTRAST_ID:
            unexpected_contrasts.add(row.get("contrast_id", ""))
        if row.get("status") == "tested":
            toolkit_tested[gene_id] = row
    gene_sets_equal = set(oracle) == set(toolkit_tested) and not unexpected_contrasts

    coefficient_passed, coefficient_metrics = _coefficient_parity(
        oracle_directory, toolkit_directory
    )
    fields = {
        "logfc": ("logFC", "logFC"),
        "statistic": ("statistic", "statistic"),
        "pvalue": ("PValue", "PValue"),
        "fdr": ("FDR", "FDR"),
    }
    numeric = {
        name: _numeric_parity(
            oracle,
            toolkit_tested,
            oracle_field=oracle_field,
            toolkit_field=toolkit_field,
        )
        for name, (oracle_field, toolkit_field) in fields.items()
    }

    analysis = read_json(toolkit_directory / "analysis.json")
    expected_semantics = {
        "statistic_type": "LRT" if mode == "lrt" else "Wald",
        "statistic_hypothesis": (
            "full_vs_reduced_omnibus"
            if mode == "lrt"
            else "contrast_equals_zero"
        ),
        "fdr_basis": (
            "omnibus_pvalue_BH"
            if mode == "lrt"
            else "contrast_pvalue_BH_after_independent_filtering"
        ),
        "test_method": (
            "DESeq2::results_LRT" if mode == "lrt" else "DESeq2::results_Wald"
        ),
    }
    result_semantics = all(
        all(row[field] == expected for row in toolkit_rows.values())
        for field, expected in expected_semantics.items()
    ) and all(row["lfc_threshold"] == "0" for row in toolkit_rows.values())
    reporting_effect = analysis.get("reporting_effect")
    analysis_semantics = (
        analysis.get("test", {}).get("mode") == mode
        and isinstance(reporting_effect, list)
        and len(reporting_effect) == 1
        and reporting_effect[0].get("contrast_id") == CONTRAST_ID
        and reporting_effect[0].get("weights") == {"dextrt": 1}
        and reporting_effect[0].get("role")
        == (
            "reported_effect_not_lrt_hypothesis"
            if mode == "lrt"
            else "tested_contrast"
        )
    )
    semantic_passed = result_semantics and analysis_semantics

    prefix = mode
    assertions: list[dict[str, object]] = [
        {
            "name": f"{prefix}_tested_gene_set_parity",
            "passed": gene_sets_equal,
            "observed": {
                "oracle_tested": len(oracle),
                "toolkit_tested": len(toolkit_tested),
                "unexpected_contrasts": sorted(unexpected_contrasts),
            },
            "requirement": "identical tested genes and one reporting contrast",
        },
        {
            "name": f"{prefix}_coefficient_numeric_parity",
            "passed": coefficient_passed,
            "observed": coefficient_metrics,
            "requirement": "all fitted coefficients satisfy rtol=1e-6 and atol=1e-10",
        },
        {
            "name": f"{prefix}_result_semantics",
            "passed": semantic_passed,
            "observed": {
                "result_semantics": result_semantics,
                "analysis_semantics": analysis_semantics,
                "expected": expected_semantics,
            },
            "requirement": (
                "Wald rows encode contrast inference; LRT rows encode omnibus "
                "inference and a reporting-only effect"
            ),
        },
    ]
    for name, comparison in numeric.items():
        assertions.append(
            {
                "name": f"{prefix}_{name}_numeric_parity",
                "passed": gene_sets_equal and comparison["violations"] == 0,
                "observed": comparison,
                "requirement": (
                    f"all tested {name} values satisfy rtol=1e-6 and atol=1e-10"
                ),
            }
        )
    metrics: dict[str, object] = {
        "input_gene_count": len(toolkit_rows),
        "tested_gene_count": len(toolkit_tested),
        "non_tested_gene_count": len(toolkit_rows) - len(toolkit_tested),
        "coefficient_parity": coefficient_metrics,
        **numeric,
    }
    return assertions, metrics


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    rscript = rscript_from_runtime(arguments.rscript)
    r_library = arguments.r_library.resolve(strict=True)
    runner_script = Path(__file__).resolve()
    direct_script = runner_script.with_name("run_deseq2_airway_direct_oracle.R")
    prepare_script = runner_script.with_name("prepare_airway_fixture.R")
    common_script = runner_script.with_name("common.py")
    schema_path = runner_script.with_name("benchmark-report-v1.schema.json")
    implementation_paths = [
        runner_script,
        direct_script,
        prepare_script,
        common_script,
        schema_path,
        PROJECT_ROOT / "rnaseq_downstream/input_semantics.py",
        PROJECT_ROOT / "rnaseq_downstream/validation_bundle.py",
        PROJECT_ROOT / "rnaseq_downstream/analysis_contract_v12.py",
        PROJECT_ROOT / "rnaseq_downstream/deseq2_backend.py",
        PROJECT_ROOT / "rnaseq_downstream/r_scripts/deseq2.R",
        PROJECT_ROOT / "rnaseq_downstream/run_summary.py",
        PROJECT_ROOT / "rnaseq_downstream/cli.py",
        PROJECT_ROOT / "conda-lock.yml",
        PROJECT_ROOT / "renv.lock",
        PROJECT_ROOT / "environment.p0.yml",
        PROJECT_ROOT / "environment/r-sources.lock",
        PROJECT_ROOT / "environment/verify.R",
    ]
    airway_description = r_library / "airway/DESCRIPTION"
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
        "implementation": [file_evidence(path) for path in implementation_paths],
        "inputs": [file_evidence(airway_description, name="airway/DESCRIPTION")],
        "artifacts": [],
        "thresholds": {
            "relative_tolerance": RELATIVE_TOLERANCE,
            "absolute_tolerance": ABSOLUTE_TOLERANCE,
        },
        "metrics": {
            "scope": "public_validate_run_summarize_same_engine_oracle"
        },
        "assertions": [
            {
                "name": "benchmark_execution",
                "passed": False,
                "observed": "not_completed",
                "requirement": "fixture, direct oracle, and both public routes complete",
            }
        ],
    }
    try:
        with _workspace(arguments.workspace) as workspace_value:
            workspace = Path(workspace_value)
            fixture = workspace / "fixture"
            direct = workspace / "direct-oracle"
            run_rscript(rscript, prepare_script, [fixture], r_library=r_library)
            run_rscript(
                rscript,
                direct_script,
                [fixture / "counts.tsv", fixture / "metadata.tsv", direct],
                r_library=r_library,
                timeout_seconds=1200,
            )
            receipts = _run_public_toolkit_route(
                workspace, rscript=rscript, r_library=r_library
            )

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
                "DESeq2": "1.52.0",
                "apeglm": "1.34.0",
                "tximport": "1.40.0",
            }
            all_assertions: list[dict[str, object]] = [
                {
                    "name": "benchmark_execution",
                    "passed": True,
                    "observed": "completed",
                    "requirement": (
                        "fixture, independent direct oracle, validate, both runs, "
                        "and both summaries complete"
                    ),
                },
                {
                    "name": "locked_fixture_identity",
                    "passed": fixture_metadata == expected_fixture,
                    "observed": fixture_metadata,
                    "requirement": "exact airway 1.32.0 integer assay and frozen design",
                },
                {
                    "name": "featurecounts_origin_not_certified",
                    "passed": True,
                    "observed": {
                        "route_fixture": "typed_combined_integer_matrix",
                        "origin_authentication": "self_attested_not_proven",
                        "claim_scope": "numerical_parity_and_public_routing_only",
                    },
                    "requirement": (
                        "the benchmark must disclose that serialization of airway "
                        "does not prove featureCounts producer provenance"
                    ),
                },
            ]
            mode_metrics: dict[str, object] = {}
            runtime_observed: dict[str, object] = {}
            for mode in ("wald", "lrt"):
                diagnostics = read_json(direct / mode / "diagnostics.json")
                analysis = read_json(workspace / f"toolkit-{mode}/analysis.json")
                runtime_observed[mode] = {
                    "r": diagnostics.get("r_version"),
                    "BiocVersion_package": diagnostics.get(
                        "bioconductor_version"
                    ),
                    "DESeq2": diagnostics.get("deseq2_version"),
                    "apeglm": diagnostics.get("apeglm_version"),
                    "tximport": diagnostics.get("tximport_version"),
                    "toolkit_runtime": analysis.get("runtime_identity"),
                }
                runtime_passed = (
                    diagnostics.get("r_version") == "4.6.1"
                    and diagnostics.get("bioconductor_version") == "3.23.1"
                    and diagnostics.get("deseq2_version") == "1.52.0"
                    and diagnostics.get("apeglm_version") == "1.34.0"
                    and diagnostics.get("tximport_version") == "1.40.0"
                    and analysis.get("runtime_identity") == expected_runtime
                    and analysis.get("backend") == "DESeq2"
                    and analysis.get("execution_scope")
                    == "validated_p1_deseq2_input"
                )
                all_assertions.append(
                    {
                        "name": f"{mode}_locked_backend_identity",
                        "passed": runtime_passed,
                        "observed": runtime_observed[mode],
                        "requirement": (
                            "R 4.6.1, Bioconductor 3.23, DESeq2 1.52.0, "
                            "and the validated public backend scope"
                        ),
                    }
                )
                assertions, metrics = _mode_comparison(
                    direct / mode,
                    workspace / f"toolkit-{mode}",
                    mode=mode,
                )
                all_assertions.extend(assertions)
                mode_metrics[mode] = metrics

            receipt_summary = {
                name: {
                    "command": value.get("command"),
                    "status": value.get("status"),
                    "verified_status": (
                        value.get("data", {}).get("status")
                        if name.startswith("summarize_")
                        else None
                    ),
                }
                for name, value in receipts.items()
            }
            receipt_passed = (
                set(receipts)
                == {
                    "inspect",
                    "validate",
                    "run_wald",
                    "summarize_wald",
                    "run_lrt",
                    "summarize_lrt",
                }
                and all(value["status"] == "success" for value in receipt_summary.values())
                and all(
                    receipt_summary[f"summarize_{mode}"]["verified_status"]
                    == "verified_complete"
                    for mode in ("wald", "lrt")
                )
            )
            all_assertions.append(
                {
                    "name": "public_workflow_receipts",
                    "passed": receipt_passed,
                    "observed": receipt_summary,
                    "requirement": (
                        "public inspect, validate, Wald/LRT run, and summarize "
                        "commands return verified success"
                    ),
                }
            )
            report["runtime"].update(
                {
                    "r": "4.6.1",
                    "Bioconductor": "3.23",
                    "DESeq2": "1.52.0",
                    "apeglm": "1.34.0",
                    "tximport": "1.40.0",
                    "per_mode_observed": runtime_observed,
                }
            )
            report["inputs"] = [
                file_evidence(fixture / "counts.tsv"),
                file_evidence(fixture / "metadata.tsv"),
                file_evidence(fixture / "fixture.json"),
                file_evidence(fixture / "benchmark-reference.fa"),
                file_evidence(fixture / "counts.manifest.json"),
                file_evidence(fixture / "input.request.json"),
                file_evidence(workspace / "wald.analysis-request.json"),
                file_evidence(workspace / "lrt.analysis-request.json"),
            ]
            artifacts: list[dict[str, object]] = []
            for mode in ("wald", "lrt"):
                artifacts.extend(
                    [
                        file_evidence(
                            direct / mode / "results.tsv",
                            name=f"direct-oracle/{mode}/results.tsv",
                        ),
                        file_evidence(
                            direct / mode / "coefficients.tsv",
                            name=f"direct-oracle/{mode}/coefficients.tsv",
                        ),
                        file_evidence(
                            workspace / f"toolkit-{mode}/results.tsv",
                            name=f"toolkit/{mode}/results.tsv",
                        ),
                        file_evidence(
                            workspace / f"toolkit-{mode}/coefficients.tsv",
                            name=f"toolkit/{mode}/coefficients.tsv",
                        ),
                        file_evidence(
                            workspace / f"toolkit-{mode}/design.tsv",
                            name=f"toolkit/{mode}/design.tsv",
                        ),
                    ]
                )
            report["artifacts"] = artifacts
            report["metrics"] = {
                "scope": "public_validate_run_summarize_same_engine_oracle",
                "featurecounts_origin_authentication": "self_attested_not_proven",
                "certified_claim": "numerical_parity_and_public_routing_only",
                "modes": mode_metrics,
            }
            report["assertions"] = all_assertions
            report["status"] = (
                "pass" if all(item["passed"] for item in all_assertions) else "fail"
            )
    except Exception as exc:
        report["assertions"] = [
            {
                "name": "benchmark_execution",
                "passed": False,
                "observed": {"error_type": type(exc).__name__, "message": str(exc)},
                "requirement": (
                    "fixture, independent direct oracle, validate, both runs, "
                    "and both summaries complete"
                ),
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
