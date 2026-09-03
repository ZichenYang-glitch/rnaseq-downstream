#!/usr/bin/env python3
"""Run the locked airway edgeR/limma pathway same-engine oracle gate."""

from __future__ import annotations

import argparse
from contextlib import nullcontext
import csv
import math
from pathlib import Path
import platform
import statistics
import tempfile
from typing import Any, Iterable, Iterator, Mapping, Sequence

from common import (
    BenchmarkError,
    PROJECT_ROOT,
    SCHEMA_VERSION,
    assert_report_shape,
    file_evidence,
    finite_float,
    read_json,
    read_tsv_records,
    rscript_from_runtime,
    run_rscript,
    sha256_file,
    write_json,
)

BENCHMARK_ID = "airway-limma-gene-set-same-engine-v1"
RELATIVE_TOLERANCE = 1e-6
ABSOLUTE_TOLERANCE = 1e-10
BH_ABSOLUTE_TOLERANCE = 1e-12
MINIMUM_TESTED_GENES = 10
CONTRAST_ID = "trt_vs_untrt"
METHOD_HYPOTHESES = {
    "limma_camera": ("directional",),
    "limma_fry": ("directional", "mixed"),
    "limma_mroast": ("directional", "mixed"),
}
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


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument("--rscript")
    parser.add_argument("--workspace", type=Path)
    return parser.parse_args()


def _workspace(path: Path | None) -> Iterator[Path]:
    if path is None:
        return tempfile.TemporaryDirectory(prefix="rnaseq-airway-pathway-")  # type: ignore[return-value]
    path.mkdir(parents=True, exist_ok=False)
    return nullcontext(path)  # type: ignore[return-value]


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
            lineterminator="\n",
            escapechar=None,
        )
        writer.writerow(header)
        writer.writerows(rows)


def _write_airway_gene_set_sources(fixture: Path) -> dict[str, Any]:
    """Build deterministic, local-only synthetic symbols and oracle gene sets."""

    counts_header, counts = read_tsv_records(fixture / "counts.tsv")
    metadata_header, metadata = read_tsv_records(fixture / "metadata.tsv")
    if counts_header[0] != "gene_id" or metadata_header != ["sample_id", "cell", "dex"]:
        raise BenchmarkError("The airway source fixture has an incompatible schema.")
    sample_ids = counts_header[1:]
    if [row["sample_id"] for row in metadata] != sample_ids:
        raise BenchmarkError("The airway source fixture sample order is not aligned.")
    control_indices = [i for i, row in enumerate(metadata) if row["dex"] == "untrt"]
    treated_indices = [i for i, row in enumerate(metadata) if row["dex"] == "trt"]
    if len(control_indices) != 4 or len(treated_indices) != 4:
        raise BenchmarkError("The airway source fixture is not the frozen 4v4 design.")

    symbols: dict[str, str] = {}
    candidates: list[tuple[float, str]] = []
    for ordinal, row in enumerate(counts, start=1):
        gene_id = row["gene_id"]
        symbols[gene_id] = f"AWSYM{ordinal:06d}"
        try:
            values = [int(row[sample]) for sample in sample_ids]
        except ValueError as exc:
            raise BenchmarkError(f"Non-integer airway count for {gene_id}.") from exc
        if min(values) < 10 or sum(values) < 1000:
            continue
        control_mean = statistics.fmean(values[index] for index in control_indices)
        treated_mean = statistics.fmean(values[index] for index in treated_indices)
        crude_logfc = math.log2((treated_mean + 0.5) / (control_mean + 0.5))
        candidates.append((crude_logfc, gene_id))
    if len(candidates) < 200:
        raise BenchmarkError("Too few robustly expressed airway genes for oracle sets.")
    candidates.sort(key=lambda item: (item[0], item[1]))
    down = [gene_id for _, gene_id in candidates[:50]]
    up = [gene_id for _, gene_id in candidates[-50:]]
    ordinary = [gene_id for _, gene_id in sorted(candidates, key=lambda x: x[1])[:50]]
    mixed = down[:25] + up[-25:]
    definitions: dict[str, tuple[str, list[str]]] = {
        "AIRWAY_DOWN": ("deterministic_crude_down_fixture_set", down),
        "AIRWAY_MIXED": ("deterministic_balanced_fixture_set", mixed),
        "AIRWAY_ORDINARY": ("deterministic_identifier_order_fixture_set", ordinary),
        "AIRWAY_UP": ("deterministic_crude_up_fixture_set", up),
        "BELOW_MINIMUM": ("mapped_but_below_minimum_negative_control", ordinary[:3]),
        "UNMAPPED_ONLY": (
            "unmapped_symbols_negative_control",
            [f"ABSENT_SYMBOL_{index:02d}" for index in range(1, 13)],
        ),
    }
    annotation_path = fixture / "annotation.tsv"
    gmt_path = fixture / "sets.gmt"
    _write_tsv(
        annotation_path,
        ["gene_id", "symbol"],
        ((gene_id, symbol) for gene_id, symbol in symbols.items()),
    )
    with gmt_path.open("w", encoding="utf-8", newline="") as handle:
        for identifier in sorted(definitions):
            description, members = definitions[identifier]
            translated = [symbols.get(gene_id, gene_id) for gene_id in members]
            if len(translated) != len(set(translated)):
                raise BenchmarkError(f"Duplicate member in airway set {identifier}.")
            handle.write("\t".join([identifier, description, *translated]) + "\n")
    source_manifest = {
        "fixture_id": "airway-1.32.0-c2-pathway-synthetic-symbols",
        "construction": (
            "synthetic one-to-one symbols; candidate genes require every count >=10 "
            "and total count >=1000; ordinary sets have 50 members; negative controls "
            "exercise below-minimum and unmapped status paths"
        ),
        "minimum_tested_genes": MINIMUM_TESTED_GENES,
        "sets": {
            identifier: {
                "description": description,
                "raw_member_count": len(members),
            }
            for identifier, (description, members) in sorted(definitions.items())
        },
    }
    manifest_path = fixture / "gene-set-source.json"
    write_json(manifest_path, source_manifest)
    return {
        "request": {
            "gmt": {
                "path": str(gmt_path.resolve()),
                "sha256": sha256_file(gmt_path),
                "collection": "airway_c2_oracle_fixture",
                "version": "1.0.0",
                "identifier_type": "symbol",
            },
            "annotation": {
                "path": str(annotation_path.resolve()),
                "sha256": sha256_file(annotation_path),
                "name": "airway_c2_synthetic_symbols",
                "version": "1.0.0",
                "gene_id_column": "gene_id",
                "symbol_column": "symbol",
            },
            "minimum_tested_genes": MINIMUM_TESTED_GENES,
        },
        "definitions": definitions,
        "manifest": source_manifest,
    }


def _private_kernel(
    counts_path: Path,
    metadata_path: Path,
    output_dir: Path,
    gene_sets: Mapping[str, Any],
    *,
    rscript: Path,
    r_library: Path,
) -> Any:
    try:
        from rnaseq_downstream.edger_backend import _run_edger_ql_benchmark_kernel
    except (ImportError, AttributeError) as exc:
        raise BenchmarkError("The private pathway benchmark adapter is unavailable.") from exc
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
                "contrast_id": CONTRAST_ID,
                "weights": {"dextrt": 1.0},
                "lfc_threshold": 0.0,
            }
        ],
        gene_sets=gene_sets,
        rscript=str(rscript),
        r_library=r_library,
    )


def _strict_pathway_rows(path: Path) -> list[dict[str, str]]:
    header, rows = read_tsv_records(path)
    if header != PATHWAY_HEADER:
        raise BenchmarkError(
            f"Pathway TSV header differs from the certified schema: {header!r}"
        )
    expected_order = sorted(
        rows,
        key=lambda row: (
            row["contrast_id"],
            row["gene_set_id"],
            row["method_id"],
            row["hypothesis"],
        ),
    )
    if rows != expected_order:
        raise BenchmarkError("Pathway result rows are not in canonical lexical order.")
    keys = [
        (row["contrast_id"], row["gene_set_id"], row["method_id"], row["hypothesis"])
        for row in rows
    ]
    if len(keys) != len(set(keys)):
        raise BenchmarkError("Pathway result rows contain duplicate inference keys.")
    return rows


def _bh(p_values: Sequence[float]) -> list[float]:
    count = len(p_values)
    adjusted = [0.0] * count
    running = 1.0
    ordered = sorted(range(count), key=lambda index: (p_values[index], index))
    for position in range(count - 1, -1, -1):
        index = ordered[position]
        running = min(running, min(1.0, p_values[index] * count / (position + 1)))
        adjusted[index] = running
    return adjusted


def _bh_parity(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    families: dict[tuple[str, str, str], list[Mapping[str, str]]] = {}
    for row in rows:
        if row["status"] == "tested":
            key = (row["contrast_id"], row["method_id"], row["hypothesis"])
            families.setdefault(key, []).append(row)
    violations = 0
    maximum = 0.0
    family_ids_valid = True
    for family_rows in families.values():
        if len({row["fdr_family_id"] for row in family_rows}) != 1:
            family_ids_valid = False
        p_values = [
            finite_float(row["p_value"], field="p_value", gene_id=row["gene_set_id"])
            for row in family_rows
        ]
        observed = [
            finite_float(row["fdr"], field="fdr", gene_id=row["gene_set_id"])
            for row in family_rows
        ]
        for expected, actual in zip(_bh(p_values), observed, strict=True):
            difference = abs(expected - actual)
            maximum = max(maximum, difference)
            violations += int(difference > BH_ABSOLUTE_TOLERANCE)
    return {
        "family_count": len(families),
        "violations": violations,
        "maximum_absolute_difference": maximum,
        "family_ids_valid": family_ids_valid,
    }


def _numeric_comparison(
    oracle: Mapping[tuple[str, str, str, str], Mapping[str, str]],
    toolkit: Mapping[tuple[str, str, str, str], Mapping[str, str]],
    field: str,
    *,
    predicate: Any | None = None,
) -> dict[str, object]:
    violations = 0
    maximum_absolute = 0.0
    maximum_relative = 0.0
    value_count = 0
    for key, expected_row in oracle.items():
        if predicate is not None and not predicate(expected_row):
            continue
        expected = finite_float(expected_row[field], field=field, gene_id=key[1])
        observed = finite_float(toolkit[key][field], field=field, gene_id=key[1])
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


def _comparison(
    direct_path: Path,
    toolkit_path: Path,
    definitions: Mapping[str, tuple[str, list[str]]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    direct_header, direct_rows = read_tsv_records(direct_path)
    expected_direct_header = [
        "contrast_id", "gene_set_id", "method_id", "hypothesis", "direction",
        "proportion_down", "proportion_up", "p_value", "fdr", "method_ngenes",
        "correlation_estimate_raw", "correlation_effective", "vif_used",
    ]
    if direct_header != expected_direct_header:
        raise BenchmarkError("Independent pathway oracle has an incompatible schema.")
    toolkit_rows = _strict_pathway_rows(toolkit_path)
    all_set_ids = sorted(definitions)
    expected_grid = [
        (CONTRAST_ID, set_id, method, hypothesis)
        for set_id in all_set_ids
        for method in sorted(METHOD_HYPOTHESES)
        for hypothesis in METHOD_HYPOTHESES[method]
    ]
    observed_grid = [
        (row["contrast_id"], row["gene_set_id"], row["method_id"], row["hypothesis"])
        for row in toolkit_rows
    ]
    expected_not_tested = {"BELOW_MINIMUM", "UNMAPPED_ONLY"}
    status_valid = all(
        row["status"] == ("not_tested" if row["gene_set_id"] in expected_not_tested else "tested")
        and row["status_reason"]
        == ("tested_gene_count_below_minimum" if row["gene_set_id"] in expected_not_tested else "")
        for row in toolkit_rows
    )
    expected_method_metadata = {
        "limma_mroast": (
            "self_contained", "corroborative", "not_applicable",
            "performed_fixed_seed_9999_rotations",
        ),
        "limma_fry": (
            "self_contained", "primary", "not_applicable",
            "not_applicable_analytic_approximation",
        ),
        "limma_camera": (
            "competitive", "supplementary", "estimated_set_specific",
            "not_applicable",
        ),
    }
    method_metadata_valid = True
    mapping_metadata_valid = True
    universe_counts: set[str] = set()
    for row in toolkit_rows:
        method = row["method_id"]
        test_class, role, tested_correlation, tested_rotation = expected_method_metadata[method]
        expected_correlation = (
            tested_correlation
            if row["status"] == "tested"
            else ("not_estimated_not_tested" if method == "limma_camera" else "not_applicable")
        )
        expected_rotation = (
            tested_rotation
            if row["status"] == "tested"
            else (
                "not_performed_not_tested"
                if method == "limma_mroast"
                else tested_rotation
            )
        )
        method_metadata_valid &= (
            row["test_class"] == test_class
            and row["inference_role"] == role
            and row["fdr_family_id"]
            == f"{CONTRAST_ID}|{method}|{row['hypothesis']}"
            and row["correlation_status"] == expected_correlation
            and row["rotation_status"] == expected_rotation
        )
        set_id = row["gene_set_id"]
        raw_count = len(definitions[set_id][1])
        if set_id == "UNMAPPED_ONLY":
            mapped_symbols = mapped_genes = tested_genes = 0
            unmapped_symbols = raw_count
            mapping_rate = 0.0
        else:
            mapped_symbols = mapped_genes = tested_genes = raw_count
            unmapped_symbols = 0
            mapping_rate = 1.0
        mapping_metadata_valid &= (
            row["gene_set_description"] == definitions[set_id][0]
            and int(row["gmt_member_count_raw"]) == raw_count
            and int(row["gmt_symbol_count_unique"]) == raw_count
            and int(row["mapped_symbol_count_unique"]) == mapped_symbols
            and int(row["ambiguous_symbol_count_unique"]) == 0
            and int(row["unmapped_symbol_count_unique"]) == unmapped_symbols
            and float(row["mapping_rate"]) == mapping_rate
            and int(row["mapped_gene_id_count_unique"]) == mapped_genes
            and int(row["tested_gene_count"]) == tested_genes
            and int(row["filtered_gene_count"]) == 0
        )
        universe_counts.add(row["tested_universe_gene_count"])
    mapping_metadata_valid &= len(universe_counts) == 1
    tested_toolkit_rows = [row for row in toolkit_rows if row["status"] == "tested"]
    toolkit = {
        (row["contrast_id"], row["gene_set_id"], row["method_id"], row["hypothesis"]): row
        for row in tested_toolkit_rows
    }
    direct = {
        (row["contrast_id"], row["gene_set_id"], row["method_id"], row["hypothesis"]): row
        for row in direct_rows
    }
    tested_grid_equal = set(direct) == set(toolkit)
    categorical_equal = tested_grid_equal and all(
        direct[key]["direction"] == toolkit[key]["direction"]
        and direct[key]["method_ngenes"] == toolkit[key]["method_ngenes"]
        for key in direct
    )
    numeric: dict[str, dict[str, object]] = {
        "p_value": _numeric_comparison(direct, toolkit, "p_value"),
        "fdr": _numeric_comparison(direct, toolkit, "fdr"),
        "proportion_down": _numeric_comparison(
            direct, toolkit, "proportion_down",
            predicate=lambda row: row["method_id"] == "limma_mroast",
        ),
        "proportion_up": _numeric_comparison(
            direct, toolkit, "proportion_up",
            predicate=lambda row: row["method_id"] == "limma_mroast",
        ),
        "correlation_estimate_raw": _numeric_comparison(
            direct, toolkit, "correlation_estimate_raw",
            predicate=lambda row: row["method_id"] == "limma_camera",
        ),
        "correlation_effective": _numeric_comparison(
            direct, toolkit, "correlation_effective",
            predicate=lambda row: row["method_id"] == "limma_camera",
        ),
        "vif_used": _numeric_comparison(
            direct, toolkit, "vif_used",
            predicate=lambda row: row["method_id"] == "limma_camera",
        ),
    } if tested_grid_equal else {}
    bh = _bh_parity(toolkit_rows)
    assertions = [
        {
            "name": "exact_pathway_grid_status_and_order",
            "passed": observed_grid == expected_grid and status_valid,
            "observed": {
                "expected_rows": len(expected_grid),
                "observed_rows": len(observed_grid),
                "grid_equal": observed_grid == expected_grid,
                "status_semantics_valid": status_valid,
            },
            "requirement": "exact canonical 6-set/5-hypothesis grid with explicit negative-control statuses",
        },
        {
            "name": "tested_method_shape_and_direction_parity",
            "passed": categorical_equal,
            "observed": {
                "direct_rows": len(direct),
                "toolkit_tested_rows": len(toolkit),
                "key_sets_equal": tested_grid_equal,
            },
            "requirement": "identical tested keys, method NGenes, and direction labels",
        },
        {
            "name": "mapping_and_method_metadata_parity",
            "passed": bool(mapping_metadata_valid and method_metadata_valid),
            "observed": {
                "mapping_metadata_valid": bool(mapping_metadata_valid),
                "method_metadata_valid": bool(method_metadata_valid),
                "tested_universe_counts": sorted(universe_counts),
            },
            "requirement": "mapping counts/rates and self-contained/competitive method labels match frozen semantics",
        },
        *[
            {
                "name": f"{field}_numeric_parity",
                "passed": comparison["violations"] == 0,
                "observed": comparison,
                "requirement": f"all independent-oracle {field} values satisfy rtol=1e-6 and atol=1e-10",
            }
            for field, comparison in numeric.items()
        ],
        {
            "name": "python_bh_parity",
            "passed": bh["violations"] == 0 and bh["family_ids_valid"],
            "observed": bh,
            "requirement": "independent Python BH per contrast/method/hypothesis agrees within 1e-12",
        },
    ]
    return assertions, {
        "scope": "backend_kernel_pathway",
        "gene_set_count": len(all_set_ids),
        "pathway_row_count": len(toolkit_rows),
        "tested_pathway_row_count": len(tested_toolkit_rows),
        "not_tested_pathway_row_count": len(toolkit_rows) - len(tested_toolkit_rows),
        "numeric_parity": numeric,
        "bh_parity": bh,
    }


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    rscript = rscript_from_runtime(arguments.rscript)
    r_library = arguments.r_library.resolve(strict=True)
    runner = Path(__file__).resolve()
    prepare = PROJECT_ROOT / "scripts/benchmark/prepare_airway_fixture.R"
    direct_script = runner.with_name("run_airway_pathway_direct_oracle.R")
    implementation_paths = [
        runner,
        runner.with_name("common.py"),
        runner.with_name("benchmark-report-v1.schema.json"),
        prepare,
        direct_script,
        PROJECT_ROOT / "rnaseq_downstream/edger_backend.py",
        PROJECT_ROOT / "rnaseq_downstream/analysis_contract.py",
        PROJECT_ROOT / "rnaseq_downstream/gene_sets.py",
        PROJECT_ROOT / "rnaseq_downstream/r_scripts/edger_ql.R",
        PROJECT_ROOT / "rnaseq_downstream/r_scripts/pathway_tests.R",
        PROJECT_ROOT / "conda-lock.yml",
        PROJECT_ROOT / "renv.lock",
        PROJECT_ROOT / "environment.p0.yml",
        PROJECT_ROOT / "environment/r-sources.lock",
        PROJECT_ROOT / "environment/verify.R",
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
        "implementation": [file_evidence(path) for path in implementation_paths],
        "inputs": [file_evidence(r_library / "airway/DESCRIPTION", name="airway/DESCRIPTION")],
        "artifacts": [],
        "thresholds": {
            "relative_tolerance": RELATIVE_TOLERANCE,
            "absolute_tolerance": ABSOLUTE_TOLERANCE,
            "bh_absolute_tolerance": BH_ABSOLUTE_TOLERANCE,
            "minimum_tested_genes": MINIMUM_TESTED_GENES,
            "mroast_nrot": 9999,
            "mroast_seed": 1729,
        },
        "metrics": {"scope": "backend_kernel_pathway"},
        "assertions": [{
            "name": "benchmark_execution", "passed": False,
            "observed": "not_completed",
            "requirement": "fixture, frozen sources, independent oracle, and toolkit route all complete",
        }],
    }
    try:
        with _workspace(arguments.workspace) as workspace_value:
            workspace = Path(workspace_value)
            fixture = workspace / "fixture"
            direct = workspace / "direct-oracle"
            toolkit = workspace / "toolkit"
            run_rscript(rscript, prepare, [fixture], r_library=r_library)
            sources = _write_airway_gene_set_sources(fixture)
            request = sources["request"]
            run_rscript(
                rscript,
                direct_script,
                [
                    fixture / "counts.tsv", fixture / "metadata.tsv",
                    request["gmt"]["path"], request["annotation"]["path"],
                    MINIMUM_TESTED_GENES, direct,
                ],
                r_library=r_library,
                timeout_seconds=900,
            )
            _private_kernel(
                fixture / "counts.tsv", fixture / "metadata.tsv", toolkit, request,
                rscript=rscript, r_library=r_library,
            )
            diagnostics = read_json(direct / "diagnostics.json")
            analysis = read_json(toolkit / "analysis.json")
            assertions, metrics = _comparison(
                direct / "pathway_oracle.tsv",
                toolkit / "pathway_results.tsv",
                sources["definitions"],
            )
            expected_runtime = {
                "R": "4.6.1", "Bioconductor": "3.23",
                "BiocVersion_package": "3.23.1", "edgeR": "4.10.0",
                "tximport": "1.40.0", "limma": "3.68.0",
            }
            identity_passed = (
                diagnostics.get("r_version") == "4.6.1"
                and diagnostics.get("edgeR_version") == "4.10.0"
                and diagnostics.get("limma_version") == "3.68.0"
                and analysis.get("runtime_identity") == expected_runtime
                and analysis.get("execution_scope") == "backend_kernel_only"
            )
            report["inputs"] = [
                file_evidence(
                    r_library / "airway/DESCRIPTION", name="airway/DESCRIPTION"
                ),
                *(
                    file_evidence(fixture / name, name=f"fixture/{name}")
                    for name in (
                        "counts.tsv", "metadata.tsv", "fixture.json", "sets.gmt",
                        "annotation.tsv", "gene-set-source.json",
                    )
                ),
            ]
            report["artifacts"] = [
                file_evidence(direct / "pathway_oracle.tsv", name="direct-oracle/pathway_oracle.tsv"),
                file_evidence(direct / "diagnostics.json", name="direct-oracle/diagnostics.json"),
                file_evidence(toolkit / "pathway_results.tsv", name="toolkit/pathway_results.tsv"),
            ]
            report["runtime"].update({
                "r": diagnostics.get("r_version"),
                "edgeR": diagnostics.get("edgeR_version"),
                "limma": diagnostics.get("limma_version"),
                "airway": "1.32.0",
                "toolkit_runtime_identity": analysis.get("runtime_identity"),
            })
            report["metrics"] = metrics
            report["assertions"] = [
                {
                    "name": "benchmark_execution", "passed": True,
                    "observed": "completed",
                    "requirement": "fixture, frozen sources, independent oracle, and toolkit route all complete",
                },
                {
                    "name": "locked_runtime_identity", "passed": identity_passed,
                    "observed": {
                        "direct_r": diagnostics.get("r_version"),
                        "direct_edgeR": diagnostics.get("edgeR_version"),
                        "direct_limma": diagnostics.get("limma_version"),
                        "toolkit": analysis.get("runtime_identity"),
                    },
                    "requirement": "R 4.6.1, edgeR 4.10.0, limma 3.68.0 in private backend scope",
                },
                *assertions,
            ]
            report["status"] = "pass" if all(item["passed"] for item in report["assertions"]) else "fail"
    except Exception as exc:
        report["assertions"] = [{
            "name": "benchmark_execution", "passed": False,
            "observed": {"error_type": type(exc).__name__, "message": str(exc)},
            "requirement": "fixture, frozen sources, independent oracle, and toolkit route all complete",
        }]
    assert_report_shape(report, benchmark_id=BENCHMARK_ID)
    write_json(arguments.report, report)
    return report


def main() -> int:
    report = run(_arguments())
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
