#!/usr/bin/env python3
"""Run the locked compcodeR C2 pathway FDR/TPR regression gate."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from functools import partial
from pathlib import Path
import platform
import statistics
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
    read_tsv_records,
    rscript_from_runtime,
    run_rscript,
    sha256_file,
    write_json,
)

BENCHMARK_ID = "compcoder-limma-self-contained-fdr-tpr-v1"
MIXED_REPLICATES = tuple(range(1, 21))
COMPLETE_NULL_REPLICATES = tuple(range(1, 41))
METHODS = ("limma_mroast", "limma_fry")
NOMINAL_FDR = 0.05
MAX_MEAN_FDP = 0.10
MAX_WORST_FDP = 0.25
MIN_MEAN_TPR = 0.80
MIN_WORST_TPR = 0.60
MAX_COMPLETE_NULL_FAMILY_REJECTIONS = 4
BH_ABSOLUTE_TOLERANCE = 1e-12
MINIMUM_TESTED_GENES = 10
CONTRAST_ID = "treated_vs_control"


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument("--rscript")
    parser.add_argument("--workspace", type=Path)
    return parser.parse_args()


def _workspace(path: Path | None) -> Iterator[Path]:
    if path is None:
        return tempfile.TemporaryDirectory(prefix="rnaseq-pathway-compcoder-")  # type: ignore[return-value]
    path.mkdir(parents=True, exist_ok=False)
    return nullcontext(path)  # type: ignore[return-value]


def _gene_sets_request(fixture: Path, scenario: str) -> dict[str, Any]:
    gmt = fixture / "sets.gmt"
    annotation = fixture / "annotation.tsv"
    return {
        "gmt": {
            "path": str(gmt.resolve()),
            "sha256": sha256_file(gmt),
            "collection": f"compcodeR_c2_{scenario}_fixture",
            "version": "1.0.0",
            "identifier_type": "symbol",
        },
        "annotation": {
            "path": str(annotation.resolve()),
            "sha256": sha256_file(annotation),
            "name": "compcodeR_c2_synthetic_symbols",
            "version": "1.0.0",
            "gene_id_column": "gene_id",
            "symbol_column": "symbol",
        },
        "minimum_tested_genes": MINIMUM_TESTED_GENES,
    }


def _private_kernel(
    fixture: Path,
    output: Path,
    scenario: str,
    *,
    rscript: Path,
    r_library: Path,
) -> Any:
    try:
        from rnaseq_downstream.edger_backend import _run_edger_ql_benchmark_kernel
    except (ImportError, AttributeError) as exc:
        raise BenchmarkError("The private pathway benchmark adapter is unavailable.") from exc
    return _run_edger_ql_benchmark_kernel(
        fixture / "counts.tsv",
        fixture / "metadata.tsv",
        output,
        design={
            "intercept": True,
            "terms": ["condition"],
            "variables": {
                "condition": {"type": "factor", "levels": ["control", "treated"]}
            },
        },
        contrasts=[{
            "contrast_id": CONTRAST_ID,
            "weights": {"conditiontreated": 1.0},
            "lfc_threshold": 0.0,
        }],
        gene_sets=_gene_sets_request(fixture, scenario),
        rscript=str(rscript),
        r_library=r_library,
    )


def _expected_set_ids(scenario: str) -> list[str]:
    result = [
        f"NULL_{size:03d}_{index:02d}"
        for size in (20, 40, 80, 160)
        for index in range(1, 26)
    ]
    if scenario == "mixed":
        result.extend(f"POSITIVE_UP_{index:02d}" for index in range(1, 6))
        result.extend(f"POSITIVE_DOWN_{index:02d}" for index in range(1, 6))
    return sorted(result)


def _assert_fixture(scenario: str, replicate: int, fixture: Path) -> dict[str, Any]:
    metadata = read_json(fixture / "fixture.json")
    simulation_seed = (5000 if scenario == "mixed" else 15000) + replicate
    membership_seed = (25000 if scenario == "mixed" else 45000) + replicate
    expected_scalar = {
        "scenario": scenario,
        "compcodeR_version": "1.48.0",
        "simulation_seed": simulation_seed,
        "membership_seed": membership_seed,
        "replicate": replicate,
        "n_genes": 5000,
        "samples_per_condition": 6,
        "n_differential": 500 if scenario == "mixed" else 0,
        "n_upregulated": 250 if scenario == "mixed" else 0,
        "n_downregulated": 250 if scenario == "mixed" else 0,
        "effect_size": 2,
        "filter_threshold_total": 0,
        "null_set_count": 100,
        "positive_set_count": 10 if scenario == "mixed" else 0,
        "positive_set_size": 40 if scenario == "mixed" else 0,
        "minimum_tested_genes": 10,
    }
    observed_scalar = {key: metadata.get(key) for key in expected_scalar}
    if observed_scalar != expected_scalar:
        raise BenchmarkError(
            f"{scenario} replicate {replicate} fixture identity differs: {observed_scalar!r}"
        )
    if metadata.get("rng_kind") != {
        "generator": "Mersenne-Twister", "normal": "Inversion", "sample": "Rejection"
    }:
        raise BenchmarkError(f"{scenario} replicate {replicate} RNG kind differs.")
    expected_policy = (
        "membership_rng_is_separate_from_simulation_rng;"
        "null_sets_sample_truth_null_genes_without_replacement_within_set;"
        "positive_sets_are_disjoint_within_direction_and_sample_truth_DE_genes;"
        "no_observed_count_statistic_pvalue_or_gene_rank_is_used"
    )
    if metadata.get("construction_policy") != expected_policy:
        raise BenchmarkError(f"{scenario} replicate {replicate} policy differs.")

    counts = read_tsv(fixture / "counts.tsv")
    truth = read_tsv(fixture / "truth.tsv")
    annotation = read_tsv(fixture / "annotation.tsv")
    sample_header, sample_rows = read_tsv_records(fixture / "metadata.tsv")
    if (
        set(counts) != set(truth)
        or set(annotation) != set(truth)
        or len(counts) != 5000
        or sample_header != ["sample_id", "condition"]
        or len(sample_rows) != 12
        or [row["condition"] for row in sample_rows].count("control") != 6
        or [row["condition"] for row in sample_rows].count("treated") != 6
    ):
        raise BenchmarkError(f"{scenario} replicate {replicate} dimensions differ.")
    if any(
        annotation[gene_id]["symbol"] != truth_row["symbol"]
        for gene_id, truth_row in truth.items()
    ):
        raise BenchmarkError("Annotation is not the declared one-to-one truth mapping.")
    membership_header, membership_rows = read_tsv_records(fixture / "membership.tsv")
    if membership_header != [
        "gene_set_id", "truth_class", "true_direction", "nominal_size", "symbol"
    ]:
        raise BenchmarkError("Membership TSV schema differs from the frozen contract.")
    truth_by_symbol = {row["symbol"]: row for row in truth.values()}
    sets: dict[str, list[dict[str, str]]] = {}
    observed_pairs: set[tuple[str, str]] = set()
    for row in membership_rows:
        pair = (row["gene_set_id"], row["symbol"])
        if pair in observed_pairs:
            raise BenchmarkError(f"Duplicate within-set membership {pair!r}.")
        observed_pairs.add(pair)
        sets.setdefault(row["gene_set_id"], []).append(row)
    if sorted(sets) != _expected_set_ids(scenario):
        raise BenchmarkError(f"{scenario} replicate {replicate} set inventory differs.")
    gmt_lines = (fixture / "sets.gmt").read_text(encoding="utf-8").splitlines()
    gmt_fields = [line.split("\t") for line in gmt_lines]
    if (
        any(len(fields) < 3 for fields in gmt_fields)
        or [fields[0] for fields in gmt_fields] != _expected_set_ids(scenario)
    ):
        raise BenchmarkError(f"{scenario} replicate {replicate} GMT order differs.")
    for fields in gmt_fields:
        set_id, description, *members = fields
        expected_members = [row["symbol"] for row in sets[set_id]]
        if description != f"frozen_{scenario}_{set_id}" or members != expected_members:
            raise BenchmarkError(f"GMT and membership inventory differ for {set_id}.")
    positive_by_direction: dict[str, set[str]] = {"Up": set(), "Down": set()}
    for set_id, rows in sets.items():
        nominal_sizes = {int(row["nominal_size"]) for row in rows}
        truth_classes = {row["truth_class"] for row in rows}
        directions = {row["true_direction"] for row in rows}
        if len(nominal_sizes) != 1 or nominal_sizes != {len(rows)}:
            raise BenchmarkError(f"Nominal membership size differs for {set_id}.")
        if len(truth_classes) != 1 or len(directions) != 1:
            raise BenchmarkError(f"Truth labels vary within {set_id}.")
        truth_class = next(iter(truth_classes))
        direction = next(iter(directions))
        for row in rows:
            truth_row = truth_by_symbol.get(row["symbol"])
            if truth_row is None:
                raise BenchmarkError(f"Membership symbol is absent from truth: {row['symbol']}")
            is_de = truth_row["differential_expression"] == "1"
            lfc = finite_float(
                truth_row["true_log2_fold_change"],
                field="true_log2_fold_change", gene_id=row["symbol"],
            )
            if truth_class == "null" and (is_de or direction != "None"):
                raise BenchmarkError(f"Null membership is contaminated in {set_id}.")
            if truth_class == "positive":
                expected_sign = 1 if direction == "Up" else -1
                if not is_de or lfc * expected_sign <= 0:
                    raise BenchmarkError(f"Positive membership has wrong truth in {set_id}.")
                if row["symbol"] in positive_by_direction[direction]:
                    raise BenchmarkError(
                        f"Positive sets are not disjoint within {direction}: {row['symbol']}"
                    )
                positive_by_direction[direction].add(row["symbol"])
    return metadata


def _read_directional_rows(output: Path, scenario: str) -> list[dict[str, str]]:
    header, all_rows = read_tsv_records(output / "pathway_results.tsv")
    required = {
        "contrast_id", "gene_set_id", "method_id", "hypothesis", "status",
        "direction", "p_value", "fdr", "fdr_family_id", "tested_gene_count",
        "method_ngenes",
    }
    if not required <= set(header):
        raise BenchmarkError("Pathway output is missing simulation-gate fields.")
    expected_sets = _expected_set_ids(scenario)
    observed_sets = sorted({row["gene_set_id"] for row in all_rows})
    if observed_sets != expected_sets:
        raise BenchmarkError(f"{scenario} output set inventory differs.")
    method_hypotheses = {
        "limma_camera": ("directional",),
        "limma_fry": ("directional", "mixed"),
        "limma_mroast": ("directional", "mixed"),
    }
    expected_grid = [
        (CONTRAST_ID, set_id, method, hypothesis)
        for set_id in expected_sets
        for method in sorted(method_hypotheses)
        for hypothesis in method_hypotheses[method]
    ]
    observed_grid = [
        (row["contrast_id"], row["gene_set_id"], row["method_id"], row["hypothesis"])
        for row in all_rows
    ]
    if observed_grid != expected_grid:
        raise BenchmarkError(
            f"{scenario} output does not have the exact canonical five-cell set grid."
        )
    for row in all_rows:
        if row["contrast_id"] != CONTRAST_ID or row["status"] != "tested":
            raise BenchmarkError("A simulation set was filtered/not tested or used another contrast.")
        tested = int(row["tested_gene_count"])
        method_n = int(row["method_ngenes"])
        if tested < MINIMUM_TESTED_GENES or method_n != tested:
            raise BenchmarkError("A simulation set failed the no-replacement tested-size rule.")
    result = [
        row for row in all_rows
        if row["method_id"] in METHODS and row["hypothesis"] == "directional"
    ]
    if len(result) != len(expected_sets) * len(METHODS):
        raise BenchmarkError("Directional self-contained result grid is incomplete.")
    return result


def _bh(p_values: Sequence[float]) -> list[float]:
    count = len(p_values)
    result = [0.0] * count
    running = 1.0
    order = sorted(range(count), key=lambda index: (p_values[index], index))
    for position in range(count - 1, -1, -1):
        index = order[position]
        running = min(running, min(1.0, p_values[index] * count / (position + 1)))
        result[index] = running
    return result


def _bh_parity(output: Path) -> dict[str, object]:
    _, rows = read_tsv_records(output / "pathway_results.tsv")
    families: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        if row["status"] == "tested":
            key = (row["contrast_id"], row["method_id"], row["hypothesis"])
            families.setdefault(key, []).append(row)
    violations = 0
    maximum = 0.0
    identifiers_valid = True
    for family_rows in families.values():
        identifiers_valid &= len({row["fdr_family_id"] for row in family_rows}) == 1
        p_values = [
            finite_float(row["p_value"], field="p_value", gene_id=row["gene_set_id"])
            for row in family_rows
        ]
        fdr_values = [
            finite_float(row["fdr"], field="fdr", gene_id=row["gene_set_id"])
            for row in family_rows
        ]
        for expected, observed in zip(_bh(p_values), fdr_values, strict=True):
            difference = abs(expected - observed)
            maximum = max(maximum, difference)
            violations += int(difference > BH_ABSOLUTE_TOLERANCE)
    return {
        "family_count": len(families),
        "violations": violations,
        "maximum_absolute_difference": maximum,
        "family_ids_valid": bool(identifiers_valid),
    }


def _score_mixed(replicate: int, output: Path) -> dict[str, object]:
    rows = _read_directional_rows(output, "mixed")
    method_metrics: dict[str, dict[str, object]] = {}
    for method in METHODS:
        method_rows = [row for row in rows if row["method_id"] == method]
        false_positives = 0
        true_positives = 0
        wrong_direction = 0
        for row in method_rows:
            called = finite_float(
                row["fdr"], field="fdr", gene_id=row["gene_set_id"]
            ) <= NOMINAL_FDR
            if row["gene_set_id"].startswith("NULL_"):
                false_positives += int(called)
            else:
                expected_direction = "Up" if "_UP_" in row["gene_set_id"] else "Down"
                true_positives += int(called and row["direction"] == expected_direction)
                wrong_direction += int(called and row["direction"] != expected_direction)
        discoveries = false_positives + true_positives + wrong_direction
        method_metrics[method] = {
            "null_set_count": 100,
            "positive_set_count": 10,
            "discoveries": discoveries,
            "false_positives": false_positives,
            "true_positives": true_positives,
            "wrong_direction_significant_positives": wrong_direction,
            "false_discovery_proportion": false_positives / discoveries if discoveries else 0.0,
            "true_positive_rate": true_positives / 10,
        }
    return {
        "scenario": "mixed",
        "replicate": replicate,
        "simulation_seed": 5000 + replicate,
        "membership_seed": 25000 + replicate,
        "mroast_seed": 1729,
        "methods": method_metrics,
        "bh_parity": _bh_parity(output),
    }


def _score_complete_null(replicate: int, output: Path) -> dict[str, object]:
    rows = _read_directional_rows(output, "complete_null")
    family_rejections = {
        method: int(any(
            finite_float(row["fdr"], field="fdr", gene_id=row["gene_set_id"])
            <= NOMINAL_FDR
            for row in rows if row["method_id"] == method
        ))
        for method in METHODS
    }
    return {
        "scenario": "complete_null",
        "replicate": replicate,
        "simulation_seed": 15000 + replicate,
        "membership_seed": 45000 + replicate,
        "mroast_seed": 1729,
        "family_rejections": family_rejections,
        "bh_parity": _bh_parity(output),
    }


def _run_replicate(
    job: tuple[str, int],
    *,
    workspace: Path,
    generator: Path,
    rscript: Path,
    r_library: Path,
) -> dict[str, Any]:
    """Generate and score one independent replicate in a private directory."""

    scenario, replicate = job
    label = f"{scenario}-replicate-{replicate:02d}"
    fixture = workspace / label / "fixture"
    output = workspace / label / "toolkit"
    run_rscript(
        rscript,
        generator,
        [scenario, replicate, fixture],
        r_library=r_library,
    )
    _assert_fixture(scenario, replicate, fixture)
    _private_kernel(
        fixture,
        output,
        scenario,
        rscript=rscript,
        r_library=r_library,
    )
    analysis = read_json(output / "analysis.json")
    if (
        analysis.get("execution_scope") != "backend_kernel_only"
        or analysis.get("input_semantics") != "benchmark_kernel_integer_counts"
    ):
        raise BenchmarkError(f"{label} did not use private backend scope.")
    return {
        "scenario": scenario,
        "inputs": [
            file_evidence(fixture / name, name=f"{label}/{name}")
            for name in (
                "counts.tsv", "metadata.tsv", "truth.tsv", "fixture.json",
                "sets.gmt", "annotation.tsv", "membership.tsv",
            )
        ],
        "artifact": file_evidence(
            output / "pathway_results.tsv",
            name=f"{label}/pathway_results.tsv",
        ),
        "metrics": (
            _score_mixed(replicate, output)
            if scenario == "mixed"
            else _score_complete_null(replicate, output)
        ),
        "runtime_identity": analysis.get("runtime_identity"),
    }


def _aggregate(
    mixed: Sequence[Mapping[str, Any]],
    complete_null: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    per_method: dict[str, dict[str, object]] = {}
    assertions: list[dict[str, object]] = [
        {
            "name": "exact_replicate_and_set_counts",
            "passed": len(mixed) == 20 and len(complete_null) == 40,
            "observed": {
                "mixed_replicates": len(mixed), "mixed_sets_each": 110,
                "complete_null_replicates": len(complete_null),
                "complete_null_sets_each": 100,
            },
            "requirement": "exactly 20 mixed 110-set and 40 complete-null 100-set replicates",
        }
    ]
    for method in METHODS:
        fdps = [float(item["methods"][method]["false_discovery_proportion"]) for item in mixed]
        tprs = [float(item["methods"][method]["true_positive_rate"]) for item in mixed]
        wrong = sum(int(item["methods"][method]["wrong_direction_significant_positives"]) for item in mixed)
        complete_rejections = sum(int(item["family_rejections"][method]) for item in complete_null)
        method_metrics = {
            "mean_mixed_fdp": statistics.fmean(fdps),
            "worst_mixed_fdp": max(fdps),
            "mean_mixed_tpr": statistics.fmean(tprs),
            "worst_mixed_tpr": min(tprs),
            "wrong_direction_significant_positives": wrong,
            "complete_null_family_rejections": complete_rejections,
            "complete_null_empirical_family_fdr": complete_rejections / 40,
        }
        per_method[method] = method_metrics
        checks = (
            ("mean_mixed_fdp", method_metrics["mean_mixed_fdp"] <= MAX_MEAN_FDP, "mean FDP <= 0.10"),
            ("worst_mixed_fdp", method_metrics["worst_mixed_fdp"] <= MAX_WORST_FDP, "worst FDP <= 0.25"),
            ("mean_mixed_tpr", method_metrics["mean_mixed_tpr"] >= MIN_MEAN_TPR, "mean TPR >= 0.80"),
            ("worst_mixed_tpr", method_metrics["worst_mixed_tpr"] >= MIN_WORST_TPR, "worst TPR >= 0.60"),
            ("wrong_direction", wrong == 0, "zero wrong-direction significant positive sets"),
            ("complete_null_family_fdr", complete_rejections <= MAX_COMPLETE_NULL_FAMILY_REJECTIONS, "at most 4 of 40 complete-null families reject"),
        )
        for suffix, passed, requirement in checks:
            assertions.append({
                "name": f"{method}_{suffix}_gate",
                "passed": bool(passed),
                "observed": method_metrics.get(suffix, wrong if suffix == "wrong_direction" else complete_rejections),
                "requirement": requirement,
            })
    bh_records = [item["bh_parity"] for item in [*mixed, *complete_null]]
    bh_violations = sum(int(item["violations"]) for item in bh_records)
    bh_maximum = max(float(item["maximum_absolute_difference"]) for item in bh_records)
    bh_ids = all(bool(item["family_ids_valid"]) for item in bh_records)
    assertions.append({
        "name": "python_bh_parity",
        "passed": bh_violations == 0 and bh_ids,
        "observed": {
            "replicates": len(bh_records), "violations": bh_violations,
            "maximum_absolute_difference": bh_maximum,
            "family_ids_valid": bh_ids,
        },
        "requirement": "Python BH per contrast/method/hypothesis agrees within 1e-12",
    })
    return assertions, {
        "scope": "backend_kernel_pathway",
        "mixed_replicate_count": len(mixed),
        "complete_null_replicate_count": len(complete_null),
        "method_metrics": per_method,
        "mixed_replicates": list(mixed),
        "complete_null_replicates": list(complete_null),
        "gross_regression_limitation": (
            "Frozen finite Monte Carlo scenarios detect gross calibration/power regressions; "
            "they do not establish universal pathway-test error control or biological power."
        ),
        "membership_construction": (
            "A separate frozen RNG samples truth-null genes without replacement within "
            "each null set; positive sets are disjoint within direction and sample pure "
            "truth-DE genes without using counts, fitted statistics, p-values, or ranks."
        ),
        "null_family_metric": (
            "Per complete-null replicate and method, the family outcome is one if any "
            "of 100 BH-adjusted directional tests rejects, otherwise zero."
        ),
    }


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    rscript = rscript_from_runtime(arguments.rscript)
    r_library = arguments.r_library.resolve(strict=True)
    runner = Path(__file__).resolve()
    generator = runner.with_name("generate_pathway_compcoder_fixture.R")
    implementation_paths = [
        runner, runner.with_name("common.py"),
        runner.with_name("benchmark-report-v1.schema.json"), generator,
        PROJECT_ROOT / "rnaseq_downstream/edger_backend.py",
        PROJECT_ROOT / "rnaseq_downstream/analysis_contract.py",
        PROJECT_ROOT / "rnaseq_downstream/gene_sets.py",
        PROJECT_ROOT / "rnaseq_downstream/r_scripts/edger_ql.R",
        PROJECT_ROOT / "rnaseq_downstream/r_scripts/pathway_tests.R",
        PROJECT_ROOT / "conda-lock.yml", PROJECT_ROOT / "renv.lock",
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
        "inputs": [file_evidence(r_library / "compcodeR/DESCRIPTION", name="compcodeR/DESCRIPTION")],
        "artifacts": [],
        "thresholds": {
            "nominal_bh_fdr": NOMINAL_FDR,
            "maximum_mean_mixed_fdp": MAX_MEAN_FDP,
            "maximum_worst_mixed_fdp": MAX_WORST_FDP,
            "minimum_mean_mixed_tpr": MIN_MEAN_TPR,
            "minimum_worst_mixed_tpr": MIN_WORST_TPR,
            "maximum_complete_null_family_rejections": MAX_COMPLETE_NULL_FAMILY_REJECTIONS,
            "complete_null_replicates": 40,
            "bh_absolute_tolerance": BH_ABSOLUTE_TOLERANCE,
            "selection_context": "frozen_after_disclosed_exploratory_design_pilot",
            "mroast_production_seed": 1729,
            "mroast_seed_policy": (
                "The production seed is reset per contrast and is intentionally independent "
                "of simulation and membership RNG streams; sharing it across replicates can "
                "induce Monte Carlo dependence, so complete-null evidence is a gross-regression gate."
            ),
        },
        "metrics": {"scope": "backend_kernel_pathway"},
        "assertions": [{
            "name": "benchmark_execution", "passed": False,
            "observed": "not_completed",
            "requirement": "all 20 mixed and 40 complete-null pathway simulations complete",
        }],
    }
    try:
        with _workspace(arguments.workspace) as workspace_value:
            workspace = Path(workspace_value)
            inputs: list[dict[str, object]] = []
            artifacts: list[dict[str, object]] = []
            mixed_metrics: list[dict[str, object]] = []
            null_metrics: list[dict[str, object]] = []
            runtime_identities: list[object] = []
            jobs = [
                *(('mixed', replicate) for replicate in MIXED_REPLICATES),
                *(('complete_null', replicate) for replicate in COMPLETE_NULL_REPLICATES),
            ]
            worker = partial(
                _run_replicate,
                workspace=workspace,
                generator=generator,
                rscript=rscript,
                r_library=r_library,
            )
            with ThreadPoolExecutor(max_workers=4) as executor:
                # executor.map preserves the frozen job order even if workers
                # finish out of order. Any worker exception propagates here and
                # is serialized by the outer fail-report boundary.
                replicate_results = list(executor.map(worker, jobs))
            for result in replicate_results:
                inputs.extend(result["inputs"])
                artifacts.append(result["artifact"])
                runtime_identities.append(result["runtime_identity"])
                if result["scenario"] == "mixed":
                    mixed_metrics.append(result["metrics"])
                else:
                    null_metrics.append(result["metrics"])

            # Regenerate representative fixtures without rerunning inference and
            # demand byte-identical count, membership, and source streams.
            determinism_checks: dict[str, bool] = {}
            for scenario in ("mixed", "complete_null"):
                original = workspace / f"{scenario}-replicate-01" / "fixture"
                rerun = workspace / f"determinism-{scenario}" / "fixture"
                run_rscript(rscript, generator, [scenario, 1, rerun], r_library=r_library)
                for name in (
                    "counts.tsv", "metadata.tsv", "truth.tsv", "fixture.json",
                    "sets.gmt", "annotation.tsv", "membership.tsv",
                ):
                    determinism_checks[f"{scenario}/{name}"] = (
                        sha256_file(original / name) == sha256_file(rerun / name)
                    )

            aggregate_assertions, metrics = _aggregate(mixed_metrics, null_metrics)
            expected_runtime = {
                "R": "4.6.1", "Bioconductor": "3.23",
                "BiocVersion_package": "3.23.1", "edgeR": "4.10.0",
                "tximport": "1.40.0", "limma": "3.68.0",
            }
            runtime_passed = all(identity == expected_runtime for identity in runtime_identities)
            assertions = [
                {
                    "name": "benchmark_execution", "passed": True,
                    "observed": {"mixed": len(mixed_metrics), "complete_null": len(null_metrics)},
                    "requirement": "all 20 mixed and 40 complete-null pathway simulations complete",
                },
                {
                    "name": "locked_runtime_identity", "passed": runtime_passed,
                    "observed": {"replicate_identities_equal": runtime_passed, "expected": expected_runtime},
                    "requirement": "all replicates use R 4.6.1, edgeR 4.10.0, and limma 3.68.0",
                },
                {
                    "name": "fixture_byte_determinism", "passed": all(determinism_checks.values()),
                    "observed": determinism_checks,
                    "requirement": "representative simulation, truth, membership, GMT, and annotation bytes rerun identically",
                },
                *aggregate_assertions,
            ]
            report["inputs"] = [
                file_evidence(
                    r_library / "compcodeR/DESCRIPTION",
                    name="compcodeR/DESCRIPTION",
                ),
                *inputs,
            ]
            report["artifacts"] = artifacts
            report["runtime"].update({
                "r": "4.6.1", "edgeR": "4.10.0", "limma": "3.68.0",
                "compcodeR": "1.48.0", "replicate_workers": 4,
            })
            report["metrics"] = metrics
            report["assertions"] = assertions
            report["status"] = "pass" if all(item["passed"] for item in assertions) else "fail"
    except Exception as exc:
        report["assertions"] = [{
            "name": "benchmark_execution", "passed": False,
            "observed": {"error_type": type(exc).__name__, "message": str(exc)},
            "requirement": "all 20 mixed and 40 complete-null pathway simulations complete",
        }]
    assert_report_shape(report, benchmark_id=BENCHMARK_ID)
    write_json(arguments.report, report)
    return report


def main() -> int:
    report = run(_arguments())
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
