#!/usr/bin/env python3
"""Diagnose the disclosed DESeq2 compcodeR calibration failure read-only."""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import csv
import html
import math
from pathlib import Path
import platform
import statistics
import sys
import tempfile
from typing import Any, Iterable, Iterator, Mapping, Sequence

from common import (
    BenchmarkError,
    PROJECT_ROOT,
    file_evidence,
    read_json,
    read_tsv_records,
    rscript_from_runtime,
    run_rscript,
    sha256_file,
    write_json,
)

DIAGNOSTIC_ID = "compcoder-deseq2-calibration-mechanism-exploratory-v1"
SCHEMA_VERSION = "rnaseq-downstream-diagnostic-report-v1"
EXPLORATORY_REPORT_ID = "compcoder-deseq2-nb-exploratory-v1"
EXPLORATORY_SEEDS = tuple(range(61001, 61021))
REPLICATES = tuple(range(1, 21))
WORKERS = 2
NOMINAL_FDR = 0.05
DISPERSION_LABELS = (
    "Q1_lowest",
    "Q2_low",
    "Q3_middle",
    "Q4_high",
    "Q5_highest",
)
EXPECTED_RUNTIME = {
    "R": "4.6.1",
    "DESeq2": "1.52.0",
    "edgeR": "4.10.0",
    "compcodeR": "1.48.0",
}
FIXTURE_FILES = ("counts.tsv", "metadata.tsv", "truth.tsv", "fixture.json")
GENE_FIELDS = (
    "gene_id",
    "truth_class",
    "true_dispersion",
    "base_mean",
    "dispersion_gene_estimate",
    "dispersion_fitted_trend",
    "dispersion_final",
    "dispersion_outlier",
    "true_dispersion_above_deseq2_max_disp",
    "wald_pvalue",
    "wald_fdr_if_on",
    "wald_fdr_if_off",
    "lrt_pvalue",
    "lrt_fdr_if_on",
    "edger_tested",
    "edger_pvalue",
    "edger_fdr_native",
)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--r-library", type=Path, required=True)
    parser.add_argument("--rscript")
    parser.add_argument("--workspace", type=Path)
    return parser.parse_args()


def _workspace(path: Path | None) -> Iterator[Path]:
    if path is None:
        return tempfile.TemporaryDirectory(  # type: ignore[return-value]
            prefix="rnaseq-deseq2-mechanism-exploratory-"
        )
    path.mkdir(parents=True, exist_ok=False)
    return nullcontext(path)  # type: ignore[return-value]


def _optional_float(value: str, *, field: str, gene_id: str) -> float | None:
    if value == "":
        return None
    try:
        parsed = float(value)
    except ValueError as error:
        raise BenchmarkError(
            f"Invalid {field} for gene {gene_id!r}: {value!r}"
        ) from error
    if not math.isfinite(parsed):
        raise BenchmarkError(f"Non-finite {field} for gene {gene_id!r}.")
    return parsed


def _quantile_type7(values: Sequence[float], probability: float) -> float:
    if not values or not 0 <= probability <= 1:
        raise BenchmarkError("Cannot calculate the requested diagnostic quantile.")
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def _dispersion_strata(
    values: Sequence[float], gene_ids: Sequence[str]
) -> tuple[list[str], list[dict[str, Any]]]:
    if len(values) != len(gene_ids) or len(values) % len(DISPERSION_LABELS) != 0:
        raise BenchmarkError("Dispersion ranks cannot be split into equal quintiles.")
    ordered_indices = sorted(
        range(len(values)), key=lambda index: (values[index], gene_ids[index])
    )
    group_size = len(values) // len(DISPERSION_LABELS)
    labels = [""] * len(values)
    for rank, original_index in enumerate(ordered_indices):
        labels[original_index] = DISPERSION_LABELS[rank // group_size]
    boundaries: list[dict[str, Any]] = []
    for boundary in range(1, len(DISPERSION_LABELS)):
        lower_index = ordered_indices[boundary * group_size - 1]
        upper_index = ordered_indices[boundary * group_size]
        boundaries.append(
            {
                "lower_rank": boundary * group_size,
                "lower_max_true_dispersion": values[lower_index],
                "lower_max_gene_id": gene_ids[lower_index],
                "upper_min_true_dispersion": values[upper_index],
                "upper_min_gene_id": gene_ids[upper_index],
            }
        )
    return labels, boundaries


def _bh_adjust(pairs: Sequence[tuple[str, float]]) -> dict[str, float]:
    ordered = sorted(pairs, key=lambda item: (item[1], item[0]))
    count = len(ordered)
    adjusted: dict[str, float] = {}
    running = 1.0
    for reverse_index in range(count - 1, -1, -1):
        gene_id, pvalue = ordered[reverse_index]
        rank = reverse_index + 1
        running = min(running, pvalue * count / rank, 1.0)
        adjusted[gene_id] = running
    return adjusted


def _add_aligned_bh(rows: Sequence[dict[str, Any]]) -> int:
    for row in rows:
        row["aligned_tested"] = (
            row["wald_pvalue"] is not None and row["edger_pvalue"] is not None
        )
    common_rows = [row for row in rows if row["aligned_tested"]]
    wald_adjusted = _bh_adjust(
        [(str(row["gene_id"]), float(row["wald_pvalue"])) for row in common_rows]
    )
    edger_adjusted = _bh_adjust(
        [(str(row["gene_id"]), float(row["edger_pvalue"])) for row in common_rows]
    )
    for row in rows:
        gene_id = str(row["gene_id"])
        row["wald_fdr_aligned"] = wald_adjusted.get(gene_id)
        row["edger_fdr_aligned"] = edger_adjusted.get(gene_id)
    return len(common_rows)


def _score(
    rows: Sequence[Mapping[str, Any]],
    *,
    fdr_field: str,
    tested_field: str | None = None,
) -> dict[str, Any]:
    calls = 0
    true_positives = 0
    false_positives = 0
    tested = 0
    truth_de_total = 0
    truth_de_in_tested_family = 0
    for row in rows:
        is_de = row["truth_class"] == "DE"
        truth_de_total += int(is_de)
        fdr = row.get(fdr_field)
        is_tested = fdr is not None
        if tested_field is not None:
            is_tested = is_tested and bool(row[tested_field])
        tested += int(is_tested)
        truth_de_in_tested_family += int(is_tested and is_de)
        called = is_tested and float(fdr) <= NOMINAL_FDR
        calls += int(called)
        true_positives += int(called and is_de)
        false_positives += int(called and not is_de)
    return {
        "tested_genes": tested,
        "discoveries": calls,
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_discovery_proportion": false_positives / calls if calls else 0.0,
        "truth_de_total": truth_de_total,
        "truth_de_in_tested_family": truth_de_in_tested_family,
        "all_truth_de_tpr": (
            true_positives / truth_de_total if truth_de_total else 0.0
        ),
        "conditional_tpr_within_tested_family": (
            true_positives / truth_de_in_tested_family
            if truth_de_in_tested_family
            else 0.0
        ),
    }


def _parse_gene_table(path: Path, *, replicate: int) -> list[dict[str, Any]]:
    fields, raw_rows = read_tsv_records(path)
    if tuple(fields) != GENE_FIELDS or len(raw_rows) != 5000:
        raise BenchmarkError(
            f"Replicate {replicate} diagnostic table has the wrong shape."
        )
    rows: list[dict[str, Any]] = []
    identifiers: set[str] = set()
    numeric_fields = GENE_FIELDS[2:]
    for raw in raw_rows:
        gene_id = raw["gene_id"]
        if not gene_id or gene_id in identifiers:
            raise BenchmarkError(
                f"Replicate {replicate} has a blank or duplicate gene identifier."
            )
        identifiers.add(gene_id)
        if raw["truth_class"] not in {"null", "DE"}:
            raise BenchmarkError(
                f"Replicate {replicate} has invalid truth class for {gene_id}."
            )
        row: dict[str, Any] = {
            "replicate": replicate,
            "seed": 61000 + replicate,
            "gene_id": gene_id,
            "truth_class": raw["truth_class"],
        }
        for field in numeric_fields:
            if field in {
                "dispersion_outlier",
                "true_dispersion_above_deseq2_max_disp",
                "edger_tested",
            }:
                if raw[field] not in {"0", "1"}:
                    raise BenchmarkError(
                        f"Replicate {replicate} has invalid {field} state for {gene_id}."
                    )
                row[field] = raw[field] == "1"
            else:
                row[field] = _optional_float(raw[field], field=field, gene_id=gene_id)
        if row["true_dispersion"] is None or row["true_dispersion"] <= 0:
            raise BenchmarkError(
                f"Replicate {replicate} has invalid true dispersion for {gene_id}."
            )
        rows.append(row)
    labels, boundaries = _dispersion_strata(
        [float(row["true_dispersion"]) for row in rows],
        [str(row["gene_id"]) for row in rows],
    )
    for row, label in zip(rows, labels, strict=True):
        row["dispersion_stratum"] = label
    _add_aligned_bh(rows)
    rows[0]["_stratum_boundaries"] = boundaries
    return rows


def _validate_exploratory_report(
    path: Path,
) -> tuple[dict[str, Any], dict[str, Mapping[str, Any]]]:
    report = read_json(path)
    if (
        report.get("benchmark_id") != EXPLORATORY_REPORT_ID
        or report.get("status") != "pass"
    ):
        raise BenchmarkError("The disclosed exploratory report is not complete.")
    metrics = report.get("metrics")
    if not isinstance(metrics, Mapping) or metrics.get("seed_grid") != list(
        EXPLORATORY_SEEDS
    ):
        raise BenchmarkError("The disclosed report has the wrong exploratory grid.")
    if metrics.get("evidence_role") != "threshold_selection_only_not_certification":
        raise BenchmarkError("The disclosed report has the wrong evidence role.")
    indexed: dict[str, Mapping[str, Any]] = {}
    for record in report.get("inputs", []):
        if not isinstance(record, Mapping) or not isinstance(record.get("name"), str):
            raise BenchmarkError("The disclosed report has malformed input evidence.")
        indexed[str(record["name"])] = record
    return report, indexed


def _verify_fixture_evidence(
    fixture: Path,
    *,
    replicate: int,
    archived: Mapping[str, Mapping[str, Any]],
) -> list[dict[str, object]]:
    evidence: list[dict[str, object]] = []
    for filename in FIXTURE_FILES:
        name = f"replicate-{replicate:02d}/{filename}"
        expected = archived.get(name)
        observed = file_evidence(fixture / filename, name=name)
        if expected is None or any(
            observed[field] != expected.get(field)
            for field in ("name", "sha256", "size_bytes")
        ):
            raise BenchmarkError(
                f"Regenerated exploratory fixture differs from archived evidence: {name}."
            )
        evidence.append(observed)
    metadata = read_json(fixture / "fixture.json")
    if (
        metadata.get("mode") != "exploratory"
        or metadata.get("seed") != 61000 + replicate
        or metadata.get("replicate") != replicate
    ):
        raise BenchmarkError(
            f"Replicate {replicate} fixture metadata is not exploratory."
        )
    return evidence


def _run_replicate(
    replicate: int,
    *,
    workspace: Path,
    generator: Path,
    diagnostic: Path,
    rscript: Path,
    r_library: Path,
    archived_inputs: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    root = workspace / f"replicate-{replicate:02d}"
    fixture = root / "fixture"
    output = root / "diagnostic"
    run_rscript(
        rscript,
        generator,
        ["exploratory", replicate, fixture],
        r_library=r_library,
    )
    evidence = _verify_fixture_evidence(
        fixture, replicate=replicate, archived=archived_inputs
    )
    run_rscript(
        rscript,
        diagnostic,
        [replicate, fixture, output],
        r_library=r_library,
        timeout_seconds=600,
    )
    fit_summary = read_json(output / "fit-summary.json")
    expected_seed = 61000 + replicate
    if (
        fit_summary.get("replicate") != replicate
        or fit_summary.get("seed") != expected_seed
        or fit_summary.get("fixture_regeneration_match") is not True
        or fit_summary.get("runtime") != EXPECTED_RUNTIME
    ):
        raise BenchmarkError(f"Replicate {replicate} fit summary is inconsistent.")
    for mode, type_field, coefficient_field in (
        ("Wald", "resolved_fit_type", "fit_coefficients"),
        ("LRT", "lrt_resolved_fit_type", "lrt_fit_coefficients"),
    ):
        resolved = fit_summary.get(type_field)
        coefficients = fit_summary.get(coefficient_field)
        if resolved not in {"parametric", "local", "mean"}:
            raise BenchmarkError(
                f"Replicate {replicate} has invalid {mode} dispersion fit type."
            )
        if resolved == "parametric":
            if not isinstance(coefficients, Mapping) or set(coefficients) != {
                "asymptDisp",
                "extraPois",
            }:
                raise BenchmarkError(
                    f"Replicate {replicate} lacks {mode} parametric coefficients."
                )
            values = [float(value) for value in coefficients.values()]
            if any(not math.isfinite(value) or value <= 0 for value in values):
                raise BenchmarkError(
                    f"Replicate {replicate} has invalid {mode} fit coefficients."
                )
    rows = _parse_gene_table(output / "gene-diagnostics.tsv", replicate=replicate)
    boundaries = rows[0].pop("_stratum_boundaries")
    scores = {
        "deseq2_wald_independent_filtering_on": _score(
            rows, fdr_field="wald_fdr_if_on"
        ),
        "deseq2_wald_independent_filtering_off": _score(
            rows, fdr_field="wald_fdr_if_off"
        ),
        "deseq2_lrt_independent_filtering_on": _score(rows, fdr_field="lrt_fdr_if_on"),
        "edger_ql_native": _score(
            rows, fdr_field="edger_fdr_native", tested_field="edger_tested"
        ),
        "deseq2_wald_aligned_common_tested_family": _score(
            rows,
            fdr_field="wald_fdr_aligned",
            tested_field="aligned_tested",
        ),
        "edger_ql_aligned_common_tested_family": _score(
            rows,
            fdr_field="edger_fdr_aligned",
            tested_field="aligned_tested",
        ),
    }
    return {
        "replicate": replicate,
        "seed": expected_seed,
        "rows": rows,
        "fit_summary": fit_summary,
        "fixture_inputs": evidence,
        "stratum_boundaries": boundaries,
        "scores": scores,
    }


def _aggregate_scores(jobs: Sequence[Mapping[str, Any]], method: str) -> dict[str, Any]:
    per_replicate = []
    for job in jobs:
        score = dict(job["scores"][method])
        score.update({"replicate": job["replicate"], "seed": job["seed"]})
        per_replicate.append(score)
    discoveries = sum(item["discoveries"] for item in per_replicate)
    false_positives = sum(item["false_positives"] for item in per_replicate)
    true_positives = sum(item["true_positives"] for item in per_replicate)
    return {
        "method": method,
        "replicate_count": len(per_replicate),
        "mean_tested_genes": statistics.fmean(
            item["tested_genes"] for item in per_replicate
        ),
        "minimum_tested_genes": min(item["tested_genes"] for item in per_replicate),
        "maximum_tested_genes": max(item["tested_genes"] for item in per_replicate),
        "mean_replicate_fdp": statistics.fmean(
            item["false_discovery_proportion"] for item in per_replicate
        ),
        "minimum_replicate_fdp": min(
            item["false_discovery_proportion"] for item in per_replicate
        ),
        "maximum_replicate_fdp": max(
            item["false_discovery_proportion"] for item in per_replicate
        ),
        "pooled_fdp_descriptive_only": (
            false_positives / discoveries if discoveries else 0.0
        ),
        "mean_all_truth_de_tpr": statistics.fmean(
            item["all_truth_de_tpr"] for item in per_replicate
        ),
        "minimum_all_truth_de_tpr": min(
            item["all_truth_de_tpr"] for item in per_replicate
        ),
        "maximum_all_truth_de_tpr": max(
            item["all_truth_de_tpr"] for item in per_replicate
        ),
        "mean_conditional_tpr_within_tested_family": statistics.fmean(
            item["conditional_tpr_within_tested_family"] for item in per_replicate
        ),
        "minimum_conditional_tpr_within_tested_family": min(
            item["conditional_tpr_within_tested_family"] for item in per_replicate
        ),
        "maximum_conditional_tpr_within_tested_family": max(
            item["conditional_tpr_within_tested_family"] for item in per_replicate
        ),
        "mean_truth_de_in_tested_family": statistics.fmean(
            item["truth_de_in_tested_family"] for item in per_replicate
        ),
        "total_discoveries": discoveries,
        "total_true_positives": true_positives,
        "total_false_positives": false_positives,
        "per_replicate": per_replicate,
    }


def _ratio_summary(values: Iterable[float]) -> dict[str, Any]:
    materialized = [value for value in values if math.isfinite(value) and value > 0]
    if not materialized:
        return {
            "count": 0,
            "median": None,
            "q25": None,
            "q75": None,
            "below_one_fraction": None,
        }
    return {
        "count": len(materialized),
        "median": statistics.median(materialized),
        "q25": _quantile_type7(materialized, 0.25),
        "q75": _quantile_type7(materialized, 0.75),
        "below_one_fraction": sum(value < 1 for value in materialized)
        / len(materialized),
    }


def _dispersion_estimation(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    groups: dict[str, list[Mapping[str, Any]]] = {
        "all_genes": list(rows),
        "truth_null": [row for row in rows if row["truth_class"] == "null"],
        "truth_de": [row for row in rows if row["truth_class"] == "DE"],
        "null_false_positive": [
            row
            for row in rows
            if row["truth_class"] == "null"
            and row["wald_fdr_if_on"] is not None
            and row["wald_fdr_if_on"] <= NOMINAL_FDR
        ],
        "null_not_called": [
            row
            for row in rows
            if row["truth_class"] == "null"
            and not (
                row["wald_fdr_if_on"] is not None
                and row["wald_fdr_if_on"] <= NOMINAL_FDR
            )
        ],
    }
    summaries: dict[str, Any] = {}
    for name, group in groups.items():
        summaries[name] = {
            "final_to_true_ratio": _ratio_summary(
                float(row["dispersion_final"]) / float(row["true_dispersion"])
                for row in group
                if row["dispersion_final"] is not None
            ),
            "gene_estimate_to_true_ratio": _ratio_summary(
                float(row["dispersion_gene_estimate"]) / float(row["true_dispersion"])
                for row in group
                if row["dispersion_gene_estimate"] is not None
                and row["dispersion_gene_estimate"] > 0
            ),
            "fitted_trend_to_true_ratio": _ratio_summary(
                float(row["dispersion_fitted_trend"]) / float(row["true_dispersion"])
                for row in group
                if row["dispersion_fitted_trend"] is not None
                and row["dispersion_fitted_trend"] > 0
            ),
        }
    return summaries


def _extreme_dispersion_tail(
    rows: Sequence[Mapping[str, Any]], jobs: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    bounds = {float(job["fit_summary"]["max_dispersion_bound"]) for job in jobs}
    if len(bounds) != 1:
        raise BenchmarkError("DESeq2 max-dispersion bounds differ across replicates.")
    bound = next(iter(bounds))
    tail = [row for row in rows if row["true_dispersion_above_deseq2_max_disp"]]
    null_tail = [row for row in tail if row["truth_class"] == "null"]
    de_tail = [row for row in tail if row["truth_class"] == "DE"]
    raw_testable_null = [row for row in null_tail if row["wald_pvalue"] is not None]
    published_tested_null = [
        row for row in null_tail if row["wald_fdr_if_on"] is not None
    ]
    false_positives = [
        row for row in published_tested_null if row["wald_fdr_if_on"] <= NOMINAL_FDR
    ]
    all_false_positives = sum(
        row["truth_class"] == "null"
        and row["wald_fdr_if_on"] is not None
        and row["wald_fdr_if_on"] <= NOMINAL_FDR
        for row in rows
    )
    all_discoveries = sum(
        row["wald_fdr_if_on"] is not None and row["wald_fdr_if_on"] <= NOMINAL_FDR
        for row in rows
    )
    dispersion_outliers = [row for row in rows if row["dispersion_outlier"]]
    return {
        "deseq2_max_dispersion_bound": bound,
        "definition": "simulator_true_dispersion_strictly_above_DESeq2_maxDisp",
        "genes": len(tail),
        "truth_null_genes": len(null_tail),
        "truth_de_genes": len(de_tail),
        "raw_testable_truth_null_genes_after_cooks": len(raw_testable_null),
        "published_tested_truth_null_genes_after_independent_filtering": len(
            published_tested_null
        ),
        "false_positives": len(false_positives),
        "null_false_positive_rate_among_published_tested": (
            len(false_positives) / len(published_tested_null)
            if published_tested_null
            else None
        ),
        "false_positive_share": (
            len(false_positives) / all_false_positives if all_false_positives else 0.0
        ),
        "pooled_fdp_contribution": (
            len(false_positives) / all_discoveries if all_discoveries else 0.0
        ),
        "deseq2_dispersion_outlier_flag_count": len(dispersion_outliers),
        "deseq2_dispersion_outlier_truth_null_count": sum(
            row["truth_class"] == "null" for row in dispersion_outliers
        ),
        "deseq2_dispersion_outlier_false_positive_count": sum(
            row["truth_class"] == "null"
            and row["wald_fdr_if_on"] is not None
            and row["wald_fdr_if_on"] <= NOMINAL_FDR
            for row in dispersion_outliers
        ),
    }


def _stratified_table(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    total_false_positives = sum(
        row["truth_class"] == "null"
        and row["wald_fdr_if_on"] is not None
        and row["wald_fdr_if_on"] <= NOMINAL_FDR
        for row in rows
    )
    total_discoveries = sum(
        row["wald_fdr_if_on"] is not None and row["wald_fdr_if_on"] <= NOMINAL_FDR
        for row in rows
    )
    by_truth: list[dict[str, Any]] = []
    combined: list[dict[str, Any]] = []
    for label in DISPERSION_LABELS:
        stratum = [row for row in rows if row["dispersion_stratum"] == label]
        stratum_fp = sum(
            row["truth_class"] == "null"
            and row["wald_fdr_if_on"] is not None
            and row["wald_fdr_if_on"] <= NOMINAL_FDR
            for row in stratum
        )
        stratum_tp = sum(
            row["truth_class"] == "DE"
            and row["wald_fdr_if_on"] is not None
            and row["wald_fdr_if_on"] <= NOMINAL_FDR
            for row in stratum
        )
        discoveries = stratum_fp + stratum_tp
        combined.append(
            {
                "dispersion_stratum": label,
                "genes": len(stratum),
                "null_genes": sum(row["truth_class"] == "null" for row in stratum),
                "de_genes": sum(row["truth_class"] == "DE" for row in stratum),
                "discoveries": discoveries,
                "false_positives": stratum_fp,
                "true_positives": stratum_tp,
                "within_stratum_fdp": stratum_fp / discoveries if discoveries else 0.0,
                "false_positive_share": (
                    stratum_fp / total_false_positives if total_false_positives else 0.0
                ),
                "pooled_fdp_contribution": (
                    stratum_fp / total_discoveries if total_discoveries else 0.0
                ),
            }
        )
        for truth_class in ("null", "DE"):
            group = [row for row in stratum if row["truth_class"] == truth_class]
            raw_testable = sum(row["wald_pvalue"] is not None for row in group)
            published_tested = sum(row["wald_fdr_if_on"] is not None for row in group)
            rejected = sum(
                row["wald_fdr_if_on"] is not None
                and row["wald_fdr_if_on"] <= NOMINAL_FDR
                for row in group
            )
            by_truth.append(
                {
                    "dispersion_stratum": label,
                    "truth_class": truth_class,
                    "genes": len(group),
                    "raw_testable_genes_after_cooks": raw_testable,
                    "published_tested_genes_after_independent_filtering": (
                        published_tested
                    ),
                    "discoveries": rejected,
                    "false_positives": rejected if truth_class == "null" else 0,
                    "true_positives": rejected if truth_class == "DE" else 0,
                    "null_false_positive_rate": (
                        rejected / published_tested
                        if truth_class == "null" and published_tested
                        else None
                    ),
                    "de_true_positive_rate": (
                        rejected / len(group) if truth_class == "DE" and group else None
                    ),
                }
            )
    return by_truth, combined


def _per_seed_stratification(
    jobs: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    records: list[dict[str, Any]] = []
    for job in jobs:
        truth_rows, combined_rows = _stratified_table(job["rows"])
        null_by_label = {
            row["dispersion_stratum"]: row
            for row in truth_rows
            if row["truth_class"] == "null"
        }
        for combined in combined_rows:
            null = null_by_label[combined["dispersion_stratum"]]
            records.append(
                {
                    "replicate": job["replicate"],
                    "seed": job["seed"],
                    **combined,
                    "null_raw_testable_genes_after_cooks": null[
                        "raw_testable_genes_after_cooks"
                    ],
                    "null_published_tested_genes_after_independent_filtering": (
                        null["published_tested_genes_after_independent_filtering"]
                    ),
                    "null_false_positive_rate": null["null_false_positive_rate"],
                }
            )
    ranges: list[dict[str, Any]] = []
    for label in DISPERSION_LABELS:
        group = [row for row in records if row["dispersion_stratum"] == label]
        ranges.append(
            {
                "dispersion_stratum": label,
                "replicate_count": len(group),
                "mean_false_positives": statistics.fmean(
                    row["false_positives"] for row in group
                ),
                "minimum_false_positives": min(row["false_positives"] for row in group),
                "maximum_false_positives": max(row["false_positives"] for row in group),
                "mean_null_false_positive_rate": statistics.fmean(
                    row["null_false_positive_rate"] for row in group
                ),
                "minimum_null_false_positive_rate": min(
                    row["null_false_positive_rate"] for row in group
                ),
                "maximum_null_false_positive_rate": max(
                    row["null_false_positive_rate"] for row in group
                ),
                "mean_replicate_false_positive_share": statistics.fmean(
                    row["false_positive_share"] for row in group
                ),
                "minimum_replicate_false_positive_share": min(
                    row["false_positive_share"] for row in group
                ),
                "maximum_replicate_false_positive_share": max(
                    row["false_positive_share"] for row in group
                ),
                "mean_replicate_fdp_contribution": statistics.fmean(
                    row["pooled_fdp_contribution"] for row in group
                ),
                "minimum_replicate_fdp_contribution": min(
                    row["pooled_fdp_contribution"] for row in group
                ),
                "maximum_replicate_fdp_contribution": max(
                    row["pooled_fdp_contribution"] for row in group
                ),
            }
        )
    return records, ranges


def _wald_lrt_overlap(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for truth_scope in ("all", "null", "DE"):
        group = (
            list(rows)
            if truth_scope == "all"
            else [row for row in rows if row["truth_class"] == truth_scope]
        )
        both = 0
        wald_only = 0
        lrt_only = 0
        neither = 0
        for row in group:
            wald_called = (
                row["wald_fdr_if_on"] is not None
                and row["wald_fdr_if_on"] <= NOMINAL_FDR
            )
            lrt_called = (
                row["lrt_fdr_if_on"] is not None and row["lrt_fdr_if_on"] <= NOMINAL_FDR
            )
            both += int(wald_called and lrt_called)
            wald_only += int(wald_called and not lrt_called)
            lrt_only += int(lrt_called and not wald_called)
            neither += int(not wald_called and not lrt_called)
        union = both + wald_only + lrt_only
        records.append(
            {
                "truth_scope": truth_scope,
                "genes": len(group),
                "both_called": both,
                "wald_only": wald_only,
                "lrt_only": lrt_only,
                "neither_called": neither,
                "call_jaccard": both / union if union else 1.0,
            }
        )
    return records


def _archived_wald_parity(
    jobs: Sequence[Mapping[str, Any]], exploratory_report: Mapping[str, Any]
) -> dict[str, Any]:
    metrics = exploratory_report.get("metrics")
    if not isinstance(metrics, Mapping) or not isinstance(
        metrics.get("replicates"), list
    ):
        raise BenchmarkError("The exploratory report lacks replicate metrics.")
    archived = {
        int(item["seed"]): item
        for item in metrics["replicates"]
        if isinstance(item, Mapping)
    }
    per_seed: list[dict[str, Any]] = []
    for job in jobs:
        seed = int(job["seed"])
        expected = archived.get(seed)
        observed = job["scores"]["deseq2_wald_independent_filtering_on"]
        if expected is None:
            raise BenchmarkError(f"Archived exploratory metrics omit seed {seed}.")
        expected_tested = int(expected["status_counts"].get("tested", 0))
        checks = {
            "tested_genes": observed["tested_genes"] == expected_tested,
            "discoveries": observed["discoveries"] == expected["discoveries"],
            "true_positives": observed["true_positives"] == expected["true_positives"],
            "false_positives": observed["false_positives"]
            == expected["false_positives"],
            "false_discovery_proportion": math.isclose(
                observed["false_discovery_proportion"],
                float(expected["false_discovery_proportion"]),
                rel_tol=1e-14,
                abs_tol=1e-15,
            ),
            "all_truth_de_tpr": math.isclose(
                observed["all_truth_de_tpr"],
                float(expected["tpr"]),
                rel_tol=1e-14,
                abs_tol=1e-15,
            ),
        }
        per_seed.append(
            {
                "seed": seed,
                "all_fields_match": all(checks.values()),
                "field_checks": checks,
            }
        )
    if len(per_seed) != 20 or not all(item["all_fields_match"] for item in per_seed):
        raise BenchmarkError(
            "Direct Wald diagnostics do not reproduce the disclosed exploratory failure."
        )
    return {
        "matched_replicates": 20,
        "required_replicates": 20,
        "all_match": True,
        "fields": [
            "tested_genes",
            "discoveries",
            "true_positives",
            "false_positives",
            "false_discovery_proportion",
            "all_truth_de_tpr",
        ],
        "per_seed": per_seed,
    }


def _svg_bar_chart(
    path: Path,
    *,
    title: str,
    labels: Sequence[str],
    series: Sequence[tuple[str, Sequence[float], str]],
    y_label: str,
    y_max: float,
) -> None:
    width, height = 900, 520
    left, right, top, bottom = 90, 30, 70, 110
    plot_width = width - left - right
    plot_height = height - top - bottom
    group_width = plot_width / len(labels)
    bar_width = group_width * 0.72 / len(series)
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520" viewBox="0 0 900 520">',
        '<rect width="900" height="520" fill="white"/>',
        f'<text x="450" y="34" text-anchor="middle" font-family="sans-serif" font-size="20">{html.escape(title)}</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#222"/>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#222"/>',
    ]
    for tick in range(6):
        value = y_max * tick / 5
        y = top + plot_height - plot_height * value / y_max
        lines.extend(
            [
                f'<line x1="{left}" y1="{y:.3f}" x2="{left + plot_width}" y2="{y:.3f}" stroke="#ddd"/>',
                f'<text x="{left - 10}" y="{y + 4:.3f}" text-anchor="end" font-family="sans-serif" font-size="12">{value:.2f}</text>',
            ]
        )
    for group_index, label in enumerate(labels):
        group_x = left + group_width * group_index
        start_x = group_x + group_width * 0.14
        for series_index, (_, values, color) in enumerate(series):
            value = values[group_index]
            bar_height = plot_height * value / y_max
            x = start_x + series_index * bar_width
            y = top + plot_height - bar_height
            lines.append(
                f'<rect x="{x:.3f}" y="{y:.3f}" width="{bar_width:.3f}" height="{bar_height:.3f}" fill="{color}"/>'
            )
        center = group_x + group_width / 2
        lines.append(
            f'<text x="{center:.3f}" y="{top + plot_height + 24}" text-anchor="middle" font-family="sans-serif" font-size="12">{html.escape(label)}</text>'
        )
    lines.append(
        f'<text x="18" y="{top + plot_height / 2:.3f}" transform="rotate(-90 18 {top + plot_height / 2:.3f})" text-anchor="middle" font-family="sans-serif" font-size="14">{html.escape(y_label)}</text>'
    )
    legend_x = left
    legend_y = height - 34
    for index, (name, _, color) in enumerate(series):
        x = legend_x + index * 260
        lines.extend(
            [
                f'<rect x="{x}" y="{legend_y - 12}" width="16" height="16" fill="{color}"/>',
                f'<text x="{x + 23}" y="{legend_y + 1}" font-family="sans-serif" font-size="13">{html.escape(name)}</text>',
            ]
        )
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_plots(
    artifact_dir: Path,
    *,
    strata: Sequence[Mapping[str, Any]],
    method_summaries: Mapping[str, Mapping[str, Any]],
) -> list[Path]:
    artifact_dir.mkdir(parents=True, exist_ok=True)
    if any(artifact_dir.iterdir()):
        raise BenchmarkError("The diagnostic artifact directory must be empty.")
    strata_plot = artifact_dir / "false-positive-by-true-dispersion.svg"
    _svg_bar_chart(
        strata_plot,
        title="DESeq2 Wald false positives by within-seed true-dispersion quintile",
        labels=[str(item["dispersion_stratum"]) for item in strata],
        series=[
            (
                "share of all false positives",
                [float(item["false_positive_share"]) for item in strata],
                "#b33c3c",
            ),
            (
                "contribution to pooled FDP",
                [float(item["pooled_fdp_contribution"]) for item in strata],
                "#4666a5",
            ),
        ],
        y_label="fraction",
        y_max=0.5,
    )
    method_order = (
        "deseq2_wald_independent_filtering_on",
        "deseq2_wald_independent_filtering_off",
        "deseq2_lrt_independent_filtering_on",
        "edger_ql_native",
        "deseq2_wald_aligned_common_tested_family",
        "edger_ql_aligned_common_tested_family",
    )
    method_plot = artifact_dir / "method-fdp-tpr-comparison.svg"
    _svg_bar_chart(
        method_plot,
        title="Exploratory mean replicate FDP and TPR",
        labels=(
            "Wald IF on",
            "Wald IF off",
            "LRT IF on",
            "edgeR native",
            "Wald aligned",
            "edgeR aligned",
        ),
        series=[
            (
                "mean replicate FDP",
                [
                    float(method_summaries[name]["mean_replicate_fdp"])
                    for name in method_order
                ],
                "#b33c3c",
            ),
            (
                "mean all-truth-DE TPR",
                [
                    float(method_summaries[name]["mean_all_truth_de_tpr"])
                    for name in method_order
                ],
                "#3b7f5f",
            ),
        ],
        y_label="proportion",
        y_max=0.8,
    )
    return [strata_plot, method_plot]


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise BenchmarkError(f"Cannot write an empty diagnostic table: {path.name}.")
    fields = list(rows[0])
    if any(list(row) != fields for row in rows):
        raise BenchmarkError(
            f"Diagnostic table rows have inconsistent fields: {path.name}."
        )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: ""
                    if value is None
                    else str(value).lower()
                    if isinstance(value, bool)
                    else value
                    for key, value in row.items()
                }
            )


def _write_tables(
    artifact_dir: Path,
    *,
    truth_table: Sequence[Mapping[str, Any]],
    combined_table: Sequence[Mapping[str, Any]],
    per_seed_strata: Sequence[Mapping[str, Any]],
    method_summaries: Mapping[str, Mapping[str, Any]],
    overlap_per_seed: Sequence[Mapping[str, Any]],
    overlap_pooled: Sequence[Mapping[str, Any]],
    jobs: Sequence[Mapping[str, Any]],
) -> list[Path]:
    null_by_label = {
        row["dispersion_stratum"]: row
        for row in truth_table
        if row["truth_class"] == "null"
    }
    pooled_strata = []
    for combined in combined_table:
        null = null_by_label[combined["dispersion_stratum"]]
        pooled_strata.append(
            {
                **combined,
                "null_raw_testable_genes_after_cooks": null[
                    "raw_testable_genes_after_cooks"
                ],
                "null_published_tested_genes_after_independent_filtering": null[
                    "published_tested_genes_after_independent_filtering"
                ],
                "null_false_positive_rate": null["null_false_positive_rate"],
            }
        )
    summary_fields = (
        "method",
        "replicate_count",
        "mean_tested_genes",
        "minimum_tested_genes",
        "maximum_tested_genes",
        "mean_replicate_fdp",
        "minimum_replicate_fdp",
        "maximum_replicate_fdp",
        "pooled_fdp_descriptive_only",
        "mean_all_truth_de_tpr",
        "minimum_all_truth_de_tpr",
        "maximum_all_truth_de_tpr",
        "mean_conditional_tpr_within_tested_family",
        "minimum_conditional_tpr_within_tested_family",
        "maximum_conditional_tpr_within_tested_family",
        "mean_truth_de_in_tested_family",
        "total_discoveries",
        "total_true_positives",
        "total_false_positives",
    )
    method_rows = [
        {field: summary[field] for field in summary_fields}
        for summary in method_summaries.values()
    ]
    replicate_method_rows = [
        {"method": name, **record}
        for name, summary in method_summaries.items()
        for record in summary["per_replicate"]
    ]
    overlap_rows = [
        {"aggregation": "per_seed", **record} for record in overlap_per_seed
    ] + [
        {"aggregation": "pooled", "replicate": None, "seed": None, **record}
        for record in overlap_pooled
    ]
    fit_rows = []
    for job in jobs:
        fit = job["fit_summary"]
        fit_rows.append(
            {
                "replicate": job["replicate"],
                "seed": job["seed"],
                "wald_resolved_fit_type": fit["resolved_fit_type"],
                "lrt_resolved_fit_type": fit["lrt_resolved_fit_type"],
                "automatic_fallback": fit["automatic_fallback"],
                "wald_asymptDisp": fit["fit_coefficients"]["asymptDisp"],
                "wald_extraPois": fit["fit_coefficients"]["extraPois"],
                "lrt_asymptDisp": fit["lrt_fit_coefficients"]["asymptDisp"],
                "lrt_extraPois": fit["lrt_fit_coefficients"]["extraPois"],
                "max_dispersion_bound": fit["max_dispersion_bound"],
                "dispersion_outlier_count": fit["dispersion_outlier_count"],
                "true_dispersion_above_max_count": fit[
                    "true_dispersion_above_max_count"
                ],
            }
        )
    tables = {
        "dispersion-strata-pooled.tsv": pooled_strata,
        "dispersion-strata-by-seed.tsv": list(per_seed_strata),
        "method-summary.tsv": method_rows,
        "method-summary-by-seed.tsv": replicate_method_rows,
        "wald-lrt-call-overlap.tsv": overlap_rows,
        "dispersion-fit-by-seed.tsv": fit_rows,
    }
    paths: list[Path] = []
    for name, rows in tables.items():
        path = artifact_dir / name
        _write_tsv(path, rows)
        paths.append(path)
    return paths


def _implementation_inventory(
    runner: Path, diagnostic: Path, generator: Path
) -> list[dict[str, object]]:
    paths = (runner, diagnostic, generator, runner.with_name("common.py"))
    return [file_evidence(path, name=path.name) for path in paths]


def _paired_differences(
    left: Mapping[str, Any], right: Mapping[str, Any]
) -> dict[str, Any]:
    left_by_seed = {item["seed"]: item for item in left["per_replicate"]}
    right_by_seed = {item["seed"]: item for item in right["per_replicate"]}
    if left_by_seed.keys() != right_by_seed.keys():
        raise BenchmarkError("Method comparison seed grids differ.")
    fdp = [
        left_by_seed[seed]["false_discovery_proportion"]
        - right_by_seed[seed]["false_discovery_proportion"]
        for seed in sorted(left_by_seed)
    ]
    all_truth_tpr = [
        left_by_seed[seed]["all_truth_de_tpr"] - right_by_seed[seed]["all_truth_de_tpr"]
        for seed in sorted(left_by_seed)
    ]
    conditional_tpr = [
        left_by_seed[seed]["conditional_tpr_within_tested_family"]
        - right_by_seed[seed]["conditional_tpr_within_tested_family"]
        for seed in sorted(left_by_seed)
    ]
    return {
        "left_minus_right": True,
        "mean_paired_fdp_difference": statistics.fmean(fdp),
        "minimum_paired_fdp_difference": min(fdp),
        "maximum_paired_fdp_difference": max(fdp),
        "mean_paired_all_truth_de_tpr_difference": statistics.fmean(all_truth_tpr),
        "minimum_paired_all_truth_de_tpr_difference": min(all_truth_tpr),
        "maximum_paired_all_truth_de_tpr_difference": max(all_truth_tpr),
        "mean_paired_conditional_tpr_difference": statistics.fmean(conditional_tpr),
        "minimum_paired_conditional_tpr_difference": min(conditional_tpr),
        "maximum_paired_conditional_tpr_difference": max(conditional_tpr),
    }


def _format_number(value: Any, digits: int = 5) -> str:
    if value is None:
        return "not applicable"
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.{digits}f}"


def _render_markdown(report: Mapping[str, Any]) -> str:
    methods = report["method_comparisons"]
    dispersion = report["dispersion_fit"]
    strata = report["truth_dispersion_stratification"]
    exploratory = report["inputs"][0]
    resolved_fit_counts = "; ".join(
        f"{mode.upper()} {dict(sorted(dispersion['resolved_fit_type_counts'][mode].items()))}"
        for mode in ("wald", "lrt")
    )
    lines = [
        "# DESeq2 compcodeR calibration-failure mechanism diagnostic",
        "",
        "## Scope",
        "",
        "This is a read-only, hypothesis-generating diagnostic of the disclosed exploratory seeds 61001–61020. It is not certification evidence and does not change any gate threshold or DESeq2 setting. Its declared inputs and executable seed enumeration contain only the exploratory grid; because the simulation seeds are public and deterministic, historical non-viewing of other seeds cannot be cryptographically proved.",
        "",
        f"The diagnostic is bound to exploratory report SHA-256 `{exploratory['sha256']}`. Each fixture was regenerated and matched to that report before analysis. The adjacent JSON is the authoritative machine record; this Markdown is a deterministic derived view.",
        "",
        "## Dispersion fit",
        "",
        f"DESeq2 requested the parametric dispersion trend in 20/20 replicates. The resolved fit counts were `{resolved_fit_counts}`; automatic fallback occurred in {dispersion['automatic_fallback_count']} replicates. Because no fallback was observed, a fallback-effect contrast is not estimable from this seed grid.",
        "",
        dispersion["numerical_success_caveat"],
        "",
        f"The simulator true-dispersion tail above DESeq2's maxDisp={_format_number(dispersion['extreme_true_dispersion_tail']['deseq2_max_dispersion_bound'], 1)} contains {dispersion['extreme_true_dispersion_tail']['truth_null_genes']} truth-null genes; {dispersion['extreme_true_dispersion_tail']['published_tested_truth_null_genes_after_independent_filtering']} were published-tested and {dispersion['extreme_true_dispersion_tail']['false_positives']} were false positives.",
        "",
        "Final-to-true dispersion ratios:",
        "",
        "| Group | Genes | Median | Q25 | Q75 | Fraction below 1 |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for name in (
        "all_genes",
        "truth_null",
        "truth_de",
        "null_false_positive",
        "null_not_called",
    ):
        values = dispersion["estimation_ratio_summaries"][name]["final_to_true_ratio"]
        lines.append(
            f"| {name.replace('_', ' ')} | {values['count']} | {_format_number(values['median'])} | {_format_number(values['q25'])} | {_format_number(values['q75'])} | {_format_number(values['below_one_fraction'])} |"
        )
    lines.extend(
        [
            "",
            "The false-positive row conditions on the observed outcome and is therefore selection-biased. Its contrast with uncalled null genes is hypothesis-generating association, not a causal estimate.",
            "",
            "## False positives by true dispersion",
            "",
            "Equal-count quintiles are formed independently within each seed by ranking true dispersion with gene_id as a deterministic tie-break. They never use fitted statistics, calls, P values, or observed-result ranks. `Null false-positive rate` is FP / published-tested truth-null genes after independent filtering in the stratum. `FP share` is FP / all false positives. `Pooled-FDP contribution` is FP / all discoveries and is not mislabeled as a within-stratum FDP.",
            "",
            "| Within-seed quintile | Null raw-testable | Null published-tested | False positives | Null false-positive rate | FP share | Pooled-FDP contribution |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    truth_rows = {
        item["dispersion_stratum"]: item
        for item in strata["truth_by_stratum"]
        if item["truth_class"] == "null"
    }
    for item in strata["combined_by_stratum"]:
        null = truth_rows[item["dispersion_stratum"]]
        lines.append(
            f"| {item['dispersion_stratum']} | {null['raw_testable_genes_after_cooks']} | {null['published_tested_genes_after_independent_filtering']} | {item['false_positives']} | {_format_number(null['null_false_positive_rate'])} | {_format_number(item['false_positive_share'])} | {_format_number(item['pooled_fdp_contribution'])} |"
        )
    lines.extend(
        [
            "",
            "Per-seed concentration ranges:",
            "",
            "| Quintile | Mean FP | FP range | Mean null FPR | Null FPR range |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
    )
    for item in strata["replicate_ranges_by_stratum"]:
        lines.append(
            f"| {item['dispersion_stratum']} | {_format_number(item['mean_false_positives'], 2)} | {item['minimum_false_positives']}–{item['maximum_false_positives']} | {_format_number(item['mean_null_false_positive_rate'])} | {_format_number(item['minimum_null_false_positive_rate'])}–{_format_number(item['maximum_null_false_positive_rate'])} |"
        )
    lines.extend(
        [
            "",
            "![False positives by true dispersion](deseq2-compcoder-mechanism-artifacts/false-positive-by-true-dispersion.svg)",
            "",
            "## Method diagnostics",
            "",
            "| Analysis | Mean tested | Mean replicate FDP | Pooled FDP | Mean all-truth-DE TPR | Mean conditional TPR in tested family |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    order = (
        "deseq2_wald_independent_filtering_on",
        "deseq2_wald_independent_filtering_off",
        "deseq2_lrt_independent_filtering_on",
        "edger_ql_native",
        "deseq2_wald_aligned_common_tested_family",
        "edger_ql_aligned_common_tested_family",
    )
    for name in order:
        values = methods["summaries"][name]
        lines.append(
            f"| {name.replace('_', ' ')} | {_format_number(values['mean_tested_genes'], 1)} | {_format_number(values['mean_replicate_fdp'])} | {_format_number(values['pooled_fdp_descriptive_only'])} | {_format_number(values['mean_all_truth_de_tpr'])} | {_format_number(values['mean_conditional_tpr_within_tested_family'])} |"
        )
    pooled_overlap = {
        item["truth_scope"]: item
        for item in methods["wald_lrt_call_overlap"]["pooled_descriptive_only"]
    }
    lines.extend(
        [
            "",
            "![Method comparison](deseq2-compcoder-mechanism-artifacts/method-fdp-tpr-comparison.svg)",
            "",
            "Independent filtering is compared using the same Wald fit and the same Cook's-cutoff behavior; only `independentFiltering` changes. Wald and LRT use the same counts, design, size-factor/dispersion settings, and seed. The LRT tests the omnibus full `~ condition` versus reduced `~ 1` hypothesis.",
            "",
            "Pooled Wald/LRT call overlap (descriptive only):",
            "",
            "| Truth scope | Both | Wald only | LRT only | Neither | Jaccard |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
            *[
                f"| {scope} | {pooled_overlap[scope]['both_called']} | {pooled_overlap[scope]['wald_only']} | {pooled_overlap[scope]['lrt_only']} | {pooled_overlap[scope]['neither_called']} | {_format_number(pooled_overlap[scope]['call_jaccard'])} |"
                for scope in ("all", "null", "DE")
            ],
            "",
            "For the aligned edgeR comparison, both raw-P families are restricted to the intersection of finite DESeq2 Wald P values and genes retained by edgeR `filterByExpr`; BH is then recomputed separately for each method on those identical gene IDs. Native results are retained only as unaligned context.",
            "",
            "## Interpretation boundary",
            "",
            "These 20 seeds may generate hypotheses about the observed calibration failure. They cannot establish a new applicability boundary, support changing a gate, or certify either method. The self-attested synthetic matrix also does not authenticate featureCounts provenance. No held-out artifacts are declared or consumed by this diagnostic; opening a held-out route still requires a separate user decision.",
            "",
        ]
    )
    return "\n".join(lines)


def run(arguments: argparse.Namespace) -> dict[str, Any]:
    runner = Path(__file__).resolve()
    diagnostic = runner.with_name("run_deseq2_compcoder_mechanism_diagnostic.R")
    generator = runner.with_name("generate_deseq2_compcoder_fixture.R")
    exploratory_path = (
        PROJECT_ROOT / "tests/simulation/deseq2-compcoder-exploratory-report.json"
    )
    rscript = rscript_from_runtime(arguments.rscript)
    r_library = arguments.r_library.resolve(strict=True)
    exploratory_report, archived_inputs = _validate_exploratory_report(exploratory_path)
    if platform.python_version() != "3.12.12":
        raise BenchmarkError("The diagnostic requires locked Python 3.12.12.")
    if Path(sys.executable).resolve().parent != rscript.parent:
        raise BenchmarkError("Python and Rscript must be from the same locked prefix.")
    with _workspace(arguments.workspace) as workspace_value:
        workspace = Path(workspace_value)
        worker = lambda replicate: _run_replicate(  # noqa: E731
            replicate,
            workspace=workspace,
            generator=generator,
            diagnostic=diagnostic,
            rscript=rscript,
            r_library=r_library,
            archived_inputs=archived_inputs,
        )
        with ThreadPoolExecutor(max_workers=WORKERS) as executor:
            jobs = list(executor.map(worker, REPLICATES))

    jobs.sort(key=lambda item: item["replicate"])
    if [job["seed"] for job in jobs] != list(EXPLORATORY_SEEDS):
        raise BenchmarkError("The completed diagnostic seed grid is not exploratory.")
    all_rows = [row for job in jobs for row in job["rows"]]
    fit_types: dict[str, dict[str, int]] = {
        "wald": defaultdict(int),
        "lrt": defaultdict(int),
    }
    for job in jobs:
        fit_types["wald"][str(job["fit_summary"]["resolved_fit_type"])] += 1
        fit_types["lrt"][str(job["fit_summary"]["lrt_resolved_fit_type"])] += 1
    fallback_count = sum(bool(job["fit_summary"]["automatic_fallback"]) for job in jobs)
    fit_coefficients = {
        mode: {
            name: [
                float(job["fit_summary"][field][name])
                for job in jobs
                if job["fit_summary"][field] is not None
            ]
            for name in ("asymptDisp", "extraPois")
        }
        for mode, field in (
            ("wald", "fit_coefficients"),
            ("lrt", "lrt_fit_coefficients"),
        )
    }
    fit_pair_groups: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for job in jobs:
        fit_pair_groups[
            (
                str(job["fit_summary"]["resolved_fit_type"]),
                str(job["fit_summary"]["lrt_resolved_fit_type"]),
            )
        ].append(job)
    fit_impact_groups = [
        {
            "wald_resolved_fit_type": pair[0],
            "lrt_resolved_fit_type": pair[1],
            "replicates": [job["replicate"] for job in group],
            "wald_scores": _aggregate_scores(
                group, "deseq2_wald_independent_filtering_on"
            ),
            "lrt_scores": _aggregate_scores(
                group, "deseq2_lrt_independent_filtering_on"
            ),
        }
        for pair, group in sorted(fit_pair_groups.items())
    ]
    method_names = tuple(jobs[0]["scores"])
    method_summaries = {
        method: _aggregate_scores(jobs, method) for method in method_names
    }
    truth_table, combined_table = _stratified_table(all_rows)
    per_seed_strata, replicate_strata_ranges = _per_seed_stratification(jobs)
    overlap_per_seed = [
        {"replicate": job["replicate"], "seed": job["seed"], **record}
        for job in jobs
        for record in _wald_lrt_overlap(job["rows"])
    ]
    overlap_pooled = _wald_lrt_overlap(all_rows)
    archived_parity = _archived_wald_parity(jobs, exploratory_report)
    artifact_paths = _write_plots(
        arguments.artifact_dir,
        strata=combined_table,
        method_summaries=method_summaries,
    )
    artifact_paths.extend(
        _write_tables(
            arguments.artifact_dir,
            truth_table=truth_table,
            combined_table=combined_table,
            per_seed_strata=per_seed_strata,
            method_summaries=method_summaries,
            overlap_per_seed=overlap_per_seed,
            overlap_pooled=overlap_pooled,
            jobs=jobs,
        )
    )
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "diagnostic_id": DIAGNOSTIC_ID,
        "status": "complete",
        "evidence_role": "hypothesis_generation_only_not_certification",
        "scope": {
            "seed_grid": list(EXPLORATORY_SEEDS),
            "seed_role": "disclosed_exploratory_only",
            "held_out_seeds_consumed_by_this_diagnostic": False,
            "gate_thresholds_changed": False,
            "fixture_changed": False,
            "backend_semantics_changed": False,
            "analysis_is_read_only": True,
        },
        "runtime": {
            "python": platform.python_version(),
            "python_executable": Path(sys.executable).name,
            "rscript_executable": rscript.name,
            "rscript_sha256": sha256_file(rscript),
            "r_library_identity": EXPECTED_RUNTIME,
            "workers": WORKERS,
        },
        "implementation": _implementation_inventory(runner, diagnostic, generator),
        "inputs": [
            file_evidence(
                exploratory_path, name="deseq2-compcoder-exploratory-report.json"
            ),
            file_evidence(r_library / "DESeq2/DESCRIPTION", name="DESeq2/DESCRIPTION"),
            file_evidence(r_library / "edgeR/DESCRIPTION", name="edgeR/DESCRIPTION"),
            file_evidence(
                r_library / "compcodeR/DESCRIPTION", name="compcodeR/DESCRIPTION"
            ),
            *(record for job in jobs for record in job["fixture_inputs"]),
        ],
        "artifacts": [
            file_evidence(
                path, name=f"deseq2-compcoder-mechanism-artifacts/{path.name}"
            )
            for path in artifact_paths
        ],
        "methods": {
            "deseq2_common": {
                "version": "1.52.0",
                "design": "~ condition",
                "sfType": "ratio",
                "fitType_requested": "parametric",
                "betaPrior": False,
                "minReplicatesForReplace": 7,
                "useT": False,
                "minmu": 0.5,
                "lfcThreshold": 0,
                "altHypothesis": "greaterAbs",
                "alpha_for_independent_filter_optimization": 0.1,
                "pAdjustMethod": "BH",
                "cooksCutoff": "automatic",
            },
            "lrt": {
                "full": "~ condition",
                "reduced": "~ 1",
                "hypothesis": "omnibus_full_vs_reduced",
            },
            "edger_ql": {
                "version": "4.10.0",
                "steps": ["filterByExpr", "normLibSizes_TMM", "glmQLFit", "glmQLFTest"],
                "abundance_trend": True,
                "robust": True,
                "winsor_tail_p": [0.05, 0.1],
                "legacy": False,
                "top_proportion": None,
            },
            "aligned_comparison": {
                "family": "intersection_of_finite_DESeq2_Wald_raw_P_and_edgeR_filterByExpr_genes",
                "adjustment": "separate_BH_per_method_on_identical_gene_ids",
            },
            "independent_filtering_comparison": {
                "isolated_argument": "independentFiltering_TRUE_vs_FALSE",
                "unchanged": [
                    "same_fitted_Wald_DESeqDataSet",
                    "same_contrast",
                    "same_alpha_0.1",
                    "same_BH_adjustment",
                    "same_automatic_Cooks_cutoff",
                ],
            },
        },
        "dispersion_fit": {
            "requested_fit_type": "parametric",
            "resolved_fit_type_counts": {
                mode: dict(sorted(counts.items())) for mode, counts in fit_types.items()
            },
            "automatic_fallback_count": fallback_count,
            "fallback_impact": {
                "comparison_status": (
                    "not_estimable_no_automatic_fallback_observed"
                    if fallback_count == 0
                    else "fallback_observed_requires_explicit_group_comparison"
                ),
                "fallback_replicates": [
                    job["replicate"]
                    for job in jobs
                    if job["fit_summary"]["automatic_fallback"]
                ],
                "score_groups_by_resolved_fit_pair": fit_impact_groups,
            },
            "numerical_success_caveat": (
                "A resolved parametric curve with finite positive coefficients "
                "establishes numerical completion, not dispersion-model adequacy "
                "or causality for false discoveries."
            ),
            "parametric_coefficient_ranges": {
                mode: {
                    name: {"minimum": min(values), "maximum": max(values)}
                    for name, values in coefficients.items()
                }
                for mode, coefficients in fit_coefficients.items()
            },
            "per_replicate": [
                {
                    "replicate": job["replicate"],
                    "seed": job["seed"],
                    "requested_fit_type": job["fit_summary"]["requested_fit_type"],
                    "resolved_fit_type": job["fit_summary"]["resolved_fit_type"],
                    "lrt_resolved_fit_type": job["fit_summary"][
                        "lrt_resolved_fit_type"
                    ],
                    "automatic_fallback": job["fit_summary"]["automatic_fallback"],
                    "fit_coefficients": job["fit_summary"]["fit_coefficients"],
                    "lrt_fit_coefficients": job["fit_summary"]["lrt_fit_coefficients"],
                    "dispersion_outlier_count": job["fit_summary"][
                        "dispersion_outlier_count"
                    ],
                    "wald_notices": job["fit_summary"]["wald_notices"],
                    "lrt_notices": job["fit_summary"]["lrt_notices"],
                }
                for job in jobs
            ],
            "estimation_ratio_summaries": _dispersion_estimation(all_rows),
            "extreme_true_dispersion_tail": _extreme_dispersion_tail(all_rows, jobs),
        },
        "truth_dispersion_stratification": {
            "definition": (
                "within_seed_equal_count_rank_quintiles_of_true_dispersion_"
                "with_gene_id_tie_break"
            ),
            "labels": list(DISPERSION_LABELS),
            "per_seed_rank_boundaries": [
                {
                    "replicate": job["replicate"],
                    "seed": job["seed"],
                    "boundaries": job["stratum_boundaries"],
                }
                for job in jobs
            ],
            "denominators": {
                "null_false_positive_rate": (
                    "false_positives_in_stratum/published_tested_truth_null_"
                    "genes_after_independent_filtering_in_stratum"
                ),
                "false_positive_share": "false_positives_in_stratum/all_false_positives",
                "pooled_fdp_contribution": "false_positives_in_stratum/all_discoveries",
                "within_stratum_fdp": "false_positives_in_stratum/discoveries_in_stratum",
            },
            "truth_by_stratum": truth_table,
            "combined_by_stratum": combined_table,
            "per_seed_by_stratum": per_seed_strata,
            "replicate_ranges_by_stratum": replicate_strata_ranges,
        },
        "method_comparisons": {
            "summaries": method_summaries,
            "archived_exploratory_wald_parity": archived_parity,
            "wald_lrt_call_overlap": {
                "per_seed": overlap_per_seed,
                "pooled_descriptive_only": overlap_pooled,
            },
            "independent_filtering_on_minus_off": _paired_differences(
                method_summaries["deseq2_wald_independent_filtering_on"],
                method_summaries["deseq2_wald_independent_filtering_off"],
            ),
            "wald_minus_lrt": _paired_differences(
                method_summaries["deseq2_wald_independent_filtering_on"],
                method_summaries["deseq2_lrt_independent_filtering_on"],
            ),
            "native_deseq2_minus_edger": _paired_differences(
                method_summaries["deseq2_wald_independent_filtering_on"],
                method_summaries["edger_ql_native"],
            ),
            "aligned_deseq2_minus_edger": _paired_differences(
                method_summaries["deseq2_wald_aligned_common_tested_family"],
                method_summaries["edger_ql_aligned_common_tested_family"],
            ),
        },
        "limitations": [
            "exploratory seeds are reused for hypothesis generation and cannot certify a revised claim",
            (
                "declared diagnostic inputs contain no held-out artifacts; because "
                "simulation seeds are public and deterministic, historical non-viewing "
                "cannot be cryptographically proved"
            ),
            "the finite balanced 6-vs-6 independent-gene simulation does not represent all bulk RNA-seq settings",
            "the self-attested featureCounts manifest does not authenticate producer origin",
            "no fallback-effect comparison is possible when every replicate resolves the parametric fit",
            "dispersion ratios conditioned on false-positive status are selection-biased associations, not causal estimates",
        ],
    }
    write_json(arguments.report, report)
    arguments.markdown.parent.mkdir(parents=True, exist_ok=True)
    arguments.markdown.write_text(_render_markdown(report), encoding="utf-8")
    return report


def main() -> int:
    try:
        report = run(_arguments())
    except (BenchmarkError, OSError, ValueError) as error:
        print(f"DESeq2 mechanism diagnostic failed: {error}", file=sys.stderr)
        return 1
    print(
        f"DESeq2 mechanism diagnostic complete: {report['diagnostic_id']}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
