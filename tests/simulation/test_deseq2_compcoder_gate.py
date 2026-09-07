"""DESeq2 compcodeR exploratory evidence and fail-closed gate mechanics."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType

import pytest

from scripts.benchmark.evidence_resolver import verify_archived_implementation

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUNNER = PROJECT_ROOT / "scripts/benchmark/run_deseq2_compcoder_gate.py"
EXPLORATORY_REPORT = Path(__file__).with_name(
    "deseq2-compcoder-exploratory-report.json"
)
EXPLORATORY_ID = "compcoder-deseq2-nb-exploratory-v1"
IMPLEMENTATION_PATHS = {
    "run_deseq2_compcoder_gate.py": RUNNER,
    "common.py": PROJECT_ROOT / "scripts/benchmark/common.py",
    "benchmark-report-v1.schema.json": (
        PROJECT_ROOT / "scripts/benchmark/benchmark-report-v1.schema.json"
    ),
    "generate_deseq2_compcoder_fixture.R": (
        PROJECT_ROOT / "scripts/benchmark/generate_deseq2_compcoder_fixture.R"
    ),
    "cli.py": PROJECT_ROOT / "rnaseq_downstream/cli.py",
    "input_semantics.py": PROJECT_ROOT / "rnaseq_downstream/input_semantics.py",
    "validation_bundle.py": PROJECT_ROOT / "rnaseq_downstream/validation_bundle.py",
    "analysis_contract_v12.py": (
        PROJECT_ROOT / "rnaseq_downstream/analysis_contract_v12.py"
    ),
    "deseq2_backend.py": PROJECT_ROOT / "rnaseq_downstream/deseq2_backend.py",
    "run_summary.py": PROJECT_ROOT / "rnaseq_downstream/run_summary.py",
    "deseq2.R": PROJECT_ROOT / "rnaseq_downstream/r_scripts/deseq2.R",
    "conda-lock.yml": PROJECT_ROOT / "conda-lock.yml",
    "renv.lock": PROJECT_ROOT / "renv.lock",
    "environment.p0.yml": PROJECT_ROOT / "environment.p0.yml",
    "environment/r-sources.lock": PROJECT_ROOT / "environment/r-sources.lock",
    "environment/verify.R": PROJECT_ROOT / "environment/verify.R",
}


def _load_runner() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "deseq2_compcoder_runner", RUNNER
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.path.insert(0, str(RUNNER.parent))
    try:
        specification.loader.exec_module(module)
    finally:
        sys.path.remove(str(RUNNER.parent))
    return module


def _strict_json(path: Path) -> dict[str, object]:
    def reject_constant(value: str) -> None:
        raise ValueError(f"non-finite JSON constant: {value}")

    def unique_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    value = json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=reject_constant,
        object_pairs_hook=unique_pairs,
    )
    assert isinstance(value, dict)
    return value


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


@pytest.mark.simulation
def test_archived_exploratory_report_is_disclosed_but_not_certification() -> None:
    runner = _load_runner()
    report = _strict_json(EXPLORATORY_REPORT)
    runner.assert_report_shape(report, benchmark_id=EXPLORATORY_ID)
    assert report["status"] == "pass"
    assert report["metrics"]["scope"] == (  # type: ignore[index]
        "public_cli_inspect_validate_run_summarize"
    )
    assert report["metrics"]["evidence_role"] == (  # type: ignore[index]
        "threshold_selection_only_not_certification"
    )
    assert report["metrics"]["seed_grid"] == list(range(61001, 61021))  # type: ignore[index]
    assert report["thresholds"] == {  # type: ignore[comparison-overlap]
        "applied_to_this_report": False,
        "held_out_policy": (
            "certification_seeds_are_disjoint_and_thresholds_are_not_revised_"
            "after_certification_results_are_observed"
        ),
        "maximum_mean_replicate_fdp": 0.065,
        "maximum_replicate_fdp": 0.12,
        "minimum_mean_tpr": 0.45,
        "minimum_replicate_tpr": 0.35,
        "nominal_bh_fdr": 0.05,
        "selection_context": (
            "candidate_limits_recorded_before_held_out_certification"
        ),
        "threshold_basis": (
            "nominal_0.05_BH_plus_finite_simulation_tolerance_and_"
            "conservative_DESeq2_power_regression_floor"
        ),
    }
    assertion_names = {item["name"] for item in report["assertions"]}  # type: ignore[union-attr]
    assert all(item["passed"] for item in report["assertions"])  # type: ignore[union-attr]
    assert not assertion_names & {
        "mean_replicate_fdp_gate",
        "worst_replicate_fdp_gate",
        "mean_tpr_gate",
        "worst_replicate_tpr_gate",
    }
    assert len(report["inputs"]) == 2 + 20 * 8  # type: ignore[arg-type]
    assert len(report["artifacts"]) == 20 * 3  # type: ignore[arg-type]
    serialized = EXPLORATORY_REPORT.read_text(encoding="utf-8")
    assert "/tmp/" not in serialized
    assert "/home/" not in serialized


@pytest.mark.simulation
def test_exploratory_report_binds_the_exact_executed_implementation() -> None:
    report = _strict_json(EXPLORATORY_REPORT)
    verify_archived_implementation(report["implementation"], IMPLEMENTATION_PATHS)


def test_exploratory_and_held_out_seed_grids_are_disjoint() -> None:
    runner = _load_runner()
    exploratory = {
        runner.SEED_BASES["exploratory"] + replicate for replicate in runner.REPLICATES
    }
    certification = {
        runner.SEED_BASES["certification"] + replicate
        for replicate in runner.REPLICATES
    }
    assert exploratory == set(range(61001, 61021))
    assert certification == set(range(62001, 62021))
    assert exploratory.isdisjoint(certification)


def test_scoring_uses_truth_nulls_and_explicit_statuses(tmp_path: Path) -> None:
    runner = _load_runner()
    fixture = tmp_path / "fixture"
    output = tmp_path / "output"
    _write_tsv(
        fixture / "truth.tsv",
        [
            {
                "gene_id": "true_called",
                "differential_expression": 1,
                "true_log2_fold_change": 2,
            },
            {
                "gene_id": "null_called",
                "differential_expression": 0,
                "true_log2_fold_change": 0,
            },
            {
                "gene_id": "null_filtered",
                "differential_expression": 0,
                "true_log2_fold_change": 0,
            },
            {
                "gene_id": "true_failed",
                "differential_expression": 1,
                "true_log2_fold_change": -2,
            },
        ],
    )
    common = {
        "contrast_id": "treated_vs_control",
        "statistic_type": "Wald",
        "statistic_hypothesis": "contrast_equals_zero",
        "fdr_basis": "contrast_pvalue_BH_after_independent_filtering",
        "test_method": "DESeq2::results_Wald",
        "lfc_threshold": 0,
        "shrinkage_method": "none",
    }
    _write_tsv(
        output / "results.tsv",
        [
            {
                "gene_id": "true_called",
                **common,
                "status": "tested",
                "status_reason": "tested",
                "PValue": 0.001,
                "FDR": 0.01,
            },
            {
                "gene_id": "null_called",
                **common,
                "status": "tested",
                "status_reason": "tested",
                "PValue": 0.002,
                "FDR": 0.02,
            },
            {
                "gene_id": "null_filtered",
                **common,
                "status": "filtered",
                "status_reason": "independent_filtering",
                "PValue": 0.8,
                "FDR": "",
            },
            {
                "gene_id": "true_failed",
                **common,
                "status": "failed",
                "status_reason": "cooks_outlier",
                "PValue": "",
                "FDR": "",
            },
        ],
    )
    stages = {
        "inspect": "success",
        "validate": "success",
        "run": "success",
        "summarize": "verified_complete",
    }

    score = runner._score_replicate(
        "exploratory", 1, fixture, output, cli_stages=stages
    )

    assert score["discoveries"] == 2
    assert score["true_positives"] == 1
    assert score["false_positives"] == 1
    assert score["false_discovery_proportion"] == 0.5
    assert score["tpr"] == 0.5
    assert score["status_counts"] == {
        "failed": 1,
        "filtered": 1,
        "tested": 2,
    }


def test_aggregate_keeps_predeclared_bad_metrics_failed() -> None:
    runner = _load_runner()
    stages = {
        "inspect": "success",
        "validate": "success",
        "run": "success",
        "summarize": "verified_complete",
    }
    metrics = [
        {
            "replicate": replicate,
            "seed": 62000 + replicate,
            "cli_stages": stages,
            "true_differential": 500,
            "discoveries": 400,
            "true_positives": 200,
            "false_positives": 200,
            "false_discovery_proportion": 0.5,
            "tpr": 0.4,
        }
        for replicate in range(1, 21)
    ]

    assertions, _ = runner._aggregate(
        "certification",
        metrics,
        runtime_identities=[dict(runner.EXPECTED_RUNTIME) for _ in range(20)],
        fixture_determinism={name: True for name in runner.FIXTURE_MEMBERS},
    )
    by_name = {item["name"]: item for item in assertions}

    assert not by_name["mean_replicate_fdp_gate"]["passed"]
    assert not by_name["worst_replicate_fdp_gate"]["passed"]
    assert not by_name["mean_tpr_gate"]["passed"]
    assert by_name["worst_replicate_tpr_gate"]["passed"]


def test_certification_refuses_tampered_exploratory_thresholds(
    tmp_path: Path,
) -> None:
    runner = _load_runner()
    document = _strict_json(EXPLORATORY_REPORT)
    document["thresholds"]["maximum_mean_replicate_fdp"] = 0.2  # type: ignore[index]
    tampered = tmp_path / "exploratory.json"
    runner.write_json(tampered, document)

    with pytest.raises(runner.BenchmarkError, match="differ from code"):
        runner._validate_exploratory_report(tampered)


def test_runner_has_no_private_backend_shortcut() -> None:
    source = RUNNER.read_text(encoding="utf-8")
    assert "_run_edger_ql_benchmark_kernel" not in source
    assert "_run_validated_deseq2" not in source
    assert '[sys.executable, "-m", "rnaseq_downstream", *arguments]' in source
