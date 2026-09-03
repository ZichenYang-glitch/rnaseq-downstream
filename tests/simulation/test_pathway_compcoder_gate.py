"""Locked compcodeR pathway-level FDR/TPR regression gate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import ModuleType

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUNNER = PROJECT_ROOT / "scripts/benchmark/run_pathway_compcoder_gate.py"
ARCHIVED_REPORT = Path(__file__).with_name("pathway-compcoder-benchmark-report.json")
BENCHMARK_ID = "compcoder-limma-self-contained-fdr-tpr-v1"
EVIDENCE_INVENTORY_SHA256 = (
    "8cf97c955739fcd5e962ceedf5ced92a01e32786a2161639f8988eeac6441141"
)
IMPLEMENTATION_PATHS = {
    "run_pathway_compcoder_gate.py": RUNNER,
    "common.py": PROJECT_ROOT / "scripts/benchmark/common.py",
    "benchmark-report-v1.schema.json": (
        PROJECT_ROOT / "scripts/benchmark/benchmark-report-v1.schema.json"
    ),
    "generate_pathway_compcoder_fixture.R": (
        PROJECT_ROOT / "scripts/benchmark/generate_pathway_compcoder_fixture.R"
    ),
    "edger_backend.py": PROJECT_ROOT / "rnaseq_downstream/edger_backend.py",
    "analysis_contract.py": PROJECT_ROOT / "rnaseq_downstream/analysis_contract.py",
    "gene_sets.py": PROJECT_ROOT / "rnaseq_downstream/gene_sets.py",
    "edger_ql.R": PROJECT_ROOT / "rnaseq_downstream/r_scripts/edger_ql.R",
    "pathway_tests.R": PROJECT_ROOT / "rnaseq_downstream/r_scripts/pathway_tests.R",
    "conda-lock.yml": PROJECT_ROOT / "conda-lock.yml",
    "renv.lock": PROJECT_ROOT / "renv.lock",
    "environment.p0.yml": PROJECT_ROOT / "environment.p0.yml",
    "r-sources.lock": PROJECT_ROOT / "environment/r-sources.lock",
    "verify.R": PROJECT_ROOT / "environment/verify.R",
}


def _load_common() -> ModuleType:
    path = PROJECT_ROOT / "scripts/benchmark/common.py"
    specification = importlib.util.spec_from_file_location(
        "pathway_compcoder_common", path
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _load_runner() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "pathway_compcoder_runner", RUNNER
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.path.insert(0, str(RUNNER.parent))
    try:
        specification.loader.exec_module(module)
    finally:
        sys.path.remove(str(RUNNER.parent))
    return module


def _write_scoring_fixture(
    output: Path,
    *,
    scenario: str,
    significant_null_sets: int = 0,
) -> None:
    runner = _load_runner()
    set_ids = runner._expected_set_ids(scenario)
    header = [
        "contrast_id", "gene_set_id", "method_id", "hypothesis", "status",
        "direction", "p_value", "fdr", "fdr_family_id", "tested_gene_count",
        "method_ngenes",
    ]
    rows: list[list[object]] = []
    for set_index, set_id in enumerate(set_ids):
        if "_UP_" in set_id:
            direction = "Up"
        else:
            direction = "Down"
        is_positive = set_id.startswith("POSITIVE_")
        is_called = is_positive or (
            set_id.startswith("NULL_") and set_index < significant_null_sets
        )
        for method, hypotheses in (
            ("limma_camera", ("directional",)),
            ("limma_fry", ("directional", "mixed")),
            ("limma_mroast", ("directional", "mixed")),
        ):
            for hypothesis in hypotheses:
                rows.append([
                    "treated_vs_control", set_id, method, hypothesis, "tested",
                    direction if hypothesis == "directional" else "Mixed",
                    0.001 if is_called else 0.5,
                    0.01 if is_called else 0.5,
                    f"treated_vs_control|{method}|{hypothesis}", 20, 20,
                ])
    output.mkdir()
    with (output / "pathway_results.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


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


@pytest.mark.simulation
def test_archived_pathway_compcoder_report_is_passing_evidence() -> None:
    report = _strict_json(ARCHIVED_REPORT)
    _load_common().assert_report_shape(report, benchmark_id=BENCHMARK_ID)
    assert report["schema_version"] == "rnaseq-downstream-benchmark-report-v1"
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    thresholds = report["thresholds"]
    assert thresholds["nominal_bh_fdr"] == 0.05  # type: ignore[index]
    assert thresholds["maximum_mean_mixed_fdp"] == 0.10  # type: ignore[index]
    assert thresholds["maximum_worst_mixed_fdp"] == 0.25  # type: ignore[index]
    assert thresholds["minimum_mean_mixed_tpr"] == 0.80  # type: ignore[index]
    assert thresholds["minimum_worst_mixed_tpr"] == 0.60  # type: ignore[index]
    assert thresholds["maximum_complete_null_family_rejections"] == 4  # type: ignore[index]
    assert thresholds["complete_null_replicates"] == 40  # type: ignore[index]
    assert thresholds["bh_absolute_tolerance"] == 1e-12  # type: ignore[index]
    assert thresholds["selection_context"] == (  # type: ignore[index]
        "frozen_after_disclosed_exploratory_design_pilot"
    )
    assert thresholds["mroast_production_seed"] == 1729  # type: ignore[index]

    metrics = report["metrics"]
    assert metrics["scope"] == "backend_kernel_pathway"  # type: ignore[index]
    assert metrics["mixed_replicate_count"] == 20  # type: ignore[index]
    assert metrics["complete_null_replicate_count"] == 40  # type: ignore[index]
    assert len(metrics["mixed_replicates"]) == 20  # type: ignore[index]
    assert len(metrics["complete_null_replicates"]) == 40  # type: ignore[index]
    assert "gross calibration/power regressions" in metrics["gross_regression_limitation"]  # type: ignore[index]
    for method in ("limma_mroast", "limma_fry"):
        method_metrics = metrics["method_metrics"][method]  # type: ignore[index]
        assert method_metrics["mean_mixed_fdp"] <= 0.10
        assert method_metrics["worst_mixed_fdp"] <= 0.25
        assert method_metrics["mean_mixed_tpr"] >= 0.80
        assert method_metrics["worst_mixed_tpr"] >= 0.60
        assert method_metrics["wrong_direction_significant_positives"] == 0
        assert method_metrics["complete_null_family_rejections"] <= 4

    assertions = {item["name"]: item for item in report["assertions"]}  # type: ignore[union-attr]
    assert all(item["passed"] for item in assertions.values())
    assert "exact_replicate_and_set_counts" in assertions
    assert "python_bh_parity" in assertions
    for method in ("limma_mroast", "limma_fry"):
        for suffix in (
            "mean_mixed_fdp", "worst_mixed_fdp", "mean_mixed_tpr",
            "worst_mixed_tpr", "wrong_direction", "complete_null_family_fdr",
        ):
            assert f"{method}_{suffix}_gate" in assertions

    recorded = {item["name"]: item for item in report["implementation"]}  # type: ignore[union-attr]
    assert set(recorded) == set(IMPLEMENTATION_PATHS)
    for name, path in IMPLEMENTATION_PATHS.items():
        payload = path.read_bytes()
        assert recorded[name]["sha256"] == hashlib.sha256(payload).hexdigest()
        assert recorded[name]["size_bytes"] == len(payload)
    assert len(report["inputs"]) == 1 + 60 * 7  # type: ignore[arg-type]
    assert len(report["artifacts"]) == 60  # type: ignore[arg-type]
    jobs = [
        *(("mixed", replicate) for replicate in range(1, 21)),
        *(("complete_null", replicate) for replicate in range(1, 41)),
    ]
    fixture_members = (
        "counts.tsv", "metadata.tsv", "truth.tsv", "fixture.json",
        "sets.gmt", "annotation.tsv", "membership.tsv",
    )
    expected_input_names = [
        "compcodeR/DESCRIPTION",
        *(
            f"{scenario}-replicate-{replicate:02d}/{name}"
            for scenario, replicate in jobs
            for name in fixture_members
        ),
    ]
    expected_artifact_names = [
        f"{scenario}-replicate-{replicate:02d}/pathway_results.tsv"
        for scenario, replicate in jobs
    ]
    assert [item["name"] for item in report["inputs"]] == expected_input_names  # type: ignore[union-attr]
    assert [item["name"] for item in report["artifacts"]] == expected_artifact_names  # type: ignore[union-attr]
    evidence_payload = json.dumps(
        {"inputs": report["inputs"], "artifacts": report["artifacts"]},
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    assert hashlib.sha256(evidence_payload).hexdigest() == EVIDENCE_INVENTORY_SHA256
    serialized = ARCHIVED_REPORT.read_text(encoding="utf-8")
    assert "/tmp/" not in serialized
    assert "/home/" not in serialized


def test_complete_null_score_is_family_indicator_not_set_proportion(
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    _write_scoring_fixture(
        output, scenario="complete_null", significant_null_sets=17
    )

    score = _load_runner()._score_complete_null(1, output)

    assert score["family_rejections"] == {"limma_mroast": 1, "limma_fry": 1}


def test_scoring_rejects_a_duplicate_and_missing_non_scored_grid_cell(
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    _write_scoring_fixture(output, scenario="complete_null")
    path = output / "pathway_results.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    fields = [line.split("\t") for line in lines]
    fry_mixed = next(
        index
        for index, row in enumerate(fields[1:], start=1)
        if row[2:4] == ["limma_fry", "mixed"]
    )
    camera = next(
        index
        for index, row in enumerate(fields[1:], start=1)
        if row[1] == fields[fry_mixed][1]
        and row[2:4] == ["limma_camera", "directional"]
    )
    lines[fry_mixed] = lines[camera]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    runner = _load_runner()
    with pytest.raises(runner.BenchmarkError, match="exact canonical five-cell"):
        runner._read_directional_rows(output, "complete_null")


def test_mixed_score_uses_all_ten_pure_de_sets_as_tpr_denominator(
    tmp_path: Path,
) -> None:
    output = tmp_path / "output"
    _write_scoring_fixture(output, scenario="mixed", significant_null_sets=2)

    score = _load_runner()._score_mixed(1, output)

    for method in ("limma_mroast", "limma_fry"):
        metrics = score["methods"][method]
        assert metrics["true_positives"] == 10
        assert metrics["true_positive_rate"] == 1.0
        assert metrics["false_positives"] == 2
        assert metrics["false_discovery_proportion"] == pytest.approx(2 / 12)
        assert metrics["wrong_direction_significant_positives"] == 0


def test_aggregate_turns_each_bad_regression_metric_into_a_failed_gate() -> None:
    runner = _load_runner()
    bad_method = {
        "false_discovery_proportion": 0.5,
        "true_positive_rate": 0.1,
        "wrong_direction_significant_positives": 1,
    }
    mixed = [
        {
            "methods": {
                "limma_mroast": dict(bad_method),
                "limma_fry": dict(bad_method),
            },
            "bh_parity": {
                "violations": 0,
                "maximum_absolute_difference": 0.0,
                "family_ids_valid": True,
            },
        }
        for _ in range(20)
    ]
    complete_null = [
        {
            "family_rejections": {"limma_mroast": 1, "limma_fry": 1},
            "bh_parity": {
                "violations": 0,
                "maximum_absolute_difference": 0.0,
                "family_ids_valid": True,
            },
        }
        for _ in range(40)
    ]

    assertions, _ = runner._aggregate(mixed, complete_null)
    by_name = {item["name"]: item for item in assertions}

    for method in ("limma_mroast", "limma_fry"):
        for suffix in (
            "mean_mixed_fdp", "worst_mixed_fdp", "mean_mixed_tpr",
            "worst_mixed_tpr", "wrong_direction", "complete_null_family_fdr",
        ):
            assert not by_name[f"{method}_{suffix}_gate"]["passed"]


@pytest.mark.simulation
def test_pathway_compcoder_live_gate(tmp_path: Path) -> None:
    r_library_value = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not r_library_value:
        if os.environ.get("RNASEQ_P0_REQUIRE_BENCHMARKS") == "1":
            pytest.fail(
                "Certification mode requires RNASEQ_P0_R_LIBRARY; skipping is forbidden"
            )
        pytest.skip("RNASEQ_P0_R_LIBRARY is required for the locked pathway simulation gate")
    r_library = Path(r_library_value)
    rscript = Path(sys.executable).with_name("Rscript")
    if not rscript.is_file() or not r_library.is_dir():
        pytest.fail("The declared locked P0 R runtime is unavailable")

    report_directory_value = os.environ.get("RNASEQ_P0_BENCHMARK_REPORT_DIR")
    report_directory = Path(report_directory_value) if report_directory_value else tmp_path
    if report_directory.is_symlink():
        pytest.fail("RNASEQ_P0_BENCHMARK_REPORT_DIR must not be a symbolic link")
    report_directory.mkdir(parents=True, exist_ok=True)
    if report_directory.is_symlink() or not report_directory.is_dir():
        pytest.fail("RNASEQ_P0_BENCHMARK_REPORT_DIR is not a directory")
    report_path = report_directory / "pathway-compcoder-benchmark-report.json"
    completed = subprocess.run(
        [
            sys.executable, str(RUNNER), "--rscript", str(rscript),
            "--r-library", str(r_library), "--report", str(report_path),
        ],
        cwd=PROJECT_ROOT,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=7200,
    )
    assert report_path.is_file(), (completed.stdout, completed.stderr)
    report = _strict_json(report_path)
    assert completed.returncode == 0, report
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assert all(item["passed"] for item in report["assertions"])  # type: ignore[union-attr]
