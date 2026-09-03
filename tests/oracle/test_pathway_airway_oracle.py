"""Locked airway edgeR/limma pathway oracle parity gate."""

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
RUNNER = PROJECT_ROOT / "scripts/benchmark/run_airway_pathway_oracle.py"
ARCHIVED_REPORT = Path(__file__).with_name("pathway-airway-benchmark-report.json")
BENCHMARK_ID = "airway-limma-gene-set-same-engine-v1"
EVIDENCE_INVENTORY_SHA256 = (
    "2b60f13beb6e9dd06f6cb39f5639d1a376012208c26af2b9abf6b5e57fb16c62"
)
IMPLEMENTATION_PATHS = {
    "run_airway_pathway_oracle.py": RUNNER,
    "common.py": PROJECT_ROOT / "scripts/benchmark/common.py",
    "benchmark-report-v1.schema.json": (
        PROJECT_ROOT / "scripts/benchmark/benchmark-report-v1.schema.json"
    ),
    "prepare_airway_fixture.R": (
        PROJECT_ROOT / "scripts/benchmark/prepare_airway_fixture.R"
    ),
    "run_airway_pathway_direct_oracle.R": (
        PROJECT_ROOT / "scripts/benchmark/run_airway_pathway_direct_oracle.R"
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
REQUIRED_ASSERTIONS = {
    "benchmark_execution",
    "locked_runtime_identity",
    "exact_pathway_grid_status_and_order",
    "tested_method_shape_and_direction_parity",
    "mapping_and_method_metadata_parity",
    "p_value_numeric_parity",
    "fdr_numeric_parity",
    "proportion_down_numeric_parity",
    "proportion_up_numeric_parity",
    "correlation_estimate_raw_numeric_parity",
    "correlation_effective_numeric_parity",
    "vif_used_numeric_parity",
    "python_bh_parity",
}


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


def _load_runner() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "airway_pathway_oracle_runner", RUNNER
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.path.insert(0, str(RUNNER.parent))
    try:
        specification.loader.exec_module(module)
    finally:
        sys.path.remove(str(RUNNER.parent))
    return module


@pytest.mark.oracle
def test_archived_pathway_airway_report_is_passing_evidence() -> None:
    report = _strict_json(ARCHIVED_REPORT)
    runner = _load_runner()
    runner.assert_report_shape(report, benchmark_id=BENCHMARK_ID)
    assert report["schema_version"] == "rnaseq-downstream-benchmark-report-v1"
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assert report["metrics"]["scope"] == "backend_kernel_pathway"  # type: ignore[index]
    assertions = {item["name"]: item for item in report["assertions"]}  # type: ignore[union-attr]
    assert set(assertions) == REQUIRED_ASSERTIONS
    assert all(item["passed"] for item in assertions.values())
    assert report["thresholds"] == {  # type: ignore[comparison-overlap]
        "absolute_tolerance": 1e-10,
        "bh_absolute_tolerance": 1e-12,
        "minimum_tested_genes": 10,
        "mroast_nrot": 9999,
        "mroast_seed": 1729,
        "relative_tolerance": 1e-6,
    }
    recorded = {item["name"]: item for item in report["implementation"]}  # type: ignore[union-attr]
    assert set(recorded) == set(IMPLEMENTATION_PATHS)
    for name, path in IMPLEMENTATION_PATHS.items():
        payload = path.read_bytes()
        assert recorded[name]["sha256"] == hashlib.sha256(payload).hexdigest()
        assert recorded[name]["size_bytes"] == len(payload)
    inputs = {item["name"] for item in report["inputs"]}  # type: ignore[union-attr]
    assert inputs == {
        "airway/DESCRIPTION",
        "fixture/counts.tsv", "fixture/metadata.tsv", "fixture/fixture.json",
        "fixture/sets.gmt", "fixture/annotation.tsv", "fixture/gene-set-source.json",
    }
    artifacts = {item["name"] for item in report["artifacts"]}  # type: ignore[union-attr]
    assert artifacts == {
        "direct-oracle/pathway_oracle.tsv",
        "direct-oracle/diagnostics.json",
        "toolkit/pathway_results.tsv",
    }
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


def test_pathway_oracle_bh_implementation_handles_ties_and_unsorted_values() -> None:
    runner = _load_runner()
    observed = runner._bh([0.01, 0.04, 0.01, 0.9])
    assert observed == pytest.approx([0.02, 0.05333333333333334, 0.02, 0.9])


@pytest.mark.oracle
def test_airway_pathway_oracle_live_gate(tmp_path: Path) -> None:
    r_library_value = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not r_library_value:
        if os.environ.get("RNASEQ_P0_REQUIRE_BENCHMARKS") == "1":
            pytest.fail(
                "Certification mode requires RNASEQ_P0_R_LIBRARY; skipping is forbidden"
            )
        pytest.skip("RNASEQ_P0_R_LIBRARY is required for the locked pathway oracle")
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
    report_path = report_directory / "pathway-airway-benchmark-report.json"
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
        timeout=1800,
    )
    assert report_path.is_file(), (completed.stdout, completed.stderr)
    report = _strict_json(report_path)
    assert completed.returncode == 0, report
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assert all(item["passed"] for item in report["assertions"])  # type: ignore[union-attr]
