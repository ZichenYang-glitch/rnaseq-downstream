"""Locked compcodeR negative-binomial FDR/TPR gate."""

from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import subprocess
import sys

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUNNER = PROJECT_ROOT / "scripts/benchmark/run_compcoder_gate.py"
ARCHIVED_REPORT = Path(__file__).with_name("compcoder-benchmark-report.json")
BENCHMARK_ID = "compcoder-edger-ql-nb-fdr-tpr-v1"
IMPLEMENTATION_PATHS = {
    "run_compcoder_gate.py": RUNNER,
    "common.py": PROJECT_ROOT / "scripts/benchmark/common.py",
    "benchmark-report-v1.schema.json": (
        PROJECT_ROOT / "scripts/benchmark/benchmark-report-v1.schema.json"
    ),
    "generate_compcoder_fixture.R": (
        PROJECT_ROOT / "scripts/benchmark/generate_compcoder_fixture.R"
    ),
    "edger_backend.py": PROJECT_ROOT / "rnaseq_downstream/edger_backend.py",
    "analysis_contract.py": PROJECT_ROOT / "rnaseq_downstream/analysis_contract.py",
    "edger_ql.R": PROJECT_ROOT / "rnaseq_downstream/r_scripts/edger_ql.R",
    "conda-lock.yml": PROJECT_ROOT / "conda-lock.yml",
    "renv.lock": PROJECT_ROOT / "renv.lock",
    "environment.p0.yml": PROJECT_ROOT / "environment.p0.yml",
    "r-sources.lock": PROJECT_ROOT / "environment/r-sources.lock",
    "verify.R": PROJECT_ROOT / "environment/verify.R",
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


@pytest.mark.simulation
def test_archived_compcoder_report_is_passing_evidence() -> None:
    report = _strict_json(ARCHIVED_REPORT)
    assert report["schema_version"] == "rnaseq-downstream-benchmark-report-v1"
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assert report["metrics"]["scope"] == "backend_kernel"  # type: ignore[index]
    assert all(item["passed"] for item in report["assertions"])  # type: ignore[union-attr]
    recorded = {item["name"]: item for item in report["implementation"]}  # type: ignore[union-attr]
    assert set(recorded) == set(IMPLEMENTATION_PATHS)
    for name, path in IMPLEMENTATION_PATHS.items():
        payload = path.read_bytes()
        assert recorded[name]["sha256"] == hashlib.sha256(payload).hexdigest()
        assert recorded[name]["size_bytes"] == len(payload)
    serialized = ARCHIVED_REPORT.read_text(encoding="utf-8")
    assert "/tmp/" not in serialized
    assert "/home/" not in serialized


@pytest.mark.simulation
def test_compcoder_negative_binomial_fdr_tpr_gate(tmp_path: Path) -> None:
    r_library_value = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not r_library_value:
        if os.environ.get("RNASEQ_P0_REQUIRE_BENCHMARKS") == "1":
            pytest.fail(
                "Certification mode requires RNASEQ_P0_R_LIBRARY; skipping is forbidden"
            )
        pytest.skip("RNASEQ_P0_R_LIBRARY is required for the locked simulation gate")
    r_library = Path(r_library_value)
    rscript = Path(sys.executable).with_name("Rscript")
    if not rscript.is_file() or not r_library.is_dir():
        pytest.fail("The declared locked P0 R runtime is unavailable")

    report_directory_value = os.environ.get("RNASEQ_P0_BENCHMARK_REPORT_DIR")
    report_directory = (
        Path(report_directory_value) if report_directory_value else tmp_path
    )
    if report_directory.is_symlink():
        pytest.fail("RNASEQ_P0_BENCHMARK_REPORT_DIR must not be a symbolic link")
    report_directory.mkdir(parents=True, exist_ok=True)
    if report_directory.is_symlink() or not report_directory.is_dir():
        pytest.fail("RNASEQ_P0_BENCHMARK_REPORT_DIR is not a directory")
    report_path = report_directory / "compcoder-benchmark-report.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(RUNNER),
            "--rscript",
            str(rscript),
            "--r-library",
            str(r_library),
            "--report",
            str(report_path),
        ],
        cwd=PROJECT_ROOT,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert report_path.is_file(), (completed.stdout, completed.stderr)
    report = _strict_json(report_path)
    assert completed.returncode == 0, report
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assert all(item["passed"] for item in report["assertions"])  # type: ignore[union-attr]
