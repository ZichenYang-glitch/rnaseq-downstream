"""Read-only CI diagnostics must expose dispatch and reject artifact drift."""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path

import pytest

from scripts.ci import check_benchmark_artifacts as artifacts
from scripts.ci import report_blas_runtime as runtime


pytestmark = pytest.mark.unit


def report() -> dict:
    return {
        "benchmark_id": "airway-example",
        "status": "pass",
        "artifacts": [
            {"name": "toolkit/results.tsv", "sha256": "a" * 64, "size_bytes": 120},
            {"name": "toolkit/design.tsv", "sha256": "b" * 64, "size_bytes": 10},
        ],
    }


def write_report(path: Path, value: object) -> Path:
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_cpu_identity_reports_model_and_relevant_flags() -> None:
    lines = runtime.cpu_identity(
        "Architecture: x86_64\nModel name: Example CPU\n"
        "Flags: sse2 avx2 avx512f avx512vl fma avx2\n"
    )
    assert lines == [
        "CPU Model name: Example CPU",
        "CPU BLAS ISA flags: avx2 avx512f avx512vl fma",
    ]
    assert "unavailable" in runtime.cpu_identity("")[0]


def test_symlink_chain_reports_relative_hops_missing_and_cycles(tmp_path: Path) -> None:
    library = tmp_path / "libopenblas.so"
    library.write_bytes(b"fixture")
    blas = tmp_path / "libblas.so.3"
    blas.symlink_to("libopenblas.so")
    assert runtime.link_chain(blas) == f"{blas} -> {library} -> [file]"
    missing = tmp_path / "libRblas.so"
    assert runtime.link_chain(missing).endswith("[missing]")
    missing.symlink_to("libRblas.so")
    assert runtime.link_chain(missing).endswith("[symlink cycle]")


def test_openblas_uses_loaded_c_api(monkeypatch: pytest.MonkeyPatch) -> None:
    class Library:
        pass

    library = Library()
    library.openblas_get_config = lambda: b"OpenBLAS DYNAMIC_ARCH NO_AFFINITY"
    library.openblas_get_corename = lambda: b"Haswell"
    library.openblas_get_num_threads = lambda: 7
    loaded = []

    def load(path):
        loaded.append(path)
        return library

    monkeypatch.setattr(runtime.ctypes, "CDLL", load)
    assert runtime.openblas_identity(Path("/locked/lib/libblas.so.3")) == (
        "OpenBLAS DYNAMIC_ARCH NO_AFFINITY",
        "Haswell",
        7,
    )
    assert loaded == ["/locked/lib/libblas.so.3"]


@pytest.mark.parametrize(
    ("config", "resolved", "requested", "expected", "message"),
    [
        ("OpenBLAS", "Haswell", "Haswell", None, "DYNAMIC_ARCH"),
        ("DYNAMIC_ARCH", "SkylakeX", "Haswell", "Haswell", "resolved core"),
        ("DYNAMIC_ARCH", "Haswell", "", "Haswell", "OPENBLAS_CORETYPE"),
    ],
)
def test_dispatch_validation_fails_closed(
    monkeypatch, config, resolved, requested, expected, message
) -> None:
    monkeypatch.setenv("OPENBLAS_CORETYPE", requested)
    with pytest.raises(ValueError, match=message):
        runtime.validate_identity(config, resolved, expected)


def test_runtime_logs_to_summary_without_mutating_thread_environment(
    tmp_path, monkeypatch, capsys
) -> None:
    summary = tmp_path / "summary.md"
    summary.write_text("Previous step\n", encoding="utf-8")
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    monkeypatch.setenv("OPENBLAS_CORETYPE", "Haswell")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "7")
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    monkeypatch.setattr(runtime, "command_log", lambda command: "Model name: Mock CPU")
    monkeypatch.setattr(
        runtime, "openblas_identity", lambda path: ("DYNAMIC_ARCH", "Haswell", 7)
    )
    before = dict(os.environ)
    assert runtime.main(["--prefix", str(tmp_path), "--expect-core", "Haswell"]) == 0
    assert dict(os.environ) == before
    stdout = capsys.readouterr().out
    assert "OpenBLAS resolved core: Haswell" in stdout
    assert "OPENBLAS_NUM_THREADS: 7" in stdout
    assert summary.read_text().startswith("Previous step\n")
    assert stdout.strip() in summary.read_text()


def test_runtime_failure_is_reported_in_log_and_summary(tmp_path, monkeypatch, capsys):
    summary = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    monkeypatch.setattr(runtime, "command_log", lambda command: "Model name: Mock CPU")
    monkeypatch.setattr(runtime, "openblas_identity", lambda path: ("static", "CPU", 1))
    assert runtime.main(["--prefix", str(tmp_path)]) == 1
    assert "ERROR:" in capsys.readouterr().out
    assert "DYNAMIC_ARCH" in summary.read_text()


def test_pre_restore_cpu_probe_does_not_load_a_library(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(tmp_path / "summary.md"))
    monkeypatch.setattr(runtime, "command_log", lambda command: "Model name: Mock CPU")

    def unexpected_load(path):
        pytest.fail("CPU-only probe must not load BLAS")

    monkeypatch.setattr(runtime, "openblas_identity", unexpected_load)
    assert runtime.main([]) == 0
    assert "CPU Model name: Mock CPU" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        runtime.main(["--expect-core", "Haswell"])


def test_exact_inventory_ignores_record_order_and_does_not_modify_reports(tmp_path):
    original = report()
    reordered = copy.deepcopy(original)
    reordered["artifacts"].reverse()
    baseline = write_report(tmp_path / "baseline.json", original)
    current = write_report(tmp_path / "current.json", reordered)
    before = baseline.read_bytes(), current.read_bytes()
    assert artifacts.compare_inventories(baseline, current)[0]
    assert before == (baseline.read_bytes(), current.read_bytes())


@pytest.mark.parametrize("change", ["sha256", "size", "missing", "extra"])
def test_inventory_rejects_changed_missing_or_extra_artifacts(tmp_path, change):
    original = report()
    changed = copy.deepcopy(original)
    if change == "sha256":
        changed["artifacts"][0]["sha256"] = "c" * 64
    elif change == "size":
        changed["artifacts"][0]["size_bytes"] += 1
    elif change == "missing":
        changed["artifacts"].pop()
    else:
        changed["artifacts"].append(
            {"name": "extra.tsv", "sha256": "d" * 64, "size_bytes": 1}
        )
    passed, lines = artifacts.compare_inventories(
        write_report(tmp_path / "baseline.json", original),
        write_report(tmp_path / "current.json", changed),
    )
    assert not passed
    assert lines[-1] == "Artifact byte identity: FAIL"


@pytest.mark.parametrize(
    "change", ["status", "id", "empty", "duplicate", "digest", "size", "shape"]
)
def test_inventory_rejects_malformed_or_wrong_benchmark_reports(tmp_path, change):
    original = report()
    changed = copy.deepcopy(original)
    if change == "status":
        changed["status"] = "fail"
    elif change == "id":
        changed["benchmark_id"] = "other"
    elif change == "empty":
        changed["artifacts"] = []
    elif change == "duplicate":
        changed["artifacts"].append(changed["artifacts"][0])
    elif change == "digest":
        changed["artifacts"][0]["sha256"] = "bad"
    elif change == "size":
        changed["artifacts"][0]["size_bytes"] = True
    else:
        changed["artifacts"][0]["new_field"] = 1
    with pytest.raises(ValueError):
        artifacts.compare_inventories(
            write_report(tmp_path / "baseline.json", original),
            write_report(tmp_path / "current.json", changed),
        )


def test_inventory_cli_rejects_invalid_json_and_summarizes(
    tmp_path, monkeypatch, capsys
):
    baseline = write_report(tmp_path / "baseline.json", report())
    current = tmp_path / "current.json"
    current.write_text('{"status":"pass","status":"fail"}', encoding="utf-8")
    summary = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    assert artifacts.main(["--baseline", str(baseline), "--current", str(current)]) == 2
    assert "Duplicate JSON object key" in capsys.readouterr().out
    assert "Duplicate JSON object key" in summary.read_text()


def test_committed_airway_inventory_comparisons_are_read_only():
    root = Path(__file__).resolve().parents[2]
    for name in ("airway-benchmark-report.json", "deseq2-airway-benchmark-report.json"):
        path = root / "tests/oracle" / name
        before = path.read_bytes()
        passed, _ = artifacts.compare_inventories(path, path)
        assert passed
        assert path.read_bytes() == before
