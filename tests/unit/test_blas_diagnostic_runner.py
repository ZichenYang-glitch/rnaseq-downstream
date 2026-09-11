"""Mocked diagnostic control-flow tests; no R process or BLAS library is loaded."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.ci import run_blas_diagnostic as diagnostic


@pytest.mark.unit
def test_cpu_flags_handles_linux_sse3_alias_and_requires_flags() -> None:
    assert diagnostic.cpu_flags("Model name: Example\nFlags: pni avx2 fma\n") == {
        "pni",
        "sse3",
        "avx2",
        "fma",
    }
    with pytest.raises(ValueError, match="explicit core probes are unsafe"):
        diagnostic.cpu_flags("Model name: Example\n")


@pytest.mark.unit
def test_candidate_environment_preserves_all_thread_controls(tmp_path: Path) -> None:
    original = {
        "OPENBLAS_CORETYPE": "Haswell",
        "OPENBLAS_VERBOSE": "9",
        "OPENBLAS_NUM_THREADS": "7",
        "OPENBLAS_DEFAULT_NUM_THREADS": "9",
        "OMP_NUM_THREADS": "3",
        "GOTO_NUM_THREADS": "4",
        "RNASEQ_P0_R_LIBRARY": "/locked/library",
    }
    automatic = diagnostic.candidate_environment(original, "automatic", tmp_path)
    explicit = diagnostic.candidate_environment(original, "Core2", tmp_path)
    assert "OPENBLAS_CORETYPE" not in automatic
    assert explicit["OPENBLAS_CORETYPE"] == "Core2"
    assert original["OPENBLAS_CORETYPE"] == "Haswell"
    for key, value in original.items():
        if key not in {"OPENBLAS_CORETYPE", "OPENBLAS_VERBOSE"}:
            assert automatic[key] == explicit[key] == value
    assert "OPENBLAS_VERBOSE" not in automatic
    assert "OPENBLAS_VERBOSE" not in explicit
    assert original["OPENBLAS_VERBOSE"] == "9"
    assert automatic["RNASEQ_P0_REQUIRE_BENCHMARKS"] == "1"
    assert automatic["RNASEQ_P0_BENCHMARK_REPORT_DIR"] == str(tmp_path)


def _mock_commands(monkeypatch, outcomes=None):
    calls = []
    outcomes = outcomes or {}

    def run(command, environment):
        core = environment.get("OPENBLAS_CORETYPE", "automatic")
        stage = (
            "probe"
            if "report_blas_runtime.py" in command[1]
            else "gate"
            if command[1] == "-m"
            else "compare"
        )
        calls.append((core, stage, command, environment.copy()))
        return outcomes.get((core, stage), 0)

    def r_probe(prefix, environment):
        core = environment.get("OPENBLAS_CORETYPE", "automatic")
        calls.append((core, "r_probe", [str(prefix)], environment.copy()))
        return outcomes.get((core, "r_probe"), 0)

    monkeypatch.setattr(diagnostic, "run_command", run)
    monkeypatch.setattr(diagnostic, "run_r_probe", r_probe)
    return calls


@pytest.mark.unit
def test_unsafe_candidates_are_skipped_before_probe_or_gate(
    monkeypatch, tmp_path
) -> None:
    calls = _mock_commands(monkeypatch)
    result = diagnostic.run_diagnostic(tmp_path, tmp_path / "results", {}, {"sse3"})
    assert result == 0
    assert {call[0] for call in calls} == {"automatic", "Prescott", "Barcelona"}
    assert [call[1] for call in calls[:4]] == ["probe", "r_probe", "gate", "compare"]


@pytest.mark.unit
def test_full_grid_runs_only_airway_with_unique_report_directories(
    monkeypatch,
    tmp_path,
) -> None:
    calls = _mock_commands(monkeypatch)
    flags = set().union(*diagnostic.CORE_REQUIREMENTS.values())
    result = diagnostic.run_diagnostic(
        tmp_path, tmp_path / "results", {"OPENBLAS_VERBOSE": "9"}, flags
    )
    assert result == 0
    gates = [call for call in calls if call[1] == "gate"]
    assert [call[0] for call in gates] == list(diagnostic.CORE_REQUIREMENTS)
    assert len({call[3]["RNASEQ_P0_BENCHMARK_REPORT_DIR"] for call in gates}) == 10
    for _, _, command, environment in gates:
        assert command == [
            str(tmp_path / "bin/python"),
            "-m",
            "pytest",
            "tests/oracle/test_airway_oracle.py",
            "-v",
        ]
        assert environment["RNASEQ_P0_REQUIRE_BENCHMARKS"] == "1"
    probes = [call for call in calls if call[1] == "probe"]
    assert "--expect-core" not in probes[0][2]
    for core, _, command, _ in probes[1:]:
        assert command[-2:] == ["--expect-core", core]
    for _, stage, _, environment in calls:
        if stage in {"probe", "r_probe"}:
            assert environment["OPENBLAS_VERBOSE"] == "2"
        else:
            assert "OPENBLAS_VERBOSE" not in environment


@pytest.mark.unit
def test_mismatch_does_not_prevent_later_candidate_matching(
    monkeypatch, tmp_path
) -> None:
    calls = _mock_commands(monkeypatch, {("automatic", "compare"): 1})
    assert diagnostic.run_diagnostic(tmp_path, tmp_path / "results", {}, {"sse3"}) == 0
    assert ("Prescott", "compare") in [(call[0], call[1]) for call in calls]


@pytest.mark.unit
def test_no_exact_core_fails(monkeypatch, tmp_path) -> None:
    _mock_commands(monkeypatch, {("automatic", "compare"): 1})
    assert diagnostic.run_diagnostic(tmp_path, tmp_path / "results", {}, set()) == 1


@pytest.mark.unit
@pytest.mark.parametrize("stage", ["probe", "r_probe", "gate", "compare"])
def test_execution_failure_is_fatal_even_when_another_core_matches(
    monkeypatch,
    tmp_path,
    stage,
) -> None:
    calls = _mock_commands(monkeypatch, {("automatic", stage): 2})
    assert diagnostic.run_diagnostic(tmp_path, tmp_path / "results", {}, {"sse3"}) == 1
    automatic_stages = [call[1] for call in calls if call[0] == "automatic"]
    assert automatic_stages[-1] == stage
    assert ("Prescott", "compare") in [(call[0], call[1]) for call in calls]


@pytest.mark.unit
def test_candidate_outcomes_are_logged_and_summarized(
    monkeypatch, tmp_path, capsys
) -> None:
    _mock_commands(monkeypatch, {("automatic", "compare"): 1})
    summary = tmp_path / "summary.md"
    diagnostic.run_diagnostic(
        tmp_path,
        tmp_path / "results",
        {"GITHUB_STEP_SUMMARY": str(summary)},
        {"sse3"},
    )
    for text in (capsys.readouterr().out, summary.read_text()):
        assert "automatic: nonmatch" in text
        assert "Prescott: match" in text
        assert "Haswell: skipped: missing ISA flags avx2, fma" in text


@pytest.mark.unit
def test_skylakex_and_cooperlake_require_complete_conservative_flags(
    monkeypatch,
    tmp_path,
) -> None:
    calls = _mock_commands(monkeypatch)
    flags = set(diagnostic.CORE_REQUIREMENTS["Cooperlake"]) - {"avx512cd"}
    diagnostic.run_diagnostic(tmp_path, tmp_path / "results", {}, flags)
    assert not {"SkylakeX", "Cooperlake"} & {call[0] for call in calls}


@pytest.mark.unit
def test_standalone_r_probe_forwards_both_output_channels_with_timeout(
    monkeypatch,
    tmp_path,
    capsys,
) -> None:
    calls = []

    def run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(
            stdout="Core: Haswell\nR BLAS: /locked/libblas.so\n",
            returncode=0,
        )

    monkeypatch.setattr(diagnostic.subprocess, "run", run)
    summary = tmp_path / "summary.md"
    environment = {"GITHUB_STEP_SUMMARY": str(summary), "OPENBLAS_VERBOSE": "2"}
    assert diagnostic.run_r_probe(tmp_path, environment) == 0
    command, kwargs = calls[0]
    assert command[:3] == [str(tmp_path / "bin/Rscript"), "--vanilla", "-e"]
    assert "crossprod" in command[3]
    assert "extSoftVersion" in command[3]
    assert kwargs["stderr"] == diagnostic.subprocess.STDOUT
    assert kwargs["timeout"] == 30
    for text in (capsys.readouterr().out, summary.read_text()):
        assert "Core: Haswell" in text
        assert "R BLAS: /locked/libblas.so" in text
