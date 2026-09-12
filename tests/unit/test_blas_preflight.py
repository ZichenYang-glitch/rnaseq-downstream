"""Certification preflight tests; no R process or BLAS library is loaded."""

from __future__ import annotations

import ast
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from scripts.ci import blas_isa, preflight_blas as preflight
from scripts.ci import run_blas_diagnostic as diagnostic


pytestmark = pytest.mark.unit
REQUIRED_FLAGS = frozenset(
    {"avx2", "fma", "avx512f", "avx512dq", "avx512bw", "avx512vl", "avx512cd"}
)
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _mock_lscpu(monkeypatch, output: str):
    calls = []

    def run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(stdout=output, returncode=0)

    monkeypatch.setattr(preflight.subprocess, "run", run)
    return calls


def test_preflight_and_diagnostic_share_requirements_and_parser() -> None:
    assert preflight.CERTIFICATION_CORE == "SkylakeX"
    assert preflight.CORE_REQUIREMENTS is diagnostic.CORE_REQUIREMENTS
    assert preflight.CORE_REQUIREMENTS is blas_isa.CORE_REQUIREMENTS
    assert preflight.cpu_flags is diagnostic.cpu_flags is blas_isa.cpu_flags
    assert preflight.CORE_REQUIREMENTS["SkylakeX"] == REQUIRED_FLAGS


@pytest.mark.parametrize("extra", [set(), {"sse", "ssse3", "avx512_bf16"}])
def test_required_subset_and_superset_are_admitted(extra: set[str]) -> None:
    assert preflight.require_certification_isa(set(REQUIRED_FLAGS) | extra) is None


@pytest.mark.parametrize("missing", sorted(REQUIRED_FLAGS))
def test_each_missing_flag_is_a_blocking_failure(missing: str) -> None:
    with pytest.raises(ValueError) as error:
        preflight.require_certification_isa(set(REQUIRED_FLAGS) - {missing})
    message = str(error.value)
    assert "runner does not expose the required AVX-512" in message
    assert f"missing CPU flags: {missing}" in message
    assert "rerun to obtain a compatible machine" in message
    assert not any(
        word in message for word in ("--skip", "override", "force", "warning")
    )


@pytest.mark.parametrize(
    "output",
    [
        "",
        "Model name: Example\n",
        "Flags:\n",
        "Flags:   \t\n",
        "Flags: avx2,fma\n",
        "Flags: avx2 fma;exit\n",
        "Flags: avx2\nFlags: fma\n",
        "Flags: avx2\nFlags: avx2\n",
        "Flags: avx2 FMA\n",
    ],
)
def test_absent_empty_malformed_or_duplicate_flags_are_rejected(output: str) -> None:
    with pytest.raises(ValueError, match="lscpu"):
        blas_isa.cpu_flags(output)


def test_success_uses_c_locale_and_records_logs_and_summary_without_mutating_env(
    monkeypatch, tmp_path, capsys
) -> None:
    summary = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    monkeypatch.setenv("LC_ALL", "zh_CN.UTF-8")
    for key, value in {
        "OPENBLAS_CORETYPE": "SkylakeX",
        "OPENBLAS_NUM_THREADS": "7",
        "OPENBLAS_DEFAULT_NUM_THREADS": "8",
        "GOTO_NUM_THREADS": "9",
        "OMP_NUM_THREADS": "10",
    }.items():
        monkeypatch.setenv(key, value)
    original = os.environ.copy()
    calls = _mock_lscpu(monkeypatch, "Flags: " + " ".join(sorted(REQUIRED_FLAGS)))
    assert preflight.main([]) == 0
    assert os.environ == original
    assert len(calls) == 1
    command, kwargs = calls[0]
    assert command == ["lscpu"]
    assert kwargs["env"] == {**original, "LC_ALL": "C"}
    assert kwargs["check"] is True
    assert kwargs["timeout"] == 30
    assert kwargs["stdin"] == subprocess.DEVNULL
    captured = capsys.readouterr()
    assert not captured.err
    for message in (captured.out, summary.read_text()):
        assert "PASS: runner exposes the required SkylakeX ISA flags" in message
        assert all(flag in message for flag in REQUIRED_FLAGS)


def test_failure_is_logged_and_summarized_with_missing_flags(
    monkeypatch, tmp_path, capsys
) -> None:
    summary = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))
    _mock_lscpu(monkeypatch, "Flags: avx2 fma\n")
    assert preflight.main([]) == 1
    captured = capsys.readouterr()
    assert not captured.out
    for message in (captured.err, summary.read_text()):
        assert "ERROR: Certification BLAS ISA preflight failed" in message
        assert "runner does not expose the required AVX-512" in message
        assert "rerun to obtain a compatible machine" in message
        assert (
            "missing CPU flags: avx512bw, avx512cd, avx512dq, avx512f, avx512vl"
            in message
        )


@pytest.mark.parametrize(
    "error",
    [
        FileNotFoundError("lscpu unavailable"),
        subprocess.CalledProcessError(1, ["lscpu"]),
        subprocess.TimeoutExpired(["lscpu"], 30),
    ],
)
def test_unreadable_lscpu_fails_closed(monkeypatch, tmp_path, capsys, error) -> None:
    summary = tmp_path / "summary.md"
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(summary))

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(preflight.subprocess, "run", fail)
    assert preflight.main([]) == 1
    for message in (capsys.readouterr().err, summary.read_text()):
        assert "ERROR: Certification BLAS ISA preflight failed" in message
        assert "PASS" not in message


def test_malformed_lscpu_cannot_pass_main(monkeypatch, tmp_path, capsys) -> None:
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(tmp_path / "summary.md"))
    _mock_lscpu(monkeypatch, "Flags:\n")
    assert preflight.main([]) == 1
    assert "nonempty CPU Flags field" in capsys.readouterr().err


@pytest.mark.parametrize("core", ["automatic", "Haswell", "SkylakeX", ""])
def test_environment_cannot_bypass_isa_check(monkeypatch, tmp_path, core) -> None:
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(tmp_path / "summary.md"))
    monkeypatch.setenv("OPENBLAS_CORETYPE", core)
    monkeypatch.setenv("SKIP_BLAS_PREFLIGHT", "1")
    monkeypatch.setenv("BLAS_PREFLIGHT_SKIP", "true")
    monkeypatch.setenv("CPU_FLAGS", " ".join(REQUIRED_FLAGS))
    original = os.environ.copy()
    calls = _mock_lscpu(monkeypatch, "Flags: avx2 fma\n")
    assert preflight.main([]) == 1
    assert os.environ == original
    assert [call[0] for call in calls] == [["lscpu"]]


@pytest.mark.parametrize(
    "arguments",
    [
        ["--skip"],
        ["--skip-preflight"],
        ["--force"],
        ["--core", "Haswell"],
        ["--flags", " ".join(sorted(REQUIRED_FLAGS))],
        ["--allow-unsupported"],
    ],
)
def test_cli_rejects_bypass_and_alternate_core_options(monkeypatch, arguments) -> None:
    calls = _mock_lscpu(monkeypatch, "Flags: " + " ".join(REQUIRED_FLAGS))
    with pytest.raises(SystemExit) as error:
        preflight.main(arguments)
    assert error.value.code == 2
    assert not calls


def test_unwritable_summary_is_a_blocking_failure(
    monkeypatch, tmp_path, capsys
) -> None:
    monkeypatch.setenv("GITHUB_STEP_SUMMARY", str(tmp_path))
    _mock_lscpu(monkeypatch, "Flags: " + " ".join(REQUIRED_FLAGS))
    assert preflight.main([]) == 1
    assert "Cannot write the BLAS preflight step summary" in capsys.readouterr().err


def test_preflight_is_standalone_without_site_packages_or_blas_imports(
    tmp_path,
) -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-S",
            str(PROJECT_ROOT / "scripts/ci/preflight_blas.py"),
            "--help",
        ],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert "SkylakeX ISA flags" in result.stdout
    allowed_imports = {
        "__future__",
        "argparse",
        "os",
        "pathlib",
        "subprocess",
        "sys",
        "blas_isa",
        "re",
    }
    for module in (preflight, blas_isa):
        source = Path(module.__file__).read_text()
        imports = set()
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.Import):
                imports.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.add(node.module)
        assert imports <= allowed_imports
