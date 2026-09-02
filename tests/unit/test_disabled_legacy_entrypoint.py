"""The duplicated unsafe legacy entry point must remain fail-closed."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.unit
def test_legacy_module_help_resolves_without_optional_runtime() -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "legacy", "--help"],
        cwd=PROJECT_ROOT,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=10,
    )

    assert completed.returncode == 0
    assert "RNA-seq Downstream Analysis Pipeline" in completed.stdout
    assert completed.stderr == ""


@pytest.mark.unit
def test_integrated_legacy_script_refuses_to_bypass_checkpoint_a(
    tmp_path: Path,
) -> None:
    script = PROJECT_ROOT / "legacy" / "run_integrated_pydeseq2.py"

    completed = subprocess.run(
        [sys.executable, str(script)],
        cwd=tmp_path,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=10,
    )

    assert completed.returncode == 2
    assert completed.stdout == ""
    assert "disabled" in completed.stderr
    assert "Neither path is P0-certified" in completed.stderr
    assert list(tmp_path.iterdir()) == []


@pytest.mark.unit
def test_integrated_legacy_script_contains_no_unsafe_duplicate_analysis() -> None:
    source = (
        PROJECT_ROOT / "legacy" / "run_integrated_pydeseq2.py"
    ).read_text(encoding="utf-8")

    for prohibited in (
        "StandardScaler",
        ".round().astype(int)",
        "groupby(idx_col)",
        "common_samples",
        "DeseqStats",
    ):
        assert prohibited not in source
