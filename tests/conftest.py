"""Shared pytest configuration for the P0 evidence lanes."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CLI_TIMEOUT_SECONDS = 10


def pytest_configure(config: pytest.Config) -> None:
    """Register the evidence-lane markers without relying on local plugins."""

    markers = {
        "unit": "contract, input-semantics, data-integrity, and QC-math tests",
        "integration": "cross-process CLI and evidence-publication tests",
        "oracle": "locked same-engine airway oracle parity tests",
        "simulation": "locked compcodeR FDR/TPR simulation gates",
    }
    for name, description in markers.items():
        config.addinivalue_line("markers", f"{name}: {description}")


class CliResult:
    """Small assertion-friendly wrapper around a CLI subprocess result."""

    def __init__(self, completed: subprocess.CompletedProcess[str]) -> None:
        self.returncode = completed.returncode
        self.stdout = completed.stdout
        self.stderr = completed.stderr

    def json(self) -> dict[str, Any]:
        """Parse stdout as exactly one JSON document."""

        assert self.stdout.endswith("\n"), "CLI stdout must end with one newline"
        assert self.stdout.count("\n") == 1, (
            "CLI stdout must contain exactly one newline-terminated JSON line"
        )
        line = self.stdout[:-1]
        assert line, "CLI stdout must contain a JSON document"
        document = json.loads(line)
        assert isinstance(document, dict), "CLI response must be a JSON object"
        return document


@pytest.fixture
def run_module_cli():
    """Run the installed module entry point from any requested directory.

    Tests intentionally remove ``PYTHONPATH``. This ensures subprocesses use
    the installed distribution instead of succeeding only because the checkout
    happens to be their current working directory.
    """

    def _run(*arguments: str, cwd: Path | None = None) -> CliResult:
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        completed = subprocess.run(
            [sys.executable, "-m", "rnaseq_downstream", *arguments],
            cwd=str(cwd or PROJECT_ROOT),
            env=env,
            check=False,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=CLI_TIMEOUT_SECONDS,
        )
        return CliResult(completed)

    return _run
