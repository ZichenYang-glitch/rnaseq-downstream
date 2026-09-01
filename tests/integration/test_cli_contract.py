"""Cross-process tests for the agent-safe JSON CLI foundation."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CLI_TIMEOUT_SECONDS = 10
EXPECTED_KEYS = [
    "schema_version",
    "command",
    "status",
    "data",
    "warnings",
    "errors",
    "artifacts",
]


def assert_envelope_shape(document: dict, *, command: str, status: str) -> None:
    """Assert the stable, machine-oriented top-level response contract."""

    assert list(document) == EXPECTED_KEYS
    assert document["schema_version"] == "1.0"
    assert document["command"] == command
    assert document["status"] == status
    assert isinstance(document["warnings"], list)
    assert isinstance(document["errors"], list)
    assert isinstance(document["artifacts"], list)


@pytest.mark.integration
def test_capabilities_is_one_json_document_on_stdout(run_module_cli) -> None:
    result = run_module_cli("capabilities")

    assert result.returncode == 0
    document = result.json()
    assert_envelope_shape(document, command="capabilities", status="success")
    assert isinstance(document["data"], dict)
    assert document["errors"] == []
    assert document["data"]["toolkit"]["maturity"] == "research_preview"
    assert document["data"]["certified_analysis_paths"] == []
    assert "NaN" not in result.stdout
    assert "Infinity" not in result.stdout
    assert result.stderr == ""


@pytest.mark.integration
def test_module_entrypoint_does_not_depend_on_checkout_cwd(
    run_module_cli, tmp_path: Path
) -> None:
    result = run_module_cli("capabilities", cwd=tmp_path)

    assert result.returncode == 0
    assert_envelope_shape(
        result.json(), command="capabilities", status="success"
    )


@pytest.mark.integration
def test_console_entrypoint_is_installed_and_json_only(tmp_path: Path) -> None:
    executable = shutil.which("rnaseq-downstream")
    assert executable is not None, (
        "The console script is absent; install the project with "
        "`python -m pip install --no-deps -e .` before running integration tests"
    )
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)

    completed = subprocess.run(
        [executable, "capabilities"],
        cwd=tmp_path,
        env=env,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=CLI_TIMEOUT_SECONDS,
    )

    assert completed.returncode == 0
    assert completed.stdout.endswith("\n")
    assert completed.stdout.count("\n") == 1
    document = json.loads(completed.stdout[:-1])
    assert_envelope_shape(document, command="capabilities", status="success")
    assert completed.stderr == ""


@pytest.mark.integration
def test_editable_install_resolves_to_current_checkout(tmp_path: Path) -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    source = (
        "import json, pathlib, rnaseq_downstream; "
        "print(json.dumps({'origin': str(pathlib.Path("
        "rnaseq_downstream.__file__).resolve())}))"
    )

    completed = subprocess.run(
        [sys.executable, "-c", source],
        cwd=tmp_path,
        env=env,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=CLI_TIMEOUT_SECONDS,
    )

    assert completed.returncode == 0
    assert completed.stderr == ""
    assert completed.stdout.count("\n") == 1
    origin = Path(json.loads(completed.stdout)["origin"])
    assert origin == (PROJECT_ROOT / "rnaseq_downstream" / "__init__.py").resolve()


@pytest.mark.integration
def test_legacy_checkout_resolves_flat_error_package_without_site_packages() -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    source = (
        "import json, pathlib, sys, types; "
        "sys.modules.update({name: types.ModuleType(name) "
        "for name in ('pandas', 'yaml')}); "
        "from modules import data; "
        "import rnaseq_downstream.errors as errors; "
        "from rnaseq_downstream.errors import InputReadError; "
        "print(json.dumps({'same_type': data.InputReadError is InputReadError, "
        "'origin': str(pathlib.Path(errors.__file__).resolve())}))"
    )

    completed = subprocess.run(
        [sys.executable, "-S", "-c", source],
        cwd=PROJECT_ROOT,
        env=env,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=CLI_TIMEOUT_SECONDS,
    )

    assert completed.returncode == 0
    assert completed.stderr == ""
    assert completed.stdout.count("\n") == 1
    document = json.loads(completed.stdout)
    assert document["same_type"] is True
    assert Path(document["origin"]) == (
        PROJECT_ROOT / "rnaseq_downstream" / "errors.py"
    ).resolve()


@pytest.mark.integration
@pytest.mark.parametrize("command", ["inspect", "validate", "run", "summarize"])
def test_unimplemented_p0_stages_fail_explicitly(
    run_module_cli, command: str
) -> None:
    result = run_module_cli(command)

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command=command, status="error")
    assert document["errors"]
    assert document["errors"][0]["code"] == "FEATURE_NOT_IMPLEMENTED"
    assert result.stderr == ""


@pytest.mark.integration
def test_invalid_arguments_are_structured_json(run_module_cli) -> None:
    result = run_module_cli("not-a-command")

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command="not-a-command", status="error")
    assert document["errors"][0]["code"] == "INVALID_REQUEST"
    assert document["errors"][0]["details"]["usage"].startswith("usage:")
    assert result.stderr == ""


@pytest.mark.integration
def test_missing_command_returns_usage_inside_json(run_module_cli) -> None:
    result = run_module_cli()

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command="cli", status="error")
    assert document["errors"][0]["code"] == "INVALID_REQUEST"
    assert "usage" in document["errors"][0]["details"]
    assert result.stderr == ""


@pytest.mark.integration
def test_help_is_machine_readable_json(run_module_cli) -> None:
    result = run_module_cli("--help")

    assert result.returncode == 0
    document = result.json()
    assert_envelope_shape(document, command="cli", status="success")
    assert isinstance(document["data"].get("help"), str)
    assert document["data"]["help"]
    assert "usage:" not in result.stderr.lower()
