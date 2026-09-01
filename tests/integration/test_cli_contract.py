"""Cross-process tests for the agent-safe JSON CLI foundation."""

from __future__ import annotations

import hashlib
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


def _write_json(path: Path, document: object) -> None:
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _featurecounts_request(
    root: Path,
    *,
    invalid_count: str | None = None,
    metadata_order: tuple[str, str] = ("s1", "s2"),
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    metadata = root / "metadata.tsv"
    metadata_rows = ["sample_id\tcondition"]
    metadata_rows.extend(
        f"{sample_id}\t{'control' if sample_id == 's1' else 'treated'}"
        for sample_id in metadata_order
    )
    metadata.write_text("\n".join(metadata_rows) + "\n", encoding="utf-8")
    reference = root / "reference.fa"
    reference.write_text(">gene1\nACGT\n", encoding="utf-8")
    matrix = root / "counts.tsv"
    first = invalid_count if invalid_count is not None else "10"
    matrix.write_text(
        f"gene_id\tgene_name\ts1\ts2\nENSG1\tGENE1\t{first}\t20\n",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": "1.0",
        "artifact_type": "featurecounts_integer_matrix",
        "artifact": {"path": "counts.tsv", "sha256": _sha256(matrix)},
        "gene_id_column": "gene_id",
        "display_columns": ["gene_name"],
        "sample_columns": ["s1", "s2"],
        "producer": {"name": "featureCounts", "version": "2.0.6"},
        "reference": {
            "name": "synthetic",
            "version": "v1",
            "source": "reference.fa",
            "sha256": _sha256(reference),
        },
    }
    _write_json(root / "counts.manifest.json", manifest)
    request = {
        "schema_version": "1.0",
        "input_semantics": "featurecounts_integer",
        "metadata": {"path": "metadata.tsv", "sample_id_column": "sample_id"},
        "producer": {"name": "featureCounts", "version": "2.0.6"},
        "reference": {
            "name": "synthetic",
            "version": "v1",
            "source": "reference.fa",
            "sha256": _sha256(reference),
        },
        "gene_id": {"strip_version": False},
        "samples": [{"sample_id": "s1"}, {"sample_id": "s2"}],
        "featurecounts": {
            "layout": "combined_matrix",
            "matrix": "counts.tsv",
            "manifest": "counts.manifest.json",
        },
    }
    request_path = root / "request.json"
    _write_json(request_path, request)
    return request_path


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
    assert_envelope_shape(result.json(), command="capabilities", status="success")


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
    assert (
        Path(document["origin"])
        == (PROJECT_ROOT / "rnaseq_downstream" / "errors.py").resolve()
    )


@pytest.mark.integration
@pytest.mark.parametrize("command", ["run", "summarize"])
def test_unimplemented_p0_stages_fail_explicitly(run_module_cli, command: str) -> None:
    result = run_module_cli(command)

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command=command, status="error")
    assert document["errors"]
    assert document["errors"][0]["code"] == "FEATURE_NOT_IMPLEMENTED"
    assert result.stderr == ""


@pytest.mark.integration
def test_inspect_is_read_only_and_reports_input_only_scope(
    run_module_cli,
    tmp_path: Path,
) -> None:
    request = _featurecounts_request(tmp_path / "inputs")
    unrelated_cwd = tmp_path / "elsewhere"
    unrelated_cwd.mkdir()

    result = run_module_cli(
        "inspect",
        "--request",
        str(request),
        cwd=unrelated_cwd,
    )

    assert result.returncode == 0
    document = result.json()
    assert_envelope_shape(document, command="inspect", status="success")
    assert document["data"]["scope"]["validation_scope"] == "input_semantics_only"
    assert document["data"]["scope"]["input_semantics"] == "inspected"
    assert document["data"]["scope"]["full_numeric_validation"] == "not_run"
    assert document["data"]["scope"]["design"] == "not_run"
    assert document["data"]["scope"]["backend"] == "not_run"
    assert document["data"]["scope"]["runnable"] is False
    assert document["artifacts"]
    assert not (tmp_path / "evidence").exists()
    assert result.stderr == ""


@pytest.mark.integration
def test_inspect_never_claims_full_validation_for_invalid_counts(
    run_module_cli,
    tmp_path: Path,
) -> None:
    request = _featurecounts_request(tmp_path / "inputs", invalid_count="1.5")

    result = run_module_cli("inspect", "--request", str(request))

    assert result.returncode == 0
    document = result.json()
    assert_envelope_shape(document, command="inspect", status="success")
    assert document["data"]["scope"]["input_semantics"] == "inspected"
    assert document["data"]["scope"]["full_numeric_validation"] == "not_run"
    assert document["data"]["input"]["input_certification_eligible"] is None
    assert result.stderr == ""


@pytest.mark.integration
def test_structured_input_warning_propagates_to_cli_envelope(
    run_module_cli,
    tmp_path: Path,
) -> None:
    request = _featurecounts_request(
        tmp_path / "inputs",
        metadata_order=("s2", "s1"),
    )

    result = run_module_cli("inspect", "--request", str(request))

    assert result.returncode == 0
    document = result.json()
    assert [warning["code"] for warning in document["warnings"]] == [
        "METADATA_ORDER_NORMALIZED"
    ]
    assert document["warnings"][0]["severity"] == "warning"
    assert result.stderr == ""


@pytest.mark.integration
def test_validate_publishes_complete_noncertifying_bundle(
    run_module_cli,
    tmp_path: Path,
) -> None:
    request = _featurecounts_request(tmp_path / "inputs")
    output = tmp_path / "evidence"

    result = run_module_cli(
        "validate",
        "--request",
        str(request),
        "--output-dir",
        str(output),
        cwd=tmp_path,
    )

    assert result.returncode == 0
    document = result.json()
    assert_envelope_shape(document, command="validate", status="success")
    assert document["data"]["scope"]["analysis_path_certified"] is False
    assert document["data"]["scope"]["runnable"] is False
    assert document["data"]["bundle"]["validation"]["design"] == "not_run"
    assert document["data"]["bundle"]["publication_status"] == ("durability_confirmed")
    assert document["data"]["bundle"]["warnings"] == []
    assert sorted(path.name for path in output.iterdir()) == [
        "bundle_manifest.json",
        "input_plan.json",
        "provenance.json",
        "validated_request.json",
    ]
    generated = [
        artifact
        for artifact in document["artifacts"]
        if artifact["kind"] == "generated_validation_artifact"
    ]
    assert len(generated) == 4
    assert all(Path(artifact["path"]).is_file() for artifact in generated)
    assert result.stderr == ""


@pytest.mark.integration
def test_validate_never_overwrites_existing_evidence(
    run_module_cli,
    tmp_path: Path,
) -> None:
    request = _featurecounts_request(tmp_path / "inputs")
    output = tmp_path / "evidence"
    output.mkdir()
    marker = output / "archived.txt"
    marker.write_text("keep", encoding="utf-8")

    result = run_module_cli(
        "validate",
        "--request",
        str(request),
        "--output-dir",
        str(output),
    )

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command="validate", status="error")
    assert document["errors"][0]["code"] == "INVALID_REQUEST"
    assert marker.read_text(encoding="utf-8") == "keep"
    assert result.stderr == ""


@pytest.mark.integration
def test_validation_failure_creates_no_bundle(
    run_module_cli,
    tmp_path: Path,
) -> None:
    request = _featurecounts_request(tmp_path / "inputs", invalid_count="1.5")
    output = tmp_path / "evidence"

    result = run_module_cli(
        "validate",
        "--request",
        str(request),
        "--output-dir",
        str(output),
    )

    assert result.returncode == 3
    document = result.json()
    assert_envelope_shape(document, command="validate", status="error")
    assert document["errors"][0]["code"] == "COUNT_VALUES_INVALID"
    assert not output.exists()
    assert result.stderr == ""


@pytest.mark.integration
def test_validation_output_failure_uses_stable_json_error(
    run_module_cli,
    tmp_path: Path,
) -> None:
    request = _featurecounts_request(tmp_path / "inputs")
    not_a_directory = tmp_path / "occupied"
    not_a_directory.write_text("file", encoding="utf-8")
    output = not_a_directory / "evidence"

    result = run_module_cli(
        "validate",
        "--request",
        str(request),
        "--output-dir",
        str(output),
    )

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command="validate", status="error")
    assert document["errors"][0]["code"] == "OUTPUT_WRITE_FAILED"
    assert document["errors"][0]["details"]["operation"] == "create_output_parent"
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
def test_raw_os_argument_bytes_still_produce_utf8_json() -> None:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [
            os.fsencode(sys.executable),
            b"-m",
            b"rnaseq_downstream",
            b"inspect",
            b"--request",
            b"\xff",
        ],
        cwd=PROJECT_ROOT,
        env=environment,
        check=False,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=CLI_TIMEOUT_SECONDS,
    )

    decoded = completed.stdout.decode("utf-8")
    document = json.loads(decoded)
    assert completed.returncode == 2
    assert_envelope_shape(document, command="inspect", status="error")
    assert document["errors"][0]["code"] == "INPUT_READ_FAILED"
    assert completed.stderr == b""


@pytest.mark.integration
@pytest.mark.parametrize(
    "arguments",
    [
        ("inspect",),
        ("validate",),
        ("validate", "--request", "request.json"),
    ],
)
def test_input_commands_require_explicit_noninteractive_arguments(
    run_module_cli,
    arguments: tuple[str, ...],
) -> None:
    result = run_module_cli(*arguments)

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command=arguments[0], status="error")
    assert document["errors"][0]["code"] == "INVALID_REQUEST"
    assert "usage" in document["errors"][0]["details"]
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
