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


def _runnable_featurecounts_request(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    samples = [f"s{index + 1}" for index in range(6)]
    (root / "metadata.tsv").write_text(
        "sample_id\tcondition\n"
        + "\n".join(
            f"{sample}\t{'control' if index < 3 else 'treated'}"
            for index, sample in enumerate(samples)
        )
        + "\n",
        encoding="utf-8",
    )
    reference = root / "reference.fa"
    reference.write_text(">synthetic\nACGT\n", encoding="utf-8")
    matrix = root / "counts.tsv"
    lines = ["\t".join(["gene_id", "gene_name", *samples])]
    for index in range(40):
        baseline = 40 + (index % 7)
        values = [baseline] * 3 + [baseline * 4 if index < 8 else baseline] * 3
        lines.append(
            "\t".join(
                [f"gene_{index + 1}", f"GENE_{index + 1}", *(str(v) for v in values)]
            )
        )
    matrix.write_text("\n".join(lines) + "\n", encoding="utf-8")
    _write_json(
        root / "counts.manifest.json",
        {
            "schema_version": "1.0",
            "artifact_type": "featurecounts_integer_matrix",
            "artifact": {"path": "counts.tsv", "sha256": _sha256(matrix)},
            "gene_id_column": "gene_id",
            "display_columns": ["gene_name"],
            "sample_columns": samples,
            "producer": {"name": "featureCounts", "version": "2.0.6"},
            "reference": {
                "name": "synthetic",
                "version": "1",
                "source": "reference.fa",
                "sha256": _sha256(reference),
            },
        },
    )
    request = root / "input-request.json"
    _write_json(
        request,
        {
            "schema_version": "1.0",
            "input_semantics": "featurecounts_integer",
            "metadata": {"path": "metadata.tsv", "sample_id_column": "sample_id"},
            "producer": {"name": "featureCounts", "version": "2.0.6"},
            "reference": {
                "name": "synthetic",
                "version": "1",
                "source": "reference.fa",
                "sha256": _sha256(reference),
            },
            "gene_id": {"strip_version": False},
            "samples": [{"sample_id": sample} for sample in samples],
            "featurecounts": {
                "layout": "combined_matrix",
                "matrix": "counts.tsv",
                "manifest": "counts.manifest.json",
            },
        },
    )
    return request


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
    assert document["data"]["analysis_request_schema_versions"] == ["1.0", "1.1"]
    assert document["data"]["certified_analysis_paths"] == []
    paths = document["data"]["evidence_gated_analysis_paths"]
    assert [path["path_id"] for path in paths] == [
        "edger_ql_p0_v1",
        "edger_ql_p0_v1_limma_gene_sets_v1",
    ]
    assert paths[0]["maturity"] == "research_preview"
    assert paths[0]["benchmark_scope"] == "backend_kernel_only"
    assert paths[0]["combined_manifest_origin_authentication"] == (
        "self_attested_not_proven"
    )
    assert paths[0]["publication_grade_claim"] is False
    assert paths[1]["parent_path_id"] == "edger_ql_p0_v1"
    assert paths[1]["self_contained"] == {
        "primary": "limma_fry",
        "corroborative": "limma_mroast",
    }
    assert paths[1]["competitive"] == {"supplementary": "limma_camera"}
    assert paths[1]["benchmark_evidence"] == [
        "airway-limma-gene-set-same-engine-v1",
        "compcoder-limma-self-contained-fdr-tpr-v1",
    ]
    assert paths[1]["publication_grade_claim"] is False
    displays = document["data"]["non_statistical_display_capabilities"]
    assert [item["capability_id"] for item in displays] == [
        "edger_ql_p0_v1_static_svg_display_v1"
    ]
    assert displays[0]["analysis_path_id"] == "edger_ql_p0_v1"
    assert displays[0]["analysis_request_schema_version"] == "1.1"
    assert displays[0]["statistical_role"] == "display_only_no_inference"
    assert displays[0]["plot_types"] == {
        "volcano": "one_per_contrast",
        "ma": "one_per_contrast",
        "pca": "one_per_analysis",
    }
    assert displays[0]["verification"] == "summarize_source_reproduction"
    assert displays[0]["publication_grade_claim"] is False
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
def test_legacy_checkout_resolves_core_errors_without_site_packages() -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    source = (
        "import json, pathlib, sys, types; "
        "sys.modules.update({name: types.ModuleType(name) "
        "for name in ('pandas', 'yaml')}); "
        "from legacy.modules import data; "
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
@pytest.mark.parametrize(
    ("arguments", "command"),
    [
        (("run",), "run"),
        (("run", "--request", "analysis-request.json"), "run"),
        (("summarize",), "summarize"),
    ],
)
def test_analysis_commands_require_explicit_noninteractive_arguments(
    run_module_cli,
    arguments: tuple[str, ...],
    command: str,
) -> None:
    result = run_module_cli(*arguments)

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command=command, status="error")
    assert document["errors"][0]["code"] == "INVALID_REQUEST"
    assert "usage" in document["errors"][0]["details"]
    assert result.stderr == ""


@pytest.mark.integration
def test_public_summarize_missing_bundle_is_structured_json(
    run_module_cli, tmp_path: Path
) -> None:
    result = run_module_cli(
        "summarize", "--run-dir", str(tmp_path / "missing-results"), cwd=tmp_path
    )

    assert result.returncode == 2
    document = result.json()
    assert_envelope_shape(document, command="summarize", status="error")
    assert document["errors"][0]["code"] == "INPUT_READ_FAILED"
    assert document["errors"][0]["details"]["operation"] == "open_result_bundle"
    assert result.stderr == ""


@pytest.mark.integration
def test_packaged_r_backend_is_readable_outside_checkout_cwd(tmp_path: Path) -> None:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    source = (
        "import hashlib, importlib.resources, json; "
        "resource = importlib.resources.files('rnaseq_downstream').joinpath("
        "'r_scripts', 'edger_ql.R'); "
        "content = resource.read_bytes(); "
        "print(json.dumps({'name': resource.name, 'size': len(content), "
        "'sha256': hashlib.sha256(content).hexdigest()}))"
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
    document = json.loads(completed.stdout)
    assert document["name"] == "edger_ql.R"
    assert document["size"] > 1000
    assert len(document["sha256"]) == 64


@pytest.mark.integration
@pytest.mark.parametrize("schema_version", ["1.0", "1.1"])
def test_locked_cli_chain_validate_run_and_summarize(
    run_module_cli, tmp_path: Path, schema_version: str
) -> None:
    rscript = os.environ.get("RNASEQ_P0_RSCRIPT")
    r_library = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not rscript or not r_library:
        if os.environ.get("RNASEQ_P0_REQUIRE_BENCHMARKS") == "1":
            pytest.fail(
                "Certification mode requires RNASEQ_P0_RSCRIPT and "
                "RNASEQ_P0_R_LIBRARY; the public CLI chain cannot skip"
            )
        pytest.skip("the locked R runtime is required for the public CLI chain")

    input_request = _runnable_featurecounts_request(tmp_path / "inputs")
    inspect_result = run_module_cli(
        "inspect", "--request", str(input_request), cwd=tmp_path
    )
    assert inspect_result.returncode == 0
    assert_envelope_shape(inspect_result.json(), command="inspect", status="success")

    evidence = tmp_path / "evidence"
    validate_result = run_module_cli(
        "validate",
        "--request",
        str(input_request),
        "--output-dir",
        str(evidence),
        cwd=tmp_path,
    )
    assert validate_result.returncode == 0
    assert_envelope_shape(validate_result.json(), command="validate", status="success")

    analysis_request = tmp_path / "analysis-request.json"
    request: dict[str, object] = {
        "schema_version": schema_version,
        "validated_input_bundle": str(evidence),
        "design": {
            "intercept": True,
            "terms": ["condition"],
            "variables": {
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                }
            },
        },
        "contrasts": [
            {
                "contrast_id": "treated_vs_control",
                "weights": {"conditiontreated": 1},
                "lfc_threshold": 0,
            }
        ],
    }
    if schema_version == "1.1":
        request["display"] = {
            "fdr_threshold": 0.05,
            "pca_top_n": 20,
            "pca_components": [1, 2],
        }
    _write_json(analysis_request, request)
    run_dir = tmp_path / "results"
    run_result = run_module_cli(
        "run",
        "--request",
        str(analysis_request),
        "--output-dir",
        str(run_dir),
        "--rscript",
        rscript,
        "--r-library",
        r_library,
        cwd=tmp_path,
    )
    assert run_result.returncode == 0, (run_result.stdout, run_result.stderr)
    run_document = run_result.json()
    assert_envelope_shape(run_document, command="run", status="success")
    assert run_document["data"]["scope"]["execution_scope"] == ("validated_p0_input")
    assert "backend_stderr" not in run_document["data"]
    assert run_document["data"]["ql_fit_parameters"]["abundance.trend"] is True
    expected_entries = [
        "analysis.json",
        "backend_manifest.json",
        "coefficients.tsv",
        "design.tsv",
        "results.tsv",
    ]
    if schema_version == "1.1":
        expected_entries.append("display")
        assert run_document["data"]["display_export"]["purpose"] == (
            "display_only_not_for_inference"
        )
    else:
        assert run_document["data"]["display_export"] is None
    assert sorted(path.name for path in run_dir.iterdir()) == sorted(expected_entries)
    expected_artifact_count = 12 if schema_version == "1.1" else 5
    assert len(run_document["artifacts"]) == expected_artifact_count

    summary_result = run_module_cli(
        "summarize", "--run-dir", str(run_dir), cwd=tmp_path
    )
    assert summary_result.returncode == 0, (
        summary_result.stdout,
        summary_result.stderr,
    )
    summary_document = summary_result.json()
    assert_envelope_shape(summary_document, command="summarize", status="success")
    assert summary_document["data"]["status"] == "verified_complete"
    assert summary_document["data"]["execution_scope"] == "validated_p0_input"
    assert summary_document["data"]["gene_count"] == 40
    assert summary_document["data"]["result_row_count"] == 40
    assert len(summary_document["artifacts"]) == expected_artifact_count
    if schema_version == "1.1":
        assert summary_document["data"]["display"]["plot_types"] == {
            "pca": 1,
            "volcano": 1,
            "ma": 1,
        }
        # The byte-identical legacy core intentionally carries no display
        # intent; this documented compatibility boundary must remain explicit.
        (run_dir / "display").rename(tmp_path / "detached-display")
        detached_result = run_module_cli(
            "summarize", "--run-dir", str(run_dir), cwd=tmp_path
        )
        assert detached_result.returncode == 0
        detached_document = detached_result.json()
        assert "display" not in detached_document["data"]
        assert len(detached_document["artifacts"]) == 5
    else:
        assert "display" not in summary_document["data"]
    assert summary_result.stderr == ""


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
