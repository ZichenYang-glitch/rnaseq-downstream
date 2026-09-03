"""Unit tests for the real CLI execution and exit-code boundary."""

from __future__ import annotations

import json
import logging
import math
from types import SimpleNamespace

import pytest

from rnaseq_downstream import cli
from rnaseq_downstream import validation_bundle
from rnaseq_downstream.contracts import Status, build_envelope
from rnaseq_downstream.errors import (
    BackendFailedError,
    DesignRankDeficientError,
    ExitCode,
    PartialRunError,
    ToolkitError,
)


@pytest.mark.unit
def test_capabilities_exposes_v11_display_without_changing_analysis_path() -> None:
    data = cli._capabilities(SimpleNamespace())

    assert data["analysis_request_schema_versions"] == ["1.0", "1.1"]
    assert [path["path_id"] for path in data["evidence_gated_analysis_paths"]] == [
        "edger_ql_p0_v1"
    ]
    assert data["certified_analysis_paths"] == []
    assert data["non_statistical_display_capabilities"] == [
        {
            "capability_id": "edger_ql_p0_v1_static_svg_display_v1",
            "maturity": "research_preview",
            "analysis_path_id": "edger_ql_p0_v1",
            "analysis_request_schema_version": "1.1",
            "invocation": "optional_same_run",
            "statistical_role": "display_only_no_inference",
            "output_location": "display/",
            "output_format": "svg",
            "plot_types": {
                "volcano": "one_per_contrast",
                "ma": "one_per_contrast",
                "pca": "one_per_analysis",
            },
            "pca_input": "post_filter_post_tmm_edger_logcpm",
            "pca_scaling": "centered_unscaled",
            "determinism_scope": "locked_runtime",
            "verification": "summarize_source_reproduction",
            "publication_grade_claim": False,
        }
    ]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("command", "failure", "expected_exit", "expected_status", "expected_code"),
    [
        (
            "validate",
            DesignRankDeficientError("The design is not full rank."),
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
            "DESIGN_RANK_DEFICIENT",
        ),
        (
            "run",
            BackendFailedError("edgeR did not produce a result."),
            ExitCode.BACKEND_ERROR,
            "error",
            "BACKEND_FAILED",
        ),
        (
            "run",
            PartialRunError("At least one requested contrast failed."),
            ExitCode.PARTIAL,
            "partial",
            "PARTIAL_RUN",
        ),
    ],
)
def test_execute_maps_typed_failures_to_public_contract(
    monkeypatch: pytest.MonkeyPatch,
    command: str,
    failure: ToolkitError,
    expected_exit: ExitCode,
    expected_status: str,
    expected_code: str,
) -> None:
    def fail_command(_arguments) -> dict:
        raise failure

    if command == "validate":
        monkeypatch.setattr(cli, "_validate", fail_command)
        arguments = [
            command,
            "--request",
            "request.json",
            "--output-dir",
            "evidence",
        ]
    else:
        monkeypatch.setattr(cli, "_run", fail_command)
        arguments = [
            command,
            "--request",
            "analysis-request.json",
            "--output-dir",
            "results",
        ]

    envelope, exit_code = cli._execute(arguments)

    assert exit_code is expected_exit
    assert envelope["command"] == command
    assert envelope["status"] == expected_status
    assert envelope["errors"][0]["code"] == expected_code


@pytest.mark.unit
def test_execute_maps_unexpected_failure_to_internal_error(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def fail_unexpectedly(_arguments) -> dict:
        raise RuntimeError("private implementation detail")

    monkeypatch.setattr(cli, "_run", fail_unexpectedly)

    with caplog.at_level(logging.ERROR, logger=cli.__name__):
        envelope, exit_code = cli._execute(
            [
                "run",
                "--request",
                "analysis-request.json",
                "--output-dir",
                "results",
            ]
        )

    assert exit_code is ExitCode.INTERNAL_ERROR
    assert envelope["status"] == "error"
    assert envelope["errors"] == [
        {
            "code": "INTERNAL_ERROR",
            "message": "An unexpected internal error occurred.",
            "details": {},
        }
    ]
    assert "Unhandled CLI failure" in caplog.text


@pytest.mark.unit
def test_main_returns_internal_exit_when_serialization_falls_back(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    unsafe = build_envelope(
        "capabilities",
        status=Status.SUCCESS,
        data={"unsafe": math.nan},
    )
    monkeypatch.setattr(
        cli,
        "_execute",
        lambda _arguments: (unsafe, ExitCode.SUCCESS),
    )

    exit_code = cli.main(["capabilities"])
    captured = capsys.readouterr()
    document = json.loads(captured.out)

    assert exit_code == int(ExitCode.INTERNAL_ERROR)
    assert document["status"] == "error"
    assert document["errors"][0]["code"] == "INTERNAL_ERROR"
    assert captured.out.count("\n") == 1
    assert captured.err == ""


@pytest.mark.unit
def test_validate_scope_preserves_ineligible_input_status(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    validation_result = {
        "data": {
            "input_certification_status": "ineligible",
            "input_certification_eligible": False,
        },
        "warnings": [{"code": "HIGH_RISK_OVERRIDE", "severity": "high"}],
        "artifacts": [],
    }
    monkeypatch.setattr(
        validation_bundle,
        "validate_request_to_bundle",
        lambda _request, _output: {
            "validation_result": validation_result,
            "bundle": {
                "artifacts": [],
                "validation": {},
                "warnings": [
                    {"code": "OUTPUT_DURABILITY_UNCONFIRMED", "severity": "high"}
                ],
            },
        },
    )

    result = cli._validate(
        SimpleNamespace(request="request.json", output_dir="evidence")
    )

    assert result.data["scope"]["input_semantics"] == "ineligible"
    assert result.data["scope"]["full_numeric_validation"] == "passed"
    assert [warning["code"] for warning in result.warnings] == [
        "HIGH_RISK_OVERRIDE",
        "OUTPUT_DURABILITY_UNCONFIRMED",
    ]


@pytest.mark.unit
def test_run_forwards_backend_logs_to_stderr_and_returns_artifacts(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from rnaseq_downstream import edger_backend

    monkeypatch.setattr(
        edger_backend,
        "run_edger_ql",
        lambda *_args, **_kwargs: {
            "status": "success",
            "backend": "edgeR_QL",
            "execution_scope": "validated_p0_input",
            "output_dir": "/results",
            "publication_status": "durability_confirmed",
            "plan_id": "a" * 64,
            "analysis_request_sha256": "b" * 64,
            "data": {
                "runtime_identity": {"edgeR": "4.10.0"},
                "input_semantics": "featurecounts_integer",
                "route_observed": {"constructor": "edgeR::DGEList"},
                "design_columns": ["(Intercept)", "conditiontreated"],
                "design_rank": 2,
                "residual_df": 4,
                "gene_count": 100,
                "tested_gene_count": 80,
                "filtered_gene_count": 20,
                "ql_fit_parameters": {
                    "abundance.trend": True,
                    "robust": True,
                    "winsor.tail.p": [0.05, 0.1],
                    "legacy": False,
                    "top.proportion": None,
                    "resolved_top.proportion": 0.1,
                    "keep.unit.mat": False,
                },
                "display_export": None,
                "contrasts": [{"contrast_id": "treated_vs_control"}],
            },
            "warnings": [],
            "artifacts": [{"role": "results", "path": "/results/results.tsv"}],
            "backend_stderr": "edgeR diagnostic\n",
        },
    )

    envelope, exit_code = cli._execute(
        [
            "run",
            "--request",
            "analysis-request.json",
            "--output-dir",
            "results",
            "--rscript",
            "/locked/Rscript",
            "--r-library",
            "/locked/library",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code is ExitCode.SUCCESS
    assert envelope["status"] == "success"
    assert envelope["data"]["scope"]["analysis_path"] == "edger_ql_p0_v1"
    assert envelope["data"]["scope"]["benchmark_scope"] == "backend_kernel_only"
    assert envelope["data"]["scope"]["combined_manifest_origin_authentication"] == (
        "self_attested_not_proven"
    )
    assert envelope["data"]["scope"]["publication_grade_claim"] is False
    assert envelope["data"]["ql_fit_parameters"]["abundance.trend"] is True
    assert envelope["data"]["ql_fit_parameters"]["winsor.tail.p"] == [0.05, 0.1]
    assert envelope["data"]["display_export"] is None
    assert envelope["artifacts"][0]["role"] == "results"
    assert captured.out == ""
    assert captured.err == "edgeR diagnostic\n"


@pytest.mark.unit
def test_summarize_handler_exposes_verified_data_and_consumed_artifacts(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from rnaseq_downstream import run_summary

    monkeypatch.setattr(
        run_summary,
        "summarize_run",
        lambda _path: {
            "schema_version": "1.0",
            "status": "verified_complete",
            "contrasts": [],
            "artifacts": [{"role": "backend_manifest"}],
        },
    )

    envelope, exit_code = cli._execute(["summarize", "--run-dir", "results"])

    assert exit_code is ExitCode.SUCCESS
    assert envelope["status"] == "success"
    assert envelope["data"]["status"] == "verified_complete"
    assert "artifacts" not in envelope["data"]
    assert envelope["artifacts"] == [{"role": "backend_manifest"}]
