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
    def fail_reserved_command(_arguments) -> dict:
        raise failure

    if command == "validate":
        monkeypatch.setattr(cli, "_validate", fail_reserved_command)
        arguments = [
            command,
            "--request",
            "request.json",
            "--output-dir",
            "evidence",
        ]
    else:
        monkeypatch.setattr(cli, "_reserved_command", fail_reserved_command)
        arguments = [command]

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

    monkeypatch.setattr(cli, "_reserved_command", fail_unexpectedly)

    with caplog.at_level(logging.ERROR, logger=cli.__name__):
        envelope, exit_code = cli._execute(["run"])

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
