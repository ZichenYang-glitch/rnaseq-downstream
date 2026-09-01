"""Unit tests for the complete and fail-closed JSON response envelope."""

from __future__ import annotations

from io import StringIO
import json
import math

import pytest

from rnaseq_downstream.contracts import (
    SCHEMA_VERSION,
    Status,
    build_envelope,
    serialize_envelope,
    write_json_document,
)
from rnaseq_downstream.errors import PartialRunError

EXPECTED_KEYS = [
    "schema_version",
    "command",
    "status",
    "data",
    "warnings",
    "errors",
    "artifacts",
]


@pytest.mark.unit
def test_build_envelope_never_omits_contract_fields() -> None:
    envelope = build_envelope("inspect", status=Status.SUCCESS)

    assert list(envelope) == EXPECTED_KEYS
    assert envelope == {
        "schema_version": SCHEMA_VERSION,
        "command": "inspect",
        "status": "success",
        "data": None,
        "warnings": [],
        "errors": [],
        "artifacts": [],
    }


@pytest.mark.unit
def test_envelope_rejects_unknown_status() -> None:
    with pytest.raises(ValueError, match="Unsupported response status"):
        build_envelope("run", status="ambiguous")


@pytest.mark.unit
@pytest.mark.parametrize("status", [Status.ERROR, Status.PARTIAL])
def test_non_success_envelope_requires_an_error(status: Status) -> None:
    with pytest.raises(ValueError, match="must contain at least one error"):
        build_envelope("run", status=status)


@pytest.mark.unit
def test_success_envelope_rejects_errors() -> None:
    with pytest.raises(ValueError, match="successful response cannot contain errors"):
        build_envelope(
            "validate",
            status=Status.SUCCESS,
            errors=[{"code": "IMPOSSIBLE", "message": "contradiction"}],
        )


@pytest.mark.unit
def test_partial_response_is_explicit_and_preserves_available_evidence() -> None:
    failure = PartialRunError(
        "One requested contrast failed.",
        details={"failed_contrasts": ["treated_vs_control"]},
    )
    envelope = build_envelope(
        "run",
        status=failure.status,
        data={"completed_contrasts": ["dose_vs_control"]},
        warnings=[{"code": "INCOMPLETE_OUTPUT", "message": "Do not certify."}],
        errors=[failure.to_dict()],
        artifacts=[{"kind": "provenance", "path": "run/provenance.json"}],
    )

    assert envelope["status"] == "partial"
    assert envelope["errors"][0]["code"] == "PARTIAL_RUN"
    assert envelope["data"]["completed_contrasts"] == ["dose_vs_control"]
    assert envelope["artifacts"] == [
        {"kind": "provenance", "path": "run/provenance.json"}
    ]


@pytest.mark.unit
def test_serializer_is_compact_valid_json() -> None:
    envelope = build_envelope(
        "capabilities",
        status=Status.SUCCESS,
        data={"unicode": "基因"},
    )

    document = serialize_envelope(envelope)

    assert json.loads(document) == envelope
    assert document.isascii()
    assert "\\u57fa\\u56e0" in document
    assert "\n" not in document
    assert ": " not in document


@pytest.mark.unit
def test_nonfinite_number_fails_closed_to_internal_error() -> None:
    unsafe = build_envelope(
        "summarize",
        status=Status.SUCCESS,
        data={"fdr": math.nan},
    )

    serialized = serialize_envelope(unsafe)
    fallback = json.loads(serialized)

    assert "NaN" not in serialized
    assert fallback["command"] == "summarize"
    assert fallback["status"] == "error"
    assert fallback["data"] is None
    assert fallback["errors"][0]["code"] == "INTERNAL_ERROR"
    assert fallback["errors"][0]["details"]["cause_type"] == "ValueError"


@pytest.mark.unit
def test_writer_reports_fallback_and_still_emits_one_document() -> None:
    stream = StringIO()
    unsafe = build_envelope(
        "run",
        status=Status.SUCCESS,
        data={"coefficient": math.inf},
    )

    serialized_safely = write_json_document(unsafe, stream=stream)

    assert serialized_safely is False
    assert stream.getvalue().endswith("\n")
    assert stream.getvalue().count("\n") == 1
    assert json.loads(stream.getvalue())["errors"][0]["code"] == "INTERNAL_ERROR"


@pytest.mark.unit
def test_writer_reports_primary_serialization_success() -> None:
    stream = StringIO()
    safe = build_envelope("validate", status=Status.SUCCESS, data={"valid": True})

    serialized_safely = write_json_document(safe, stream=stream)

    assert serialized_safely is True
    assert json.loads(stream.getvalue()) == safe
