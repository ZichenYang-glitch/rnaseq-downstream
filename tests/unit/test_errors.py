"""Unit tests for stable public error and process-exit semantics."""

from __future__ import annotations

from pathlib import Path

import pytest

from rnaseq_downstream.errors import (
    BackendFailedError,
    AssayProtocolRequiredError,
    CountValuesInvalidError,
    ContrastNotEstimableError,
    CovariateConfoundedError,
    DesignRankDeficientError,
    ErrorCode,
    ExitCode,
    FeatureNotImplementedError,
    GeneIdentifierError,
    InputEvidenceRequiredError,
    InputIntegrityError,
    InputReadError,
    InputValidationError,
    InternalToolkitError,
    InvalidRequestError,
    OutputWriteError,
    PartialRunError,
    QCValidationError,
    SampleSetMismatchError,
    SalmonOffsetRequiredError,
    ToolkitError,
)


@pytest.mark.unit
def test_exit_code_values_are_stable() -> None:
    assert {member.name: int(member) for member in ExitCode} == {
        "SUCCESS": 0,
        "REQUEST_ERROR": 2,
        "SCIENTIFIC_VALIDATION_ERROR": 3,
        "BACKEND_ERROR": 4,
        "PARTIAL": 5,
        "INTERNAL_ERROR": 70,
    }


@pytest.mark.unit
@pytest.mark.parametrize(
    ("error_type", "code", "exit_code", "status"),
    [
        (
            InvalidRequestError,
            ErrorCode.INVALID_REQUEST,
            ExitCode.REQUEST_ERROR,
            "error",
        ),
        (
            InputValidationError,
            ErrorCode.INPUT_VALIDATION_FAILED,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            InputEvidenceRequiredError,
            ErrorCode.INPUT_EVIDENCE_REQUIRED,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            InputIntegrityError,
            ErrorCode.INPUT_INTEGRITY_FAILED,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            CountValuesInvalidError,
            ErrorCode.COUNT_VALUES_INVALID,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            SampleSetMismatchError,
            ErrorCode.SAMPLE_SET_MISMATCH,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            GeneIdentifierError,
            ErrorCode.GENE_ID_INVALID,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            AssayProtocolRequiredError,
            ErrorCode.ASSAY_PROTOCOL_REQUIRED,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            QCValidationError,
            ErrorCode.QC_VALIDATION_FAILED,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            SalmonOffsetRequiredError,
            ErrorCode.SALMON_OFFSET_REQUIRED,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            DesignRankDeficientError,
            ErrorCode.DESIGN_RANK_DEFICIENT,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            CovariateConfoundedError,
            ErrorCode.COVARIATE_CONFOUNDED,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (
            ContrastNotEstimableError,
            ErrorCode.CONTRAST_NOT_ESTIMABLE,
            ExitCode.SCIENTIFIC_VALIDATION_ERROR,
            "error",
        ),
        (BackendFailedError, ErrorCode.BACKEND_FAILED, ExitCode.BACKEND_ERROR, "error"),
        (PartialRunError, ErrorCode.PARTIAL_RUN, ExitCode.PARTIAL, "partial"),
        (
            FeatureNotImplementedError,
            ErrorCode.FEATURE_NOT_IMPLEMENTED,
            ExitCode.REQUEST_ERROR,
            "error",
        ),
        (
            InternalToolkitError,
            ErrorCode.INTERNAL_ERROR,
            ExitCode.INTERNAL_ERROR,
            "error",
        ),
    ],
)
def test_typed_errors_have_fixed_machine_semantics(
    error_type: type[ToolkitError],
    code: ErrorCode,
    exit_code: ExitCode,
    status: str,
) -> None:
    error = error_type("Evidence-bearing failure.", details={"stage": "test"})

    assert error.code is code
    assert error.exit_code is exit_code
    assert error.status == status
    assert error.to_dict() == {
        "code": code.value,
        "message": "Evidence-bearing failure.",
        "details": {"stage": "test"},
    }


@pytest.mark.unit
def test_input_read_error_preserves_source_without_leaking_exception_objects() -> None:
    cause = OSError("permission denied")
    error = InputReadError(
        "Could not read metadata.",
        path=Path("inputs/metadata.tsv"),
        operation="load_metadata",
        cause=cause,
        details={"design_column": "condition"},
    )

    assert error.code is ErrorCode.INPUT_READ_FAILED
    assert error.exit_code is ExitCode.REQUEST_ERROR
    assert error.cause is cause
    assert error.details == {
        "path": "inputs/metadata.tsv",
        "operation": "load_metadata",
        "cause_type": "OSError",
        "cause_message": "permission denied",
        "design_column": "condition",
    }
    assert error.to_dict()["details"] == error.details


@pytest.mark.unit
def test_output_write_error_preserves_stable_request_failure_context() -> None:
    cause = OSError("read-only filesystem")
    error = OutputWriteError(
        "Could not publish validation evidence.",
        path=Path("evidence"),
        operation="publish_validation_bundle",
        cause=cause,
    )

    assert error.code is ErrorCode.OUTPUT_WRITE_FAILED
    assert error.exit_code is ExitCode.REQUEST_ERROR
    assert error.cause is cause
    assert error.details == {
        "path": "evidence",
        "operation": "publish_validation_bundle",
        "cause_type": "OSError",
        "cause_message": "read-only filesystem",
    }


@pytest.mark.unit
def test_error_details_are_copied_at_input_and_output_boundaries() -> None:
    original = {"sample": "s1"}
    error = InvalidRequestError("Invalid sample.", details=original)
    original["sample"] = "mutated"
    serialized = error.to_dict()
    serialized["details"]["sample"] = "also-mutated"

    assert error.details == {"sample": "s1"}
