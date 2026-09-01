"""Stable error types shared by library and command-line boundaries."""

from __future__ import annotations

from enum import Enum, IntEnum
from pathlib import Path
from typing import Any, Mapping


class ErrorCode(str, Enum):
    """Machine-readable error codes in CLI contract version 1.0."""

    INVALID_REQUEST = "INVALID_REQUEST"
    INPUT_READ_FAILED = "INPUT_READ_FAILED"
    SALMON_OFFSET_REQUIRED = "SALMON_OFFSET_REQUIRED"
    DESIGN_RANK_DEFICIENT = "DESIGN_RANK_DEFICIENT"
    CONTRAST_NOT_ESTIMABLE = "CONTRAST_NOT_ESTIMABLE"
    BACKEND_FAILED = "BACKEND_FAILED"
    PARTIAL_RUN = "PARTIAL_RUN"
    FEATURE_NOT_IMPLEMENTED = "FEATURE_NOT_IMPLEMENTED"
    INTERNAL_ERROR = "INTERNAL_ERROR"


class ExitCode(IntEnum):
    """Process exit codes in CLI contract version 1.0."""

    SUCCESS = 0
    REQUEST_ERROR = 2
    SCIENTIFIC_VALIDATION_ERROR = 3
    BACKEND_ERROR = 4
    PARTIAL = 5
    INTERNAL_ERROR = 70


class ToolkitError(Exception):
    """Base class for expected, machine-readable toolkit failures.

    Subclasses fix the error and exit codes at class definition time. ``cause`` is
    retained for Python callers and exception chaining; it is never serialized
    implicitly.
    """

    code = ErrorCode.INTERNAL_ERROR
    exit_code = ExitCode.INTERNAL_ERROR
    status = "error"

    def __init__(
        self,
        message: str,
        *,
        details: Mapping[str, Any] | None = None,
        cause: BaseException | None = None,
    ) -> None:
        super().__init__(message)
        self.message = str(message)
        self.details = dict(details or {})
        self.cause = cause

    def to_dict(self) -> dict[str, Any]:
        """Return the stable JSON error representation."""

        return {
            "code": self.code.value,
            "message": self.message,
            "details": dict(self.details),
        }


class InvalidRequestError(ToolkitError):
    """The request cannot be parsed or does not satisfy the CLI schema."""

    code = ErrorCode.INVALID_REQUEST
    exit_code = ExitCode.REQUEST_ERROR


class InputReadError(ToolkitError):
    """An explicitly named input could not be read or parsed."""

    code = ErrorCode.INPUT_READ_FAILED
    exit_code = ExitCode.REQUEST_ERROR

    def __init__(
        self,
        message: str,
        *,
        path: str | Path,
        operation: str,
        cause: BaseException,
        details: Mapping[str, Any] | None = None,
    ) -> None:
        error_details: dict[str, Any] = {
            "path": str(path),
            "operation": operation,
            "cause_type": type(cause).__name__,
            "cause_message": str(cause),
        }
        error_details.update(details or {})
        super().__init__(message, details=error_details, cause=cause)


class SalmonOffsetRequiredError(ToolkitError):
    """A full-length Salmon request lacks the required length-offset evidence."""

    code = ErrorCode.SALMON_OFFSET_REQUIRED
    exit_code = ExitCode.SCIENTIFIC_VALIDATION_ERROR


class DesignRankDeficientError(ToolkitError):
    """The declared model design is not full rank."""

    code = ErrorCode.DESIGN_RANK_DEFICIENT
    exit_code = ExitCode.SCIENTIFIC_VALIDATION_ERROR


class ContrastNotEstimableError(ToolkitError):
    """The requested contrast cannot be estimated from the model design."""

    code = ErrorCode.CONTRAST_NOT_ESTIMABLE
    exit_code = ExitCode.SCIENTIFIC_VALIDATION_ERROR


class BackendFailedError(ToolkitError):
    """The selected statistical backend failed to complete."""

    code = ErrorCode.BACKEND_FAILED
    exit_code = ExitCode.BACKEND_ERROR


class PartialRunError(ToolkitError):
    """A run produced explicitly incomplete output."""

    code = ErrorCode.PARTIAL_RUN
    exit_code = ExitCode.PARTIAL
    status = "partial"


class FeatureNotImplementedError(ToolkitError):
    """A reserved command surface is not available in this checkpoint."""

    code = ErrorCode.FEATURE_NOT_IMPLEMENTED
    exit_code = ExitCode.REQUEST_ERROR


class InternalToolkitError(ToolkitError):
    """An unexpected internal failure hidden behind a stable public code."""

    code = ErrorCode.INTERNAL_ERROR
    exit_code = ExitCode.INTERNAL_ERROR
