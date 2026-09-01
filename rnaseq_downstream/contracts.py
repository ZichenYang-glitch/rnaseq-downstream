"""JSON document construction and fail-closed serialization."""

from __future__ import annotations

import json
import sys
from enum import Enum
from typing import Any, Mapping, Sequence, TextIO

from .errors import ErrorCode


SCHEMA_VERSION = "1.0"


class Status(str, Enum):
    """Allowed top-level command statuses."""

    SUCCESS = "success"
    ERROR = "error"
    PARTIAL = "partial"


def build_envelope(
    command: str,
    *,
    status: Status | str,
    data: Any = None,
    warnings: Sequence[Mapping[str, Any]] | None = None,
    errors: Sequence[Mapping[str, Any]] | None = None,
    artifacts: Sequence[Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Build a complete, internally consistent CLI response envelope."""

    normalized_status = status.value if isinstance(status, Status) else str(status)
    allowed = {item.value for item in Status}
    if normalized_status not in allowed:
        raise ValueError(f"Unsupported response status: {normalized_status}")

    normalized_errors = list(errors or [])
    if normalized_status == Status.SUCCESS.value and normalized_errors:
        raise ValueError("A successful response cannot contain errors")
    if normalized_status in {Status.ERROR.value, Status.PARTIAL.value}:
        if not normalized_errors:
            raise ValueError(
                f"A '{normalized_status}' response must contain at least one error"
            )

    return {
        "schema_version": SCHEMA_VERSION,
        "command": str(command),
        "status": normalized_status,
        "data": data,
        "warnings": list(warnings or []),
        "errors": normalized_errors,
        "artifacts": list(artifacts or []),
    }


def _fallback_envelope(envelope: object, error: BaseException) -> dict[str, Any]:
    command = "unknown"
    if isinstance(envelope, Mapping) and isinstance(envelope.get("command"), str):
        command = envelope["command"]
    return build_envelope(
        command,
        status=Status.ERROR,
        errors=[
            {
                "code": ErrorCode.INTERNAL_ERROR.value,
                "message": "The response could not be serialized safely.",
                "details": {"cause_type": type(error).__name__},
            }
        ],
    )


def _serialize_with_fallback(envelope: Mapping[str, Any]) -> tuple[str, bool]:
    try:
        document = json.dumps(
            envelope,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
        )
        return document, False
    except (TypeError, ValueError, OverflowError) as error:
        fallback = _fallback_envelope(envelope, error)
        document = json.dumps(
            fallback,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
        )
        return document, True


def serialize_envelope(envelope: Mapping[str, Any]) -> str:
    """Serialize exactly one compact JSON document, with a safe fallback."""

    document, _ = _serialize_with_fallback(envelope)
    return document


def write_json_document(
    envelope: Mapping[str, Any],
    *,
    stream: TextIO | None = None,
) -> bool:
    """Write one JSON document and return whether fallback was unnecessary."""

    output = stream if stream is not None else sys.stdout
    document, used_fallback = _serialize_with_fallback(envelope)
    output.write(document)
    output.write("\n")
    output.flush()
    return not used_fallback
