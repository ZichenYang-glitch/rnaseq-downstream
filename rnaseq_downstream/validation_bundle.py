"""Crash-durable, non-overwriting input-validation evidence bundles."""

from __future__ import annotations

import ctypes
import errno
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import tempfile
from typing import Any, Mapping, Sequence

from . import __version__
from .errors import (
    InternalToolkitError,
    InvalidRequestError,
    OutputWriteError,
)

BUNDLE_SCHEMA_VERSION = "1.0"
_AT_FDCWD = -100
_RENAME_NOREPLACE = 1


def _canonical_bytes(document: Any) -> bytes:
    try:
        serialized = json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        return serialized.encode("utf-8")
    except (TypeError, ValueError, OverflowError) as error:
        raise InternalToolkitError(
            "Validated input data could not be represented as strict JSON.",
            details={"cause_type": type(error).__name__},
            cause=error,
        ) from error


def compute_validation_plan_id(normalized: Mapping[str, Any]) -> str:
    """Return the deterministic identity of normalized validation evidence.

    ``normalized`` is the exact ``data``/``warnings``/``artifacts`` object
    written into ``validated_request.json``.  Keeping this computation in one
    place lets evidence consumers independently derive the identity instead of
    trusting the digest repeated by bundle members.
    """

    expected_keys = {"data", "warnings", "artifacts"}
    if not isinstance(normalized, Mapping) or set(normalized) != expected_keys:
        raise InternalToolkitError(
            "Normalized validation evidence has an incompatible identity schema.",
            details={
                "expected_keys": sorted(expected_keys),
                "observed_keys": (
                    sorted(str(key) for key in normalized)
                    if isinstance(normalized, Mapping)
                    else None
                ),
                "observed_type": type(normalized).__name__,
            },
        )
    return hashlib.sha256(_canonical_bytes(normalized)).hexdigest()


def canonical_json_equal(left: Any, right: Any) -> bool:
    """Compare JSON values without Python's ``True == 1`` coercion semantics."""

    return _canonical_bytes({"value": left}) == _canonical_bytes({"value": right})


def _write_json(path: Path, document: Mapping[str, Any]) -> tuple[str, int]:
    payload = _canonical_bytes(document) + b"\n"
    try:
        with path.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except OSError as error:
        raise OutputWriteError(
            "A validation artifact could not be written.",
            path=path,
            operation="write_validation_artifact",
            cause=error,
        ) from error
    return hashlib.sha256(payload).hexdigest(), len(payload)


def _fsync_directory(path: Path, *, operation: str) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    try:
        descriptor = os.open(path, flags)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except OSError as error:
        raise OutputWriteError(
            "A validation evidence directory could not be synchronized.",
            path=path,
            operation=operation,
            cause=error,
        ) from error


def _resolve_output_directory(path: str | Path) -> Path:
    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise InvalidRequestError(
            "The validation output directory must be a non-empty path.",
            details={"field": "output_dir"},
        )
    candidate = Path(path)
    if candidate.is_symlink():
        raise InvalidRequestError(
            "The validation output path must not be a symbolic link.",
            details={"output_dir": str(candidate), "reason": "output_symlink"},
        )
    try:
        return candidate.resolve(strict=False)
    except (OSError, RuntimeError) as error:
        raise OutputWriteError(
            "The validation output directory could not be resolved.",
            path=candidate,
            operation="resolve_output_directory",
            cause=error,
        ) from error


def _publish_directory_noreplace(source: Path, target: Path) -> None:
    """Atomically rename ``source`` to absent ``target`` without replacement."""

    if sys.platform.startswith("linux"):
        library = ctypes.CDLL(None, use_errno=True)
        renameat2 = getattr(library, "renameat2", None)
        if renameat2 is None:
            cause = OSError(errno.ENOSYS, "renameat2 is unavailable")
            raise OutputWriteError(
                "This platform cannot publish a directory without replacement.",
                path=target,
                operation="publish_validation_bundle_noreplace",
                cause=cause,
            )
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        result = renameat2(
            _AT_FDCWD,
            os.fsencode(source),
            _AT_FDCWD,
            os.fsencode(target),
            _RENAME_NOREPLACE,
        )
        if result == 0:
            return
        error_number = ctypes.get_errno()
        if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
            raise InvalidRequestError(
                "The validation output path appeared during staging; no evidence was overwritten.",
                details={"output_dir": str(target), "reason": "output_race"},
            )
        cause = OSError(error_number, os.strerror(error_number), str(target))
        raise OutputWriteError(
            "The completed validation bundle could not be published.",
            path=target,
            operation="publish_validation_bundle_noreplace",
            cause=cause,
        )

    if os.name == "nt":
        try:
            os.rename(source, target)
        except FileExistsError as error:
            raise InvalidRequestError(
                "The validation output path appeared during staging; no evidence was overwritten.",
                details={"output_dir": str(target), "reason": "output_race"},
                cause=error,
            ) from error
        except OSError as error:
            raise OutputWriteError(
                "The completed validation bundle could not be published.",
                path=target,
                operation="publish_validation_bundle_noreplace",
                cause=error,
            ) from error
        return

    cause = OSError(errno.ENOTSUP, "atomic no-replace directory rename unavailable")
    raise OutputWriteError(
        "This platform cannot publish a directory without replacement.",
        path=target,
        operation="publish_validation_bundle_noreplace",
        cause=cause,
    )


def _assert_validated_input_data(input_data: Mapping[str, Any]) -> None:
    if not isinstance(input_data, Mapping):
        raise InternalToolkitError(
            "The validator did not return a mapping for validated input data.",
            details={"observed_type": type(input_data).__name__},
        )
    required = {
        "validation_level": "validate",
        "certification_scope": "input_semantics_only",
    }
    mismatches = {
        key: {"expected": expected, "observed": input_data.get(key)}
        for key, expected in required.items()
        if input_data.get(key) != expected
    }
    eligible = input_data.get("input_certification_eligible")
    if not isinstance(eligible, bool):
        mismatches["input_certification_eligible"] = {
            "expected": "boolean",
            "observed": eligible,
        }
    else:
        expected_status = "passed" if eligible else "ineligible"
        if input_data.get("input_certification_status") != expected_status:
            mismatches["input_certification_status"] = {
                "expected": expected_status,
                "observed": input_data.get("input_certification_status"),
            }
    for key in ("input_semantics", "request", "route"):
        if key not in input_data:
            mismatches[key] = {"expected": "present", "observed": "missing"}
    if mismatches:
        raise InternalToolkitError(
            "The validator did not return a complete validated input receipt.",
            details={"mismatches": mismatches},
        )


def _write_validation_bundle(
    output_dir: str | Path,
    *,
    input_data: Mapping[str, Any],
    warnings: Sequence[Mapping[str, Any]],
    consumed_artifacts: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Stage and publish an already validated receipt without replacement."""

    _assert_validated_input_data(input_data)
    target = _resolve_output_directory(output_dir)
    if target.exists() or target.is_symlink():
        raise InvalidRequestError(
            "The validation output path already exists; existing evidence is never overwritten.",
            details={"output_dir": str(target), "reason": "output_exists"},
        )

    parent = target.parent
    try:
        parent.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise OutputWriteError(
            "The parent of the validation output directory could not be created.",
            path=parent,
            operation="create_output_parent",
            cause=error,
        ) from error

    normalized = {
        "data": dict(input_data),
        "warnings": [dict(item) for item in warnings],
        "artifacts": [dict(item) for item in consumed_artifacts],
    }
    plan_id = compute_validation_plan_id(normalized)
    scope = {
        "input_semantics": input_data["input_certification_status"],
        "design": "not_run",
        "backend": "not_run",
        "runnable": False,
        "analysis_path_certified": False,
    }
    documents: dict[str, dict[str, Any]] = {
        "validated_request.json": {
            "schema_version": BUNDLE_SCHEMA_VERSION,
            "kind": "validated_input_request",
            "plan_id": plan_id,
            "validation_scope": "input_semantics_only",
            **normalized,
        },
        "input_plan.json": {
            "schema_version": BUNDLE_SCHEMA_VERSION,
            "kind": "input_plan",
            "plan_id": plan_id,
            "validation_scope": "input_semantics_only",
            "validation": scope,
            "input": dict(input_data),
        },
        "provenance.json": {
            "schema_version": BUNDLE_SCHEMA_VERSION,
            "kind": "input_provenance",
            "plan_id": plan_id,
            "toolkit": {"name": "rnaseq-downstream", "version": __version__},
            "runtime": {
                "python": platform.python_version(),
                "implementation": platform.python_implementation(),
                "platform": platform.platform(),
            },
            "consumed_artifacts": [dict(item) for item in consumed_artifacts],
            "warnings": [dict(item) for item in warnings],
        },
    }

    try:
        staging = Path(
            tempfile.mkdtemp(prefix=f".{target.name}.staging-", dir=str(parent))
        )
    except OSError as error:
        raise OutputWriteError(
            "A validation staging directory could not be created.",
            path=parent,
            operation="create_staging_directory",
            cause=error,
        ) from error

    written: list[dict[str, Any]] = []
    published = False
    publication_status = "durability_confirmed"
    publication_warnings: list[dict[str, Any]] = []
    try:
        manifest_members: list[dict[str, Any]] = []
        for filename in sorted(documents):
            digest, size = _write_json(staging / filename, documents[filename])
            manifest_members.append(
                {"path": filename, "sha256": digest, "size_bytes": size}
            )
            written.append(
                _generated_artifact(target, filename, digest=digest, size=size)
            )

        manifest = {
            "schema_version": BUNDLE_SCHEMA_VERSION,
            "kind": "validation_bundle_manifest",
            "plan_id": plan_id,
            "members": manifest_members,
        }
        manifest_name = "bundle_manifest.json"
        digest, size = _write_json(staging / manifest_name, manifest)
        written.append(
            _generated_artifact(target, manifest_name, digest=digest, size=size)
        )
        _fsync_directory(staging, operation="synchronize_staging_directory")
        _publish_directory_noreplace(staging, target)
        published = True
        try:
            _fsync_directory(parent, operation="synchronize_output_parent")
        except OutputWriteError as error:
            publication_status = "published_durability_unconfirmed"
            publication_warnings.append(
                {
                    "code": "OUTPUT_DURABILITY_UNCONFIRMED",
                    "severity": "high",
                    "message": (
                        "The complete bundle is visible, but synchronizing its parent "
                        "directory failed; crash durability is not confirmed."
                    ),
                    "details": {
                        "output_dir": str(target),
                        "plan_id": plan_id,
                        "manifest": str(target / "bundle_manifest.json"),
                        "operation": error.details.get("operation"),
                        "cause_type": error.details.get("cause_type"),
                    },
                }
            )
    except Exception:
        if not published and staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
        raise

    return {
        "output_dir": str(target),
        "plan_id": plan_id,
        "artifacts": written,
        "validation": scope,
        "publication_status": publication_status,
        "warnings": publication_warnings,
    }


def _generated_artifact(
    target: Path,
    filename: str,
    *,
    digest: str,
    size: int,
) -> dict[str, Any]:
    return {
        "kind": "generated_validation_artifact",
        "role": filename.removesuffix(".json"),
        "path": str(target / filename),
        "sha256": digest,
        "size_bytes": size,
        "media_type": "application/json",
        "schema_version": BUNDLE_SCHEMA_VERSION,
    }


def validate_request_to_bundle(
    request_path: str | Path,
    output_dir: str | Path,
) -> dict[str, Any]:
    """Validate one request and publish its evidence without an inspect shortcut."""

    from .input_semantics import validate_request

    result = validate_request(request_path)
    if not isinstance(result, Mapping):
        raise InternalToolkitError(
            "The input validator returned an invalid result type."
        )
    try:
        input_data = result["data"]
        warnings = result["warnings"]
        artifacts = result["artifacts"]
    except KeyError as error:
        raise InternalToolkitError(
            "The input validator returned an incomplete result.",
            details={"missing_key": str(error)},
            cause=error,
        ) from error
    bundle = _write_validation_bundle(
        output_dir,
        input_data=input_data,
        warnings=warnings,
        consumed_artifacts=artifacts,
    )
    return {"validation_result": result, "bundle": bundle}


__all__ = [
    "BUNDLE_SCHEMA_VERSION",
    "canonical_json_equal",
    "compute_validation_plan_id",
    "validate_request_to_bundle",
]
