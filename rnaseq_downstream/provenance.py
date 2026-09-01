"""Deterministic, read-only provenance helpers for declared input files."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping
import unicodedata

from .errors import InputIntegrityError, InputReadError, InvalidRequestError

_SHA256_PATTERN = re.compile(r"^[0-9a-fA-F]{64}$")
_MAX_EXACT_JSON_INTEGER = (2**53) - 1


class _DuplicateJsonKeyError(ValueError):
    def __init__(self, key: str) -> None:
        super().__init__(key)
        self.key = key


class _NonStandardJsonConstantError(ValueError):
    def __init__(self, value: str) -> None:
        super().__init__(value)
        self.value = value


class _OversizedJsonIntegerError(ValueError):
    def __init__(self, literal: str) -> None:
        super().__init__(literal)
        self.literal = literal


def _strict_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    document: dict[str, Any] = {}
    for key, value in pairs:
        if key in document:
            raise _DuplicateJsonKeyError(key)
        document[key] = value
    return document


def _reject_json_constant(value: str) -> None:
    raise _NonStandardJsonConstantError(value)


def _strict_json_integer(literal: str) -> int:
    unsigned = literal.removeprefix("-")
    normalized = unsigned.lstrip("0") or "0"
    maximum = str(_MAX_EXACT_JSON_INTEGER)
    if len(normalized) > len(maximum) or (
        len(normalized) == len(maximum) and normalized > maximum
    ):
        raise _OversizedJsonIntegerError(literal)
    return int(literal)


@dataclass(frozen=True)
class FileFingerprint:
    """A stable description of one consumed regular file."""

    role: str
    declared_path: str
    path: str
    sha256: str
    size_bytes: int

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable artifact record."""

        return {
            "kind": "consumed_input",
            "role": self.role,
            "declared_path": self.declared_path,
            "path": self.path,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
        }


@dataclass(frozen=True)
class _CapturedFile:
    """Immutable bytes used for both provenance and subsequent parsing."""

    content: bytes
    sha256: str
    size_bytes: int


def resolve_declared_path(declared_path: str, *, relative_to: Path) -> Path:
    """Resolve a JSON path relative to the JSON document that declares it."""

    if not isinstance(declared_path, str) or not declared_path.strip():
        raise InvalidRequestError(
            "A declared path must be a non-empty string.",
            details={"declared_path": declared_path},
        )
    if has_control_characters(declared_path):
        raise InvalidRequestError(
            "A declared path contains a Unicode control character.",
            details={"declared_path": declared_path},
        )
    candidate = Path(declared_path)
    if not candidate.is_absolute():
        candidate = relative_to.parent / candidate
    try:
        return candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise InputReadError(
            "A declared input path could not be resolved.",
            path=candidate,
            operation="resolve_input_path",
            cause=error,
            details={"declared_path": declared_path},
        ) from error


def require_directory(path: Path, *, declared_path: str, role: str) -> Path:
    """Require a resolved path to be a directory."""

    if not path.is_dir():
        error = NotADirectoryError(str(path))
        raise InputReadError(
            "A declared input directory is not a directory.",
            path=path,
            operation="read_input_directory",
            cause=error,
            details={"declared_path": declared_path, "role": role},
        )
    return path


def _capture_file(path: Path) -> _CapturedFile:
    """Capture a regular file once so the digest and parser see identical bytes."""

    if not path.is_file():
        error = (
            IsADirectoryError(str(path))
            if path.is_dir()
            else FileNotFoundError(str(path))
        )
        raise InputReadError(
            "A declared input file is not a regular file.",
            path=path,
            operation="hash_input_file",
            cause=error,
        )

    try:
        before = path.stat()
        with path.open("rb") as handle:
            content = handle.read()
        after = path.stat()
    except OSError as error:
        raise InputReadError(
            "A declared input file could not be hashed.",
            path=path,
            operation="hash_input_file",
            cause=error,
        ) from error

    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if identity_before != identity_after:
        raise InputIntegrityError(
            "An input file changed while its digest was being calculated.",
            details={"path": str(path)},
        )
    digest = hashlib.sha256(content).hexdigest()
    if len(content) != after.st_size:
        raise InputIntegrityError(
            "An input file size did not match its captured byte snapshot.",
            details={
                "path": str(path),
                "stat_size_bytes": after.st_size,
                "captured_size_bytes": len(content),
            },
        )
    return _CapturedFile(
        content=content,
        sha256=digest,
        size_bytes=len(content),
    )


def _fingerprint(path: Path) -> tuple[str, int]:
    if not path.is_file():
        error = (
            IsADirectoryError(str(path))
            if path.is_dir()
            else FileNotFoundError(str(path))
        )
        raise InputReadError(
            "A declared input file is not a regular file.",
            path=path,
            operation="hash_input_file",
            cause=error,
        )
    digest = hashlib.sha256()
    try:
        before = path.stat()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        after = path.stat()
    except OSError as error:
        raise InputReadError(
            "A declared input file could not be hashed.",
            path=path,
            operation="hash_input_file",
            cause=error,
        ) from error
    identity_before = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    identity_after = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if identity_before != identity_after:
        raise InputIntegrityError(
            "An input file changed while its digest was being calculated.",
            details={"path": str(path)},
        )
    return digest.hexdigest(), after.st_size


class ArtifactInventory:
    """Collect file fingerprints and emit them in deterministic order."""

    def __init__(self) -> None:
        self._fingerprints: list[FileFingerprint] = []
        self._digest_cache: dict[Path, tuple[str, int]] = {}
        self._snapshot_cache: dict[Path, bytes] = {}

    def add(
        self,
        path: Path,
        *,
        role: str,
        declared_path: str,
        expected_sha256: str | None = None,
        capture_for_parse: bool = False,
    ) -> FileFingerprint:
        """Fingerprint a file, optionally enforcing a predeclared digest."""

        if expected_sha256 is not None:
            if not isinstance(expected_sha256, str) or not _SHA256_PATTERN.fullmatch(
                expected_sha256
            ):
                raise InvalidRequestError(
                    "A declared SHA-256 digest must contain exactly 64 hexadecimal characters.",
                    details={"role": role, "sha256": expected_sha256},
                )
            expected_sha256 = expected_sha256.lower()

        try:
            resolved = path.resolve(strict=True)
        except (OSError, RuntimeError) as error:
            raise InputReadError(
                "A declared input file disappeared before it could be fingerprinted.",
                path=path,
                operation="resolve_input_for_hash",
                cause=error,
                details={"role": role},
            ) from error
        cached = self._digest_cache.get(resolved)
        if capture_for_parse:
            snapshot = _capture_file(resolved)
            if cached is not None and cached != (
                snapshot.sha256,
                snapshot.size_bytes,
            ):
                raise InputIntegrityError(
                    "An input file changed before its parse snapshot was captured.",
                    details={
                        "path": str(resolved),
                        "expected_sha256": cached[0],
                        "observed_sha256": snapshot.sha256,
                    },
                )
            digest, size = snapshot.sha256, snapshot.size_bytes
            self._snapshot_cache[resolved] = snapshot.content
        elif cached is None:
            digest, size = _fingerprint(resolved)
        else:
            digest, size = cached
        self._digest_cache[resolved] = (digest, size)

        if expected_sha256 is not None and digest != expected_sha256:
            raise InputIntegrityError(
                "An input file does not match its declared SHA-256 digest.",
                details={
                    "role": role,
                    "path": str(resolved),
                    "expected_sha256": expected_sha256,
                    "observed_sha256": digest,
                },
            )

        fingerprint = FileFingerprint(
            role=role,
            declared_path=declared_path,
            path=str(resolved),
            sha256=digest,
            size_bytes=size,
        )
        self._fingerprints.append(fingerprint)
        return fingerprint

    def consume_bytes(self, path: Path) -> bytes:
        """Transfer and release the immutable bytes captured for one parser."""

        content = self._snapshot_cache.pop(path, None)
        if content is None:
            raise InputIntegrityError(
                "A parser requested bytes that were not captured for provenance.",
                details={"path": str(path)},
            )
        return content

    def verify_unchanged(self) -> None:
        """Fail if any consumed path no longer contains the fingerprinted bytes."""

        if self._snapshot_cache:
            raise InputIntegrityError(
                "One or more captured parse snapshots were not released.",
                details={"paths": sorted(str(path) for path in self._snapshot_cache)},
            )
        for path, (expected_digest, expected_size) in sorted(
            self._digest_cache.items(), key=lambda item: str(item[0])
        ):
            observed_digest, observed_size = _fingerprint(path)
            if observed_digest != expected_digest or observed_size != expected_size:
                raise InputIntegrityError(
                    "An input file changed between provenance capture and parsing.",
                    details={
                        "path": str(path),
                        "expected_sha256": expected_digest,
                        "observed_sha256": observed_digest,
                        "expected_size_bytes": expected_size,
                        "observed_size_bytes": observed_size,
                    },
                )

    def records(self) -> list[dict[str, Any]]:
        """Return deterministic JSON artifact records."""

        return [
            item.to_dict()
            for item in sorted(
                self._fingerprints,
                key=lambda item: (item.role, item.path, item.declared_path),
            )
        ]


def read_json_object(
    path: Path,
    *,
    document_role: str,
    content: bytes | None = None,
) -> dict[str, Any]:
    """Read a UTF-8 JSON object and map parse failures to public errors."""

    try:
        text = (
            path.read_text(encoding="utf-8")
            if content is None
            else content.decode("utf-8")
        )
        document = json.loads(
            text,
            object_pairs_hook=_strict_json_object,
            parse_constant=_reject_json_constant,
            parse_int=_strict_json_integer,
        )
    except _DuplicateJsonKeyError as error:
        raise InvalidRequestError(
            "A JSON input document contains a duplicate object key.",
            details={
                "path": str(path),
                "document_role": document_role,
                "duplicate_key": error.key,
            },
            cause=error,
        ) from error
    except _NonStandardJsonConstantError as error:
        raise InvalidRequestError(
            "A JSON input document contains a non-standard numeric constant.",
            details={
                "path": str(path),
                "document_role": document_role,
                "value": error.value,
            },
            cause=error,
        ) from error
    except _OversizedJsonIntegerError as error:
        raise InvalidRequestError(
            "A JSON integer exceeds the exact interoperable numeric range.",
            details={
                "path": str(path),
                "document_role": document_role,
                "maximum_absolute_integer": _MAX_EXACT_JSON_INTEGER,
                "digit_count": len(error.literal.removeprefix("-").lstrip("0") or "0"),
            },
            cause=error,
        ) from error
    except json.JSONDecodeError as error:
        raise InvalidRequestError(
            "A JSON input document is not valid JSON.",
            details={
                "path": str(path),
                "document_role": document_role,
                "line": error.lineno,
                "column": error.colno,
            },
            cause=error,
        ) from error
    except (OSError, UnicodeError) as error:
        raise InputReadError(
            "A JSON input document could not be read as UTF-8.",
            path=path,
            operation="read_json",
            cause=error,
            details={"document_role": document_role},
        ) from error
    if not isinstance(document, dict):
        raise InvalidRequestError(
            "A JSON input document must contain an object at its root.",
            details={"path": str(path), "document_role": document_role},
        )
    return document


def require_expected_keys(
    document: Mapping[str, Any],
    *,
    allowed: set[str],
    required: set[str],
    context: str,
) -> None:
    """Reject missing and unknown JSON keys rather than guessing intent."""

    missing = sorted(required - document.keys())
    unknown = sorted(document.keys() - allowed)
    if missing or unknown:
        raise InvalidRequestError(
            f"The {context} object does not satisfy the input schema.",
            details={
                "context": context,
                "missing_keys": missing,
                "unknown_keys": unknown,
            },
        )


def require_nonempty_string(value: Any, *, field: str) -> str:
    """Return an exact required string without coercion or normalization."""

    if not isinstance(value, str) or not value.strip():
        raise InvalidRequestError(
            f"'{field}' must be a non-empty string.",
            details={"field": field},
        )
    if value != value.strip():
        raise InvalidRequestError(
            f"'{field}' must not contain surrounding whitespace.",
            details={"field": field, "value": value},
        )
    if has_control_characters(value):
        raise InvalidRequestError(
            f"'{field}' must not contain Unicode control characters.",
            details={"field": field},
        )
    return value


def has_control_characters(value: str) -> bool:
    """Return true for C0, C1, or any other Unicode ``Cc`` character."""

    return any(unicodedata.category(character) == "Cc" for character in value)


__all__ = [
    "ArtifactInventory",
    "FileFingerprint",
    "has_control_characters",
    "read_json_object",
    "require_directory",
    "require_expected_keys",
    "require_nonempty_string",
    "resolve_declared_path",
]
