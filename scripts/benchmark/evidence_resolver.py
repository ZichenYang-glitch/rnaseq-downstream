"""Resolve archived benchmark evidence against immutable environment snapshots.

This module is intentionally separate from :mod:`common`: the four historical
P0/C1/C2 reports bind ``common.py`` by digest, so changing that file would make
the evidence-evolution mechanism circular.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path, PurePosixPath
import re
from typing import Any, Mapping

try:  # Package import in tests; direct-script import in benchmark commands.
    from .common import BenchmarkError, PROJECT_ROOT, sha256_file
except ImportError:  # pragma: no cover - exercised by standalone scripts
    from common import BenchmarkError, PROJECT_ROOT, sha256_file


SNAPSHOT_SCHEMA_VERSION = "rnaseq-downstream-environment-snapshot-v1"
SNAPSHOT_ROOT = PROJECT_ROOT / "environment" / "snapshots"
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")


@dataclass(frozen=True)
class SnapshotFile:
    """One verified file in an append-only environment snapshot."""

    snapshot_id: str
    source_revision: str
    source_path: str
    snapshot_path: str
    path: Path
    sha256: str
    size_bytes: int


def _strict_json(path: Path) -> dict[str, Any]:
    def reject_constant(value: str) -> None:
        raise BenchmarkError(f"Non-finite JSON constant in {path}: {value}")

    def unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise BenchmarkError(f"Duplicate JSON key in {path}: {key}")
            result[key] = value
        return result

    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=unique_pairs,
            parse_constant=reject_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise BenchmarkError(f"Could not read snapshot manifest {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise BenchmarkError(f"Snapshot manifest must be an object: {path}")
    return value


def _portable_relative_path(value: object, *, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise BenchmarkError(f"{context} must be a non-empty relative POSIX path")
    parsed = PurePosixPath(value)
    if parsed.is_absolute() or parsed.as_posix() != value or ".." in parsed.parts:
        raise BenchmarkError(
            f"{context} is not a canonical relative POSIX path: {value!r}"
        )
    return value


def _digest(value: object, *, context: str) -> str:
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise BenchmarkError(f"{context} is not a canonical SHA-256 digest")
    return value


def _size(value: object, *, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise BenchmarkError(f"{context} is not a non-negative integer")
    return value


def _snapshot_files(snapshot_directory: Path) -> tuple[SnapshotFile, ...]:
    if snapshot_directory.is_symlink() or not snapshot_directory.is_dir():
        raise BenchmarkError(
            f"Environment snapshot is not a real directory: {snapshot_directory}"
        )
    manifest_path = snapshot_directory / "manifest.json"
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise BenchmarkError(
            f"Environment snapshot manifest is missing: {manifest_path}"
        )
    manifest = _strict_json(manifest_path)
    expected_keys = {
        "schema_version",
        "snapshot_id",
        "source_revision",
        "immutability",
        "files",
    }
    if set(manifest) != expected_keys:
        raise BenchmarkError(
            f"Snapshot manifest has an unexpected root shape: {manifest_path}"
        )
    if manifest["schema_version"] != SNAPSHOT_SCHEMA_VERSION:
        raise BenchmarkError(
            f"Snapshot manifest schema is unsupported: {manifest_path}"
        )
    snapshot_id = manifest["snapshot_id"]
    if not isinstance(snapshot_id, str) or snapshot_id != snapshot_directory.name:
        raise BenchmarkError(
            f"Snapshot id does not match its directory: {manifest_path}"
        )
    source_revision = manifest["source_revision"]
    if (
        not isinstance(source_revision, str)
        or _REVISION_PATTERN.fullmatch(source_revision) is None
    ):
        raise BenchmarkError(f"Snapshot source revision is invalid: {manifest_path}")
    if manifest["immutability"] != "append_only":
        raise BenchmarkError(f"Snapshot is not marked append-only: {manifest_path}")
    file_records = manifest["files"]
    if not isinstance(file_records, list) or not file_records:
        raise BenchmarkError(
            f"Snapshot files must be a non-empty array: {manifest_path}"
        )

    observed_sources: set[str] = set()
    observed_snapshots: set[str] = set()
    verified: list[SnapshotFile] = []
    for index, record in enumerate(file_records):
        if not isinstance(record, dict) or set(record) != {
            "source_path",
            "snapshot_path",
            "sha256",
            "size_bytes",
        }:
            raise BenchmarkError(
                f"Snapshot file record {index} has an unexpected shape: {manifest_path}"
            )
        source_path = _portable_relative_path(
            record["source_path"], context=f"snapshot files[{index}].source_path"
        )
        snapshot_path = _portable_relative_path(
            record["snapshot_path"], context=f"snapshot files[{index}].snapshot_path"
        )
        if snapshot_path == "manifest.json":
            raise BenchmarkError("A snapshot manifest cannot list itself as a payload")
        if source_path in observed_sources or snapshot_path in observed_snapshots:
            raise BenchmarkError(
                f"Snapshot manifest contains a duplicate file: {manifest_path}"
            )
        observed_sources.add(source_path)
        observed_snapshots.add(snapshot_path)
        digest = _digest(record["sha256"], context=f"snapshot files[{index}].sha256")
        size = _size(
            record["size_bytes"], context=f"snapshot files[{index}].size_bytes"
        )
        payload_path = snapshot_directory.joinpath(*PurePosixPath(snapshot_path).parts)
        if payload_path.is_symlink() or not payload_path.is_file():
            raise BenchmarkError(f"Snapshot payload is missing: {payload_path}")
        if payload_path.stat().st_size != size or sha256_file(payload_path) != digest:
            raise BenchmarkError(
                f"Snapshot payload disagrees with its manifest: {payload_path}"
            )
        verified.append(
            SnapshotFile(
                snapshot_id=snapshot_id,
                source_revision=source_revision,
                source_path=source_path,
                snapshot_path=snapshot_path,
                path=payload_path,
                sha256=digest,
                size_bytes=size,
            )
        )

    actual_files: set[str] = set()
    for candidate in snapshot_directory.rglob("*"):
        if candidate.is_symlink():
            raise BenchmarkError(f"Snapshot content cannot be a symlink: {candidate}")
        if candidate.is_file():
            actual_files.add(candidate.relative_to(snapshot_directory).as_posix())
    expected_files = {"manifest.json", *observed_snapshots}
    if actual_files != expected_files:
        raise BenchmarkError(
            "Snapshot file inventory disagrees with its manifest: "
            f"missing={sorted(expected_files - actual_files)}, "
            f"unexpected={sorted(actual_files - expected_files)}"
        )
    return tuple(verified)


def load_environment_snapshots(
    snapshot_root: Path = SNAPSHOT_ROOT,
) -> tuple[SnapshotFile, ...]:
    """Load and verify every registered append-only environment snapshot."""

    if snapshot_root.is_symlink() or not snapshot_root.is_dir():
        raise BenchmarkError(
            f"Environment snapshot root is unavailable: {snapshot_root}"
        )
    children = sorted(snapshot_root.iterdir(), key=lambda path: path.name)
    if not children:
        raise BenchmarkError(f"Environment snapshot root is empty: {snapshot_root}")
    unexpected = [path for path in children if path.is_symlink() or not path.is_dir()]
    if unexpected:
        raise BenchmarkError(
            f"Environment snapshot root contains unexpected entries: {unexpected}"
        )
    snapshots = tuple(
        record for directory in children for record in _snapshot_files(directory)
    )
    identities = [
        (record.snapshot_id, record.source_path, record.sha256, record.size_bytes)
        for record in snapshots
    ]
    if len(identities) != len(set(identities)):
        raise BenchmarkError(
            "Environment snapshots contain duplicate evidence identities"
        )
    return snapshots


def resolve_archived_implementation_path(
    current_path: Path,
    *,
    expected_sha256: str,
    expected_size: int,
    project_root: Path = PROJECT_ROOT,
    snapshot_root: Path = SNAPSHOT_ROOT,
) -> Path:
    """Resolve current implementation evidence, falling back by hash to a snapshot."""

    digest = _digest(expected_sha256, context="archived implementation sha256")
    size = _size(expected_size, context="archived implementation size_bytes")
    root = project_root.resolve(strict=True)
    lexical_path = current_path if current_path.is_absolute() else root / current_path
    if lexical_path.is_symlink():
        raise BenchmarkError(
            f"Implementation evidence cannot be a symlink: {lexical_path}"
        )
    resolved_current = lexical_path.resolve(strict=False)
    try:
        source_path = resolved_current.relative_to(root).as_posix()
    except ValueError as exc:
        raise BenchmarkError(
            f"Implementation evidence is outside the project root: {current_path}"
        ) from exc
    if resolved_current.is_file():
        current_size = resolved_current.stat().st_size
        current_digest = sha256_file(resolved_current)
        if current_size == size and current_digest == digest:
            return resolved_current

    matches = [
        record
        for record in load_environment_snapshots(snapshot_root)
        if record.source_path == source_path
        and record.sha256 == digest
        and record.size_bytes == size
    ]
    if not matches:
        raise BenchmarkError(
            "No current or immutable snapshot file matches archived implementation "
            f"evidence for {source_path!r} ({digest}, {size} bytes)"
        )
    return matches[0].path


def verify_archived_implementation(
    recorded_items: object,
    current_paths: Mapping[str, Path],
    *,
    project_root: Path = PROJECT_ROOT,
    snapshot_root: Path = SNAPSHOT_ROOT,
) -> dict[str, Path]:
    """Verify a report's full implementation inventory and return resolved paths."""

    if not isinstance(recorded_items, list):
        raise BenchmarkError("Archived implementation evidence must be an array")
    recorded: dict[str, Mapping[str, object]] = {}
    for index, item in enumerate(recorded_items):
        if not isinstance(item, Mapping) or set(item) != {
            "name",
            "sha256",
            "size_bytes",
        }:
            raise BenchmarkError(
                f"Archived implementation record {index} has an unexpected shape"
            )
        name = item["name"]
        if not isinstance(name, str) or not name or name in recorded:
            raise BenchmarkError(
                f"Archived implementation record {index} has an invalid name"
            )
        recorded[name] = item
    if set(recorded) != set(current_paths):
        raise BenchmarkError(
            "Archived implementation inventory names disagree with the expected paths"
        )

    resolved: dict[str, Path] = {}
    for name, current_path in current_paths.items():
        item = recorded[name]
        resolved[name] = resolve_archived_implementation_path(
            current_path,
            expected_sha256=item["sha256"],  # type: ignore[arg-type]
            expected_size=item["size_bytes"],  # type: ignore[arg-type]
            project_root=project_root,
            snapshot_root=snapshot_root,
        )
    return resolved


def manifest_evidence(
    snapshot_id: str,
    snapshot_root: Path = SNAPSHOT_ROOT,
) -> dict[str, object]:
    """Return hash evidence for one fully verified snapshot manifest."""

    snapshots = load_environment_snapshots(snapshot_root)
    matching = {
        record.snapshot_id for record in snapshots if record.snapshot_id == snapshot_id
    }
    if matching != {snapshot_id}:
        raise BenchmarkError(f"Environment snapshot is unavailable: {snapshot_id}")
    manifest = snapshot_root / snapshot_id / "manifest.json"
    return {
        "name": f"environment/snapshots/{snapshot_id}/manifest.json",
        "sha256": sha256_file(manifest),
        "size_bytes": manifest.stat().st_size,
    }
