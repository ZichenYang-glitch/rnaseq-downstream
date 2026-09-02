"""Shared, dependency-free helpers for locked P0 certification benchmarks."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable

SCHEMA_VERSION = "rnaseq-downstream-benchmark-report-v1"
PROJECT_ROOT = Path(__file__).resolve().parents[2]


class BenchmarkError(RuntimeError):
    """A certification benchmark could not produce a valid passing report."""


def canonical_json_bytes(value: object) -> bytes:
    """Serialize deterministic, strict JSON suitable for archived evidence."""

    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")


def write_json(path: Path, value: object) -> None:
    """Atomically write strict canonical JSON to *path*."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_bytes(canonical_json_bytes(value))
    os.replace(temporary, path)


def read_json(path: Path) -> dict[str, Any]:
    """Read one JSON object while rejecting duplicate keys and non-finite values."""

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
        raise BenchmarkError(f"Could not read benchmark JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise BenchmarkError(f"Benchmark JSON must be an object: {path}")
    return value


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_evidence(path: Path, *, name: str | None = None) -> dict[str, object]:
    """Describe an input or output file without embedding its machine path."""

    return {
        "name": name or path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def rscript_from_runtime(explicit: str | None = None) -> Path:
    """Resolve Rscript from an explicit path or the active locked Python prefix."""

    candidate = (
        Path(explicit) if explicit else Path(sys.executable).with_name("Rscript")
    )
    if not candidate.is_file():
        raise BenchmarkError(
            "Rscript is not available in the active Python prefix; run this benchmark "
            "with the locked P0 Conda interpreter or pass --rscript."
        )
    return candidate.resolve()


def run_rscript(
    rscript: Path,
    script: Path,
    arguments: Iterable[object],
    *,
    r_library: Path,
    timeout_seconds: int = 300,
) -> subprocess.CompletedProcess[str]:
    """Run one locked R script without inheriting a user R library."""

    if not r_library.is_dir():
        raise BenchmarkError(f"Locked R library does not exist: {r_library}")
    environment = os.environ.copy()
    environment["R_LIBS_USER"] = str(r_library.resolve())
    environment.pop("R_PROFILE_USER", None)
    environment.pop("R_ENVIRON_USER", None)
    command = [
        str(rscript),
        "--vanilla",
        str(script),
        str(r_library.resolve()),
        *(str(argument) for argument in arguments),
    ]
    completed = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        env=environment,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    if completed.returncode != 0:
        raise BenchmarkError(
            f"R benchmark step failed ({script.name}, exit {completed.returncode}).\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed


def read_tsv(path: Path, *, key: str = "gene_id") -> dict[str, dict[str, str]]:
    """Read a strict, rectangular TSV keyed by a unique identifier column."""

    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
            if reader.fieldnames is None or key not in reader.fieldnames:
                raise BenchmarkError(f"TSV {path} is missing required column {key!r}")
            rows: dict[str, dict[str, str]] = {}
            for line_number, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise BenchmarkError(
                        f"Non-rectangular TSV row {line_number}: {path}"
                    )
                identifier = row[key]
                if not identifier or identifier in rows:
                    raise BenchmarkError(
                        f"Blank or duplicate {key!r} at TSV row {line_number}: {path}"
                    )
                rows[identifier] = row
    except (OSError, UnicodeError, csv.Error) as exc:
        raise BenchmarkError(f"Could not read benchmark TSV {path}: {exc}") from exc
    return rows


def read_tsv_records(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read all records from a strict rectangular TSV without choosing a key."""

    fieldnames: list[str]
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
            if reader.fieldnames is None or any(not name for name in reader.fieldnames):
                raise BenchmarkError(f"TSV {path} has an invalid header")
            if len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise BenchmarkError(f"TSV {path} has duplicate header columns")
            fieldnames = list(reader.fieldnames)
            rows: list[dict[str, str]] = []
            for line_number, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise BenchmarkError(
                        f"Non-rectangular TSV row {line_number}: {path}"
                    )
                rows.append(row)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise BenchmarkError(f"Could not read benchmark TSV {path}: {exc}") from exc
    return fieldnames, rows


def finite_float(value: str, *, field: str, gene_id: str) -> float:
    """Parse a finite floating-point benchmark result."""

    try:
        parsed = float(value)
    except ValueError as exc:
        raise BenchmarkError(
            f"Invalid numeric {field!r} for gene {gene_id!r}: {value!r}"
        ) from exc
    if not math.isfinite(parsed):
        raise BenchmarkError(f"Non-finite {field!r} for gene {gene_id!r}")
    return parsed


def assert_report_shape(report: dict[str, Any], *, benchmark_id: str) -> None:
    """Validate dependency-free benchmark report schema v1 invariants."""

    expected_keys = {
        "schema_version",
        "benchmark_id",
        "status",
        "runtime",
        "implementation",
        "inputs",
        "artifacts",
        "thresholds",
        "metrics",
        "assertions",
    }
    if set(report) != expected_keys:
        raise BenchmarkError(
            "Benchmark report has an unexpected root shape: "
            f"missing={sorted(expected_keys - report.keys())}, "
            f"unexpected={sorted(report.keys() - expected_keys)}"
        )
    if report["schema_version"] != SCHEMA_VERSION:
        raise BenchmarkError(
            f"Unexpected benchmark schema: {report['schema_version']!r}"
        )
    if report["benchmark_id"] != benchmark_id:
        raise BenchmarkError(f"Unexpected benchmark id: {report['benchmark_id']!r}")
    if report["status"] not in {"pass", "fail"}:
        raise BenchmarkError(f"Unexpected benchmark status: {report['status']!r}")
    for field in ("runtime", "thresholds", "metrics"):
        if not isinstance(report[field], dict):
            raise BenchmarkError(f"Benchmark report field {field!r} must be an object")
    if not report["runtime"]:
        raise BenchmarkError("Benchmark report runtime must be non-empty")
    for field in ("implementation", "inputs"):
        if not isinstance(report[field], list) or not report[field]:
            raise BenchmarkError(f"Benchmark report {field} must be a non-empty array")
    if not isinstance(report["artifacts"], list):
        raise BenchmarkError("Benchmark report artifacts must be an array")
    for field in ("implementation", "inputs", "artifacts"):
        observed_names: set[str] = set()
        for index, record in enumerate(report[field]):
            if not isinstance(record, dict) or set(record) != {
                "name",
                "sha256",
                "size_bytes",
            }:
                raise BenchmarkError(
                    f"Benchmark report {field}[{index}] has an unexpected shape"
                )
            name = record["name"]
            digest = record["sha256"]
            size = record["size_bytes"]
            if not isinstance(name, str) or not name:
                raise BenchmarkError(
                    f"Benchmark report {field}[{index}].name must be non-empty"
                )
            if name in observed_names:
                raise BenchmarkError(
                    f"Benchmark report {field} contains duplicate name {name!r}"
                )
            observed_names.add(name)
            if (
                not isinstance(digest, str)
                or len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
            ):
                raise BenchmarkError(
                    f"Benchmark report {field}[{index}].sha256 is not canonical"
                )
            if isinstance(size, bool) or not isinstance(size, int) or size < 0:
                raise BenchmarkError(
                    f"Benchmark report {field}[{index}].size_bytes is invalid"
                )
    if not isinstance(report["assertions"], list) or not report["assertions"]:
        raise BenchmarkError("Benchmark report assertions must be a non-empty array")
    for assertion in report["assertions"]:
        if not isinstance(assertion, dict):
            raise BenchmarkError("Every benchmark assertion must be an object")
        if set(assertion) != {"name", "passed", "observed", "requirement"}:
            raise BenchmarkError("Benchmark assertion has an unexpected shape")
        if not isinstance(assertion["name"], str) or not assertion["name"]:
            raise BenchmarkError("Benchmark assertion name must be non-empty")
        if (
            not isinstance(assertion["requirement"], str)
            or not assertion["requirement"]
        ):
            raise BenchmarkError("Benchmark assertion requirement must be non-empty")
        if not isinstance(assertion["passed"], bool):
            raise BenchmarkError("Benchmark assertion passed must be boolean")
    expected_status = (
        "pass"
        if all(assertion["passed"] for assertion in report["assertions"])
        else "fail"
    )
    if report["status"] != expected_status:
        raise BenchmarkError(
            "Benchmark report status disagrees with its assertion outcomes"
        )
