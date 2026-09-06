#!/usr/bin/env python3
"""Archive proof that frozen P0/C1/C2 gates survive an environment expansion."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
from pathlib import Path
from typing import Any, Mapping

try:  # Package import in tests; direct-script import on the command line.
    from .common import (
        BenchmarkError,
        PROJECT_ROOT,
        assert_report_shape,
        canonical_json_bytes,
        file_evidence,
        read_json,
        rscript_from_runtime,
        sha256_file,
        write_json,
    )
    from .evidence_resolver import (
        SNAPSHOT_ROOT,
        load_environment_snapshots,
        manifest_evidence,
        resolve_archived_implementation_path,
    )
except ImportError:  # pragma: no cover - exercised by standalone invocation
    from common import (
        BenchmarkError,
        PROJECT_ROOT,
        assert_report_shape,
        canonical_json_bytes,
        file_evidence,
        read_json,
        rscript_from_runtime,
        sha256_file,
        write_json,
    )
    from evidence_resolver import (
        SNAPSHOT_ROOT,
        load_environment_snapshots,
        manifest_evidence,
        resolve_archived_implementation_path,
    )


SCHEMA_VERSION = "rnaseq-downstream-environment-compatibility-report-v1"
COMPATIBILITY_ID = "p1-deseq2-environment-expansion-v1"
BASELINE_SNAPSHOT_ID = "p0-c1-c2-c6b6cd9"
SOURCE_REVISION = "c6b6cd970c1a0c6807145dd01033505efae60215"
ENVIRONMENT_IMPLEMENTATION_PATHS = {
    "renv.lock": PROJECT_ROOT / "renv.lock",
    "r-sources.lock": PROJECT_ROOT / "environment/r-sources.lock",
    "verify.R": PROJECT_ROOT / "environment/verify.R",
}
CONDA_LOCK_PATH = PROJECT_ROOT / "conda-lock.yml"
RSCRIPT_RUNTIME_FIELD = "rscript_sha256"
RSCRIPT_RELOCATION_REASON = "locked_conda_prefix_relocation"
BENCHMARKS = (
    {
        "benchmark_id": "airway-edger-ql-same-engine-v1",
        "filename": "airway-benchmark-report.json",
        "archived_path": PROJECT_ROOT / "tests/oracle/airway-benchmark-report.json",
        "archived_sha256": "446955638d367dc85b3b6440ddce7b5299717417aae2e0394991681f5df4c23f",
        "archived_size": 6734,
        "rebuild_metadata_input": None,
        "source_package": None,
    },
    {
        "benchmark_id": "airway-limma-gene-set-same-engine-v1",
        "filename": "pathway-airway-benchmark-report.json",
        "archived_path": PROJECT_ROOT
        / "tests/oracle/pathway-airway-benchmark-report.json",
        "archived_sha256": "c14c1db1c94f53abee7c5d116b9b627d75041f5a8eb4edf7caafa5e9a059ca62",
        "archived_size": 8236,
        "rebuild_metadata_input": "airway/DESCRIPTION",
        "source_package": "airway",
    },
    {
        "benchmark_id": "compcoder-edger-ql-nb-fdr-tpr-v1",
        "filename": "compcoder-benchmark-report.json",
        "archived_path": PROJECT_ROOT
        / "tests/simulation/compcoder-benchmark-report.json",
        "archived_sha256": "a95f50ecbc96a48847d15d2587b421603222f1f975cc4bfdbc123e11f16ecbce",
        "archived_size": 16624,
        "rebuild_metadata_input": None,
        "source_package": None,
    },
    {
        "benchmark_id": "compcoder-limma-self-contained-fdr-tpr-v1",
        "filename": "pathway-compcoder-benchmark-report.json",
        "archived_path": PROJECT_ROOT
        / "tests/simulation/pathway-compcoder-benchmark-report.json",
        "archived_sha256": "0f34458fa1fe43d31e3fcf5b6909362043229c1ea98df14e2d2015f696f6bc52",
        "archived_size": 101207,
        "rebuild_metadata_input": "compcodeR/DESCRIPTION",
        "source_package": "compcodeR",
    },
)


def _implementation_map(
    report: Mapping[str, Any],
    *,
    context: str,
) -> dict[str, dict[str, object]]:
    items = report.get("implementation")
    if not isinstance(items, list):
        raise BenchmarkError(f"{context} implementation must be an array")
    result: dict[str, dict[str, object]] = {}
    for index, item in enumerate(items):
        if not isinstance(item, dict) or set(item) != {
            "name",
            "sha256",
            "size_bytes",
        }:
            raise BenchmarkError(
                f"{context} implementation[{index}] has an unexpected shape"
            )
        name = item["name"]
        if not isinstance(name, str) or not name or name in result:
            raise BenchmarkError(
                f"{context} implementation[{index}] has an invalid name"
            )
        result[name] = item
    return result


def _source_manifest_records(path: Path) -> dict[str, dict[str, str]]:
    expected_fields = [
        "package",
        "version",
        "repository",
        "role",
        "url",
        "sha256",
    ]
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
            if reader.fieldnames != expected_fields:
                raise BenchmarkError(f"R source manifest header is invalid: {path}")
            records: dict[str, dict[str, str]] = {}
            for row in reader:
                if None in row or any(
                    value is None or value == "" for value in row.values()
                ):
                    raise BenchmarkError(f"R source manifest row is invalid: {path}")
                package = row["package"]
                if package in records:
                    raise BenchmarkError(
                        f"R source manifest package is duplicated: {path}/{package}"
                    )
                records[package] = dict(row)  # type: ignore[arg-type]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise BenchmarkError(f"Could not read R source manifest {path}: {exc}") from exc
    return records


def _unchanged_source_archives(
    *,
    expanded_source_lock: Path = ENVIRONMENT_IMPLEMENTATION_PATHS["r-sources.lock"],
    snapshot_root: Path = SNAPSHOT_ROOT,
) -> dict[str, dict[str, str]]:
    snapshot_records = load_environment_snapshots(snapshot_root)
    source_locks = [
        record.path
        for record in snapshot_records
        if record.snapshot_id == BASELINE_SNAPSHOT_ID
        and record.source_path == "environment/r-sources.lock"
    ]
    if len(source_locks) != 1:
        raise BenchmarkError("The baseline R source manifest snapshot is unavailable")
    baseline = _source_manifest_records(source_locks[0])
    current = _source_manifest_records(expanded_source_lock)
    packages = {
        package
        for specification in BENCHMARKS
        if isinstance((package := specification["source_package"]), str)
    }
    unchanged: dict[str, dict[str, str]] = {}
    for package in packages:
        if package not in baseline or package not in current:
            raise BenchmarkError(f"Pinned source package is missing: {package}")
        if baseline[package] != current[package]:
            raise BenchmarkError(
                f"Pinned source archive changed across environments: {package}"
            )
        unchanged[package] = current[package]
    return unchanged


def _assert_frozen_archive(specification: Mapping[str, Any]) -> dict[str, Any]:
    path = specification["archived_path"]
    if not isinstance(path, Path):
        raise BenchmarkError("Internal archived benchmark path is invalid")
    if (
        not path.is_file()
        or path.stat().st_size != specification["archived_size"]
        or sha256_file(path) != specification["archived_sha256"]
    ):
        raise BenchmarkError(f"Frozen archived benchmark report changed: {path}")
    report = read_json(path)
    assert_report_shape(report, benchmark_id=specification["benchmark_id"])
    if report["status"] != "pass":
        raise BenchmarkError(f"Frozen archived benchmark is not passing: {path}")
    return report


def _file_record(value: object, *, context: str) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != {"name", "sha256", "size_bytes"}:
        raise BenchmarkError(f"{context} has an unexpected file-record shape")
    name = value["name"]
    digest = value["sha256"]
    size = value["size_bytes"]
    if not isinstance(name, str) or not name:
        raise BenchmarkError(f"{context}.name is invalid")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        raise BenchmarkError(f"{context}.sha256 is invalid")
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        raise BenchmarkError(f"{context}.size_bytes is invalid")
    return value


def assert_compatibility_report_shape(
    report: Mapping[str, Any],
    *,
    environment_paths: Mapping[str, Path] = ENVIRONMENT_IMPLEMENTATION_PATHS,
    project_root: Path = PROJECT_ROOT,
    snapshot_root: Path = SNAPSHOT_ROOT,
) -> None:
    """Independently validate an archived environment-compatibility report."""

    expected_keys = {
        "schema_version",
        "compatibility_id",
        "status",
        "baseline_snapshot",
        "expanded_environment",
        "permitted_implementation_changes",
        "benchmarks",
        "assertions",
    }
    if set(report) != expected_keys:
        raise BenchmarkError("Compatibility report has an unexpected root shape")
    if report["schema_version"] != SCHEMA_VERSION:
        raise BenchmarkError("Compatibility report schema is unsupported")
    if report["compatibility_id"] != COMPATIBILITY_ID or report["status"] != "pass":
        raise BenchmarkError("Compatibility report identity or status is invalid")
    baseline = report["baseline_snapshot"]
    if not isinstance(baseline, dict) or set(baseline) != {
        "snapshot_id",
        "source_revision",
        "manifest",
    }:
        raise BenchmarkError("Compatibility baseline snapshot has an unexpected shape")
    if (
        baseline["snapshot_id"] != BASELINE_SNAPSHOT_ID
        or baseline["source_revision"] != SOURCE_REVISION
        or baseline["manifest"]
        != manifest_evidence(BASELINE_SNAPSHOT_ID, snapshot_root)
    ):
        raise BenchmarkError("Compatibility baseline snapshot identity is invalid")

    permitted = report["permitted_implementation_changes"]
    if (
        not isinstance(permitted, list)
        or len(permitted) != len(environment_paths)
        or set(permitted) != set(environment_paths)
    ):
        raise BenchmarkError("Compatibility permitted-change inventory is invalid")
    expanded_environment = report["expanded_environment"]
    if not isinstance(expanded_environment, dict) or set(expanded_environment) != {
        "files"
    }:
        raise BenchmarkError("Expanded environment evidence has an unexpected shape")
    current_files = expanded_environment["files"]
    if not isinstance(current_files, list):
        raise BenchmarkError("Expanded environment files must be an array")
    observed_current = {
        _file_record(item, context="expanded environment file")["name"]: item
        for item in current_files
    }
    if set(observed_current) != set(environment_paths):
        raise BenchmarkError("Expanded environment file inventory is invalid")
    resolved_environment: dict[str, Path] = {}
    for name, current_path in environment_paths.items():
        evidence = observed_current[name]
        resolved_environment[name] = resolve_archived_implementation_path(
            current_path,
            expected_sha256=evidence["sha256"],  # type: ignore[arg-type]
            expected_size=evidence["size_bytes"],  # type: ignore[arg-type]
            project_root=project_root,
            snapshot_root=snapshot_root,
        )
    source_archives = _unchanged_source_archives(
        expanded_source_lock=resolved_environment["r-sources.lock"],
        snapshot_root=snapshot_root,
    )

    benchmark_records = report["benchmarks"]
    if not isinstance(benchmark_records, list) or len(benchmark_records) != len(
        BENCHMARKS
    ):
        raise BenchmarkError(
            "Compatibility report must contain exactly four benchmarks"
        )
    observed_benchmarks: dict[str, dict[str, Any]] = {}
    for index, record in enumerate(benchmark_records):
        if not isinstance(record, dict) or set(record) != {
            "benchmark_id",
            "archived_report",
            "expanded_environment_report",
            "artifact_inventory",
            "implementation_changes",
            "runtime_metadata_changes",
            "rebuild_metadata_changes",
        }:
            raise BenchmarkError(
                f"Compatibility benchmark record {index} has an unexpected shape"
            )
        benchmark_id = record["benchmark_id"]
        if (
            not isinstance(benchmark_id, str)
            or not benchmark_id
            or benchmark_id in observed_benchmarks
        ):
            raise BenchmarkError(
                f"Compatibility benchmark record {index} has an invalid id"
            )
        observed_benchmarks[benchmark_id] = record

    observed_runtime_before: set[str] = set()
    observed_runtime_after: set[str] = set()
    for specification in BENCHMARKS:
        benchmark_id = specification["benchmark_id"]
        assert isinstance(benchmark_id, str)
        if benchmark_id not in observed_benchmarks:
            raise BenchmarkError(f"Compatibility benchmark is missing: {benchmark_id}")
        record = observed_benchmarks[benchmark_id]
        archived = _assert_frozen_archive(specification)
        archived_path = specification["archived_path"]
        filename = specification["filename"]
        rebuild_metadata_input = specification["rebuild_metadata_input"]
        source_package = specification["source_package"]
        assert isinstance(archived_path, Path)
        assert isinstance(filename, str)
        assert rebuild_metadata_input is None or isinstance(rebuild_metadata_input, str)
        assert source_package is None or isinstance(source_package, str)
        expected_archived_evidence = file_evidence(
            archived_path,
            name=archived_path.relative_to(PROJECT_ROOT).as_posix(),
        )
        if record["archived_report"] != expected_archived_evidence:
            raise BenchmarkError(f"Archived report evidence changed: {benchmark_id}")
        expanded_report = _file_record(
            record["expanded_environment_report"],
            context=f"{benchmark_id} expanded report",
        )
        if expanded_report["name"] != f"expanded-environment/{filename}":
            raise BenchmarkError(f"Expanded report name is invalid: {benchmark_id}")
        inventory = record["artifact_inventory"]
        expected_inventory = {
            "count": len(archived["artifacts"]),
            "sha256": hashlib.sha256(
                canonical_json_bytes(archived["artifacts"])
            ).hexdigest(),
        }
        if inventory != expected_inventory:
            raise BenchmarkError(f"Artifact inventory identity changed: {benchmark_id}")
        changes = record["implementation_changes"]
        if not isinstance(changes, list) or len(changes) != len(environment_paths):
            raise BenchmarkError(
                f"Environment change records are invalid: {benchmark_id}"
            )
        archived_implementation = _implementation_map(
            archived, context=f"{benchmark_id} archived report"
        )
        observed_changes: dict[str, dict[str, Any]] = {}
        for change in changes:
            if not isinstance(change, dict) or set(change) != {
                "name",
                "before",
                "after",
            }:
                raise BenchmarkError(
                    f"Environment change shape is invalid: {benchmark_id}"
                )
            name = change["name"]
            if not isinstance(name, str) or name in observed_changes:
                raise BenchmarkError(
                    f"Environment change name is invalid: {benchmark_id}"
                )
            observed_changes[name] = change
        if set(observed_changes) != set(environment_paths):
            raise BenchmarkError(
                f"Environment change inventory is invalid: {benchmark_id}"
            )
        for name in environment_paths:
            if (
                observed_changes[name]["before"] != archived_implementation[name]
                or observed_changes[name]["after"] != observed_current[name]
            ):
                raise BenchmarkError(
                    f"Environment transition evidence is invalid: {benchmark_id}/{name}"
                )
        runtime_changes = record["runtime_metadata_changes"]
        if not isinstance(runtime_changes, list) or len(runtime_changes) != 1:
            raise BenchmarkError(
                f"Runtime metadata evidence is invalid: {benchmark_id}"
            )
        runtime_change = runtime_changes[0]
        if not isinstance(runtime_change, dict) or set(runtime_change) != {
            "field",
            "reason",
            "before",
            "after",
            "conda_lock",
        }:
            raise BenchmarkError(
                f"Runtime metadata change shape is invalid: {benchmark_id}"
            )
        archived_runtime = archived["runtime"]
        archived_conda_lock = archived_implementation.get("conda-lock.yml")
        before_digest = runtime_change["before"]
        after_digest = runtime_change["after"]
        if (
            runtime_change["field"] != RSCRIPT_RUNTIME_FIELD
            or runtime_change["reason"] != RSCRIPT_RELOCATION_REASON
            or before_digest != archived_runtime.get(RSCRIPT_RUNTIME_FIELD)
            or runtime_change["conda_lock"] != archived_conda_lock
            or not isinstance(before_digest, str)
            or not isinstance(after_digest, str)
            or len(before_digest) != 64
            or len(after_digest) != 64
            or any(character not in "0123456789abcdef" for character in before_digest)
            or any(character not in "0123456789abcdef" for character in after_digest)
            or before_digest == after_digest
        ):
            raise BenchmarkError(
                f"Rscript runtime transition is invalid: {benchmark_id}"
            )
        observed_runtime_before.add(before_digest)
        observed_runtime_after.add(after_digest)
        rebuild_changes = record["rebuild_metadata_changes"]
        expected_rebuild_count = 1 if rebuild_metadata_input is not None else 0
        if not isinstance(rebuild_changes, list) or len(rebuild_changes) != (
            expected_rebuild_count
        ):
            raise BenchmarkError(
                f"Rebuild metadata evidence is invalid: {benchmark_id}"
            )
        if rebuild_metadata_input is not None:
            assert source_package is not None
            rebuild_change = rebuild_changes[0]
            if not isinstance(rebuild_change, dict) or set(rebuild_change) != {
                "input_name",
                "reason",
                "before",
                "after",
                "source_archive",
            }:
                raise BenchmarkError(
                    f"Rebuild metadata change shape is invalid: {benchmark_id}"
                )
            archived_inputs = {item["name"]: item for item in archived["inputs"]}
            if (
                rebuild_change["input_name"] != rebuild_metadata_input
                or rebuild_change["reason"] != "installed_DESCRIPTION_build_metadata"
                or rebuild_change["before"] != archived_inputs[rebuild_metadata_input]
                or rebuild_change["source_archive"] != source_archives[source_package]
            ):
                raise BenchmarkError(
                    f"Rebuild metadata transition is invalid: {benchmark_id}"
                )
            rebuilt_after = _file_record(
                rebuild_change["after"],
                context=f"{benchmark_id} rebuilt DESCRIPTION",
            )
            if rebuilt_after["name"] != rebuild_metadata_input:
                raise BenchmarkError(
                    f"Rebuilt DESCRIPTION name is invalid: {benchmark_id}"
                )

    if len(observed_runtime_before) != 1 or len(observed_runtime_after) != 1:
        raise BenchmarkError(
            "All compatibility benchmarks must share one Rscript transition"
        )

    assertions = report["assertions"]
    expected_assertion_names = {
        "all_four_expanded_environment_gates_pass",
        "numeric_artifact_sha256_and_size_identity",
        "rscript_relocation_is_conda_lock_bound",
        "installed_description_drift_is_source_archive_bound",
        "non_environment_evidence_is_byte_identical",
        "environment_transition_is_snapshot_bound",
    }
    if not isinstance(assertions, list) or len(assertions) != len(
        expected_assertion_names
    ):
        raise BenchmarkError("Compatibility assertions have an unexpected shape")
    observed_assertions: dict[str, dict[str, object]] = {}
    for assertion in assertions:
        if not isinstance(assertion, dict) or set(assertion) != {
            "name",
            "passed",
            "observed",
            "requirement",
        }:
            raise BenchmarkError("Compatibility assertion has an unexpected shape")
        name = assertion["name"]
        if (
            not isinstance(name, str)
            or name in observed_assertions
            or assertion["passed"] is not True
            or not isinstance(assertion["requirement"], str)
            or not assertion["requirement"]
        ):
            raise BenchmarkError("Compatibility assertion is invalid")
        observed_assertions[name] = assertion
    if set(observed_assertions) != expected_assertion_names:
        raise BenchmarkError("Compatibility assertion identities are invalid")
    expected_rscript_after = next(iter(observed_runtime_after))
    if (
        observed_assertions["rscript_relocation_is_conda_lock_bound"]["observed"]
        != expected_rscript_after
    ):
        raise BenchmarkError("Compatibility Rscript assertion evidence is invalid")


def _validate_environment_transition(
    archived: dict[str, Any],
    expanded: dict[str, Any],
    *,
    benchmark_id: str,
    rebuild_metadata_input: str | None,
    source_package: str | None,
    unchanged_source_archives: Mapping[str, dict[str, str]],
    live_rscript_sha256: str,
) -> tuple[
    list[dict[str, object]],
    dict[str, object],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    archived_implementation = _implementation_map(
        archived, context=f"{benchmark_id} archived report"
    )
    expanded_implementation = _implementation_map(
        expanded, context=f"{benchmark_id} expanded report"
    )
    if set(archived_implementation) != set(expanded_implementation):
        raise BenchmarkError(f"{benchmark_id} implementation inventory names changed")
    changed = {
        name
        for name in archived_implementation
        if archived_implementation[name] != expanded_implementation[name]
    }
    expected_changed = set(ENVIRONMENT_IMPLEMENTATION_PATHS)
    if changed != expected_changed:
        raise BenchmarkError(
            f"{benchmark_id} implementation changes are not environment-only: "
            f"observed={sorted(changed)}, expected={sorted(expected_changed)}"
        )

    conda_lock_evidence = file_evidence(CONDA_LOCK_PATH, name="conda-lock.yml")
    if (
        archived_implementation.get("conda-lock.yml") != conda_lock_evidence
        or expanded_implementation.get("conda-lock.yml") != conda_lock_evidence
    ):
        raise BenchmarkError(
            f"{benchmark_id} Rscript transition is not bound to current conda-lock.yml"
        )
    archived_runtime = archived.get("runtime")
    expanded_runtime = expanded.get("runtime")
    if not isinstance(archived_runtime, dict) or not isinstance(expanded_runtime, dict):
        raise BenchmarkError(f"{benchmark_id} runtime metadata must be objects")
    if set(archived_runtime) != set(expanded_runtime):
        raise BenchmarkError(f"{benchmark_id} runtime metadata fields changed")
    changed_runtime = {
        field
        for field in archived_runtime
        if archived_runtime[field] != expanded_runtime[field]
    }
    if changed_runtime != {RSCRIPT_RUNTIME_FIELD}:
        raise BenchmarkError(
            f"{benchmark_id} runtime changes are not the approved Rscript relocation: "
            f"observed={sorted(changed_runtime)}"
        )
    before_rscript = archived_runtime[RSCRIPT_RUNTIME_FIELD]
    after_rscript = expanded_runtime[RSCRIPT_RUNTIME_FIELD]
    if (
        not isinstance(before_rscript, str)
        or len(before_rscript) != 64
        or any(character not in "0123456789abcdef" for character in before_rscript)
        or not isinstance(after_rscript, str)
        or len(after_rscript) != 64
        or any(character not in "0123456789abcdef" for character in after_rscript)
        or after_rscript != live_rscript_sha256
    ):
        raise BenchmarkError(
            f"{benchmark_id} Rscript digest does not match the live locked runtime"
        )
    runtime_metadata_changes = [
        {
            "field": RSCRIPT_RUNTIME_FIELD,
            "reason": RSCRIPT_RELOCATION_REASON,
            "before": before_rscript,
            "after": after_rscript,
            "conda_lock": conda_lock_evidence,
        }
    ]

    if expanded["artifacts"] != archived["artifacts"]:
        raise BenchmarkError(
            f"{benchmark_id} numeric artifact SHA-256/size inventory changed"
        )
    artifact_inventory = archived["artifacts"]
    artifact_identity = {
        "count": len(artifact_inventory),
        "sha256": hashlib.sha256(canonical_json_bytes(artifact_inventory)).hexdigest(),
    }

    archived_inputs = {
        item["name"]: item
        for item in archived["inputs"]
        if isinstance(item, dict) and isinstance(item.get("name"), str)
    }
    expanded_inputs = {
        item["name"]: item
        for item in expanded["inputs"]
        if isinstance(item, dict) and isinstance(item.get("name"), str)
    }
    if (
        len(archived_inputs) != len(archived["inputs"])
        or len(expanded_inputs) != len(expanded["inputs"])
        or set(archived_inputs) != set(expanded_inputs)
    ):
        raise BenchmarkError(f"{benchmark_id} input inventory shape changed")
    changed_inputs = {
        name
        for name in archived_inputs
        if archived_inputs[name] != expanded_inputs[name]
    }
    expected_input_changes = (
        {rebuild_metadata_input} if rebuild_metadata_input is not None else set()
    )
    if changed_inputs != expected_input_changes:
        raise BenchmarkError(
            f"{benchmark_id} input changes are not the approved rebuild metadata: "
            f"observed={sorted(changed_inputs)}, expected={sorted(expected_input_changes)}"
        )
    rebuild_metadata_changes: list[dict[str, object]] = []
    if rebuild_metadata_input is not None:
        if source_package is None or source_package not in unchanged_source_archives:
            raise BenchmarkError(
                f"{benchmark_id} lacks unchanged source-archive evidence"
            )
        source_archive = unchanged_source_archives[source_package]
        if (
            archived["runtime"].get(source_package) != source_archive["version"]
            or expanded["runtime"].get(source_package) != source_archive["version"]
        ):
            raise BenchmarkError(
                f"{benchmark_id} runtime package version disagrees with its source archive"
            )
        before_input = _file_record(
            archived_inputs[rebuild_metadata_input],
            context=f"{benchmark_id} archived rebuild metadata",
        )
        after_input = _file_record(
            expanded_inputs[rebuild_metadata_input],
            context=f"{benchmark_id} expanded rebuild metadata",
        )
        rebuild_metadata_changes.append(
            {
                "input_name": rebuild_metadata_input,
                "reason": "installed_DESCRIPTION_build_metadata",
                "before": before_input,
                "after": after_input,
                "source_archive": source_archive,
            }
        )

    changes: list[dict[str, object]] = []
    for name, current_path in ENVIRONMENT_IMPLEMENTATION_PATHS.items():
        before = archived_implementation[name]
        resolve_archived_implementation_path(
            current_path,
            expected_sha256=before["sha256"],  # type: ignore[arg-type]
            expected_size=before["size_bytes"],  # type: ignore[arg-type]
        )
        after = expanded_implementation[name]
        expected_after = file_evidence(current_path, name=name)
        if after != expected_after:
            raise BenchmarkError(
                f"{benchmark_id} does not bind current environment file {name!r}"
            )
        changes.append({"name": name, "before": before, "after": after})

    normalized_expanded = copy.deepcopy(expanded)
    normalized_expanded["implementation"] = [
        archived_implementation[item["name"]]
        if item["name"] in expected_changed
        else item
        for item in normalized_expanded["implementation"]
    ]
    if rebuild_metadata_input is not None:
        normalized_expanded["inputs"] = [
            archived_inputs[item["name"]]
            if item["name"] == rebuild_metadata_input
            else item
            for item in normalized_expanded["inputs"]
        ]
    normalized_expanded["runtime"][RSCRIPT_RUNTIME_FIELD] = before_rscript
    if normalized_expanded != archived:
        raise BenchmarkError(
            f"{benchmark_id} changed outside approved environment/rebuild metadata records"
        )
    return (
        changes,
        artifact_identity,
        runtime_metadata_changes,
        rebuild_metadata_changes,
    )


def build_compatibility_report(
    report_directory: Path,
    *,
    rscript_path: Path | None = None,
) -> dict[str, object]:
    """Validate four real reruns and return their canonical compatibility report."""

    if report_directory.is_symlink() or not report_directory.is_dir():
        raise BenchmarkError(
            f"Expanded-environment report directory is unavailable: {report_directory}"
        )
    snapshot_records = load_environment_snapshots(SNAPSHOT_ROOT)
    baseline_records = [
        record
        for record in snapshot_records
        if record.snapshot_id == BASELINE_SNAPSHOT_ID
    ]
    if not baseline_records or {
        record.source_revision for record in baseline_records
    } != {SOURCE_REVISION}:
        raise BenchmarkError("The required P0/C1/C2 baseline snapshot is unavailable")
    source_archives = _unchanged_source_archives()
    live_rscript = rscript_from_runtime(
        str(rscript_path) if rscript_path is not None else None
    )
    live_rscript_sha256 = sha256_file(live_rscript)

    benchmark_records: list[dict[str, object]] = []
    for specification in BENCHMARKS:
        benchmark_id = specification["benchmark_id"]
        filename = specification["filename"]
        assert isinstance(benchmark_id, str)
        assert isinstance(filename, str)
        archived = _assert_frozen_archive(specification)
        expanded_path = report_directory / filename
        if expanded_path.is_symlink() or not expanded_path.is_file():
            raise BenchmarkError(
                f"Expanded-environment report is missing: {expanded_path}"
            )
        expanded = read_json(expanded_path)
        assert_report_shape(expanded, benchmark_id=benchmark_id)
        if expanded["status"] != "pass":
            raise BenchmarkError(
                f"Expanded-environment benchmark is not passing: {expanded_path}"
            )
        rebuild_metadata_input = specification["rebuild_metadata_input"]
        source_package = specification["source_package"]
        assert rebuild_metadata_input is None or isinstance(rebuild_metadata_input, str)
        assert source_package is None or isinstance(source_package, str)
        (
            changes,
            artifact_identity,
            runtime_metadata_changes,
            rebuild_metadata_changes,
        ) = _validate_environment_transition(
            archived,
            expanded,
            benchmark_id=benchmark_id,
            rebuild_metadata_input=rebuild_metadata_input,
            source_package=source_package,
            unchanged_source_archives=source_archives,
            live_rscript_sha256=live_rscript_sha256,
        )
        archived_path = specification["archived_path"]
        assert isinstance(archived_path, Path)
        benchmark_records.append(
            {
                "benchmark_id": benchmark_id,
                "archived_report": file_evidence(
                    archived_path,
                    name=archived_path.relative_to(PROJECT_ROOT).as_posix(),
                ),
                "expanded_environment_report": file_evidence(
                    expanded_path,
                    name=f"expanded-environment/{filename}",
                ),
                "artifact_inventory": artifact_identity,
                "implementation_changes": changes,
                "runtime_metadata_changes": runtime_metadata_changes,
                "rebuild_metadata_changes": rebuild_metadata_changes,
            }
        )

    observed_rscript_transitions = {
        (
            record["runtime_metadata_changes"][0]["before"],
            record["runtime_metadata_changes"][0]["after"],
        )
        for record in benchmark_records
    }
    if len(observed_rscript_transitions) != 1:
        raise BenchmarkError(
            "All four expanded-environment reports must share one Rscript transition"
        )

    current_environment = [
        file_evidence(path, name=name)
        for name, path in ENVIRONMENT_IMPLEMENTATION_PATHS.items()
    ]
    assertions = [
        {
            "name": "all_four_expanded_environment_gates_pass",
            "passed": True,
            "observed": [record["benchmark_id"] for record in benchmark_records],
            "requirement": "all four frozen P0/C1/C2 benchmarks rerun and pass",
        },
        {
            "name": "numeric_artifact_sha256_and_size_identity",
            "passed": True,
            "observed": sum(
                record["artifact_inventory"]["count"] for record in benchmark_records
            ),
            "requirement": (
                "every archived numeric artifact name, SHA-256 digest, and byte size "
                "is identical in the expanded-environment rerun"
            ),
        },
        {
            "name": "rscript_relocation_is_conda_lock_bound",
            "passed": True,
            "observed": live_rscript_sha256,
            "requirement": (
                "all four reports use the same live Rscript digest, and only that "
                "runtime field changes under the unchanged conda-lock.yml"
            ),
        },
        {
            "name": "installed_description_drift_is_source_archive_bound",
            "passed": True,
            "observed": sum(
                len(record["rebuild_metadata_changes"]) for record in benchmark_records
            ),
            "requirement": (
                "the only rebuilt installed-package DESCRIPTION records are disclosed "
                "and bind to unchanged versioned source archives"
            ),
        },
        {
            "name": "non_environment_evidence_is_byte_identical",
            "passed": True,
            "observed": len(benchmark_records),
            "requirement": (
                "each rerun report equals its archive after restoring only the approved "
                "environment and installed DESCRIPTION build-metadata records"
            ),
        },
        {
            "name": "environment_transition_is_snapshot_bound",
            "passed": True,
            "observed": sorted(ENVIRONMENT_IMPLEMENTATION_PATHS),
            "requirement": (
                "old hashes resolve through the immutable baseline snapshot and new "
                "hashes match the current locked environment"
            ),
        },
    ]
    report = {
        "schema_version": SCHEMA_VERSION,
        "compatibility_id": COMPATIBILITY_ID,
        "status": "pass",
        "baseline_snapshot": {
            "snapshot_id": BASELINE_SNAPSHOT_ID,
            "source_revision": SOURCE_REVISION,
            "manifest": manifest_evidence(BASELINE_SNAPSHOT_ID),
        },
        "expanded_environment": {"files": current_environment},
        "permitted_implementation_changes": sorted(ENVIRONMENT_IMPLEMENTATION_PATHS),
        "benchmarks": benchmark_records,
        "assertions": assertions,
    }
    assert_compatibility_report_shape(report)
    return report


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reports",
        type=Path,
        required=True,
        help="Directory containing the four real expanded-environment benchmark reports.",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--rscript",
        type=Path,
        help=(
            "Rscript executable used by all four reruns; defaults to the executable "
            "beside the active locked Python interpreter."
        ),
    )
    return parser.parse_args()


def main() -> int:
    arguments = _arguments()
    if arguments.output.exists() or arguments.output.is_symlink():
        raise BenchmarkError(
            f"Refusing to overwrite compatibility evidence: {arguments.output}"
        )
    report = build_compatibility_report(
        arguments.reports,
        rscript_path=arguments.rscript,
    )
    write_json(arguments.output, report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
