"""Evidence-coupled adapter for the locked edgeR v4 QL backend."""

from __future__ import annotations

import copy
import csv
import ctypes
import errno
import hashlib
from importlib import resources
import io
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tempfile
from typing import Any, Mapping, Sequence

from .analysis_contract import (
    ANALYSIS_SCHEMA_VERSION,
    ValidatedAnalysisRequest,
    _parse_contrasts,
    _parse_design,
    load_analysis_request,
)
from .errors import (
    BackendFailedError,
    ContrastNotEstimableError,
    DesignRankDeficientError,
    InputIntegrityError,
    InputReadError,
    InvalidRequestError,
    OutputWriteError,
    ToolkitError,
)
from .provenance import has_control_characters, read_json_object

BACKEND_NAME = "edgeR_QL"
BACKEND_OUTPUT_SCHEMA_VERSION = "1.0"
PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION = "1.1"
_EXPECTED_RUNTIME = {
    "R": "4.6.1",
    "Bioconductor": "3.23",
    "BiocVersion_package": "3.23.1",
    "edgeR": "4.10.0",
    "tximport": "1.40.0",
    "limma": "3.68.0",
}
_OUTPUT_MEMBERS = {
    "analysis.json",
    "coefficients.tsv",
    "design.tsv",
    "results.tsv",
}
_PATHWAY_OUTPUT_MEMBERS = _OUTPUT_MEMBERS | {"pathway_results.tsv"}
_AT_FDCWD = -100
_RENAME_NOREPLACE = 1


def _canonical_bytes(document: Mapping[str, Any]) -> bytes:
    try:
        return json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    except (TypeError, ValueError, OverflowError) as error:
        raise InvalidRequestError(
            "The normalized edgeR request is not strict JSON.",
            details={"cause_type": type(error).__name__},
            cause=error,
        ) from error


def _write_private_json(path: Path, document: Mapping[str, Any]) -> None:
    payload = _canonical_bytes(document) + b"\n"
    try:
        with path.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        path.chmod(0o400)
    except OSError as error:
        raise OutputWriteError(
            "The private edgeR request could not be written.",
            path=path,
            operation="write_private_backend_request",
            cause=error,
        ) from error


def _resolve_output(path: str | Path) -> Path:
    if not isinstance(path, (str, Path)) or not str(path):
        raise InvalidRequestError("The backend output directory must be non-empty.")
    if has_control_characters(str(path)):
        raise InvalidRequestError(
            "The backend output directory contains a control character."
        )
    candidate = Path(path)
    if candidate.is_symlink():
        raise InvalidRequestError(
            "The backend output path must not be a symbolic link.",
            details={"output_dir": str(candidate)},
        )
    try:
        resolved = candidate.resolve(strict=False)
    except (OSError, RuntimeError) as error:
        raise OutputWriteError(
            "The backend output path could not be resolved.",
            path=candidate,
            operation="resolve_backend_output",
            cause=error,
        ) from error
    if resolved.exists() or resolved.is_symlink():
        raise InvalidRequestError(
            "The backend output path already exists; results are never overwritten.",
            details={"output_dir": str(resolved)},
        )
    return resolved


def _artifact_for_role(
    artifacts: Sequence[Mapping[str, Any]],
    role: str,
    *,
    expected_path: str | Path,
) -> dict[str, Any]:
    matches = [dict(item) for item in artifacts if item.get("role") == role]
    if len(matches) != 1:
        raise InputIntegrityError(
            "Validated input evidence does not identify one required artifact.",
            details={"role": role, "match_count": len(matches)},
        )
    record = matches[0]
    try:
        observed = Path(record["path"]).resolve(strict=True)
        expected = Path(expected_path).resolve(strict=True)
    except (KeyError, OSError, RuntimeError, TypeError) as error:
        raise InputIntegrityError(
            "A validated artifact path cannot be resolved for execution.",
            details={"role": role},
            cause=error,
        ) from error
    if observed != expected:
        raise InputIntegrityError(
            "Normalized input data and validated artifact evidence disagree.",
            details={
                "role": role,
                "normalized_path": str(expected),
                "evidence_path": str(observed),
            },
        )
    return record


def _copy_evidence_snapshot(
    artifact: Mapping[str, Any],
    destination: Path,
    *,
    capture_bytes: bool = False,
) -> tuple[dict[str, Any], bytes | None]:
    """Copy and hash one source stream, coupling R bytes to validated evidence."""

    source = Path(str(artifact["path"]))
    expected_digest = artifact.get("sha256")
    expected_size = artifact.get("size_bytes")
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    digest = hashlib.sha256()
    captured = bytearray() if capture_bytes else None
    try:
        resolved = source.resolve(strict=True)
        if not resolved.is_file():
            raise IsADirectoryError(str(resolved))
        before = resolved.stat()
        with resolved.open("rb") as source_handle, destination.open("xb") as target:
            observed_size = 0
            for chunk in iter(lambda: source_handle.read(1024 * 1024), b""):
                digest.update(chunk)
                target.write(chunk)
                observed_size += len(chunk)
                if captured is not None:
                    captured.extend(chunk)
            target.flush()
            os.fsync(target.fileno())
        after = resolved.stat()
        destination.chmod(0o400)
    except OSError as error:
        raise InputReadError(
            "A validated input could not be materialized for the R backend.",
            path=source,
            operation="materialize_verified_backend_input",
            cause=error,
            details={"role": artifact.get("role")},
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
    observed_digest = digest.hexdigest()
    if (
        identity_before != identity_after
        or observed_digest != expected_digest
        or observed_size != expected_size
    ):
        raise InputIntegrityError(
            "An input no longer matches the bytes accepted by validation.",
            details={
                "role": artifact.get("role"),
                "path": str(resolved),
                "expected_sha256": expected_digest,
                "observed_sha256": observed_digest,
                "expected_size_bytes": expected_size,
                "observed_size_bytes": observed_size,
                "source_identity_stable": identity_before == identity_after,
            },
        )
    snapshot = {
        "role": artifact["role"],
        "source_path": str(resolved),
        "sha256": observed_digest,
        "size_bytes": observed_size,
        "private_relative_path": destination.name,
        "coupling": "copied_and_hashed_from_one_source_stream",
    }
    return snapshot, bytes(captured) if captured is not None else None


def _parse_metadata_snapshot(
    content: bytes,
    *,
    sample_id_column: str,
    sample_order: Sequence[str],
    design_terms: Sequence[str],
) -> dict[str, list[str]]:
    try:
        text = content.decode("utf-8")
        rows = list(
            csv.reader(io.StringIO(text, newline=""), delimiter="\t", strict=True)
        )
    except (UnicodeError, csv.Error) as error:
        raise InputIntegrityError(
            "The private metadata snapshot cannot be parsed identically to validation.",
            details={"cause_type": type(error).__name__},
            cause=error,
        ) from error
    if not rows or not rows[0]:
        raise InputIntegrityError("The private metadata snapshot is empty.")
    header = rows[0]
    required = [sample_id_column, *design_terms]
    missing = sorted(set(required) - set(header))
    if missing:
        raise InvalidRequestError(
            "The statistical design names columns absent from validated metadata.",
            details={"missing_metadata_columns": missing},
        )
    if len(set(header)) != len(header):
        raise InputIntegrityError(
            "The private metadata snapshot has duplicate columns."
        )
    by_sample: dict[str, list[str]] = {}
    sample_index = header.index(sample_id_column)
    for row_number, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise InputIntegrityError(
                "The private metadata snapshot has an invalid row width.",
                details={"row": row_number},
            )
        sample_id = row[sample_index]
        if sample_id in by_sample:
            raise InputIntegrityError(
                "The private metadata snapshot contains duplicate samples.",
                details={"sample_id": sample_id},
            )
        by_sample[sample_id] = row
    if set(by_sample) != set(sample_order):
        raise InputIntegrityError(
            "The private metadata snapshot no longer matches the assay sample set."
        )
    selected = [sample_id_column, *design_terms]
    return {
        column: [by_sample[sample][header.index(column)] for sample in sample_order]
        for column in selected
    }


def _assert_p0_route(input_data: Mapping[str, Any]) -> None:
    semantics = input_data.get("input_semantics")
    route = input_data.get("route")
    if not isinstance(route, Mapping):
        raise InputIntegrityError("The validated input route is missing.")
    if semantics == "featurecounts_integer":
        expected = {
            "backend_input": "edgeR::DGEList",
            "count_semantics": "integer",
            "transcript_length_offset": False,
        }
    elif semantics == "salmon_quant_dirs_full_length":
        expected = {
            "edgeR_constructor": "edgeR::DGEListFromTximport",
            "count_source": "txi$counts",
            "transcript_length_offset": True,
        }
        tximport = route.get("tximport")
        if not isinstance(tximport, Mapping) or (
            tximport.get("countsFromAbundance") != "no"
            or tximport.get("dropInfReps") is not False
        ):
            raise InputIntegrityError(
                "The validated full-length tximport route is incompatible."
            )
        options = route.get("edgeR_options")
        if not isinstance(options, Mapping) or not isinstance(
            options.get("divide"), bool
        ):
            raise InputIntegrityError(
                "The validated full-length edgeR divide policy is missing."
            )
    elif semantics == "salmon_quant_dirs_three_prime":
        expected = {
            "edgeR_constructor": "edgeR::DGEList",
            "count_source": "txi$counts",
            "transcript_length_offset": False,
            "gene_length_correction": False,
            "certified_path_execution_permitted": True,
        }
        tximport = route.get("tximport")
        if not isinstance(tximport, Mapping) or (
            tximport.get("countsFromAbundance") != "no"
            or tximport.get("dropInfReps") is not False
        ):
            raise InputIntegrityError(
                "The validated three-prime tximport route is incompatible."
            )
    else:
        raise InvalidRequestError(
            "Only the three validated P0 input semantics can reach edgeR.",
            details={"input_semantics": semantics},
        )
    mismatches = {
        key: {"expected": value, "observed": route.get(key)}
        for key, value in expected.items()
        if route.get(key) != value
    }
    if mismatches:
        raise InputIntegrityError(
            "The validated input route does not match the fixed edgeR P0 policy.",
            details={"mismatches": mismatches},
        )


def _materialize_validated_inputs(
    validated: ValidatedAnalysisRequest,
    input_root: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, list[str]]]:
    rewritten = copy.deepcopy(validated.input_data)
    _assert_p0_route(rewritten)
    artifacts = validated.consumed_artifacts
    snapshots: list[dict[str, Any]] = []

    metadata_path = rewritten["metadata"]["path"]
    metadata_artifact = _artifact_for_role(
        artifacts, "metadata", expected_path=metadata_path
    )
    metadata_destination = input_root / "metadata.tsv"
    snapshot, metadata_content = _copy_evidence_snapshot(
        metadata_artifact,
        metadata_destination,
        capture_bytes=True,
    )
    snapshot["private_relative_path"] = "metadata.tsv"
    snapshots.append(snapshot)
    rewritten["metadata"]["path"] = str(metadata_destination)

    semantics = rewritten["input_semantics"]
    if semantics == "featurecounts_integer":
        featurecounts = rewritten["featurecounts"]
        if featurecounts["layout"] == "combined_matrix":
            artifact = _artifact_for_role(
                artifacts,
                "featurecounts.matrix",
                expected_path=featurecounts["matrix_path"],
            )
            destination = input_root / "featurecounts" / "counts.tsv"
            record, _ = _copy_evidence_snapshot(artifact, destination)
            record["private_relative_path"] = "featurecounts/counts.tsv"
            snapshots.append(record)
            featurecounts["matrix_path"] = str(destination)
        elif featurecounts["layout"] == "per_sample_files":
            for index, file_record in enumerate(featurecounts["files"]):
                role = f"featurecounts.sample.{file_record['sample_id']}"
                artifact = _artifact_for_role(
                    artifacts, role, expected_path=file_record["path"]
                )
                destination = input_root / "featurecounts" / f"sample-{index + 1}.txt"
                record, _ = _copy_evidence_snapshot(artifact, destination)
                record["private_relative_path"] = (
                    f"featurecounts/sample-{index + 1}.txt"
                )
                snapshots.append(record)
                file_record["path"] = str(destination)
        else:
            raise InputIntegrityError(
                "The normalized featureCounts layout is unsupported."
            )
    else:
        salmon = rewritten["salmon"]
        tx2gene_artifact = _artifact_for_role(
            artifacts, "salmon.tx2gene", expected_path=salmon["tx2gene_path"]
        )
        tx2gene_destination = input_root / "salmon" / "tx2gene.tsv"
        record, _ = _copy_evidence_snapshot(tx2gene_artifact, tx2gene_destination)
        record["private_relative_path"] = "salmon/tx2gene.tsv"
        snapshots.append(record)
        salmon["tx2gene_path"] = str(tx2gene_destination)
        for index, sample in enumerate(salmon["samples"]):
            sample_id = sample["sample_id"]
            quant_root = input_root / "salmon" / f"sample-{index + 1}"
            sources = {
                f"salmon.quant.{sample_id}": (
                    Path(sample["quant_dir"]) / "quant.sf",
                    quant_root / "quant.sf",
                    f"salmon/sample-{index + 1}/quant.sf",
                ),
                f"salmon.cmd_info.{sample_id}": (
                    Path(sample["quant_dir"]) / "cmd_info.json",
                    quant_root / "cmd_info.json",
                    f"salmon/sample-{index + 1}/cmd_info.json",
                ),
                f"salmon.meta_info.{sample_id}": (
                    Path(sample["aux_dir"]) / "meta_info.json",
                    quant_root
                    / PurePosixPath(sample["aux_dir_declared"])
                    / "meta_info.json",
                    f"salmon/sample-{index + 1}/{sample['aux_dir_declared']}/meta_info.json",
                ),
            }
            replicate_prefix = f"salmon.inferential_replicates.{sample_id}"
            available_roles = {item["role"] for item in artifacts}
            if f"{replicate_prefix}.names.tsv.gz" in available_roles:
                bootstrap_source = Path(sample["aux_dir"]) / "bootstrap"
                sources[f"{replicate_prefix}.names.tsv.gz"] = (
                    bootstrap_source / "names.tsv.gz",
                    quant_root
                    / PurePosixPath(sample["aux_dir_declared"])
                    / "bootstrap"
                    / "names.tsv.gz",
                    f"salmon/sample-{index + 1}/{sample['aux_dir_declared']}/bootstrap/names.tsv.gz",
                )
                sources[f"{replicate_prefix}.bootstraps.gz"] = (
                    bootstrap_source / "bootstraps.gz",
                    quant_root
                    / PurePosixPath(sample["aux_dir_declared"])
                    / "bootstrap"
                    / "bootstraps.gz",
                    f"salmon/sample-{index + 1}/{sample['aux_dir_declared']}/bootstrap/bootstraps.gz",
                )
            for role, (source, destination, relative) in sources.items():
                artifact = _artifact_for_role(artifacts, role, expected_path=source)
                copied, _ = _copy_evidence_snapshot(artifact, destination)
                copied["private_relative_path"] = relative
                snapshots.append(copied)
            sample["quant_dir"] = str(quant_root)
            sample["aux_dir"] = str(
                quant_root / PurePosixPath(sample["aux_dir_declared"])
            )

    assert metadata_content is not None
    metadata_values = _parse_metadata_snapshot(
        metadata_content,
        sample_id_column=rewritten["metadata"]["sample_id_column"],
        sample_order=rewritten["sample_order"],
        design_terms=validated.design["terms"],
    )
    rewritten["metadata_values"] = metadata_values
    return rewritten, snapshots, metadata_values


def _materialize_gene_set_inputs(
    specification: Mapping[str, Any],
    input_root: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Recapture frozen pathway inputs into the private R workspace.

    The public request parser has already captured and verified each declared
    digest.  This second pass deliberately hashes the same source stream that
    is copied for R, closing the request-load-to-execution race.
    """

    rewritten = copy.deepcopy(dict(specification))
    snapshots: list[dict[str, Any]] = []
    destinations = {
        "gmt": ("pathways.gmt", input_root / "gene-sets" / "sets.gmt"),
        "annotation": (
            "pathways.annotation",
            input_root / "gene-sets" / "annotation.tsv",
        ),
    }
    for source_kind, (role, destination) in destinations.items():
        source = rewritten.get(source_kind)
        if not isinstance(source, Mapping):
            raise InputIntegrityError(
                "The normalized pathway source evidence is incomplete.",
                details={"source_kind": source_kind},
            )
        artifact = {
            "role": role,
            "path": source.get("path"),
            "sha256": source.get("sha256"),
            "size_bytes": source.get("size_bytes"),
        }
        record, _ = _copy_evidence_snapshot(artifact, destination)
        record["private_relative_path"] = str(
            PurePosixPath("gene-sets") / destination.name
        )
        snapshots.append(record)
        rewritten[source_kind]["path"] = str(destination)
    return rewritten, snapshots


def _capture_declared_benchmark_file(
    source: str | Path,
    destination: Path,
    *,
    role: str,
) -> tuple[dict[str, Any], bytes]:
    candidate = Path(source)
    try:
        resolved = candidate.resolve(strict=True)
        if not resolved.is_file():
            raise IsADirectoryError(str(resolved))
        before = resolved.stat()
        digest = hashlib.sha256()
        captured = bytearray()
        destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        with resolved.open("rb") as source_handle, destination.open("xb") as target:
            for chunk in iter(lambda: source_handle.read(1024 * 1024), b""):
                digest.update(chunk)
                captured.extend(chunk)
                target.write(chunk)
            target.flush()
            os.fsync(target.fileno())
        after = resolved.stat()
        destination.chmod(0o400)
    except OSError as error:
        raise InputReadError(
            "A benchmark-kernel input could not be captured.",
            path=candidate,
            operation="capture_benchmark_kernel_input",
            cause=error,
            details={"role": role},
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
    if identity_before != identity_after or len(captured) != after.st_size:
        raise InputIntegrityError(
            "A benchmark-kernel input changed while it was captured.",
            details={"role": role, "path": str(resolved)},
        )
    record = {
        "role": role,
        "source_path": str(resolved),
        "sha256": digest.hexdigest(),
        "size_bytes": len(captured),
        "private_relative_path": destination.name,
        "coupling": "copied_and_hashed_from_one_source_stream",
    }
    return record, bytes(captured)


def _parse_benchmark_counts_header(content: bytes) -> list[str]:
    try:
        rows = list(
            csv.reader(
                io.StringIO(content.decode("utf-8"), newline=""),
                delimiter="\t",
                strict=True,
            )
        )
    except (UnicodeError, csv.Error) as error:
        raise InvalidRequestError(
            "The benchmark count matrix is not strict UTF-8 TSV.",
            details={"cause_type": type(error).__name__},
            cause=error,
        ) from error
    if not rows or len(rows[0]) < 3 or rows[0][0] != "gene_id":
        raise InvalidRequestError(
            "The benchmark count header must be gene_id followed by at least two samples."
        )
    header = rows[0]
    if len(set(header)) != len(header) or any(not value for value in header):
        raise InvalidRequestError(
            "The benchmark count header must contain unique non-empty columns."
        )
    for row_number, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise InvalidRequestError(
                "A benchmark count row has the wrong width.",
                details={"row": row_number},
            )
    if len(rows) == 1:
        raise InvalidRequestError("The benchmark count matrix contains no genes.")
    return header[1:]


def _r_script_path() -> Path:
    resource = resources.files("rnaseq_downstream").joinpath("r_scripts", "edger_ql.R")
    path = Path(str(resource))
    if not path.is_file():
        raise BackendFailedError(
            "The packaged edgeR backend resource is missing.",
            details={"resource": "r_scripts/edger_ql.R"},
        )
    return path


def _parse_backend_stdout(stdout: str, *, returncode: int) -> dict[str, Any]:
    lines = stdout.splitlines()
    if len(lines) != 1 or not lines[0]:
        raise BackendFailedError(
            "The R backend violated its one-document stdout contract.",
            details={"returncode": returncode, "stdout_line_count": len(lines)},
        )
    try:
        response = json.loads(
            lines[0],
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-standard constant {value}")
            ),
        )
    except (json.JSONDecodeError, ValueError) as error:
        raise BackendFailedError(
            "The R backend returned invalid JSON.",
            details={"returncode": returncode, "cause_type": type(error).__name__},
            cause=error,
        ) from error
    if not isinstance(response, dict):
        raise BackendFailedError(
            "The R backend response root is not an object.",
            details={"returncode": returncode},
        )
    return response


def _raise_backend_response(
    response: Mapping[str, Any],
    *,
    returncode: int,
    stderr: str,
) -> None:
    errors = response.get("errors")
    if (
        not isinstance(errors, list)
        or len(errors) != 1
        or not isinstance(errors[0], Mapping)
    ):
        raise BackendFailedError(
            "The failed R backend response lacks one structured error.",
            details={"returncode": returncode, "stderr": stderr[-8000:]},
        )
    error = errors[0]
    code = error.get("code")
    message = str(error.get("message", "The edgeR backend failed."))
    details = error.get("details")
    if not isinstance(details, Mapping):
        details = {}
    normalized_details = {
        **dict(details),
        "backend_returncode": returncode,
        "backend_stderr": stderr[-8000:],
    }
    if code == "DESIGN_RANK_DEFICIENT":
        raise DesignRankDeficientError(message, details=normalized_details)
    if code == "CONTRAST_NOT_ESTIMABLE":
        raise ContrastNotEstimableError(message, details=normalized_details)
    raise BackendFailedError(message, details=normalized_details)


def _invoke_r(
    backend_request: Path,
    result_stage: Path,
    *,
    rscript: str | Path,
    r_library: str | Path | None,
    expected_schema_version: str = BACKEND_OUTPUT_SCHEMA_VERSION,
) -> tuple[dict[str, Any], str]:
    environment = os.environ.copy()
    if r_library is not None:
        try:
            library_path = Path(r_library).resolve(strict=True)
        except (OSError, RuntimeError, TypeError) as error:
            raise InvalidRequestError(
                "The requested R package library cannot be resolved.",
                details={
                    "r_library": str(r_library),
                    "cause_type": type(error).__name__,
                },
                cause=error,
            ) from error
        if not library_path.is_dir():
            raise InvalidRequestError(
                "The requested R package library is not a directory.",
                details={"r_library": str(library_path)},
            )
        environment["R_LIBS_USER"] = str(library_path)
    command = [
        str(rscript),
        "--vanilla",
        str(_r_script_path()),
        str(backend_request),
        str(result_stage),
    ]
    try:
        completed = subprocess.run(
            command,
            check=False,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="strict",
            env=environment,
        )
    except (OSError, UnicodeError) as error:
        raise BackendFailedError(
            "The independent R backend process could not be executed.",
            details={
                "cause_type": type(error).__name__,
                "rscript": str(rscript),
            },
            cause=error,
        ) from error
    response = _parse_backend_stdout(completed.stdout, returncode=completed.returncode)
    if completed.returncode != 0 or response.get("status") != "success":
        _raise_backend_response(
            response,
            returncode=completed.returncode,
            stderr=completed.stderr,
        )
    if (
        response.get("backend") != BACKEND_NAME
        or response.get("schema_version") != expected_schema_version
    ):
        raise BackendFailedError(
            "The R backend response has an unexpected identity.",
            details={
                "observed_backend": response.get("backend"),
                "observed_schema_version": response.get("schema_version"),
                "expected_schema_version": expected_schema_version,
            },
        )
    return response, completed.stderr


def _capture_output_json(path: Path, *, role: str) -> tuple[dict[str, Any], str, int]:
    try:
        before = path.stat()
        with path.open("rb") as handle:
            content = handle.read()
        after = path.stat()
    except OSError as error:
        raise BackendFailedError(
            "A backend JSON artifact could not be captured.",
            details={"path": str(path), "role": role},
            cause=error,
        ) from error
    identity_before = (before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after or len(content) != after.st_size:
        raise BackendFailedError(
            "A backend JSON artifact changed during verification.",
            details={"path": str(path), "role": role},
        )
    try:
        document = read_json_object(path, document_role=role, content=content)
    except (InvalidRequestError, InputReadError) as error:
        raise BackendFailedError(
            "A backend JSON artifact does not contain a valid strict JSON object.",
            details={
                "path": str(path),
                "role": role,
                "parse_error_code": error.code.value,
                "parse_error_details": dict(error.details),
            },
            cause=error,
        ) from error
    return document, hashlib.sha256(content).hexdigest(), len(content)


def _hash_and_fsync_output(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
            os.fsync(handle.fileno())
        size = path.stat().st_size
    except OSError as error:
        raise BackendFailedError(
            "A backend output artifact could not be verified.",
            details={"path": str(path)},
            cause=error,
        ) from error
    return digest.hexdigest(), size


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _verify_result_stage(
    result_stage: Path,
    *,
    execution_scope: str,
    expected_schema_version: str = BACKEND_OUTPUT_SCHEMA_VERSION,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if expected_schema_version == BACKEND_OUTPUT_SCHEMA_VERSION:
        expected_members = _OUTPUT_MEMBERS
    elif expected_schema_version == PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION:
        expected_members = _PATHWAY_OUTPUT_MEMBERS
    else:
        raise BackendFailedError(
            "The expected backend output schema is unsupported.",
            details={"expected_schema_version": expected_schema_version},
        )
    expected_files = expected_members | {"backend_manifest.json"}
    try:
        entries = list(result_stage.iterdir())
    except OSError as error:
        raise BackendFailedError(
            "The R backend did not create a readable result stage.",
            details={"result_stage": str(result_stage)},
            cause=error,
        ) from error
    names = {entry.name for entry in entries}
    unsafe = sorted(
        entry.name for entry in entries if entry.is_symlink() or not entry.is_file()
    )
    if names != expected_files or unsafe:
        raise BackendFailedError(
            "The R backend output inventory is incomplete or unsafe.",
            details={
                "missing_files": sorted(expected_files - names),
                "unexpected_files": sorted(names - expected_files),
                "unsafe_files": unsafe,
            },
        )
    manifest_path = result_stage / "backend_manifest.json"
    manifest, manifest_digest, manifest_size = _capture_output_json(
        manifest_path, role="backend_manifest"
    )
    if (
        manifest.get("schema_version") != expected_schema_version
        or manifest.get("kind") != "edger_ql_backend_manifest"
        or manifest.get("backend") != BACKEND_NAME
        or manifest.get("execution_scope") != execution_scope
        or manifest.get("runtime_identity") != _EXPECTED_RUNTIME
    ):
        raise BackendFailedError(
            "The backend manifest identity is incompatible.",
            details={"manifest": str(manifest_path)},
        )
    members = manifest.get("members")
    if not isinstance(members, list):
        raise BackendFailedError("The backend manifest member inventory is invalid.")
    observed_members: set[str] = set()
    artifacts: list[dict[str, Any]] = []
    for index, member in enumerate(members):
        if not isinstance(member, Mapping):
            raise BackendFailedError(
                "A backend manifest member is not an object.",
                details={"member_index": index},
            )
        name = member.get("path")
        if name not in expected_members or name in observed_members:
            raise BackendFailedError(
                "A backend manifest member path is invalid.",
                details={"member_index": index, "path": name},
            )
        observed_members.add(name)
        digest, size = _hash_and_fsync_output(result_stage / name)
        if digest != member.get("sha256") or size != member.get("size_bytes"):
            raise BackendFailedError(
                "A backend output does not match its manifest.",
                details={"path": name},
            )
        artifacts.append(
            {
                "kind": "generated_analysis_artifact",
                "role": name.removesuffix(".json").removesuffix(".tsv"),
                "relative_path": name,
                "sha256": digest,
                "size_bytes": size,
                "media_type": (
                    "application/json"
                    if name.endswith(".json")
                    else "text/tab-separated-values"
                ),
                "schema_version": expected_schema_version,
            }
        )
    if observed_members != expected_members:
        raise BackendFailedError(
            "The backend manifest omits required outputs.",
            details={"missing_members": sorted(expected_members - observed_members)},
        )
    _hash_and_fsync_output(manifest_path)
    artifacts.append(
        {
            "kind": "generated_analysis_artifact",
            "role": "backend_manifest",
            "relative_path": "backend_manifest.json",
            "sha256": manifest_digest,
            "size_bytes": manifest_size,
            "media_type": "application/json",
            "schema_version": expected_schema_version,
        }
    )
    analysis, _, _ = _capture_output_json(
        result_stage / "analysis.json", role="edger_analysis"
    )
    if analysis.get("schema_version") != expected_schema_version:
        raise BackendFailedError(
            "The analysis document schema does not match the requested backend schema.",
            details={
                "observed_schema_version": analysis.get("schema_version"),
                "expected_schema_version": expected_schema_version,
            },
        )
    _fsync_directory(result_stage)
    return analysis, sorted(artifacts, key=lambda item: item["relative_path"])


def _verify_complete_public_stage(result_stage: Path) -> None:
    """Run the public fail-closed reader before a C1 bundle is published."""

    try:
        from .run_summary import summarize_run

        summary = summarize_run(result_stage)
    except Exception as error:
        details: dict[str, Any] = {
            "reason": "staged_bundle_verification_failed",
            "cause_type": type(error).__name__,
        }
        if isinstance(error, ToolkitError):
            details.update(
                {
                    "cause_code": error.code.value,
                    "cause_details": dict(error.details),
                }
            )
        raise BackendFailedError(
            "The complete staged result bundle failed independent verification.",
            details=details,
            cause=error,
        ) from error
    if not isinstance(summary, Mapping) or summary.get("status") != "verified_complete":
        raise BackendFailedError(
            "The complete staged result bundle failed independent verification.",
            details={
                "reason": "staged_bundle_verification_incomplete",
                "observed_status": (
                    summary.get("status") if isinstance(summary, Mapping) else None
                ),
            },
        )


def _publish_noreplace(source: Path, target: Path) -> None:
    if not sys_platform_linux():
        raise OutputWriteError(
            "Atomic no-replace backend publication requires the locked Linux runtime.",
            path=target,
            operation="publish_backend_results_noreplace",
            cause=OSError(errno.ENOTSUP, "renameat2 unavailable"),
        )
    library = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(library, "renameat2", None)
    if renameat2 is None:
        raise OutputWriteError(
            "Atomic no-replace backend publication is unavailable.",
            path=target,
            operation="publish_backend_results_noreplace",
            cause=OSError(errno.ENOSYS, "renameat2 unavailable"),
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
            "The backend output path appeared during execution; nothing was overwritten.",
            details={"output_dir": str(target)},
        )
    raise OutputWriteError(
        "The verified backend result bundle could not be published.",
        path=target,
        operation="publish_backend_results_noreplace",
        cause=OSError(error_number, os.strerror(error_number), str(target)),
    )


def sys_platform_linux() -> bool:
    return os.name == "posix" and Path("/proc/sys/kernel").exists()


def _execute_document(
    document: dict[str, Any],
    target: Path,
    workspace: Path,
    *,
    rscript: str | Path,
    r_library: str | Path | None,
    display_configuration: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    expected_schema_version = str(document.get("schema_version", ""))
    if expected_schema_version not in {
        BACKEND_OUTPUT_SCHEMA_VERSION,
        PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION,
    }:
        raise InvalidRequestError(
            "The private backend request schema is unsupported.",
            details={"schema_version": expected_schema_version},
        )
    has_gene_sets = "gene_sets" in document
    if has_gene_sets != (
        expected_schema_version == PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION
    ):
        raise InvalidRequestError(
            "The private backend schema and pathway request presence disagree.",
            details={
                "schema_version": expected_schema_version,
                "gene_sets_present": has_gene_sets,
            },
        )
    request_path = workspace / "backend_request.json"
    result_stage = workspace / "results"
    display_logcpm = workspace / "display-logcpm.tsv"
    if display_configuration is not None:
        document["display_export"] = {"path": str(display_logcpm)}
    _write_private_json(request_path, document)
    response, stderr = _invoke_r(
        request_path,
        result_stage,
        rscript=rscript,
        r_library=r_library,
        expected_schema_version=expected_schema_version,
    )
    execution_scope = document["execution_scope"]
    analysis, artifacts = _verify_result_stage(
        result_stage,
        execution_scope=execution_scope,
        expected_schema_version=expected_schema_version,
    )
    response_data = response.get("data")
    if not isinstance(response_data, Mapping):
        raise BackendFailedError("The successful R response lacks backend data.")
    if response_data.get("runtime_identity") != _EXPECTED_RUNTIME:
        raise BackendFailedError(
            "The R response does not carry the exact locked runtime identity."
        )
    if display_configuration is not None:
        from .display_bundle import build_display_bundle

        display_artifacts = build_display_bundle(
            display_dir=result_stage / "display",
            logcpm_path=display_logcpm,
            core_dir=result_stage,
            core_artifacts=[
                item
                for item in artifacts
                if item.get("relative_path")
                in (_OUTPUT_MEMBERS | {"backend_manifest.json"})
            ],
            backend_document=document,
            backend_data=response_data,
            configuration=display_configuration,
        )
        artifacts.extend(display_artifacts)
        _fsync_directory(result_stage)
        _verify_complete_public_stage(result_stage)
    _publish_noreplace(result_stage, target)
    publication_status = "durability_confirmed"
    warnings: list[dict[str, Any]] = []
    try:
        _fsync_directory(target.parent)
    except OSError as error:
        publication_status = "published_durability_unconfirmed"
        warnings.append(
            {
                "code": "OUTPUT_DURABILITY_UNCONFIRMED",
                "severity": "high",
                "message": (
                    "The complete backend bundle is visible, but parent-directory "
                    "synchronization failed."
                ),
                "details": {
                    "output_dir": str(target),
                    "cause_type": type(error).__name__,
                },
            }
        )
    for artifact in artifacts:
        artifact["path"] = str(target / artifact.pop("relative_path"))
    return {
        "schema_version": expected_schema_version,
        "status": "success",
        "backend": BACKEND_NAME,
        "execution_scope": execution_scope,
        "output_dir": str(target),
        "publication_status": publication_status,
        "data": dict(response_data),
        "analysis": analysis,
        "warnings": warnings,
        "errors": [],
        "artifacts": artifacts,
        "backend_stderr": stderr,
    }


def run_edger_ql(
    request_path: str | Path,
    output_dir: str | Path,
    *,
    rscript: str | Path = "Rscript",
    r_library: str | Path | None = None,
) -> dict[str, Any]:
    """Run the certified-input edgeR QL route through an independent R process.

    The public production adapter accepts only a complete, eligible validation
    bundle.  Every byte consumed by R is a private copy whose digest was checked
    while that same source stream was copied.
    """

    validated = load_analysis_request(request_path)
    target = _resolve_output(output_dir)
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        workspace = Path(
            tempfile.mkdtemp(prefix=f".{target.name}.edger-", dir=target.parent)
        )
        workspace.chmod(0o700)
    except OSError as error:
        raise OutputWriteError(
            "A private edgeR workspace could not be created.",
            path=target.parent,
            operation="create_backend_workspace",
            cause=error,
        ) from error
    try:
        input_root = workspace / "inputs"
        input_root.mkdir(mode=0o700)
        rewritten_input, snapshots, _ = _materialize_validated_inputs(
            validated, input_root
        )
        document = validated.to_backend_document()
        if validated.gene_sets is not None:
            rewritten_gene_sets, pathway_snapshots = _materialize_gene_set_inputs(
                validated.gene_sets, input_root
            )
            document["schema_version"] = PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION
            document["gene_sets"] = rewritten_gene_sets
            snapshots.extend(pathway_snapshots)
        document.update(
            {
                "execution_scope": "validated_p0_input",
                "input": rewritten_input,
                "input_evidence": {
                    "kind": "validated_input_bundle",
                    "plan_id": validated.plan_id,
                    "bundle_path": str(validated.bundle_path),
                    "validation_bundle_artifacts": [
                        dict(item) for item in validated.validation_bundle_artifacts
                    ],
                    "r_input_snapshots": snapshots,
                    "digest_coupling": "private_copy_and_hash_same_source_stream",
                },
            }
        )
        result = _execute_document(
            document,
            target,
            workspace,
            rscript=rscript,
            r_library=r_library,
            display_configuration=validated.display,
        )
        result["plan_id"] = validated.plan_id
        result["analysis_request_sha256"] = validated.request_sha256
        return result
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def _run_edger_ql_benchmark_kernel(
    counts_path: str | Path,
    metadata_path: str | Path,
    output_dir: str | Path,
    *,
    design: Mapping[str, Any],
    contrasts: Sequence[Mapping[str, Any]],
    gene_sets: Mapping[str, Any] | None = None,
    rscript: str | Path = "Rscript",
    r_library: str | Path | None = None,
) -> dict[str, Any]:
    """Run only the statistical kernel for oracle/simulation certification.

    This intentionally private entry point makes no featureCounts provenance
    claim and is not a certified analysis route.  It exists so packaged data
    such as ``airway`` and simulated counts can test the exact same R kernel.
    """

    normalized_design = _parse_design(design)
    normalized_contrasts = _parse_contrasts(list(contrasts))
    normalized_gene_sets: dict[str, Any] | None = None
    if gene_sets is not None:
        from .gene_sets import parse_gene_sets_request

        normalized_gene_sets = parse_gene_sets_request(
            gene_sets,
            request_path=Path.cwd() / "private-benchmark-request.json",
        )
    target = _resolve_output(output_dir)
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        workspace = Path(
            tempfile.mkdtemp(prefix=f".{target.name}.kernel-", dir=target.parent)
        )
        workspace.chmod(0o700)
    except OSError as error:
        raise OutputWriteError(
            "A private benchmark-kernel workspace could not be created.",
            path=target.parent,
            operation="create_backend_workspace",
            cause=error,
        ) from error
    try:
        input_root = workspace / "inputs"
        input_root.mkdir(mode=0o700)
        counts_record, counts_content = _capture_declared_benchmark_file(
            counts_path,
            input_root / "counts.tsv",
            role="benchmark.counts",
        )
        metadata_record, metadata_content = _capture_declared_benchmark_file(
            metadata_path,
            input_root / "metadata.tsv",
            role="benchmark.metadata",
        )
        counts_record["private_relative_path"] = "counts.tsv"
        metadata_record["private_relative_path"] = "metadata.tsv"
        sample_order = _parse_benchmark_counts_header(counts_content)
        metadata_values = _parse_metadata_snapshot(
            metadata_content,
            sample_id_column="sample_id",
            sample_order=sample_order,
            design_terms=normalized_design["terms"],
        )
        gene_set_identity = None
        if normalized_gene_sets is not None:
            gene_set_identity = {
                "gmt": {
                    key: value
                    for key, value in normalized_gene_sets["gmt"].items()
                    if key not in {"path", "declared_path"}
                },
                "annotation": {
                    key: value
                    for key, value in normalized_gene_sets["annotation"].items()
                    if key not in {"path", "declared_path"}
                },
                "minimum_tested_genes": normalized_gene_sets[
                    "minimum_tested_genes"
                ],
            }
        invocation_identity = hashlib.sha256(
            _canonical_bytes(
                {
                    "counts_sha256": counts_record["sha256"],
                    "metadata_sha256": metadata_record["sha256"],
                    "design": normalized_design,
                    "contrasts": list(normalized_contrasts),
                    "gene_sets": gene_set_identity,
                }
            )
        ).hexdigest()
        document = {
            "schema_version": (
                PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION
                if normalized_gene_sets is not None
                else ANALYSIS_SCHEMA_VERSION
            ),
            "kind": "edger_ql_backend_request",
            "execution_scope": "backend_kernel_only",
            "analysis_request": {
                "kind": "private_benchmark_kernel_invocation",
                "sha256": invocation_identity,
            },
            "input_evidence": {
                "kind": "benchmark_kernel_inputs",
                "benchmark_id": invocation_identity,
                "certified_analysis_input": False,
                "r_input_snapshots": [counts_record, metadata_record],
                "digest_coupling": "private_copy_and_hash_same_source_stream",
            },
            "input": {
                "input_semantics": "benchmark_kernel_integer_counts",
                "sample_order": sample_order,
                "metadata": {
                    "path": str(input_root / "metadata.tsv"),
                    "sample_id_column": "sample_id",
                },
                "metadata_values": metadata_values,
                "gene_id_policy": {
                    "internal_key": "gene_id",
                    "symbol_role": "display_only",
                    "strip_version": False,
                },
                "benchmark_counts": {"matrix_path": str(input_root / "counts.tsv")},
            },
            "design": normalized_design,
            "contrasts": [dict(item) for item in normalized_contrasts],
        }
        if normalized_gene_sets is not None:
            rewritten_gene_sets, pathway_snapshots = _materialize_gene_set_inputs(
                normalized_gene_sets, input_root
            )
            document["gene_sets"] = rewritten_gene_sets
            document["input_evidence"]["r_input_snapshots"].extend(
                pathway_snapshots
            )
        result = _execute_document(
            document,
            target,
            workspace,
            rscript=rscript,
            r_library=r_library,
        )
        result["benchmark_id"] = invocation_identity
        return result
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


__all__ = ["BACKEND_NAME", "run_edger_ql"]
