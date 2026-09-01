"""Strict, read-only validation for the three P0 input semantics.

The public functions in this module do not write analysis artifacts.  They
return plain JSON-compatible dictionaries that can be placed directly into the
CLI response envelope.
"""

from __future__ import annotations

import csv
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import io
import math
from pathlib import Path, PurePosixPath
import re
import struct
from typing import Any, Mapping, Sequence
import zlib

from .errors import (
    AssayProtocolRequiredError,
    CountValuesInvalidError,
    GeneIdentifierError,
    InputEvidenceRequiredError,
    InputIntegrityError,
    InputReadError,
    InputValidationError,
    InvalidRequestError,
    SalmonOffsetRequiredError,
    SampleSetMismatchError,
)
from .provenance import (
    ArtifactInventory,
    has_control_characters,
    read_json_object,
    require_directory,
    require_expected_keys,
    require_nonempty_string,
    resolve_declared_path,
)

REQUEST_SCHEMA_VERSION = "1.0"
MAX_EXACT_INTEGER_COUNT = (2**53) - 1
FEATURECOUNTS_INTEGER = "featurecounts_integer"
SALMON_FULL_LENGTH = "salmon_quant_dirs_full_length"
SALMON_THREE_PRIME = "salmon_quant_dirs_three_prime"
SUPPORTED_INPUT_SEMANTICS = (
    FEATURECOUNTS_INTEGER,
    SALMON_FULL_LENGTH,
    SALMON_THREE_PRIME,
)

_TOP_LEVEL_ALLOWED = {
    "schema_version",
    "input_semantics",
    "metadata",
    "producer",
    "reference",
    "gene_id",
    "samples",
    "featurecounts",
    "salmon",
    "assay_protocol",
    "analysis_options",
}
_TOP_LEVEL_REQUIRED = {
    "schema_version",
    "input_semantics",
    "metadata",
    "producer",
    "reference",
    "gene_id",
    "samples",
}
_FEATURECOUNTS_ANNOTATION_COLUMNS = (
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length",
)
_INTEGER_COUNT_PATTERN = re.compile(r"^[0-9]+$")
_DECIMAL_PATTERN = re.compile(
    r"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
)
_GENE_VERSION_PATTERN = re.compile(r"\.([0-9]+)$")
_SALMON_INDEX_HASH_LENGTHS = {
    "index_seq_hash": 64,
    "index_name_hash": 64,
    "index_decoy_seq_hash": 64,
    "index_decoy_name_hash": 64,
    "index_seq_hash512": 128,
    "index_name_hash512": 128,
    "index_decoy_seq_hash512": 128,
    "index_decoy_name_hash512": 128,
}
_FEATURECOUNTS_VERSION_PATTERN = re.compile(
    r"Program\s*:\s*featureCounts\s+v?([^;\s]+)", re.IGNORECASE
)


def inspect_request(path: str | Path) -> dict[str, Any]:
    """Inspect declared inputs, path resolution, and provenance without writes."""

    return _process_request(path, validation_level="inspect")


def validate_request(path: str | Path) -> dict[str, Any]:
    """Validate all P0 input-semantic invariants without running statistics."""

    return _process_request(path, validation_level="validate")


def _process_request(path: str | Path, *, validation_level: str) -> dict[str, Any]:
    request_path = _resolve_request_path(path)
    inventory = ArtifactInventory()
    request_fingerprint = inventory.add(
        request_path,
        role="request",
        declared_path=str(path),
        capture_for_parse=True,
    )
    document = read_json_object(
        request_path,
        document_role="analysis_request",
        content=inventory.consume_bytes(request_path),
    )
    require_expected_keys(
        document,
        allowed=_TOP_LEVEL_ALLOWED,
        required=_TOP_LEVEL_REQUIRED,
        context="analysis request",
    )

    if document["schema_version"] != REQUEST_SCHEMA_VERSION:
        raise InvalidRequestError(
            "The request schema version is not supported.",
            details={
                "observed": document["schema_version"],
                "supported": REQUEST_SCHEMA_VERSION,
            },
        )
    semantics = require_nonempty_string(
        document["input_semantics"], field="input_semantics"
    )
    if semantics not in SUPPORTED_INPUT_SEMANTICS:
        raise InvalidRequestError(
            "The declared input semantics are not supported.",
            details={
                "input_semantics": semantics,
                "supported": list(SUPPORTED_INPUT_SEMANTICS),
            },
        )

    producer = _parse_producer(document["producer"], semantics=semantics)
    reference = _parse_reference(
        document["reference"],
        request_path=request_path,
        inventory=inventory,
        role="reference.source",
    )
    gene_policy = _parse_gene_policy(document["gene_id"])
    sample_entries = _parse_sample_entries(document["samples"], semantics=semantics)
    sample_order = [sample["sample_id"] for sample in sample_entries]
    warnings: list[dict[str, Any]] = []
    metadata = _parse_metadata(
        document["metadata"],
        request_path=request_path,
        expected_sample_order=sample_order,
        inventory=inventory,
        warnings=warnings,
    )

    common_data: dict[str, Any] = {
        "validation_level": validation_level,
        "input_semantics": semantics,
        "request": {
            "schema_version": REQUEST_SCHEMA_VERSION,
            "path": str(request_path),
            "sha256": request_fingerprint.sha256,
        },
        "sample_order": sample_order,
        "sample_count": len(sample_order),
        "metadata": metadata,
        "producer": producer,
        "reference": reference,
        "gene_id_policy": {
            "internal_key": "gene_id",
            "symbol_role": "display_only",
            "strip_version": gene_policy["strip_version"],
        },
    }

    if semantics == FEATURECOUNTS_INTEGER:
        mode_data, certification_eligible = _process_featurecounts(
            document,
            request_path=request_path,
            sample_entries=sample_entries,
            producer=producer,
            reference=reference,
            strip_gene_version=gene_policy["strip_version"],
            inventory=inventory,
            validation_level=validation_level,
        )
    else:
        reference["salmon_index_linkage_status"] = (
            "declared_source_recorded_separately_from_salmon_index_hash"
        )
        mode_data, certification_eligible = _process_salmon(
            document,
            semantics=semantics,
            request_path=request_path,
            sample_entries=sample_entries,
            producer=producer,
            strip_gene_version=gene_policy["strip_version"],
            inventory=inventory,
            warnings=warnings,
            validation_level=validation_level,
        )

    common_data.update(mode_data)
    common_data["certification_scope"] = "input_semantics_only"
    high_risk_reasons = [
        warning["code"] for warning in warnings if warning.get("severity") == "high"
    ]
    if validation_level == "inspect" and certification_eligible:
        common_data["input_certification_eligible"] = None
        common_data["input_certification_status"] = "not_evaluated"
        common_data["input_certification_reasons"] = ["NUMERIC_DOMAIN_NOT_EVALUATED"]
    else:
        common_data["input_certification_eligible"] = certification_eligible
        common_data["input_certification_status"] = (
            "passed" if certification_eligible else "ineligible"
        )
        common_data["input_certification_reasons"] = high_risk_reasons
    inventory.verify_unchanged()
    return {
        "data": common_data,
        "warnings": warnings,
        "artifacts": inventory.records(),
    }


def _resolve_request_path(path: str | Path) -> Path:
    if not isinstance(path, (str, Path)) or not str(path):
        raise InvalidRequestError(
            "The request path must be a non-empty string or Path.",
            details={"path_type": type(path).__name__},
        )
    if has_control_characters(str(path)):
        raise InvalidRequestError(
            "The request path must not contain Unicode control characters."
        )
    candidate = Path(path)
    try:
        resolved = candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise InputReadError(
            "The request file could not be resolved.",
            path=candidate,
            operation="resolve_request",
            cause=error,
        ) from error
    if not resolved.is_file():
        error = IsADirectoryError(str(resolved))
        raise InputReadError(
            "The request path is not a regular file.",
            path=resolved,
            operation="read_request",
            cause=error,
        )
    return resolved


def _parse_producer(value: Any, *, semantics: str) -> dict[str, str]:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'producer' must be an object.")
    require_expected_keys(
        value,
        allowed={"name", "version"},
        required={"name", "version"},
        context="producer",
    )
    name = require_nonempty_string(value["name"], field="producer.name")
    version = require_nonempty_string(value["version"], field="producer.version")
    expected_name = "featureCounts" if semantics == FEATURECOUNTS_INTEGER else "Salmon"
    if name != expected_name:
        raise InputIntegrityError(
            "The declared producer is incompatible with the input semantics.",
            details={"expected": expected_name, "observed": name},
        )
    return {"name": name, "version": version}


def _parse_reference(
    value: Any,
    *,
    request_path: Path,
    inventory: ArtifactInventory,
    role: str,
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'reference' must be an object.")
    require_expected_keys(
        value,
        allowed={"name", "version", "source", "sha256"},
        required={"name", "version", "source"},
        context="reference",
    )
    name = require_nonempty_string(value["name"], field="reference.name")
    version = require_nonempty_string(value["version"], field="reference.version")
    declared_source = require_nonempty_string(value["source"], field="reference.source")
    source_path = resolve_declared_path(declared_source, relative_to=request_path)
    fingerprint = inventory.add(
        source_path,
        role=role,
        declared_path=declared_source,
        expected_sha256=value.get("sha256"),
    )
    return {
        "name": name,
        "version": version,
        "source_path": fingerprint.path,
        "source_sha256": fingerprint.sha256,
    }


def _parse_gene_policy(value: Any) -> dict[str, bool]:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'gene_id' must be an object.")
    require_expected_keys(
        value,
        allowed={"strip_version"},
        required={"strip_version"},
        context="gene_id",
    )
    strip_version = value["strip_version"]
    if not isinstance(strip_version, bool):
        raise InvalidRequestError(
            "'gene_id.strip_version' must be a JSON boolean.",
            details={"observed_type": type(strip_version).__name__},
        )
    return {"strip_version": strip_version}


def _parse_sample_entries(value: Any, *, semantics: str) -> list[dict[str, str]]:
    if not isinstance(value, list) or not value:
        raise InvalidRequestError("'samples' must be a non-empty array.")
    if semantics == FEATURECOUNTS_INTEGER:
        allowed = {"sample_id", "counts_file", "count_column"}
    else:
        allowed = {"sample_id", "quant_dir"}

    entries: list[dict[str, str]] = []
    seen: set[str] = set()
    for index, raw_entry in enumerate(value):
        if not isinstance(raw_entry, Mapping):
            raise InvalidRequestError(
                "Every sample entry must be an object.",
                details={"sample_index": index},
            )
        require_expected_keys(
            raw_entry,
            allowed=allowed,
            required={"sample_id"},
            context=f"samples[{index}]",
        )
        sample_id = _require_identifier(raw_entry["sample_id"], field="sample_id")
        if sample_id in seen:
            raise SampleSetMismatchError(
                "Sample identifiers must be unique.",
                details={"duplicate_sample_id": sample_id},
            )
        seen.add(sample_id)
        entry = {"sample_id": sample_id}
        for key in allowed - {"sample_id"}:
            if key in raw_entry:
                entry[key] = require_nonempty_string(
                    raw_entry[key], field=f"samples[{index}].{key}"
                )
        entries.append(entry)
    return entries


def _require_identifier(value: Any, *, field: str) -> str:
    identifier = require_nonempty_string(value, field=field)
    if identifier != value or any(character.isspace() for character in identifier):
        raise InvalidRequestError(
            f"'{field}' must not contain leading, trailing, or embedded whitespace.",
            details={"field": field, "value": value},
        )
    return identifier


def _parse_metadata(
    value: Any,
    *,
    request_path: Path,
    expected_sample_order: Sequence[str],
    inventory: ArtifactInventory,
    warnings: list[dict[str, Any]],
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'metadata' must be an object.")
    require_expected_keys(
        value,
        allowed={"path", "sample_id_column"},
        required={"path", "sample_id_column"},
        context="metadata",
    )
    declared_path = require_nonempty_string(value["path"], field="metadata.path")
    id_column = require_nonempty_string(
        value["sample_id_column"], field="metadata.sample_id_column"
    )
    metadata_path = resolve_declared_path(declared_path, relative_to=request_path)
    fingerprint = inventory.add(
        metadata_path,
        role="metadata",
        declared_path=declared_path,
        capture_for_parse=True,
    )
    header, rows = _read_tsv(
        metadata_path,
        role="metadata",
        content=inventory.consume_bytes(metadata_path),
    )
    _require_unique_header(header, path=metadata_path)
    if id_column not in header:
        raise SampleSetMismatchError(
            "The metadata sample identifier column is absent.",
            details={"path": str(metadata_path), "sample_id_column": id_column},
        )
    id_index = header.index(id_column)
    observed_order: list[str] = []
    seen: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        _require_row_width(row, header, path=metadata_path, row_number=row_number)
        sample_id = _require_identifier(
            row[id_index], field=f"metadata row {row_number} sample ID"
        )
        if sample_id in seen:
            raise SampleSetMismatchError(
                "Metadata sample identifiers must be unique.",
                details={"path": str(metadata_path), "duplicate_sample_id": sample_id},
            )
        if any(cell == "" for cell in row):
            raise InputValidationError(
                "Metadata rows must not contain empty fields.",
                details={"path": str(metadata_path), "row": row_number},
            )
        seen.add(sample_id)
        observed_order.append(sample_id)

    expected_set = set(expected_sample_order)
    observed_set = set(observed_order)
    if expected_set != observed_set:
        raise SampleSetMismatchError(
            "Metadata and assay inputs must contain exactly the same sample set.",
            details={
                "missing_from_metadata": sorted(expected_set - observed_set),
                "unexpected_in_metadata": sorted(observed_set - expected_set),
                "assay_sample_order": list(expected_sample_order),
                "metadata_sample_order": observed_order,
            },
        )
    order_matches = observed_order == list(expected_sample_order)
    if not order_matches:
        warnings.append(
            {
                "code": "METADATA_ORDER_NORMALIZED",
                "severity": "warning",
                "message": "Metadata rows will be aligned to the declared assay sample order.",
                "details": {
                    "metadata_sample_order": observed_order,
                    "analysis_sample_order": list(expected_sample_order),
                },
            }
        )
    return {
        "path": fingerprint.path,
        "sha256": fingerprint.sha256,
        "sample_id_column": id_column,
        "columns": header,
        "row_count": len(rows),
        "sample_order": observed_order,
        "order_matches_analysis": order_matches,
    }


def _process_featurecounts(
    document: Mapping[str, Any],
    *,
    request_path: Path,
    sample_entries: Sequence[Mapping[str, str]],
    producer: Mapping[str, str],
    reference: Mapping[str, Any],
    strip_gene_version: bool,
    inventory: ArtifactInventory,
    validation_level: str,
) -> tuple[dict[str, Any], bool]:
    if "featurecounts" not in document:
        raise InputEvidenceRequiredError(
            "A featureCounts request requires an explicit input layout declaration."
        )
    prohibited = [
        key
        for key in ("salmon", "assay_protocol", "analysis_options")
        if key in document
    ]
    if prohibited:
        raise InvalidRequestError(
            "A featureCounts request contains fields for another input route.",
            details={"prohibited_keys": prohibited},
        )
    config = document["featurecounts"]
    if not isinstance(config, Mapping):
        raise InvalidRequestError("'featurecounts' must be an object.")
    require_expected_keys(
        config,
        allowed={"layout", "matrix", "manifest"},
        required={"layout"},
        context="featurecounts",
    )
    layout = require_nonempty_string(config["layout"], field="featurecounts.layout")
    if layout == "combined_matrix":
        required = {"layout", "matrix", "manifest"}
        missing_evidence = required - set(config)
        if missing_evidence:
            raise InputEvidenceRequiredError(
                "A combined featureCounts matrix requires its typed manifest and matrix path.",
                details={"missing_keys": sorted(missing_evidence)},
            )
        require_expected_keys(
            config,
            allowed=required,
            required=required,
            context="featurecounts combined_matrix",
        )
        if any(set(sample) != {"sample_id"} for sample in sample_entries):
            raise InvalidRequestError(
                "Combined-matrix sample entries may contain only 'sample_id'."
            )
        result = _validate_featurecounts_matrix(
            config,
            request_path=request_path,
            sample_order=[sample["sample_id"] for sample in sample_entries],
            producer=producer,
            reference=reference,
            strip_gene_version=strip_gene_version,
            inventory=inventory,
            validation_level=validation_level,
        )
    elif layout == "per_sample_files":
        require_expected_keys(
            config,
            allowed={"layout"},
            required={"layout"},
            context="featurecounts per_sample_files",
        )
        for sample in sample_entries:
            missing = {"counts_file", "count_column"} - set(sample)
            if missing:
                raise InputEvidenceRequiredError(
                    "Every per-sample featureCounts entry requires its original "
                    "file and count column.",
                    details={
                        "sample_id": sample["sample_id"],
                        "missing_keys": sorted(missing),
                    },
                )
        result = _validate_featurecounts_sample_files(
            sample_entries,
            request_path=request_path,
            producer=producer,
            strip_gene_version=strip_gene_version,
            inventory=inventory,
            validation_level=validation_level,
        )
    else:
        raise InvalidRequestError(
            "The featureCounts layout is not supported.",
            details={
                "observed": layout,
                "supported": ["combined_matrix", "per_sample_files"],
            },
        )

    return (
        {
            "featurecounts": result,
            "route": {
                "backend_input": "edgeR::DGEList",
                "count_semantics": "integer",
                "maximum_exact_integer": MAX_EXACT_INTEGER_COUNT,
                "transcript_length_offset": False,
            },
        },
        True,
    )


def _validate_featurecounts_matrix(
    config: Mapping[str, Any],
    *,
    request_path: Path,
    sample_order: list[str],
    producer: Mapping[str, str],
    reference: Mapping[str, Any],
    strip_gene_version: bool,
    inventory: ArtifactInventory,
    validation_level: str,
) -> dict[str, Any]:
    matrix_declared = require_nonempty_string(
        config["matrix"], field="featurecounts.matrix"
    )
    manifest_declared = require_nonempty_string(
        config["manifest"], field="featurecounts.manifest"
    )
    matrix_path = resolve_declared_path(matrix_declared, relative_to=request_path)
    manifest_path = resolve_declared_path(manifest_declared, relative_to=request_path)
    manifest_fingerprint = inventory.add(
        manifest_path,
        role="featurecounts.manifest",
        declared_path=manifest_declared,
        capture_for_parse=True,
    )
    manifest = read_json_object(
        manifest_path,
        document_role="featurecounts_manifest",
        content=inventory.consume_bytes(manifest_path),
    )
    require_expected_keys(
        manifest,
        allowed={
            "schema_version",
            "artifact_type",
            "artifact",
            "gene_id_column",
            "display_columns",
            "sample_columns",
            "producer",
            "reference",
        },
        required={
            "schema_version",
            "artifact_type",
            "artifact",
            "gene_id_column",
            "sample_columns",
            "producer",
            "reference",
        },
        context="featureCounts manifest",
    )
    if manifest["schema_version"] != REQUEST_SCHEMA_VERSION:
        raise InputIntegrityError(
            "The featureCounts manifest schema version is incompatible.",
            details={"observed": manifest["schema_version"]},
        )
    if manifest["artifact_type"] != "featurecounts_integer_matrix":
        raise InputEvidenceRequiredError(
            "The manifest does not declare a featureCounts integer matrix.",
            details={"artifact_type": manifest["artifact_type"]},
        )
    _compare_manifest_producer(manifest["producer"], producer)
    _compare_manifest_reference(
        manifest["reference"],
        reference,
        manifest_path=manifest_path,
        inventory=inventory,
    )

    artifact = manifest["artifact"]
    if not isinstance(artifact, Mapping):
        raise InvalidRequestError("'manifest.artifact' must be an object.")
    require_expected_keys(
        artifact,
        allowed={"path", "sha256"},
        required={"path", "sha256"},
        context="manifest.artifact",
    )
    manifest_matrix_declared = require_nonempty_string(
        artifact["path"], field="manifest.artifact.path"
    )
    manifest_matrix_path = resolve_declared_path(
        manifest_matrix_declared, relative_to=manifest_path
    )
    if manifest_matrix_path != matrix_path:
        raise InputIntegrityError(
            "The request and manifest identify different count matrices.",
            details={
                "request_matrix": str(matrix_path),
                "manifest_matrix": str(manifest_matrix_path),
            },
        )
    matrix_fingerprint = inventory.add(
        matrix_path,
        role="featurecounts.matrix",
        declared_path=matrix_declared,
        expected_sha256=artifact["sha256"],
        capture_for_parse=True,
    )

    gene_id_column = require_nonempty_string(
        manifest["gene_id_column"], field="manifest.gene_id_column"
    )
    if gene_id_column.lower() in {"gene_name", "symbol", "gene_symbol"}:
        raise GeneIdentifierError(
            "A display symbol cannot be used as the internal gene key.",
            details={"gene_id_column": gene_id_column},
        )
    display_columns = manifest.get("display_columns", [])
    if not isinstance(display_columns, list) or any(
        not isinstance(item, str) or not item for item in display_columns
    ):
        raise InvalidRequestError(
            "'manifest.display_columns' must be an array of strings."
        )
    display_columns = [
        require_nonempty_string(item, field=f"manifest.display_columns[{index}]")
        for index, item in enumerate(display_columns)
    ]
    if len(set(display_columns)) != len(display_columns):
        raise InvalidRequestError("'manifest.display_columns' must be unique.")
    declared_samples = manifest["sample_columns"]
    if not isinstance(declared_samples, list) or any(
        not isinstance(item, str) or not item for item in declared_samples
    ):
        raise InvalidRequestError(
            "'manifest.sample_columns' must be an array of strings."
        )
    declared_samples = [
        _require_identifier(item, field=f"manifest.sample_columns[{index}]")
        for index, item in enumerate(declared_samples)
    ]
    if declared_samples != sample_order:
        raise SampleSetMismatchError(
            "The manifest sample columns do not exactly match the request order.",
            details={"request_order": sample_order, "manifest_order": declared_samples},
        )

    header, rows = _read_tsv(
        matrix_path,
        role="featurecounts_matrix",
        content=inventory.consume_bytes(matrix_path),
    )
    expected_header = [gene_id_column, *display_columns, *sample_order]
    _require_unique_header(expected_header, path=manifest_path)
    _require_unique_header(header, path=matrix_path)
    if header != expected_header:
        raise InputIntegrityError(
            "The matrix header does not exactly match its typed manifest.",
            details={"expected_header": expected_header, "observed_header": header},
        )
    gene_ids: list[str] = []
    count_start = 1 + len(display_columns)
    for row_number, row in enumerate(rows, start=2):
        _require_row_width(row, header, path=matrix_path, row_number=row_number)
        gene_ids.append(row[0])
        if validation_level == "validate":
            for column_index, sample_id in enumerate(sample_order, start=count_start):
                _validate_integer_count(
                    row[column_index],
                    path=matrix_path,
                    row_number=row_number,
                    sample_id=sample_id,
                )
    normalized = _normalize_gene_ids(
        gene_ids,
        strip_version=strip_gene_version,
        context="featureCounts matrix",
        allow_repeated_original=False,
    )
    return {
        "layout": "combined_matrix",
        "matrix_path": matrix_fingerprint.path,
        "matrix_sha256": matrix_fingerprint.sha256,
        "manifest_path": manifest_fingerprint.path,
        "manifest_sha256": manifest_fingerprint.sha256,
        "gene_id_column": gene_id_column,
        "display_columns": display_columns,
        "sample_columns": sample_order,
        "gene_count": len(normalized),
        "gene_id_inventory_sha256": _sequence_digest(normalized),
    }


def _compare_manifest_producer(value: Any, expected: Mapping[str, str]) -> None:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'manifest.producer' must be an object.")
    require_expected_keys(
        value,
        allowed={"name", "version"},
        required={"name", "version"},
        context="manifest.producer",
    )
    observed = {
        "name": require_nonempty_string(value["name"], field="manifest.producer.name"),
        "version": require_nonempty_string(
            value["version"], field="manifest.producer.version"
        ),
    }
    if observed != dict(expected):
        raise InputIntegrityError(
            "The request and manifest producer evidence disagree.",
            details={"request": dict(expected), "manifest": observed},
        )


def _compare_manifest_reference(
    value: Any,
    expected: Mapping[str, Any],
    *,
    manifest_path: Path,
    inventory: ArtifactInventory,
) -> None:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'manifest.reference' must be an object.")
    require_expected_keys(
        value,
        allowed={"name", "version", "source", "sha256"},
        required={"name", "version", "source", "sha256"},
        context="manifest.reference",
    )
    source_declared = require_nonempty_string(
        value["source"], field="manifest.reference.source"
    )
    source_path = resolve_declared_path(source_declared, relative_to=manifest_path)
    fingerprint = inventory.add(
        source_path,
        role="featurecounts.manifest_reference",
        declared_path=source_declared,
        expected_sha256=value["sha256"],
    )
    observed = {
        "name": require_nonempty_string(value["name"], field="manifest.reference.name"),
        "version": require_nonempty_string(
            value["version"], field="manifest.reference.version"
        ),
        "source_path": fingerprint.path,
        "source_sha256": fingerprint.sha256,
    }
    if observed != dict(expected):
        raise InputIntegrityError(
            "The request and manifest reference evidence disagree.",
            details={"request": dict(expected), "manifest": observed},
        )


def _validate_featurecounts_sample_files(
    sample_entries: Sequence[Mapping[str, str]],
    *,
    request_path: Path,
    producer: Mapping[str, str],
    strip_gene_version: bool,
    inventory: ArtifactInventory,
    validation_level: str,
) -> dict[str, Any]:
    expected_gene_ids: list[str] | None = None
    expected_annotations: list[str] | None = None
    files: list[dict[str, Any]] = []
    for sample in sample_entries:
        sample_id = sample["sample_id"]
        declared_path = sample["counts_file"]
        count_column = sample["count_column"]
        counts_path = resolve_declared_path(declared_path, relative_to=request_path)
        fingerprint = inventory.add(
            counts_path,
            role=f"featurecounts.sample.{sample_id}",
            declared_path=declared_path,
            capture_for_parse=True,
        )
        comments, header, rows = _read_featurecounts_tsv(
            counts_path,
            content=inventory.consume_bytes(counts_path),
        )
        observed_version = _featurecounts_version(comments, path=counts_path)
        if observed_version != producer["version"]:
            raise InputIntegrityError(
                "A featureCounts file producer version differs from the request.",
                details={
                    "sample_id": sample_id,
                    "expected_version": producer["version"],
                    "observed_version": observed_version,
                },
            )
        expected_header = [*_FEATURECOUNTS_ANNOTATION_COLUMNS, count_column]
        _require_unique_header(expected_header, path=counts_path)
        _require_unique_header(header, path=counts_path)
        if header != expected_header:
            raise InputIntegrityError(
                "A per-sample featureCounts file must contain exactly one declared count column.",
                details={
                    "sample_id": sample_id,
                    "expected_header": expected_header,
                    "observed_header": header,
                },
            )
        gene_ids: list[str] = []
        annotations: list[str] = []
        for row_number, row in enumerate(rows, start=len(comments) + 2):
            _require_row_width(row, header, path=counts_path, row_number=row_number)
            gene_ids.append(row[0])
            annotations.append("\t".join(row[: len(_FEATURECOUNTS_ANNOTATION_COLUMNS)]))
            if validation_level == "validate":
                _validate_integer_count(
                    row[-1],
                    path=counts_path,
                    row_number=row_number,
                    sample_id=sample_id,
                )
        normalized = _normalize_gene_ids(
            gene_ids,
            strip_version=strip_gene_version,
            context=f"featureCounts sample '{sample_id}'",
            allow_repeated_original=False,
        )
        if expected_gene_ids is None:
            expected_gene_ids = normalized
            expected_annotations = annotations
        elif normalized != expected_gene_ids:
            raise GeneIdentifierError(
                "Per-sample featureCounts files must contain the same genes in the same order.",
                details={"sample_id": sample_id},
            )
        elif annotations != expected_annotations:
            raise InputIntegrityError(
                "Per-sample featureCounts annotation fields must match exactly.",
                details={"sample_id": sample_id},
            )
        files.append(
            {
                "sample_id": sample_id,
                "path": fingerprint.path,
                "sha256": fingerprint.sha256,
                "count_column": count_column,
                "producer_version": observed_version,
            }
        )
    normalized_ids = expected_gene_ids or []
    return {
        "layout": "per_sample_files",
        "files": files,
        "sample_columns": [sample["sample_id"] for sample in sample_entries],
        "gene_count": len(normalized_ids),
        "gene_id_inventory_sha256": _sequence_digest(normalized_ids),
        "annotation_inventory_sha256": _sequence_digest(expected_annotations or []),
    }


def _process_salmon(
    document: Mapping[str, Any],
    *,
    semantics: str,
    request_path: Path,
    sample_entries: Sequence[Mapping[str, str]],
    producer: Mapping[str, str],
    strip_gene_version: bool,
    inventory: ArtifactInventory,
    warnings: list[dict[str, Any]],
    validation_level: str,
) -> tuple[dict[str, Any], bool]:
    if "salmon" not in document:
        raise InputEvidenceRequiredError(
            "A Salmon request requires an explicit tx2gene declaration."
        )
    if "featurecounts" in document:
        raise InvalidRequestError(
            "A Salmon request cannot contain a featureCounts declaration."
        )
    for sample in sample_entries:
        if "quant_dir" not in sample:
            raise InputEvidenceRequiredError(
                "Every Salmon sample requires its original quantification directory.",
                details={"sample_id": sample["sample_id"]},
            )

    if semantics == SALMON_THREE_PRIME:
        if document.get("assay_protocol") != "three_prime":
            raise AssayProtocolRequiredError(
                "Three-prime Salmon inputs require assay_protocol='three_prime'.",
                details={"observed": document.get("assay_protocol")},
            )
        override, certification_eligible = _parse_three_prime_override(
            document.get("analysis_options"), warnings=warnings
        )
    else:
        if "assay_protocol" in document and document["assay_protocol"] != "full_length":
            raise AssayProtocolRequiredError(
                "Full-length Salmon inputs may only declare assay_protocol='full_length'.",
                details={"observed": document["assay_protocol"]},
            )
        if "analysis_options" in document:
            raise InvalidRequestError(
                "Three-prime analysis options are not valid for a full-length request."
            )
        override = None
        certification_eligible = True

    salmon_config = document["salmon"]
    if not isinstance(salmon_config, Mapping):
        raise InvalidRequestError("'salmon' must be an object.")
    require_expected_keys(
        salmon_config,
        allowed={"tx2gene", "tx2gene_sha256"},
        required={"tx2gene"},
        context="salmon",
    )
    tx2gene_declared = require_nonempty_string(
        salmon_config["tx2gene"], field="salmon.tx2gene"
    )
    tx2gene_path = resolve_declared_path(tx2gene_declared, relative_to=request_path)

    sample_records: list[dict[str, Any]] = []
    transcript_order: list[str] | None = None
    replicate_records: list[dict[str, Any]] = []
    index_hashes: dict[str, str] | None = None
    index_identity: str | None = None
    for sample in sample_entries:
        sample_id = sample["sample_id"]
        quant_declared = sample["quant_dir"]
        quant_dir = require_directory(
            resolve_declared_path(quant_declared, relative_to=request_path),
            declared_path=quant_declared,
            role=f"salmon.quant_dir.{sample_id}",
        )
        quant_path = _required_child(quant_dir, "quant.sf", sample_id=sample_id)
        cmd_path = _required_child(quant_dir, "cmd_info.json", sample_id=sample_id)
        cmd_fingerprint = inventory.add(
            cmd_path,
            role=f"salmon.cmd_info.{sample_id}",
            declared_path=f"{quant_declared}/cmd_info.json",
            capture_for_parse=True,
        )
        cmd_info = read_json_object(
            cmd_path,
            document_role="salmon_cmd_info",
            content=inventory.consume_bytes(cmd_path),
        )
        aux_declared, aux_dir = _resolve_salmon_aux_dir(
            quant_dir,
            cmd_info=cmd_info,
            sample_id=sample_id,
        )
        meta_path = _required_child(aux_dir, "meta_info.json", sample_id=sample_id)
        meta_fingerprint = inventory.add(
            meta_path,
            role=f"salmon.meta_info.{sample_id}",
            declared_path=f"{quant_declared}/{aux_declared}/meta_info.json",
            capture_for_parse=True,
        )
        meta_info = read_json_object(
            meta_path,
            document_role="salmon_meta_info",
            content=inventory.consume_bytes(meta_path),
        )
        salmon_metadata = _validate_salmon_metadata(
            cmd_info,
            meta_info,
            sample_id=sample_id,
            expected_version=producer["version"],
        )
        if index_hashes is None:
            index_hashes = salmon_metadata["index_hashes"]
        elif salmon_metadata["index_hashes"] != index_hashes:
            raise InputIntegrityError(
                "Salmon samples were quantified against different index identities.",
                details={"sample_id": sample_id},
            )
        if index_identity is None:
            index_identity = salmon_metadata["index"]
        elif salmon_metadata["index"] != index_identity:
            warnings.append(
                {
                    "code": "SALMON_INDEX_PATHS_DIFFER",
                    "severity": "warning",
                    "message": (
                        "Salmon index paths differ across samples, while their "
                        "required content-hash identity remains equal."
                    ),
                    "details": {
                        "sample_id": sample_id,
                        "first_index_path": index_identity,
                        "observed_index_path": salmon_metadata["index"],
                    },
                },
            )
        quant_fingerprint = inventory.add(
            quant_path,
            role=f"salmon.quant.{sample_id}",
            declared_path=f"{quant_declared}/quant.sf",
            capture_for_parse=True,
        )
        observed_transcripts = _read_salmon_quant(
            quant_path,
            sample_id=sample_id,
            require_offset=semantics == SALMON_FULL_LENGTH,
            validation_level=validation_level,
            content=inventory.consume_bytes(quant_path),
        )
        if transcript_order is None:
            transcript_order = observed_transcripts
        elif observed_transcripts != transcript_order:
            raise InputIntegrityError(
                "Salmon quantifications must contain the same transcripts in the same order.",
                details={"sample_id": sample_id},
            )
        replicate = _inspect_inferential_replicates(
            aux_dir,
            aux_declared=aux_declared,
            cmd_info=cmd_info,
            meta_info=meta_info,
            transcript_order=observed_transcripts,
            sample_id=sample_id,
            quant_declared=quant_declared,
            inventory=inventory,
        )
        replicate_records.append(replicate)
        sample_records.append(
            {
                "sample_id": sample_id,
                "quant_dir": str(quant_dir),
                "aux_dir": str(aux_dir),
                "aux_dir_declared": aux_declared,
                "quant_sha256": quant_fingerprint.sha256,
                "cmd_info_sha256": cmd_fingerprint.sha256,
                "meta_info_sha256": meta_fingerprint.sha256,
                "salmon_version": salmon_metadata["salmon_version"],
                "index": salmon_metadata["index"],
                "index_hashes": salmon_metadata["index_hashes"],
            }
        )

    tx2gene_fingerprint = inventory.add(
        tx2gene_path,
        role="salmon.tx2gene",
        declared_path=tx2gene_declared,
        expected_sha256=salmon_config.get("tx2gene_sha256"),
        capture_for_parse=True,
    )
    transcripts = transcript_order or []
    mapping = _read_tx2gene(
        tx2gene_path,
        strip_gene_version=strip_gene_version,
        content=inventory.consume_bytes(tx2gene_path),
    )
    quant_set = set(transcripts)
    mapping_set = set(mapping["transcript_to_gene"])
    if quant_set != mapping_set:
        raise GeneIdentifierError(
            "tx2gene must map every quantified transcript exactly once and no others.",
            details={
                "missing_from_tx2gene": sorted(quant_set - mapping_set)[:20],
                "extra_in_tx2gene": sorted(mapping_set - quant_set)[:20],
                "missing_count": len(quant_set - mapping_set),
                "extra_count": len(mapping_set - quant_set),
            },
        )

    replicate_summary, replicate_certifiable = _summarize_replicates(
        replicate_records,
        validation_level=validation_level,
        warnings=warnings,
    )
    certification_eligible = certification_eligible and replicate_certifiable
    if semantics == SALMON_FULL_LENGTH:
        use_replicate_overdispersion = (
            replicate_summary["state"] == "all"
            and replicate_summary["consistent_method_and_count"]
        )
        route = {
            "tximport": {
                "countsFromAbundance": "no",
                "dropInfReps": False,
            },
            "edgeR_constructor": "edgeR::DGEListFromTximport",
            "edgeR_options": {"divide": use_replicate_overdispersion},
            "count_source": "txi$counts",
            "transcript_length_offset": True,
            "inferential_overdispersion": {
                "enabled": use_replicate_overdispersion,
                "source": "tximport.infReps",
                "relative_abundance_adjustment": use_replicate_overdispersion,
            },
        }
    else:
        route = {
            "tximport": {
                "countsFromAbundance": "no",
                "dropInfReps": False,
            },
            "edgeR_constructor": "edgeR::DGEList",
            "count_source": "txi$counts",
            "transcript_length_offset": False,
            "gene_length_correction": False,
            "route_interpretation": "certified_three_prime_default",
            "certified_path_execution_permitted": True,
        }
        if override is not None:
            route["route_interpretation"] = (
                "certified_default_with_blocked_override_request"
            )
            route["certified_path_execution_permitted"] = False
            route["high_risk_override"] = override

    return (
        {
            "salmon": {
                "samples": sample_records,
                "tx2gene_path": tx2gene_fingerprint.path,
                "tx2gene_sha256": tx2gene_fingerprint.sha256,
                "transcript_count": len(transcripts),
                "transcript_id_inventory_sha256": _sequence_digest(transcripts),
                "gene_count": len(mapping["gene_ids"]),
                "gene_id_inventory_sha256": _sequence_digest(mapping["gene_ids"]),
                "estimated_counts_may_be_fractional": True,
                "inferential_replicates": replicate_summary,
                "index_hash_identity": index_hashes,
            },
            "route": route,
        },
        certification_eligible,
    )


def _parse_three_prime_override(
    value: Any,
    *,
    warnings: list[dict[str, Any]],
) -> tuple[dict[str, Any] | None, bool]:
    if value is None:
        return None, True
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'analysis_options' must be an object.")
    require_expected_keys(
        value,
        allowed={"three_prime_length_correction_override"},
        required={"three_prime_length_correction_override"},
        context="analysis_options",
    )
    override = value["three_prime_length_correction_override"]
    if not isinstance(override, Mapping):
        raise InvalidRequestError(
            "'three_prime_length_correction_override' must be an object."
        )
    require_expected_keys(
        override,
        allowed={"enabled", "reason"},
        required={"enabled"},
        context="three_prime_length_correction_override",
    )
    enabled = override["enabled"]
    if not isinstance(enabled, bool):
        raise InvalidRequestError(
            "The three-prime override 'enabled' field must be boolean."
        )
    if not enabled:
        return None, True
    if "reason" not in override:
        raise InputEvidenceRequiredError(
            "A high-risk three-prime length-correction override requires a reason."
        )
    reason = require_nonempty_string(
        override["reason"], field="three_prime_length_correction_override.reason"
    )
    warning = {
        "code": "HIGH_RISK_THREE_PRIME_LENGTH_CORRECTION_OVERRIDE",
        "severity": "high",
        "message": (
            "Gene-length correction was explicitly requested for three-prime RNA-seq; "
            "this run is not eligible for P0 certification."
        ),
        "details": {"reason": reason},
    }
    warnings.append(warning)
    return (
        {
            "gene_length_correction_requested": True,
            "reason": reason,
            "execution_policy": "not_executable_in_certified_path",
        },
        False,
    )


def _required_child(directory: Path, relative: str, *, sample_id: str) -> Path:
    candidate = directory / relative
    try:
        return candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise InputEvidenceRequiredError(
            "A Salmon quantification directory is missing required evidence.",
            details={
                "sample_id": sample_id,
                "quant_dir": str(directory),
                "missing": relative,
            },
            cause=error,
        ) from error


def _resolve_salmon_aux_dir(
    quant_dir: Path,
    *,
    cmd_info: Mapping[str, Any],
    sample_id: str,
) -> tuple[str, Path]:
    """Resolve the exact relative auxDir consumed by locked tximport."""

    raw_aux_dir = cmd_info.get("auxDir", "aux_info")
    if (
        not isinstance(raw_aux_dir, str)
        or not raw_aux_dir
        or raw_aux_dir != raw_aux_dir.strip()
        or has_control_characters(raw_aux_dir)
    ):
        raise InputIntegrityError(
            "Salmon cmd_info.auxDir must be a non-empty exact relative path.",
            details={"sample_id": sample_id},
        )
    if "\\" in raw_aux_dir or re.match(r"^[A-Za-z]:", raw_aux_dir):
        raise InputIntegrityError(
            "Salmon cmd_info.auxDir must use a portable relative POSIX path.",
            details={"sample_id": sample_id, "auxDir": raw_aux_dir},
        )
    pure_path = PurePosixPath(raw_aux_dir)
    parts = raw_aux_dir.split("/")
    if pure_path.is_absolute() or any(part in {"", ".", ".."} for part in parts):
        raise InputIntegrityError(
            "Salmon cmd_info.auxDir must not be absolute or contain unsafe path segments.",
            details={"sample_id": sample_id, "auxDir": raw_aux_dir},
        )
    candidate = quant_dir.joinpath(*parts)
    try:
        resolved = candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise InputEvidenceRequiredError(
            "Salmon cmd_info.auxDir does not identify an available evidence directory.",
            details={
                "sample_id": sample_id,
                "quant_dir": str(quant_dir),
                "auxDir": raw_aux_dir,
            },
            cause=error,
        ) from error
    if not resolved.is_relative_to(quant_dir):
        raise InputIntegrityError(
            "Salmon cmd_info.auxDir resolves outside its quantification directory.",
            details={
                "sample_id": sample_id,
                "quant_dir": str(quant_dir),
                "auxDir": raw_aux_dir,
                "resolved_aux_dir": str(resolved),
            },
        )
    if not resolved.is_dir():
        raise InputEvidenceRequiredError(
            "Salmon cmd_info.auxDir is not an evidence directory.",
            details={
                "sample_id": sample_id,
                "auxDir": raw_aux_dir,
                "resolved_aux_dir": str(resolved),
            },
        )
    return raw_aux_dir, resolved


def _validate_salmon_metadata(
    cmd_info: Mapping[str, Any],
    meta_info: Mapping[str, Any],
    *,
    sample_id: str,
    expected_version: str,
) -> dict[str, Any]:
    versions: list[str] = []
    for role, document in (("cmd_info", cmd_info), ("meta_info", meta_info)):
        if "salmon_version" not in document:
            raise InputEvidenceRequiredError(
                "Every Salmon metadata document must record the producer version.",
                details={"sample_id": sample_id, "document": role},
            )
        versions.append(
            require_nonempty_string(
                document["salmon_version"], field=f"{role}.salmon_version"
            )
        )
    if len(set(versions)) != 1 or versions[0] != expected_version:
        raise InputIntegrityError(
            "Salmon metadata and the requested producer version disagree.",
            details={
                "sample_id": sample_id,
                "metadata_versions": versions,
                "request_version": expected_version,
            },
        )
    index = None
    for key in ("index", "indexDirectory", "index_name"):
        if key in cmd_info and isinstance(cmd_info[key], str) and cmd_info[key].strip():
            index = require_nonempty_string(cmd_info[key], field=f"cmd_info.{key}")
            break
    if index is None:
        raise InputEvidenceRequiredError(
            "Salmon cmd_info.json does not identify the quantification index.",
            details={"sample_id": sample_id},
        )
    quant_errors = meta_info.get("quant_errors", [])
    if quant_errors not in (None, []) and quant_errors != "":
        raise InputIntegrityError(
            "Salmon reports quantification errors for a sample.",
            details={"sample_id": sample_id, "quant_errors": quant_errors},
        )
    index_hashes: dict[str, str] = {}
    for document in (cmd_info, meta_info):
        for key, expected_length in _SALMON_INDEX_HASH_LENGTHS.items():
            if key not in document:
                continue
            raw_value = document[key]
            if (
                not isinstance(raw_value, str)
                or len(raw_value) != expected_length
                or not all(
                    character in "0123456789abcdefABCDEF" for character in raw_value
                )
            ):
                raise InputIntegrityError(
                    "A Salmon index hash is not valid hexadecimal evidence.",
                    details={
                        "sample_id": sample_id,
                        "field": key,
                        "expected_hex_length": expected_length,
                    },
                )
            normalized_value = raw_value.lower()
            previous = index_hashes.get(key)
            if previous is not None and previous != normalized_value:
                raise InputIntegrityError(
                    "Salmon metadata contains conflicting index hashes.",
                    details={"sample_id": sample_id, "field": key},
                )
            index_hashes[key] = normalized_value
    if "index_seq_hash" not in index_hashes:
        raise InputEvidenceRequiredError(
            "Salmon metadata must contain the canonical index_seq_hash identity.",
            details={"sample_id": sample_id},
        )
    return {
        "salmon_version": versions[0],
        "index": index,
        "index_hashes": dict(sorted(index_hashes.items())),
    }


def _read_salmon_quant(
    path: Path,
    *,
    sample_id: str,
    require_offset: bool,
    validation_level: str,
    content: bytes,
) -> list[str]:
    header, rows = _read_tsv(path, role="salmon_quant", content=content)
    expected_header = ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]
    if header != expected_header:
        error_type = (
            SalmonOffsetRequiredError if require_offset else InputIntegrityError
        )
        raise error_type(
            "Salmon quant.sf must contain the canonical length and count columns.",
            details={"sample_id": sample_id, "observed_header": header},
        )
    transcripts: list[str] = []
    seen: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        _require_row_width(row, header, path=path, row_number=row_number)
        transcript = _require_identifier(
            row[0], field=f"quant.sf transcript at row {row_number}"
        )
        if transcript in seen:
            raise GeneIdentifierError(
                "Salmon transcript identifiers must be unique.",
                details={"sample_id": sample_id, "transcript_id": transcript},
            )
        seen.add(transcript)
        transcripts.append(transcript)
        if validation_level == "validate":
            _validate_decimal(
                row[1],
                path=path,
                row_number=row_number,
                field="Length",
                minimum=Decimal(0),
                strictly_greater=True,
            )
            effective_length = _validate_decimal(
                row[2],
                path=path,
                row_number=row_number,
                field="EffectiveLength",
                minimum=Decimal(0),
            )
            _validate_decimal(
                row[3],
                path=path,
                row_number=row_number,
                field="TPM",
                minimum=Decimal(0),
            )
            _validate_decimal(
                row[4],
                path=path,
                row_number=row_number,
                field="NumReads",
                minimum=Decimal(0),
            )
            if require_offset and effective_length <= 0:
                raise SalmonOffsetRequiredError(
                    "Full-length Salmon inputs require positive effective lengths.",
                    details={"sample_id": sample_id, "row": row_number},
                )
    if not transcripts:
        raise InputValidationError(
            "Salmon quant.sf must contain at least one transcript.",
            details={"sample_id": sample_id},
        )
    return transcripts


def _read_tx2gene(
    path: Path,
    *,
    strip_gene_version: bool,
    content: bytes,
) -> dict[str, Any]:
    header, rows = _read_tsv(path, role="tx2gene", content=content)
    if header not in (
        ["transcript_id", "gene_id"],
        ["transcript_id", "gene_id", "gene_name"],
    ):
        raise GeneIdentifierError(
            "tx2gene must have transcript_id and gene_id columns, with optional gene_name.",
            details={"observed_header": header},
        )
    transcript_to_original_gene: dict[str, str] = {}
    symbols_by_original_gene: dict[str, set[str]] = {}
    original_genes: list[str] = []
    for row_number, row in enumerate(rows, start=2):
        _require_row_width(row, header, path=path, row_number=row_number)
        transcript = _require_identifier(
            row[0], field=f"tx2gene transcript at row {row_number}"
        )
        gene = _require_gene_id(row[1], context=f"tx2gene row {row_number}")
        if transcript in transcript_to_original_gene:
            raise GeneIdentifierError(
                "Each tx2gene transcript must appear exactly once.",
                details={"transcript_id": transcript, "row": row_number},
            )
        transcript_to_original_gene[transcript] = gene
        original_genes.append(gene)
        if len(header) == 3 and row[2]:
            symbol = require_nonempty_string(
                row[2], field=f"tx2gene gene_name at row {row_number}"
            )
            symbols_by_original_gene.setdefault(gene, set()).add(symbol)
    if not transcript_to_original_gene:
        raise GeneIdentifierError("tx2gene must contain at least one mapping.")
    conflicts = {
        gene: sorted(symbols)
        for gene, symbols in symbols_by_original_gene.items()
        if len(symbols) > 1
    }
    if conflicts:
        raise GeneIdentifierError(
            "One stable gene_id maps to conflicting display symbols.",
            details={"conflicts": conflicts},
        )
    normalized_genes = _normalize_gene_ids(
        original_genes,
        strip_version=strip_gene_version,
        context="tx2gene",
        allow_repeated_original=True,
    )
    transcript_to_gene = {
        transcript: normalized
        for transcript, normalized in zip(
            transcript_to_original_gene,
            normalized_genes,
            strict=True,
        )
    }
    gene_ids = list(dict.fromkeys(transcript_to_gene.values()))
    return {"transcript_to_gene": transcript_to_gene, "gene_ids": gene_ids}


def _inspect_inferential_replicates(
    aux_dir: Path,
    *,
    aux_declared: str,
    cmd_info: Mapping[str, Any],
    meta_info: Mapping[str, Any],
    transcript_order: Sequence[str],
    sample_id: str,
    quant_declared: str,
    inventory: ArtifactInventory,
) -> dict[str, Any]:
    bootstrap_count = _metadata_nonnegative_integer(
        cmd_info,
        keys=("numBootstraps", "num_bootstraps"),
        sample_id=sample_id,
    )
    gibbs_count = _metadata_nonnegative_integer(
        cmd_info,
        keys=("numGibbsSamples", "num_gibbs_samples", "numGibbs"),
        sample_id=sample_id,
    )
    meta_count = _metadata_nonnegative_integer(
        meta_info,
        keys=("num_bootstraps", "numBootstraps"),
        sample_id=sample_id,
    )
    sample_type = _salmon_sample_type(meta_info, sample_id=sample_id)
    bootstrap_value = bootstrap_count if bootstrap_count is not None else 0
    gibbs_value = gibbs_count if gibbs_count is not None else 0
    if bootstrap_value > 0 and gibbs_value > 0:
        raise InputIntegrityError(
            "Salmon metadata declares both bootstrap and Gibbs replicates.",
            details={"sample_id": sample_id},
        )
    if meta_count is None:
        raise InputEvidenceRequiredError(
            "Salmon meta_info.json must declare num_bootstraps, including an explicit zero.",
            details={"sample_id": sample_id},
        )
    if meta_count == 0:
        if sample_type != "none" or bootstrap_value > 0 or gibbs_value > 0:
            raise InputIntegrityError(
                "Salmon zero-replicate metadata contains a contradictory sampling declaration.",
                details={
                    "sample_id": sample_id,
                    "samp_type": sample_type,
                    "cmd_bootstrap_count": bootstrap_count,
                    "cmd_gibbs_count": gibbs_count,
                    "meta_count": meta_count,
                },
            )
        method = None
        count = 0
    else:
        if sample_type == "none":
            raise InputIntegrityError(
                "Positive Salmon replicate metadata cannot declare samp_type='none'.",
                details={"sample_id": sample_id, "meta_count": meta_count},
            )
        method = sample_type
        relevant_command_count = (
            bootstrap_count if method == "bootstrap" else gibbs_count
        )
        irrelevant_command_count = (
            gibbs_count if method == "bootstrap" else bootstrap_count
        )
        if relevant_command_count is None:
            raise InputEvidenceRequiredError(
                "Positive Salmon replicate metadata requires its matching command-level count.",
                details={
                    "sample_id": sample_id,
                    "method": method,
                    "required_command_fields": (
                        ["numBootstraps", "num_bootstraps"]
                        if method == "bootstrap"
                        else ["numGibbsSamples", "num_gibbs_samples", "numGibbs"]
                    ),
                },
            )
        if (
            irrelevant_command_count not in (None, 0)
            or relevant_command_count != meta_count
        ):
            raise InputIntegrityError(
                "Salmon inferential-replicate counts disagree between metadata files.",
                details={
                    "sample_id": sample_id,
                    "method": method,
                    "cmd_bootstrap_count": bootstrap_count,
                    "cmd_gibbs_count": gibbs_count,
                    "meta_count": meta_count,
                },
            )
        count = meta_count

    bootstrap_path = aux_dir / "bootstrap" / "bootstraps.gz"
    names_path = aux_dir / "bootstrap" / "names.tsv.gz"
    evidence_exists = (bootstrap_path.exists(), names_path.exists())
    if count > 0 and not all(evidence_exists):
        raise InputEvidenceRequiredError(
            "Salmon declares inferential replicates but its archive evidence is incomplete.",
            details={
                "sample_id": sample_id,
                "replicate_count": count,
                "bootstraps_gz_present": evidence_exists[0],
                "names_tsv_gz_present": evidence_exists[1],
            },
        )
    if count == 0 and any(evidence_exists):
        raise InputIntegrityError(
            "Salmon replicate files exist but metadata declares no replicates.",
            details={"sample_id": sample_id},
        )
    file_records: list[dict[str, Any]] = []
    if count > 0:
        target_count = _salmon_target_count(meta_info, sample_id=sample_id)
        if target_count != len(transcript_order):
            raise InputIntegrityError(
                "Salmon target-count metadata does not match quant.sf.",
                details={
                    "sample_id": sample_id,
                    "metadata_target_count": target_count,
                    "quant_transcript_count": len(transcript_order),
                },
            )
        names_fingerprint = inventory.add(
            names_path,
            role=f"salmon.inferential_replicates.{sample_id}.names.tsv.gz",
            declared_path=(f"{quant_declared}/{aux_declared}/bootstrap/names.tsv.gz"),
            capture_for_parse=True,
        )
        _validate_replicate_names(
            names_path=names_path,
            transcript_order=transcript_order,
            sample_id=sample_id,
            names_content=inventory.consume_bytes(Path(names_fingerprint.path)),
        )
        bootstrap_fingerprint = inventory.add(
            bootstrap_path,
            role=f"salmon.inferential_replicates.{sample_id}.bootstraps.gz",
            declared_path=(f"{quant_declared}/{aux_declared}/bootstrap/bootstraps.gz"),
            capture_for_parse=True,
        )
        archive_encoding = _validate_replicate_values(
            bootstrap_path=bootstrap_path,
            replicate_count=count,
            target_count=target_count,
            sample_id=sample_id,
            bootstrap_content=inventory.consume_bytes(Path(bootstrap_fingerprint.path)),
        )
        file_records.extend(
            {
                "path": fingerprint.path,
                "sha256": fingerprint.sha256,
            }
            for fingerprint in (names_fingerprint, bootstrap_fingerprint)
        )
    else:
        archive_encoding = None
    return {
        "sample_id": sample_id,
        "present": count > 0,
        "method": method,
        "count": count,
        "archive_numeric_encoding": archive_encoding,
        "files": file_records,
    }


def _salmon_sample_type(meta_info: Mapping[str, Any], *, sample_id: str) -> str:
    if "samp_type" not in meta_info:
        raise InputEvidenceRequiredError(
            "Salmon meta_info.json must explicitly declare samp_type.",
            details={"sample_id": sample_id},
        )
    raw_sample_type = meta_info["samp_type"]
    if (
        not isinstance(raw_sample_type, str)
        or raw_sample_type != raw_sample_type.strip()
        or has_control_characters(raw_sample_type)
        or raw_sample_type not in {"none", "bootstrap", "gibbs"}
    ):
        raise InputIntegrityError(
            "Salmon samp_type must be exactly 'none', 'bootstrap', or 'gibbs'.",
            details={"sample_id": sample_id, "samp_type": raw_sample_type},
        )
    return raw_sample_type


def _metadata_nonnegative_integer(
    document: Mapping[str, Any],
    *,
    keys: Sequence[str],
    sample_id: str,
) -> int | None:
    values: list[int] = []
    for key in keys:
        if key not in document:
            continue
        raw = document[key]
        if isinstance(raw, bool):
            valid = False
        elif isinstance(raw, int):
            valid = 0 <= raw <= MAX_EXACT_INTEGER_COUNT
            parsed = raw
        elif isinstance(raw, str) and _INTEGER_COUNT_PATTERN.fullmatch(raw):
            valid = not _integer_literal_exceeds(raw, MAX_EXACT_INTEGER_COUNT)
            parsed = int(raw) if valid else 0
        else:
            valid = False
        if not valid:
            raise InputIntegrityError(
                "Salmon replicate metadata must be a bounded nonnegative integer.",
                details={
                    "sample_id": sample_id,
                    "field": key,
                    "value": document[key],
                    "maximum_exact_integer": MAX_EXACT_INTEGER_COUNT,
                },
            )
        values.append(parsed)
    if len(set(values)) > 1:
        raise InputIntegrityError(
            "Salmon replicate metadata fields disagree.",
            details={"sample_id": sample_id, "fields": list(keys), "values": values},
        )
    return values[0] if values else None


def _salmon_target_count(meta_info: Mapping[str, Any], *, sample_id: str) -> int:
    target_key = (
        "num_valid_targets" if "num_valid_targets" in meta_info else "num_targets"
    )
    target_count = _metadata_nonnegative_integer(
        meta_info,
        keys=(target_key,),
        sample_id=sample_id,
    )
    if target_count is None or target_count == 0:
        raise InputEvidenceRequiredError(
            "Salmon replicate metadata must declare a positive target count.",
            details={"sample_id": sample_id},
        )
    return target_count


def _validate_replicate_names(
    *,
    names_path: Path,
    transcript_order: Sequence[str],
    sample_id: str,
    names_content: bytes,
) -> None:
    expected_names_payload = ("\t".join(transcript_order) + "\n").encode("utf-8")
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(names_content), mode="rb") as handle:
            observed_names_payload = handle.read(len(expected_names_payload) + 1)
    except (OSError, EOFError, zlib.error) as error:
        raise InputIntegrityError(
            "Salmon names.tsv.gz is not a valid gzip stream.",
            details={"sample_id": sample_id, "path": str(names_path)},
            cause=error,
        ) from error
    if observed_names_payload != expected_names_payload:
        raise InputIntegrityError(
            "Salmon replicate names must exactly match quant.sf transcript order.",
            details={
                "sample_id": sample_id,
                "expected_name_count": len(transcript_order),
            },
        )


def _validate_replicate_values(
    *,
    bootstrap_path: Path,
    replicate_count: int,
    target_count: int,
    sample_id: str,
    bootstrap_content: bytes,
) -> str:
    expected_values = replicate_count * target_count
    allowed_sizes = {
        expected_values * 8: ("float64_le", "<d", 8),
        expected_values * 4: ("int32_le", "<i", 4),
    }
    maximum_bytes = max(allowed_sizes)
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(bootstrap_content), mode="rb") as handle:
            observed_bytes = 0
            while chunk := handle.read(1024 * 1024):
                observed_bytes += len(chunk)
                if observed_bytes > maximum_bytes:
                    raise InputIntegrityError(
                        "Salmon bootstraps.gz contains more values than metadata declares.",
                        details={
                            "sample_id": sample_id,
                            "accepted_uncompressed_bytes": sorted(allowed_sizes),
                            "observed_uncompressed_bytes_at_failure": observed_bytes,
                        },
                    )
    except InputIntegrityError:
        raise
    except (OSError, EOFError, zlib.error) as error:
        raise InputIntegrityError(
            "Salmon bootstraps.gz is not a valid inferential-replicate stream.",
            details={"sample_id": sample_id, "path": str(bootstrap_path)},
            cause=error,
        ) from error
    if observed_bytes not in allowed_sizes:
        raise InputIntegrityError(
            "Salmon bootstraps.gz dimensions do not match metadata.",
            details={
                "sample_id": sample_id,
                "accepted_uncompressed_bytes": sorted(allowed_sizes),
                "observed_uncompressed_bytes": observed_bytes,
            },
        )
    encoding, unpack_format, value_width = allowed_sizes[observed_bytes]
    remainder = b""
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(bootstrap_content), mode="rb") as handle:
            while chunk := handle.read(1024 * 1024):
                values = remainder + chunk
                complete_bytes = len(values) - (len(values) % value_width)
                for (value,) in struct.iter_unpack(
                    unpack_format, values[:complete_bytes]
                ):
                    if not math.isfinite(value) or value < 0:
                        raise InputIntegrityError(
                            "Salmon inferential-replicate values must be finite and nonnegative.",
                            details={"sample_id": sample_id, "encoding": encoding},
                        )
                remainder = values[complete_bytes:]
    except InputIntegrityError:
        raise
    except (OSError, EOFError, struct.error, zlib.error) as error:
        raise InputIntegrityError(
            "Salmon bootstraps.gz could not be decoded safely.",
            details={"sample_id": sample_id, "encoding": encoding},
            cause=error,
        ) from error
    if remainder:
        raise InputIntegrityError(
            "Salmon bootstraps.gz ends with an incomplete numeric value.",
            details={"sample_id": sample_id, "encoding": encoding},
        )
    return encoding


def _summarize_replicates(
    records: Sequence[Mapping[str, Any]],
    *,
    validation_level: str,
    warnings: list[dict[str, Any]],
) -> tuple[dict[str, Any], bool]:
    present = [bool(record["present"]) for record in records]
    if all(present):
        state = "all"
    elif any(present):
        state = "mixed"
    else:
        state = "none"
    methods = {record["method"] for record in records if record["present"]}
    counts = {record["count"] for record in records if record["present"]}
    consistent = len(methods) <= 1 and len(counts) <= 1
    if state == "mixed" or not consistent:
        details = {
            "state": state,
            "consistent_method_and_count": consistent,
            "per_sample": [dict(record) for record in records],
        }
        if validation_level == "validate":
            raise InputValidationError(
                "Inferential replicates must be present consistently for every sample.",
                details=details,
            )
        warnings.append(
            {
                "code": "INFERENTIAL_REPLICATES_INCONSISTENT",
                "severity": "high",
                "message": "Inferential-replicate evidence is inconsistent across samples.",
                "details": details,
            }
        )
    summary = {
        "state": state,
        "consistent_method_and_count": consistent,
        "method": next(iter(methods)) if len(methods) == 1 else None,
        "replicate_count": next(iter(counts)) if len(counts) == 1 else None,
        "per_sample": [dict(record) for record in records],
    }
    return summary, state != "mixed" and consistent


def _read_tsv(
    path: Path,
    *,
    role: str,
    content: bytes,
) -> tuple[list[str], list[list[str]]]:
    try:
        text = content.decode("utf-8")
        rows = list(
            csv.reader(io.StringIO(text, newline=""), delimiter="\t", strict=True)
        )
    except (UnicodeError, csv.Error) as error:
        raise InputReadError(
            "A tabular input could not be parsed as strict UTF-8 TSV.",
            path=path,
            operation="read_tsv",
            cause=error,
            details={"role": role},
        ) from error
    if not rows or not rows[0]:
        raise InputValidationError(
            "A tabular input must contain a non-empty header.",
            details={"path": str(path), "role": role},
        )
    if any(not row for row in rows):
        raise InputValidationError(
            "Blank rows are not permitted in tabular inputs.",
            details={"path": str(path), "role": role},
        )
    return rows[0], rows[1:]


def _read_featurecounts_tsv(
    path: Path,
    *,
    content: bytes,
) -> tuple[list[str], list[str], list[list[str]]]:
    try:
        raw_lines = content.decode("utf-8").splitlines(keepends=True)
    except UnicodeError as error:
        raise InputReadError(
            "A featureCounts file could not be read as UTF-8.",
            path=path,
            operation="read_featurecounts",
            cause=error,
        ) from error
    comments: list[str] = []
    data_lines: list[str] = []
    for line in raw_lines:
        if not data_lines and line.startswith("#"):
            comments.append(line.rstrip("\r\n"))
        else:
            data_lines.append(line)
    try:
        rows = list(csv.reader(data_lines, delimiter="\t", strict=True))
    except csv.Error as error:
        raise InputReadError(
            "A featureCounts file is not valid TSV.",
            path=path,
            operation="parse_featurecounts",
            cause=error,
        ) from error
    if not rows or not rows[0]:
        raise InputValidationError(
            "A featureCounts file must contain a header.", details={"path": str(path)}
        )
    return comments, rows[0], rows[1:]


def _featurecounts_version(comments: Sequence[str], *, path: Path) -> str:
    versions = {
        match.group(1)
        for comment in comments
        if (match := _FEATURECOUNTS_VERSION_PATTERN.search(comment)) is not None
    }
    if len(versions) != 1:
        raise InputEvidenceRequiredError(
            "A per-sample featureCounts file must contain one unambiguous producer-version header.",
            details={"path": str(path), "observed_versions": sorted(versions)},
        )
    return next(iter(versions))


def _require_unique_header(header: Sequence[str], *, path: Path) -> None:
    duplicates = sorted({column for column in header if header.count(column) > 1})
    control_columns = [
        index for index, column in enumerate(header) if has_control_characters(column)
    ]
    if any(column == "" for column in header) or duplicates or control_columns:
        raise InputValidationError(
            "Tabular headers must contain unique, non-empty, control-free column names.",
            details={
                "path": str(path),
                "duplicate_columns": duplicates,
                "control_character_columns": control_columns,
            },
        )


def _require_row_width(
    row: Sequence[str],
    header: Sequence[str],
    *,
    path: Path,
    row_number: int,
) -> None:
    if len(row) != len(header):
        raise InputValidationError(
            "A tabular row does not match the header width.",
            details={
                "path": str(path),
                "row": row_number,
                "expected_columns": len(header),
                "observed_columns": len(row),
            },
        )


def _validate_integer_count(
    value: str,
    *,
    path: Path,
    row_number: int,
    sample_id: str,
) -> None:
    if not _INTEGER_COUNT_PATTERN.fullmatch(value):
        raise CountValuesInvalidError(
            "featureCounts values must be literal nonnegative integers; no coercion is performed.",
            details={
                "path": str(path),
                "row": row_number,
                "sample_id": sample_id,
                "value": value,
                "maximum_exact_integer": MAX_EXACT_INTEGER_COUNT,
            },
        )
    if _integer_literal_exceeds(value, MAX_EXACT_INTEGER_COUNT):
        raise CountValuesInvalidError(
            "A featureCounts value exceeds the largest integer represented exactly downstream.",
            details={
                "path": str(path),
                "row": row_number,
                "sample_id": sample_id,
                "value": value,
                "maximum_exact_integer": MAX_EXACT_INTEGER_COUNT,
            },
        )


def _integer_literal_exceeds(value: str, maximum: int) -> bool:
    normalized = value.lstrip("0") or "0"
    maximum_text = str(maximum)
    return len(normalized) > len(maximum_text) or (
        len(normalized) == len(maximum_text) and normalized > maximum_text
    )


def _validate_decimal(
    value: str,
    *,
    path: Path,
    row_number: int,
    field: str,
    minimum: Decimal,
    strictly_greater: bool = False,
) -> Decimal:
    if not _DECIMAL_PATTERN.fullmatch(value):
        raise CountValuesInvalidError(
            "A Salmon numeric field does not use the accepted decimal grammar.",
            details={
                "path": str(path),
                "row": row_number,
                "field": field,
                "value": value,
            },
        )
    try:
        number = Decimal(value)
    except InvalidOperation as error:
        raise CountValuesInvalidError(
            "A Salmon numeric field is not a decimal value.",
            details={
                "path": str(path),
                "row": row_number,
                "field": field,
                "value": value,
            },
            cause=error,
        ) from error
    machine_value = float(number) if number.is_finite() else float("nan")
    if (
        not number.is_finite()
        or not math.isfinite(machine_value)
        or (number != 0 and machine_value == 0)
    ):
        invalid_bound = True
    else:
        invalid_bound = number <= minimum if strictly_greater else number < minimum
    if invalid_bound:
        raise CountValuesInvalidError(
            "Salmon numeric fields must be finite and within their valid range.",
            details={
                "path": str(path),
                "row": row_number,
                "field": field,
                "value": value,
            },
        )
    return number


def _require_gene_id(value: Any, *, context: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise GeneIdentifierError(
            "gene_id values must be non-empty strings without surrounding whitespace.",
            details={"context": context, "value": value},
        )
    if any(character.isspace() for character in value):
        raise GeneIdentifierError(
            "gene_id values must not contain whitespace.",
            details={"context": context, "value": value},
        )
    if has_control_characters(value):
        raise GeneIdentifierError(
            "gene_id values must not contain Unicode control characters.",
            details={"context": context},
        )
    return value


def _normalize_gene_ids(
    values: Sequence[str],
    *,
    strip_version: bool,
    context: str,
    allow_repeated_original: bool,
) -> list[str]:
    normalized: list[str] = []
    normalized_to_original: dict[str, str] = {}
    seen_original: set[str] = set()
    for row_index, raw_value in enumerate(values, start=1):
        original = _require_gene_id(raw_value, context=f"{context} item {row_index}")
        if not allow_repeated_original and original in seen_original:
            raise GeneIdentifierError(
                "Gene-level inputs must contain each stable gene_id exactly once.",
                details={"context": context, "duplicate_gene_id": original},
            )
        seen_original.add(original)
        candidate = (
            _GENE_VERSION_PATTERN.sub("", original) if strip_version else original
        )
        previous = normalized_to_original.get(candidate)
        if previous is not None and previous != original:
            raise GeneIdentifierError(
                "Version stripping would collapse distinct stable gene identifiers.",
                details={
                    "context": context,
                    "normalized_gene_id": candidate,
                    "original_gene_ids": sorted({previous, original}),
                },
            )
        normalized_to_original[candidate] = original
        normalized.append(candidate)
    if not normalized:
        raise GeneIdentifierError(
            "A gene-level input must contain at least one gene_id.",
            details={"context": context},
        )
    return normalized


def _sequence_digest(values: Sequence[str]) -> str:
    digest = hashlib.sha256()
    digest.update(b"rnaseq-downstream.sequence.v1\0")
    digest.update(struct.pack(">Q", len(values)))
    for value in values:
        encoded = value.encode("utf-8")
        digest.update(struct.pack(">Q", len(encoded)))
        digest.update(encoded)
    return digest.hexdigest()


__all__ = [
    "FEATURECOUNTS_INTEGER",
    "MAX_EXACT_INTEGER_COUNT",
    "REQUEST_SCHEMA_VERSION",
    "SALMON_FULL_LENGTH",
    "SALMON_THREE_PRIME",
    "SUPPORTED_INPUT_SEMANTICS",
    "inspect_request",
    "validate_request",
]
