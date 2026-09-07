"""Fail-closed reader for published evidence-gated result bundles."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from .errors import InputIntegrityError, InputReadError, InvalidRequestError
from .provenance import has_control_characters

SUMMARY_SCHEMA_VERSION = "1.0"
PATHWAY_SUMMARY_SCHEMA_VERSION = "1.1"
_CORE_FILES = {
    "analysis.json",
    "backend_manifest.json",
    "coefficients.tsv",
    "design.tsv",
    "results.tsv",
}
_MANIFESTED_FILES = _CORE_FILES - {"backend_manifest.json"}
_PATHWAY_RESULT_FILE = "pathway_results.tsv"
_PATHWAY_FILES = _CORE_FILES | {_PATHWAY_RESULT_FILE}
_PATHWAY_MANIFESTED_FILES = _PATHWAY_FILES - {"backend_manifest.json"}
_DESEQ2_FILES = _CORE_FILES
_DESEQ2_MANIFESTED_FILES = _DESEQ2_FILES - {"backend_manifest.json"}
_EXPECTED_RUNTIME = {
    "R": "4.6.1",
    "Bioconductor": "3.23",
    "BiocVersion_package": "3.23.1",
    "edgeR": "4.10.0",
    "tximport": "1.40.0",
    "limma": "3.68.0",
}
_DESEQ2_EXPECTED_RUNTIME = {
    "R": "4.6.1",
    "Bioconductor": "3.23",
    "BiocVersion_package": "3.23.1",
    "DESeq2": "1.52.0",
    "apeglm": "1.34.0",
    "tximport": "1.40.0",
}
_DESEQ2_SCHEMA_VERSION = "1.0"
_DESEQ2_EXECUTION_SCOPE = "validated_p1_deseq2_input"
_RESULT_HEADER = [
    "gene_id",
    "contrast_id",
    "status",
    "logFC",
    "unshrunk_logFC",
    "logCPM",
    "statistic",
    "statistic_type",
    "statistic_status",
    "PValue",
    "FDR",
    "test_method",
    "lfc_threshold",
]
_COEFFICIENT_HEADER = ["gene_id", "status", "coefficient", "estimate", "scale"]
_DESEQ2_RESULT_HEADER = [
    "gene_id",
    "contrast_id",
    "status",
    "status_reason",
    "baseMean",
    "logFC",
    "unshrunk_logFC",
    "lfcSE",
    "statistic",
    "statistic_type",
    "statistic_hypothesis",
    "PValue",
    "FDR",
    "fdr_basis",
    "test_method",
    "lfc_threshold",
    "shrinkage_method",
]
_DESEQ2_COEFFICIENT_HEADER = [
    "gene_id",
    "status",
    "status_reason",
    "coefficient",
    "estimate",
    "scale",
]
_DESIGN_HEADER = ["sample_id", "coefficient", "value"]
_PATHWAY_HEADER = [
    "contrast_id",
    "gene_set_id",
    "gene_set_description",
    "method_id",
    "test_class",
    "hypothesis",
    "inference_role",
    "status",
    "status_reason",
    "direction",
    "proportion_down",
    "proportion_up",
    "p_value",
    "fdr",
    "fdr_family_id",
    "gmt_member_count_raw",
    "gmt_symbol_count_unique",
    "mapped_symbol_count_unique",
    "ambiguous_symbol_count_unique",
    "unmapped_symbol_count_unique",
    "mapping_rate",
    "mapped_gene_id_count_unique",
    "tested_gene_count",
    "filtered_gene_count",
    "tested_universe_gene_count",
    "method_ngenes",
    "correlation_status",
    "correlation_estimate_raw",
    "correlation_effective",
    "vif_used",
    "rotation_status",
]
_STATUSES = {"filtered", "not_tested", "not_estimable", "failed", "tested"}
_PATHWAY_STATUSES = {"tested", "not_tested"}
_PATHWAY_METHODS = {
    "limma_mroast": {
        "test_class": "self_contained",
        "inference_role": "corroborative",
        "hypotheses": ("directional", "mixed"),
    },
    "limma_fry": {
        "test_class": "self_contained",
        "inference_role": "primary",
        "hypotheses": ("directional", "mixed"),
    },
    "limma_camera": {
        "test_class": "competitive",
        "inference_role": "supplementary",
        "hypotheses": ("directional",),
    },
}
_NUMBER_PATTERN = re.compile(
    r"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
)
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _integrity(message: str, **details: Any) -> InputIntegrityError:
    return InputIntegrityError(message, details=details)


def _resolve_run_dir(value: str | Path) -> Path:
    if not isinstance(value, (str, Path)) or not str(value):
        raise InvalidRequestError("The result-bundle directory must be non-empty.")
    if has_control_characters(str(value)):
        raise InvalidRequestError(
            "The result-bundle directory contains a control character."
        )
    candidate = Path(value)
    try:
        if candidate.is_symlink():
            raise OSError("symbolic links are not accepted")
        resolved = candidate.resolve(strict=True)
        if not resolved.is_dir():
            raise NotADirectoryError(str(resolved))
    except (OSError, RuntimeError) as error:
        raise InputReadError(
            "The result bundle could not be opened as a regular directory.",
            path=candidate,
            operation="open_result_bundle",
            cause=error,
        ) from error
    return resolved


def _capture(path: Path) -> tuple[bytes, str, int]:
    try:
        if path.is_symlink():
            raise OSError("symbolic links are not accepted")
        resolved = path.resolve(strict=True)
        if not resolved.is_file():
            raise IsADirectoryError(str(resolved))
        before = resolved.stat()
        with resolved.open("rb") as handle:
            content = handle.read()
        after = resolved.stat()
    except OSError as error:
        raise InputReadError(
            "A result-bundle member could not be captured.",
            path=path,
            operation="capture_result_member",
            cause=error,
        ) from error
    before_identity = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if before_identity != after_identity or len(content) != after.st_size:
        raise _integrity(
            "A result-bundle member changed while it was captured.",
            path=path.name,
        )
    return content, hashlib.sha256(content).hexdigest(), len(content)


def _capture_bundle(
    run_dir: Path,
) -> tuple[dict[str, tuple[bytes, str, int]], Path | None]:
    try:
        before = run_dir.stat()
        entries = list(run_dir.iterdir())
    except OSError as error:
        raise InputReadError(
            "The result bundle could not be enumerated.",
            path=run_dir,
            operation="enumerate_result_bundle",
            cause=error,
        ) from error
    names = {entry.name for entry in entries}
    is_pathway = _PATHWAY_RESULT_FILE in names
    expected_files = _PATHWAY_FILES if is_pathway else _CORE_FILES
    if is_pathway:
        expected_names = expected_files | {"display"}
        inventory_message = (
            "A schema 1.1 pathway bundle must contain exactly the six root "
            "artifacts and its display directory."
        )
    else:
        expected_names = expected_files | ({"display"} if "display" in names else set())
        inventory_message = (
            "The result bundle must contain the five core files and optionally "
            "one display directory."
        )
    unsafe = sorted(
        entry.name
        for entry in entries
        if entry.is_symlink()
        or (entry.name == "display" and not entry.is_dir())
        or (entry.name != "display" and not entry.is_file())
    )
    if names != expected_names or unsafe:
        raise _integrity(
            inventory_message,
            missing_files=sorted(expected_names - names),
            unexpected_files=sorted(names - expected_names),
            unsafe_files=unsafe,
        )
    captured = {name: _capture(run_dir / name) for name in sorted(expected_files)}
    try:
        after = run_dir.stat()
        final_names = {entry.name for entry in run_dir.iterdir()}
    except OSError as error:
        raise InputReadError(
            "The result bundle could not be rechecked after capture.",
            path=run_dir,
            operation="recheck_result_bundle",
            cause=error,
        ) from error
    before_identity = (before.st_dev, before.st_ino)
    after_identity = (after.st_dev, after.st_ino)
    if before_identity != after_identity or final_names != expected_names:
        raise _integrity("The result bundle changed while it was captured.")
    display_dir = run_dir / "display" if "display" in names else None
    return captured, display_dir


def _json_object(content: bytes, *, role: str) -> dict[str, Any]:
    def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError(f"duplicate JSON key: {key}")
            value[key] = item
        return value

    try:
        document = json.loads(
            content.decode("utf-8"),
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-standard numeric constant {value}")
            ),
            object_pairs_hook=reject_duplicate_keys,
        )
    except (UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise _integrity(
            "A result-bundle JSON member is not strict UTF-8 JSON.",
            role=role,
            cause_type=type(error).__name__,
        ) from error
    if not isinstance(document, dict):
        raise _integrity(
            "A result-bundle JSON member must contain an object.", role=role
        )
    return document


def _exact_keys(value: Mapping[str, Any], expected: set[str], *, context: str) -> None:
    observed = set(value)
    if observed != expected:
        raise _integrity(
            f"The {context} schema is incompatible.",
            missing_fields=sorted(expected - observed),
            unexpected_fields=sorted(observed - expected),
        )


def _finite_json_number(
    value: Any, *, context: str, nonnegative: bool = False
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise _integrity(f"The {context} is not a finite JSON number.")
    normalized = float(value)
    if not math.isfinite(normalized) or (nonnegative and normalized < 0):
        raise _integrity(f"The {context} is outside its finite numeric domain.")
    return normalized


def _verify_input_evidence(value: Any) -> str:
    if not isinstance(value, Mapping):
        raise _integrity("The backend manifest lacks validated-input evidence.")
    _exact_keys(
        value,
        {
            "kind",
            "plan_id",
            "bundle_path",
            "validation_bundle_artifacts",
            "r_input_snapshots",
            "digest_coupling",
        },
        context="backend input evidence",
    )
    plan_id = value.get("plan_id")
    if (
        value.get("kind") != "validated_input_bundle"
        or not isinstance(plan_id, str)
        or _SHA256_PATTERN.fullmatch(plan_id) is None
        or not isinstance(value.get("bundle_path"), str)
        or not value["bundle_path"]
        or value.get("digest_coupling") != "private_copy_and_hash_same_source_stream"
    ):
        raise _integrity("The backend input-evidence identity is invalid.")
    inventories = (
        ("validation_bundle_artifacts", "path"),
        ("r_input_snapshots", "source_path"),
    )
    for field, path_field in inventories:
        records = value.get(field)
        if not isinstance(records, list) or not records:
            raise _integrity(
                "A backend input-evidence inventory is empty or invalid.", field=field
            )
        roles: set[str] = set()
        for index, record in enumerate(records):
            if not isinstance(record, Mapping):
                raise _integrity(
                    "A backend input-evidence record is invalid.",
                    field=field,
                    record_index=index,
                )
            role = record.get("role")
            digest = record.get("sha256")
            size = record.get("size_bytes")
            path = record.get(path_field)
            if (
                not isinstance(role, str)
                or not role
                or role in roles
                or not isinstance(path, str)
                or not path
                or not isinstance(digest, str)
                or _SHA256_PATTERN.fullmatch(digest) is None
                or isinstance(size, bool)
                or not isinstance(size, int)
                or size < 0
            ):
                raise _integrity(
                    "A backend input-evidence record has invalid identity fields.",
                    field=field,
                    record_index=index,
                )
            roles.add(role)
    return plan_id


def _verify_manifest(
    manifest: Mapping[str, Any],
    captured: Mapping[str, tuple[bytes, str, int]],
) -> None:
    _exact_keys(
        manifest,
        {
            "schema_version",
            "kind",
            "backend",
            "runtime_identity",
            "execution_scope",
            "input_evidence",
            "members",
        },
        context="backend manifest",
    )
    if manifest.get("execution_scope") == "backend_kernel_only":
        raise InvalidRequestError(
            "Public summarize refuses private benchmark-kernel output.",
            details={"execution_scope": "backend_kernel_only"},
        )
    identity = {
        "schema_version": manifest.get("schema_version"),
        "kind": manifest.get("kind"),
        "backend": manifest.get("backend"),
        "execution_scope": manifest.get("execution_scope"),
        "runtime_identity": manifest.get("runtime_identity"),
    }
    pathway_bundle = _PATHWAY_RESULT_FILE in captured
    schema_version = "1.1" if pathway_bundle else "1.0"
    manifested_files = (
        _PATHWAY_MANIFESTED_FILES if pathway_bundle else _MANIFESTED_FILES
    )
    expected = {
        "schema_version": schema_version,
        "kind": "edger_ql_backend_manifest",
        "backend": "edgeR_QL",
        "execution_scope": "validated_p0_input",
        "runtime_identity": _EXPECTED_RUNTIME,
    }
    if identity != expected:
        raise _integrity(
            "The backend manifest identity is not the locked public P0 identity.",
            observed=identity,
        )
    _verify_input_evidence(manifest.get("input_evidence"))
    members = manifest.get("members")
    if not isinstance(members, list):
        raise _integrity("The backend manifest member inventory is not an array.")
    observed: set[str] = set()
    for index, member in enumerate(members):
        if not isinstance(member, Mapping):
            raise _integrity(
                "A backend manifest member is not an object.", member_index=index
            )
        _exact_keys(
            member,
            {"path", "sha256", "size_bytes"},
            context=f"backend manifest member {index}",
        )
        name = member.get("path")
        digest = member.get("sha256")
        size = member.get("size_bytes")
        if name not in manifested_files or name in observed:
            raise _integrity(
                "A backend manifest member path is invalid.",
                member_index=index,
                path=name,
            )
        if not isinstance(digest, str) or _SHA256_PATTERN.fullmatch(digest) is None:
            raise _integrity("A backend manifest member digest is invalid.", path=name)
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise _integrity("A backend manifest member size is invalid.", path=name)
        observed.add(name)
        _, actual_digest, actual_size = captured[name]
        if digest != actual_digest or size != actual_size:
            raise _integrity(
                "A result member does not match the backend manifest.", path=name
            )
    if observed != manifested_files:
        raise _integrity(
            "The backend manifest omits required result members.",
            missing_members=sorted(manifested_files - observed),
        )


def _verify_deseq2_manifest(
    manifest: Mapping[str, Any],
    captured: Mapping[str, tuple[bytes, str, int]],
) -> None:
    """Verify the independent DESeq2 bundle identity and member commitments."""

    _exact_keys(
        manifest,
        {
            "schema_version",
            "kind",
            "backend",
            "runtime_identity",
            "execution_scope",
            "input_evidence",
            "members",
        },
        context="DESeq2 backend manifest",
    )
    if manifest.get("execution_scope") == "backend_kernel_only":
        raise InvalidRequestError(
            "Public summarize refuses private benchmark-kernel output.",
            details={"execution_scope": "backend_kernel_only"},
        )
    identity = {
        "schema_version": manifest.get("schema_version"),
        "kind": manifest.get("kind"),
        "backend": manifest.get("backend"),
        "execution_scope": manifest.get("execution_scope"),
        "runtime_identity": manifest.get("runtime_identity"),
    }
    expected = {
        "schema_version": _DESEQ2_SCHEMA_VERSION,
        "kind": "deseq2_backend_manifest",
        "backend": "DESeq2",
        "execution_scope": _DESEQ2_EXECUTION_SCOPE,
        "runtime_identity": _DESEQ2_EXPECTED_RUNTIME,
    }
    if identity != expected:
        raise _integrity(
            "The DESeq2 backend manifest identity is incompatible.",
            observed=identity,
        )
    _verify_input_evidence(manifest.get("input_evidence"))
    members = manifest.get("members")
    if not isinstance(members, list):
        raise _integrity("The DESeq2 manifest member inventory is not an array.")
    observed: set[str] = set()
    for index, member in enumerate(members):
        if not isinstance(member, Mapping):
            raise _integrity(
                "A DESeq2 manifest member is not an object.", member_index=index
            )
        _exact_keys(
            member,
            {"path", "sha256", "size_bytes"},
            context=f"DESeq2 backend manifest member {index}",
        )
        name = member.get("path")
        digest = member.get("sha256")
        size = member.get("size_bytes")
        if name not in _DESEQ2_MANIFESTED_FILES or name in observed:
            raise _integrity(
                "A DESeq2 manifest member path is invalid.",
                member_index=index,
                path=name,
            )
        if not isinstance(digest, str) or _SHA256_PATTERN.fullmatch(digest) is None:
            raise _integrity("A DESeq2 manifest member digest is invalid.", path=name)
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise _integrity("A DESeq2 manifest member size is invalid.", path=name)
        observed.add(name)
        _, actual_digest, actual_size = captured[name]
        if digest != actual_digest or size != actual_size:
            raise _integrity(
                "A DESeq2 result member does not match its manifest.", path=name
            )
    if observed != _DESEQ2_MANIFESTED_FILES:
        raise _integrity(
            "The DESeq2 manifest omits required result members.",
            missing_members=sorted(_DESEQ2_MANIFESTED_FILES - observed),
        )


def _tsv_rows(
    content: bytes, *, expected_header: Sequence[str], role: str
) -> list[list[str]]:
    try:
        text = content.decode("utf-8")
        rows = list(
            csv.reader(
                io.StringIO(text, newline=""),
                delimiter="\t",
                quoting=csv.QUOTE_NONE,
                strict=True,
            )
        )
    except (UnicodeError, csv.Error) as error:
        raise _integrity(
            "A result member is not strict UTF-8 TSV.",
            role=role,
            cause_type=type(error).__name__,
        ) from error
    if not rows or rows[0] != list(expected_header):
        raise _integrity(
            "A result TSV has an incompatible header.",
            role=role,
            observed_header=rows[0] if rows else None,
        )
    result = rows[1:]
    for row_number, row in enumerate(result, start=2):
        if len(row) != len(expected_header):
            raise _integrity(
                "A result TSV row has an incompatible width.",
                role=role,
                row=row_number,
            )
    return result


def _number(value: str, *, role: str, row: int, field: str) -> float:
    if _NUMBER_PATTERN.fullmatch(value) is None:
        raise _integrity(
            "A result TSV numeric field is not a strict finite decimal.",
            role=role,
            row=row,
            field=field,
            observed=value,
        )
    parsed = float(value)
    if not math.isfinite(parsed):
        raise _integrity(
            "A result TSV numeric field is non-finite.",
            role=role,
            row=row,
            field=field,
        )
    return parsed


def _benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    """Recompute the within-contrast BH adjustment recorded by the backend."""

    count = len(p_values)
    adjusted = [0.0] * count
    running_minimum = 1.0
    order = sorted(range(count), key=lambda index: p_values[index])
    for reverse_rank in range(count - 1, -1, -1):
        index = order[reverse_rank]
        rank = reverse_rank + 1
        candidate = min(1.0, p_values[index] * count / rank)
        running_minimum = min(running_minimum, candidate)
        adjusted[index] = running_minimum
    return adjusted


def _analysis_contrasts(analysis: Mapping[str, Any]) -> list[dict[str, Any]]:
    value = analysis.get("contrasts")
    if not isinstance(value, list) or not value:
        raise _integrity("The analysis contrast inventory is invalid.")
    normalized: list[dict[str, Any]] = []
    identifiers: set[str] = set()
    for index, item in enumerate(value):
        if not isinstance(item, Mapping):
            raise _integrity(
                "An analysis contrast record is invalid.", contrast_index=index
            )
        required = {
            "contrast_id",
            "weights",
            "lfc_threshold",
            "test_method",
            "treat_null",
            "estimability_residual",
            "estimability_tolerance",
        }
        _exact_keys(item, required, context=f"analysis contrast {index}")
        identifier = item.get("contrast_id")
        threshold = item.get("lfc_threshold")
        method = item.get("test_method")
        if (
            not isinstance(identifier, str)
            or not identifier
            or identifier in identifiers
        ):
            raise _integrity("Analysis contrast identifiers are invalid or duplicated.")
        if isinstance(threshold, bool) or not isinstance(threshold, (int, float)):
            raise _integrity(
                "An analysis contrast threshold is invalid.", contrast_id=identifier
            )
        threshold = float(threshold)
        if not math.isfinite(threshold) or threshold < 0:
            raise _integrity(
                "An analysis contrast threshold is invalid.", contrast_id=identifier
            )
        expected_method = "glmTreat" if threshold > 0 else "glmQLFTest"
        expected_null = "interval" if threshold > 0 else None
        if method != expected_method or item.get("treat_null") != expected_null:
            raise _integrity(
                "An analysis contrast does not match the fixed threshold-test policy.",
                contrast_id=identifier,
            )
        weights = item.get("weights")
        if not isinstance(weights, Mapping) or not weights:
            raise _integrity(
                "An analysis contrast weight inventory is invalid.",
                contrast_id=identifier,
            )
        normalized_weights = []
        for coefficient, weight in weights.items():
            if not isinstance(coefficient, str) or not coefficient:
                raise _integrity(
                    "An analysis contrast coefficient name is invalid.",
                    contrast_id=identifier,
                )
            normalized_weights.append(
                _finite_json_number(
                    weight,
                    context=(f"analysis contrast {identifier} weight {coefficient}"),
                )
            )
        if not any(normalized_weights):
            raise _integrity(
                "An analysis contrast has identically zero weights.",
                contrast_id=identifier,
            )
        residual = _finite_json_number(
            item.get("estimability_residual"),
            context=f"analysis contrast {identifier} estimability residual",
            nonnegative=True,
        )
        tolerance = _finite_json_number(
            item.get("estimability_tolerance"),
            context=f"analysis contrast {identifier} estimability tolerance",
            nonnegative=True,
        )
        if tolerance <= 0:
            raise _integrity(
                "An analysis contrast estimability tolerance is not positive.",
                contrast_id=identifier,
            )
        if residual > tolerance:
            raise _integrity(
                "An analysis contrast exceeds its estimability tolerance.",
                contrast_id=identifier,
            )
        identifiers.add(identifier)
        normalized.append(
            {
                "contrast_id": identifier,
                "lfc_threshold": threshold,
                "test_method": method,
            }
        )
    return normalized


def _verify_analysis(
    analysis: Mapping[str, Any], manifest: Mapping[str, Any]
) -> list[dict[str, Any]]:
    expected_keys = {
        "schema_version",
        "kind",
        "backend",
        "execution_scope",
        "analysis_request",
        "input_evidence",
        "runtime_identity",
        "input_semantics",
        "route_observed",
        "pipeline",
        "design",
        "contrasts",
        "genes",
        "status_vocabulary",
        "result_logFC_scale",
        "coefficient_scale",
        "multiple_testing",
    }
    manifest_schema = manifest.get("schema_version")
    if manifest_schema == "1.1":
        expected_keys.add("pathway_analysis")
    _exact_keys(analysis, expected_keys, context="analysis document")
    identity = {
        "schema_version": analysis.get("schema_version"),
        "kind": analysis.get("kind"),
        "backend": analysis.get("backend"),
        "execution_scope": analysis.get("execution_scope"),
        "runtime_identity": analysis.get("runtime_identity"),
    }
    expected = {
        "schema_version": manifest_schema,
        "kind": "edger_ql_analysis",
        "backend": "edgeR_QL",
        "execution_scope": "validated_p0_input",
        "runtime_identity": _EXPECTED_RUNTIME,
    }
    if identity != expected:
        raise _integrity("The analysis identity is incompatible.", observed=identity)
    if analysis.get("input_evidence") != manifest.get("input_evidence"):
        raise _integrity("Analysis and manifest input evidence disagree.")
    request = analysis.get("analysis_request")
    if not isinstance(request, Mapping):
        raise _integrity("The analysis request identity is invalid.")
    _exact_keys(request, {"path", "sha256"}, context="analysis request identity")
    if (
        not isinstance(request.get("path"), str)
        or not request["path"]
        or not isinstance(request.get("sha256"), str)
        or _SHA256_PATTERN.fullmatch(request["sha256"]) is None
    ):
        raise _integrity("The analysis request identity is invalid.")
    if analysis.get("status_vocabulary") != [
        "filtered",
        "not_tested",
        "not_estimable",
        "failed",
        "tested",
    ]:
        raise _integrity("The analysis status vocabulary is incompatible.")
    if (
        analysis.get("result_logFC_scale") != "log2"
        or analysis.get("coefficient_scale") != "natural_log"
        or analysis.get("multiple_testing") != "Benjamini-Hochberg within each contrast"
    ):
        raise _integrity("The analysis scale or multiple-testing identity is invalid.")
    semantics = analysis.get("input_semantics")
    if semantics not in {
        "featurecounts_integer",
        "salmon_quant_dirs_full_length",
        "salmon_quant_dirs_three_prime",
    }:
        raise _integrity("The analysis input semantics are not a public P0 route.")
    expected_pipeline = [
        {"step": "filterByExpr", "arguments": {"design": True}},
        {"step": "normLibSizes", "arguments": {"method": "TMM"}},
        {
            "step": "glmQLFit",
            "arguments": {"legacy": False, "robust": True},
        },
        {
            "step": "contrast_test",
            "dispatch": ("lfc_threshold == 0: glmQLFTest; lfc_threshold > 0: glmTreat"),
        },
    ]
    if analysis.get("pipeline") != expected_pipeline:
        raise _integrity("The analysis pipeline identity is incompatible.")
    route = analysis.get("route_observed")
    if not isinstance(route, Mapping):
        raise _integrity("The observed input route is invalid.")
    if semantics == "featurecounts_integer":
        if route != {
            "constructor": "edgeR::DGEList",
            "count_semantics": "integer",
            "transcript_length_offset": False,
        }:
            raise _integrity("The observed featureCounts route is incompatible.")
    elif semantics == "salmon_quant_dirs_full_length":
        fixed_route = {
            "constructor": "edgeR::DGEListFromTximport",
            "transcript_length_offset": True,
            "countsFromAbundance": "no",
            "dropInfReps": False,
        }
        if any(route.get(key) != value for key, value in fixed_route.items()) or (
            set(route)
            != {
                *fixed_route,
                "divide",
                "inferential_replicates_imported",
            }
        ):
            raise _integrity("The observed full-length Salmon route is incompatible.")
        if not isinstance(route.get("divide"), bool) or not isinstance(
            route.get("inferential_replicates_imported"), bool
        ):
            raise _integrity("The observed full-length Salmon route is invalid.")
        if route["divide"] != route["inferential_replicates_imported"]:
            raise _integrity(
                "The full-length divide policy disagrees with imported "
                "inferential-replicate evidence."
            )
    else:
        fixed_route = {
            "constructor": "edgeR::DGEList",
            "count_source": "txi$counts",
            "transcript_length_offset": False,
            "gene_length_correction": False,
            "countsFromAbundance": "no",
            "dropInfReps": False,
            "divide": False,
        }
        if any(route.get(key) != value for key, value in fixed_route.items()) or (
            set(route) != {*fixed_route, "inferential_replicates_imported"}
        ):
            raise _integrity("The observed three-prime Salmon route is incompatible.")
        if not isinstance(route.get("inferential_replicates_imported"), bool):
            raise _integrity("The observed three-prime Salmon route is invalid.")
    return _analysis_contrasts(analysis)


def _verify_design(
    rows: Sequence[Sequence[str]], analysis: Mapping[str, Any]
) -> tuple[list[str], list[str]]:
    design = analysis.get("design")
    if not isinstance(design, Mapping):
        raise _integrity("The analysis design evidence is invalid.")
    _exact_keys(
        design,
        {
            "intercept",
            "terms",
            "variable_types",
            "factor_levels",
            "columns",
            "sample_count",
            "rank",
            "residual_df",
            "qr_tolerance",
        },
        context="analysis design evidence",
    )
    raw_columns = design.get("columns")
    columns = [raw_columns] if isinstance(raw_columns, str) else raw_columns
    if (
        not isinstance(columns, list)
        or not columns
        or not all(isinstance(value, str) and value for value in columns)
    ):
        raise _integrity("The analysis design-column inventory is invalid.")
    if len(set(columns)) != len(columns):
        raise _integrity("The analysis design columns are duplicated.")
    sample_count = design.get("sample_count")
    rank = design.get("rank")
    residual_df = design.get("residual_df")
    if any(
        isinstance(value, bool) or not isinstance(value, int)
        for value in (sample_count, rank, residual_df)
    ):
        raise _integrity("The analysis design dimensions are invalid.")
    if (
        sample_count <= 0
        or rank <= 0
        or residual_df <= 0
        or rank != len(columns)
        or residual_df != sample_count - rank
    ):
        raise _integrity("The analysis design dimensions are inconsistent.")
    if not isinstance(design.get("intercept"), bool):
        raise _integrity("The analysis design intercept flag is invalid.")
    qr_tolerance = _finite_json_number(
        design.get("qr_tolerance"),
        context="analysis design QR tolerance",
        nonnegative=True,
    )
    if qr_tolerance <= 0:
        raise _integrity("The analysis design QR tolerance is not positive.")
    observed: set[tuple[str, str]] = set()
    samples: list[str] = []
    sample_set: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        sample, coefficient, raw = row
        if not sample or coefficient not in columns:
            raise _integrity("A design TSV identity field is invalid.", row=row_number)
        key = (sample, coefficient)
        if key in observed:
            raise _integrity(
                "The design TSV contains a duplicate cell.", row=row_number
            )
        observed.add(key)
        _number(raw, role="design.tsv", row=row_number, field="value")
        if sample not in sample_set:
            sample_set.add(sample)
            samples.append(sample)
    expected = {(sample, column) for sample in samples for column in columns}
    if len(samples) != sample_count or observed != expected:
        raise _integrity("The design TSV is not the declared complete matrix.")
    return samples, list(columns)


def _verify_coefficients(
    rows: Sequence[Sequence[str]], *, columns: Sequence[str], expected_genes: set[str]
) -> dict[str, str]:
    observed: set[tuple[str, str]] = set()
    gene_statuses: dict[str, str] = {}
    for row_number, row in enumerate(rows, start=2):
        gene, status, coefficient, estimate, scale = row
        if (
            not gene
            or coefficient not in columns
            or status not in {"filtered", "tested"}
        ):
            raise _integrity(
                "A coefficient TSV identity or status is invalid.", row=row_number
            )
        key = (gene, coefficient)
        if key in observed:
            raise _integrity(
                "The coefficient TSV contains a duplicate gene/coefficient row.",
                row=row_number,
            )
        observed.add(key)
        previous = gene_statuses.setdefault(gene, status)
        if previous != status or scale != "natural_log":
            raise _integrity(
                "A coefficient row disagrees with its gene status or scale.",
                row=row_number,
            )
        if status == "tested":
            _number(estimate, role="coefficients.tsv", row=row_number, field="estimate")
        elif estimate:
            raise _integrity(
                "A filtered coefficient estimate must be empty.", row=row_number
            )
    expected = {
        (gene, coefficient) for gene in expected_genes for coefficient in columns
    }
    if set(gene_statuses) != expected_genes or observed != expected:
        raise _integrity("The coefficient TSV is not the complete declared matrix.")
    return gene_statuses


def _verify_results(
    rows: Sequence[Sequence[str]],
    *,
    contrasts: Sequence[Mapping[str, Any]],
    analysis: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], set[str], dict[str, str]]:
    contrast_by_id = {item["contrast_id"]: item for item in contrasts}
    observed: set[tuple[str, str]] = set()
    genes_by_contrast: dict[str, set[str]] = {
        identifier: set() for identifier in contrast_by_id
    }
    gene_filter_status: dict[str, str] = {}
    status_counts = {
        identifier: {status: 0 for status in sorted(_STATUSES)}
        for identifier in contrast_by_id
    }
    significance = {
        identifier: {"fdr_le_0_05": 0, "fdr_gt_0_05": 0, "not_tested": 0}
        for identifier in contrast_by_id
    }
    tested_probabilities: dict[str, list[tuple[float, float, int]]] = {
        identifier: [] for identifier in contrast_by_id
    }
    for row_number, row in enumerate(rows, start=2):
        values = dict(zip(_RESULT_HEADER, row, strict=True))
        gene = values["gene_id"]
        identifier = values["contrast_id"]
        status = values["status"]
        if not gene or identifier not in contrast_by_id or status not in _STATUSES:
            raise _integrity(
                "A result row identity or status is invalid.", row=row_number
            )
        if status not in {"filtered", "tested"}:
            raise _integrity(
                "A verified-complete backend result may contain only filtered or tested rows.",
                row=row_number,
                status=status,
            )
        key = (gene, identifier)
        if key in observed:
            raise _integrity(
                "The results TSV contains a duplicate gene/contrast row.",
                row=row_number,
            )
        observed.add(key)
        genes_by_contrast[identifier].add(gene)
        specification = contrast_by_id[identifier]
        threshold = _number(
            values["lfc_threshold"],
            role="results.tsv",
            row=row_number,
            field="lfc_threshold",
        )
        if (
            threshold != specification["lfc_threshold"]
            or values["test_method"] != (specification["test_method"])
        ):
            raise _integrity(
                "A result row disagrees with its declared contrast test.",
                row=row_number,
            )
        is_treat = values["test_method"] == "glmTreat"
        expected_statistic_type = "not_reported_by_glmTreat" if is_treat else "F"
        if values["statistic_type"] != expected_statistic_type:
            raise _integrity(
                "A result statistic type disagrees with its method.", row=row_number
            )
        numeric_outcomes = [
            "logFC",
            "unshrunk_logFC",
            "logCPM",
            "statistic",
            "PValue",
            "FDR",
        ]
        if status == "tested":
            required = ["logFC", "logCPM", "PValue", "FDR"]
            parsed = {
                field: _number(
                    values[field], role="results.tsv", row=row_number, field=field
                )
                for field in required
            }
            for probability in ("PValue", "FDR"):
                if not 0 <= parsed[probability] <= 1:
                    raise _integrity(
                        "A tested probability is outside [0, 1].",
                        row=row_number,
                        field=probability,
                    )
            if values["unshrunk_logFC"]:
                _number(
                    values["unshrunk_logFC"],
                    role="results.tsv",
                    row=row_number,
                    field="unshrunk_logFC",
                )
            if is_treat:
                if values["statistic"] or values["statistic_status"] != "not_reported":
                    raise _integrity(
                        "glmTreat must explicitly carry an unreported statistic.",
                        row=row_number,
                    )
            else:
                statistic = _number(
                    values["statistic"],
                    role="results.tsv",
                    row=row_number,
                    field="statistic",
                )
                if statistic < 0 or values["statistic_status"] != "reported":
                    raise _integrity(
                        "glmQLFTest must carry a reported non-negative F statistic.",
                        row=row_number,
                    )
            significance[identifier][
                "fdr_le_0_05" if parsed["FDR"] <= 0.05 else "fdr_gt_0_05"
            ] += 1
            tested_probabilities[identifier].append(
                (parsed["PValue"], parsed["FDR"], row_number)
            )
        else:
            if any(values[field] for field in numeric_outcomes):
                raise _integrity(
                    "A non-tested result row must not carry numeric outcomes.",
                    row=row_number,
                    status=status,
                )
            expected_status = (
                "not_applicable_filtered"
                if status == "filtered"
                else f"not_applicable_{status}"
            )
            if values["statistic_status"] != expected_status:
                raise _integrity(
                    "A non-tested result has an incompatible statistic status.",
                    row=row_number,
                )
            significance[identifier]["not_tested"] += 1
        status_counts[identifier][status] += 1
        filter_class = "tested" if status == "tested" else "filtered"
        previous = gene_filter_status.setdefault(gene, filter_class)
        if previous != filter_class:
            raise _integrity(
                "A gene changes filter/test status across contrasts.", gene_id=gene
            )
    if not rows:
        raise _integrity("The results TSV contains no result rows.")
    inventories = list(genes_by_contrast.values())
    genes = inventories[0]
    if any(inventory != genes for inventory in inventories[1:]):
        raise _integrity("Contrasts do not contain the same complete gene inventory.")
    expected_pairs = {
        (gene, identifier) for gene in genes for identifier in contrast_by_id
    }
    if observed != expected_pairs:
        raise _integrity("The results TSV is not the complete gene/contrast matrix.")
    for identifier, probabilities in tested_probabilities.items():
        expected_fdr = _benjamini_hochberg([p_value for p_value, _, _ in probabilities])
        for (_, observed_fdr, row_number), expected_value in zip(
            probabilities, expected_fdr, strict=True
        ):
            if not math.isclose(
                observed_fdr,
                expected_value,
                rel_tol=1e-10,
                abs_tol=1e-12,
            ):
                raise _integrity(
                    "A tested FDR does not equal the declared within-contrast "
                    "Benjamini-Hochberg adjustment.",
                    contrast_id=identifier,
                    row=row_number,
                    observed_fdr=observed_fdr,
                    expected_fdr=expected_value,
                )
    gene_evidence = analysis.get("genes")
    if not isinstance(gene_evidence, Mapping):
        raise _integrity("The analysis gene-count evidence is invalid.")
    if set(gene_evidence) != {"total", "tested", "filtered"} or any(
        isinstance(value, bool) or not isinstance(value, int) or value < 0
        for value in gene_evidence.values()
    ):
        raise _integrity("The analysis gene-count fields are invalid.")
    expected_gene_counts = {
        "total": len(genes),
        "tested": sum(status == "tested" for status in gene_filter_status.values()),
        "filtered": sum(status == "filtered" for status in gene_filter_status.values()),
    }
    if gene_evidence != expected_gene_counts:
        raise _integrity(
            "The analysis gene counts disagree with results.tsv.",
            expected=expected_gene_counts,
            observed=gene_evidence,
        )
    summaries = []
    for specification in contrasts:
        identifier = specification["contrast_id"]
        summaries.append(
            {
                **dict(specification),
                "status_counts": status_counts[identifier],
                "significance_counts": significance[identifier],
            }
        )
    return summaries, genes, gene_filter_status


def _verify_deseq2_defaults(
    value: Any, *, design: Mapping[str, Any]
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise _integrity("The DESeq2 defaults evidence is invalid.")
    expected_keys = {
        "fit_type_requested",
        "fit_type_resolved",
        "size_factor_type",
        "beta_prior",
        "min_replicates_for_replace",
        "independent_filtering",
        "cooks_cutoff",
        "alpha",
        "p_adjust_method",
        "results_alt_hypothesis",
        "use_t",
        "minmu",
        "parallel",
        "outlier_replacement_count",
    }
    _exact_keys(value, expected_keys, context="DESeq2 defaults evidence")
    expected_fixed = {
        "fit_type_requested": "parametric",
        "size_factor_type": "ratio",
        "beta_prior": False,
        "min_replicates_for_replace": 7,
        "independent_filtering": True,
        "alpha": 0.1,
        "p_adjust_method": "BH",
        "results_alt_hypothesis": "greaterAbs",
        "use_t": False,
        "minmu": 0.5,
        "parallel": False,
    }
    if any(value.get(key) != expected for key, expected in expected_fixed.items()):
        raise _integrity("The DESeq2 defaults evidence is incompatible.")
    if value.get("fit_type_resolved") not in ("parametric", "local"):
        raise _integrity("The resolved DESeq2 fit type is invalid.")
    replacement_count = value.get("outlier_replacement_count")
    if (
        isinstance(replacement_count, bool)
        or not isinstance(replacement_count, int)
        or replacement_count < 0
    ):
        raise _integrity("The DESeq2 outlier-replacement count is invalid.")
    cooks = value.get("cooks_cutoff")
    if not isinstance(cooks, Mapping):
        raise _integrity("The DESeq2 Cook's cutoff evidence is invalid.")
    _exact_keys(
        cooks,
        {
            "requested",
            "resolved_f_quantile",
            "resolved_value",
            "numerator_df",
            "denominator_df",
        },
        context="DESeq2 Cook's cutoff evidence",
    )
    numerator = cooks.get("numerator_df")
    denominator = cooks.get("denominator_df")
    if (
        cooks.get("requested") != "automatic"
        or _finite_json_number(
            cooks.get("resolved_f_quantile"),
            context="DESeq2 Cook's F quantile",
            nonnegative=True,
        )
        != 0.99
        or _finite_json_number(
            cooks.get("resolved_value"),
            context="DESeq2 resolved Cook's cutoff",
            nonnegative=True,
        )
        <= 0
        or isinstance(numerator, bool)
        or not isinstance(numerator, int)
        or numerator != design.get("rank")
        or isinstance(denominator, bool)
        or not isinstance(denominator, int)
        or denominator != design.get("residual_df")
    ):
        raise _integrity("The DESeq2 Cook's cutoff evidence is inconsistent.")
    return dict(value)


def _verify_deseq2_test(value: Any, *, design: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise _integrity("The DESeq2 test evidence is invalid.")
    mode = value.get("mode")
    shrinkage = value.get("shrinkage")
    if shrinkage not in ("none", "apeglm"):
        raise _integrity("The DESeq2 shrinkage policy is invalid.")
    if mode == "wald":
        _exact_keys(
            value, {"mode", "shrinkage", "reduced"}, context="DESeq2 Wald test"
        )
        if value.get("reduced") is not None:
            raise _integrity("A DESeq2 Wald test cannot carry a reduced design.")
        return {"mode": mode, "shrinkage": shrinkage, "reduced": None}
    if mode != "lrt":
        raise _integrity("The DESeq2 test mode is invalid.")
    _exact_keys(
        value,
        {
            "mode",
            "shrinkage",
            "reduced",
        },
        context="DESeq2 LRT test",
    )
    reduced = value.get("reduced")
    if not isinstance(reduced, Mapping):
        raise _integrity("The DESeq2 LRT reduced-design evidence is invalid.")
    _exact_keys(
        reduced,
        {
            "intercept",
            "terms",
            "columns",
            "rank",
            "residual_df",
            "nesting_residual",
            "nesting_tolerance",
        },
        context="DESeq2 LRT reduced design",
    )
    if reduced.get("intercept") is not design.get("intercept"):
        raise _integrity("The DESeq2 reduced design changes the intercept policy.")
    raw_terms = reduced.get("terms")
    terms = [raw_terms] if isinstance(raw_terms, str) else raw_terms
    raw_full_terms = design.get("terms")
    full_terms = [raw_full_terms] if isinstance(raw_full_terms, str) else raw_full_terms
    if (
        not isinstance(terms, list)
        or not all(isinstance(term, str) and term for term in terms)
        or len(set(terms)) != len(terms)
        or not isinstance(full_terms, list)
        or not set(terms) < set(full_terms)
    ):
        raise _integrity("The DESeq2 reduced-design terms are not a strict subset.")
    raw_columns = reduced.get("columns")
    columns = [raw_columns] if isinstance(raw_columns, str) else raw_columns
    full_columns = design.get("columns")
    if isinstance(full_columns, str):
        full_columns = [full_columns]
    if (
        not isinstance(columns, list)
        or not columns
        or not all(isinstance(column, str) and column for column in columns)
        or len(set(columns)) != len(columns)
        or not isinstance(full_columns, list)
        or not set(columns).issubset(set(full_columns))
    ):
        raise _integrity("The DESeq2 reduced-design columns are invalid.")
    reduced_rank = reduced.get("rank")
    full_rank = design.get("rank")
    reduced_residual_df = reduced.get("residual_df")
    sample_count = design.get("sample_count")
    if (
        isinstance(reduced_rank, bool)
        or not isinstance(reduced_rank, int)
        or isinstance(full_rank, bool)
        or not isinstance(full_rank, int)
        or reduced_rank != len(columns)
        or reduced_rank <= 0
        or reduced_rank >= full_rank
        or isinstance(reduced_residual_df, bool)
        or not isinstance(reduced_residual_df, int)
        or isinstance(sample_count, bool)
        or not isinstance(sample_count, int)
        or reduced_residual_df != sample_count - reduced_rank
    ):
        raise _integrity("The DESeq2 LRT ranks or degrees of freedom are inconsistent.")
    residual = _finite_json_number(
        reduced.get("nesting_residual"),
        context="DESeq2 reduced-design nesting residual",
        nonnegative=True,
    )
    tolerance = _finite_json_number(
        reduced.get("nesting_tolerance"),
        context="DESeq2 reduced-design nesting tolerance",
        nonnegative=True,
    )
    if tolerance <= 0 or residual > tolerance:
        raise _integrity("The DESeq2 reduced design is not demonstrably nested.")
    return {
        "mode": mode,
        "shrinkage": shrinkage,
        "reduced": {
            **dict(reduced),
            "terms": list(terms),
            "nesting_residual": residual,
            "nesting_tolerance": tolerance,
        },
        "degrees_of_freedom": full_rank - reduced_rank,
    }


def _verify_deseq2_contrasts(
    value: Any,
    *,
    design_columns: Sequence[str],
    test: Mapping[str, Any],
    defaults: Mapping[str, Any],
) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not value:
        raise _integrity("The DESeq2 contrast inventory is invalid.")
    if test["mode"] == "lrt" and len(value) != 1:
        raise _integrity("A DESeq2 LRT bundle must contain exactly one reporting contrast.")
    expected_keys = {
        "contrast_id",
        "weights",
        "lfc_threshold",
        "test_method",
        "alternative_hypothesis",
        "independent_filter_threshold",
        "independent_filter_theta",
        "independent_filter_alpha",
        "cooks_filter_applied",
        "cooks_cutoff",
        "coefficient_name",
        "shrinkage_method",
        "shrinkage_nonconvergence_count",
        "estimability_residual",
        "estimability_tolerance",
    }
    identifiers: set[str] = set()
    normalized: list[dict[str, Any]] = []
    for index, item in enumerate(value):
        if not isinstance(item, Mapping):
            raise _integrity("A DESeq2 contrast record is invalid.", contrast_index=index)
        _exact_keys(item, expected_keys, context=f"DESeq2 contrast {index}")
        identifier = item.get("contrast_id")
        if (
            not isinstance(identifier, str)
            or not identifier
            or has_control_characters(identifier)
            or identifier in identifiers
        ):
            raise _integrity("DESeq2 contrast identifiers are invalid or duplicated.")
        weights = item.get("weights")
        if not isinstance(weights, Mapping) or not weights:
            raise _integrity("A DESeq2 contrast has no coefficient weights.", contrast_id=identifier)
        parsed_weights: dict[str, float] = {}
        for coefficient, weight in weights.items():
            if not isinstance(coefficient, str) or coefficient not in design_columns:
                raise _integrity(
                    "A DESeq2 contrast names an unknown design coefficient.",
                    contrast_id=identifier,
                    coefficient=coefficient,
                )
            parsed_weights[coefficient] = _finite_json_number(
                weight, context=f"DESeq2 contrast {identifier} weight {coefficient}"
            )
        if not any(parsed_weights.values()):
            raise _integrity("A DESeq2 contrast has identically zero weights.", contrast_id=identifier)
        threshold = _finite_json_number(
            item.get("lfc_threshold"),
            context=f"DESeq2 contrast {identifier} LFC threshold",
            nonnegative=True,
        )
        if test["mode"] == "lrt":
            expected_method = "DESeq2::results_LRT"
            expected_alternative = "full_vs_reduced_omnibus"
            if threshold != 0:
                raise _integrity("A DESeq2 LRT cannot carry an LFC threshold.")
        elif threshold > 0:
            expected_method = "DESeq2::results_Wald_greaterAbs"
            expected_alternative = "greaterAbs"
        else:
            expected_method = "DESeq2::results_Wald"
            expected_alternative = "greaterAbs_at_zero_equivalent_two_sided"
        if (
            item.get("test_method") != expected_method
            or item.get("alternative_hypothesis") != expected_alternative
        ):
            raise _integrity("A DESeq2 contrast test identity is incompatible.", contrast_id=identifier)
        filter_threshold = _finite_json_number(
            item.get("independent_filter_threshold"),
            context=f"DESeq2 contrast {identifier} filter threshold",
            nonnegative=True,
        )
        filter_theta = _finite_json_number(
            item.get("independent_filter_theta"),
            context=f"DESeq2 contrast {identifier} filter theta",
            nonnegative=True,
        )
        filter_alpha = _finite_json_number(
            item.get("independent_filter_alpha"),
            context=f"DESeq2 contrast {identifier} independent-filter alpha",
            nonnegative=True,
        )
        cooks_cutoff = _finite_json_number(
            item.get("cooks_cutoff"),
            context=f"DESeq2 contrast {identifier} Cook's cutoff",
            nonnegative=True,
        )
        shrinkage = item.get("shrinkage_method")
        coefficient_name = item.get("coefficient_name")
        if shrinkage not in ("none", "apeglm"):
            raise _integrity("A DESeq2 contrast has an invalid shrinkage method.")
        if shrinkage == "none" and coefficient_name is not None:
            raise _integrity("An unshrunk DESeq2 contrast cannot name an apeglm coefficient.")
        if shrinkage == "apeglm" and (
            not isinstance(coefficient_name, str) or not coefficient_name
        ):
            raise _integrity("An apeglm contrast must name one DESeq2 coefficient.")
        shrinkage_nonconvergence_count = item.get(
            "shrinkage_nonconvergence_count"
        )
        if (
            isinstance(shrinkage_nonconvergence_count, bool)
            or not isinstance(shrinkage_nonconvergence_count, int)
            or shrinkage_nonconvergence_count < 0
            or (shrinkage == "none" and shrinkage_nonconvergence_count != 0)
        ):
            raise _integrity(
                "A DESeq2 shrinkage convergence count is invalid.",
                contrast_id=identifier,
            )
        residual = _finite_json_number(
            item.get("estimability_residual"),
            context=f"DESeq2 contrast {identifier} estimability residual",
            nonnegative=True,
        )
        tolerance = _finite_json_number(
            item.get("estimability_tolerance"),
            context=f"DESeq2 contrast {identifier} estimability tolerance",
            nonnegative=True,
        )
        if (
            filter_theta > 1
            or filter_alpha != 0.1
            or item.get("cooks_filter_applied") is not True
            or not math.isclose(
                cooks_cutoff,
                defaults["cooks_cutoff"]["resolved_value"],
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
            or tolerance <= 0
            or residual > tolerance
        ):
            raise _integrity("A DESeq2 contrast carries invalid diagnostic thresholds.")
        identifiers.add(identifier)
        normalized.append(
            {
                "contrast_id": identifier,
                "weights": parsed_weights,
                "lfc_threshold": threshold,
                "test_method": expected_method,
                "alternative_hypothesis": expected_alternative,
                "independent_filter_threshold": filter_threshold,
                "independent_filter_theta": filter_theta,
                "independent_filter_alpha": filter_alpha,
                "cooks_filter_applied": True,
                "cooks_cutoff": cooks_cutoff,
                "coefficient_name": coefficient_name,
                "shrinkage_method": shrinkage,
                "shrinkage_nonconvergence_count": (
                    shrinkage_nonconvergence_count
                ),
                "estimability_residual": residual,
                "estimability_tolerance": tolerance,
            }
        )
    if len({item["shrinkage_method"] for item in normalized}) != 1:
        raise _integrity("DESeq2 shrinkage policy changes across reporting contrasts.")
    if normalized[0]["shrinkage_method"] != test["shrinkage"]:
        raise _integrity("DESeq2 contrasts disagree with the declared shrinkage policy.")
    return normalized


def _verify_deseq2_reporting_effects(
    value: Any,
    *,
    contrasts: Sequence[Mapping[str, Any]],
    coefficient_mapping: Mapping[str, str],
    test: Mapping[str, Any],
) -> None:
    if not isinstance(value, list) or len(value) != len(contrasts):
        raise _integrity("The DESeq2 reporting-effect inventory is invalid.")
    by_id = {item["contrast_id"]: item for item in contrasts}
    observed: set[str] = set()
    for index, item in enumerate(value):
        if not isinstance(item, Mapping):
            raise _integrity("A DESeq2 reporting-effect record is invalid.", effect_index=index)
        _exact_keys(
            item,
            {"contrast_id", "role", "weights", "coefficient_name"},
            context=f"DESeq2 reporting effect {index}",
        )
        identifier = item.get("contrast_id")
        if (
            not isinstance(identifier, str)
            or identifier not in by_id
            or identifier in observed
        ):
            raise _integrity("A DESeq2 reporting effect has an invalid contrast identity.")
        expected_role = (
            "reported_effect_not_lrt_hypothesis"
            if test["mode"] == "lrt"
            else "tested_contrast"
        )
        if item.get("role") != expected_role or item.get("weights") != by_id[identifier]["weights"]:
            raise _integrity("A DESeq2 reporting effect does not match its contrast.")
        coefficient_name = item.get("coefficient_name")
        contrast = by_id[identifier]
        if coefficient_name != contrast["coefficient_name"]:
            raise _integrity("A DESeq2 reporting effect changes its coefficient name.")
        if contrast["shrinkage_method"] == "apeglm":
            weights = contrast["weights"]
            if (
                len(weights) != 1
                or next(iter(weights.values())) != 1
                or coefficient_name != coefficient_mapping.get(next(iter(weights)))
            ):
                raise _integrity("An apeglm reporting effect is not a single +1 coefficient.")
        observed.add(identifier)
    if observed != set(by_id):
        raise _integrity("The DESeq2 reporting-effect inventory is incomplete.")


def _verify_deseq2_design(
    rows: Sequence[Sequence[str]], analysis: Mapping[str, Any]
) -> tuple[list[str], list[str], dict[str, str]]:
    design = analysis.get("design")
    if not isinstance(design, Mapping):
        raise _integrity("The DESeq2 design evidence is invalid.")
    _exact_keys(
        design,
        {
            "intercept",
            "terms",
            "variable_types",
            "factor_levels",
            "columns",
            "coefficient_mapping",
            "sample_count",
            "rank",
            "residual_df",
            "qr_tolerance",
        },
        context="DESeq2 design evidence",
    )
    edge_compatible_analysis = dict(analysis)
    edge_compatible_design = dict(design)
    mapping_records = edge_compatible_design.pop("coefficient_mapping")
    edge_compatible_analysis["design"] = edge_compatible_design
    sample_ids, columns = _verify_design(rows, edge_compatible_analysis)
    if not isinstance(mapping_records, list) or len(mapping_records) != len(columns):
        raise _integrity("The DESeq2 coefficient-name mapping is incomplete.")
    mapping: dict[str, str] = {}
    result_names: set[str] = set()
    for index, record in enumerate(mapping_records, start=1):
        if not isinstance(record, Mapping):
            raise _integrity("A DESeq2 coefficient-name mapping record is invalid.")
        _exact_keys(
            record,
            {"design_coefficient", "deseq2_result_name", "position"},
            context=f"DESeq2 coefficient-name mapping {index}",
        )
        design_name = record.get("design_coefficient")
        result_name = record.get("deseq2_result_name")
        position = record.get("position")
        if (
            design_name != columns[index - 1]
            or not isinstance(result_name, str)
            or not result_name
            or result_name in result_names
            or isinstance(position, bool)
            or not isinstance(position, int)
            or position != index
        ):
            raise _integrity("The DESeq2 coefficient-name mapping is inconsistent.")
        mapping[design_name] = result_name
        result_names.add(result_name)
    return sample_ids, columns, mapping


def _verify_deseq2_analysis(
    analysis: Mapping[str, Any], manifest: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    _exact_keys(
        analysis,
        {
            "schema_version",
            "kind",
            "backend",
            "execution_scope",
            "analysis_request",
            "input_evidence",
            "runtime_identity",
            "input_semantics",
            "route_observed",
            "pipeline",
            "defaults",
            "design",
            "test",
            "contrasts",
            "reporting_effect",
            "genes",
            "status_vocabulary",
            "result_logFC_scale",
            "coefficient_scale",
            "multiple_testing",
        },
        context="DESeq2 analysis document",
    )
    identity = {
        "schema_version": analysis.get("schema_version"),
        "kind": analysis.get("kind"),
        "backend": analysis.get("backend"),
        "execution_scope": analysis.get("execution_scope"),
        "runtime_identity": analysis.get("runtime_identity"),
    }
    expected_identity = {
        "schema_version": _DESEQ2_SCHEMA_VERSION,
        "kind": "deseq2_analysis",
        "backend": "DESeq2",
        "execution_scope": _DESEQ2_EXECUTION_SCOPE,
        "runtime_identity": _DESEQ2_EXPECTED_RUNTIME,
    }
    if identity != expected_identity:
        raise _integrity("The DESeq2 analysis identity is incompatible.", observed=identity)
    if analysis.get("input_evidence") != manifest.get("input_evidence"):
        raise _integrity("DESeq2 analysis and manifest input evidence disagree.")
    request = analysis.get("analysis_request")
    if not isinstance(request, Mapping):
        raise _integrity("The DESeq2 analysis-request identity is invalid.")
    _exact_keys(
        request,
        {"path", "sha256", "schema_version", "backend"},
        context="DESeq2 analysis-request identity",
    )
    if (
        not isinstance(request.get("path"), str)
        or not request["path"]
        or not isinstance(request.get("sha256"), str)
        or _SHA256_PATTERN.fullmatch(request["sha256"]) is None
        or request.get("schema_version") != "1.2"
        or request.get("backend") != "deseq2"
    ):
        raise _integrity("The DESeq2 analysis-request identity is invalid.")
    if analysis.get("status_vocabulary") != [
        "filtered",
        "not_tested",
        "not_estimable",
        "failed",
        "tested",
    ]:
        raise _integrity("The DESeq2 status vocabulary is incompatible.")
    if (
        analysis.get("result_logFC_scale") != "log2"
        or analysis.get("coefficient_scale") != "log2"
        or analysis.get("multiple_testing")
        != (
            "Benjamini-Hochberg within each reporting contrast after "
            "DESeq2 independent filtering"
        )
    ):
        raise _integrity("The DESeq2 scale or multiple-testing identity is invalid.")
    semantics = analysis.get("input_semantics")
    if semantics not in (
        "featurecounts_integer",
        "salmon_quant_dirs_full_length",
        "salmon_quant_dirs_three_prime",
    ):
        raise _integrity("The DESeq2 input semantics are unsupported.")
    design = analysis.get("design")
    if not isinstance(design, Mapping):
        raise _integrity("The DESeq2 design evidence is invalid.")
    defaults = _verify_deseq2_defaults(analysis.get("defaults"), design=design)
    raw_columns = design.get("columns")
    design_columns = [raw_columns] if isinstance(raw_columns, str) else raw_columns
    if not isinstance(design_columns, list):
        raise _integrity("The DESeq2 design-column inventory is invalid.")
    test = _verify_deseq2_test(analysis.get("test"), design=design)
    contrasts = _verify_deseq2_contrasts(
        analysis.get("contrasts"),
        design_columns=design_columns,
        test=test,
        defaults=defaults,
    )
    expected_pipeline = [
        {
            "step": "construct_DESeqDataSet",
            "constructor": (
                analysis.get("route_observed", {}).get("constructor")
                if isinstance(analysis.get("route_observed"), Mapping)
                else None
            ),
        },
        {
            "step": "DESeq",
            "arguments": {
                "test": "Wald" if test["mode"] == "wald" else "LRT",
                "fitType": "parametric",
                "sfType": "ratio",
                "betaPrior": False,
                "minReplicatesForReplace": 7,
                "useT": False,
                "minmu": 0.5,
                "parallel": False,
            },
        },
        {
            "step": "results",
            "arguments": {
                "independentFiltering": True,
                "cooksCutoff": "automatic",
                "alpha": 0.1,
                "pAdjustMethod": "BH",
                "altHypothesis": "greaterAbs",
                "parallel": False,
            },
        },
        {
            "step": "lfcShrink",
            "method": test["shrinkage"],
            "arguments": (
                {"lfcThreshold": 0, "svalue": False}
                if test["shrinkage"] == "apeglm"
                else None
            ),
            "inferential_columns_changed": False,
        },
    ]
    if analysis.get("pipeline") != expected_pipeline:
        raise _integrity("The DESeq2 pipeline evidence is incompatible.")
    genes = analysis.get("genes")
    if not isinstance(genes, Mapping):
        raise _integrity("The DESeq2 gene evidence is invalid.")
    _exact_keys(
        genes,
        {"total", "result_rows", "status_counts"},
        context="DESeq2 gene evidence",
    )
    total = genes.get("total")
    result_rows = genes.get("result_rows")
    if (
        isinstance(total, bool)
        or not isinstance(total, int)
        or total <= 0
        or isinstance(result_rows, bool)
        or not isinstance(result_rows, int)
        or result_rows != total * len(contrasts)
    ):
        raise _integrity("The DESeq2 gene dimensions are invalid.")
    status_counts = genes.get("status_counts")
    if (
        not isinstance(status_counts, Mapping)
        or set(status_counts) != _STATUSES
        or any(
            isinstance(count, bool) or not isinstance(count, int) or count < 0
            for count in status_counts.values()
        )
        or sum(status_counts.values()) != result_rows
    ):
        raise _integrity("The aggregate DESeq2 status-count inventory is invalid.")
    return contrasts, {"test": test, "defaults": defaults}


def _verify_deseq2_rounding_audit(
    value: Any,
    *,
    expected_mode: str,
    sample_ids: Sequence[str],
    gene_count: int,
) -> None:
    if not isinstance(value, Mapping):
        raise _integrity("The DESeq2 count-rounding audit is invalid.")
    _exact_keys(
        value,
        {
            "mode",
            "round_function",
            "rule",
            "hash_serialization",
            "source_matrix_sha256",
            "rounded_matrix_sha256",
            "gene_count",
            "sample_count",
            "cell_count",
            "changed_cell_count",
            "max_absolute_delta",
            "absolute_delta_sum",
            "per_sample",
        },
        context="DESeq2 count-rounding audit",
    )
    if value.get("mode") != expected_mode:
        raise _integrity("The DESeq2 count-rounding mode is incompatible.")
    expected_rule = (
        "not_applicable_exact_integer_input"
        if expected_mode == "none_integer_input"
        else "R_base_round_IEC_60559_ties_to_even"
    )
    expected_function = "none" if expected_mode == "none_integer_input" else "base::round"
    if (
        value.get("rule") != expected_rule
        or value.get("round_function") != expected_function
        or value.get("hash_serialization")
        != "R_serialize_version_3_locked_runtime"
    ):
        raise _integrity("The DESeq2 count-rounding rule is incompatible.")
    for field in ("source_matrix_sha256", "rounded_matrix_sha256"):
        digest = value.get(field)
        if not isinstance(digest, str) or _SHA256_PATTERN.fullmatch(digest) is None:
            raise _integrity("A DESeq2 count-rounding matrix digest is invalid.", field=field)
    sample_count = value.get("sample_count")
    observed_gene_count = value.get("gene_count")
    cell_count = value.get("cell_count")
    changed_count = value.get("changed_cell_count")
    if any(
        isinstance(item, bool) or not isinstance(item, int) or item < 0
        for item in (sample_count, observed_gene_count, cell_count, changed_count)
    ) or (
        sample_count != len(sample_ids)
        or observed_gene_count != gene_count
        or cell_count != sample_count * gene_count
        or changed_count > cell_count
    ):
        raise _integrity("The DESeq2 count-rounding dimensions are inconsistent.")
    maximum = _finite_json_number(
        value.get("max_absolute_delta"),
        context="DESeq2 maximum rounding delta",
        nonnegative=True,
    )
    absolute_sum = _finite_json_number(
        value.get("absolute_delta_sum"),
        context="DESeq2 absolute rounding-delta sum",
        nonnegative=True,
    )
    if maximum > 0.5 + 1e-12:
        raise _integrity("The DESeq2 count-rounding delta exceeds one half.")
    per_sample = value.get("per_sample")
    if not isinstance(per_sample, list) or len(per_sample) != len(sample_ids):
        raise _integrity("The DESeq2 per-sample rounding audit is incomplete.")
    observed_samples: list[str] = []
    summed_changed = 0
    summed_absolute = 0.0
    sample_maxima: list[float] = []
    for index, record in enumerate(per_sample):
        if not isinstance(record, Mapping):
            raise _integrity("A DESeq2 per-sample rounding record is invalid.")
        _exact_keys(
            record,
            {
                "sample_id",
                "cell_count",
                "changed_cell_count",
                "max_absolute_delta",
                "absolute_delta_sum",
                "total_count_before",
                "total_count_after",
                "total_count_delta",
            },
            context=f"DESeq2 per-sample rounding record {index}",
        )
        sample_id = record.get("sample_id")
        record_cells = record.get("cell_count")
        record_changed = record.get("changed_cell_count")
        if (
            not isinstance(sample_id, str)
            or not sample_id
            or isinstance(record_cells, bool)
            or not isinstance(record_cells, int)
            or record_cells != gene_count
            or isinstance(record_changed, bool)
            or not isinstance(record_changed, int)
            or not 0 <= record_changed <= gene_count
        ):
            raise _integrity("A DESeq2 per-sample rounding record is inconsistent.")
        record_maximum = _finite_json_number(
            record.get("max_absolute_delta"),
            context=f"DESeq2 sample {sample_id} maximum rounding delta",
            nonnegative=True,
        )
        record_absolute = _finite_json_number(
            record.get("absolute_delta_sum"),
            context=f"DESeq2 sample {sample_id} absolute rounding-delta sum",
            nonnegative=True,
        )
        before = _finite_json_number(
            record.get("total_count_before"),
            context=f"DESeq2 sample {sample_id} total before rounding",
            nonnegative=True,
        )
        after = _finite_json_number(
            record.get("total_count_after"),
            context=f"DESeq2 sample {sample_id} total after rounding",
            nonnegative=True,
        )
        delta = _finite_json_number(
            record.get("total_count_delta"),
            context=f"DESeq2 sample {sample_id} signed total-count delta",
        )
        if (
            record_maximum > 0.5 + 1e-12
            or record_absolute + 1e-12 < abs(delta)
            or not math.isclose(after - before, delta, rel_tol=1e-12, abs_tol=1e-9)
            or not float(after).is_integer()
        ):
            raise _integrity("A DESeq2 per-sample rounding delta is inconsistent.")
        observed_samples.append(sample_id)
        summed_changed += record_changed
        summed_absolute += record_absolute
        sample_maxima.append(record_maximum)
    if observed_samples != list(sample_ids):
        raise _integrity("The DESeq2 rounding audit sample order is incompatible.")
    expected_maximum = max(sample_maxima, default=0.0)
    if (
        summed_changed != changed_count
        or not math.isclose(summed_absolute, absolute_sum, rel_tol=1e-12, abs_tol=1e-9)
        or not math.isclose(expected_maximum, maximum, rel_tol=1e-12, abs_tol=1e-12)
    ):
        raise _integrity("The DESeq2 aggregate rounding audit is inconsistent.")
    if expected_mode == "none_integer_input" and (
        value["source_matrix_sha256"] != value["rounded_matrix_sha256"]
        or changed_count != 0
        or maximum != 0
        or absolute_sum != 0
        or any(
            record["total_count_delta"] != 0
            for record in per_sample
            if isinstance(record, Mapping)
        )
    ):
        raise _integrity("Integer featureCounts input cannot carry a rounding mutation.")


def _verify_deseq2_route(
    analysis: Mapping[str, Any],
    *,
    sample_ids: Sequence[str],
    gene_count: int,
) -> None:
    route = analysis.get("route_observed")
    if not isinstance(route, Mapping):
        raise _integrity("The observed DESeq2 input route is invalid.")
    semantics = analysis.get("input_semantics")
    if semantics == "featurecounts_integer":
        _exact_keys(
            route,
            {
                "constructor",
                "count_source",
                "count_semantics",
                "transcript_length_offset",
                "gene_length_correction",
                "rounding_audit",
                "inferential_replicates_imported",
                "inferential_replicates_used_for_inference",
                "inferential_replicates_unused_reason",
            },
            context="DESeq2 featureCounts route",
        )
        fixed = {
            "constructor": "DESeq2::DESeqDataSetFromMatrix",
            "count_source": "validated_featureCounts_integer_counts",
            "count_semantics": "integer",
            "transcript_length_offset": False,
            "gene_length_correction": False,
            "inferential_replicates_imported": False,
            "inferential_replicates_used_for_inference": False,
            "inferential_replicates_unused_reason": "not_applicable",
        }
        if any(route.get(key) != expected for key, expected in fixed.items()):
            raise _integrity("The observed DESeq2 featureCounts route is incompatible.")
        expected_rounding_mode = "none_integer_input"
    elif semantics == "salmon_quant_dirs_full_length":
        fixed = {
            "constructor": "DESeq2::DESeqDataSetFromTximport",
            "count_source": "txi$counts",
            "count_semantics": "salmon_estimated_counts_rounded_for_DESeq2",
            "transcript_length_offset": True,
            "gene_length_correction": True,
            "countsFromAbundance": "no",
            "dropInfReps": False,
            "inferential_replicates_used_for_inference": False,
        }
        _exact_keys(
            route,
            {
                *fixed,
                "inferential_replicates_imported",
                "inferential_replicate_state",
                "inferential_replicates_per_sample",
                "inferential_replicate_method",
                "inferential_replicates_unused_reason",
                "rounding_audit",
                "rounding_disclosure",
            },
            context="DESeq2 full-length Salmon route",
        )
        if any(route.get(key) != expected for key, expected in fixed.items()):
            raise _integrity("The observed DESeq2 full-length Salmon route is incompatible.")
        expected_rounding_mode = "deseq2_constructor"
    else:
        fixed = {
            "constructor": "DESeq2::DESeqDataSetFromMatrix",
            "count_source": "txi$counts",
            "count_semantics": "salmon_estimated_counts_rounded_for_DESeq2",
            "transcript_length_offset": False,
            "gene_length_correction": False,
            "countsFromAbundance": "no",
            "dropInfReps": False,
            "inferential_replicates_used_for_inference": False,
        }
        _exact_keys(
            route,
            {
                *fixed,
                "inferential_replicates_imported",
                "inferential_replicate_state",
                "inferential_replicates_per_sample",
                "inferential_replicate_method",
                "inferential_replicates_unused_reason",
                "rounding_audit",
                "rounding_disclosure",
            },
            context="DESeq2 three-prime Salmon route",
        )
        if any(route.get(key) != expected for key, expected in fixed.items()):
            raise _integrity("The observed DESeq2 three-prime Salmon route is incompatible.")
        expected_rounding_mode = "explicit_before_matrix_constructor"

    if semantics != "featurecounts_integer":
        count = route.get("inferential_replicates_per_sample")
        state = route.get("inferential_replicate_state")
        method = route.get("inferential_replicate_method")
        imported = route.get("inferential_replicates_imported")
        if (
            isinstance(count, bool)
            or not isinstance(count, int)
            or count < 0
            or state not in ("none", "all")
            or not isinstance(imported, bool)
        ):
            raise _integrity("The DESeq2 inferential-replicate counts are invalid.")
        if (
            (state == "none" and (count != 0 or method is not None or imported))
            or (
                state == "all"
                and (
                    count < 1
                    or not isinstance(method, str)
                    or not method
                    or not imported
                )
            )
        ):
            raise _integrity("The DESeq2 inferential-replicate import flag is invalid.")
        if count == 1:
            raise _integrity(
                "Salmon input must have zero or at least two inferential "
                "replicates per sample."
            )
        expected_unused = (
            "DESeq2 1.52.0 does not consume tximport infReps in this backend; "
            "they are imported only to verify declared evidence."
        )
        expected_disclosure = (
            "DESeq2::DESeqDataSetFromTximport internally calls base::round; "
            "the pre/post conversion is audited here."
            if semantics == "salmon_quant_dirs_full_length"
            else (
                "The toolkit explicitly calls base::round before "
                "DESeqDataSetFromMatrix; no length offset is attached."
            )
        )
        if (
            route.get("inferential_replicates_unused_reason") != expected_unused
            or route.get("rounding_disclosure") != expected_disclosure
        ):
            raise _integrity("The DESeq2 Salmon disclosure is incompatible.")
    _verify_deseq2_rounding_audit(
        route.get("rounding_audit"),
        expected_mode=expected_rounding_mode,
        sample_ids=sample_ids,
        gene_count=gene_count,
    )


_DESEQ2_MODEL_FAILURE_REASONS = {
    "wald_model_nonconvergence",
    "lrt_model_nonconvergence",
    "lrt_full_model_nonconvergence",
    "lrt_reduced_model_nonconvergence",
    "lrt_full_and_reduced_model_nonconvergence",
}
_DESEQ2_RESULT_FAILURE_REASONS = _DESEQ2_MODEL_FAILURE_REASONS | {
    "cooks_outlier",
    "apeglm_nonconvergence",
}


def _verify_deseq2_coefficients(
    rows: Sequence[Sequence[str]],
    *,
    columns: Sequence[str],
    expected_genes: set[str],
) -> tuple[
    dict[str, dict[str, float]],
    dict[str, set[tuple[str, str]]],
]:
    observed: set[tuple[str, str]] = set()
    observed_genes: set[str] = set()
    estimates: dict[str, dict[str, float]] = {}
    gene_states: dict[str, set[tuple[str, str]]] = {}
    for row_number, row in enumerate(rows, start=2):
        gene, status, reason, coefficient, estimate, scale = row
        if not gene or coefficient not in columns or status not in _STATUSES:
            raise _integrity(
                "A DESeq2 coefficient identity or status is invalid.", row=row_number
            )
        if status not in {"tested", "not_tested", "failed"}:
            raise _integrity(
                "A completed DESeq2 coefficient row has an impossible status.",
                row=row_number,
                status=status,
            )
        key = (gene, coefficient)
        if key in observed:
            raise _integrity(
                "The DESeq2 coefficient TSV contains a duplicate gene/coefficient row.",
                row=row_number,
            )
        observed.add(key)
        observed_genes.add(gene)
        gene_states.setdefault(gene, set()).add((status, reason))
        if scale != "log2":
            raise _integrity(
                "A DESeq2 coefficient row has an incompatible scale.",
                row=row_number,
            )
        if status == "tested":
            if reason != "fitted":
                raise _integrity("A fitted DESeq2 coefficient has an invalid reason.")
            value = _number(
                estimate, role="coefficients.tsv", row=row_number, field="estimate"
            )
        else:
            expected_reason = "all_zero" if status == "not_tested" else None
            if (
                (expected_reason is not None and reason != expected_reason)
                or (
                    status == "failed"
                    and reason
                    not in (_DESEQ2_MODEL_FAILURE_REASONS | {"coefficient_unavailable"})
                )
                or (status == "not_tested" and estimate)
            ):
                raise _integrity(
                    "An unavailable DESeq2 coefficient has incompatible evidence.",
                    row=row_number,
                )
            value = (
                _number(
                    estimate,
                    role="coefficients.tsv",
                    row=row_number,
                    field="estimate",
                )
                if estimate
                else None
            )
            if reason == "coefficient_unavailable" and value is not None:
                raise _integrity("An unavailable DESeq2 coefficient cannot carry a value.")
        if value is not None:
            estimates.setdefault(gene, {})[coefficient] = value
    expected = {
        (gene, coefficient) for gene in expected_genes for coefficient in columns
    }
    if observed_genes != expected_genes or observed != expected:
        raise _integrity("The DESeq2 coefficient TSV is not the complete declared matrix.")
    for gene, states in gene_states.items():
        if any(status == "not_tested" for status, _ in states):
            if states != {("not_tested", "all_zero")}:
                raise _integrity(
                    "An all-zero DESeq2 gene has inconsistent coefficient states.",
                    gene_id=gene,
                )
            continue
        model_failures = {
            reason
            for status, reason in states
            if status == "failed" and reason in _DESEQ2_MODEL_FAILURE_REASONS
        }
        if model_failures and (
            len(model_failures) != 1
            or states != {("failed", next(iter(model_failures)))}
        ):
            raise _integrity(
                "A model-failed DESeq2 gene has inconsistent coefficient states.",
                gene_id=gene,
            )
    return estimates, gene_states


def _optional_number(
    value: str, *, role: str, row: int, field: str
) -> float | None:
    if not value:
        return None
    return _number(value, role=role, row=row, field=field)


def _deseq2_expected_result_identity(
    contrast: Mapping[str, Any], *, test: Mapping[str, Any]
) -> tuple[str, str, str, str]:
    if test["mode"] == "lrt":
        return (
            "DESeq2::results_LRT",
            "LRT",
            "full_vs_reduced_omnibus",
            "omnibus_pvalue_BH",
        )
    if contrast["lfc_threshold"] > 0:
        return (
            "DESeq2::results_Wald_greaterAbs",
            "Wald",
            "abs_log2_fold_change_greater_than_threshold",
            "contrast_threshold_pvalue_BH_after_independent_filtering",
        )
    return (
        "DESeq2::results_Wald",
        "Wald",
        "contrast_equals_zero",
        "contrast_pvalue_BH_after_independent_filtering",
    )


def _verify_deseq2_results(
    rows: Sequence[Sequence[str]],
    *,
    contrasts: Sequence[Mapping[str, Any]],
    analysis: Mapping[str, Any],
    test: Mapping[str, Any],
    coefficient_estimates: Mapping[str, Mapping[str, float]],
    coefficient_states: Mapping[str, set[tuple[str, str]]],
) -> tuple[list[dict[str, Any]], set[str]]:
    if not rows:
        raise _integrity("The DESeq2 results TSV contains no result rows.")
    contrast_by_id = {item["contrast_id"]: item for item in contrasts}
    observed: set[tuple[str, str]] = set()
    genes_by_contrast: dict[str, set[str]] = {
        identifier: set() for identifier in contrast_by_id
    }
    status_counts = {
        identifier: {status: 0 for status in sorted(_STATUSES)}
        for identifier in contrast_by_id
    }
    significance = {
        identifier: {"fdr_le_0_05": 0, "fdr_gt_0_05": 0, "not_tested": 0}
        for identifier in contrast_by_id
    }
    adjusted_probabilities: dict[str, list[tuple[float, float, int]]] = {
        identifier: [] for identifier in contrast_by_id
    }
    apeglm_nonconvergence_counts = {
        identifier: 0 for identifier in contrast_by_id
    }
    for row_number, row in enumerate(rows, start=2):
        values = dict(zip(_DESEQ2_RESULT_HEADER, row, strict=True))
        gene = values["gene_id"]
        identifier = values["contrast_id"]
        status = values["status"]
        reason = values["status_reason"]
        if (
            not gene
            or has_control_characters(gene)
            or identifier not in contrast_by_id
            or status not in _STATUSES
        ):
            raise _integrity("A DESeq2 result identity or status is invalid.", row=row_number)
        if status == "not_estimable":
            raise _integrity(
                "A completed DESeq2 bundle cannot contain a non-estimable row.",
                row=row_number,
            )
        key = (gene, identifier)
        if key in observed:
            raise _integrity(
                "The DESeq2 results contain a duplicate gene/contrast row.",
                row=row_number,
            )
        observed.add(key)
        genes_by_contrast[identifier].add(gene)
        contrast = contrast_by_id[identifier]
        threshold = _number(
            values["lfc_threshold"],
            role="results.tsv",
            row=row_number,
            field="lfc_threshold",
        )
        expected_method, expected_statistic, expected_hypothesis, expected_fdr_basis = (
            _deseq2_expected_result_identity(contrast, test=test)
        )
        if (
            threshold != contrast["lfc_threshold"]
            or values["test_method"] != expected_method
            or values["statistic_type"] != expected_statistic
            or values["statistic_hypothesis"] != expected_hypothesis
            or values["fdr_basis"] != expected_fdr_basis
            or values["shrinkage_method"] != contrast["shrinkage_method"]
        ):
            raise _integrity(
                "A DESeq2 result row disagrees with its declared test.", row=row_number
            )

        parsed = {
            field: _optional_number(
                values[field], role="results.tsv", row=row_number, field=field
            )
            for field in (
                "baseMean",
                "logFC",
                "unshrunk_logFC",
                "lfcSE",
                "statistic",
                "PValue",
                "FDR",
            )
        }
        base_mean = parsed["baseMean"]
        if base_mean is None or base_mean < 0:
            raise _integrity("A DESeq2 baseMean is missing or negative.", row=row_number)
        for probability in ("PValue", "FDR"):
            number = parsed[probability]
            if number is not None and not 0 <= number <= 1:
                raise _integrity(
                    "A DESeq2 probability is outside [0, 1].",
                    row=row_number,
                    field=probability,
                )
        if parsed["lfcSE"] is not None and parsed["lfcSE"] < 0:
            raise _integrity("A DESeq2 LFC standard error is negative.", row=row_number)
        if (
            test["mode"] == "lrt"
            and parsed["statistic"] is not None
            and parsed["statistic"] < 0
            and (
                parsed["PValue"] != 1
                or parsed["FDR"] not in {None, 1}
            )
        ):
            raise _integrity(
                "A negative raw DESeq2 LRT statistic must carry PValue=1 "
                "and either missing or unit FDR.",
                row=row_number,
            )

        if test["shrinkage"] == "none":
            if (parsed["unshrunk_logFC"] is None) != (parsed["logFC"] is None) or (
                parsed["unshrunk_logFC"] is not None
                and not math.isclose(
                    parsed["unshrunk_logFC"],
                    parsed["logFC"],
                    rel_tol=1e-12,
                    abs_tol=1e-12,
                )
            ):
                raise _integrity(
                    "An unshrunk DESeq2 LFC must equal logFC when shrinkage is disabled.",
                    row=row_number,
                )
        elif status in {"tested", "filtered"} or reason in {
            "cooks_outlier",
            "apeglm_nonconvergence",
        }:
            if parsed["unshrunk_logFC"] is None:
                raise _integrity(
                    "An apeglm result must preserve the unshrunk DESeq2 LFC.",
                    row=row_number,
                )
        if reason == "apeglm_nonconvergence" and parsed["logFC"] is not None:
            raise _integrity(
                "An apeglm failure cannot silently fall back to an unshrunk LFC.",
                row=row_number,
            )
        if reason == "apeglm_nonconvergence":
            apeglm_nonconvergence_counts[identifier] += 1

        required_model_outcomes = ("logFC", "lfcSE", "statistic", "PValue")
        if status == "tested":
            if reason != "tested" or any(parsed[field] is None for field in required_model_outcomes):
                raise _integrity("A tested DESeq2 row is incomplete.", row=row_number)
        elif status == "filtered":
            if (
                reason != "independent_filtering"
                or any(parsed[field] is None for field in required_model_outcomes)
                or parsed["FDR"] is not None
                or base_mean >= contrast["independent_filter_threshold"]
            ):
                raise _integrity("A filtered DESeq2 row is inconsistent.", row=row_number)
        elif status == "not_tested":
            if (
                reason != "all_zero"
                or base_mean != 0
                or any(
                    parsed[field] is not None
                    for field in (
                        "logFC",
                        "unshrunk_logFC",
                        "lfcSE",
                        "statistic",
                        "PValue",
                        "FDR",
                    )
                )
            ):
                raise _integrity("An all-zero DESeq2 row is inconsistent.", row=row_number)
        else:
            if reason not in _DESEQ2_RESULT_FAILURE_REASONS:
                raise _integrity("A failed DESeq2 row has an unknown reason.", row=row_number)
            if reason == "cooks_outlier" and (
                parsed["PValue"] is not None or parsed["FDR"] is not None
            ):
                raise _integrity("A Cook's-outlier row cannot carry PValue or FDR.")

        gene_coefficient_states = coefficient_states.get(gene, set())
        if status == "not_tested":
            expected_coefficient_states = {("not_tested", "all_zero")}
            if gene_coefficient_states != expected_coefficient_states:
                raise _integrity(
                    "An all-zero DESeq2 result disagrees with coefficients.tsv.",
                    row=row_number,
                )
        elif reason in _DESEQ2_MODEL_FAILURE_REASONS:
            if gene_coefficient_states != {("failed", reason)}:
                raise _integrity(
                    "A model-failed DESeq2 result disagrees with coefficients.tsv.",
                    row=row_number,
                )
        elif any(
            coefficient_reason == "all_zero"
            or coefficient_reason in _DESEQ2_MODEL_FAILURE_REASONS
            for _, coefficient_reason in gene_coefficient_states
        ):
            raise _integrity(
                "A DESeq2 result status disagrees with the gene model state.",
                row=row_number,
            )

        p_value = parsed["PValue"]
        fdr = parsed["FDR"]
        passes_filter = base_mean >= contrast["independent_filter_threshold"]
        if p_value is None:
            if fdr is not None:
                raise _integrity("A DESeq2 FDR cannot exist without a finite PValue.")
        elif passes_filter:
            if fdr is None:
                raise _integrity(
                    "A finite, filter-passing DESeq2 PValue must carry an FDR.",
                    row=row_number,
                )
        elif fdr is not None:
            raise _integrity(
                "An independently filtered DESeq2 PValue cannot carry an FDR.",
                row=row_number,
            )
        if fdr is not None:
            adjusted_probabilities[identifier].append((p_value, fdr, row_number))
            significance[identifier][
                "fdr_le_0_05" if fdr <= 0.05 else "fdr_gt_0_05"
            ] += 1
        else:
            significance[identifier]["not_tested"] += 1

        gene_coefficients = coefficient_estimates.get(gene)
        reported_unshrunk = parsed["unshrunk_logFC"]
        if reported_unshrunk is not None:
            available_coefficients = gene_coefficients or {}
            missing_coefficients = sorted(
                set(contrast["weights"]) - set(available_coefficients)
            )
            if missing_coefficients:
                raise _integrity(
                    "A DESeq2 reporting effect cannot be reproduced from incomplete "
                    "coefficients.",
                    row=row_number,
                    missing_coefficients=missing_coefficients,
                )
            expected_effect = sum(
                available_coefficients[coefficient] * weight
                for coefficient, weight in contrast["weights"].items()
            )
            clean_contrast_all_zero = (
                any(weight < 0 for weight in contrast["weights"].values())
                and any(weight > 0 for weight in contrast["weights"].values())
                and reported_unshrunk == 0
                and parsed["statistic"] == 0
                and parsed["PValue"] == 1
            )
            if not clean_contrast_all_zero and not math.isclose(
                reported_unshrunk,
                expected_effect,
                rel_tol=1e-8,
                abs_tol=1e-10,
            ):
                raise _integrity(
                    "A DESeq2 reporting effect does not reproduce from coefficients.",
                    row=row_number,
                    expected=expected_effect,
                    observed=reported_unshrunk,
                )
        status_counts[identifier][status] += 1

    inventories = list(genes_by_contrast.values())
    genes = inventories[0]
    if any(inventory != genes for inventory in inventories[1:]):
        raise _integrity("DESeq2 contrasts do not contain the same gene inventory.")
    expected_pairs = {
        (gene, identifier) for gene in genes for identifier in contrast_by_id
    }
    if observed != expected_pairs:
        raise _integrity("The DESeq2 results are not the complete gene/contrast matrix.")
    for identifier, probabilities in adjusted_probabilities.items():
        expected_fdr = _benjamini_hochberg(
            [p_value for p_value, _, _ in probabilities]
        )
        for (_, observed_fdr, row_number), expected_value in zip(
            probabilities, expected_fdr, strict=True
        ):
            if not math.isclose(
                observed_fdr,
                expected_value,
                rel_tol=1e-10,
                abs_tol=1e-12,
            ):
                raise _integrity(
                    "A DESeq2 FDR does not equal the independently recomputed "
                    "within-contrast Benjamini-Hochberg adjustment.",
                    contrast_id=identifier,
                    row=row_number,
                    observed_fdr=observed_fdr,
                    expected_fdr=expected_value,
                )
    for identifier, observed_count in apeglm_nonconvergence_counts.items():
        expected_count = contrast_by_id[identifier][
            "shrinkage_nonconvergence_count"
        ]
        if observed_count != expected_count:
            raise _integrity(
                "The DESeq2 apeglm convergence count disagrees with results.tsv.",
                contrast_id=identifier,
                observed=observed_count,
                expected=expected_count,
            )
    gene_evidence = analysis["genes"]
    expected_status_counts = {
        status: sum(counts[status] for counts in status_counts.values())
        for status in sorted(_STATUSES)
    }
    if (
        gene_evidence["total"] != len(genes)
        or gene_evidence["result_rows"] != len(rows)
        or gene_evidence["status_counts"] != expected_status_counts
        or analysis["defaults"]["outlier_replacement_count"] > len(genes)
    ):
        raise _integrity("The DESeq2 gene evidence disagrees with results.tsv.")
    summaries: list[dict[str, Any]] = []
    for contrast in contrasts:
        identifier = contrast["contrast_id"]
        summaries.append(
            {
                "contrast_id": identifier,
                "lfc_threshold": contrast["lfc_threshold"],
                "test_method": contrast["test_method"],
                "alternative_hypothesis": contrast["alternative_hypothesis"],
                "status_counts": status_counts[identifier],
                "significance_counts": significance[identifier],
            }
        )
    return summaries, genes


def _pathway_integer(value: str, *, row: int, field: str) -> int:
    if re.fullmatch(r"0|[1-9][0-9]*", value) is None:
        raise _integrity(
            "A pathway result count is not a canonical non-negative integer.",
            role=_PATHWAY_RESULT_FILE,
            row=row,
            field=field,
            observed=value,
        )
    return int(value)


def _json_exact_equal(observed: Any, expected: Any) -> bool:
    if isinstance(expected, Mapping):
        return (
            isinstance(observed, Mapping)
            and set(observed) == set(expected)
            and all(_json_exact_equal(observed[key], expected[key]) for key in expected)
        )
    if isinstance(expected, list):
        return (
            isinstance(observed, list)
            and len(observed) == len(expected)
            and all(
                _json_exact_equal(item, expected_item)
                for item, expected_item in zip(observed, expected, strict=True)
            )
        )
    if expected is None:
        return observed is None
    if isinstance(expected, bool):
        return type(observed) is bool and observed is expected
    if isinstance(expected, int):
        return type(observed) is int and observed == expected
    if isinstance(expected, float):
        return type(observed) is float and observed == expected
    if isinstance(expected, str):
        return type(observed) is str and observed == expected
    return type(observed) is type(expected) and observed == expected


def _pathway_json_string(value: Any, *, context: str) -> str:
    if not isinstance(value, str) or not value or has_control_characters(value):
        raise _integrity(f"The {context} is not a valid non-empty string.")
    return value


def _pathway_json_integer(value: Any, *, context: str, positive: bool = False) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < (1 if positive else 0)
    ):
        raise _integrity(f"The {context} is not a valid integer count.")
    return value


def _pathway_json_digest(value: Any, *, context: str) -> str:
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise _integrity(f"The {context} is not a canonical SHA-256 digest.")
    return value


def _pathway_string_list(value: Any, *, context: str) -> list[str]:
    normalized = [value] if isinstance(value, str) else value
    if (
        not isinstance(normalized, list)
        or any(
            not isinstance(item, str) or not item or has_control_characters(item)
            for item in normalized
        )
        or len(set(normalized)) != len(normalized)
    ):
        raise _integrity(f"The {context} is not a unique string array.")
    return normalized


def _verify_pathway_analysis(
    analysis: Mapping[str, Any], contrasts: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    value = analysis.get("pathway_analysis")
    if not isinstance(value, Mapping):
        raise _integrity("The pathway analysis evidence is missing or invalid.")
    _exact_keys(
        value,
        {
            "gene_sets",
            "mapping_policy",
            "tested_universe_gene_count",
            "methods",
            "multiple_testing",
            "contrasts",
        },
        context="pathway analysis evidence",
    )

    tested_universe = _pathway_json_integer(
        value["tested_universe_gene_count"],
        context="pathway tested-universe size",
        positive=True,
    )
    if tested_universe != analysis["genes"]["tested"]:
        raise _integrity(
            "The pathway tested universe disagrees with the gene-level tested universe."
        )

    gene_sets = value["gene_sets"]
    if not isinstance(gene_sets, Mapping):
        raise _integrity("The pathway gene-set source evidence is invalid.")
    _exact_keys(
        gene_sets,
        {"gmt", "annotation", "minimum_tested_genes", "sets"},
        context="pathway gene-set evidence",
    )
    gmt = gene_sets["gmt"]
    if not isinstance(gmt, Mapping):
        raise _integrity("The pathway GMT evidence is invalid.")
    _exact_keys(
        gmt,
        {
            "collection",
            "version",
            "identifier_type",
            "sha256",
            "size_bytes",
            "gene_set_count",
        },
        context="pathway GMT evidence",
    )
    _pathway_json_string(gmt["collection"], context="pathway GMT collection")
    _pathway_json_string(gmt["version"], context="pathway GMT version")
    if gmt["identifier_type"] != "symbol":
        raise _integrity("The pathway GMT identifier type is incompatible.")
    _pathway_json_digest(gmt["sha256"], context="pathway GMT digest")
    _pathway_json_integer(gmt["size_bytes"], context="pathway GMT size", positive=True)
    gene_set_count = _pathway_json_integer(
        gmt["gene_set_count"], context="pathway gene-set count", positive=True
    )

    annotation = gene_sets["annotation"]
    if not isinstance(annotation, Mapping):
        raise _integrity("The pathway annotation evidence is invalid.")
    _exact_keys(
        annotation,
        {
            "name",
            "version",
            "gene_id_column",
            "symbol_column",
            "sha256",
            "size_bytes",
            "row_count",
        },
        context="pathway annotation evidence",
    )
    _pathway_json_string(annotation["name"], context="pathway annotation name")
    _pathway_json_string(annotation["version"], context="pathway annotation version")
    if (
        annotation["gene_id_column"] != "gene_id"
        or annotation["symbol_column"] != "symbol"
    ):
        raise _integrity("The pathway annotation column identity is incompatible.")
    _pathway_json_digest(annotation["sha256"], context="pathway annotation digest")
    _pathway_json_integer(
        annotation["size_bytes"], context="pathway annotation size", positive=True
    )
    _pathway_json_integer(
        annotation["row_count"], context="pathway annotation row count", positive=True
    )
    input_evidence = analysis["input_evidence"]
    snapshots = input_evidence.get("r_input_snapshots")
    if not isinstance(snapshots, list):
        raise _integrity("The pathway input snapshot inventory is invalid.")
    snapshots_by_role = {
        item.get("role"): item for item in snapshots if isinstance(item, Mapping)
    }
    if set(snapshots_by_role) & {"pathways.gmt", "pathways.annotation"} != {
        "pathways.gmt",
        "pathways.annotation",
    }:
        raise _integrity(
            "The pathway source snapshots are missing from input evidence."
        )
    for role, source in (
        ("pathways.gmt", gmt),
        ("pathways.annotation", annotation),
    ):
        snapshot = snapshots_by_role[role]
        if (
            snapshot.get("sha256") != source["sha256"]
            or snapshot.get("size_bytes") != source["size_bytes"]
        ):
            raise _integrity(
                "A pathway source identity disagrees with its captured input snapshot.",
                role=role,
            )
    minimum = _pathway_json_integer(
        gene_sets["minimum_tested_genes"],
        context="minimum tested genes per pathway",
        positive=True,
    )

    sets = gene_sets["sets"]
    if not isinstance(sets, list) or not sets:
        raise _integrity("The declared pathway set inventory is invalid or empty.")
    expected_set_keys = {
        "gene_set_id",
        "gene_set_description",
        "gmt_member_count_raw",
        "gmt_symbol_count_unique",
        "mapped_symbol_count_unique",
        "ambiguous_symbol_count_unique",
        "unmapped_symbol_count_unique",
        "mapping_rate",
        "mapped_gene_id_count_unique",
        "tested_gene_count",
        "filtered_gene_count",
        "absent_from_assay_gene_count",
    }
    normalized_sets: list[dict[str, Any]] = []
    set_ids: list[str] = []
    for index, item in enumerate(sets):
        if not isinstance(item, Mapping):
            raise _integrity(
                "A declared pathway set record is invalid.", set_index=index
            )
        _exact_keys(item, expected_set_keys, context=f"pathway set {index}")
        gene_set_id = _pathway_json_string(
            item["gene_set_id"], context=f"pathway set {index} identifier"
        )
        description = _pathway_json_string(
            item["gene_set_description"],
            context=f"pathway set {index} description",
        )
        count_fields = expected_set_keys - {
            "gene_set_id",
            "gene_set_description",
            "mapping_rate",
        }
        counts = {
            field: _pathway_json_integer(
                item[field], context=f"pathway set {index} {field}"
            )
            for field in count_fields
        }
        mapping_rate = _finite_json_number(
            item["mapping_rate"],
            context=f"pathway set {index} mapping rate",
            nonnegative=True,
        )
        unique_count = counts["gmt_symbol_count_unique"]
        expected_mapping_rate = (
            counts["mapped_symbol_count_unique"] / unique_count if unique_count else 0.0
        )
        if (
            counts["gmt_member_count_raw"] <= 0
            or counts["gmt_member_count_raw"] != unique_count
            or unique_count
            != counts["mapped_symbol_count_unique"]
            + counts["ambiguous_symbol_count_unique"]
            + counts["unmapped_symbol_count_unique"]
            or counts["mapped_gene_id_count_unique"]
            != counts["tested_gene_count"]
            + counts["filtered_gene_count"]
            + counts["absent_from_assay_gene_count"]
            or counts["mapped_gene_id_count_unique"]
            != counts["mapped_symbol_count_unique"]
            or counts["tested_gene_count"] > tested_universe
            or counts["filtered_gene_count"] > analysis["genes"]["filtered"]
            or mapping_rate > 1
            or not math.isclose(
                mapping_rate,
                expected_mapping_rate,
                rel_tol=1e-10,
                abs_tol=1e-12,
            )
        ):
            raise _integrity(
                "A declared pathway set has inconsistent mapping counts.",
                gene_set_id=gene_set_id,
            )
        set_ids.append(gene_set_id)
        normalized_sets.append(
            {
                "gene_set_id": gene_set_id,
                "gene_set_description": description,
                **counts,
                "mapping_rate": mapping_rate,
            }
        )
    if (
        len(set(set_ids)) != len(set_ids)
        or set_ids != sorted(set_ids)
        or gene_set_count != len(set_ids)
    ):
        raise _integrity(
            "The pathway set inventory is duplicated, unordered, or disagrees with its count."
        )

    mapping_policy = value["mapping_policy"]
    if not isinstance(mapping_policy, Mapping):
        raise _integrity("The pathway mapping policy is invalid.")
    _exact_keys(
        mapping_policy,
        {
            "source_identifier",
            "target_identifier",
            "annotation_gene_id_version_stripping",
            "one_to_many_symbols",
            "duplicate_gmt_members",
            "mapping_rate",
            "tested_membership",
            "not_tested_policy",
        },
        context="pathway mapping policy",
    )
    expected_mapping_policy = {
        "source_identifier": "symbol",
        "target_identifier": "stable_gene_id",
        "one_to_many_symbols": "ambiguous_excluded",
        "duplicate_gmt_members": "hard_fail",
        "mapping_rate": (
            "uniquely_mapped_unique_symbols divided by unique_GMT_symbols"
        ),
        "tested_membership": "intersection_with_filterByExpr_tested_universe",
        "not_tested_policy": "tested_gene_count_below_minimum",
    }
    if not isinstance(
        mapping_policy["annotation_gene_id_version_stripping"], bool
    ) or any(
        mapping_policy[key] != expected
        for key, expected in expected_mapping_policy.items()
    ):
        raise _integrity("The pathway mapping policy is incompatible.")

    expected_methods = {
        "limma_mroast": {
            "generic": "limma::mroast",
            "dispatch": "edgeR::mroast.DGEGLM",
            "test_class": "self_contained",
            "inference_role": "corroborative",
            "statistical_null": "zero_effect",
            "parameters": {
                "set.statistic": "mean",
                "gene.weights": None,
                "nrot": 9999,
                "adjust.method": "BH",
                "midp": False,
                "sort": "none",
            },
        },
        "limma_fry": {
            "generic": "limma::fry",
            "dispatch": "edgeR::fry.DGEGLM",
            "test_class": "self_contained",
            "inference_role": "primary",
            "statistical_null": "zero_effect",
            "parameters": {"sort": "none"},
        },
        "limma_camera": {
            "generic": "limma::camera",
            "dispatch": "edgeR::camera.DGEGLM",
            "test_class": "competitive",
            "inference_role": "supplementary",
            "statistical_null": "no_enrichment_relative_to_complement",
            "parameters": {
                "weights": None,
                "use.ranks": False,
                "allow.neg.cor": False,
                "inter.gene.cor": "estimated_per_gene_set",
                "sort": False,
            },
        },
    }
    if not _json_exact_equal(value["methods"], expected_methods):
        raise _integrity("The pathway method policy is incompatible.")
    expected_multiple_testing = {
        "method": "Benjamini-Hochberg",
        "scope": (
            "separately within contrast x method x hypothesis across tested "
            "gene sets only"
        ),
        "family_id_format": "contrast_id|method_id|hypothesis",
        "python_independent_recalculation_required": True,
    }
    if not _json_exact_equal(value["multiple_testing"], expected_multiple_testing):
        raise _integrity("The pathway multiple-testing policy is incompatible.")

    pathway_contrasts = value["contrasts"]
    if not isinstance(pathway_contrasts, list) or len(pathway_contrasts) != len(
        contrasts
    ):
        raise _integrity("The pathway contrast evidence is incomplete.")
    expected_self_sets = [
        item["gene_set_id"]
        for item in normalized_sets
        if item["tested_gene_count"] >= minimum
    ]
    expected_competitive_sets = [
        item["gene_set_id"]
        for item in normalized_sets
        if item["tested_gene_count"] >= minimum
        and item["tested_gene_count"] >= 2
        and item["tested_gene_count"] < tested_universe
    ]
    ordered_membership_digests: dict[str, str] = {}
    for index, (item, contrast) in enumerate(
        zip(pathway_contrasts, contrasts, strict=True)
    ):
        if not isinstance(item, Mapping):
            raise _integrity(
                "A pathway contrast evidence record is invalid.", contrast_index=index
            )
        _exact_keys(
            item,
            {
                "contrast_id",
                "gene_level_lfc_threshold",
                "self_contained_statistical_null",
                "gene_level_lfc_threshold_applied_to_pathways",
                "ordered_set_lists",
                "rotation",
            },
            context=f"pathway contrast {index}",
        )
        threshold = _finite_json_number(
            item["gene_level_lfc_threshold"],
            context=f"pathway contrast {index} gene-level LFC threshold",
            nonnegative=True,
        )
        if (
            item["contrast_id"] != contrast["contrast_id"]
            or threshold != contrast["lfc_threshold"]
            or item["self_contained_statistical_null"] != "zero_effect"
            or item["gene_level_lfc_threshold_applied_to_pathways"] is not False
        ):
            raise _integrity(
                "A pathway contrast disagrees with its gene-level contrast identity.",
                contrast_index=index,
            )
        ordered = item["ordered_set_lists"]
        if not isinstance(ordered, Mapping):
            raise _integrity("A pathway ordered-set inventory is invalid.")
        _exact_keys(
            ordered, {"self_contained", "competitive"}, context="ordered set lists"
        )
        for test_class, expected_ids in (
            ("self_contained", expected_self_sets),
            ("competitive", expected_competitive_sets),
        ):
            record = ordered[test_class]
            if not isinstance(record, Mapping):
                raise _integrity("A pathway ordered-set record is invalid.")
            _exact_keys(
                record,
                {"gene_set_ids", "sha256", "verification_scope"},
                context=f"{test_class} ordered set list",
            )
            observed_ids = (
                _pathway_string_list(
                    record["gene_set_ids"], context=f"{test_class} ordered set IDs"
                )
                if expected_ids
                else ([] if record["gene_set_ids"] == [] else None)
            )
            if observed_ids is None or observed_ids != expected_ids:
                raise _integrity(
                    "A pathway ordered-set list disagrees with set eligibility.",
                    test_class=test_class,
                    contrast_id=contrast["contrast_id"],
                )
            digest = _pathway_json_digest(
                record["sha256"], context=f"{test_class} ordered-set digest"
            )
            if record["verification_scope"] != (
                "backend_execution_evidence_syntax_and_"
                "cross_contrast_consistency_only"
            ):
                raise _integrity(
                    "A pathway ordered-set digest overstates its verification scope.",
                    test_class=test_class,
                    contrast_id=contrast["contrast_id"],
                )
            previous_digest = ordered_membership_digests.setdefault(
                test_class, digest
            )
            if previous_digest != digest:
                raise _integrity(
                    "A pathway ordered-set membership digest changes across contrasts.",
                    test_class=test_class,
                    contrast_id=contrast["contrast_id"],
                )
        rotation = item["rotation"]
        expected_rotation = {
            "method_id": "limma_mroast",
            "seed": 1729 + index,
            "seed_policy": "base_1729_plus_zero_based_declared_contrast_index",
            "rng": {
                "kind": "Mersenne-Twister",
                "normal.kind": "Inversion",
                "sample.kind": "Rejection",
            },
            "nrot": 9999,
            "reset_before_each_contrast": True,
        }
        if not _json_exact_equal(rotation, expected_rotation):
            raise _integrity("The pathway rotation policy is incompatible.")

    return {
        "sets": normalized_sets,
        "minimum_tested_genes": minimum,
    }


def _pathway_probability(value: str, *, row: int, field: str) -> float:
    parsed = _number(value, role=_PATHWAY_RESULT_FILE, row=row, field=field)
    if not 0 <= parsed <= 1:
        raise _integrity(
            "A pathway probability is outside [0, 1].",
            row=row,
            field=field,
        )
    return parsed


def _pathway_expected_status(
    *,
    method_id: str,
    tested_gene_count: int,
    tested_universe_gene_count: int,
    minimum_tested_genes: int,
) -> tuple[str, str]:
    if tested_gene_count < minimum_tested_genes:
        return "not_tested", "tested_gene_count_below_minimum"
    if method_id == "limma_camera" and tested_gene_count < 2:
        return "not_tested", "camera_requires_at_least_two_tested_genes"
    if method_id == "limma_camera" and tested_gene_count == tested_universe_gene_count:
        return "not_tested", "competitive_test_requires_background_genes"
    return "tested", ""


def _pathway_applicability(*, method_id: str, status: str) -> tuple[str, str]:
    if method_id == "limma_mroast":
        rotation = (
            "performed_fixed_seed_9999_rotations"
            if status == "tested"
            else "not_performed_not_tested"
        )
        return "not_applicable", rotation
    if method_id == "limma_fry":
        return "not_applicable", "not_applicable_analytic_approximation"
    correlation = (
        "estimated_set_specific" if status == "tested" else "not_estimated_not_tested"
    )
    return correlation, "not_applicable"


def _empty_pathway_method_summary() -> dict[str, Any]:
    return {
        "status_counts": {"tested": 0, "not_tested": 0},
        "significance_counts": {
            "fdr_le_0_05": 0,
            "fdr_gt_0_05": 0,
            "not_tested": 0,
        },
    }


def _verify_pathway_results(
    rows: Sequence[Sequence[str]],
    *,
    contrasts: Sequence[Mapping[str, Any]],
    analysis: Mapping[str, Any],
    pathway_specification: Mapping[str, Any],
) -> dict[str, Any]:
    contrast_ids = [str(item["contrast_id"]) for item in contrasts]
    declared_sets_raw = pathway_specification["sets"]
    declared_sets = {str(item["gene_set_id"]): dict(item) for item in declared_sets_raw}
    minimum_tested_genes = int(pathway_specification["minimum_tested_genes"])
    gene_evidence = analysis["genes"]
    expected_tested_universe = int(gene_evidence["tested"])
    expected_filtered_universe = int(gene_evidence["filtered"])

    expected_keys = {
        (contrast_id, gene_set_id, method_id, hypothesis)
        for contrast_id in contrast_ids
        for gene_set_id in declared_sets
        for method_id, specification in _PATHWAY_METHODS.items()
        for hypothesis in specification["hypotheses"]
    }
    observed_keys: set[tuple[str, str, str, str]] = set()
    audit_by_set: dict[str, tuple[Any, ...]] = {}
    mroast_proportions: dict[tuple[str, str], tuple[float, float]] = {}
    fdr_groups: dict[tuple[str, str, str], list[tuple[float, float, int]]] = {}
    summaries = {
        "self_contained": {
            "limma_mroast": {
                hypothesis: _empty_pathway_method_summary()
                for hypothesis in ("directional", "mixed")
            },
            "limma_fry": {
                hypothesis: _empty_pathway_method_summary()
                for hypothesis in ("directional", "mixed")
            },
        },
        "competitive": {
            "limma_camera": {"directional": _empty_pathway_method_summary()}
        },
    }

    for row_number, row in enumerate(rows, start=2):
        values = dict(zip(_PATHWAY_HEADER, row, strict=True))
        contrast_id = values["contrast_id"]
        gene_set_id = values["gene_set_id"]
        description = values["gene_set_description"]
        method_id = values["method_id"]
        hypothesis = values["hypothesis"]
        status = values["status"]
        if (
            contrast_id not in contrast_ids
            or gene_set_id not in declared_sets
            or description
            != (
                declared_sets[gene_set_id]["gene_set_description"]
                if gene_set_id in declared_sets
                else None
            )
            or not description
            or has_control_characters(description)
            or method_id not in _PATHWAY_METHODS
            or status not in _PATHWAY_STATUSES
        ):
            raise _integrity(
                "A pathway row identity, description, or status is invalid.",
                row=row_number,
            )
        method = _PATHWAY_METHODS[method_id]
        if (
            hypothesis not in method["hypotheses"]
            or values["test_class"] != method["test_class"]
            or values["inference_role"] != method["inference_role"]
        ):
            raise _integrity(
                "A pathway row mislabels its method, test class, hypothesis, or inference role.",
                row=row_number,
            )
        key = (contrast_id, gene_set_id, method_id, hypothesis)
        if key in observed_keys:
            raise _integrity(
                "The pathway results contain a duplicate result cell.", row=row_number
            )
        observed_keys.add(key)
        expected_family = f"{contrast_id}|{method_id}|{hypothesis}"
        if values["fdr_family_id"] != expected_family:
            raise _integrity(
                "A pathway row has an invalid FDR family identity.",
                row=row_number,
                expected_fdr_family_id=expected_family,
            )

        count_fields = (
            "gmt_member_count_raw",
            "gmt_symbol_count_unique",
            "mapped_symbol_count_unique",
            "ambiguous_symbol_count_unique",
            "unmapped_symbol_count_unique",
            "mapped_gene_id_count_unique",
            "tested_gene_count",
            "filtered_gene_count",
            "tested_universe_gene_count",
        )
        counts = {
            field: _pathway_integer(values[field], row=row_number, field=field)
            for field in count_fields
        }
        mapping_rate = _pathway_probability(
            values["mapping_rate"], row=row_number, field="mapping_rate"
        )
        raw_count = counts["gmt_member_count_raw"]
        unique_count = counts["gmt_symbol_count_unique"]
        mapped_symbols = counts["mapped_symbol_count_unique"]
        ambiguous_symbols = counts["ambiguous_symbol_count_unique"]
        unmapped_symbols = counts["unmapped_symbol_count_unique"]
        mapped_gene_ids = counts["mapped_gene_id_count_unique"]
        tested_genes = counts["tested_gene_count"]
        filtered_genes = counts["filtered_gene_count"]
        tested_universe = counts["tested_universe_gene_count"]
        expected_mapping_rate = mapped_symbols / unique_count if unique_count else 0.0
        declared = declared_sets[gene_set_id]
        declared_count_fields = count_fields[:-1]
        if (
            raw_count <= 0
            or raw_count != unique_count
            or unique_count != mapped_symbols + ambiguous_symbols + unmapped_symbols
            or mapped_gene_ids < tested_genes + filtered_genes
            or tested_universe != expected_tested_universe
            or tested_universe <= 0
            or tested_genes > tested_universe
            or filtered_genes > expected_filtered_universe
            or any(counts[field] != declared[field] for field in declared_count_fields)
            or not math.isclose(
                mapping_rate,
                float(declared["mapping_rate"]),
                rel_tol=1e-10,
                abs_tol=1e-12,
            )
            or not math.isclose(
                mapping_rate,
                expected_mapping_rate,
                rel_tol=1e-10,
                abs_tol=1e-12,
            )
        ):
            raise _integrity(
                "A pathway row has inconsistent mapping or tested-universe counts.",
                row=row_number,
            )
        audit = (
            description,
            *(counts[field] for field in count_fields),
            mapping_rate,
        )
        previous_audit = audit_by_set.setdefault(gene_set_id, audit)
        if previous_audit != audit:
            raise _integrity(
                "A gene-set mapping audit changes across pathway result rows.",
                row=row_number,
                gene_set_id=gene_set_id,
            )

        expected_status, expected_reason = _pathway_expected_status(
            method_id=method_id,
            tested_gene_count=tested_genes,
            tested_universe_gene_count=tested_universe,
            minimum_tested_genes=minimum_tested_genes,
        )
        if status != expected_status or values["status_reason"] != expected_reason:
            raise _integrity(
                "A pathway status or reason disagrees with the fixed eligibility policy.",
                row=row_number,
                expected_status=expected_status,
                expected_status_reason=expected_reason,
            )
        expected_correlation_status, expected_rotation_status = _pathway_applicability(
            method_id=method_id, status=status
        )
        if (
            values["correlation_status"] != expected_correlation_status
            or values["rotation_status"] != expected_rotation_status
        ):
            raise _integrity(
                "A pathway row has incompatible correlation or rotation metadata.",
                row=row_number,
            )

        outcome_fields = (
            "direction",
            "proportion_down",
            "proportion_up",
            "p_value",
            "fdr",
            "method_ngenes",
            "correlation_estimate_raw",
            "correlation_effective",
            "vif_used",
        )
        summary = summaries[method["test_class"]][method_id][hypothesis]
        summary["status_counts"][status] += 1
        if status == "not_tested":
            if any(values[field] for field in outcome_fields):
                raise _integrity(
                    "A not-tested pathway row must have blank inferential outcomes.",
                    row=row_number,
                )
            summary["significance_counts"]["not_tested"] += 1
            continue

        expected_direction = "Mixed" if hypothesis == "mixed" else None
        if expected_direction is None:
            if values["direction"] not in {"Up", "Down"}:
                raise _integrity(
                    "A tested directional pathway row has an invalid direction.",
                    row=row_number,
                )
        elif values["direction"] != expected_direction:
            raise _integrity(
                "A tested mixed pathway row must use direction 'Mixed'.",
                row=row_number,
            )
        p_value = _pathway_probability(
            values["p_value"], row=row_number, field="p_value"
        )
        fdr = _pathway_probability(values["fdr"], row=row_number, field="fdr")
        method_ngenes = _pathway_integer(
            values["method_ngenes"], row=row_number, field="method_ngenes"
        )
        if method_ngenes != tested_genes:
            raise _integrity(
                "A pathway method gene count disagrees with its tested gene count.",
                row=row_number,
            )
        if method_id == "limma_mroast":
            proportion_down = _pathway_probability(
                values["proportion_down"],
                row=row_number,
                field="proportion_down",
            )
            proportion_up = _pathway_probability(
                values["proportion_up"], row=row_number, field="proportion_up"
            )
            if proportion_down + proportion_up > 1 + 1e-12:
                raise _integrity(
                    "The mroast active proportions sum to more than one.",
                    row=row_number,
                )
            proportion_key = (contrast_id, gene_set_id)
            proportions = (proportion_down, proportion_up)
            previous_proportions = mroast_proportions.setdefault(
                proportion_key, proportions
            )
            if previous_proportions != proportions:
                raise _integrity(
                    "The mroast directional and mixed rows disagree on rotation proportions.",
                    row=row_number,
                )
        elif values["proportion_down"] or values["proportion_up"]:
            raise _integrity(
                "A non-mroast pathway row carries mroast-only rotation proportions.",
                row=row_number,
            )
        if method_id == "limma_camera":
            raw_correlation = _number(
                values["correlation_estimate_raw"],
                role=_PATHWAY_RESULT_FILE,
                row=row_number,
                field="correlation_estimate_raw",
            )
            effective_correlation = _number(
                values["correlation_effective"],
                role=_PATHWAY_RESULT_FILE,
                row=row_number,
                field="correlation_effective",
            )
            vif_used = _number(
                values["vif_used"],
                role=_PATHWAY_RESULT_FILE,
                row=row_number,
                field="vif_used",
            )
            expected_effective = max(0.0, raw_correlation)
            expected_vif = max(1.0, 1.0 + (method_ngenes - 1) * raw_correlation)
            lower_bound = -1.0 / (method_ngenes - 1)
            if (
                raw_correlation < lower_bound - 1e-10
                or raw_correlation > 1 + 1e-10
                or not math.isclose(
                    effective_correlation,
                    expected_effective,
                    rel_tol=1e-10,
                    abs_tol=1e-12,
                )
                or not math.isclose(
                    vif_used,
                    expected_vif,
                    rel_tol=1e-10,
                    abs_tol=1e-12,
                )
            ):
                raise _integrity(
                    "A camera correlation or variance-inflation audit is inconsistent.",
                    row=row_number,
                )
        elif any(
            values[field]
            for field in (
                "correlation_estimate_raw",
                "correlation_effective",
                "vif_used",
            )
        ):
            raise _integrity(
                "A non-camera pathway row carries camera-only correlation outcomes.",
                row=row_number,
            )

        group = (contrast_id, method_id, hypothesis)
        fdr_groups.setdefault(group, []).append((p_value, fdr, row_number))
        significance_key = "fdr_le_0_05" if fdr <= 0.05 else "fdr_gt_0_05"
        summary["significance_counts"][significance_key] += 1

    if observed_keys != expected_keys:
        raise _integrity(
            "The pathway results are not the complete declared contrast/set/method grid.",
            missing_cells=[
                "|".join(cell) for cell in sorted(expected_keys - observed_keys)
            ],
            unexpected_cells=[
                "|".join(cell) for cell in sorted(observed_keys - expected_keys)
            ],
        )
    for (contrast_id, method_id, hypothesis), probabilities in fdr_groups.items():
        expected_fdr = _benjamini_hochberg([p_value for p_value, _, _ in probabilities])
        for (_, observed_fdr, row_number), expected_value in zip(
            probabilities, expected_fdr, strict=True
        ):
            if not math.isclose(
                observed_fdr,
                expected_value,
                rel_tol=1e-10,
                abs_tol=1e-12,
            ):
                raise _integrity(
                    "A pathway FDR does not equal the independently recomputed "
                    "within-contrast/method/hypothesis Benjamini-Hochberg adjustment.",
                    contrast_id=contrast_id,
                    method_id=method_id,
                    hypothesis=hypothesis,
                    row=row_number,
                    observed_fdr=observed_fdr,
                    expected_fdr=expected_value,
                )
    return {
        "gene_set_count": len(declared_sets),
        "pathway_result_row_count": len(rows),
        **summaries,
    }


def _summarize_deseq2(
    resolved: Path,
    captured: Mapping[str, tuple[bytes, str, int]],
    display_dir: Path | None,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    if display_dir is not None:
        raise _integrity("A D1 DESeq2 bundle cannot contain an edgeR display sidecar.")
    if _PATHWAY_RESULT_FILE in captured:
        raise _integrity("A D1 DESeq2 bundle cannot contain edgeR pathway results.")
    _verify_deseq2_manifest(manifest, captured)
    analysis = _json_object(captured["analysis.json"][0], role="DESeq2 analysis")
    contrasts, configuration = _verify_deseq2_analysis(analysis, manifest)
    design_rows = _tsv_rows(
        captured["design.tsv"][0],
        expected_header=_DESIGN_HEADER,
        role="design.tsv",
    )
    design_samples, design_columns, coefficient_mapping = _verify_deseq2_design(
        design_rows, analysis
    )
    _verify_deseq2_reporting_effects(
        analysis.get("reporting_effect"),
        contrasts=contrasts,
        coefficient_mapping=coefficient_mapping,
        test=configuration["test"],
    )
    result_rows = _tsv_rows(
        captured["results.tsv"][0],
        expected_header=_DESEQ2_RESULT_HEADER,
        role="results.tsv",
    )
    provisional_genes = {row[0] for row in result_rows if row and row[0]}
    coefficient_rows = _tsv_rows(
        captured["coefficients.tsv"][0],
        expected_header=_DESEQ2_COEFFICIENT_HEADER,
        role="coefficients.tsv",
    )
    coefficient_estimates, coefficient_states = _verify_deseq2_coefficients(
        coefficient_rows,
        columns=design_columns,
        expected_genes=provisional_genes,
    )
    contrast_summaries, genes = _verify_deseq2_results(
        result_rows,
        contrasts=contrasts,
        analysis=analysis,
        test=configuration["test"],
        coefficient_estimates=coefficient_estimates,
        coefficient_states=coefficient_states,
    )
    if genes != provisional_genes:
        raise _integrity("The DESeq2 result gene inventory is inconsistent.")
    _verify_deseq2_route(
        analysis,
        sample_ids=design_samples,
        gene_count=len(genes),
    )
    artifacts = []
    for name in sorted(captured):
        _, digest, size = captured[name]
        artifacts.append(
            {
                "kind": "consumed_analysis_artifact",
                "role": name.removesuffix(".json").removesuffix(".tsv"),
                "path": str(resolved / name),
                "sha256": digest,
                "size_bytes": size,
            }
        )
    evidence = manifest["input_evidence"]
    return {
        "schema_version": _DESEQ2_SCHEMA_VERSION,
        "status": "verified_complete",
        "run_dir": str(resolved),
        "backend": "DESeq2",
        "execution_scope": _DESEQ2_EXECUTION_SCOPE,
        "runtime_identity": dict(_DESEQ2_EXPECTED_RUNTIME),
        "plan_id": evidence["plan_id"],
        "input_semantics": analysis["input_semantics"],
        "test": dict(configuration["test"]),
        "gene_count": len(genes),
        "result_row_count": len(result_rows),
        "contrasts": contrast_summaries,
        "artifacts": artifacts,
    }


def summarize_run(run_dir: str | Path) -> dict[str, Any]:
    """Verify and summarize one supported immutable public result bundle.

    The summary is computed only from bytes captured and verified against the
    backend manifest. It never trusts file presence as evidence of success.
    """

    resolved = _resolve_run_dir(run_dir)
    captured, display_dir = _capture_bundle(resolved)
    manifest = _json_object(
        captured["backend_manifest.json"][0], role="backend_manifest"
    )
    manifest_identity = (manifest.get("kind"), manifest.get("backend"))
    if manifest_identity == ("deseq2_backend_manifest", "DESeq2"):
        return _summarize_deseq2(
            resolved,
            captured,
            display_dir,
            manifest,
        )
    if manifest_identity != ("edger_ql_backend_manifest", "edgeR_QL"):
        raise _integrity(
            "The result bundle declares an unsupported backend identity.",
            kind=manifest.get("kind"),
            backend=manifest.get("backend"),
        )
    _verify_manifest(manifest, captured)
    analysis = _json_object(captured["analysis.json"][0], role="analysis")
    contrasts = _verify_analysis(analysis, manifest)
    design_rows = _tsv_rows(
        captured["design.tsv"][0], expected_header=_DESIGN_HEADER, role="design.tsv"
    )
    design_samples, design_columns = _verify_design(design_rows, analysis)
    result_rows = _tsv_rows(
        captured["results.tsv"][0],
        expected_header=_RESULT_HEADER,
        role="results.tsv",
    )
    contrast_summaries, genes, result_gene_statuses = _verify_results(
        result_rows, contrasts=contrasts, analysis=analysis
    )
    coefficient_rows = _tsv_rows(
        captured["coefficients.tsv"][0],
        expected_header=_COEFFICIENT_HEADER,
        role="coefficients.tsv",
    )
    coefficient_statuses = _verify_coefficients(
        coefficient_rows, columns=design_columns, expected_genes=genes
    )
    if coefficient_statuses != result_gene_statuses:
        raise _integrity("Coefficient and result gene statuses disagree.")

    pathway_summary = None
    if _PATHWAY_RESULT_FILE in captured:
        pathway_specification = _verify_pathway_analysis(analysis, contrasts)
        pathway_rows = _tsv_rows(
            captured[_PATHWAY_RESULT_FILE][0],
            expected_header=_PATHWAY_HEADER,
            role=_PATHWAY_RESULT_FILE,
        )
        pathway_summary = _verify_pathway_results(
            pathway_rows,
            contrasts=contrasts,
            analysis=analysis,
            pathway_specification=pathway_specification,
        )

    artifacts = []
    for name in sorted(captured):
        _, digest, size = captured[name]
        artifacts.append(
            {
                "kind": "consumed_analysis_artifact",
                "role": name.removesuffix(".json").removesuffix(".tsv"),
                "path": str(resolved / name),
                "sha256": digest,
                "size_bytes": size,
            }
        )
    evidence = manifest["input_evidence"]
    summary = {
        "schema_version": (
            PATHWAY_SUMMARY_SCHEMA_VERSION
            if pathway_summary is not None
            else SUMMARY_SCHEMA_VERSION
        ),
        "status": "verified_complete",
        "run_dir": str(resolved),
        "backend": "edgeR_QL",
        "execution_scope": "validated_p0_input",
        "runtime_identity": dict(_EXPECTED_RUNTIME),
        "plan_id": evidence["plan_id"],
        "input_semantics": analysis["input_semantics"],
        "gene_count": len(genes),
        "result_row_count": len(result_rows),
        "contrasts": contrast_summaries,
        "artifacts": artifacts,
    }
    if pathway_summary is not None:
        summary["pathways"] = pathway_summary
    if display_dir is not None:
        from .display_bundle import verify_display_bundle

        display = verify_display_bundle(
            display_dir=display_dir,
            core_dir=resolved,
            core_captured=captured,
            analysis=analysis,
            backend_manifest=manifest,
            expected_sample_ids=design_samples,
        )
        summary["display"] = display["summary"]
        summary["artifacts"].extend(display["artifacts"])
    return summary


__all__ = ["SUMMARY_SCHEMA_VERSION", "summarize_run"]
