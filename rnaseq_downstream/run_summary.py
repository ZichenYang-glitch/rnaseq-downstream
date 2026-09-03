"""Fail-closed reader for a published edgeR P0 result bundle."""

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
_EXPECTED_RUNTIME = {
    "R": "4.6.1",
    "Bioconductor": "3.23",
    "BiocVersion_package": "3.23.1",
    "edgeR": "4.10.0",
    "tximport": "1.40.0",
    "limma": "3.68.0",
}
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


def _tsv_rows(
    content: bytes, *, expected_header: Sequence[str], role: str
) -> list[list[str]]:
    try:
        text = content.decode("utf-8")
        rows = list(
            csv.reader(io.StringIO(text, newline=""), delimiter="\t", strict=True)
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
            > counts["mapped_symbol_count_unique"]
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
                "pathway_statistical_null",
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
            or item["pathway_statistical_null"] != "zero_effect"
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
                {"gene_set_ids", "sha256"},
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
            _pathway_json_digest(
                record["sha256"], context=f"{test_class} ordered-set digest"
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


def summarize_run(run_dir: str | Path) -> dict[str, Any]:
    """Verify and summarize one immutable public P0 result bundle.

    The summary is computed only from bytes captured and verified against the
    backend manifest. It never trusts file presence as evidence of success.
    """

    resolved = _resolve_run_dir(run_dir)
    captured, display_dir = _capture_bundle(resolved)
    manifest = _json_object(
        captured["backend_manifest.json"][0], role="backend_manifest"
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
