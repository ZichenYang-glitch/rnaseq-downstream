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
_CORE_FILES = {
    "analysis.json",
    "backend_manifest.json",
    "coefficients.tsv",
    "design.tsv",
    "results.tsv",
}
_MANIFESTED_FILES = _CORE_FILES - {"backend_manifest.json"}
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
_STATUSES = {"filtered", "not_tested", "not_estimable", "failed", "tested"}
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
    expected_names = _CORE_FILES | ({"display"} if "display" in names else set())
    unsafe = sorted(
        entry.name
        for entry in entries
        if entry.is_symlink()
        or (entry.name == "display" and not entry.is_dir())
        or (entry.name != "display" and not entry.is_file())
    )
    if names != expected_names or unsafe:
        raise _integrity(
            "The result bundle must contain the five core files and optionally one display directory.",
            missing_files=sorted(_CORE_FILES - names),
            unexpected_files=sorted(names - (_CORE_FILES | {"display"})),
            unsafe_files=unsafe,
        )
    captured = {name: _capture(run_dir / name) for name in sorted(_CORE_FILES)}
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
    expected = {
        "schema_version": "1.0",
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
        if name not in _MANIFESTED_FILES or name in observed:
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
    if observed != _MANIFESTED_FILES:
        raise _integrity(
            "The backend manifest omits required result members.",
            missing_members=sorted(_MANIFESTED_FILES - observed),
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
    _exact_keys(analysis, expected_keys, context="analysis document")
    identity = {
        "schema_version": analysis.get("schema_version"),
        "kind": analysis.get("kind"),
        "backend": analysis.get("backend"),
        "execution_scope": analysis.get("execution_scope"),
        "runtime_identity": analysis.get("runtime_identity"),
    }
    expected = {
        "schema_version": "1.0",
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

    artifacts = []
    for name in sorted(_CORE_FILES):
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
        "schema_version": SUMMARY_SCHEMA_VERSION,
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
