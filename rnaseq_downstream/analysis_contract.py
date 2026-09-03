"""Strict analysis-request and validated-bundle contract for edgeR P0.

The statistical backend never accepts raw input paths directly.  An analysis
request points to a complete checkpoint-A validation bundle, whose manifest,
plan identity, eligibility, and member digests are verified before any R
process can be started.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import re
from typing import Any, Mapping

from .errors import (
    ContrastNotEstimableError,
    InputIntegrityError,
    InputReadError,
    InvalidRequestError,
)
from .gene_sets import parse_gene_sets_request
from .provenance import (
    has_control_characters,
    read_json_object,
    require_expected_keys,
    require_nonempty_string,
)
from .validation_bundle import canonical_json_equal, compute_validation_plan_id

# Private R-backend protocol/output schema. Public analysis requests are
# versioned independently so display-only additions cannot silently alter the
# statistical protocol.
ANALYSIS_SCHEMA_VERSION = "1.0"
ANALYSIS_REQUEST_SCHEMA_VERSIONS = ("1.0", "1.1")
VALIDATION_BUNDLE_SCHEMA_VERSION = "1.0"
_BUNDLE_MEMBERS = {
    "input_plan.json",
    "provenance.json",
    "validated_request.json",
}
_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.]*$")
_CONTRAST_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_R_RESERVED_WORDS = frozenset(
    {
        "if",
        "else",
        "repeat",
        "while",
        "function",
        "for",
        "in",
        "next",
        "break",
        "TRUE",
        "FALSE",
        "NULL",
        "Inf",
        "NaN",
        "NA",
        "NA_integer_",
        "NA_real_",
        "NA_complex_",
        "NA_character_",
    }
)

_MANIFEST_KEYS = {"schema_version", "kind", "plan_id", "members"}
_MEMBER_RECORD_KEYS = {"path", "sha256", "size_bytes"}
_PLAN_KEYS = {
    "schema_version",
    "kind",
    "plan_id",
    "validation_scope",
    "validation",
    "input",
}
_RECEIPT_KEYS = {
    "schema_version",
    "kind",
    "plan_id",
    "validation_scope",
    "data",
    "warnings",
    "artifacts",
}
_PROVENANCE_KEYS = {
    "schema_version",
    "kind",
    "plan_id",
    "toolkit",
    "runtime",
    "consumed_artifacts",
    "warnings",
}
_EXPECTED_PLAN_VALIDATION_STATE = {
    "design": "not_run",
    "backend": "not_run",
    "runnable": False,
    "analysis_path_certified": False,
}


@dataclass(frozen=True)
class ValidatedAnalysisRequest:
    """A normalized analysis request plus its verified input evidence."""

    request_path: Path
    request_sha256: str
    request_schema_version: str
    bundle_path: Path
    plan_id: str
    input_data: dict[str, Any]
    consumed_artifacts: tuple[dict[str, Any], ...]
    validation_bundle_artifacts: tuple[dict[str, Any], ...]
    design: dict[str, Any]
    contrasts: tuple[dict[str, Any], ...]
    display: dict[str, Any] | None
    gene_sets: dict[str, Any] | None

    def to_backend_document(self) -> dict[str, Any]:
        """Return the JSON payload completed later with private input snapshots."""

        document = {
            "schema_version": ANALYSIS_SCHEMA_VERSION,
            "kind": "edger_ql_backend_request",
            "analysis_request": {
                "path": str(self.request_path),
                "sha256": self.request_sha256,
            },
            "validated_input_bundle": {
                "path": str(self.bundle_path),
                "plan_id": self.plan_id,
                "artifacts": [dict(item) for item in self.validation_bundle_artifacts],
            },
            "input": dict(self.input_data),
            "design": dict(self.design),
            "contrasts": [dict(item) for item in self.contrasts],
        }
        if self.gene_sets is not None:
            document["gene_sets"] = copy.deepcopy(self.gene_sets)
        return document


def _capture(path: Path) -> tuple[bytes, str, int]:
    """Capture one stable JSON file for digesting and parsing from identical bytes."""

    try:
        resolved = path.resolve(strict=True)
        if path.is_symlink() or not resolved.is_file():
            raise OSError("expected a non-symlink regular file")
        before = resolved.stat()
        with resolved.open("rb") as handle:
            content = handle.read()
        after = resolved.stat()
    except OSError as error:
        raise InputReadError(
            "An analysis JSON evidence file could not be captured.",
            path=path,
            operation="capture_analysis_evidence",
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
    if identity_before != identity_after or len(content) != after.st_size:
        raise InputIntegrityError(
            "An analysis JSON evidence file changed while it was captured.",
            details={"path": str(resolved)},
        )
    return content, hashlib.sha256(content).hexdigest(), len(content)


def _resolve_analysis_request(path: str | Path) -> Path:
    if not isinstance(path, (str, Path)) or not str(path):
        raise InvalidRequestError("The analysis request path must be non-empty.")
    if has_control_characters(str(path)):
        raise InvalidRequestError(
            "The analysis request path contains a control character."
        )
    candidate = Path(path)
    try:
        resolved = candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise InputReadError(
            "The analysis request could not be resolved.",
            path=candidate,
            operation="resolve_analysis_request",
            cause=error,
        ) from error
    if not resolved.is_file():
        error = IsADirectoryError(str(resolved))
        raise InputReadError(
            "The analysis request must be a regular file.",
            path=resolved,
            operation="read_analysis_request",
            cause=error,
        ) from error
    return resolved


def _resolve_bundle(value: Any, *, request_path: Path) -> Path:
    declared = require_nonempty_string(value, field="validated_input_bundle")
    candidate = Path(declared)
    if not candidate.is_absolute():
        candidate = request_path.parent / candidate
    try:
        resolved = candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise InputReadError(
            "The validated input bundle could not be resolved.",
            path=candidate,
            operation="resolve_validation_bundle",
            cause=error,
        ) from error
    if not resolved.is_dir():
        error = NotADirectoryError(str(resolved))
        raise InputReadError(
            "The validated input bundle must be a directory.",
            path=resolved,
            operation="read_validation_bundle",
            cause=error,
        ) from error
    return resolved


def _require_sha256(value: Any, *, field: str) -> str:
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise InputIntegrityError(
            "A validation-bundle digest is not canonical SHA-256.",
            details={"field": field, "observed": value},
        )
    return value


def _require_exact_evidence_schema(
    value: Any,
    *,
    expected_keys: set[str],
    context: str,
) -> Mapping[str, Any]:
    """Require an exact object schema for already-generated bundle evidence."""

    if not isinstance(value, Mapping):
        raise InputIntegrityError(
            "A validation-bundle evidence record is not an object.",
            details={"context": context, "observed_type": type(value).__name__},
        )
    observed_keys = set(value)
    if observed_keys != expected_keys:
        raise InputIntegrityError(
            "A validation-bundle evidence object has an incompatible schema.",
            details={
                "context": context,
                "missing_keys": sorted(expected_keys - observed_keys),
                "unknown_keys": sorted(observed_keys - expected_keys),
            },
        )
    return value


def _require_warning_inventory(value: Any, *, context: str) -> list[Any]:
    if not isinstance(value, list):
        raise InputIntegrityError(
            "A validation warning inventory must be an array.",
            details={"context": context, "observed_type": type(value).__name__},
        )
    for index, warning in enumerate(value):
        warning_context = f"{context}[{index}]"
        record = _require_exact_evidence_schema(
            warning,
            expected_keys={"code", "severity", "message", "details"},
            context=warning_context,
        )
        for field in ("code", "severity", "message"):
            if not isinstance(record[field], str) or not record[field]:
                raise InputIntegrityError(
                    "A validation warning contains an invalid string field.",
                    details={"context": warning_context, "field": field},
                )
        if not isinstance(record["details"], Mapping):
            raise InputIntegrityError(
                "A validation warning details field must be an object.",
                details={"context": warning_context},
            )
    return value


def _read_verified_bundle(
    bundle_path: Path,
) -> tuple[
    str,
    dict[str, Any],
    tuple[dict[str, Any], ...],
    tuple[dict[str, Any], ...],
]:
    expected_entries = _BUNDLE_MEMBERS | {"bundle_manifest.json"}
    try:
        entries = list(bundle_path.iterdir())
    except OSError as error:
        raise InputReadError(
            "The validation bundle directory could not be enumerated.",
            path=bundle_path,
            operation="enumerate_validation_bundle",
            cause=error,
        ) from error
    observed_entries = {entry.name for entry in entries}
    unsafe_entries = sorted(
        entry.name for entry in entries if entry.is_symlink() or not entry.is_file()
    )
    if observed_entries != expected_entries or unsafe_entries:
        raise InputIntegrityError(
            "The validation bundle contains unmanifested or unsafe entries.",
            details={
                "missing_entries": sorted(expected_entries - observed_entries),
                "unexpected_entries": sorted(observed_entries - expected_entries),
                "unsafe_entries": unsafe_entries,
            },
        )
    manifest_path = bundle_path / "bundle_manifest.json"
    manifest_content, manifest_digest, manifest_size = _capture(manifest_path)
    manifest = read_json_object(
        manifest_path,
        document_role="validation_bundle_manifest",
        content=manifest_content,
    )
    _require_exact_evidence_schema(
        manifest,
        expected_keys=_MANIFEST_KEYS,
        context="validation bundle manifest",
    )
    if (
        manifest["schema_version"] != VALIDATION_BUNDLE_SCHEMA_VERSION
        or manifest["kind"] != "validation_bundle_manifest"
    ):
        raise InputIntegrityError(
            "The validation-bundle manifest identity is incompatible.",
            details={
                "schema_version": manifest["schema_version"],
                "kind": manifest["kind"],
            },
        )
    plan_id = _require_sha256(manifest["plan_id"], field="manifest.plan_id")
    members = manifest["members"]
    if not isinstance(members, list):
        raise InputIntegrityError(
            "The validation-bundle member inventory must be an array."
        )

    observed_names: set[str] = set()
    member_documents: dict[str, dict[str, Any]] = {}
    bundle_artifacts: list[dict[str, Any]] = [
        {
            "kind": "consumed_validation_evidence",
            "role": "bundle_manifest",
            "path": str(manifest_path.resolve(strict=True)),
            "sha256": manifest_digest,
            "size_bytes": manifest_size,
        }
    ]
    for index, member in enumerate(members):
        if not isinstance(member, Mapping):
            raise InputIntegrityError(
                "A validation-bundle member record must be an object.",
                details={"member_index": index},
            )
        _require_exact_evidence_schema(
            member,
            expected_keys=_MEMBER_RECORD_KEYS,
            context=f"validation bundle member {index}",
        )
        name = member["path"]
        if not isinstance(name, str) or name not in _BUNDLE_MEMBERS:
            raise InputIntegrityError(
                "The validation bundle contains an unexpected member.",
                details={"member_index": index, "path": name},
            )
        if name in observed_names:
            raise InputIntegrityError(
                "The validation bundle contains a duplicate member record.",
                details={"path": name},
            )
        observed_names.add(name)
        expected_digest = _require_sha256(
            member["sha256"], field=f"manifest.members[{index}].sha256"
        )
        expected_size = member["size_bytes"]
        if (
            isinstance(expected_size, bool)
            or not isinstance(expected_size, int)
            or expected_size < 0
        ):
            raise InputIntegrityError(
                "A validation-bundle member size is invalid.",
                details={"path": name, "size_bytes": expected_size},
            )
        member_path = bundle_path / name
        member_content, observed_digest, observed_size = _capture(member_path)
        if observed_digest != expected_digest or observed_size != expected_size:
            raise InputIntegrityError(
                "A validation-bundle member does not match its manifest.",
                details={
                    "path": str(member_path),
                    "expected_sha256": expected_digest,
                    "observed_sha256": observed_digest,
                    "expected_size_bytes": expected_size,
                    "observed_size_bytes": observed_size,
                },
            )
        member_documents[name] = read_json_object(
            member_path,
            document_role=name.removesuffix(".json"),
            content=member_content,
        )
        bundle_artifacts.append(
            {
                "kind": "consumed_validation_evidence",
                "role": name.removesuffix(".json"),
                "path": str(member_path.resolve(strict=True)),
                "sha256": observed_digest,
                "size_bytes": observed_size,
            }
        )
    if observed_names != _BUNDLE_MEMBERS:
        raise InputIntegrityError(
            "The validation bundle is incomplete.",
            details={
                "missing_members": sorted(_BUNDLE_MEMBERS - observed_names),
                "unexpected_members": sorted(observed_names - _BUNDLE_MEMBERS),
            },
        )

    plan = member_documents["input_plan.json"]
    receipt = member_documents["validated_request.json"]
    provenance = member_documents["provenance.json"]
    _require_exact_evidence_schema(
        plan,
        expected_keys=_PLAN_KEYS,
        context="input_plan.json",
    )
    _require_exact_evidence_schema(
        receipt,
        expected_keys=_RECEIPT_KEYS,
        context="validated_request.json",
    )
    _require_exact_evidence_schema(
        provenance,
        expected_keys=_PROVENANCE_KEYS,
        context="provenance.json",
    )
    expected_identities = {
        "input_plan.json": ("input_plan", "input_semantics_only"),
        "validated_request.json": (
            "validated_input_request",
            "input_semantics_only",
        ),
        "provenance.json": ("input_provenance", None),
    }
    for name, (kind, scope) in expected_identities.items():
        document = member_documents[name]
        if document.get("schema_version") != VALIDATION_BUNDLE_SCHEMA_VERSION:
            raise InputIntegrityError(
                "A validation-bundle member has an incompatible schema.",
                details={"path": name},
            )
        if document.get("kind") != kind or document.get("plan_id") != plan_id:
            raise InputIntegrityError(
                "A validation-bundle member has a mismatched identity.",
                details={"path": name},
            )
        if scope is not None and document.get("validation_scope") != scope:
            raise InputIntegrityError(
                "A validation-bundle member has a mismatched validation scope.",
                details={"path": name},
            )

    receipt_data = receipt["data"]
    receipt_artifacts = receipt["artifacts"]
    receipt_warnings = _require_warning_inventory(
        receipt["warnings"], context="validated_request.warnings"
    )
    if not isinstance(receipt_data, Mapping):
        raise InputIntegrityError("The validated receipt data field must be an object.")
    if not isinstance(receipt_artifacts, list):
        raise InputIntegrityError(
            "The validated receipt artifact inventory must be an array."
        )
    recomputed_plan_id = compute_validation_plan_id(
        {
            "data": receipt_data,
            "warnings": receipt_warnings,
            "artifacts": receipt_artifacts,
        }
    )
    if recomputed_plan_id != plan_id:
        raise InputIntegrityError(
            "The declared validation plan identity cannot be derived from its receipt.",
            details={
                "reason": "plan_id_recomputed_mismatch",
                "declared_plan_id": plan_id,
                "recomputed_plan_id": recomputed_plan_id,
            },
        )

    input_data = plan.get("input")
    validation = plan.get("validation")
    if not isinstance(input_data, Mapping) or not isinstance(validation, Mapping):
        raise InputIntegrityError(
            "The validated input plan is structurally incomplete."
        )
    if not canonical_json_equal(receipt_data, input_data):
        raise InputIntegrityError(
            "The validated receipt and input plan contain different normalized inputs.",
            details={"reason": "normalized_input_mismatch"},
        )
    expected_validation = {
        **_EXPECTED_PLAN_VALIDATION_STATE,
        "input_semantics": receipt_data.get("input_certification_status"),
    }
    if not canonical_json_equal(validation, expected_validation):
        raise InputIntegrityError(
            "The input plan validation state has an incompatible schema or value.",
            details={"reason": "plan_validation_state_mismatch"},
        )
    eligibility_checks = {
        "plan.input_semantics": validation.get("input_semantics") == "passed",
        "input.validation_level": input_data.get("validation_level") == "validate",
        "input.certification_scope": input_data.get("certification_scope")
        == "input_semantics_only",
        "input.input_certification_eligible": input_data.get(
            "input_certification_eligible"
        )
        is True,
        "input.input_certification_status": input_data.get("input_certification_status")
        == "passed",
    }
    if not all(eligibility_checks.values()):
        raise InvalidRequestError(
            "Only an input-eligible completed validation bundle can be analyzed.",
            details={"checks": eligibility_checks},
        )
    if input_data.get("input_semantics") == "salmon_quant_dirs_three_prime":
        route = input_data.get("route")
        permitted = (
            isinstance(route, Mapping)
            and route.get("certified_path_execution_permitted") is True
        )
        if not permitted:
            raise InvalidRequestError(
                "The three-prime input route is not permitted for certified execution."
            )
    consumed = receipt_artifacts
    if not isinstance(consumed, list) or not consumed:
        raise InputIntegrityError(
            "The validated receipt lacks its consumed-artifact inventory."
        )
    provenance_warnings = _require_warning_inventory(
        provenance["warnings"], context="provenance.warnings"
    )
    if not canonical_json_equal(provenance_warnings, receipt_warnings):
        raise InputIntegrityError(
            "The validation provenance and receipt warning inventories disagree.",
            details={"reason": "validation_warnings_mismatch"},
        )
    if not canonical_json_equal(provenance["consumed_artifacts"], consumed):
        raise InputIntegrityError(
            "The validation provenance and receipt artifact inventories disagree.",
            details={"reason": "validation_artifacts_mismatch"},
        )
    toolkit = _require_exact_evidence_schema(
        provenance["toolkit"],
        expected_keys={"name", "version"},
        context="provenance.toolkit",
    )
    runtime = _require_exact_evidence_schema(
        provenance["runtime"],
        expected_keys={"python", "implementation", "platform"},
        context="provenance.runtime",
    )
    if (
        toolkit["name"] != "rnaseq-downstream"
        or not isinstance(toolkit["version"], str)
        or not toolkit["version"]
    ):
        raise InputIntegrityError(
            "The validation provenance toolkit identity is incompatible."
        )
    if any(
        not isinstance(runtime[field], str) or not runtime[field] for field in runtime
    ):
        raise InputIntegrityError(
            "The validation provenance runtime identity is incomplete."
        )
    normalized_artifacts: list[dict[str, Any]] = []
    for index, artifact in enumerate(consumed):
        if not isinstance(artifact, Mapping):
            raise InputIntegrityError(
                "A consumed-artifact record is not an object.",
                details={"artifact_index": index},
            )
        required = {"kind", "role", "declared_path", "path", "sha256", "size_bytes"}
        if set(artifact) != required or artifact.get("kind") != "consumed_input":
            raise InputIntegrityError(
                "A consumed-artifact record is incompatible.",
                details={"artifact_index": index},
            )
        for field in ("role", "declared_path", "path"):
            if not isinstance(artifact.get(field), str) or not artifact[field]:
                raise InputIntegrityError(
                    "A consumed-artifact string field is invalid.",
                    details={"artifact_index": index, "field": field},
                )
        _require_sha256(artifact.get("sha256"), field=f"artifacts[{index}].sha256")
        if (
            isinstance(artifact.get("size_bytes"), bool)
            or not isinstance(artifact.get("size_bytes"), int)
            or artifact["size_bytes"] < 0
        ):
            raise InputIntegrityError(
                "A consumed-artifact size is invalid.",
                details={"artifact_index": index},
            )
        normalized_artifacts.append(dict(artifact))
    return (
        plan_id,
        dict(input_data),
        tuple(normalized_artifacts),
        tuple(sorted(bundle_artifacts, key=lambda item: item["role"])),
    )


def _parse_design(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'design' must be an object.")
    require_expected_keys(
        value,
        allowed={"intercept", "terms", "variables"},
        required={"intercept", "terms", "variables"},
        context="analysis design",
    )
    intercept = value["intercept"]
    if not isinstance(intercept, bool):
        raise InvalidRequestError("'design.intercept' must be a boolean.")
    terms = value["terms"]
    variables = value["variables"]
    if not isinstance(terms, list) or not terms:
        raise InvalidRequestError("'design.terms' must be a non-empty array.")
    if not isinstance(variables, Mapping):
        raise InvalidRequestError("'design.variables' must be an object.")
    normalized_terms: list[str] = []
    for index, term in enumerate(terms):
        name = require_nonempty_string(term, field=f"design.terms[{index}]")
        if _IDENTIFIER_PATTERN.fullmatch(name) is None:
            raise InvalidRequestError(
                "Design term names must use a conservative R-safe identifier grammar.",
                details={"term": name},
            )
        if name in _R_RESERVED_WORDS:
            raise InvalidRequestError(
                "Design term names must not be R reserved words.",
                details={"term": name, "reason": "r_reserved_word"},
            )
        normalized_terms.append(name)
    if len(set(normalized_terms)) != len(normalized_terms):
        raise InvalidRequestError("'design.terms' must not contain duplicates.")
    if set(variables) != set(normalized_terms):
        raise InvalidRequestError(
            "'design.variables' must describe exactly the declared design terms.",
            details={
                "missing_variables": sorted(set(normalized_terms) - set(variables)),
                "unexpected_variables": sorted(set(variables) - set(normalized_terms)),
            },
        )
    normalized_variables: dict[str, dict[str, Any]] = {}
    for term in normalized_terms:
        specification = variables[term]
        if not isinstance(specification, Mapping):
            raise InvalidRequestError(
                "Every design variable specification must be an object.",
                details={"term": term},
            )
        variable_type = specification.get("type")
        if variable_type == "continuous":
            require_expected_keys(
                specification,
                allowed={"type"},
                required={"type"},
                context=f"continuous design variable '{term}'",
            )
            normalized_variables[term] = {"type": "continuous"}
        elif variable_type == "factor":
            require_expected_keys(
                specification,
                allowed={"type", "levels"},
                required={"type", "levels"},
                context=f"factor design variable '{term}'",
            )
            levels = specification["levels"]
            if not isinstance(levels, list) or len(levels) < 2:
                raise InvalidRequestError(
                    "A factor design variable requires at least two ordered levels.",
                    details={"term": term},
                )
            normalized_levels = [
                require_nonempty_string(
                    level, field=f"design.variables.{term}.levels[{index}]"
                )
                for index, level in enumerate(levels)
            ]
            if len(set(normalized_levels)) != len(normalized_levels):
                raise InvalidRequestError(
                    "Factor levels must be unique.", details={"term": term}
                )
            normalized_variables[term] = {
                "type": "factor",
                "levels": normalized_levels,
            }
        else:
            raise InvalidRequestError(
                "A design variable type must be 'factor' or 'continuous'.",
                details={"term": term, "observed": variable_type},
            )
    return {
        "intercept": intercept,
        "terms": normalized_terms,
        "variables": normalized_variables,
    }


def _finite_number(value: Any, *, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidRequestError(
            f"'{field}' must be a finite JSON number.", details={"field": field}
        )
    normalized = float(value)
    if not math.isfinite(normalized):
        raise InvalidRequestError(
            f"'{field}' must be finite.", details={"field": field}
        )
    return normalized


def _parse_contrasts(value: Any) -> tuple[dict[str, Any], ...]:
    if not isinstance(value, list) or not value:
        raise InvalidRequestError("'contrasts' must be a non-empty array.")
    normalized: list[dict[str, Any]] = []
    identifiers: set[str] = set()
    for index, contrast in enumerate(value):
        if not isinstance(contrast, Mapping):
            raise InvalidRequestError(
                "Every contrast must be an object.", details={"contrast_index": index}
            )
        require_expected_keys(
            contrast,
            allowed={"contrast_id", "weights", "lfc_threshold"},
            required={"contrast_id", "weights", "lfc_threshold"},
            context=f"contrast {index}",
        )
        identifier = require_nonempty_string(
            contrast["contrast_id"], field=f"contrasts[{index}].contrast_id"
        )
        if _CONTRAST_ID_PATTERN.fullmatch(identifier) is None:
            raise InvalidRequestError(
                "A contrast identifier has unsupported characters.",
                details={"contrast_id": identifier},
            )
        if identifier in identifiers:
            raise InvalidRequestError(
                "Contrast identifiers must be unique.",
                details={"contrast_id": identifier},
            )
        identifiers.add(identifier)
        weights = contrast["weights"]
        if not isinstance(weights, Mapping) or not weights:
            raise ContrastNotEstimableError(
                "Contrast weights must be a non-empty object.",
                details={"contrast_id": identifier},
            )
        normalized_weights: dict[str, float] = {}
        for coefficient, weight in weights.items():
            coefficient_name = require_nonempty_string(
                coefficient,
                field=f"contrasts[{index}].weights coefficient",
            )
            normalized_weights[coefficient_name] = _finite_number(
                weight,
                field=f"contrasts[{index}].weights.{coefficient_name}",
            )
        if all(weight == 0 for weight in normalized_weights.values()):
            raise ContrastNotEstimableError(
                "A contrast must contain at least one non-zero weight.",
                details={"contrast_id": identifier, "reason": "contrast_zero"},
            )
        threshold = _finite_number(
            contrast["lfc_threshold"],
            field=f"contrasts[{index}].lfc_threshold",
        )
        if threshold < 0:
            raise InvalidRequestError(
                "A log2-fold-change testing threshold cannot be negative.",
                details={"contrast_id": identifier, "lfc_threshold": threshold},
            )
        normalized.append(
            {
                "contrast_id": identifier,
                "weights": normalized_weights,
                "lfc_threshold": threshold,
            }
        )
    return tuple(normalized)


def _parse_display(value: Any) -> dict[str, Any]:
    """Parse the version-1.1 display-only request.

    Display settings are deliberately kept out of the version-1.0 backend
    protocol.  They may control only post-result rendering and exploratory PCA
    display; they must not change input construction, model fitting, contrasts,
    P values, or multiple-testing adjustment.
    """

    if not isinstance(value, Mapping):
        raise InvalidRequestError("'display' must be an object.")
    require_expected_keys(
        value,
        allowed={"fdr_threshold", "pca_top_n", "pca_components"},
        required={"fdr_threshold", "pca_top_n", "pca_components"},
        context="analysis display request",
    )

    fdr_threshold = _finite_number(
        value["fdr_threshold"], field="display.fdr_threshold"
    )
    if fdr_threshold < 0 or fdr_threshold > 1:
        raise InvalidRequestError(
            "'display.fdr_threshold' must lie in the closed interval [0, 1].",
            details={"fdr_threshold": fdr_threshold},
        )

    pca_top_n = value["pca_top_n"]
    if isinstance(pca_top_n, bool) or not isinstance(pca_top_n, int) or pca_top_n <= 0:
        raise InvalidRequestError("'display.pca_top_n' must be a positive integer.")

    pca_components = value["pca_components"]
    if not isinstance(pca_components, list) or len(pca_components) != 2:
        raise InvalidRequestError(
            "'display.pca_components' must contain exactly two component numbers."
        )
    normalized_components: list[int] = []
    for index, component in enumerate(pca_components):
        if (
            isinstance(component, bool)
            or not isinstance(component, int)
            or component <= 0
        ):
            raise InvalidRequestError(
                "PCA component numbers must be positive integers.",
                details={
                    "field": f"display.pca_components[{index}]",
                    "observed": component,
                },
            )
        normalized_components.append(component)
    if len(set(normalized_components)) != len(normalized_components):
        raise InvalidRequestError(
            "'display.pca_components' must contain two distinct component numbers."
        )
    if max(normalized_components) > pca_top_n:
        raise InvalidRequestError(
            "A requested PCA component cannot exceed 'display.pca_top_n'.",
            details={
                "pca_top_n": pca_top_n,
                "pca_components": normalized_components,
            },
        )

    return {
        "fdr_threshold": fdr_threshold,
        "pca_top_n": pca_top_n,
        "pca_components": normalized_components,
    }


def load_analysis_request(path: str | Path) -> ValidatedAnalysisRequest:
    """Validate an analysis request and its complete checkpoint-A bundle."""

    request_path = _resolve_analysis_request(path)
    request_content, request_digest, _ = _capture(request_path)
    document = read_json_object(
        request_path,
        document_role="analysis_request",
        content=request_content,
    )
    base_keys = {
        "schema_version",
        "validated_input_bundle",
        "design",
        "contrasts",
    }
    if "schema_version" not in document:
        require_expected_keys(
            document,
            allowed=base_keys,
            required=base_keys,
            context="analysis request",
        )
    request_schema_version = document["schema_version"]
    if request_schema_version not in ANALYSIS_REQUEST_SCHEMA_VERSIONS:
        raise InvalidRequestError(
            "The analysis request schema version is unsupported.",
            details={
                "observed": request_schema_version,
                "supported": list(ANALYSIS_REQUEST_SCHEMA_VERSIONS),
            },
        )
    if request_schema_version == "1.0":
        require_expected_keys(
            document,
            allowed=base_keys,
            required=base_keys,
            context="analysis request",
        )
        display = None
        gene_sets = None
    else:
        require_expected_keys(
            document,
            allowed=base_keys | {"display", "gene_sets"},
            required=base_keys | {"display"},
            context="analysis request version 1.1",
        )
        display = _parse_display(document["display"])
        gene_sets = (
            parse_gene_sets_request(document["gene_sets"], request_path=request_path)
            if "gene_sets" in document
            else None
        )
    bundle_path = _resolve_bundle(
        document["validated_input_bundle"], request_path=request_path
    )
    (
        plan_id,
        input_data,
        consumed_artifacts,
        validation_bundle_artifacts,
    ) = _read_verified_bundle(bundle_path)
    return ValidatedAnalysisRequest(
        request_path=request_path,
        request_sha256=request_digest,
        request_schema_version=request_schema_version,
        bundle_path=bundle_path,
        plan_id=plan_id,
        input_data=input_data,
        consumed_artifacts=consumed_artifacts,
        validation_bundle_artifacts=validation_bundle_artifacts,
        design=_parse_design(document["design"]),
        contrasts=_parse_contrasts(document["contrasts"]),
        display=display,
        gene_sets=gene_sets,
    )


__all__ = [
    "ANALYSIS_REQUEST_SCHEMA_VERSIONS",
    "ANALYSIS_SCHEMA_VERSION",
    "ValidatedAnalysisRequest",
    "load_analysis_request",
]
