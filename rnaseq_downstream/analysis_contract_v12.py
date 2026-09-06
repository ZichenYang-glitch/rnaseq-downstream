"""Version-1.2 analysis-request routing layered over the frozen P0 contract.

The P0/C1/C2 benchmark archives bind ``analysis_contract.py`` by digest.  This
module therefore extends that contract without changing the frozen module or
its version-1.0/version-1.1 private edgeR serialization.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from . import analysis_contract as _legacy
from .errors import InvalidRequestError
from .gene_sets import parse_gene_sets_request
from .provenance import read_json_object, require_expected_keys, require_nonempty_string


DESEQ2_ANALYSIS_SCHEMA_VERSION = "1.0"
ANALYSIS_REQUEST_SCHEMA_VERSIONS = (*_legacy.ANALYSIS_REQUEST_SCHEMA_VERSIONS, "1.2")
ANALYSIS_BACKENDS = ("edger_ql", "deseq2")


@dataclass(frozen=True)
class ValidatedAnalysisRequest(_legacy.ValidatedAnalysisRequest):
    """A normalized multi-backend request plus verified checkpoint-A evidence."""

    backend: str
    deseq2: dict[str, Any] | None

    def to_backend_document(self) -> dict[str, Any]:
        """Return the frozen edgeR private document without new public fields."""

        if self.backend != "edger_ql":
            raise InvalidRequestError(
                "An edgeR backend document cannot be created for a non-edgeR request.",
                details={"backend": self.backend},
            )
        return super().to_backend_document()

    def to_deseq2_backend_document(self) -> dict[str, Any]:
        """Return the independent normalized DESeq2 private request document."""

        if self.backend != "deseq2" or self.deseq2 is None:
            raise InvalidRequestError(
                "A normalized DESeq2 request is required for DESeq2 serialization.",
                details={"backend": self.backend},
            )
        return {
            "schema_version": DESEQ2_ANALYSIS_SCHEMA_VERSION,
            "kind": "deseq2_backend_request",
            "analysis_request": {
                "path": str(self.request_path),
                "sha256": self.request_sha256,
                "schema_version": "1.2",
                "backend": "deseq2",
            },
            "validated_input_bundle": {
                "path": str(self.bundle_path),
                "plan_id": self.plan_id,
                "artifacts": [dict(item) for item in self.validation_bundle_artifacts],
            },
            "input": dict(self.input_data),
            "design": dict(self.design),
            "contrasts": [dict(item) for item in self.contrasts],
            "deseq2": copy.deepcopy(self.deseq2),
        }


def _parse_reduced_design(
    value: Any,
    *,
    full_design: Mapping[str, Any],
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'deseq2.reduced' must be an object.")
    require_expected_keys(
        value,
        allowed={"intercept", "terms"},
        required={"intercept", "terms"},
        context="DESeq2 reduced design",
    )
    intercept = value["intercept"]
    if not isinstance(intercept, bool):
        raise InvalidRequestError("'deseq2.reduced.intercept' must be a boolean.")
    if intercept != full_design["intercept"]:
        raise InvalidRequestError(
            "The DESeq2 reduced design must keep the full design intercept policy.",
            details={
                "full_intercept": full_design["intercept"],
                "reduced_intercept": intercept,
            },
        )
    terms = value["terms"]
    if not isinstance(terms, list):
        raise InvalidRequestError("'deseq2.reduced.terms' must be an array.")
    normalized_terms = [
        require_nonempty_string(term, field=f"deseq2.reduced.terms[{index}]")
        for index, term in enumerate(terms)
    ]
    if len(set(normalized_terms)) != len(normalized_terms):
        raise InvalidRequestError("'deseq2.reduced.terms' must not contain duplicates.")
    full_terms = list(full_design["terms"])
    unknown = sorted(set(normalized_terms) - set(full_terms))
    if unknown:
        raise InvalidRequestError(
            "The DESeq2 reduced design contains terms absent from the full design.",
            details={"unknown_terms": unknown},
        )
    if set(normalized_terms) == set(full_terms):
        raise InvalidRequestError(
            "The DESeq2 reduced design must be a proper additive-term subset.",
            details={"full_terms": full_terms, "reduced_terms": normalized_terms},
        )
    if not intercept and not normalized_terms:
        raise InvalidRequestError(
            "A no-intercept DESeq2 reduced design must retain at least one term."
        )
    selected = set(normalized_terms)
    return {
        "intercept": intercept,
        "terms": [term for term in full_terms if term in selected],
    }


def _parse_deseq2_configuration(
    value: Any,
    *,
    design: Mapping[str, Any],
    contrasts: tuple[dict[str, Any], ...],
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise InvalidRequestError("'deseq2' must be an object.")
    test_mode = value.get("test_mode")
    if test_mode == "wald":
        require_expected_keys(
            value,
            allowed={"test_mode", "shrinkage"},
            required={"test_mode", "shrinkage"},
            context="DESeq2 Wald configuration",
        )
    elif test_mode == "lrt":
        require_expected_keys(
            value,
            allowed={"test_mode", "shrinkage", "reduced"},
            required={"test_mode", "shrinkage", "reduced"},
            context="DESeq2 LRT configuration",
        )
    else:
        raise InvalidRequestError(
            "'deseq2.test_mode' must be exactly 'wald' or 'lrt'.",
            details={"observed": test_mode},
        )
    shrinkage = value["shrinkage"]
    if not isinstance(shrinkage, str) or shrinkage not in {"none", "apeglm"}:
        raise InvalidRequestError(
            "'deseq2.shrinkage' must be exactly 'none' or 'apeglm'.",
            details={"observed": shrinkage},
        )

    normalized: dict[str, Any] = {
        "test_mode": test_mode,
        "shrinkage": shrinkage,
    }
    if test_mode == "lrt":
        if len(contrasts) != 1:
            raise InvalidRequestError(
                "A DESeq2 LRT request must contain exactly one reporting contrast.",
                details={"contrast_count": len(contrasts)},
            )
        if contrasts[0]["lfc_threshold"] != 0:
            raise InvalidRequestError(
                "A DESeq2 LRT reporting contrast cannot use an LFC threshold.",
                details={
                    "contrast_id": contrasts[0]["contrast_id"],
                    "lfc_threshold": contrasts[0]["lfc_threshold"],
                },
            )
        normalized["reduced"] = _parse_reduced_design(
            value["reduced"], full_design=design
        )

    if shrinkage == "apeglm":
        for contrast in contrasts:
            weights = contrast["weights"]
            valid = (
                len(weights) == 1
                and next(iter(weights.values())) == 1.0
                and next(iter(weights)) not in {"(Intercept)", "Intercept"}
            )
            if not valid:
                raise InvalidRequestError(
                    "apeglm shrinkage requires every contrast to select exactly "
                    "one non-intercept coefficient with weight +1.",
                    details={
                        "contrast_id": contrast["contrast_id"],
                        "weights": dict(weights),
                    },
                )
    return normalized


def load_analysis_request(path: str | Path) -> ValidatedAnalysisRequest:
    """Validate public request versions 1.0-1.2 against one frozen input bundle."""

    request_path = _legacy._resolve_analysis_request(path)
    request_content, request_digest, _ = _legacy._capture(request_path)
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

    if request_schema_version in _legacy.ANALYSIS_REQUEST_SCHEMA_VERSIONS:
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
            display = _legacy._parse_display(document["display"])
            gene_sets = (
                parse_gene_sets_request(
                    document["gene_sets"], request_path=request_path
                )
                if "gene_sets" in document
                else None
            )
        bundle_path = _legacy._resolve_bundle(
            document["validated_input_bundle"], request_path=request_path
        )
        (
            plan_id,
            input_data,
            consumed_artifacts,
            validation_bundle_artifacts,
        ) = _legacy._read_verified_bundle(bundle_path)
        return ValidatedAnalysisRequest(
            request_path=request_path,
            request_sha256=request_digest,
            request_schema_version=request_schema_version,
            bundle_path=bundle_path,
            plan_id=plan_id,
            input_data=input_data,
            consumed_artifacts=consumed_artifacts,
            validation_bundle_artifacts=validation_bundle_artifacts,
            design=_legacy._parse_design(document["design"]),
            contrasts=_legacy._parse_contrasts(document["contrasts"]),
            display=display,
            gene_sets=gene_sets,
            backend="edger_ql",
            deseq2=None,
        )

    require_expected_keys(
        document,
        allowed=base_keys | {"backend", "display", "gene_sets", "deseq2"},
        required=base_keys,
        context="analysis request version 1.2",
    )
    backend = document.get("backend", "edger_ql")
    if backend not in ANALYSIS_BACKENDS:
        raise InvalidRequestError(
            "'backend' must be exactly 'edger_ql' or 'deseq2'.",
            details={"observed": backend},
        )

    if backend == "edger_ql":
        if "deseq2" in document:
            raise InvalidRequestError(
                "An edgeR request must not contain a DESeq2 configuration.",
                details={"backend": backend, "incompatible_fields": ["deseq2"]},
            )
        if "gene_sets" in document and "display" not in document:
            raise InvalidRequestError(
                "An edgeR gene-set request requires the atomic display extension.",
                details={"missing_fields": ["display"]},
            )
        display = (
            _legacy._parse_display(document["display"])
            if "display" in document
            else None
        )
        gene_sets = (
            parse_gene_sets_request(document["gene_sets"], request_path=request_path)
            if "gene_sets" in document
            else None
        )
        deseq2 = None
    else:
        incompatible_fields = [
            field for field in ("display", "gene_sets") if field in document
        ]
        if incompatible_fields:
            raise InvalidRequestError(
                "The D1 DESeq2 path does not support display or gene-set analysis.",
                details={
                    "backend": backend,
                    "incompatible_fields": incompatible_fields,
                },
            )
        if "deseq2" not in document:
            raise InvalidRequestError(
                "A DESeq2 request requires the 'deseq2' configuration object.",
                details={"missing_keys": ["deseq2"]},
            )
        display = None
        gene_sets = None
        deseq2 = None

    design = _legacy._parse_design(document["design"])
    contrasts = _legacy._parse_contrasts(document["contrasts"])
    if backend == "deseq2":
        deseq2 = _parse_deseq2_configuration(
            document["deseq2"], design=design, contrasts=contrasts
        )

    bundle_path = _legacy._resolve_bundle(
        document["validated_input_bundle"], request_path=request_path
    )
    (
        plan_id,
        input_data,
        consumed_artifacts,
        validation_bundle_artifacts,
    ) = _legacy._read_verified_bundle(bundle_path)
    return ValidatedAnalysisRequest(
        request_path=request_path,
        request_sha256=request_digest,
        request_schema_version=request_schema_version,
        bundle_path=bundle_path,
        plan_id=plan_id,
        input_data=input_data,
        consumed_artifacts=consumed_artifacts,
        validation_bundle_artifacts=validation_bundle_artifacts,
        design=design,
        contrasts=contrasts,
        display=display,
        gene_sets=gene_sets,
        backend=backend,
        deseq2=deseq2,
    )


__all__ = [
    "ANALYSIS_BACKENDS",
    "ANALYSIS_REQUEST_SCHEMA_VERSIONS",
    "DESEQ2_ANALYSIS_SCHEMA_VERSION",
    "ValidatedAnalysisRequest",
    "load_analysis_request",
]
