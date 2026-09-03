from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from rnaseq_downstream.analysis_contract import (
    _parse_design,
    _parse_display,
    load_analysis_request,
)
from rnaseq_downstream.edger_backend import (
    _capture_output_json,
    _invoke_r,
    _raise_backend_response,
)
from rnaseq_downstream.errors import (
    BackendFailedError,
    ContrastNotEstimableError,
    DesignRankDeficientError,
    InputIntegrityError,
    InvalidRequestError,
)
from rnaseq_downstream.validation_bundle import validate_request_to_bundle


def _write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def _write_json(path: Path, document: object) -> Path:
    return _write(path, json.dumps(document, sort_keys=True) + "\n")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rewrite_bundle(
    bundle: Path,
    documents: dict[str, dict[str, object]],
    *,
    plan_id: str,
) -> None:
    """Rewrite declared evidence and rebuild every manifest member digest."""

    members: list[dict[str, object]] = []
    for name in ("input_plan.json", "provenance.json", "validated_request.json"):
        document = documents[name]
        document["plan_id"] = plan_id
        _write_json(bundle / name, document)
        members.append(
            {
                "path": name,
                "sha256": _sha256(bundle / name),
                "size_bytes": (bundle / name).stat().st_size,
            }
        )
    _write_json(
        bundle / "bundle_manifest.json",
        {
            "schema_version": "1.0",
            "kind": "validation_bundle_manifest",
            "plan_id": plan_id,
            "members": members,
        },
    )


def _read_bundle_documents(bundle: Path) -> dict[str, dict[str, object]]:
    return {
        name: json.loads((bundle / name).read_text(encoding="utf-8"))
        for name in ("input_plan.json", "provenance.json", "validated_request.json")
    }


def _validated_bundle(root: Path) -> Path:
    reference = _write(root / "reference.fa", ">synthetic\nACGT\n")
    matrix = _write(
        root / "counts.tsv",
        "gene_id\ts1\ts2\nENSG1\t10\t20\nENSG2\t30\t40\n",
    )
    _write(root / "metadata.tsv", "sample_id\tcondition\ns1\tcontrol\ns2\ttreated\n")
    _write_json(
        root / "counts.manifest.json",
        {
            "schema_version": "1.0",
            "artifact_type": "featurecounts_integer_matrix",
            "artifact": {"path": "counts.tsv", "sha256": _sha256(matrix)},
            "gene_id_column": "gene_id",
            "display_columns": [],
            "sample_columns": ["s1", "s2"],
            "producer": {"name": "featureCounts", "version": "2.0.6"},
            "reference": {
                "name": "synthetic",
                "version": "1",
                "source": "reference.fa",
                "sha256": _sha256(reference),
            },
        },
    )
    request = _write_json(
        root / "request.json",
        {
            "schema_version": "1.0",
            "input_semantics": "featurecounts_integer",
            "metadata": {"path": "metadata.tsv", "sample_id_column": "sample_id"},
            "producer": {"name": "featureCounts", "version": "2.0.6"},
            "reference": {
                "name": "synthetic",
                "version": "1",
                "source": "reference.fa",
            },
            "gene_id": {"strip_version": False},
            "samples": [{"sample_id": "s1"}, {"sample_id": "s2"}],
            "featurecounts": {
                "layout": "combined_matrix",
                "matrix": "counts.tsv",
                "manifest": "counts.manifest.json",
            },
        },
    )
    bundle = root / "validated"
    validate_request_to_bundle(request, bundle)
    return bundle


def _analysis_request(
    root: Path,
    bundle: Path,
    *,
    weights: dict[str, float],
    schema_version: str = "1.0",
    display: object | None = None,
) -> Path:
    document: dict[str, object] = {
        "schema_version": schema_version,
        "validated_input_bundle": str(bundle),
        "design": {
            "intercept": True,
            "terms": ["condition"],
            "variables": {
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                }
            },
        },
        "contrasts": [
            {
                "contrast_id": "treated_vs_control",
                "weights": weights,
                "lfc_threshold": 0,
            }
        ],
    }
    if display is not None:
        document["display"] = display
    return _write_json(root / "analysis.json", document)


def test_analysis_contract_accepts_only_complete_eligible_bundle(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 1})

    validated = load_analysis_request(request)

    assert validated.request_schema_version == "1.0"
    assert validated.display is None
    assert validated.input_data["input_certification_eligible"] is True
    assert validated.input_data["input_semantics"] == "featurecounts_integer"
    assert {item["role"] for item in validated.validation_bundle_artifacts} == {
        "bundle_manifest",
        "input_plan",
        "provenance",
        "validated_request",
    }


def test_analysis_request_v11_accepts_explicit_display_without_changing_backend_protocol(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.1",
        display={
            "fdr_threshold": 0.05,
            "pca_top_n": 500,
            "pca_components": [1, 2],
        },
    )

    validated = load_analysis_request(request)

    assert validated.request_schema_version == "1.1"
    assert validated.display == {
        "fdr_threshold": 0.05,
        "pca_top_n": 500,
        "pca_components": [1, 2],
    }
    backend_document = validated.to_backend_document()
    assert backend_document["schema_version"] == "1.0"
    assert "display" not in backend_document


def test_analysis_request_v10_stays_strict_and_rejects_display(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        display={
            "fdr_threshold": 0.05,
            "pca_top_n": 500,
            "pca_components": [1, 2],
        },
    )

    with pytest.raises(InvalidRequestError) as caught:
        load_analysis_request(request)

    assert caught.value.details["context"] == "analysis request"
    assert caught.value.details["unknown_keys"] == ["display"]


def test_analysis_request_v11_requires_display(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.1",
    )

    with pytest.raises(InvalidRequestError) as caught:
        load_analysis_request(request)

    assert caught.value.details["context"] == "analysis request version 1.1"
    assert caught.value.details["missing_keys"] == ["display"]


def test_analysis_request_reports_both_supported_public_versions(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
    )

    with pytest.raises(InvalidRequestError) as caught:
        load_analysis_request(request)

    assert caught.value.details == {
        "observed": "1.2",
        "supported": ["1.0", "1.1"],
    }


@pytest.mark.parametrize("fdr_threshold", [0, 1])
def test_display_request_accepts_probability_boundaries_and_axis_order(
    fdr_threshold: int,
) -> None:
    assert _parse_display(
        {
            "fdr_threshold": fdr_threshold,
            "pca_top_n": 2,
            "pca_components": [2, 1],
        }
    ) == {
        "fdr_threshold": float(fdr_threshold),
        "pca_top_n": 2,
        "pca_components": [2, 1],
    }


@pytest.mark.parametrize(
    "display",
    [
        None,
        {},
        {
            "fdr_threshold": 0.05,
            "pca_top_n": 500,
            "pca_components": [1, 2],
            "unexpected": True,
        },
        {
            "fdr_threshold": True,
            "pca_top_n": 500,
            "pca_components": [1, 2],
        },
        {
            "fdr_threshold": -0.01,
            "pca_top_n": 500,
            "pca_components": [1, 2],
        },
        {
            "fdr_threshold": 1.01,
            "pca_top_n": 500,
            "pca_components": [1, 2],
        },
        {
            "fdr_threshold": 0.05,
            "pca_top_n": False,
            "pca_components": [1, 2],
        },
        {
            "fdr_threshold": 0.05,
            "pca_top_n": 0,
            "pca_components": [1, 2],
        },
        {
            "fdr_threshold": 0.05,
            "pca_top_n": 500,
            "pca_components": [1],
        },
        {
            "fdr_threshold": 0.05,
            "pca_top_n": 500,
            "pca_components": [1, 1],
        },
        {
            "fdr_threshold": 0.05,
            "pca_top_n": 500,
            "pca_components": [True, 2],
        },
        {
            "fdr_threshold": 0.05,
            "pca_top_n": 2,
            "pca_components": [1, 3],
        },
    ],
)
def test_display_request_rejects_noncanonical_or_unsafe_values(
    display: object,
) -> None:
    with pytest.raises(InvalidRequestError):
        _parse_display(display)


def test_analysis_contract_rejects_tampered_or_extra_bundle_entries(
    tmp_path: Path,
) -> None:
    tampered_bundle = _validated_bundle(tmp_path / "tampered")
    plan = tampered_bundle / "input_plan.json"
    plan.write_bytes(plan.read_bytes() + b" ")
    request = _analysis_request(
        tmp_path / "tampered", tampered_bundle, weights={"conditiontreated": 1}
    )
    with pytest.raises(InputIntegrityError, match="does not match its manifest"):
        load_analysis_request(request)

    extra_bundle = _validated_bundle(tmp_path / "extra")
    _write(extra_bundle / "unmanifested.txt", "unexpected\n")
    request = _analysis_request(
        tmp_path / "extra", extra_bundle, weights={"conditiontreated": 1}
    )
    with pytest.raises(InputIntegrityError, match="unmanifested or unsafe"):
        load_analysis_request(request)


def test_bundle_identity_is_recomputed_after_all_declared_ids_are_rebound(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    documents = _read_bundle_documents(bundle)
    plan_input = documents["input_plan.json"]["input"]
    assert isinstance(plan_input, dict)
    route = plan_input["route"]
    assert isinstance(route, dict)
    route["transcript_length_offset"] = 0

    replacement_id = "f" * 64
    _rewrite_bundle(bundle, documents, plan_id=replacement_id)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 1})

    with pytest.raises(InputIntegrityError) as caught:
        load_analysis_request(request)

    assert caught.value.details["reason"] == "plan_id_recomputed_mismatch"
    assert caught.value.details["declared_plan_id"] == replacement_id
    assert caught.value.details["recomputed_plan_id"] != replacement_id


def test_bundle_input_comparison_does_not_treat_false_as_integer_zero(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    documents = _read_bundle_documents(bundle)
    original_id = documents["validated_request.json"]["plan_id"]
    assert isinstance(original_id, str)
    plan_input = documents["input_plan.json"]["input"]
    assert isinstance(plan_input, dict)
    route = plan_input["route"]
    assert isinstance(route, dict)
    route["transcript_length_offset"] = 0
    assert route["transcript_length_offset"] == bool(0)

    _rewrite_bundle(bundle, documents, plan_id=original_id)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 1})

    with pytest.raises(InputIntegrityError) as caught:
        load_analysis_request(request)

    assert caught.value.details["reason"] == "normalized_input_mismatch"


def test_bundle_receipt_and_provenance_warnings_must_match_exactly(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    documents = _read_bundle_documents(bundle)
    plan_id = documents["validated_request.json"]["plan_id"]
    assert isinstance(plan_id, str)
    documents["provenance.json"]["warnings"] = [
        {
            "code": "ADVERSARIAL_WARNING",
            "severity": "high",
            "message": "This warning was not present in the validated receipt.",
            "details": {},
        }
    ]
    _rewrite_bundle(bundle, documents, plan_id=plan_id)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 1})

    with pytest.raises(InputIntegrityError) as caught:
        load_analysis_request(request)

    assert caught.value.details["reason"] == "validation_warnings_mismatch"


@pytest.mark.parametrize(
    "member_name",
    ["input_plan.json", "provenance.json", "validated_request.json"],
)
def test_bundle_members_require_exact_top_level_schemas(
    tmp_path: Path,
    member_name: str,
) -> None:
    bundle = _validated_bundle(tmp_path)
    documents = _read_bundle_documents(bundle)
    plan_id = documents["validated_request.json"]["plan_id"]
    assert isinstance(plan_id, str)
    documents[member_name]["unexpected"] = "not allowed"
    _rewrite_bundle(bundle, documents, plan_id=plan_id)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 1})

    with pytest.raises(InputIntegrityError, match="incompatible schema"):
        load_analysis_request(request)


def test_zero_contrast_is_a_scientific_estimability_error(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 0})

    with pytest.raises(ContrastNotEstimableError) as caught:
        load_analysis_request(request)

    assert caught.value.to_dict()["code"] == "CONTRAST_NOT_ESTIMABLE"
    assert caught.value.details["reason"] == "contrast_zero"


@pytest.mark.parametrize(
    "reserved_word",
    [
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
    ],
)
def test_design_terms_reject_every_r_reserved_word(reserved_word: str) -> None:
    with pytest.raises(InvalidRequestError) as caught:
        _parse_design(
            {
                "intercept": True,
                "terms": [reserved_word],
                "variables": {reserved_word: {"type": "continuous"}},
            }
        )

    assert caught.value.details == {
        "term": reserved_word,
        "reason": "r_reserved_word",
    }


@pytest.mark.parametrize(
    ("filename", "role"),
    [
        ("backend_manifest.json", "backend_manifest"),
        ("analysis.json", "edger_analysis"),
    ],
)
@pytest.mark.parametrize(
    "payload",
    [
        b'{"kind":',
        b'{"kind":"first","kind":"duplicate"}\n',
        b"[]\n",
    ],
    ids=["malformed", "duplicate_key", "wrong_root"],
)
def test_backend_output_json_parse_failures_are_backend_errors(
    tmp_path: Path,
    filename: str,
    role: str,
    payload: bytes,
) -> None:
    path = tmp_path / filename
    path.write_bytes(payload)

    with pytest.raises(BackendFailedError) as caught:
        _capture_output_json(path, role=role)

    assert caught.value.to_dict()["code"] == "BACKEND_FAILED"
    assert caught.value.details["role"] == role
    assert caught.value.details["parse_error_code"] == "INVALID_REQUEST"


@pytest.mark.parametrize(
    ("code", "exception_type"),
    [
        ("DESIGN_RANK_DEFICIENT", DesignRankDeficientError),
        ("CONTRAST_NOT_ESTIMABLE", ContrastNotEstimableError),
        ("BACKEND_FAILED", BackendFailedError),
    ],
)
def test_r_error_codes_map_to_stable_python_exceptions(
    code: str, exception_type: type[Exception]
) -> None:
    response = {
        "status": "error",
        "errors": [
            {
                "code": code,
                "message": "structured failure",
                "details": {"reason": "test_reason"},
            }
        ],
    }

    with pytest.raises(exception_type) as caught:
        _raise_backend_response(response, returncode=3, stderr="diagnostic")

    assert getattr(caught.value, "details")["reason"] == "test_reason"


def test_missing_r_library_is_a_structured_request_error(tmp_path: Path) -> None:
    with pytest.raises(InvalidRequestError) as caught:
        _invoke_r(
            tmp_path / "backend-request.json",
            tmp_path / "result-stage",
            rscript="Rscript",
            r_library=tmp_path / "missing-library",
        )

    assert caught.value.to_dict()["code"] == "INVALID_REQUEST"
    assert caught.value.details["cause_type"] == "FileNotFoundError"
