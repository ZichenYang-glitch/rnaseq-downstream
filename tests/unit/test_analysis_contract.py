from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from rnaseq_downstream.analysis_contract import load_analysis_request
from rnaseq_downstream.edger_backend import _invoke_r, _raise_backend_response
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


def _analysis_request(root: Path, bundle: Path, *, weights: dict[str, float]) -> Path:
    return _write_json(
        root / "analysis.json",
        {
            "schema_version": "1.0",
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
        },
    )


def test_analysis_contract_accepts_only_complete_eligible_bundle(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 1})

    validated = load_analysis_request(request)

    assert validated.input_data["input_certification_eligible"] is True
    assert validated.input_data["input_semantics"] == "featurecounts_integer"
    assert {item["role"] for item in validated.validation_bundle_artifacts} == {
        "bundle_manifest",
        "input_plan",
        "provenance",
        "validated_request",
    }


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


def test_zero_contrast_is_a_scientific_estimability_error(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 0})

    with pytest.raises(ContrastNotEstimableError) as caught:
        load_analysis_request(request)

    assert caught.value.to_dict()["code"] == "CONTRAST_NOT_ESTIMABLE"
    assert caught.value.details["reason"] == "contrast_zero"


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
