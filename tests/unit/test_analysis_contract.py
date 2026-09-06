from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from rnaseq_downstream.analysis_contract import (
    _parse_design,
    _parse_display,
    load_analysis_request as load_legacy_analysis_request,
)
from rnaseq_downstream.analysis_contract_v12 import load_analysis_request
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


_MISSING = object()


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
    gene_sets: object | None = None,
    backend: object = _MISSING,
    deseq2: object = _MISSING,
    contrasts: list[dict[str, object]] | None = None,
    design: dict[str, object] | None = None,
) -> Path:
    document: dict[str, object] = {
        "schema_version": schema_version,
        "validated_input_bundle": str(bundle),
        "design": design
        or {
            "intercept": True,
            "terms": ["condition"],
            "variables": {
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                }
            },
        },
        "contrasts": contrasts
        or [
            {
                "contrast_id": "treated_vs_control",
                "weights": weights,
                "lfc_threshold": 0,
            }
        ],
    }
    if display is not None:
        document["display"] = display
    if gene_sets is not None:
        document["gene_sets"] = gene_sets
    if backend is not _MISSING:
        document["backend"] = backend
    if deseq2 is not _MISSING:
        document["deseq2"] = deseq2
    return _write_json(root / "analysis.json", document)


def test_analysis_contract_accepts_only_complete_eligible_bundle(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(tmp_path, bundle, weights={"conditiontreated": 1})

    validated = load_analysis_request(request)

    assert validated.request_schema_version == "1.0"
    assert validated.backend == "edger_ql"
    assert validated.deseq2 is None
    assert validated.display is None
    assert validated.input_data["input_certification_eligible"] is True
    assert validated.input_data["input_semantics"] == "featurecounts_integer"
    assert {item["role"] for item in validated.validation_bundle_artifacts} == {
        "bundle_manifest",
        "input_plan",
        "provenance",
        "validated_request",
    }
    backend_document = validated.to_backend_document()
    assert set(backend_document) == {
        "schema_version",
        "kind",
        "analysis_request",
        "validated_input_bundle",
        "input",
        "design",
        "contrasts",
    }
    assert backend_document["analysis_request"] == {
        "path": str(request.resolve()),
        "sha256": _sha256(request),
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
    assert validated.gene_sets is None
    backend_document = validated.to_backend_document()
    assert backend_document["schema_version"] == "1.0"
    assert backend_document["analysis_request"] == {
        "path": str(request.resolve()),
        "sha256": _sha256(request),
    }
    assert "display" not in backend_document


@pytest.mark.parametrize("schema_version", ["1.0", "1.1"])
def test_v12_wrapper_preserves_legacy_private_edge_document_exactly(
    tmp_path: Path,
    schema_version: str,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version=schema_version,
        display=(
            {
                "fdr_threshold": 0.05,
                "pca_top_n": 500,
                "pca_components": [1, 2],
            }
            if schema_version == "1.1"
            else None
        ),
    )

    legacy = load_legacy_analysis_request(request)
    extended = load_analysis_request(request)

    assert extended.backend == "edger_ql"
    assert extended.deseq2 is None
    assert extended.to_backend_document() == legacy.to_backend_document()


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


@pytest.mark.parametrize("schema_version", ["1.1", "1.2"])
def test_edger_analysis_request_accepts_optional_frozen_gene_set_sources(
    tmp_path: Path,
    schema_version: str,
) -> None:
    bundle = _validated_bundle(tmp_path)
    gmt = _write(tmp_path / "sets.gmt", "SET_A\tdescription\tA\tB\n")
    annotation = _write(
        tmp_path / "annotation.tsv",
        "gene_id\tsymbol\nENSG1\tA\nENSG2\tB\n",
    )
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version=schema_version,
        display={
            "fdr_threshold": 0.05,
            "pca_top_n": 500,
            "pca_components": [1, 2],
        },
        gene_sets={
            "gmt": {
                "path": gmt.name,
                "sha256": _sha256(gmt),
                "collection": "synthetic_hallmarks",
                "version": "2026.1",
                "identifier_type": "symbol",
            },
            "annotation": {
                "path": annotation.name,
                "sha256": _sha256(annotation),
                "name": "synthetic_ensembl",
                "version": "1",
                "gene_id_column": "gene_id",
                "symbol_column": "symbol",
            },
            "minimum_tested_genes": 2,
        },
    )

    validated = load_analysis_request(request)

    assert validated.backend == "edger_ql"
    assert validated.gene_sets is not None
    assert validated.gene_sets["gmt"]["path"] == str(gmt.resolve())
    assert validated.gene_sets["gmt"]["declared_path"] == gmt.name
    assert validated.gene_sets["annotation"]["path"] == str(annotation.resolve())
    backend_document = validated.to_backend_document()
    assert backend_document["gene_sets"] == validated.gene_sets
    assert backend_document["gene_sets"] is not validated.gene_sets


def test_analysis_request_v10_rejects_gene_sets(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        gene_sets={},
    )

    with pytest.raises(InvalidRequestError) as caught:
        load_analysis_request(request)

    assert caught.value.details["context"] == "analysis request"
    assert caught.value.details["unknown_keys"] == ["gene_sets"]


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


def test_analysis_request_reports_all_supported_public_versions(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.3",
    )

    with pytest.raises(InvalidRequestError) as caught:
        load_analysis_request(request)

    assert caught.value.details == {
        "observed": "1.3",
        "supported": ["1.0", "1.1", "1.2"],
    }


def test_analysis_request_v12_defaults_to_edger_without_changing_private_protocol(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
    )

    validated = load_analysis_request(request)

    assert validated.request_schema_version == "1.2"
    assert validated.backend == "edger_ql"
    assert validated.deseq2 is None
    backend_document = validated.to_backend_document()
    assert backend_document["schema_version"] == "1.0"
    assert backend_document["kind"] == "edger_ql_backend_request"
    assert backend_document["analysis_request"] == {
        "path": str(request.resolve()),
        "sha256": _sha256(request),
    }
    assert "backend" not in backend_document
    assert "deseq2" not in backend_document
    assert "display" not in backend_document
    with pytest.raises(InvalidRequestError, match="normalized DESeq2"):
        validated.to_deseq2_backend_document()


def test_analysis_request_v12_edger_accepts_existing_display_extension(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="edger_ql",
        display={
            "fdr_threshold": 0.05,
            "pca_top_n": 500,
            "pca_components": [1, 2],
        },
    )

    validated = load_analysis_request(request)

    assert validated.backend == "edger_ql"
    assert validated.display == {
        "fdr_threshold": 0.05,
        "pca_top_n": 500,
        "pca_components": [1, 2],
    }


@pytest.mark.parametrize("schema_version", ["1.0", "1.1"])
@pytest.mark.parametrize("field", ["backend", "deseq2"])
def test_legacy_analysis_request_versions_reject_v12_fields(
    tmp_path: Path,
    schema_version: str,
    field: str,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version=schema_version,
        display=(
            {
                "fdr_threshold": 0.05,
                "pca_top_n": 500,
                "pca_components": [1, 2],
            }
            if schema_version == "1.1"
            else None
        ),
        backend="edger_ql" if field == "backend" else _MISSING,
        deseq2=(
            {"test_mode": "wald", "shrinkage": "none"}
            if field == "deseq2"
            else _MISSING
        ),
    )

    with pytest.raises(InvalidRequestError) as caught:
        load_analysis_request(request)

    assert caught.value.details["unknown_keys"] == [field]


@pytest.mark.parametrize(
    "backend",
    ["edgeR_ql", "edger", "deseq2 ", "DESeq2", "", None, True, 1, {}],
)
def test_analysis_request_v12_rejects_nonexact_backend_values(
    tmp_path: Path,
    backend: object,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend=backend,
    )

    with pytest.raises(InvalidRequestError, match="'backend'"):
        load_analysis_request(request)


def test_analysis_request_v12_edger_rejects_deseq2_configuration(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        deseq2={"test_mode": "wald", "shrinkage": "none"},
    )

    with pytest.raises(InvalidRequestError, match="edgeR request"):
        load_analysis_request(request)


def test_analysis_request_v12_edger_gene_sets_still_require_display(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        gene_sets={},
    )

    with pytest.raises(InvalidRequestError, match="requires the atomic display"):
        load_analysis_request(request)


def test_analysis_request_v12_accepts_deseq2_wald_and_serializes_separately(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        deseq2={"test_mode": "wald", "shrinkage": "none"},
    )

    validated = load_analysis_request(request)

    assert validated.backend == "deseq2"
    assert validated.display is None
    assert validated.gene_sets is None
    assert validated.deseq2 == {"test_mode": "wald", "shrinkage": "none"}
    with pytest.raises(InvalidRequestError, match="non-edgeR"):
        validated.to_backend_document()
    backend_document = validated.to_deseq2_backend_document()
    assert backend_document["schema_version"] == "1.0"
    assert backend_document["kind"] == "deseq2_backend_request"
    assert backend_document["analysis_request"] == {
        "path": str(request.resolve()),
        "sha256": _sha256(request),
        "schema_version": "1.2",
        "backend": "deseq2",
    }
    assert backend_document["deseq2"] == validated.deseq2
    assert backend_document["deseq2"] is not validated.deseq2
    assert "display" not in backend_document
    assert "gene_sets" not in backend_document


@pytest.mark.parametrize("field", ["display", "gene_sets"])
def test_deseq2_d1_rejects_edger_extensions(
    tmp_path: Path,
    field: str,
) -> None:
    bundle = _validated_bundle(tmp_path)
    extra: dict[str, object] = {
        field: (
            {
                "fdr_threshold": 0.05,
                "pca_top_n": 500,
                "pca_components": [1, 2],
            }
            if field == "display"
            else {}
        )
    }
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        deseq2={"test_mode": "wald", "shrinkage": "none"},
        **extra,
    )

    with pytest.raises(InvalidRequestError, match="does not support") as caught:
        load_analysis_request(request)

    assert caught.value.details["incompatible_fields"] == [field]


def test_deseq2_configuration_is_required(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
    )

    with pytest.raises(InvalidRequestError, match="requires the 'deseq2'"):
        load_analysis_request(request)


@pytest.mark.parametrize(
    "configuration",
    [
        None,
        {},
        {"test_mode": "Wald", "shrinkage": "none"},
        {"test_mode": "wald"},
        {"test_mode": "wald", "shrinkage": "normal"},
        {"test_mode": "wald", "shrinkage": {}},
        {"test_mode": "wald", "shrinkage": "none", "reduced": {}},
        {"test_mode": "wald", "shrinkage": "none", "unexpected": True},
    ],
)
def test_deseq2_configuration_rejects_noncanonical_wald_shapes(
    tmp_path: Path,
    configuration: object,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        deseq2=configuration,
    )

    with pytest.raises(InvalidRequestError):
        load_analysis_request(request)


def test_deseq2_lrt_accepts_one_zero_threshold_reporting_contrast(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    full_design = {
        "intercept": True,
        "terms": ["batch", "condition"],
        "variables": {
            "batch": {"type": "factor", "levels": ["one", "two"]},
            "condition": {"type": "factor", "levels": ["control", "treated"]},
        },
    }
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        design=full_design,
        deseq2={
            "test_mode": "lrt",
            "shrinkage": "none",
            "reduced": {"intercept": True, "terms": ["batch"]},
        },
    )

    validated = load_analysis_request(request)

    assert validated.deseq2 == {
        "test_mode": "lrt",
        "shrinkage": "none",
        "reduced": {"intercept": True, "terms": ["batch"]},
    }


def test_deseq2_lrt_requires_reduced_design(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        deseq2={"test_mode": "lrt", "shrinkage": "none"},
    )

    with pytest.raises(InvalidRequestError) as caught:
        load_analysis_request(request)

    assert caught.value.details["missing_keys"] == ["reduced"]


def test_deseq2_lrt_rejects_multiple_reporting_contrasts(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        contrasts=[
            {
                "contrast_id": "one",
                "weights": {"conditiontreated": 1},
                "lfc_threshold": 0,
            },
            {
                "contrast_id": "two",
                "weights": {"conditiontreated": -1},
                "lfc_threshold": 0,
            },
        ],
        deseq2={
            "test_mode": "lrt",
            "shrinkage": "none",
            "reduced": {"intercept": True, "terms": []},
        },
    )

    with pytest.raises(InvalidRequestError, match="exactly one"):
        load_analysis_request(request)


def test_deseq2_lrt_rejects_lfc_threshold(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        contrasts=[
            {
                "contrast_id": "treated_vs_control",
                "weights": {"conditiontreated": 1},
                "lfc_threshold": 1,
            }
        ],
        deseq2={
            "test_mode": "lrt",
            "shrinkage": "none",
            "reduced": {"intercept": True, "terms": []},
        },
    )

    with pytest.raises(InvalidRequestError, match="cannot use an LFC threshold"):
        load_analysis_request(request)


@pytest.mark.parametrize(
    "reduced",
    [
        {"intercept": False, "terms": []},
        {"intercept": True, "terms": ["batch", "condition"]},
        {"intercept": True, "terms": ["unknown"]},
        {"intercept": True, "terms": ["batch", "batch"]},
        {"intercept": True, "terms": ["batch"], "unexpected": True},
    ],
)
def test_deseq2_lrt_rejects_non_nested_or_noncanonical_reduced_designs(
    tmp_path: Path,
    reduced: object,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        design={
            "intercept": True,
            "terms": ["batch", "condition"],
            "variables": {
                "batch": {"type": "factor", "levels": ["one", "two"]},
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                },
            },
        },
        deseq2={"test_mode": "lrt", "shrinkage": "none", "reduced": reduced},
    )

    with pytest.raises(InvalidRequestError):
        load_analysis_request(request)


def test_deseq2_lrt_normalizes_reduced_terms_to_full_design_order(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        design={
            "intercept": True,
            "terms": ["batch", "sex", "condition"],
            "variables": {
                "batch": {"type": "factor", "levels": ["one", "two"]},
                "sex": {"type": "factor", "levels": ["female", "male"]},
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                },
            },
        },
        deseq2={
            "test_mode": "lrt",
            "shrinkage": "none",
            "reduced": {"intercept": True, "terms": ["condition", "batch"]},
        },
    )

    validated = load_analysis_request(request)

    assert validated.deseq2 is not None
    assert validated.deseq2["reduced"]["terms"] == ["batch", "condition"]


def test_deseq2_lrt_no_intercept_reduced_design_cannot_be_empty(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        design={
            "intercept": False,
            "terms": ["condition"],
            "variables": {
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                }
            },
        },
        deseq2={
            "test_mode": "lrt",
            "shrinkage": "none",
            "reduced": {"intercept": False, "terms": []},
        },
    )

    with pytest.raises(InvalidRequestError, match="must retain at least one term"):
        load_analysis_request(request)


def test_deseq2_apeglm_accepts_exact_positive_single_coefficient_contrasts(
    tmp_path: Path,
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        deseq2={"test_mode": "wald", "shrinkage": "apeglm"},
    )

    validated = load_analysis_request(request)

    assert validated.deseq2 == {"test_mode": "wald", "shrinkage": "apeglm"}


def test_deseq2_wald_accepts_formal_positive_lfc_threshold(tmp_path: Path) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights={"conditiontreated": 1},
        schema_version="1.2",
        backend="deseq2",
        contrasts=[
            {
                "contrast_id": "treated_vs_control_min_1",
                "weights": {"conditiontreated": 1},
                "lfc_threshold": 1,
            }
        ],
        deseq2={"test_mode": "wald", "shrinkage": "none"},
    )

    validated = load_analysis_request(request)

    assert validated.contrasts[0]["lfc_threshold"] == 1.0


@pytest.mark.parametrize(
    "weights",
    [
        {"conditiontreated": -1},
        {"conditiontreated": 2},
        {"conditiontreated": 1, "batchtwo": 0},
        {"(Intercept)": 1},
        {"Intercept": 1},
    ],
)
def test_deseq2_apeglm_rejects_any_nonexact_coefficient_selection(
    tmp_path: Path,
    weights: dict[str, float],
) -> None:
    bundle = _validated_bundle(tmp_path)
    request = _analysis_request(
        tmp_path,
        bundle,
        weights=weights,
        schema_version="1.2",
        backend="deseq2",
        deseq2={"test_mode": "wald", "shrinkage": "apeglm"},
    )

    with pytest.raises(InvalidRequestError, match="apeglm shrinkage"):
        load_analysis_request(request)


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
