"""Boundaries for the conditional C2 edgeR adapter protocol."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from rnaseq_downstream import edger_backend
from rnaseq_downstream.errors import BackendFailedError, InputIntegrityError

pytestmark = pytest.mark.unit


def _source(path: Path, content: bytes, **identity: object) -> dict[str, object]:
    path.write_bytes(content)
    return {
        "path": str(path),
        "declared_path": path.name,
        "sha256": hashlib.sha256(content).hexdigest(),
        "size_bytes": len(content),
        **identity,
    }


def _gene_sets(tmp_path: Path) -> dict[str, object]:
    return {
        "gmt": _source(
            tmp_path / "sets.gmt",
            b"set-a\tdescription\tA\tB\n",
            collection="fixture",
            version="1",
            identifier_type="symbol",
        ),
        "annotation": _source(
            tmp_path / "annotation.tsv",
            b"gene_id\tsymbol\nENSG1\tA\nENSG2\tB\n",
            name="fixture-annotation",
            version="1",
            gene_id_column="gene_id",
            symbol_column="symbol",
        ),
        "minimum_tested_genes": 2,
    }


def test_gene_set_sources_are_recaptured_into_private_workspace(
    tmp_path: Path,
) -> None:
    normalized = _gene_sets(tmp_path)
    private_root = tmp_path / "private"

    rewritten, snapshots = edger_backend._materialize_gene_set_inputs(
        normalized, private_root
    )

    assert Path(rewritten["gmt"]["path"]).read_bytes() == (
        tmp_path / "sets.gmt"
    ).read_bytes()
    assert Path(rewritten["annotation"]["path"]).read_bytes() == (
        tmp_path / "annotation.tsv"
    ).read_bytes()
    assert rewritten["gmt"]["declared_path"] == "sets.gmt"
    assert [item["role"] for item in snapshots] == [
        "pathways.gmt",
        "pathways.annotation",
    ]
    assert [item["private_relative_path"] for item in snapshots] == [
        "gene-sets/sets.gmt",
        "gene-sets/annotation.tsv",
    ]


def test_gene_set_recapture_rejects_changed_source_bytes(tmp_path: Path) -> None:
    normalized = _gene_sets(tmp_path)
    (tmp_path / "sets.gmt").write_bytes(b"set-a\tdescription\tA\tC\n")

    with pytest.raises(InputIntegrityError) as caught:
        edger_backend._materialize_gene_set_inputs(
            normalized, tmp_path / "private"
        )

    assert caught.value.details["role"] == "pathways.gmt"


def test_pathway_success_response_requires_schema_1_1(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    response = {
        "schema_version": "1.0",
        "status": "success",
        "backend": edger_backend.BACKEND_NAME,
        "data": {},
    }
    completed = SimpleNamespace(
        returncode=0,
        stdout=json.dumps(response) + "\n",
        stderr="",
    )
    monkeypatch.setattr(edger_backend, "_r_script_path", lambda: tmp_path / "edger.R")
    monkeypatch.setattr(
        edger_backend.subprocess, "run", lambda *args, **kwargs: completed
    )

    with pytest.raises(BackendFailedError) as caught:
        edger_backend._invoke_r(
            tmp_path / "request.json",
            tmp_path / "results",
            rscript="Rscript",
            r_library=None,
            expected_schema_version="1.1",
        )

    assert caught.value.details["observed_schema_version"] == "1.0"
    assert caught.value.details["expected_schema_version"] == "1.1"


def _write_stage(tmp_path: Path, *, schema_version: str) -> Path:
    stage = tmp_path / "stage"
    stage.mkdir()
    member_names = {
        "analysis.json",
        "coefficients.tsv",
        "design.tsv",
        "results.tsv",
    }
    if schema_version == "1.1":
        member_names.add("pathway_results.tsv")
    for name in member_names:
        if name == "analysis.json":
            payload = (json.dumps({"schema_version": schema_version}) + "\n").encode()
        else:
            payload = b"fixture\n"
        (stage / name).write_bytes(payload)
    members = []
    for name in sorted(member_names):
        payload = (stage / name).read_bytes()
        members.append(
            {
                "path": name,
                "sha256": hashlib.sha256(payload).hexdigest(),
                "size_bytes": len(payload),
            }
        )
    manifest = {
        "schema_version": schema_version,
        "kind": "edger_ql_backend_manifest",
        "backend": edger_backend.BACKEND_NAME,
        "runtime_identity": dict(edger_backend._EXPECTED_RUNTIME),
        "execution_scope": "validated_p0_input",
        "members": members,
    }
    (stage / "backend_manifest.json").write_text(
        json.dumps(manifest) + "\n", encoding="utf-8"
    )
    return stage


def test_pathway_stage_requires_and_binds_the_sixth_result(tmp_path: Path) -> None:
    stage = _write_stage(tmp_path, schema_version="1.1")

    analysis, artifacts = edger_backend._verify_result_stage(
        stage,
        execution_scope="validated_p0_input",
        expected_schema_version="1.1",
    )

    assert analysis["schema_version"] == "1.1"
    assert {item["relative_path"] for item in artifacts} == {
        "analysis.json",
        "backend_manifest.json",
        "coefficients.tsv",
        "design.tsv",
        "pathway_results.tsv",
        "results.tsv",
    }


def test_schema_1_0_stage_rejects_pathway_result(tmp_path: Path) -> None:
    stage = _write_stage(tmp_path, schema_version="1.1")

    with pytest.raises(BackendFailedError, match="inventory"):
        edger_backend._verify_result_stage(
            stage,
            execution_scope="validated_p0_input",
            expected_schema_version="1.0",
        )
