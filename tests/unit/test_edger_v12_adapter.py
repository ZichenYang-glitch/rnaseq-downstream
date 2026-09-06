"""Boundary tests for the schema-1.2 wrapper around the frozen edgeR adapter."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from rnaseq_downstream import edger_backend_v12


@pytest.mark.unit
def test_v12_execute_without_display_delegates_to_frozen_adapter(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    expected = {"status": "success", "source": "frozen"}
    observed: dict[str, object] = {}

    def fake_execute(*args: object, **kwargs: object) -> dict[str, object]:
        observed["args"] = args
        observed["kwargs"] = kwargs
        return expected

    monkeypatch.setattr(edger_backend_v12, "_execute_document", fake_execute)
    document = {"schema_version": "1.0"}

    result = edger_backend_v12._execute_v12_document(
        document,
        tmp_path / "result",
        tmp_path / "workspace",
        rscript="Rscript",
        r_library=None,
        display_configuration=None,
    )

    assert result is expected
    assert observed["args"] == (
        document,
        tmp_path / "result",
        tmp_path / "workspace",
    )
    assert observed["kwargs"] == {"rscript": "Rscript", "r_library": None}


@pytest.mark.unit
def test_v12_display_execution_binds_public_request_schema(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}
    runtime = dict(edger_backend_v12._EXPECTED_RUNTIME)

    monkeypatch.setattr(
        edger_backend_v12, "_write_private_json", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        edger_backend_v12,
        "_invoke_r",
        lambda *_args, **_kwargs: ({"data": {"runtime_identity": runtime}}, ""),
    )
    monkeypatch.setattr(
        edger_backend_v12,
        "_verify_result_stage",
        lambda *_args, **_kwargs: ({"kind": "analysis"}, []),
    )

    def fake_display(**kwargs: object) -> list[dict[str, object]]:
        captured.update(kwargs)
        return []

    monkeypatch.setattr(edger_backend_v12, "build_display_bundle", fake_display)
    monkeypatch.setattr(
        edger_backend_v12, "_fsync_directory", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        edger_backend_v12,
        "_verify_complete_public_stage",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        edger_backend_v12, "_publish_noreplace", lambda *_args, **_kwargs: None
    )
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    document = {
        "schema_version": "1.0",
        "execution_scope": "validated_p0_input",
    }

    result = edger_backend_v12._execute_v12_document(
        document,
        tmp_path / "result",
        workspace,
        rscript="Rscript",
        r_library=None,
        display_configuration={"pca_top_n": 500},
    )

    assert captured["analysis_request_schema_version"] == "1.2"
    assert document["display_export"] == {"path": str(workspace / "display-logcpm.tsv")}
    assert result["status"] == "success"
    assert result["backend"] == "edgeR_QL"


@pytest.mark.unit
@pytest.mark.parametrize(
    ("schema_version", "backend"),
    [("1.1", "edger_ql"), ("1.2", "deseq2")],
)
def test_v12_adapter_rejects_wrong_public_route_before_output_resolution(
    monkeypatch: pytest.MonkeyPatch,
    schema_version: str,
    backend: str,
) -> None:
    monkeypatch.setattr(
        edger_backend_v12,
        "_resolve_output",
        lambda *_args, **_kwargs: pytest.fail("output must not be resolved"),
    )
    validated = SimpleNamespace(
        request_schema_version=schema_version,
        backend=backend,
    )

    with pytest.raises(
        edger_backend_v12.InvalidRequestError,
        match="requires an edgeR request schema 1.2",
    ):
        edger_backend_v12._run_validated_edger_v12(
            validated,  # type: ignore[arg-type]
            "results",
            rscript="Rscript",
            r_library=None,
        )
