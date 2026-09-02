"""Adversarial boundaries for the evidence-coupled edgeR adapter."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from rnaseq_downstream import edger_backend
from rnaseq_downstream.errors import (
    BackendFailedError,
    InputIntegrityError,
    InvalidRequestError,
)

pytestmark = pytest.mark.unit


def test_materialization_rejects_source_bytes_changed_after_validation(
    tmp_path: Path,
) -> None:
    source = tmp_path / "validated-counts.tsv"
    validated_bytes = b"gene_id\ts1\nENSG1\t10\n"
    source.write_bytes(validated_bytes)
    artifact = {
        "role": "featurecounts.matrix",
        "path": str(source),
        "sha256": hashlib.sha256(validated_bytes).hexdigest(),
        "size_bytes": len(validated_bytes),
    }

    changed_bytes = b"gene_id\ts1\nENSG1\t99\n"
    assert len(changed_bytes) == len(validated_bytes)
    source.write_bytes(changed_bytes)

    with pytest.raises(InputIntegrityError) as caught:
        edger_backend._copy_evidence_snapshot(
            artifact,
            tmp_path / "private" / "counts.tsv",
        )

    assert caught.value.to_dict()["code"] == "INPUT_INTEGRITY_FAILED"
    assert caught.value.details["expected_sha256"] == artifact["sha256"]
    assert (
        caught.value.details["observed_sha256"]
        == hashlib.sha256(changed_bytes).hexdigest()
    )
    assert caught.value.details["source_identity_stable"] is True


def test_resolve_output_rejects_existing_directory_without_overwriting(
    tmp_path: Path,
) -> None:
    output = tmp_path / "existing-results"
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("preserve", encoding="utf-8")
    inode = output.stat().st_ino

    with pytest.raises(InvalidRequestError) as caught:
        edger_backend._resolve_output(output)

    assert caught.value.to_dict()["code"] == "INVALID_REQUEST"
    assert output.stat().st_ino == inode
    assert marker.read_text(encoding="utf-8") == "preserve"


@pytest.mark.parametrize("live_target", [False, True])
def test_resolve_output_rejects_symlink_without_touching_referent(
    tmp_path: Path,
    live_target: bool,
) -> None:
    referent = tmp_path / "referent"
    if live_target:
        referent.mkdir()
        (referent / "keep.txt").write_text("preserve", encoding="utf-8")
    output = tmp_path / "results-link"
    output.symlink_to(referent, target_is_directory=True)

    with pytest.raises(InvalidRequestError) as caught:
        edger_backend._resolve_output(output)

    assert caught.value.to_dict()["code"] == "INVALID_REQUEST"
    assert output.is_symlink()
    if live_target:
        assert (referent / "keep.txt").read_text(encoding="utf-8") == "preserve"
    else:
        assert not referent.exists()


@pytest.mark.parametrize(
    ("returncode", "response"),
    [
        pytest.param(
            4,
            {"status": "success", "backend": edger_backend.BACKEND_NAME},
            id="nonzero-return-with-success-status",
        ),
        pytest.param(
            0,
            {
                "status": "error",
                "errors": [
                    {
                        "code": "BACKEND_FAILED",
                        "message": "backend rejected the request",
                        "details": {"reason": "synthetic"},
                    }
                ],
            },
            id="zero-return-with-error-status",
        ),
        pytest.param(
            3,
            {
                "status": "error",
                "errors": [
                    {
                        "code": "UNRECOGNIZED_BACKEND_CODE",
                        "message": "unrecognized failure",
                        "details": {},
                    }
                ],
            },
            id="unrecognized-structured-error-code",
        ),
    ],
)
def test_invoke_r_contract_mismatches_are_never_success(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    returncode: int,
    response: dict[str, object],
) -> None:
    completed = SimpleNamespace(
        returncode=returncode,
        stdout=json.dumps(response, allow_nan=False) + "\n",
        stderr="synthetic backend diagnostic",
    )
    monkeypatch.setattr(edger_backend, "_r_script_path", lambda: tmp_path / "edger.R")
    monkeypatch.setattr(
        edger_backend.subprocess, "run", lambda *args, **kwargs: completed
    )

    with pytest.raises(BackendFailedError) as caught:
        edger_backend._invoke_r(
            tmp_path / "request.json",
            tmp_path / "result-stage",
            rscript="Rscript",
            r_library=None,
        )

    assert caught.value.to_dict()["code"] == "BACKEND_FAILED"


@pytest.mark.parametrize("observed_backend", [None, "DESeq2", "edgeR_QL "])
def test_success_response_requires_exact_backend_identity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    observed_backend: str | None,
) -> None:
    response = {"status": "success", "backend": observed_backend, "data": {}}
    completed = SimpleNamespace(
        returncode=0,
        stdout=json.dumps(response, allow_nan=False) + "\n",
        stderr="",
    )
    monkeypatch.setattr(edger_backend, "_r_script_path", lambda: tmp_path / "edger.R")
    monkeypatch.setattr(
        edger_backend.subprocess, "run", lambda *args, **kwargs: completed
    )

    with pytest.raises(BackendFailedError) as caught:
        edger_backend._invoke_r(
            tmp_path / "request.json",
            tmp_path / "result-stage",
            rscript="Rscript",
            r_library=None,
        )

    assert caught.value.to_dict()["code"] == "BACKEND_FAILED"
    assert caught.value.details["observed_backend"] == observed_backend
