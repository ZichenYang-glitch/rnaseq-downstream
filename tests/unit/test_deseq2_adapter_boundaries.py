"""Fail-closed process and publication boundaries for the DESeq2 adapter."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from rnaseq_downstream import deseq2_backend
from rnaseq_downstream.errors import (
    BackendFailedError,
    CountValuesInvalidError,
    InvalidRequestError,
)

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("code", "returncode", "error_type"),
    [
        ("INVALID_REQUEST", 2, InvalidRequestError),
        ("COUNT_VALUES_INVALID", 3, CountValuesInvalidError),
        ("BACKEND_FAILED", 4, BackendFailedError),
    ],
)
def test_r_error_codes_map_to_stable_python_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    code: str,
    returncode: int,
    error_type: type[Exception],
) -> None:
    response = {
        "schema_version": "1.0",
        "status": "error",
        "backend": "DESeq2",
        "errors": [
            {
                "code": code,
                "message": "synthetic failure",
                "details": {"reason": "synthetic"},
            }
        ],
    }
    completed = SimpleNamespace(
        returncode=returncode,
        stdout=json.dumps(response) + "\n",
        stderr="private diagnostic",
    )
    monkeypatch.setattr(
        deseq2_backend.subprocess, "run", lambda *args, **kwargs: completed
    )
    monkeypatch.setattr(deseq2_backend, "_r_script_path", lambda: tmp_path / "deseq2.R")

    with pytest.raises(error_type) as caught:
        deseq2_backend._invoke_r(
            tmp_path / "request.json",
            tmp_path / "results",
            rscript="Rscript",
            r_library=None,
        )

    assert caught.value.details["reason"] == "synthetic"
    assert caught.value.details["backend_returncode"] == returncode


@pytest.mark.parametrize(
    "stdout",
    [
        "",
        "{}\n{}\n",
        "not-json\n",
        "[]\n",
        '{"status":"success","backend":"edgeR_QL","schema_version":"1.0"}\n',
    ],
)
def test_r_stdout_or_identity_mismatch_never_succeeds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    stdout: str,
) -> None:
    completed = SimpleNamespace(returncode=0, stdout=stdout, stderr="")
    monkeypatch.setattr(
        deseq2_backend.subprocess, "run", lambda *args, **kwargs: completed
    )
    monkeypatch.setattr(deseq2_backend, "_r_script_path", lambda: tmp_path / "deseq2.R")

    with pytest.raises(BackendFailedError):
        deseq2_backend._invoke_r(
            tmp_path / "request.json",
            tmp_path / "results",
            rscript="Rscript",
            r_library=None,
        )


def test_deseq2_adapter_rejects_wrong_normalized_backend(tmp_path: Path) -> None:
    validated = SimpleNamespace(
        backend="edger_ql",
        display=None,
        gene_sets=None,
    )

    with pytest.raises(InvalidRequestError) as caught:
        deseq2_backend._run_validated_deseq2(
            validated,  # type: ignore[arg-type]
            tmp_path / "results",
            rscript="Rscript",
            r_library=None,
        )

    assert caught.value.details == {"backend": "edger_ql"}
    assert not (tmp_path / "results").exists()


def test_staged_summary_failure_prevents_publication(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from rnaseq_downstream import run_summary

    stage = tmp_path / "stage"
    stage.mkdir()

    def reject(_path: Path) -> dict[str, object]:
        raise InvalidRequestError("synthetic verifier failure")

    monkeypatch.setattr(run_summary, "summarize_run", reject)

    with pytest.raises(BackendFailedError) as caught:
        deseq2_backend._verify_complete_public_stage(stage)

    assert caught.value.details["reason"] == "staged_bundle_verification_failed"
    assert caught.value.details["cause_code"] == "INVALID_REQUEST"


def test_process_receipt_must_match_manifested_analysis() -> None:
    analysis = {
        "execution_scope": "validated_p1_deseq2_input",
        "runtime_identity": dict(deseq2_backend._EXPECTED_RUNTIME),
        "input_semantics": "featurecounts_integer",
        "route_observed": {"constructor": "DESeqDataSetFromMatrix"},
        "design": {"columns": ["(Intercept)"], "rank": 1, "residual_df": 3},
        "genes": {
            "total": 10,
            "status_counts": {
                "filtered": 0,
                "not_tested": 0,
                "not_estimable": 0,
                "failed": 0,
                "tested": 10,
            },
        },
        "test": {"mode": "wald", "shrinkage": "none", "reduced": None},
        "defaults": {"fit_type_requested": "parametric"},
        "contrasts": [{"contrast_id": "treated_vs_control"}],
        "reporting_effect": [{"contrast_id": "treated_vs_control"}],
    }
    receipt = {
        "execution_scope": analysis["execution_scope"],
        "runtime_identity": analysis["runtime_identity"],
        "input_semantics": analysis["input_semantics"],
        "route_observed": analysis["route_observed"],
        "design_columns": analysis["design"]["columns"],
        "design_rank": analysis["design"]["rank"],
        "residual_df": analysis["design"]["residual_df"],
        "gene_count": analysis["genes"]["total"],
        "result_status_counts": analysis["genes"]["status_counts"],
        "test": analysis["test"],
        "defaults": analysis["defaults"],
        "contrasts": analysis["contrasts"],
        "reporting_effect": analysis["reporting_effect"],
    }

    assert deseq2_backend._reconcile_response_data(receipt, analysis) == receipt
    receipt["gene_count"] = 9
    with pytest.raises(BackendFailedError) as caught:
        deseq2_backend._reconcile_response_data(receipt, analysis)

    assert caught.value.details["reason"] == "backend_receipt_analysis_mismatch"
