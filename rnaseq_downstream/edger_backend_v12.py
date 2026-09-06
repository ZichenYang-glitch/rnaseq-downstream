"""Version-1.2 edgeR adapter layered over the byte-frozen P0 implementation."""

from __future__ import annotations

from pathlib import Path
import shutil
import tempfile
from typing import Any, Mapping

from .analysis_contract_v12 import ValidatedAnalysisRequest, load_analysis_request
from .display_bundle import build_display_bundle
from .edger_backend import (
    ANALYSIS_SCHEMA_VERSION,
    BACKEND_NAME,
    PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION,
    _EXPECTED_RUNTIME,
    _OUTPUT_MEMBERS,
    _execute_document,
    _fsync_directory,
    _invoke_r,
    _materialize_gene_set_inputs,
    _materialize_validated_inputs,
    _publish_noreplace,
    _resolve_output,
    _verify_complete_public_stage,
    _verify_result_stage,
    _write_private_json,
)
from .errors import (
    BackendFailedError,
    InvalidRequestError,
    OutputWriteError,
)


def _execute_v12_document(
    document: dict[str, Any],
    target: Path,
    workspace: Path,
    *,
    rscript: str | Path,
    r_library: str | Path | None,
    display_configuration: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Execute the frozen edgeR path, adding only the v1.2 display identity."""

    if display_configuration is None:
        return _execute_document(
            document,
            target,
            workspace,
            rscript=rscript,
            r_library=r_library,
        )

    expected_schema_version = str(document.get("schema_version", ""))
    if expected_schema_version not in {
        ANALYSIS_SCHEMA_VERSION,
        PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION,
    }:
        raise InvalidRequestError(
            "The private backend request schema is unsupported.",
            details={"schema_version": expected_schema_version},
        )
    has_gene_sets = "gene_sets" in document
    if has_gene_sets != (
        expected_schema_version == PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION
    ):
        raise InvalidRequestError(
            "The private backend schema and pathway request presence disagree.",
            details={
                "schema_version": expected_schema_version,
                "gene_sets_present": has_gene_sets,
            },
        )
    request_path = workspace / "backend_request.json"
    result_stage = workspace / "results"
    display_logcpm = workspace / "display-logcpm.tsv"
    document["display_export"] = {"path": str(display_logcpm)}
    _write_private_json(request_path, document)
    response, stderr = _invoke_r(
        request_path,
        result_stage,
        rscript=rscript,
        r_library=r_library,
        expected_schema_version=expected_schema_version,
    )
    analysis, artifacts = _verify_result_stage(
        result_stage,
        execution_scope=document["execution_scope"],
        expected_schema_version=expected_schema_version,
    )
    response_data = response.get("data")
    if not isinstance(response_data, Mapping):
        raise BackendFailedError("The successful R response lacks backend data.")
    if response_data.get("runtime_identity") != _EXPECTED_RUNTIME:
        raise BackendFailedError(
            "The R response does not carry the exact locked runtime identity."
        )
    display_artifacts = build_display_bundle(
        display_dir=result_stage / "display",
        logcpm_path=display_logcpm,
        core_dir=result_stage,
        core_artifacts=[
            item
            for item in artifacts
            if item.get("relative_path")
            in (_OUTPUT_MEMBERS | {"backend_manifest.json"})
        ],
        backend_document=document,
        backend_data=response_data,
        configuration=display_configuration,
        analysis_request_schema_version="1.2",
    )
    artifacts.extend(display_artifacts)
    _fsync_directory(result_stage)
    _verify_complete_public_stage(result_stage)
    _publish_noreplace(result_stage, target)
    publication_status = "durability_confirmed"
    warnings: list[dict[str, Any]] = []
    try:
        _fsync_directory(target.parent)
    except OSError as error:
        publication_status = "published_durability_unconfirmed"
        warnings.append(
            {
                "code": "OUTPUT_DURABILITY_UNCONFIRMED",
                "severity": "high",
                "message": (
                    "The complete backend bundle is visible, but parent-directory "
                    "synchronization failed."
                ),
                "details": {
                    "output_dir": str(target),
                    "cause_type": type(error).__name__,
                },
            }
        )
    for artifact in artifacts:
        artifact["path"] = str(target / artifact.pop("relative_path"))
    return {
        "schema_version": expected_schema_version,
        "status": "success",
        "backend": BACKEND_NAME,
        "execution_scope": document["execution_scope"],
        "output_dir": str(target),
        "publication_status": publication_status,
        "data": dict(response_data),
        "analysis": analysis,
        "warnings": warnings,
        "errors": [],
        "artifacts": artifacts,
        "backend_stderr": stderr,
    }


def _run_validated_edger_v12(
    validated: ValidatedAnalysisRequest,
    output_dir: str | Path,
    *,
    rscript: str | Path,
    r_library: str | Path | None,
) -> dict[str, Any]:
    if validated.request_schema_version != "1.2" or validated.backend != "edger_ql":
        raise InvalidRequestError(
            "The version-1.2 edgeR adapter requires an edgeR request schema 1.2.",
            details={
                "request_schema_version": validated.request_schema_version,
                "backend": validated.backend,
            },
        )
    target = _resolve_output(output_dir)
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        workspace = Path(
            tempfile.mkdtemp(prefix=f".{target.name}.edger-v12-", dir=target.parent)
        )
        workspace.chmod(0o700)
    except OSError as error:
        raise OutputWriteError(
            "A private edgeR workspace could not be created.",
            path=target.parent,
            operation="create_backend_workspace",
            cause=error,
        ) from error
    try:
        input_root = workspace / "inputs"
        input_root.mkdir(mode=0o700)
        rewritten_input, snapshots, _ = _materialize_validated_inputs(
            validated, input_root
        )
        document = validated.to_backend_document()
        if validated.gene_sets is not None:
            rewritten_gene_sets, pathway_snapshots = _materialize_gene_set_inputs(
                validated.gene_sets, input_root
            )
            document["schema_version"] = PATHWAY_BACKEND_OUTPUT_SCHEMA_VERSION
            document["gene_sets"] = rewritten_gene_sets
            snapshots.extend(pathway_snapshots)
        document.update(
            {
                "execution_scope": "validated_p0_input",
                "input": rewritten_input,
                "input_evidence": {
                    "kind": "validated_input_bundle",
                    "plan_id": validated.plan_id,
                    "bundle_path": str(validated.bundle_path),
                    "validation_bundle_artifacts": [
                        dict(item) for item in validated.validation_bundle_artifacts
                    ],
                    "r_input_snapshots": snapshots,
                    "digest_coupling": "private_copy_and_hash_same_source_stream",
                },
            }
        )
        result = _execute_v12_document(
            document,
            target,
            workspace,
            rscript=rscript,
            r_library=r_library,
            display_configuration=validated.display,
        )
        result["plan_id"] = validated.plan_id
        result["analysis_request_sha256"] = validated.request_sha256
        return result
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def run_edger_ql_v12(
    request_path: str | Path,
    output_dir: str | Path,
    *,
    rscript: str | Path = "Rscript",
    r_library: str | Path | None = None,
) -> dict[str, Any]:
    """Run one public schema-1.2 edgeR request through the frozen kernel."""

    validated = load_analysis_request(request_path)
    return _run_validated_edger_v12(
        validated,
        output_dir,
        rscript=rscript,
        r_library=r_library,
    )


__all__ = ["run_edger_ql_v12"]
