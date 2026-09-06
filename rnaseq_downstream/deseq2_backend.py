"""Evidence-coupled adapter for the locked DESeq2 1.52 backend."""

from __future__ import annotations

import hashlib
from importlib import resources
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Mapping

from .analysis_contract_v12 import ValidatedAnalysisRequest, load_analysis_request
from .edger_backend import (
    _fsync_directory,
    _hash_and_fsync_output,
    _materialize_validated_inputs,
    _publish_noreplace,
    _resolve_output,
)
from .errors import (
    BackendFailedError,
    ContrastNotEstimableError,
    CountValuesInvalidError,
    DesignRankDeficientError,
    InputReadError,
    InvalidRequestError,
    OutputWriteError,
    ToolkitError,
)
from .provenance import read_json_object

BACKEND_NAME = "DESeq2"
BACKEND_OUTPUT_SCHEMA_VERSION = "1.0"
EXECUTION_SCOPE = "validated_p1_deseq2_input"
_EXPECTED_RUNTIME = {
    "R": "4.6.1",
    "Bioconductor": "3.23",
    "BiocVersion_package": "3.23.1",
    "DESeq2": "1.52.0",
    "apeglm": "1.34.0",
    "tximport": "1.40.0",
}
_OUTPUT_MEMBERS = {
    "analysis.json",
    "coefficients.tsv",
    "design.tsv",
    "results.tsv",
}


def _canonical_bytes(document: Mapping[str, Any]) -> bytes:
    try:
        return json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    except (TypeError, ValueError, OverflowError) as error:
        raise InvalidRequestError(
            "The normalized DESeq2 request is not strict JSON.",
            details={"cause_type": type(error).__name__},
            cause=error,
        ) from error


def _write_private_json(path: Path, document: Mapping[str, Any]) -> None:
    payload = _canonical_bytes(document) + b"\n"
    try:
        with path.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        path.chmod(0o400)
    except OSError as error:
        raise OutputWriteError(
            "The private DESeq2 request could not be written.",
            path=path,
            operation="write_private_deseq2_request",
            cause=error,
        ) from error


def _r_script_path() -> Path:
    resource = resources.files("rnaseq_downstream").joinpath("r_scripts", "deseq2.R")
    path = Path(str(resource))
    if not path.is_file():
        raise BackendFailedError(
            "The packaged DESeq2 backend resource is missing.",
            details={"resource": "r_scripts/deseq2.R"},
        )
    return path


def _parse_backend_stdout(stdout: str, *, returncode: int) -> dict[str, Any]:
    lines = stdout.splitlines()
    if len(lines) != 1 or not lines[0]:
        raise BackendFailedError(
            "The DESeq2 R backend violated its one-document stdout contract.",
            details={"returncode": returncode, "stdout_line_count": len(lines)},
        )
    try:
        response = json.loads(
            lines[0],
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-standard constant {value}")
            ),
        )
    except (json.JSONDecodeError, ValueError) as error:
        raise BackendFailedError(
            "The DESeq2 R backend returned invalid JSON.",
            details={"returncode": returncode, "cause_type": type(error).__name__},
            cause=error,
        ) from error
    if not isinstance(response, dict):
        raise BackendFailedError(
            "The DESeq2 R backend response root is not an object.",
            details={"returncode": returncode},
        )
    return response


def _raise_backend_response(
    response: Mapping[str, Any],
    *,
    returncode: int,
    stderr: str,
) -> None:
    errors = response.get("errors")
    if (
        not isinstance(errors, list)
        or len(errors) != 1
        or not isinstance(errors[0], Mapping)
    ):
        raise BackendFailedError(
            "The failed DESeq2 backend response lacks one structured error.",
            details={"returncode": returncode, "stderr": stderr[-8000:]},
        )
    error = errors[0]
    code = error.get("code")
    message = str(error.get("message", "The DESeq2 backend failed."))
    details = error.get("details")
    if not isinstance(details, Mapping):
        details = {}
    normalized_details = {
        **dict(details),
        "backend_returncode": returncode,
        "backend_stderr": stderr[-8000:],
    }
    if code == "DESIGN_RANK_DEFICIENT":
        raise DesignRankDeficientError(message, details=normalized_details)
    if code == "CONTRAST_NOT_ESTIMABLE":
        raise ContrastNotEstimableError(message, details=normalized_details)
    if code == "COUNT_VALUES_INVALID":
        raise CountValuesInvalidError(message, details=normalized_details)
    if code == "INVALID_REQUEST":
        raise InvalidRequestError(message, details=normalized_details)
    raise BackendFailedError(message, details=normalized_details)


def _invoke_r(
    backend_request: Path,
    result_stage: Path,
    *,
    rscript: str | Path,
    r_library: str | Path | None,
) -> tuple[dict[str, Any], str]:
    environment = os.environ.copy()
    if r_library is not None:
        try:
            library_path = Path(r_library).resolve(strict=True)
        except (OSError, RuntimeError, TypeError) as error:
            raise InvalidRequestError(
                "The requested R package library cannot be resolved.",
                details={
                    "r_library": str(r_library),
                    "cause_type": type(error).__name__,
                },
                cause=error,
            ) from error
        if not library_path.is_dir():
            raise InvalidRequestError(
                "The requested R package library is not a directory.",
                details={"r_library": str(library_path)},
            )
        environment["R_LIBS_USER"] = str(library_path)
    command = [
        str(rscript),
        "--vanilla",
        str(_r_script_path()),
        str(backend_request),
        str(result_stage),
    ]
    try:
        completed = subprocess.run(
            command,
            check=False,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="strict",
            env=environment,
        )
    except (OSError, UnicodeError) as error:
        raise BackendFailedError(
            "The independent DESeq2 R process could not be executed.",
            details={"cause_type": type(error).__name__, "rscript": str(rscript)},
            cause=error,
        ) from error
    response = _parse_backend_stdout(completed.stdout, returncode=completed.returncode)
    if completed.returncode != 0 or response.get("status") != "success":
        _raise_backend_response(
            response,
            returncode=completed.returncode,
            stderr=completed.stderr,
        )
    if (
        response.get("backend") != BACKEND_NAME
        or response.get("schema_version") != BACKEND_OUTPUT_SCHEMA_VERSION
    ):
        raise BackendFailedError(
            "The DESeq2 R response has an unexpected identity.",
            details={
                "observed_backend": response.get("backend"),
                "observed_schema_version": response.get("schema_version"),
            },
        )
    return response, completed.stderr


def _capture_output_json(path: Path, *, role: str) -> tuple[dict[str, Any], str, int]:
    try:
        before = path.stat()
        content = path.read_bytes()
        after = path.stat()
    except OSError as error:
        raise BackendFailedError(
            "A DESeq2 JSON artifact could not be captured.",
            details={"path": str(path), "role": role},
            cause=error,
        ) from error
    if (before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ) or len(content) != after.st_size:
        raise BackendFailedError(
            "A DESeq2 JSON artifact changed during verification.",
            details={"path": str(path), "role": role},
        )
    try:
        document = read_json_object(path, document_role=role, content=content)
    except (InvalidRequestError, InputReadError) as error:
        raise BackendFailedError(
            "A DESeq2 JSON artifact is not a strict JSON object.",
            details={
                "path": str(path),
                "role": role,
                "parse_error_code": error.code.value,
                "parse_error_details": dict(error.details),
            },
            cause=error,
        ) from error
    return document, hashlib.sha256(content).hexdigest(), len(content)


def _verify_result_stage(
    result_stage: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    expected_files = _OUTPUT_MEMBERS | {"backend_manifest.json"}
    try:
        entries = list(result_stage.iterdir())
    except OSError as error:
        raise BackendFailedError(
            "The DESeq2 backend did not create a readable result stage.",
            details={"result_stage": str(result_stage)},
            cause=error,
        ) from error
    names = {entry.name for entry in entries}
    unsafe = sorted(
        entry.name for entry in entries if entry.is_symlink() or not entry.is_file()
    )
    if names != expected_files or unsafe:
        raise BackendFailedError(
            "The DESeq2 output inventory is incomplete or unsafe.",
            details={
                "missing_files": sorted(expected_files - names),
                "unexpected_files": sorted(names - expected_files),
                "unsafe_files": unsafe,
            },
        )
    manifest_path = result_stage / "backend_manifest.json"
    manifest, manifest_digest, manifest_size = _capture_output_json(
        manifest_path, role="deseq2_backend_manifest"
    )
    if (
        manifest.get("schema_version") != BACKEND_OUTPUT_SCHEMA_VERSION
        or manifest.get("kind") != "deseq2_backend_manifest"
        or manifest.get("backend") != BACKEND_NAME
        or manifest.get("execution_scope") != EXECUTION_SCOPE
        or manifest.get("runtime_identity") != _EXPECTED_RUNTIME
    ):
        raise BackendFailedError(
            "The DESeq2 backend manifest identity is incompatible.",
            details={"manifest": str(manifest_path)},
        )
    members = manifest.get("members")
    if not isinstance(members, list):
        raise BackendFailedError("The DESeq2 manifest member inventory is invalid.")
    observed_members: set[str] = set()
    artifacts: list[dict[str, Any]] = []
    for index, member in enumerate(members):
        if not isinstance(member, Mapping):
            raise BackendFailedError(
                "A DESeq2 manifest member is not an object.",
                details={"member_index": index},
            )
        name = member.get("path")
        if name not in _OUTPUT_MEMBERS or name in observed_members:
            raise BackendFailedError(
                "A DESeq2 manifest member path is invalid.",
                details={"member_index": index, "path": name},
            )
        observed_members.add(str(name))
        digest, size = _hash_and_fsync_output(result_stage / str(name))
        if digest != member.get("sha256") or size != member.get("size_bytes"):
            raise BackendFailedError(
                "A DESeq2 output does not match its manifest.",
                details={"path": name},
            )
        artifacts.append(
            {
                "kind": "generated_analysis_artifact",
                "role": str(name).removesuffix(".json").removesuffix(".tsv"),
                "relative_path": name,
                "sha256": digest,
                "size_bytes": size,
                "media_type": (
                    "application/json"
                    if str(name).endswith(".json")
                    else "text/tab-separated-values"
                ),
                "schema_version": BACKEND_OUTPUT_SCHEMA_VERSION,
            }
        )
    if observed_members != _OUTPUT_MEMBERS:
        raise BackendFailedError(
            "The DESeq2 manifest omits required outputs.",
            details={"missing_members": sorted(_OUTPUT_MEMBERS - observed_members)},
        )
    # The manifest is not one of its own members, so synchronize it explicitly
    # before the containing directory can be published.
    _hash_and_fsync_output(manifest_path)
    artifacts.append(
        {
            "kind": "generated_analysis_artifact",
            "role": "backend_manifest",
            "relative_path": "backend_manifest.json",
            "sha256": manifest_digest,
            "size_bytes": manifest_size,
            "media_type": "application/json",
            "schema_version": BACKEND_OUTPUT_SCHEMA_VERSION,
        }
    )
    analysis, _, _ = _capture_output_json(
        result_stage / "analysis.json", role="deseq2_analysis"
    )
    if (
        analysis.get("schema_version") != BACKEND_OUTPUT_SCHEMA_VERSION
        or analysis.get("kind") != "deseq2_analysis"
        or analysis.get("backend") != BACKEND_NAME
        or analysis.get("runtime_identity") != _EXPECTED_RUNTIME
    ):
        raise BackendFailedError("The DESeq2 analysis identity is incompatible.")
    _fsync_directory(result_stage)
    return analysis, sorted(artifacts, key=lambda item: str(item["relative_path"]))


def _verify_complete_public_stage(result_stage: Path) -> None:
    try:
        from .run_summary import summarize_run

        summary = summarize_run(result_stage)
    except Exception as error:
        details: dict[str, Any] = {
            "reason": "staged_bundle_verification_failed",
            "cause_type": type(error).__name__,
        }
        if isinstance(error, ToolkitError):
            details.update(
                {
                    "cause_code": error.code.value,
                    "cause_details": dict(error.details),
                }
            )
        raise BackendFailedError(
            "The staged DESeq2 bundle failed independent verification.",
            details=details,
            cause=error,
        ) from error
    if not isinstance(summary, Mapping) or summary.get("status") != "verified_complete":
        raise BackendFailedError(
            "The staged DESeq2 bundle failed independent verification.",
            details={
                "reason": "staged_bundle_verification_incomplete",
                "observed_status": (
                    summary.get("status") if isinstance(summary, Mapping) else None
                ),
            },
        )


def _reconcile_response_data(
    response_data: Mapping[str, Any], analysis: Mapping[str, Any]
) -> dict[str, Any]:
    """Require the process receipt to repeat the manifested analysis evidence."""

    design = analysis.get("design")
    genes = analysis.get("genes")
    if not isinstance(design, Mapping) or not isinstance(genes, Mapping):
        raise BackendFailedError(
            "The DESeq2 analysis lacks response-reconciliation evidence."
        )
    expected = {
        "execution_scope": analysis.get("execution_scope"),
        "runtime_identity": analysis.get("runtime_identity"),
        "input_semantics": analysis.get("input_semantics"),
        "route_observed": analysis.get("route_observed"),
        "design_columns": design.get("columns"),
        "design_rank": design.get("rank"),
        "residual_df": design.get("residual_df"),
        "gene_count": genes.get("total"),
        "result_status_counts": genes.get("status_counts"),
        "test": analysis.get("test"),
        "defaults": analysis.get("defaults"),
        "contrasts": analysis.get("contrasts"),
        "reporting_effect": analysis.get("reporting_effect"),
    }
    if dict(response_data) != expected:
        raise BackendFailedError(
            "The DESeq2 process receipt disagrees with the manifested analysis.",
            details={"reason": "backend_receipt_analysis_mismatch"},
        )
    return expected


def _execute_document(
    document: dict[str, Any],
    target: Path,
    workspace: Path,
    *,
    rscript: str | Path,
    r_library: str | Path | None,
) -> dict[str, Any]:
    request_path = workspace / "backend_request.json"
    result_stage = workspace / "results"
    _write_private_json(request_path, document)
    response, stderr = _invoke_r(
        request_path,
        result_stage,
        rscript=rscript,
        r_library=r_library,
    )
    analysis, artifacts = _verify_result_stage(result_stage)
    response_data = response.get("data")
    if not isinstance(response_data, Mapping):
        raise BackendFailedError("The successful DESeq2 response lacks backend data.")
    verified_response_data = _reconcile_response_data(response_data, analysis)
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
                    "The complete DESeq2 bundle is visible, but parent-directory "
                    "synchronization failed."
                ),
                "details": {
                    "output_dir": str(target),
                    "cause_type": type(error).__name__,
                },
            }
        )
    for artifact in artifacts:
        artifact["path"] = str(target / str(artifact.pop("relative_path")))
    return {
        "schema_version": BACKEND_OUTPUT_SCHEMA_VERSION,
        "status": "success",
        "backend": BACKEND_NAME,
        "execution_scope": EXECUTION_SCOPE,
        "output_dir": str(target),
        "publication_status": publication_status,
        "data": verified_response_data,
        "analysis": analysis,
        "warnings": warnings,
        "errors": [],
        "artifacts": artifacts,
        "backend_stderr": stderr,
    }


def _run_validated_deseq2(
    validated: ValidatedAnalysisRequest,
    output_dir: str | Path,
    *,
    rscript: str | Path,
    r_library: str | Path | None,
) -> dict[str, Any]:
    if validated.backend != "deseq2":
        raise InvalidRequestError(
            "The DESeq2 adapter requires backend='deseq2'.",
            details={"backend": validated.backend},
        )
    if validated.display is not None or validated.gene_sets is not None:
        raise InvalidRequestError(
            "DESeq2 display and gene-set analysis are outside D1 scope.",
            details={"backend": validated.backend},
        )
    target = _resolve_output(output_dir)
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        workspace = Path(
            tempfile.mkdtemp(prefix=f".{target.name}.deseq2-", dir=target.parent)
        )
        workspace.chmod(0o700)
    except OSError as error:
        raise OutputWriteError(
            "A private DESeq2 workspace could not be created.",
            path=target.parent,
            operation="create_deseq2_workspace",
            cause=error,
        ) from error
    try:
        input_root = workspace / "inputs"
        input_root.mkdir(mode=0o700)
        rewritten_input, snapshots, _ = _materialize_validated_inputs(
            validated, input_root
        )
        document = validated.to_deseq2_backend_document()
        document.update(
            {
                "execution_scope": EXECUTION_SCOPE,
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
        result = _execute_document(
            document,
            target,
            workspace,
            rscript=rscript,
            r_library=r_library,
        )
        result["plan_id"] = validated.plan_id
        result["analysis_request_sha256"] = validated.request_sha256
        return result
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def run_deseq2(
    request_path: str | Path,
    output_dir: str | Path,
    *,
    rscript: str | Path = "Rscript",
    r_library: str | Path | None = None,
) -> dict[str, Any]:
    """Execute one validated-input DESeq2 request in an independent R process."""

    validated = load_analysis_request(request_path)
    return _run_validated_deseq2(
        validated,
        output_dir,
        rscript=rscript,
        r_library=r_library,
    )


__all__ = ["BACKEND_NAME", "run_deseq2"]
