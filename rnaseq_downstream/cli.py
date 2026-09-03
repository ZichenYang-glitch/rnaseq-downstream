"""Non-interactive, JSON-only command-line interface."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import logging
import sys
from collections.abc import Callable, Mapping, Sequence
from typing import Any

from . import __version__
from .contracts import Status, build_envelope, write_json_document
from .errors import (
    BackendFailedError,
    ErrorCode,
    ExitCode,
    InternalToolkitError,
    InvalidRequestError,
    ToolkitError,
)

LOGGER = logging.getLogger(__name__)
COMMANDS = ("capabilities", "inspect", "validate", "run", "summarize")


@dataclass(frozen=True)
class _CommandResult:
    """Handler result before it is wrapped in the public response envelope."""

    data: Any
    warnings: tuple[Mapping[str, Any], ...] = ()
    artifacts: tuple[Mapping[str, Any], ...] = ()


class _HelpRequested(Exception):
    def __init__(self, command: str, text: str) -> None:
        super().__init__(text)
        self.command = command
        self.text = text


class _JsonArgumentParser(argparse.ArgumentParser):
    """Capture argparse output so stdout remains a single JSON document."""

    def __init__(self, *args: Any, command_name: str = "cli", **kwargs: Any) -> None:
        self.command_name = command_name
        self._captured_messages: list[str] = []
        super().__init__(*args, **kwargs)

    def _print_message(self, message: str | None, file: Any = None) -> None:
        del file
        if message:
            self._captured_messages.append(message)

    def exit(self, status: int = 0, message: str | None = None) -> None:
        text = "".join(self._captured_messages)
        if message:
            text += message
        if status == 0:
            raise _HelpRequested(self.command_name, text.rstrip())
        raise InvalidRequestError(
            "Command-line parsing failed.",
            details={"parser_message": text.strip()},
        )

    def error(self, message: str) -> None:
        raise InvalidRequestError(
            message,
            details={"usage": self.format_usage().strip()},
        )


def _build_parser() -> _JsonArgumentParser:
    parser = _JsonArgumentParser(
        prog="rnaseq-downstream",
        command_name="cli",
        description=(
            "Evidence-traceable bulk RNA-seq analysis toolkit (research preview)."
        ),
    )
    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        metavar="COMMAND",
    )

    capabilities = subparsers.add_parser(
        "capabilities",
        help="Report implemented and reserved command surfaces.",
    )
    capabilities.command_name = "capabilities"
    capabilities.set_defaults(handler=_capabilities)

    inspect = subparsers.add_parser(
        "inspect",
        help="Inspect declared inputs and provenance without writing artifacts.",
    )
    inspect.command_name = "inspect"
    inspect.add_argument("--request", required=True, metavar="REQUEST.json")
    inspect.set_defaults(handler=_inspect)

    validate = subparsers.add_parser(
        "validate",
        help="Validate input semantics and atomically archive input evidence.",
    )
    validate.command_name = "validate"
    validate.add_argument("--request", required=True, metavar="REQUEST.json")
    validate.add_argument("--output-dir", required=True, metavar="DIRECTORY")
    validate.set_defaults(handler=_validate)

    run = subparsers.add_parser(
        "run",
        help="Execute the evidence-gated edgeR v4 QL analysis path.",
    )
    run.command_name = "run"
    run.add_argument("--request", required=True, metavar="ANALYSIS_REQUEST.json")
    run.add_argument("--output-dir", required=True, metavar="DIRECTORY")
    run.add_argument("--rscript", default="Rscript", metavar="EXECUTABLE")
    run.add_argument("--r-library", default=None, metavar="DIRECTORY")
    run.set_defaults(handler=_run)

    summarize = subparsers.add_parser(
        "summarize",
        help="Verify and summarize a complete public edgeR result bundle.",
    )
    summarize.command_name = "summarize"
    summarize.add_argument("--run-dir", required=True, metavar="DIRECTORY")
    summarize.set_defaults(handler=_summarize)

    return parser


def _capabilities(_arguments: argparse.Namespace) -> dict[str, Any]:
    from .analysis_contract import ANALYSIS_REQUEST_SCHEMA_VERSIONS

    return {
        "toolkit": {
            "name": "rnaseq-downstream",
            "version": __version__,
            "maturity": "research_preview",
        },
        "contract": {
            "schema_version": "1.0",
            "stdout": "one_json_document",
            "non_interactive": True,
            "exit_codes": {
                "success": int(ExitCode.SUCCESS),
                "request_error": int(ExitCode.REQUEST_ERROR),
                "scientific_validation_error": int(
                    ExitCode.SCIENTIFIC_VALIDATION_ERROR
                ),
                "backend_error": int(ExitCode.BACKEND_ERROR),
                "partial": int(ExitCode.PARTIAL),
                "internal_error": int(ExitCode.INTERNAL_ERROR),
            },
        },
        "commands": {
            "capabilities": "implemented",
            "inspect": "implemented_input_semantics_read_only",
            "validate": "implemented_input_semantics_only",
            "run": "implemented_edger_ql_p0",
            "summarize": "implemented_verified_result_bundle",
        },
        "validated_input_semantics": [
            "featurecounts_integer",
            "salmon_quant_dirs_full_length",
            "salmon_quant_dirs_three_prime",
        ],
        "analysis_request_schema_versions": list(ANALYSIS_REQUEST_SCHEMA_VERSIONS),
        "certified_analysis_paths": [],
        "evidence_gated_analysis_paths": [
            {
                "path_id": "edger_ql_p0_v1",
                "maturity": "research_preview",
                "execution_scope": "validated_p0_input",
                "backend": "edgeR_QL",
                "pipeline": [
                    "filterByExpr",
                    "normLibSizes:TMM",
                    "glmQLFit:legacy_false_robust_true",
                    "glmQLFTest_or_glmTreat",
                ],
                "benchmark_evidence": [
                    "airway-edger-ql-same-engine-v1",
                    "compcoder-edger-ql-nb-fdr-tpr-v1",
                ],
                "benchmark_scope": "backend_kernel_only",
                "input_route_evidence": (
                    "validation_contract_plus_locked_integration_tests"
                ),
                "combined_manifest_origin_authentication": ("self_attested_not_proven"),
                "end_to_end_publication_grade_claim": False,
                "publication_grade_claim": False,
            },
            {
                "path_id": "edger_ql_p0_v1_limma_gene_sets_v1",
                "parent_path_id": "edger_ql_p0_v1",
                "maturity": "research_preview",
                "execution_scope": "validated_p0_input",
                "analysis_request_schema_version": "1.1",
                "backend": "edgeR_QL_plus_limma_gene_set_tests",
                "input_sources": "frozen_local_gmt_plus_annotation",
                "self_contained": {
                    "primary": "limma_fry",
                    "corroborative": "limma_mroast",
                },
                "competitive": {"supplementary": "limma_camera"},
                "multiple_testing": (
                    "BH_within_contrast_method_and_hypothesis"
                ),
                "benchmark_evidence": [
                    "airway-limma-gene-set-same-engine-v1",
                    "compcoder-limma-self-contained-fdr-tpr-v1",
                ],
                "benchmark_scope": "backend_kernel_only",
                "online_resources": False,
                "end_to_end_publication_grade_claim": False,
                "publication_grade_claim": False,
            },
        ],
        "non_statistical_display_capabilities": [
            {
                "capability_id": "edger_ql_p0_v1_static_svg_display_v1",
                "maturity": "research_preview",
                "analysis_path_id": "edger_ql_p0_v1",
                "analysis_request_schema_version": "1.1",
                "invocation": "optional_same_run",
                "statistical_role": "display_only_no_inference",
                "output_location": "display/",
                "output_format": "svg",
                "plot_types": {
                    "volcano": "one_per_contrast",
                    "ma": "one_per_contrast",
                    "pca": "one_per_analysis",
                },
                "pca_input": "post_filter_post_tmm_edger_logcpm",
                "pca_scaling": "centered_unscaled",
                "determinism_scope": "locked_runtime",
                "verification": "summarize_source_reproduction",
                "publication_grade_claim": False,
            }
        ],
    }


def _input_only_scope(
    *,
    input_semantics: str,
    full_numeric_validation: str,
) -> dict[str, Any]:
    return {
        "validation_scope": "input_semantics_only",
        "input_semantics": input_semantics,
        "full_numeric_validation": full_numeric_validation,
        "design": "not_run",
        "backend": "not_run",
        "runnable": False,
        "analysis_path_certified": False,
    }


def _inspect(arguments: argparse.Namespace) -> _CommandResult:
    from .input_semantics import inspect_request

    result = inspect_request(arguments.request)
    data = {
        "scope": _input_only_scope(
            input_semantics="inspected",
            full_numeric_validation="not_run",
        ),
        "input": result["data"],
    }
    return _CommandResult(
        data=data,
        warnings=tuple(result["warnings"]),
        artifacts=tuple(result["artifacts"]),
    )


def _validate(arguments: argparse.Namespace) -> _CommandResult:
    from .validation_bundle import validate_request_to_bundle

    completed = validate_request_to_bundle(arguments.request, arguments.output_dir)
    result = completed["validation_result"]
    bundle = completed["bundle"]
    data = {
        "scope": _input_only_scope(
            input_semantics=result["data"]["input_certification_status"],
            full_numeric_validation="passed",
        ),
        "input": result["data"],
        "bundle": bundle,
    }
    return _CommandResult(
        data=data,
        warnings=tuple([*result["warnings"], *bundle["warnings"]]),
        artifacts=tuple([*result["artifacts"], *bundle["artifacts"]]),
    )


def _run(arguments: argparse.Namespace) -> _CommandResult:
    from .edger_backend import run_edger_ql

    completed = run_edger_ql(
        arguments.request,
        arguments.output_dir,
        rscript=arguments.rscript,
        r_library=arguments.r_library,
    )
    if (
        completed.get("status") != "success"
        or completed.get("execution_scope") != "validated_p0_input"
    ):
        raise BackendFailedError(
            "The edgeR adapter returned an incompatible completion record."
        )
    backend_stderr = completed.get("backend_stderr", "")
    if not isinstance(backend_stderr, str):
        raise BackendFailedError("The edgeR adapter returned invalid backend logs.")
    if backend_stderr:
        sys.stderr.write(backend_stderr)
        sys.stderr.flush()
    backend_data = completed.get("data")
    if not isinstance(backend_data, Mapping):
        raise BackendFailedError("The edgeR adapter omitted its completion data.")
    data = {
        "scope": {
            "analysis_path": "edger_ql_p0_v1",
            "execution_scope": "validated_p0_input",
            "input_semantics": "passed",
            "design": "passed",
            "backend": "passed",
            "evidence_status": "evidence_gated_research_preview",
            "benchmark_scope": "backend_kernel_only",
            "input_route_evidence": (
                "validation_contract_plus_locked_integration_tests"
            ),
            "combined_manifest_origin_authentication": "self_attested_not_proven",
            "publication_grade_claim": False,
        },
        "backend": completed["backend"],
        "output_dir": completed["output_dir"],
        "publication_status": completed["publication_status"],
        "plan_id": completed["plan_id"],
        "analysis_request_sha256": completed["analysis_request_sha256"],
        "runtime_identity": backend_data.get("runtime_identity"),
        "input_semantics": backend_data.get("input_semantics"),
        "route_observed": backend_data.get("route_observed"),
        "design": {
            "columns": backend_data.get("design_columns"),
            "rank": backend_data.get("design_rank"),
            "residual_df": backend_data.get("residual_df"),
        },
        "gene_counts": {
            "total": backend_data.get("gene_count"),
            "tested": backend_data.get("tested_gene_count"),
            "filtered": backend_data.get("filtered_gene_count"),
        },
        "ql_fit_parameters": backend_data.get("ql_fit_parameters"),
        "display_export": backend_data.get("display_export"),
        "contrasts": backend_data.get("contrasts"),
    }
    pathway_analysis = backend_data.get("pathway_analysis")
    if pathway_analysis is not None:
        if not isinstance(pathway_analysis, Mapping):
            raise BackendFailedError(
                "The edgeR adapter returned invalid pathway completion data."
            )
        data["scope"]["analysis_path"] = (
            "edger_ql_p0_v1_limma_gene_sets_v1"
        )
        data["pathways"] = dict(pathway_analysis)
    return _CommandResult(
        data=data,
        warnings=tuple(completed.get("warnings", ())),
        artifacts=tuple(completed.get("artifacts", ())),
    )


def _summarize(arguments: argparse.Namespace) -> _CommandResult:
    from .run_summary import summarize_run

    completed = summarize_run(arguments.run_dir)
    artifacts = completed.pop("artifacts")
    return _CommandResult(data=completed, artifacts=tuple(artifacts))


def _infer_command(arguments: Sequence[str]) -> str:
    for argument in arguments:
        if not argument.startswith("-"):
            return argument
    return "cli"


def _success_envelope(
    command: str,
    result: _CommandResult | Any,
) -> dict[str, Any]:
    if isinstance(result, _CommandResult):
        return build_envelope(
            command,
            status=Status.SUCCESS,
            data=result.data,
            warnings=result.warnings,
            artifacts=result.artifacts,
        )
    return build_envelope(command, status=Status.SUCCESS, data=result)


def _error_envelope(command: str, error: ToolkitError) -> dict[str, Any]:
    return build_envelope(
        command,
        status=error.status,
        errors=[error.to_dict()],
    )


def _execute(arguments: Sequence[str]) -> tuple[dict[str, Any], ExitCode]:
    command = _infer_command(arguments)
    try:
        parsed = _build_parser().parse_args(list(arguments))
        handler: Callable[[argparse.Namespace], _CommandResult | Any] = parsed.handler
        result = handler(parsed)
        return _success_envelope(parsed.command, result), ExitCode.SUCCESS
    except _HelpRequested as help_request:
        return (
            _success_envelope(help_request.command, {"help": help_request.text}),
            ExitCode.SUCCESS,
        )
    except ToolkitError as error:
        return _error_envelope(command, error), error.exit_code
    except Exception:
        LOGGER.exception("Unhandled CLI failure")
        error = InternalToolkitError("An unexpected internal error occurred.")
        return _error_envelope(command, error), ExitCode.INTERNAL_ERROR


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI and return a stable process exit code."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    envelope, exit_code = _execute(arguments)
    serialized_safely = write_json_document(envelope)
    if not serialized_safely:
        return int(ExitCode.INTERNAL_ERROR)
    return int(exit_code)


__all__ = ["main", "ErrorCode", "ExitCode"]
