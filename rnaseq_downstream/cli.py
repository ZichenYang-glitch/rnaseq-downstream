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
    ErrorCode,
    ExitCode,
    FeatureNotImplementedError,
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

    reserved_help = {
        "run": "Execute the planned evidence-gated analysis path.",
        "summarize": "Summarize a completed run from machine-readable artifacts.",
    }
    for command, help_text in reserved_help.items():
        command_parser = subparsers.add_parser(command, help=help_text)
        command_parser.command_name = command
        command_parser.set_defaults(handler=_reserved_command)

    return parser


def _capabilities(_arguments: argparse.Namespace) -> dict[str, Any]:
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
            "run": "reserved_not_implemented",
            "summarize": "reserved_not_implemented",
        },
        "validated_input_semantics": [
            "featurecounts_integer",
            "salmon_quant_dirs_full_length",
            "salmon_quant_dirs_three_prime",
        ],
        "certified_analysis_paths": [],
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


def _reserved_command(arguments: argparse.Namespace) -> dict[str, Any]:
    raise FeatureNotImplementedError(
        f"The '{arguments.command}' command is reserved but not implemented in this checkpoint.",
        details={
            "command": arguments.command,
            "checkpoint": "P0_checkpoint_A_work_items_3_4",
        },
    )


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
