"""Regression tests proving legacy library loaders never terminate callers."""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
from typing import Any

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_MODULE = PROJECT_ROOT / "modules" / "data.py"
DATA_SOURCE = DATA_MODULE.read_text(encoding="utf-8")
DATA_TREE = ast.parse(DATA_SOURCE, filename=str(DATA_MODULE))
LOADER_NAMES = {"load_metadata", "load_counts"}


def _loader_nodes() -> dict[str, ast.FunctionDef]:
    return {
        node.name: node
        for node in DATA_TREE.body
        if isinstance(node, ast.FunctionDef) and node.name in LOADER_NAMES
    }


@pytest.mark.unit
def test_legacy_data_module_has_no_process_exit_calls() -> None:
    exit_calls = [
        node
        for node in ast.walk(DATA_TREE)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "sys"
        and node.func.attr == "exit"
    ]

    assert exit_calls == []


@pytest.mark.unit
def test_loader_failure_paths_are_narrow_and_never_print() -> None:
    loaders = _loader_nodes()
    assert set(loaders) == LOADER_NAMES

    for name, function in loaders.items():
        print_calls = [
            node
            for node in ast.walk(function)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "print"
        ]
        broad_handlers = [
            node
            for node in ast.walk(function)
            if isinstance(node, ast.ExceptHandler)
            and isinstance(node.type, ast.Name)
            and node.type.id == "Exception"
        ]

        assert print_calls == [], f"{name} must not print library errors"
        assert broad_handlers == [], f"{name} must not relabel semantic failures as I/O"


@pytest.mark.unit
def test_data_module_imports_the_public_input_error_type() -> None:
    imports = [
        node
        for node in DATA_TREE.body
        if isinstance(node, ast.ImportFrom)
        and node.module == "rnaseq_downstream.errors"
    ]

    assert any(
        alias.name == "InputReadError"
        for statement in imports
        for alias in statement.names
    )


_HAS_LEGACY_RUNTIME = all(
    importlib.util.find_spec(dependency) is not None
    for dependency in ("pandas", "yaml")
)

if _HAS_LEGACY_RUNTIME:
    from modules.data import load_counts, load_metadata
    from rnaseq_downstream.errors import ErrorCode, InputReadError

    @pytest.mark.unit
    @pytest.mark.parametrize(
        ("loader", "kwargs", "operation"),
        [
            (load_metadata, {"design_col": "condition"}, "load_metadata"),
            (
                load_counts,
                {"metadata_samples": ["sample-a"]},
                "load_counts",
            ),
        ],
    )
    def test_loader_failure_is_typed_and_silent_when_legacy_deps_are_available(
        loader: Any,
        kwargs: dict[str, Any],
        operation: str,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        missing = tmp_path / "missing.tsv"

        with pytest.raises(InputReadError) as captured:
            loader(missing, **kwargs)

        error = captured.value
        assert error.code is ErrorCode.INPUT_READ_FAILED
        assert error.details["path"] == str(missing)
        assert error.details["operation"] == operation
        assert error.cause is not None
        assert capsys.readouterr() == ("", "")
