"""Dependency-free benchmark report schema regression tests."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
from types import ModuleType

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
COMMON_PATH = PROJECT_ROOT / "scripts/benchmark/common.py"
ARCHIVED_REPORT = Path(__file__).with_name("airway-benchmark-report.json")
BENCHMARK_ID = "airway-edger-ql-same-engine-v1"


def _load_common() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "rnaseq_benchmark_common_for_tests", COMMON_PATH
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


COMMON = _load_common()


def _valid_report() -> dict[str, object]:
    value = json.loads(ARCHIVED_REPORT.read_text(encoding="utf-8"))
    assert isinstance(value, dict)
    return value


@pytest.mark.oracle
def test_dependency_free_validator_accepts_archived_report() -> None:
    COMMON.assert_report_shape(_valid_report(), benchmark_id=BENCHMARK_ID)


@pytest.mark.oracle
@pytest.mark.parametrize(
    "mutation",
    [
        lambda report: report.update({"unexpected": True}),
        lambda report: report["implementation"][0].update({"sha256": "A" * 64}),
        lambda report: report["inputs"][0].update({"size_bytes": True}),
        lambda report: report.update({"status": "fail"}),
        lambda report: report["assertions"][0].update({"requirement": ""}),
    ],
    ids=[
        "unexpected-root-key",
        "noncanonical-digest",
        "boolean-size",
        "status-disagrees",
        "blank-requirement",
    ],
)
def test_dependency_free_validator_rejects_schema_or_invariant_violation(
    mutation,
) -> None:
    report = copy.deepcopy(_valid_report())
    mutation(report)
    with pytest.raises(COMMON.BenchmarkError):
        COMMON.assert_report_shape(report, benchmark_id=BENCHMARK_ID)
