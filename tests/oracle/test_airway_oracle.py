"""Locked airway same-engine oracle parity gate."""

from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import subprocess
import sys
from types import ModuleType

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUNNER = PROJECT_ROOT / "scripts/benchmark/run_airway_oracle.py"
ARCHIVED_REPORT = Path(__file__).with_name("airway-benchmark-report.json")
BENCHMARK_ID = "airway-edger-ql-same-engine-v1"
IMPLEMENTATION_PATHS = {
    "run_airway_oracle.py": RUNNER,
    "common.py": PROJECT_ROOT / "scripts/benchmark/common.py",
    "benchmark-report-v1.schema.json": (
        PROJECT_ROOT / "scripts/benchmark/benchmark-report-v1.schema.json"
    ),
    "prepare_airway_fixture.R": (
        PROJECT_ROOT / "scripts/benchmark/prepare_airway_fixture.R"
    ),
    "run_airway_direct_oracle.R": (
        PROJECT_ROOT / "scripts/benchmark/run_airway_direct_oracle.R"
    ),
    "edger_backend.py": PROJECT_ROOT / "rnaseq_downstream/edger_backend.py",
    "analysis_contract.py": PROJECT_ROOT / "rnaseq_downstream/analysis_contract.py",
    "edger_ql.R": PROJECT_ROOT / "rnaseq_downstream/r_scripts/edger_ql.R",
    "pathway_tests.R": (
        PROJECT_ROOT / "rnaseq_downstream/r_scripts/pathway_tests.R"
    ),
    "conda-lock.yml": PROJECT_ROOT / "conda-lock.yml",
    "renv.lock": PROJECT_ROOT / "renv.lock",
    "environment.p0.yml": PROJECT_ROOT / "environment.p0.yml",
    "r-sources.lock": PROJECT_ROOT / "environment/r-sources.lock",
    "verify.R": PROJECT_ROOT / "environment/verify.R",
}
PARITY_ASSERTIONS = {
    "tested_gene_set_parity",
    "coefficient_shape_parity",
    "coefficient_numeric_parity",
    "logfc_numeric_parity",
    "pvalue_numeric_parity",
    "f_statistic_numeric_parity",
    "fdr_numeric_parity",
}
P0_NUMERIC_ARTIFACTS = {
    "direct-oracle/results.tsv": (
        "d83a4dbfe00a33546fd5116369d443c6f8f664e1d66330c7d9488ee4fd0eafb8",
        1_794_936,
    ),
    "direct-oracle/coefficients.tsv": (
        "e8eb3c8e867f17d9729bf49a681af0867782a70ebf3aadcf0e457923b7ddc94f",
        1_837_974,
    ),
    "toolkit/results.tsv": (
        "986d5cdb96ffeba3ccecf1abf8af6381fe13826276625b09b56a043226b3b4a2",
        6_439_547,
    ),
    "toolkit/coefficients.tsv": (
        "91fcb08a47718fb9bd73afcb03b40af7c4476ba1215952858f4989193628a49f",
        16_852_482,
    ),
    "toolkit/design.tsv": (
        "dea28100f11c06058d783f24bc1469ff306dbed8ac5c3c33d1152f40ba08901b",
        980,
    ),
}


def _load_runner() -> ModuleType:
    import importlib.util

    spec = importlib.util.spec_from_file_location("airway_oracle_runner", RUNNER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(RUNNER.parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(RUNNER.parent))
    return module


def _write_comparison_fixture(root: Path) -> tuple[Path, Path]:
    oracle = root / "oracle"
    toolkit = root / "toolkit"
    oracle.mkdir()
    toolkit.mkdir()
    (oracle / "results.tsv").write_text(
        "gene_id\tlogFC\tlogCPM\tF\tPValue\tFDR\n"
        "gene-1\t1.25\t5.5\t12.5\t0.001\t0.002\n"
        "gene-2\t-0.75\t4.5\t6.25\t0.02\t0.02\n",
        encoding="utf-8",
    )
    (oracle / "coefficients.tsv").write_text(
        "gene_id\t(Intercept)\tdextrt\ngene-1\t2.0\t1.25\ngene-2\t1.0\t-0.75\n",
        encoding="utf-8",
    )
    (toolkit / "results.tsv").write_text(
        "gene_id\tcontrast_id\tstatus\tlogFC\tstatistic\tPValue\tFDR\n"
        "gene-1\ttrt_vs_untrt\ttested\t1.25\t12.5\t0.001\t0.002\n"
        "gene-2\ttrt_vs_untrt\ttested\t-0.75\t6.25\t0.02\t0.02\n",
        encoding="utf-8",
    )
    (toolkit / "coefficients.tsv").write_text(
        "gene_id\tstatus\tcoefficient\testimate\tscale\n"
        "gene-1\ttested\t(Intercept)\t2.0\tnatural_log\n"
        "gene-1\ttested\tdextrt\t1.25\tnatural_log\n"
        "gene-2\ttested\t(Intercept)\t1.0\tnatural_log\n"
        "gene-2\ttested\tdextrt\t-0.75\tnatural_log\n",
        encoding="utf-8",
    )
    return oracle, toolkit


def _replace_toolkit_result(toolkit: Path, column: str, value: str) -> None:
    path = toolkit / "results.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    column_index = header.index(column)
    first_row = lines[1].split("\t")
    first_row[column_index] = value
    lines[1] = "\t".join(first_row)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _strict_json(path: Path) -> dict[str, object]:
    def reject_constant(value: str) -> None:
        raise ValueError(f"non-finite JSON constant: {value}")

    def unique_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    value = json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=reject_constant,
        object_pairs_hook=unique_pairs,
    )
    assert isinstance(value, dict)
    return value


@pytest.mark.oracle
def test_archived_airway_report_is_passing_evidence() -> None:
    report = _strict_json(ARCHIVED_REPORT)
    assert report["schema_version"] == "rnaseq-downstream-benchmark-report-v1"
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assert report["metrics"]["scope"] == "backend_kernel"  # type: ignore[index]
    assert all(item["passed"] for item in report["assertions"])  # type: ignore[union-attr]
    assertion_names = {item["name"] for item in report["assertions"]}  # type: ignore[union-attr]
    assert PARITY_ASSERTIONS <= assertion_names
    artifacts = {item["name"]: item for item in report["artifacts"]}  # type: ignore[union-attr]
    assert set(artifacts) == set(P0_NUMERIC_ARTIFACTS)
    for name, (digest, size) in P0_NUMERIC_ARTIFACTS.items():
        assert artifacts[name]["sha256"] == digest
        assert artifacts[name]["size_bytes"] == size
    recorded = {item["name"]: item for item in report["implementation"]}  # type: ignore[union-attr]
    assert set(recorded) == set(IMPLEMENTATION_PATHS)
    for name, path in IMPLEMENTATION_PATHS.items():
        payload = path.read_bytes()
        assert recorded[name]["sha256"] == hashlib.sha256(payload).hexdigest()
        assert recorded[name]["size_bytes"] == len(payload)
    serialized = ARCHIVED_REPORT.read_text(encoding="utf-8")
    assert "/tmp/" not in serialized
    assert "/home/" not in serialized


def test_airway_comparison_reports_result_field_parity_metrics(
    tmp_path: Path,
) -> None:
    oracle, toolkit = _write_comparison_fixture(tmp_path)
    runner = _load_runner()

    assertions, metrics = runner._comparison(oracle, toolkit)
    by_name = {item["name"]: item for item in assertions}

    assert PARITY_ASSERTIONS == set(by_name)
    assert all(item["passed"] for item in assertions)
    for prefix in ("logfc", "pvalue", "f_statistic", "fdr"):
        assert metrics[f"{prefix}_value_count"] == 2
        assert metrics[f"{prefix}_violations"] == 0
        assert metrics[f"{prefix}_max_absolute_difference"] == 0.0
        assert metrics[f"{prefix}_max_relative_difference"] == 0.0


@pytest.mark.parametrize(
    ("column", "assertion_name", "metric_prefix"),
    [
        ("PValue", "pvalue_numeric_parity", "pvalue"),
        ("statistic", "f_statistic_numeric_parity", "f_statistic"),
        ("FDR", "fdr_numeric_parity", "fdr"),
    ],
)
def test_airway_comparison_gates_each_inferential_field_independently(
    tmp_path: Path,
    column: str,
    assertion_name: str,
    metric_prefix: str,
) -> None:
    oracle, toolkit = _write_comparison_fixture(tmp_path)
    _replace_toolkit_result(toolkit, column, "999")
    runner = _load_runner()

    assertions, metrics = runner._comparison(oracle, toolkit)
    by_name = {item["name"]: item for item in assertions}

    assert not by_name[assertion_name]["passed"]
    assert by_name[assertion_name]["observed"]["violations"] == 1
    assert metrics[f"{metric_prefix}_violations"] == 1
    other_inferential_assertions = {
        "pvalue_numeric_parity",
        "f_statistic_numeric_parity",
        "fdr_numeric_parity",
    } - {assertion_name}
    assert all(by_name[name]["passed"] for name in other_inferential_assertions)


@pytest.mark.oracle
def test_airway_edger_ql_oracle_parity_rtol_1e_6(tmp_path: Path) -> None:
    r_library_value = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not r_library_value:
        if os.environ.get("RNASEQ_P0_REQUIRE_BENCHMARKS") == "1":
            pytest.fail(
                "Certification mode requires RNASEQ_P0_R_LIBRARY; skipping is forbidden"
            )
        pytest.skip("RNASEQ_P0_R_LIBRARY is required for the locked oracle gate")
    r_library = Path(r_library_value)
    rscript = Path(sys.executable).with_name("Rscript")
    if not rscript.is_file() or not r_library.is_dir():
        pytest.fail("The declared locked P0 R runtime is unavailable")

    report_directory_value = os.environ.get("RNASEQ_P0_BENCHMARK_REPORT_DIR")
    report_directory = (
        Path(report_directory_value) if report_directory_value else tmp_path
    )
    if report_directory.is_symlink():
        pytest.fail("RNASEQ_P0_BENCHMARK_REPORT_DIR must not be a symbolic link")
    report_directory.mkdir(parents=True, exist_ok=True)
    if report_directory.is_symlink() or not report_directory.is_dir():
        pytest.fail("RNASEQ_P0_BENCHMARK_REPORT_DIR is not a directory")
    report_path = report_directory / "airway-benchmark-report.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(RUNNER),
            "--rscript",
            str(rscript),
            "--r-library",
            str(r_library),
            "--report",
            str(report_path),
        ],
        cwd=PROJECT_ROOT,
        check=False,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert report_path.is_file(), (completed.stdout, completed.stderr)
    report = _strict_json(report_path)
    assert completed.returncode == 0, report
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assert all(item["passed"] for item in report["assertions"])  # type: ignore[union-attr]
    assertion_names = {item["name"] for item in report["assertions"]}  # type: ignore[union-attr]
    assert PARITY_ASSERTIONS <= assertion_names
