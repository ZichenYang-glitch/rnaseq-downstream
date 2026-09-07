"""Locked airway same-engine oracle coverage for DESeq2 Wald and LRT."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import ModuleType

import pytest

from scripts.benchmark.evidence_resolver import verify_archived_implementation

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUNNER = PROJECT_ROOT / "scripts/benchmark/run_deseq2_airway_oracle.py"
ARCHIVED_REPORT = Path(__file__).with_name("deseq2-airway-benchmark-report.json")
BENCHMARK_ID = "airway-deseq2-wald-lrt-same-engine-v1"
IMPLEMENTATION_PATHS = {
    "run_deseq2_airway_oracle.py": RUNNER,
    "run_deseq2_airway_direct_oracle.R": (
        PROJECT_ROOT / "scripts/benchmark/run_deseq2_airway_direct_oracle.R"
    ),
    "prepare_airway_fixture.R": (
        PROJECT_ROOT / "scripts/benchmark/prepare_airway_fixture.R"
    ),
    "common.py": PROJECT_ROOT / "scripts/benchmark/common.py",
    "benchmark-report-v1.schema.json": (
        PROJECT_ROOT / "scripts/benchmark/benchmark-report-v1.schema.json"
    ),
    "input_semantics.py": PROJECT_ROOT / "rnaseq_downstream/input_semantics.py",
    "validation_bundle.py": PROJECT_ROOT / "rnaseq_downstream/validation_bundle.py",
    "analysis_contract_v12.py": (
        PROJECT_ROOT / "rnaseq_downstream/analysis_contract_v12.py"
    ),
    "deseq2_backend.py": PROJECT_ROOT / "rnaseq_downstream/deseq2_backend.py",
    "deseq2.R": PROJECT_ROOT / "rnaseq_downstream/r_scripts/deseq2.R",
    "run_summary.py": PROJECT_ROOT / "rnaseq_downstream/run_summary.py",
    "cli.py": PROJECT_ROOT / "rnaseq_downstream/cli.py",
    "conda-lock.yml": PROJECT_ROOT / "conda-lock.yml",
    "renv.lock": PROJECT_ROOT / "renv.lock",
    "environment.p0.yml": PROJECT_ROOT / "environment.p0.yml",
    "r-sources.lock": PROJECT_ROOT / "environment/r-sources.lock",
    "verify.R": PROJECT_ROOT / "environment/verify.R",
}
PARITY_ASSERTIONS = {
    f"{mode}_{field}_numeric_parity"
    for mode in ("wald", "lrt")
    for field in ("coefficient", "logfc", "statistic", "pvalue", "fdr")
}
GENE_SET_ASSERTIONS = {
    "wald_tested_gene_set_parity",
    "lrt_tested_gene_set_parity",
}
EXPECTED_ARTIFACT_NAMES = {
    f"{route}/{mode}/{artifact}"
    for mode in ("wald", "lrt")
    for route, artifacts in (
        ("direct-oracle", ("results.tsv", "coefficients.tsv")),
        ("toolkit", ("results.tsv", "coefficients.tsv", "design.tsv")),
    )
    for artifact in artifacts
}
EXPECTED_ARTIFACTS = {
    "direct-oracle/wald/results.tsv": (
        "0602bdd535f3a09f7cdcb05770af4a22da02a4273713716e85ffe550e880c0a5",
        1_613_467,
    ),
    "direct-oracle/wald/coefficients.tsv": (
        "94b5ed60af152d8630028dabffb9ae7f20ea21ffe23fc47cae406f727df791a1",
        7_529_501,
    ),
    "toolkit/wald/results.tsv": (
        "371c12a2571d96ac9e812831728cac2d590ea35dbf59534a9a10ef66c6480790",
        13_816_688,
    ),
    "toolkit/wald/coefficients.tsv": (
        "c3428f14876abcdcd50476df28d2b37cb0fb585253701761ec066ebacb4d30c2",
        18_683_995,
    ),
    "toolkit/wald/design.tsv": (
        "dea28100f11c06058d783f24bc1469ff306dbed8ac5c3c33d1152f40ba08901b",
        980,
    ),
    "direct-oracle/lrt/results.tsv": (
        "1a17c0f6e071e48d79083f81d1e4334e99e9e97ecf0955dd9c3a9da507dca0b7",
        1_607_257,
    ),
    "direct-oracle/lrt/coefficients.tsv": (
        "94b5ed60af152d8630028dabffb9ae7f20ea21ffe23fc47cae406f727df791a1",
        7_529_501,
    ),
    "toolkit/lrt/results.tsv": (
        "925a94f1fea675d5232b8292731921ef371b9234c5fd51ba41a15c6932fe8dee",
        12_022_626,
    ),
    "toolkit/lrt/coefficients.tsv": (
        "c3428f14876abcdcd50476df28d2b37cb0fb585253701761ec066ebacb4d30c2",
        18_683_995,
    ),
    "toolkit/lrt/design.tsv": (
        "dea28100f11c06058d783f24bc1469ff306dbed8ac5c3c33d1152f40ba08901b",
        980,
    ),
}


def _load_runner() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "deseq2_airway_oracle_runner", RUNNER
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.path.insert(0, str(RUNNER.parent))
    try:
        specification.loader.exec_module(module)
    finally:
        sys.path.remove(str(RUNNER.parent))
    return module


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


def _write_comparison_fixture(root: Path, *, mode: str) -> tuple[Path, Path]:
    oracle = root / "oracle"
    toolkit = root / "toolkit"
    oracle.mkdir()
    toolkit.mkdir()
    (oracle / "results.tsv").write_text(
        "gene_id\tlogFC\tstatistic\tPValue\tFDR\n"
        "gene-1\t1.25\t12.5\t0.001\t0.002\n"
        "gene-2\t-0.75\t6.25\t0.02\t0.02\n",
        encoding="utf-8",
    )
    (oracle / "coefficients.tsv").write_text(
        "gene_id\tcoefficient\testimate\n"
        "gene-1\t(Intercept)\t2.0\n"
        "gene-1\tdextrt\t1.25\n"
        "gene-2\t(Intercept)\t1.0\n"
        "gene-2\tdextrt\t-0.75\n",
        encoding="utf-8",
    )
    statistic_type = "LRT" if mode == "lrt" else "Wald"
    hypothesis = "full_vs_reduced_omnibus" if mode == "lrt" else "contrast_equals_zero"
    fdr_basis = (
        "omnibus_pvalue_BH"
        if mode == "lrt"
        else "contrast_pvalue_BH_after_independent_filtering"
    )
    method = "DESeq2::results_LRT" if mode == "lrt" else "DESeq2::results_Wald"
    header = (
        "gene_id\tcontrast_id\tstatus\tstatus_reason\tbaseMean\tlogFC\t"
        "unshrunk_logFC\tlfcSE\tstatistic\tstatistic_type\t"
        "statistic_hypothesis\tPValue\tFDR\tfdr_basis\ttest_method\t"
        "lfc_threshold\tshrinkage_method\n"
    )
    rows = [
        ("gene-1", "1.25", "12.5", "0.001", "0.002"),
        ("gene-2", "-0.75", "6.25", "0.02", "0.02"),
    ]
    body = "".join(
        "\t".join(
            [
                gene,
                "trt_vs_untrt",
                "tested",
                "tested",
                "10",
                logfc,
                logfc,
                "0.1",
                statistic,
                statistic_type,
                hypothesis,
                pvalue,
                fdr,
                fdr_basis,
                method,
                "0",
                "none",
            ]
        )
        + "\n"
        for gene, logfc, statistic, pvalue, fdr in rows
    )
    (toolkit / "results.tsv").write_text(header + body, encoding="utf-8")
    (toolkit / "coefficients.tsv").write_text(
        "gene_id\tstatus\tstatus_reason\tcoefficient\testimate\tscale\n"
        "gene-1\ttested\tfitted\t(Intercept)\t2.0\tlog2\n"
        "gene-1\ttested\tfitted\tdextrt\t1.25\tlog2\n"
        "gene-2\ttested\tfitted\t(Intercept)\t1.0\tlog2\n"
        "gene-2\ttested\tfitted\tdextrt\t-0.75\tlog2\n",
        encoding="utf-8",
    )
    role = "reported_effect_not_lrt_hypothesis" if mode == "lrt" else "tested_contrast"
    (toolkit / "analysis.json").write_text(
        json.dumps(
            {
                "test": {"mode": mode},
                "reporting_effect": [
                    {
                        "contrast_id": "trt_vs_untrt",
                        "role": role,
                        "weights": {"dextrt": 1},
                    }
                ],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    return oracle, toolkit


def _replace_result(toolkit: Path, column: str, value: str) -> None:
    path = toolkit / "results.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    fields = lines[1].split("\t")
    fields[header.index(column)] = value
    lines[1] = "\t".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.mark.oracle
def test_archived_deseq2_airway_report_is_passing_evidence() -> None:
    report = _strict_json(ARCHIVED_REPORT)
    assert report["schema_version"] == "rnaseq-downstream-benchmark-report-v1"
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assertions = {item["name"]: item for item in report["assertions"]}  # type: ignore[union-attr]
    assert PARITY_ASSERTIONS | GENE_SET_ASSERTIONS <= set(assertions)
    assert all(item["passed"] for item in assertions.values())
    assert assertions["featurecounts_origin_not_certified"]["observed"] == {
        "route_fixture": "typed_combined_integer_matrix",
        "origin_authentication": "self_attested_not_proven",
        "claim_scope": "numerical_parity_and_public_routing_only",
    }
    metrics = report["metrics"]  # type: ignore[assignment]
    assert metrics["scope"] == "public_validate_run_summarize_same_engine_oracle"
    assert metrics["featurecounts_origin_authentication"] == "self_attested_not_proven"
    artifacts = {item["name"]: item for item in report["artifacts"]}  # type: ignore[union-attr]
    assert set(artifacts) == EXPECTED_ARTIFACT_NAMES
    assert set(artifacts) == set(EXPECTED_ARTIFACTS)
    for name, (digest, size) in EXPECTED_ARTIFACTS.items():
        assert artifacts[name]["sha256"] == digest
        assert artifacts[name]["size_bytes"] == size
    verify_archived_implementation(report["implementation"], IMPLEMENTATION_PATHS)
    serialized = ARCHIVED_REPORT.read_text(encoding="utf-8")
    assert "/tmp/" not in serialized
    assert "/home/" not in serialized


@pytest.mark.parametrize("mode", ["wald", "lrt"])
def test_mode_comparison_covers_all_required_numeric_fields(
    tmp_path: Path, mode: str
) -> None:
    oracle, toolkit = _write_comparison_fixture(tmp_path, mode=mode)
    runner = _load_runner()

    assertions, metrics = runner._mode_comparison(
        oracle, toolkit, mode=mode
    )
    by_name = {item["name"]: item for item in assertions}

    expected = {
        f"{mode}_tested_gene_set_parity",
        f"{mode}_coefficient_numeric_parity",
        f"{mode}_result_semantics",
        f"{mode}_logfc_numeric_parity",
        f"{mode}_statistic_numeric_parity",
        f"{mode}_pvalue_numeric_parity",
        f"{mode}_fdr_numeric_parity",
    }
    assert set(by_name) == expected
    assert all(item["passed"] for item in assertions)
    assert metrics["tested_gene_count"] == 2
    assert metrics["coefficient_parity"]["violations"] == 0


@pytest.mark.parametrize("mode", ["wald", "lrt"])
@pytest.mark.parametrize(
    ("column", "field"),
    [
        ("logFC", "logfc"),
        ("statistic", "statistic"),
        ("PValue", "pvalue"),
        ("FDR", "fdr"),
    ],
)
def test_mode_comparison_gates_each_result_field_independently(
    tmp_path: Path, mode: str, column: str, field: str
) -> None:
    oracle, toolkit = _write_comparison_fixture(tmp_path, mode=mode)
    _replace_result(toolkit, column, "999")
    runner = _load_runner()

    assertions, _ = runner._mode_comparison(oracle, toolkit, mode=mode)
    by_name = {item["name"]: item for item in assertions}

    assert not by_name[f"{mode}_{field}_numeric_parity"]["passed"]
    other = {"logfc", "statistic", "pvalue", "fdr"} - {field}
    assert all(by_name[f"{mode}_{name}_numeric_parity"]["passed"] for name in other)


def test_mode_comparison_gates_coefficient_drift(tmp_path: Path) -> None:
    oracle, toolkit = _write_comparison_fixture(tmp_path, mode="wald")
    path = toolkit / "coefficients.tsv"
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            "gene-1\ttested\tfitted\tdextrt\t1.25\tlog2",
            "gene-1\ttested\tfitted\tdextrt\t999\tlog2",
        ),
        encoding="utf-8",
    )
    runner = _load_runner()

    assertions, _ = runner._mode_comparison(oracle, toolkit, mode="wald")
    by_name = {item["name"]: item for item in assertions}

    assert not by_name["wald_coefficient_numeric_parity"]["passed"]


def test_lrt_reporting_effect_mislabel_fails_semantic_gate(tmp_path: Path) -> None:
    oracle, toolkit = _write_comparison_fixture(tmp_path, mode="lrt")
    analysis_path = toolkit / "analysis.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    analysis["reporting_effect"][0]["role"] = "tested_contrast"
    analysis_path.write_text(json.dumps(analysis) + "\n", encoding="utf-8")
    runner = _load_runner()

    assertions, _ = runner._mode_comparison(oracle, toolkit, mode="lrt")
    by_name = {item["name"]: item for item in assertions}

    assert not by_name["lrt_result_semantics"]["passed"]
    assert by_name["lrt_pvalue_numeric_parity"]["passed"]
    assert by_name["lrt_fdr_numeric_parity"]["passed"]


@pytest.mark.oracle
def test_airway_deseq2_wald_lrt_oracle_parity_rtol_1e_6(
    tmp_path: Path,
) -> None:
    r_library_value = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not r_library_value:
        if os.environ.get("RNASEQ_P0_REQUIRE_BENCHMARKS") == "1":
            pytest.fail(
                "Certification mode requires RNASEQ_P0_R_LIBRARY; skipping is forbidden"
            )
        pytest.skip("RNASEQ_P0_R_LIBRARY is required for the locked oracle gate")
    r_library = Path(r_library_value)
    rscript_value = os.environ.get("RNASEQ_P0_RSCRIPT")
    rscript = (
        Path(rscript_value)
        if rscript_value
        else Path(sys.executable).with_name("Rscript")
    )
    if not rscript.is_file() or not r_library.is_dir():
        pytest.fail("The declared locked P1 R runtime is unavailable")

    report_directory_value = os.environ.get("RNASEQ_P0_BENCHMARK_REPORT_DIR")
    report_directory = (
        Path(report_directory_value) if report_directory_value else tmp_path
    )
    if report_directory.is_symlink():
        pytest.fail("RNASEQ_P0_BENCHMARK_REPORT_DIR must not be a symbolic link")
    report_directory.mkdir(parents=True, exist_ok=True)
    if report_directory.is_symlink() or not report_directory.is_dir():
        pytest.fail("RNASEQ_P0_BENCHMARK_REPORT_DIR is not a directory")
    report_path = report_directory / "deseq2-airway-benchmark-report.json"
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
        timeout=1800,
    )
    assert report_path.is_file(), (completed.stdout, completed.stderr)
    report = _strict_json(report_path)
    assert completed.returncode == 0, report
    assert report["benchmark_id"] == BENCHMARK_ID
    assert report["status"] == "pass"
    assertions = {item["name"]: item for item in report["assertions"]}  # type: ignore[union-attr]
    assert PARITY_ASSERTIONS | GENE_SET_ASSERTIONS <= set(assertions)
    assert all(item["passed"] for item in assertions.values())
