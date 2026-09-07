"""Archived checks for the read-only DESeq2 calibration diagnostic."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType

import pytest

from scripts.benchmark.evidence_resolver import verify_archived_implementation

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUNNER = PROJECT_ROOT / "scripts/benchmark/run_deseq2_compcoder_mechanism_diagnostic.py"
R_DIAGNOSTIC = (
    PROJECT_ROOT / "scripts/benchmark/run_deseq2_compcoder_mechanism_diagnostic.R"
)
GENERATOR = PROJECT_ROOT / "scripts/benchmark/generate_deseq2_compcoder_fixture.R"
REPORT = Path(__file__).with_name("deseq2-compcoder-mechanism-diagnostic.json")
MARKDOWN = REPORT.with_suffix(".md")
EXPLORATORY_REPORT = Path(__file__).with_name(
    "deseq2-compcoder-exploratory-report.json"
)
ARTIFACT_DIR = Path(__file__).with_name("deseq2-compcoder-mechanism-artifacts")
IMPLEMENTATION_PATHS = {
    RUNNER.name: RUNNER,
    R_DIAGNOSTIC.name: R_DIAGNOSTIC,
    GENERATOR.name: GENERATOR,
    "common.py": PROJECT_ROOT / "scripts/benchmark/common.py",
}


def _load_runner() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "deseq2_mechanism_diagnostic", RUNNER
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
        document: dict[str, object] = {}
        for key, value in pairs:
            if key in document:
                raise ValueError(f"duplicate JSON key: {key}")
            document[key] = value
        return document

    value = json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=reject_constant,
        object_pairs_hook=unique_pairs,
    )
    assert isinstance(value, dict)
    return value


def test_runner_is_hard_scoped_to_the_disclosed_exploratory_grid() -> None:
    runner = _load_runner()
    assert runner.EXPLORATORY_SEEDS == tuple(range(61001, 61021))
    assert runner.REPLICATES == tuple(range(1, 21))
    assert not hasattr(runner, "CERTIFICATION_SEEDS")
    assert not hasattr(runner, "CERTIFICATION_SEED_BASE")
    assert "620" not in RUNNER.read_text(encoding="utf-8")
    assert "620" not in R_DIAGNOSTIC.read_text(encoding="utf-8")


def test_archived_diagnostic_is_strict_read_only_hypothesis_evidence() -> None:
    report = _strict_json(REPORT)
    assert set(report) == {
        "schema_version",
        "diagnostic_id",
        "status",
        "evidence_role",
        "scope",
        "runtime",
        "implementation",
        "inputs",
        "artifacts",
        "methods",
        "dispersion_fit",
        "truth_dispersion_stratification",
        "method_comparisons",
        "limitations",
    }
    assert report["schema_version"] == "rnaseq-downstream-diagnostic-report-v1"
    assert report["diagnostic_id"] == (
        "compcoder-deseq2-calibration-mechanism-exploratory-v1"
    )
    assert report["status"] == "complete"
    assert report["evidence_role"] == "hypothesis_generation_only_not_certification"
    scope = report["scope"]
    assert isinstance(scope, dict)
    assert scope == {
        "analysis_is_read_only": True,
        "backend_semantics_changed": False,
        "fixture_changed": False,
        "gate_thresholds_changed": False,
        "held_out_seeds_consumed_by_this_diagnostic": False,
        "seed_grid": list(range(61001, 61021)),
        "seed_role": "disclosed_exploratory_only",
    }
    serialized = REPORT.read_text(encoding="utf-8")
    assert "/tmp/" not in serialized
    assert "/home/" not in serialized


def test_diagnostic_binds_current_implementation_and_final_exploratory_report() -> None:
    report = _strict_json(REPORT)
    verify_archived_implementation(report["implementation"], IMPLEMENTATION_PATHS)
    inputs = {
        item["name"]: item
        for item in report["inputs"]  # type: ignore[union-attr]
    }
    exploratory = inputs["deseq2-compcoder-exploratory-report.json"]
    runner = _load_runner()
    assert exploratory == runner.file_evidence(
        EXPLORATORY_REPORT,
        name="deseq2-compcoder-exploratory-report.json",
    )


def test_regenerated_fixture_evidence_matches_the_archived_exploratory_inputs() -> None:
    diagnostic = _strict_json(REPORT)
    exploratory = _strict_json(EXPLORATORY_REPORT)
    diagnostic_inputs = {
        item["name"]: item
        for item in diagnostic["inputs"]  # type: ignore[union-attr]
    }
    exploratory_inputs = {
        item["name"]: item
        for item in exploratory["inputs"]  # type: ignore[union-attr]
    }
    fixture_names = {
        f"replicate-{replicate:02d}/{filename}"
        for replicate in range(1, 21)
        for filename in ("counts.tsv", "metadata.tsv", "truth.tsv", "fixture.json")
    }
    assert len(fixture_names) == 80
    assert all(
        diagnostic_inputs[name] == exploratory_inputs[name] for name in fixture_names
    )
    declared_names = set(diagnostic_inputs)
    assert len(declared_names) == 84
    assert not any(
        "held-out" in name.lower() or "certification" in name.lower() or "620" in name
        for name in declared_names
    )


def test_dispersion_strata_use_unambiguous_denominators() -> None:
    report = _strict_json(REPORT)
    strata = report["truth_dispersion_stratification"]
    assert isinstance(strata, dict)
    assert strata["definition"] == (
        "within_seed_equal_count_rank_quintiles_of_true_dispersion_"
        "with_gene_id_tie_break"
    )
    assert strata["denominators"] == {
        "false_positive_share": "false_positives_in_stratum/all_false_positives",
        "null_false_positive_rate": (
            "false_positives_in_stratum/published_tested_truth_null_"
            "genes_after_independent_filtering_in_stratum"
        ),
        "pooled_fdp_contribution": ("false_positives_in_stratum/all_discoveries"),
        "within_stratum_fdp": ("false_positives_in_stratum/discoveries_in_stratum"),
    }
    combined = strata["combined_by_stratum"]
    assert isinstance(combined, list)
    assert [item["dispersion_stratum"] for item in combined] == [
        "Q1_lowest",
        "Q2_low",
        "Q3_middle",
        "Q4_high",
        "Q5_highest",
    ]
    assert [item["false_positives"] for item in combined] == [72, 108, 158, 197, 340]
    assert sum(item["false_positive_share"] for item in combined) == 1.0
    assert all("null_false_positive_rate" not in item for item in combined)
    truth_rows = strata["truth_by_stratum"]
    assert isinstance(truth_rows, list)
    null_rows = [item for item in truth_rows if item["truth_class"] == "null"]
    assert all(item["null_false_positive_rate"] is not None for item in null_rows)


def test_diagnostic_covers_fit_filter_test_and_aligned_method_questions() -> None:
    report = _strict_json(REPORT)
    fit = report["dispersion_fit"]
    assert isinstance(fit, dict)
    assert fit["requested_fit_type"] == "parametric"
    assert fit["resolved_fit_type_counts"] == {
        "lrt": {"parametric": 20},
        "wald": {"parametric": 20},
    }
    assert fit["automatic_fallback_count"] == 0
    assert fit["fallback_impact"]["comparison_status"] == (  # type: ignore[index]
        "not_estimable_no_automatic_fallback_observed"
    )
    tail = fit["extreme_true_dispersion_tail"]
    assert tail["deseq2_max_dispersion_bound"] == 12.0  # type: ignore[index]
    assert tail["truth_null_genes"] == 185  # type: ignore[index]
    assert tail["published_tested_truth_null_genes_after_independent_filtering"] == 81  # type: ignore[index]
    assert tail["false_positives"] == 8  # type: ignore[index]

    comparisons = report["method_comparisons"]
    assert isinstance(comparisons, dict)
    summaries = comparisons["summaries"]
    assert isinstance(summaries, dict)
    expected = {
        "deseq2_wald_independent_filtering_on",
        "deseq2_wald_independent_filtering_off",
        "deseq2_lrt_independent_filtering_on",
        "edger_ql_native",
        "deseq2_wald_aligned_common_tested_family",
        "edger_ql_aligned_common_tested_family",
    }
    assert set(summaries) == expected
    assert set(comparisons) == {
        "summaries",
        "archived_exploratory_wald_parity",
        "wald_lrt_call_overlap",
        "independent_filtering_on_minus_off",
        "wald_minus_lrt",
        "native_deseq2_minus_edger",
        "aligned_deseq2_minus_edger",
    }
    parity = comparisons["archived_exploratory_wald_parity"]
    assert parity["all_match"] is True  # type: ignore[index]
    assert parity["matched_replicates"] == 20  # type: ignore[index]
    assert all(  # type: ignore[union-attr]
        item["all_fields_match"]
        for item in parity["per_seed"]  # type: ignore[index]
    )
    assert (
        summaries["deseq2_wald_independent_filtering_on"][  # type: ignore[index]
            "mean_replicate_fdp"
        ]
        > summaries["edger_ql_native"]["mean_replicate_fdp"]
    )  # type: ignore[index]
    assert (
        summaries["deseq2_wald_aligned_common_tested_family"][  # type: ignore[index]
            "mean_replicate_fdp"
        ]
        > summaries["edger_ql_aligned_common_tested_family"][  # type: ignore[index]
            "mean_replicate_fdp"
        ]
    )
    wald_on = summaries["deseq2_wald_independent_filtering_on"]
    wald_off = summaries["deseq2_wald_independent_filtering_off"]
    lrt = summaries["deseq2_lrt_independent_filtering_on"]
    deseq2_aligned = summaries["deseq2_wald_aligned_common_tested_family"]
    edger_aligned = summaries["edger_ql_aligned_common_tested_family"]
    assert wald_on["mean_replicate_fdp"] == pytest.approx(  # type: ignore[index]
        0.1182056911947614
    )
    assert wald_off["mean_replicate_fdp"] == pytest.approx(  # type: ignore[index]
        0.11796783157302045
    )
    assert lrt["mean_replicate_fdp"] == pytest.approx(  # type: ignore[index]
        0.1062335913561622
    )
    assert lrt["mean_all_truth_de_tpr"] == pytest.approx(0.6385)  # type: ignore[index]
    assert deseq2_aligned["mean_replicate_fdp"] == pytest.approx(  # type: ignore[index]
        0.11171804504704605
    )
    assert edger_aligned["mean_replicate_fdp"] == pytest.approx(  # type: ignore[index]
        0.04833102013811631
    )
    assert deseq2_aligned[  # type: ignore[index]
        "mean_conditional_tpr_within_tested_family"
    ] == pytest.approx(0.7182067533376091)
    assert edger_aligned[  # type: ignore[index]
        "mean_conditional_tpr_within_tested_family"
    ] == pytest.approx(0.6417175954659752)
    assert wald_on["total_discoveries"] == 7380  # type: ignore[index]
    assert wald_on["total_true_positives"] == 6505  # type: ignore[index]
    assert wald_on["total_false_positives"] == 875  # type: ignore[index]
    assert lrt["total_discoveries"] == 7146  # type: ignore[index]
    assert lrt["total_true_positives"] == 6385  # type: ignore[index]
    assert lrt["total_false_positives"] == 761  # type: ignore[index]
    assert deseq2_aligned["mean_tested_genes"] == pytest.approx(4463.35)  # type: ignore[index]
    assert edger_aligned["mean_tested_genes"] == pytest.approx(4463.35)  # type: ignore[index]
    overlap = comparisons["wald_lrt_call_overlap"][  # type: ignore[index]
        "pooled_descriptive_only"
    ]
    overall_overlap = next(  # type: ignore[arg-type]
        item for item in overlap if item["truth_scope"] == "all"
    )
    assert overall_overlap == {
        "both_called": 7144,
        "call_jaccard": pytest.approx(0.9677594147927391),
        "genes": 100000,
        "lrt_only": 2,
        "neither_called": 92618,
        "truth_scope": "all",
        "wald_only": 236,
    }


def test_diagnostic_artifacts_and_human_report_are_bound_and_present() -> None:
    report = _strict_json(REPORT)
    artifacts = {
        item["name"]: item
        for item in report["artifacts"]  # type: ignore[union-attr]
    }
    assert set(artifacts) == {
        "deseq2-compcoder-mechanism-artifacts/false-positive-by-true-dispersion.svg",
        "deseq2-compcoder-mechanism-artifacts/method-fdp-tpr-comparison.svg",
        "deseq2-compcoder-mechanism-artifacts/dispersion-strata-pooled.tsv",
        "deseq2-compcoder-mechanism-artifacts/dispersion-strata-by-seed.tsv",
        "deseq2-compcoder-mechanism-artifacts/method-summary.tsv",
        "deseq2-compcoder-mechanism-artifacts/method-summary-by-seed.tsv",
        "deseq2-compcoder-mechanism-artifacts/wald-lrt-call-overlap.tsv",
        "deseq2-compcoder-mechanism-artifacts/dispersion-fit-by-seed.tsv",
    }
    runner = _load_runner()
    for name, expected in artifacts.items():
        path = ARTIFACT_DIR / Path(name).name
        assert expected == runner.file_evidence(path, name=name)
        if path.suffix == ".svg":
            assert path.read_text(encoding="utf-8").startswith(
                '<svg xmlns="http://www.w3.org/2000/svg"'
            )
    markdown = MARKDOWN.read_text(encoding="utf-8")
    assert markdown == runner._render_markdown(report)
    assert "authoritative machine record" in markdown
    assert "selection-biased" in markdown
    assert "hypothesis-generating diagnostic" in markdown
    assert "cannot be cryptographically proved" in markdown
    assert "FP / published-tested truth-null genes" in markdown
    assert "not mislabeled as a within-stratum FDP" in markdown
    assert "aligned edgeR comparison" in markdown


def test_equal_count_strata_use_gene_id_to_break_true_dispersion_ties() -> None:
    runner = _load_runner()
    gene_ids = ["g09", "g01", "g08", "g02", "g07", "g03", "g06", "g04", "g05", "g00"]
    labels, boundaries = runner._dispersion_strata([1.0] * 10, gene_ids)
    by_gene = dict(zip(gene_ids, labels, strict=True))
    assert [by_gene[f"g{index:02d}"] for index in range(10)] == [
        "Q1_lowest",
        "Q1_lowest",
        "Q2_low",
        "Q2_low",
        "Q3_middle",
        "Q3_middle",
        "Q4_high",
        "Q4_high",
        "Q5_highest",
        "Q5_highest",
    ]
    assert [item["lower_rank"] for item in boundaries] == [2, 4, 6, 8]


def test_aligned_family_recomputes_two_independent_bh_vectors() -> None:
    runner = _load_runner()
    rows = [
        {"gene_id": "g1", "wald_pvalue": 0.01, "edger_pvalue": 0.4},
        {"gene_id": "g2", "wald_pvalue": 0.04, "edger_pvalue": 0.001},
        {"gene_id": "wald_only", "wald_pvalue": 0.02, "edger_pvalue": None},
        {"gene_id": "edge_only", "wald_pvalue": None, "edger_pvalue": 0.02},
    ]
    assert runner._add_aligned_bh(rows) == 2
    by_gene = {row["gene_id"]: row for row in rows}
    assert by_gene["g1"]["wald_fdr_aligned"] == 0.02
    assert by_gene["g2"]["wald_fdr_aligned"] == 0.04
    assert by_gene["g1"]["edger_fdr_aligned"] == 0.4
    assert by_gene["g2"]["edger_fdr_aligned"] == 0.002
    assert by_gene["wald_only"]["wald_fdr_aligned"] is None
    assert by_gene["edge_only"]["edger_fdr_aligned"] is None


def test_score_separates_all_truth_and_tested_family_denominators() -> None:
    runner = _load_runner()
    rows = [
        {"truth_class": "DE", "fdr": 0.01, "in_family": True},
        {"truth_class": "DE", "fdr": None, "in_family": False},
        {"truth_class": "null", "fdr": 0.02, "in_family": True},
        {"truth_class": "null", "fdr": 0.8, "in_family": True},
    ]
    score = runner._score(rows, fdr_field="fdr", tested_field="in_family")
    assert score["truth_de_total"] == 2
    assert score["truth_de_in_tested_family"] == 1
    assert score["all_truth_de_tpr"] == 0.5
    assert score["conditional_tpr_within_tested_family"] == 1.0
    assert score["false_discovery_proportion"] == 0.5
