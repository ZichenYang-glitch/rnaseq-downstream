"""Locked-runtime integration coverage for the D1 DESeq2 backend."""

from __future__ import annotations

import copy
import csv
import json
import math
import os
from pathlib import Path
import shutil

import pytest

from rnaseq_downstream import deseq2_backend
from rnaseq_downstream.deseq2_backend import run_deseq2
from rnaseq_downstream.errors import (
    BackendFailedError,
    CountValuesInvalidError,
    DesignRankDeficientError,
)
from rnaseq_downstream.run_summary import summarize_run
from rnaseq_downstream.validation_bundle import validate_request_to_bundle
from tests.integration.test_edger_backend import (
    _design,
    _salmon_request,
    _sha256,
    _validated_featurecounts_bundle,
    _write_json,
)

pytestmark = pytest.mark.integration

_RESULT_HEADER = [
    "gene_id",
    "contrast_id",
    "status",
    "status_reason",
    "baseMean",
    "logFC",
    "unshrunk_logFC",
    "lfcSE",
    "statistic",
    "statistic_type",
    "statistic_hypothesis",
    "PValue",
    "FDR",
    "fdr_basis",
    "test_method",
    "lfc_threshold",
    "shrinkage_method",
]


def _locked_runtime() -> tuple[str, str]:
    rscript = os.environ.get("RNASEQ_P0_RSCRIPT")
    r_library = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not rscript or not r_library:
        pytest.skip("set RNASEQ_P0_RSCRIPT and RNASEQ_P0_R_LIBRARY")
    return rscript, r_library


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _heterogeneous_replicates(
    gene_index: int,
    *,
    treated: bool,
) -> list[int]:
    """Return a deterministic, dispersion-heterogeneous three-sample block."""

    base = 30 + 3 * gene_index
    patterns = (
        (0.95, 1.0, 1.05),
        (0.65, 1.0, 1.35),
        (0.3, 1.0, 1.9),
        (0.12, 0.8, 2.6),
        (0.05, 1.4, 3.2),
    )
    multiplier = 3.5 if treated and gene_index < 16 else 1.0
    pattern = patterns[gene_index % len(patterns)]
    return [max(1, int(round(base * multiplier * factor))) for factor in pattern]


def _validated_deseq2_featurecounts_bundle(root: Path) -> Path:
    """Build a valid input bundle whose dispersion range supports DESeq2 fitting."""

    bundle = _validated_featurecounts_bundle(root, layout="combined_matrix")
    samples = [f"s{index}" for index in range(1, 7)]
    lines = ["\t".join(["gene_id", *samples])]
    for gene_index in range(80):
        values = [
            *_heterogeneous_replicates(gene_index, treated=False),
            *_heterogeneous_replicates(gene_index, treated=True),
        ]
        lines.append(
            "\t".join([f"gene_{gene_index + 1}", *(str(value) for value in values)])
        )
    counts = root / "counts.tsv"
    counts.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest_path = root / "counts.manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["artifact"]["sha256"] = _sha256(counts)
    _write_json(manifest_path, manifest)
    shutil.rmtree(bundle)
    validate_request_to_bundle(root / "request.json", bundle)
    return bundle


def _request(
    root: Path,
    bundle: Path,
    *,
    test_mode: str = "wald",
    shrinkage: str = "none",
    threshold: float = 0,
) -> Path:
    deseq2: dict[str, object] = {
        "test_mode": test_mode,
        "shrinkage": shrinkage,
    }
    if test_mode == "lrt":
        deseq2["reduced"] = {"intercept": True, "terms": []}
    return _write_json(
        root / f"deseq2-{test_mode}-{shrinkage}.request.json",
        {
            "schema_version": "1.2",
            "backend": "deseq2",
            "validated_input_bundle": str(bundle),
            "design": _design(),
            "contrasts": [
                {
                    "contrast_id": "treated_vs_control",
                    "weights": {"conditiontreated": 1},
                    "lfc_threshold": threshold,
                }
            ],
            "deseq2": deseq2,
        },
    )


@pytest.mark.parametrize("threshold", [0, 0.5])
def test_locked_featurecounts_wald_executes_and_summarizes(
    tmp_path: Path, threshold: float
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / f"wald-{threshold}"
    bundle = _validated_deseq2_featurecounts_bundle(root)
    output = root / "results"

    completed = run_deseq2(
        _request(root, bundle, threshold=threshold),
        output,
        rscript=rscript,
        r_library=r_library,
    )

    assert completed["status"] == "success"
    assert completed["backend"] == "DESeq2"
    assert completed["execution_scope"] == "validated_p1_deseq2_input"
    assert {path.name for path in output.iterdir()} == {
        "analysis.json",
        "backend_manifest.json",
        "coefficients.tsv",
        "design.tsv",
        "results.tsv",
    }
    results = _read_tsv(output / "results.tsv")
    assert list(results[0]) == _RESULT_HEADER
    assert {row["statistic_type"] for row in results} == {"Wald"}
    assert {float(row["lfc_threshold"]) for row in results} == {threshold}
    assert {row["shrinkage_method"] for row in results} == {"none"}
    assert {row["status"] for row in results} <= {
        "filtered",
        "not_tested",
        "failed",
        "tested",
    }
    summary = summarize_run(output)
    assert summary["status"] == "verified_complete"
    assert summary["backend"] == "DESeq2"
    assert summary["test"]["mode"] == "wald"


@pytest.mark.parametrize("threshold", [0, 0.5])
def test_locked_apeglm_changes_only_published_logfc(
    tmp_path: Path,
    threshold: float,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "apeglm"
    bundle = _validated_deseq2_featurecounts_bundle(root)
    output = root / "results"

    completed = run_deseq2(
        _request(root, bundle, shrinkage="apeglm", threshold=threshold),
        output,
        rscript=rscript,
        r_library=r_library,
    )

    rows = _read_tsv(output / "results.tsv")
    assert {row["shrinkage_method"] for row in rows} == {"apeglm"}
    assert {float(row["lfc_threshold"]) for row in rows} == {threshold}
    comparable = [
        row
        for row in rows
        if row["logFC"] and row["unshrunk_logFC"] and row["status"] == "tested"
    ]
    assert comparable
    assert any(
        float(row["logFC"]) != float(row["unshrunk_logFC"]) for row in comparable
    )
    assert all(row["PValue"] and row["FDR"] for row in comparable)
    analysis = json.loads((output / "analysis.json").read_text(encoding="utf-8"))
    shrink_step = analysis["pipeline"][-1]
    assert shrink_step["arguments"] == {"lfcThreshold": 0, "svalue": False}
    assert completed["data"]["contrasts"][0]["shrinkage_nonconvergence_count"] == 0
    assert summarize_run(output)["status"] == "verified_complete"


def test_locked_lrt_is_one_omnibus_reporting_effect(
    tmp_path: Path,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "lrt"
    bundle = _validated_deseq2_featurecounts_bundle(root)
    output = root / "results"

    completed = run_deseq2(
        _request(root, bundle, test_mode="lrt"),
        output,
        rscript=rscript,
        r_library=r_library,
    )

    rows = _read_tsv(output / "results.tsv")
    assert {row["contrast_id"] for row in rows} == {"treated_vs_control"}
    assert {row["statistic_type"] for row in rows} == {"LRT"}
    assert {row["statistic_hypothesis"] for row in rows} == {"full_vs_reduced_omnibus"}
    assert {row["fdr_basis"] for row in rows} == {"omnibus_pvalue_BH"}
    assert {float(row["lfc_threshold"]) for row in rows} == {0}
    effect = completed["data"]["reporting_effect"]
    assert effect == [
        {
            "contrast_id": "treated_vs_control",
            "role": "reported_effect_not_lrt_hypothesis",
            "weights": {"conditiontreated": 1},
            "coefficient_name": None,
        }
    ]
    assert summarize_run(output)["test"]["mode"] == "lrt"


@pytest.mark.parametrize("design_case", ["paired", "multiweight"])
def test_locked_wald_supports_paired_and_vector_designs(
    tmp_path: Path, design_case: str
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / design_case
    bundle = _validated_deseq2_featurecounts_bundle(root)
    metadata = root / "metadata.tsv"
    request_path = _request(root, bundle)
    request = json.loads(request_path.read_text(encoding="utf-8"))
    if design_case == "paired":
        metadata.write_text(
            "sample_id\tcondition\tsubject\n"
            "s1\tcontrol\tA\n"
            "s2\tcontrol\tB\n"
            "s3\tcontrol\tC\n"
            "s4\ttreated\tA\n"
            "s5\ttreated\tB\n"
            "s6\ttreated\tC\n",
            encoding="utf-8",
        )
        request["design"] = {
            "intercept": True,
            "terms": ["subject", "condition"],
            "variables": {
                "subject": {"type": "factor", "levels": ["A", "B", "C"]},
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                },
            },
        }
        request["contrasts"][0]["weights"] = {"conditiontreated": 1}
    else:
        request["design"] = {
            "intercept": False,
            "terms": ["condition"],
            "variables": {
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                }
            },
        }
        request["contrasts"][0]["weights"] = {
            "conditioncontrol": -1,
            "conditiontreated": 1,
        }
    shutil.rmtree(bundle)
    validate_request_to_bundle(root / "request.json", bundle)
    _write_json(request_path, request)
    output = root / "results"

    completed = run_deseq2(
        request_path,
        output,
        rscript=rscript,
        r_library=r_library,
    )

    assert completed["status"] == "success"
    assert (
        completed["data"]["contrasts"][0]["weights"]
        == request["contrasts"][0]["weights"]
    )
    assert summarize_run(output)["status"] == "verified_complete"


@pytest.mark.parametrize("failure", ["confounded", "saturated"])
def test_locked_deseq2_design_lint_fails_closed(
    tmp_path: Path,
    failure: str,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / failure
    bundle = _validated_deseq2_featurecounts_bundle(root)
    metadata = root / "metadata.tsv"
    if failure == "confounded":
        metadata.write_text(
            "sample_id\tcondition\tbatch\n"
            "s1\tcontrol\tA\n"
            "s2\tcontrol\tA\n"
            "s3\tcontrol\tA\n"
            "s4\ttreated\tB\n"
            "s5\ttreated\tB\n"
            "s6\ttreated\tB\n",
            encoding="utf-8",
        )
        design = {
            "intercept": True,
            "terms": ["condition", "batch"],
            "variables": {
                "condition": {
                    "type": "factor",
                    "levels": ["control", "treated"],
                },
                "batch": {"type": "factor", "levels": ["A", "B"]},
            },
        }
        weight = "conditiontreated"
        expected_reason = "complete_confounding_or_redundant_columns"
    else:
        metadata.write_text(
            "sample_id\tcondition\tsubject\n"
            + "\n".join(
                f"s{index}\t{'control' if index <= 3 else 'treated'}\ts{index}"
                for index in range(1, 7)
            )
            + "\n",
            encoding="utf-8",
        )
        design = {
            "intercept": True,
            "terms": ["subject"],
            "variables": {
                "subject": {
                    "type": "factor",
                    "levels": [f"s{index}" for index in range(1, 7)],
                }
            },
        }
        weight = "subjects2"
        expected_reason = "residual_df_nonpositive"
    shutil.rmtree(bundle)
    validate_request_to_bundle(root / "request.json", bundle)
    request_path = _request(root, bundle)
    request_document = json.loads(request_path.read_text(encoding="utf-8"))
    request_document["design"] = design
    request_document["contrasts"][0]["weights"] = {weight: 1}
    _write_json(request_path, request_document)
    output = root / "must-not-exist"

    with pytest.raises(DesignRankDeficientError) as caught:
        run_deseq2(
            request_path,
            output,
            rscript=rscript,
            r_library=r_library,
        )

    assert caught.value.details["reason"] == expected_reason
    assert not output.exists()
    assert not list(root.glob(f".{output.name}.deseq2-*"))


@pytest.mark.parametrize(
    ("three_prime", "inferential_replicates", "expected_mode"),
    [
        (False, 0, "deseq2_constructor"),
        (False, 2, "deseq2_constructor"),
        (True, 0, "explicit_before_matrix_constructor"),
        (True, 2, "explicit_before_matrix_constructor"),
    ],
)
def test_locked_salmon_routes_disclose_rounding_and_unused_infreps(
    tmp_path: Path,
    three_prime: bool,
    inferential_replicates: int,
    expected_mode: str,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / f"salmon-{three_prime}-{inferential_replicates}"
    input_request = _salmon_request(
        root,
        three_prime=three_prime,
        inferential_replicates=inferential_replicates,
    )
    # Make every Salmon count genuinely fractional so the approved conversion
    # audit cannot pass vacuously.
    for sample_index, quant in enumerate(sorted((root / "salmon").glob("*/quant.sf"))):
        lines = quant.read_text(encoding="utf-8").splitlines()
        rewritten = [lines[0]]
        for gene_index, line in enumerate(lines[1:]):
            fields = line.split("\t")
            block = _heterogeneous_replicates(
                gene_index,
                treated=sample_index >= 3,
            )
            fractional_part = 0.5 if gene_index % 3 in (0, 1) else 0.25
            fields[-1] = str(float(block[sample_index % 3]) + fractional_part)
            rewritten.append("\t".join(fields))
        quant.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    bundle = root / "validated"
    validate_request_to_bundle(input_request, bundle)
    output = root / "results"

    completed = run_deseq2(
        _request(root, bundle),
        output,
        rscript=rscript,
        r_library=r_library,
    )

    route = completed["data"]["route_observed"]
    audit = route["rounding_audit"]
    assert audit["mode"] == expected_mode
    assert audit["round_function"] == "base::round"
    assert audit["rule"] == "R_base_round_IEC_60559_ties_to_even"
    assert audit["changed_cell_count"] > 0
    assert audit["max_absolute_delta"] == 0.5
    assert len(audit["source_matrix_sha256"]) == 64
    assert len(audit["rounded_matrix_sha256"]) == 64
    assert len(audit["per_sample"]) == 6
    assert all("total_count_delta" in item for item in audit["per_sample"])
    for sample_index, item in enumerate(audit["per_sample"]):
        expected_delta = 0.0
        for gene_index in range(40):
            block = _heterogeneous_replicates(
                gene_index,
                treated=sample_index >= 3,
            )
            fraction = 0.5 if gene_index % 3 in (0, 1) else 0.25
            before = block[sample_index % 3] + fraction
            expected_delta += round(before) - before
        assert math.isclose(item["total_count_delta"], expected_delta, abs_tol=1e-12)
    assert route["inferential_replicates_used_for_inference"] is False
    assert route["inferential_replicates_imported"] is (inferential_replicates > 0)
    assert route["inferential_replicate_state"] == (
        "all" if inferential_replicates else "none"
    )
    assert route["inferential_replicates_per_sample"] == inferential_replicates
    assert route["inferential_replicate_method"] == (
        "bootstrap" if inferential_replicates else None
    )
    assert (
        "does not consume tximport infReps"
        in route["inferential_replicates_unused_reason"]
    )
    assert route["transcript_length_offset"] is (not three_prime)
    assert route["gene_length_correction"] is (not three_prime)
    assert route["countsFromAbundance"] == "no"
    assert route["dropInfReps"] is False
    if three_prime:
        assert route["constructor"] == "DESeq2::DESeqDataSetFromMatrix"
        assert route["count_source"] == "txi$counts"
        assert audit["source_matrix_sha256"] != audit["rounded_matrix_sha256"]
    else:
        assert route["constructor"] == "DESeq2::DESeqDataSetFromTximport"
    assert summarize_run(output)["status"] == "verified_complete"


def test_locked_r_backend_independently_rejects_one_three_prime_replicate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "three-prime-single-replicate-defense"
    input_request = _salmon_request(
        root,
        three_prime=True,
        inferential_replicates=2,
    )
    bundle = root / "validated"
    validate_request_to_bundle(input_request, bundle)
    output = root / "must-not-exist"
    original_write = deseq2_backend._write_private_json

    def write_with_single_replicate(
        path: Path,
        document: dict[str, object],
    ) -> None:
        altered = copy.deepcopy(document)
        salmon = altered["input"]["salmon"]  # type: ignore[index]
        summary = salmon["inferential_replicates"]  # type: ignore[index]
        summary["replicate_count"] = 1
        for record in summary["per_sample"]:
            record["count"] = 1
        original_write(path, altered)

    monkeypatch.setattr(
        deseq2_backend,
        "_write_private_json",
        write_with_single_replicate,
    )

    with pytest.raises(BackendFailedError) as caught:
        run_deseq2(
            _request(root, bundle),
            output,
            rscript=rscript,
            r_library=r_library,
        )

    assert caught.value.details["reason"] == (
        "inferential_replicate_count_below_minimum"
    )
    assert caught.value.details["input_semantics"] == ("salmon_quant_dirs_three_prime")
    assert caught.value.details["observed_replicates_per_sample"] == 1
    assert caught.value.details["minimum_replicates_per_sample"] == 2
    assert caught.value.details["backend_returncode"] == 4
    assert not output.exists()
    assert not list(root.glob(f".{output.name}.deseq2-*"))


def test_locked_deseq2_rejects_counts_above_r_integer_max(
    tmp_path: Path,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "integer-max"
    bundle = _validated_deseq2_featurecounts_bundle(root)
    # Rebuild validation evidence around an otherwise valid value which exceeds
    # the narrower DESeq2 integer domain but remains exact in the input contract.
    counts = root / "counts.tsv"
    lines = counts.read_text(encoding="utf-8").splitlines()
    first_gene = lines[1].split("\t")
    first_gene[1] = "2147483648"
    lines[1] = "\t".join(first_gene)
    counts.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest_path = root / "counts.manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["artifact"]["sha256"] = _sha256(counts)
    _write_json(manifest_path, manifest)
    # Replace the now-stale bundle rather than ever bypassing input validation.
    shutil.rmtree(bundle)
    validate_request_to_bundle(root / "request.json", bundle)
    output = root / "must-not-exist"

    with pytest.raises(CountValuesInvalidError) as caught:
        run_deseq2(
            _request(root, bundle),
            output,
            rscript=rscript,
            r_library=r_library,
        )

    assert caught.value.details["reason"] == (
        "count_value_out_of_deseq2_integer_domain"
    )
    assert not output.exists()


def test_packaged_deseq2_r_backend_is_readable() -> None:
    from rnaseq_downstream.deseq2_backend import _r_script_path

    script = _r_script_path()
    assert script.is_file()
    source = script.read_text(encoding="utf-8")
    assert "DESeq2::DESeq(" in source
    assert "returnList = TRUE" in source
    assert 'shrinkage_fit$fit$diag[, "conv"]' in source


def test_locked_public_cli_executes_and_verifies_deseq2(
    tmp_path: Path,
    run_module_cli,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "public-cli"
    bundle = _validated_deseq2_featurecounts_bundle(root)
    shutil.rmtree(bundle)

    inspected = run_module_cli(
        "inspect",
        "--request",
        str(root / "request.json"),
        cwd=root,
    )
    assert inspected.returncode == 0
    assert inspected.json()["status"] == "success"
    validated = run_module_cli(
        "validate",
        "--request",
        str(root / "request.json"),
        "--output-dir",
        str(bundle),
        cwd=root,
    )
    assert validated.returncode == 0
    assert validated.json()["status"] == "success"

    request = _request(root, bundle)
    output = root / "results"
    completed = run_module_cli(
        "run",
        "--request",
        str(request),
        "--output-dir",
        str(output),
        "--rscript",
        rscript,
        "--r-library",
        r_library,
        cwd=root,
    )
    assert completed.returncode == 0, (completed.stdout, completed.stderr)
    response = completed.json()
    assert response["status"] == "success"
    assert response["data"]["backend"] == "DESeq2"
    assert response["data"]["scope"]["analysis_path"] == ("deseq2_p1_v1_gate_pending")
    assert response["data"]["scope"]["evidence_status"] == (
        "implementation_complete_gate_pending"
    )
    assert response["data"]["scope"]["publication_grade_claim"] is False
    assert len(response["artifacts"]) == 5

    summarized = run_module_cli(
        "summarize",
        "--run-dir",
        str(output),
        cwd=root,
    )
    assert summarized.returncode == 0, (summarized.stdout, summarized.stderr)
    summary = summarized.json()
    assert summary["status"] == "success"
    assert summary["data"]["status"] == "verified_complete"
    assert summary["data"]["backend"] == "DESeq2"
