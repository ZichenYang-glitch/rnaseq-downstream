from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import subprocess

import pytest

from rnaseq_downstream.edger_backend import (
    _run_edger_ql_benchmark_kernel,
    run_edger_ql,
)
from rnaseq_downstream.errors import (
    BackendFailedError,
    ContrastNotEstimableError,
    DesignRankDeficientError,
    InputIntegrityError,
    QCValidationError,
)
from rnaseq_downstream.input_semantics import inspect_request
from rnaseq_downstream.run_summary import summarize_run
from rnaseq_downstream.validation_bundle import validate_request_to_bundle


_CORE_OUTPUT_NAMES = {
    "analysis.json",
    "backend_manifest.json",
    "coefficients.tsv",
    "design.tsv",
    "results.tsv",
}


def _locked_runtime() -> tuple[str, str]:
    rscript = os.environ.get("RNASEQ_P0_RSCRIPT")
    r_library = os.environ.get("RNASEQ_P0_R_LIBRARY")
    if not rscript or not r_library:
        pytest.skip("set RNASEQ_P0_RSCRIPT and RNASEQ_P0_R_LIBRARY")
    return rscript, r_library


def _write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def _write_json(path: Path, document: object) -> Path:
    return _write(
        path,
        json.dumps(document, allow_nan=False, sort_keys=True) + "\n",
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _design() -> dict[str, object]:
    return {
        "intercept": True,
        "terms": ["condition"],
        "variables": {
            "condition": {
                "type": "factor",
                "levels": ["control", "treated"],
            }
        },
    }


def _contrasts(*, include_treat: bool = False) -> list[dict[str, object]]:
    values: list[dict[str, object]] = [
        {
            "contrast_id": "treated_vs_control",
            "weights": {"conditiontreated": 1},
            "lfc_threshold": 0,
        }
    ]
    if include_treat:
        values.append(
            {
                "contrast_id": "treated_vs_control_treat",
                "weights": {"conditiontreated": 1},
                "lfc_threshold": 0.5,
            }
        )
    return values


def _matrix_and_metadata(root: Path, *, sample_count: int = 6) -> tuple[Path, Path]:
    samples = [f"s{index + 1}" for index in range(sample_count)]
    count_lines = ["\t".join(["gene_id", *samples])]
    for gene_index in range(40):
        control = 40 + (gene_index % 7)
        treated = control * 4 if gene_index < 8 else control
        values = [control] * (sample_count // 2) + [treated] * (
            sample_count - sample_count // 2
        )
        count_lines.append(
            "\t".join([f"gene_{gene_index + 1}", *(str(value) for value in values)])
        )
    metadata_lines = ["sample_id\tcondition"]
    metadata_lines.extend(
        f"{sample}\t{'control' if index < sample_count // 2 else 'treated'}"
        for index, sample in enumerate(samples)
    )
    return (
        _write(root / "counts.tsv", "\n".join(count_lines) + "\n"),
        _write(root / "metadata.tsv", "\n".join(metadata_lines) + "\n"),
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _display_configuration(
    *, pca_components: list[int] | None = None
) -> dict[str, object]:
    return {
        "fdr_threshold": 0.05,
        "pca_top_n": 40,
        "pca_components": pca_components or [1, 2],
    }


def _assert_display_bundle(
    output: Path,
    result: dict[str, object],
    *,
    configuration: dict[str, object],
    sample_count: int = 6,
) -> None:
    analysis = result["analysis"]
    data = result["data"]
    assert isinstance(analysis, dict)
    assert isinstance(data, dict)
    contrast_ids = [item["contrast_id"] for item in analysis["contrasts"]]
    expected_display_names = {
        "display_manifest.json",
        "logcpm.tsv",
        "pca.svg",
        "pca_coordinates.tsv",
        "pca_selected_genes.tsv",
        *{f"volcano--{identifier}.svg" for identifier in contrast_ids},
        *{f"ma--{identifier}.svg" for identifier in contrast_ids},
    }
    assert {item.name for item in output.iterdir()} == _CORE_OUTPUT_NAMES | {"display"}
    display = output / "display"
    assert display.is_dir() and not display.is_symlink()
    assert {item.name for item in display.iterdir()} == expected_display_names

    logcpm_rows = _read_tsv(display / "logcpm.tsv")
    assert len(logcpm_rows) == data["tested_gene_count"]
    assert list(logcpm_rows[0]) == [
        "gene_id",
        *(f"s{index + 1}" for index in range(sample_count)),
    ]
    assert len({row["gene_id"] for row in logcpm_rows}) == len(logcpm_rows)
    assert all(
        math.isfinite(float(row[f"s{index + 1}"]))
        for row in logcpm_rows
        for index in range(sample_count)
    )

    expected_logcpm_method = {
        "method": "edgeR::cpm.DGEList",
        "source": "post_filter_post_TMM_observed_DGEList",
        "purpose": "display_only_not_for_inference",
        "arguments": {
            "normalized.lib.sizes": True,
            "log": True,
            "prior.count": 2,
        },
        "scale": "log2",
        "gene_count": len(logcpm_rows),
        "sample_count": sample_count,
    }
    assert data["display_export"] == expected_logcpm_method

    manifest = json.loads(
        (display / "display_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["schema_version"] == "1.0"
    assert manifest["kind"] == "rnaseq_downstream_display_manifest"
    assert manifest["configuration"] == configuration
    assert manifest["methods"]["logcpm"] == expected_logcpm_method
    source = manifest["source_bundle"]
    assert source["analysis_request_schema_version"] == "1.1"
    assert source["analysis_request_sha256"] == result["analysis_request_sha256"]
    assert source["plan_id"] == result["plan_id"]
    source_members = {item["path"]: item for item in source["members"]}
    assert set(source_members) == _CORE_OUTPUT_NAMES
    for name, record in source_members.items():
        assert record["sha256"] == _sha256(output / name)
        assert record["size_bytes"] == (output / name).stat().st_size
    display_members = {item["path"]: item for item in manifest["members"]}
    assert set(display_members) == expected_display_names - {"display_manifest.json"}
    logcpm_member = display_members["logcpm.tsv"]
    assert logcpm_member["role"] == "display_logcpm"
    assert logcpm_member["sha256"] == _sha256(display / "logcpm.tsv")
    assert logcpm_member["size_bytes"] == (display / "logcpm.tsv").stat().st_size

    pca_selected_rows = _read_tsv(display / "pca_selected_genes.tsv")
    summary = summarize_run(output)
    assert summary["status"] == "verified_complete"
    assert summary["display"] == {
        "schema_version": "1.0",
        "status": "verified_complete",
        "source_bundle_id": source["source_bundle_id"],
        "plot_count": 1 + 2 * len(contrast_ids),
        "plot_types": {
            "pca": 1,
            "volcano": len(contrast_ids),
            "ma": len(contrast_ids),
        },
        "tested_gene_count": len(logcpm_rows),
        "pca_gene_count": len(pca_selected_rows),
        "sample_count": sample_count,
        "configuration": configuration,
    }
    artifact_paths = {
        Path(item["path"]).relative_to(output).as_posix()
        for item in summary["artifacts"]
    }
    assert artifact_paths == _CORE_OUTPUT_NAMES | {
        f"display/{name}" for name in expected_display_names
    }


def test_locked_benchmark_kernel_runs_ql_and_treat(tmp_path: Path) -> None:
    rscript, r_library = _locked_runtime()
    counts, metadata = _matrix_and_metadata(tmp_path)
    output = tmp_path / "kernel-output"

    result = _run_edger_ql_benchmark_kernel(
        counts,
        metadata,
        output,
        design=_design(),
        contrasts=_contrasts(include_treat=True),
        rscript=rscript,
        r_library=r_library,
    )

    assert result["status"] == "success"
    assert result["execution_scope"] == "backend_kernel_only"
    assert result["data"]["runtime_identity"] == {
        "R": "4.6.1",
        "Bioconductor": "3.23",
        "BiocVersion_package": "3.23.1",
        "edgeR": "4.10.0",
        "tximport": "1.40.0",
        "limma": "3.68.0",
    }
    results = _read_tsv(output / "results.tsv")
    assert {row["test_method"] for row in results} == {"glmQLFTest", "glmTreat"}
    assert {row["status"] for row in results} <= {"filtered", "tested"}
    for row in results:
        if row["status"] == "tested":
            assert all(row[column] for column in ("logFC", "logCPM", "PValue", "FDR"))
            if row["test_method"] == "glmQLFTest":
                assert row["statistic"]
                assert row["statistic_type"] == "F"
                assert row["statistic_status"] == "reported"
            else:
                assert row["statistic"] == ""
                assert row["statistic_type"] == "not_reported_by_glmTreat"
                assert row["statistic_status"] == "not_reported"
    coefficients = _read_tsv(output / "coefficients.tsv")
    assert list(coefficients[0]) == [
        "gene_id",
        "status",
        "coefficient",
        "estimate",
        "scale",
    ]
    assert all(
        row["estimate"] and row["scale"] == "natural_log"
        for row in coefficients
        if row["status"] == "tested"
    )
    assert {item.name for item in output.iterdir()} == _CORE_OUTPUT_NAMES
    assert result["data"]["display_export"] is None


def test_locked_design_and_contrast_failures_publish_nothing(tmp_path: Path) -> None:
    rscript, r_library = _locked_runtime()
    counts, metadata = _matrix_and_metadata(tmp_path / "rank")
    rank_output = tmp_path / "rank-output"
    rank_design = {
        "intercept": True,
        "terms": ["condition", "duplicate"],
        "variables": {
            "condition": {
                "type": "factor",
                "levels": ["control", "treated"],
            },
            "duplicate": {
                "type": "factor",
                "levels": ["control", "treated"],
            },
        },
    }
    text = metadata.read_text(encoding="utf-8")
    _write(
        metadata,
        text.replace("sample_id\tcondition", "sample_id\tcondition\tduplicate")
        .replace("\tcontrol\n", "\tcontrol\tcontrol\n")
        .replace("\ttreated\n", "\ttreated\ttreated\n"),
    )
    with pytest.raises(DesignRankDeficientError) as rank_error:
        _run_edger_ql_benchmark_kernel(
            counts,
            metadata,
            rank_output,
            design=rank_design,
            contrasts=_contrasts(),
            rscript=rscript,
            r_library=r_library,
        )
    assert rank_error.value.details["reason"] == (
        "complete_confounding_or_redundant_columns"
    )
    assert not rank_output.exists()

    counts, metadata = _matrix_and_metadata(tmp_path / "contrast")
    contrast_output = tmp_path / "contrast-output"
    with pytest.raises(ContrastNotEstimableError) as contrast_error:
        _run_edger_ql_benchmark_kernel(
            counts,
            metadata,
            contrast_output,
            design=_design(),
            contrasts=[
                {
                    "contrast_id": "unknown",
                    "weights": {"missing_coefficient": 1},
                    "lfc_threshold": 0,
                }
            ],
            rscript=rscript,
            r_library=r_library,
        )
    assert contrast_error.value.details["reason"] == "contrast_coefficients_unknown"
    assert not contrast_output.exists()


def test_locked_residual_df_boundary_is_structured(tmp_path: Path) -> None:
    rscript, r_library = _locked_runtime()
    counts, metadata = _matrix_and_metadata(tmp_path, sample_count=2)
    output = tmp_path / "residual-output"
    with pytest.raises(DesignRankDeficientError) as caught:
        _run_edger_ql_benchmark_kernel(
            counts,
            metadata,
            output,
            design=_design(),
            contrasts=_contrasts(),
            rscript=rscript,
            r_library=r_library,
        )
    assert caught.value.details["reason"] == "residual_df_nonpositive"
    assert not output.exists()


def test_locked_duplicate_design_column_names_are_rejected(tmp_path: Path) -> None:
    rscript, r_library = _locked_runtime()
    counts, metadata = _matrix_and_metadata(tmp_path, sample_count=8)
    metadata_lines = ["sample_id\ta\tab"]
    combinations = [
        ("x", "x"),
        ("bc", "x"),
        ("x", "c"),
        ("bc", "c"),
    ]
    metadata_lines.extend(
        f"s{index + 1}\t{combinations[index % 4][0]}\t{combinations[index % 4][1]}"
        for index in range(8)
    )
    _write(metadata, "\n".join(metadata_lines) + "\n")
    output = tmp_path / "duplicate-column-output"
    design = {
        "intercept": True,
        "terms": ["a", "ab"],
        "variables": {
            "a": {"type": "factor", "levels": ["x", "bc"]},
            "ab": {"type": "factor", "levels": ["x", "c"]},
        },
    }

    with pytest.raises(DesignRankDeficientError) as caught:
        _run_edger_ql_benchmark_kernel(
            counts,
            metadata,
            output,
            design=design,
            contrasts=[
                {
                    "contrast_id": "ambiguous",
                    "weights": {"abc": 1},
                    "lfc_threshold": 0,
                }
            ],
            rscript=rscript,
            r_library=r_library,
        )

    assert caught.value.details["reason"] == "design_column_names_duplicated"
    assert caught.value.details["duplicated_columns"] == "abc"
    assert not output.exists()


def _common_request(
    root: Path,
    *,
    semantics: str,
    producer_name: str,
    producer_version: str,
) -> tuple[dict[str, object], list[str]]:
    samples = [f"s{index + 1}" for index in range(6)]
    _write(root / "reference.fa", ">synthetic\nACGT\n")
    _write(
        root / "metadata.tsv",
        "sample_id\tcondition\n"
        + "\n".join(
            f"{sample}\t{'control' if index < 3 else 'treated'}"
            for index, sample in enumerate(samples)
        )
        + "\n",
    )
    return (
        {
            "schema_version": "1.0",
            "input_semantics": semantics,
            "metadata": {
                "path": "metadata.tsv",
                "sample_id_column": "sample_id",
            },
            "producer": {"name": producer_name, "version": producer_version},
            "reference": {
                "name": "synthetic",
                "version": "1",
                "source": "reference.fa",
            },
            "gene_id": {"strip_version": False},
            "samples": [{"sample_id": sample} for sample in samples],
        },
        samples,
    )


def _analysis_request(
    root: Path,
    bundle: Path,
    *,
    display: dict[str, object] | None = None,
) -> Path:
    request: dict[str, object] = {
        "schema_version": "1.1" if display is not None else "1.0",
        "validated_input_bundle": str(bundle),
        "design": _design(),
        "contrasts": _contrasts(),
    }
    if display is not None:
        request["display"] = display
    return _write_json(root / "analysis-request.json", request)


def _validated_featurecounts_bundle(root: Path, *, layout: str) -> Path:
    request, samples = _common_request(
        root,
        semantics="featurecounts_integer",
        producer_name="featureCounts",
        producer_version="2.0.6",
    )
    if layout == "combined_matrix":
        matrix_lines = ["\t".join(["gene_id", *samples])]
        for gene_index in range(40):
            base = 40 + gene_index % 5
            values = [base] * 3 + [base * (4 if gene_index < 8 else 1)] * 3
            matrix_lines.append(
                "\t".join([f"gene_{gene_index + 1}", *(str(value) for value in values)])
            )
        matrix = _write(root / "counts.tsv", "\n".join(matrix_lines) + "\n")
        manifest = {
            "schema_version": "1.0",
            "artifact_type": "featurecounts_integer_matrix",
            "artifact": {"path": "counts.tsv", "sha256": _sha256(matrix)},
            "gene_id_column": "gene_id",
            "display_columns": [],
            "sample_columns": samples,
            "producer": {"name": "featureCounts", "version": "2.0.6"},
            "reference": {
                "name": "synthetic",
                "version": "1",
                "source": "reference.fa",
                "sha256": _sha256(root / "reference.fa"),
            },
        }
        _write_json(root / "counts.manifest.json", manifest)
        request["featurecounts"] = {
            "layout": layout,
            "matrix": "counts.tsv",
            "manifest": "counts.manifest.json",
        }
    else:
        sample_entries = []
        for sample_index, sample in enumerate(samples):
            count_column = f"{sample}.bam"
            lines = [
                "# Program:featureCounts v2.0.6; Command:synthetic",
                "\t".join(
                    ["Geneid", "Chr", "Start", "End", "Strand", "Length", count_column]
                ),
            ]
            for gene_index in range(40):
                base = 40 + gene_index % 5
                count = base * (4 if sample_index >= 3 and gene_index < 8 else 1)
                lines.append(
                    "\t".join(
                        [
                            f"gene_{gene_index + 1}",
                            "1",
                            str(1 + gene_index * 100),
                            str(100 + gene_index * 100),
                            "+",
                            "100",
                            str(count),
                        ]
                    )
                )
            path = _write(
                root / "featurecounts" / f"{sample}.txt",
                "\n".join(lines) + "\n",
            )
            sample_entries.append(
                {
                    "sample_id": sample,
                    "counts_file": str(path.relative_to(root)),
                    "count_column": count_column,
                }
            )
        request["samples"] = sample_entries
        request["featurecounts"] = {"layout": layout}
    request_path = _write_json(root / "request.json", request)
    bundle = root / "validated"
    validate_request_to_bundle(request_path, bundle)
    return bundle


@pytest.mark.parametrize(
    ("layout", "with_display"),
    [("combined_matrix", True), ("per_sample_files", False)],
)
def test_locked_featurecounts_routes_execute(
    tmp_path: Path, layout: str, with_display: bool
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / layout
    bundle = _validated_featurecounts_bundle(root, layout=layout)
    output = root / "results"
    display = _display_configuration() if with_display else None

    result = run_edger_ql(
        _analysis_request(root, bundle, display=display),
        output,
        rscript=rscript,
        r_library=r_library,
    )

    assert result["status"] == "success"
    assert result["execution_scope"] == "validated_p0_input"
    assert result["data"]["route_observed"] == {
        "constructor": "edgeR::DGEList",
        "count_semantics": "integer",
        "transcript_length_offset": False,
    }
    if display is None:
        assert {item.name for item in output.iterdir()} == _CORE_OUTPUT_NAMES
        assert result["data"]["display_export"] is None
        summary = summarize_run(output)
        assert summary["status"] == "verified_complete"
        assert "display" not in summary
        assert len(summary["artifacts"]) == len(_CORE_OUTPUT_NAMES)
    else:
        _assert_display_bundle(output, result, configuration=display)


def test_locked_v11_display_keeps_de_tables_byte_identical_to_v10(
    tmp_path: Path,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "display-byte-invariance"
    bundle = _validated_featurecounts_bundle(root, layout="combined_matrix")
    core_output = root / "v1.0-results"
    display_output = root / "v1.1-results"

    run_edger_ql(
        _analysis_request(root, bundle),
        core_output,
        rscript=rscript,
        r_library=r_library,
    )
    run_edger_ql(
        _analysis_request(root, bundle, display=_display_configuration()),
        display_output,
        rscript=rscript,
        r_library=r_library,
    )

    for name in ("results.tsv", "coefficients.tsv", "design.tsv"):
        assert (display_output / name).read_bytes() == (core_output / name).read_bytes()


def test_locked_featurecounts_display_failure_publishes_nothing(
    tmp_path: Path,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "display-dimension-failure"
    bundle = _validated_featurecounts_bundle(root, layout="combined_matrix")
    output = root / "must-not-exist"
    display = _display_configuration(pca_components=[1, 6])

    with pytest.raises(QCValidationError) as caught:
        run_edger_ql(
            _analysis_request(root, bundle, display=display),
            output,
            rscript=rscript,
            r_library=r_library,
        )

    assert caught.value.details["reason"] == "insufficient_pca_dimensions"
    assert caught.value.details["maximum_components"] == 5
    assert not output.exists()
    assert not list(root.glob(f".{output.name}.edger-*"))


@pytest.mark.parametrize("failure_kind", ["verifier", "internal"])
def test_locked_staged_display_verification_failure_publishes_nothing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure_kind: str,
) -> None:
    from rnaseq_downstream import run_summary

    rscript, r_library = _locked_runtime()
    root = tmp_path / f"staged-verification-{failure_kind}"
    bundle = _validated_featurecounts_bundle(root, layout="combined_matrix")
    output = root / "must-not-exist"

    def fail_staged_verification(_run_dir: str | Path) -> dict[str, object]:
        if failure_kind == "verifier":
            raise InputIntegrityError(
                "Injected staged verifier failure.", details={"injected": True}
            )
        raise RuntimeError("injected staged verifier internal failure")

    monkeypatch.setattr(run_summary, "summarize_run", fail_staged_verification)

    with pytest.raises(BackendFailedError) as caught:
        run_edger_ql(
            _analysis_request(root, bundle, display=_display_configuration()),
            output,
            rscript=rscript,
            r_library=r_library,
        )

    assert caught.value.details["reason"] == "staged_bundle_verification_failed"
    if failure_kind == "verifier":
        assert caught.value.details["cause_type"] == "InputIntegrityError"
        assert caught.value.details["cause_code"] == "INPUT_INTEGRITY_FAILED"
        assert caught.value.details["cause_details"] == {"injected": True}
    else:
        assert caught.value.details == {
            "reason": "staged_bundle_verification_failed",
            "cause_type": "RuntimeError",
        }
    assert not output.exists()
    assert not list(root.glob(f".{output.name}.edger-*"))


def test_locked_r_backend_rejects_unknown_private_request_field(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from rnaseq_downstream import edger_backend

    rscript, r_library = _locked_runtime()
    root = tmp_path / "private-envelope-defense"
    bundle = _validated_featurecounts_bundle(root, layout="combined_matrix")
    output = root / "must-not-exist"
    original_write = edger_backend._write_private_json

    def write_with_unknown_field(path: Path, document: dict[str, object]) -> None:
        altered = dict(document)
        altered["unexpected_private_field"] = True
        original_write(path, altered)

    monkeypatch.setattr(edger_backend, "_write_private_json", write_with_unknown_field)

    with pytest.raises(BackendFailedError) as caught:
        run_edger_ql(
            _analysis_request(root, bundle),
            output,
            rscript=rscript,
            r_library=r_library,
        )

    assert caught.value.details["reason"] == "backend_request_schema_invalid"
    assert caught.value.details["unexpected_fields"] == "unexpected_private_field"
    assert not output.exists()
    assert not list(root.glob(f".{output.name}.edger-*"))


def _salmon_request(
    root: Path,
    *,
    three_prime: bool,
    inferential_replicates: int,
) -> Path:
    semantics = (
        "salmon_quant_dirs_three_prime"
        if three_prime
        else "salmon_quant_dirs_full_length"
    )
    request, samples = _common_request(
        root,
        semantics=semantics,
        producer_name="Salmon",
        producer_version="1.10.3",
    )
    sample_entries = []
    transcript_count = 40
    for sample_index, sample in enumerate(samples):
        quant_dir = root / "salmon" / sample
        quant_lines = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
        sample_counts: list[float] = []
        for gene_index in range(transcript_count):
            base = 40 + gene_index % 5
            count = float(base * (4 if sample_index >= 3 and gene_index < 8 else 1))
            sample_counts.append(count)
            quant_lines.append(
                f"tx_{gene_index + 1}\t1000\t{800 + (gene_index % 7)}\t1\t{count}"
            )
        _write(quant_dir / "quant.sf", "\n".join(quant_lines) + "\n")
        _write_json(
            quant_dir / "cmd_info.json",
            {
                "salmon_version": "1.10.3",
                "index": "/portable/synthetic-index",
                "numBootstraps": inferential_replicates,
                "numGibbsSamples": 0,
                "index_seq_hash": "a" * 64,
            },
        )
        _write_json(
            quant_dir / "aux_info" / "meta_info.json",
            {
                "salmon_version": "1.10.3",
                "num_bootstraps": inferential_replicates,
                "samp_type": "bootstrap" if inferential_replicates else "none",
                "num_valid_targets": transcript_count,
                "quant_errors": [],
                "index_seq_hash": "a" * 64,
            },
        )
        if inferential_replicates:
            bootstrap_dir = quant_dir / "aux_info" / "bootstrap"
            bootstrap_dir.mkdir(parents=True, exist_ok=True)
            with gzip.open(
                bootstrap_dir / "names.tsv.gz", "wt", encoding="utf-8"
            ) as handle:
                handle.write(
                    "\t".join(f"tx_{index + 1}" for index in range(transcript_count))
                    + "\n"
                )
            values = [
                count * (1 + 0.01 * (replicate - 1))
                for replicate in range(inferential_replicates)
                for count in sample_counts
            ]
            with gzip.open(bootstrap_dir / "bootstraps.gz", "wb") as handle:
                handle.write(struct.pack(f"<{len(values)}d", *values))
        sample_entries.append({"sample_id": sample, "quant_dir": f"salmon/{sample}"})
    tx2gene = _write(
        root / "tx2gene.tsv",
        "transcript_id\tgene_id\n"
        + "\n".join(
            f"tx_{index + 1}\tgene_{index + 1}" for index in range(transcript_count)
        )
        + "\n",
    )
    request["samples"] = sample_entries
    request["salmon"] = {
        "tx2gene": "tx2gene.tsv",
        "tx2gene_sha256": _sha256(tx2gene),
    }
    if three_prime:
        request["assay_protocol"] = "three_prime"
    return _write_json(root / "request.json", request)


@pytest.mark.parametrize(
    ("three_prime", "inferential_replicates"),
    [(False, 0), (False, 2), (True, 0), (True, 1)],
)
def test_locked_salmon_routes_execute(
    tmp_path: Path,
    three_prime: bool,
    inferential_replicates: int,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / ("three-prime" if three_prime else "full-length")
    request_path = _salmon_request(
        root,
        three_prime=three_prime,
        inferential_replicates=inferential_replicates,
    )
    bundle = root / "validated"
    validate_request_to_bundle(request_path, bundle)
    output = root / "results"
    display = _display_configuration()

    result = run_edger_ql(
        _analysis_request(root, bundle, display=display),
        output,
        rscript=rscript,
        r_library=r_library,
    )

    route = result["data"]["route_observed"]
    if three_prime:
        assert route["constructor"] == "edgeR::DGEList"
        assert route["count_source"] == "txi$counts"
        assert route["transcript_length_offset"] is False
        assert route["gene_length_correction"] is False
    else:
        assert route["constructor"] == "edgeR::DGEListFromTximport"
        assert route["transcript_length_offset"] is True
        assert route["divide"] is (inferential_replicates >= 2)
        assert route["inferential_replicates_imported"] is (inferential_replicates > 0)
    _assert_display_bundle(output, result, configuration=display)


def test_locked_r_backend_independently_rejects_one_full_length_replicate(
    tmp_path: Path,
) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / "single-replicate-defense"
    request_path = _salmon_request(
        root,
        three_prime=False,
        inferential_replicates=1,
    )
    inspected = inspect_request(request_path)
    input_data = inspected["data"]
    samples = input_data["sample_order"]
    input_data["metadata_values"] = {
        "sample_id": samples,
        "condition": ["control"] * 3 + ["treated"] * 3,
    }
    backend_request = _write_json(
        root / "normalized-backend-request.json",
        {
            "schema_version": "1.0",
            "kind": "edger_ql_backend_request",
            "execution_scope": "adversarial_test_only",
            "analysis_request": {"path": str(request_path), "sha256": "a" * 64},
            "input_evidence": {"plan_id": "b" * 64},
            "input": input_data,
            "design": _design(),
            "contrasts": _contrasts(),
        },
    )
    output = root / "must-not-exist"
    environment = os.environ.copy()
    environment["R_LIBS_USER"] = r_library
    script = (
        Path(__file__).resolve().parents[2]
        / "rnaseq_downstream"
        / "r_scripts"
        / "edger_ql.R"
    )

    completed = subprocess.run(
        [rscript, "--vanilla", str(script), str(backend_request), str(output)],
        check=False,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="strict",
        env=environment,
    )

    assert completed.returncode == 4
    response = json.loads(completed.stdout)
    assert response["status"] == "error"
    assert response["errors"] == [
        {
            "code": "BACKEND_FAILED",
            "message": (
                "Full-length Salmon inferential overdispersion requires at "
                "least two replicates per sample."
            ),
            "details": {
                "reason": "inferential_replicate_count_below_minimum",
                "observed_replicates_per_sample": 1,
                "minimum_replicates_per_sample": 2,
            },
        }
    ]
    assert not output.exists()
