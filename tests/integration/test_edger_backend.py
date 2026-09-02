from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import struct

import pytest

from rnaseq_downstream.edger_backend import (
    _run_edger_ql_benchmark_kernel,
    run_edger_ql,
)
from rnaseq_downstream.errors import (
    ContrastNotEstimableError,
    DesignRankDeficientError,
)
from rnaseq_downstream.validation_bundle import validate_request_to_bundle


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


def _analysis_request(root: Path, bundle: Path) -> Path:
    return _write_json(
        root / "analysis-request.json",
        {
            "schema_version": "1.0",
            "validated_input_bundle": str(bundle),
            "design": _design(),
            "contrasts": _contrasts(),
        },
    )


@pytest.mark.parametrize("layout", ["combined_matrix", "per_sample_files"])
def test_locked_featurecounts_routes_execute(tmp_path: Path, layout: str) -> None:
    rscript, r_library = _locked_runtime()
    root = tmp_path / layout
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
    output = root / "results"

    result = run_edger_ql(
        _analysis_request(root, bundle),
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
    ("three_prime", "inferential_replicates"), [(False, 3), (True, 0)]
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

    result = run_edger_ql(
        _analysis_request(root, bundle),
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
        assert route["divide"] is True
        assert route["inferential_replicates_imported"] is True
