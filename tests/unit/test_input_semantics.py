"""Synthetic tests for the strict P0 input-semantics layer."""

from __future__ import annotations

import hashlib
import gzip
import json
from pathlib import Path
import struct
import zlib

import pytest

import rnaseq_downstream.input_semantics as input_semantics_module
from rnaseq_downstream.errors import (
    AssayProtocolRequiredError,
    CountValuesInvalidError,
    GeneIdentifierError,
    InputEvidenceRequiredError,
    InputIntegrityError,
    InputValidationError,
    InvalidRequestError,
    SalmonOffsetRequiredError,
    SampleSetMismatchError,
)
from rnaseq_downstream.input_semantics import inspect_request, validate_request
from rnaseq_downstream.input_semantics import MAX_EXACT_INTEGER_COUNT

pytestmark = pytest.mark.unit


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _write_json(path: Path, value: object) -> Path:
    return _write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _metadata(root: Path, sample_order: tuple[str, ...] = ("s1", "s2")) -> Path:
    rows = ["sample_id\tcondition\tbatch"]
    for index, sample_id in enumerate(sample_order):
        condition = "treated" if sample_id == "s2" else "control"
        rows.append(f"{sample_id}\t{condition}\tb{index + 1}")
    return _write(root / "metadata.tsv", "\n".join(rows) + "\n")


def _reference(root: Path) -> Path:
    return _write(root / "reference.fa", ">tx1\nACGT\n>tx2\nTGCA\n")


def _common_request(
    root: Path,
    *,
    semantics: str,
    producer: str,
    version: str,
    sample_order: tuple[str, ...] = ("s1", "s2"),
) -> dict[str, object]:
    reference = _reference(root)
    _metadata(root, sample_order)
    return {
        "schema_version": "1.0",
        "input_semantics": semantics,
        "metadata": {"path": "metadata.tsv", "sample_id_column": "sample_id"},
        "producer": {"name": producer, "version": version},
        "reference": {
            "name": "synthetic-transcriptome",
            "version": "v1",
            "source": "reference.fa",
            "sha256": _sha256(reference),
        },
        "gene_id": {"strip_version": False},
        "samples": [{"sample_id": sample_id} for sample_id in sample_order],
    }


def _featurecounts_combined(
    root: Path,
    *,
    counts: tuple[tuple[str, str], ...] = (("10", "20"), ("30", "40")),
    gene_ids: tuple[str, ...] = ("ENSG1.1", "ENSG2.1"),
    symbols: tuple[str, ...] = ("DUP", "DUP"),
    sample_order: tuple[str, ...] = ("s1", "s2"),
    metadata_order: tuple[str, ...] | None = None,
    strip_version: bool = False,
) -> tuple[Path, dict[str, object]]:
    request = _common_request(
        root,
        semantics="featurecounts_integer",
        producer="featureCounts",
        version="2.0.6",
        sample_order=sample_order,
    )
    if metadata_order is not None:
        _metadata(root, metadata_order)
    request["gene_id"] = {"strip_version": strip_version}
    matrix_lines = ["gene_id\tgene_name\t" + "\t".join(sample_order)]
    for gene_id, symbol, row_counts in zip(gene_ids, symbols, counts, strict=True):
        matrix_lines.append("\t".join((gene_id, symbol, *row_counts)))
    matrix = _write(root / "counts.tsv", "\n".join(matrix_lines) + "\n")
    reference = root / "reference.fa"
    manifest = {
        "schema_version": "1.0",
        "artifact_type": "featurecounts_integer_matrix",
        "artifact": {"path": "counts.tsv", "sha256": _sha256(matrix)},
        "gene_id_column": "gene_id",
        "display_columns": ["gene_name"],
        "sample_columns": list(sample_order),
        "producer": {"name": "featureCounts", "version": "2.0.6"},
        "reference": {
            "name": "synthetic-transcriptome",
            "version": "v1",
            "source": "reference.fa",
            "sha256": _sha256(reference),
        },
    }
    _write_json(root / "counts.manifest.json", manifest)
    request["featurecounts"] = {
        "layout": "combined_matrix",
        "matrix": "counts.tsv",
        "manifest": "counts.manifest.json",
    }
    request_path = _write_json(root / "request.json", request)
    return request_path, request


def _featurecounts_per_sample(root: Path) -> Path:
    request = _common_request(
        root,
        semantics="featurecounts_integer",
        producer="featureCounts",
        version="2.0.6",
    )
    entries = []
    for sample_id, counts in (("s1", ("10", "30")), ("s2", ("20", "40"))):
        count_column = f"{sample_id}.bam"
        lines = [
            '# Program:featureCounts v2.0.6; Command:"featureCounts" ...',
            "\t".join((*_FC_ANNOTATION_COLUMNS, count_column)),
            f"ENSG1\tchr1\t1\t10\t+\t10\t{counts[0]}",
            f"ENSG2\tchr1\t20\t30\t-\t11\t{counts[1]}",
        ]
        _write(root / f"{sample_id}.featureCounts.txt", "\n".join(lines) + "\n")
        entries.append(
            {
                "sample_id": sample_id,
                "counts_file": f"{sample_id}.featureCounts.txt",
                "count_column": count_column,
            }
        )
    request["samples"] = entries
    request["featurecounts"] = {"layout": "per_sample_files"}
    return _write_json(root / "request.json", request)


_FC_ANNOTATION_COLUMNS = ("Geneid", "Chr", "Start", "End", "Strand", "Length")


def _salmon_request(
    root: Path,
    *,
    three_prime: bool = False,
    replicate_counts: tuple[int, int] = (0, 0),
    numreads: tuple[str, str] = ("10.25", "20.5"),
) -> tuple[Path, dict[str, object]]:
    semantics = (
        "salmon_quant_dirs_three_prime"
        if three_prime
        else "salmon_quant_dirs_full_length"
    )
    request = _common_request(
        root,
        semantics=semantics,
        producer="Salmon",
        version="1.10.3",
    )
    samples = []
    for sample_index, (sample_id, replicate_count) in enumerate(
        zip(("s1", "s2"), replicate_counts, strict=True)
    ):
        quant_dir = root / "salmon" / sample_id
        quant_lines = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
        quant_lines.extend(
            [
                f"tx1\t100\t80.5\t5.25\t{numreads[0]}",
                f"tx2\t200\t180.25\t10.5\t{numreads[1]}",
            ]
        )
        _write(quant_dir / "quant.sf", "\n".join(quant_lines) + "\n")
        _write_json(
            quant_dir / "cmd_info.json",
            {
                "salmon_version": "1.10.3",
                "index": "/portable/location/synthetic-index",
                "numBootstraps": replicate_count,
                "numGibbsSamples": 0,
                "index_seq_hash": "a" * 64,
            },
        )
        _write_json(
            quant_dir / "aux_info" / "meta_info.json",
            {
                "salmon_version": "1.10.3",
                "num_bootstraps": replicate_count,
                "samp_type": "bootstrap" if replicate_count else "none",
                "num_valid_targets": 2,
                "quant_errors": [],
                "index_seq_hash": "a" * 64,
            },
        )
        if replicate_count:
            bootstrap_dir = quant_dir / "aux_info" / "bootstrap"
            bootstrap_dir.mkdir(parents=True, exist_ok=True)
            with gzip.open(
                bootstrap_dir / "names.tsv.gz", "wt", encoding="utf-8"
            ) as handle:
                handle.write("tx1\ttx2\n")
            values = [
                float(index + sample_index + 1) for index in range(replicate_count * 2)
            ]
            with gzip.open(bootstrap_dir / "bootstraps.gz", "wb") as handle:
                handle.write(struct.pack(f"<{len(values)}d", *values))
        samples.append({"sample_id": sample_id, "quant_dir": f"salmon/{sample_id}"})
    tx2gene = _write(
        root / "tx2gene.tsv",
        "transcript_id\tgene_id\tgene_name\n"
        "tx1\tENSG1.1\tDUP\n"
        "tx2\tENSG2.1\tDUP\n",
    )
    request["samples"] = samples
    request["salmon"] = {
        "tx2gene": "tx2gene.tsv",
        "tx2gene_sha256": _sha256(tx2gene),
    }
    if three_prime:
        request["assay_protocol"] = "three_prime"
    return _write_json(root / "request.json", request), request


def _rewrite_request(path: Path, request: dict[str, object]) -> None:
    _write_json(path, request)


def _rewrite_manifest_digest(root: Path) -> None:
    manifest_path = root / "counts.manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["artifact"]["sha256"] = _sha256(root / "counts.tsv")
    _write_json(manifest_path, manifest)


def test_combined_featurecounts_validates_without_using_symbol_as_key(
    tmp_path: Path,
) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)

    result = validate_request(request_path)

    assert json.loads(json.dumps(result, allow_nan=False)) == result
    data = result["data"]
    assert data["certification_scope"] == "input_semantics_only"
    assert data["input_certification_eligible"] is True
    assert data["featurecounts"]["gene_count"] == 2
    assert data["gene_id_policy"]["internal_key"] == "gene_id"
    assert data["gene_id_policy"]["symbol_role"] == "display_only"
    assert data["route"] == {
        "backend_input": "edgeR::DGEList",
        "count_semantics": "integer",
        "maximum_exact_integer": MAX_EXACT_INTEGER_COUNT,
        "transcript_length_offset": False,
    }


def test_documented_limitation_self_attested_manifest_cannot_prove_origin(
    tmp_path: Path,
) -> None:
    """A forged but self-consistent producer claim passes only the input gate."""

    # These integer values could be rounded Salmon estimated counts. The
    # validator has no authenticated producer record beyond the caller-supplied
    # request and matching manifest, so their true origin is unknowable here.
    request_path, _ = _featurecounts_combined(
        tmp_path,
        counts=(("10", "21"), ("30", "41")),
    )

    result = validate_request(request_path)

    assert result["data"]["input_certification_eligible"] is True
    assert result["data"]["certification_scope"] == "input_semantics_only"
    assert result["data"]["producer"] == {
        "name": "featureCounts",
        "version": "2.0.6",
    }


def test_paths_are_relative_to_request_not_process_cwd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    request_path, _ = _featurecounts_combined(tmp_path / "portable")
    monkeypatch.chdir(tmp_path)

    result = validate_request(Path("portable/request.json"))

    assert result["data"]["featurecounts"]["matrix_path"] == str(
        (tmp_path / "portable" / "counts.tsv").resolve()
    )


def test_artifact_inventory_is_deterministic_and_all_digests_match(
    tmp_path: Path,
) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)

    first = inspect_request(request_path)
    second = inspect_request(request_path)

    assert first == second
    artifacts = first["artifacts"]
    assert artifacts == sorted(
        artifacts,
        key=lambda item: (item["role"], item["path"], item["declared_path"]),
    )
    assert {item["role"] for item in artifacts} == {
        "featurecounts.manifest",
        "featurecounts.manifest_reference",
        "featurecounts.matrix",
        "metadata",
        "reference.source",
        "request",
    }
    for artifact in artifacts:
        assert artifact["sha256"] == _sha256(Path(artifact["path"]))


def test_validation_releases_each_parse_snapshot_before_capturing_the_next(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    inventory_type = input_semantics_module.ArtifactInventory
    original_add = inventory_type.add
    original_verify = inventory_type.verify_unchanged
    peak_snapshot_count = 0

    def tracked_add(self: object, *args: object, **kwargs: object):
        nonlocal peak_snapshot_count
        fingerprint = original_add(self, *args, **kwargs)
        peak_snapshot_count = max(
            peak_snapshot_count, len(self._snapshot_cache)  # type: ignore[attr-defined]
        )
        return fingerprint

    def assert_released_before_verify(self: object) -> None:
        assert not self._snapshot_cache  # type: ignore[attr-defined]
        original_verify(self)

    monkeypatch.setattr(inventory_type, "add", tracked_add)
    monkeypatch.setattr(
        inventory_type, "verify_unchanged", assert_released_before_verify
    )

    validate_request(request_path)

    assert peak_snapshot_count == 1


def test_file_mutation_between_fingerprint_and_parse_is_a_hard_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)
    original_read_tsv = input_semantics_module._read_tsv
    mutated = False

    def mutate_before_parse(path: Path, *, role: str, content: bytes):
        nonlocal mutated
        if role == "featurecounts_matrix" and not mutated:
            mutated = True
            text = path.read_text(encoding="utf-8").replace(
                "ENSG1.1\tDUP\t10\t20", "ENSG1.1\tDUP\t11\t20"
            )
            _write(path, text)
        return original_read_tsv(path, role=role, content=content)

    monkeypatch.setattr(input_semantics_module, "_read_tsv", mutate_before_parse)

    with pytest.raises(InputIntegrityError) as raised:
        validate_request(request_path)

    assert "expected_sha256" in raised.value.details
    assert "observed_sha256" in raised.value.details


def test_aba_path_mutation_cannot_change_the_captured_bytes_being_parsed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)
    matrix_path = tmp_path / "counts.tsv"
    expected_digest = _sha256(matrix_path)
    original_read_tsv = input_semantics_module._read_tsv
    mutated = False

    def mutate_restore_before_parse(path: Path, *, role: str, content: bytes):
        nonlocal mutated
        if role != "featurecounts_matrix" or mutated:
            return original_read_tsv(path, role=role, content=content)
        mutated = True
        original_bytes = path.read_bytes()
        malicious_bytes = original_bytes.replace(b"ENSG1.1", b"EVIL001")
        path.write_bytes(malicious_bytes)
        try:
            return original_read_tsv(path, role=role, content=content)
        finally:
            path.write_bytes(original_bytes)

    monkeypatch.setattr(
        input_semantics_module, "_read_tsv", mutate_restore_before_parse
    )

    result = validate_request(request_path)

    assert mutated is True
    assert result["data"]["featurecounts"]["matrix_sha256"] == expected_digest
    assert result["data"]["featurecounts"]["gene_id_inventory_sha256"] == (
        input_semantics_module._sequence_digest(["ENSG1.1", "ENSG2.1"])
    )


def test_metadata_permutation_is_explicitly_recorded(tmp_path: Path) -> None:
    request_path, _ = _featurecounts_combined(tmp_path, metadata_order=("s2", "s1"))

    result = validate_request(request_path)

    assert result["data"]["sample_order"] == ["s1", "s2"]
    assert result["data"]["metadata"]["sample_order"] == ["s2", "s1"]
    assert [warning["code"] for warning in result["warnings"]] == [
        "METADATA_ORDER_NORMALIZED"
    ]


def test_metadata_sample_set_mismatch_is_a_hard_failure(tmp_path: Path) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)
    _metadata(tmp_path, ("s1", "s3"))

    with pytest.raises(SampleSetMismatchError) as raised:
        validate_request(request_path)

    assert raised.value.details["missing_from_metadata"] == ["s2"]
    assert raised.value.details["unexpected_in_metadata"] == ["s3"]


def test_combined_matrix_without_typed_manifest_is_rejected(tmp_path: Path) -> None:
    request_path, request = _featurecounts_combined(tmp_path)
    request["featurecounts"] = {"layout": "combined_matrix", "matrix": "counts.tsv"}
    _rewrite_request(request_path, request)

    with pytest.raises(InputEvidenceRequiredError):
        validate_request(request_path)


def test_manifest_artifact_digest_mismatch_is_rejected(tmp_path: Path) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)
    with (tmp_path / "counts.tsv").open("a", encoding="utf-8") as handle:
        handle.write("ENSG3\tNEW\t1\t2\n")

    with pytest.raises(InputIntegrityError) as raised:
        validate_request(request_path)

    assert "expected_sha256" in raised.value.details
    assert "observed_sha256" in raised.value.details


@pytest.mark.parametrize("invalid_count", ["1.5", "-1", "NA", "NaN", "Inf", ""])
def test_featurecounts_never_coerces_invalid_count_values(
    tmp_path: Path, invalid_count: str
) -> None:
    request_path, _ = _featurecounts_combined(
        tmp_path, counts=((invalid_count, "20"), ("30", "40"))
    )

    with pytest.raises(CountValuesInvalidError) as raised:
        validate_request(request_path)

    assert raised.value.details["value"] == invalid_count


def test_inspect_does_not_claim_to_validate_count_values(tmp_path: Path) -> None:
    request_path, _ = _featurecounts_combined(
        tmp_path, counts=(("1.5", "20"), ("30", "40"))
    )

    inspected = inspect_request(request_path)

    assert inspected["data"]["validation_level"] == "inspect"
    assert inspected["data"]["input_certification_eligible"] is None
    assert inspected["data"]["input_certification_status"] == "not_evaluated"
    assert inspected["data"]["input_certification_reasons"] == [
        "NUMERIC_DOMAIN_NOT_EVALUATED"
    ]
    with pytest.raises(CountValuesInvalidError):
        validate_request(request_path)


def test_featurecounts_rejects_integer_above_exact_r_numeric_bound(
    tmp_path: Path,
) -> None:
    request_path, _ = _featurecounts_combined(
        tmp_path,
        counts=((str(MAX_EXACT_INTEGER_COUNT + 1), "20"), ("30", "40")),
    )

    with pytest.raises(CountValuesInvalidError) as raised:
        validate_request(request_path)

    assert raised.value.details["maximum_exact_integer"] == MAX_EXACT_INTEGER_COUNT


def test_featurecounts_oversized_digit_literal_maps_to_stable_count_error(
    tmp_path: Path,
) -> None:
    request_path, _ = _featurecounts_combined(
        tmp_path,
        counts=(("9" * 5000, "20"), ("30", "40")),
    )

    with pytest.raises(CountValuesInvalidError) as raised:
        validate_request(request_path)

    assert raised.value.code.value == "COUNT_VALUES_INVALID"
    assert raised.value.details["maximum_exact_integer"] == MAX_EXACT_INTEGER_COUNT


@pytest.mark.parametrize("gene_ids", [("", "ENSG2"), ("ENSG1", "ENSG1")])
def test_empty_or_duplicate_matrix_gene_id_is_rejected(
    tmp_path: Path, gene_ids: tuple[str, str]
) -> None:
    request_path, _ = _featurecounts_combined(tmp_path, gene_ids=gene_ids)

    with pytest.raises(GeneIdentifierError):
        validate_request(request_path)


def test_gene_identifier_with_nul_control_character_is_rejected(
    tmp_path: Path,
) -> None:
    request_path, _ = _featurecounts_combined(
        tmp_path,
        gene_ids=("ENSG1\x00hidden", "ENSG2"),
    )

    with pytest.raises(GeneIdentifierError):
        validate_request(request_path)


def test_version_stripping_collision_is_rejected(tmp_path: Path) -> None:
    request_path, _ = _featurecounts_combined(
        tmp_path,
        gene_ids=("ENSG1.1", "ENSG1.2"),
        strip_version=True,
    )

    with pytest.raises(GeneIdentifierError) as raised:
        validate_request(request_path)

    assert raised.value.details["normalized_gene_id"] == "ENSG1"


def test_per_sample_featurecounts_files_are_accepted_as_upstream_evidence(
    tmp_path: Path,
) -> None:
    request_path = _featurecounts_per_sample(tmp_path)

    result = validate_request(request_path)

    files = result["data"]["featurecounts"]["files"]
    assert [item["sample_id"] for item in files] == ["s1", "s2"]
    assert all(item["producer_version"] == "2.0.6" for item in files)


def test_per_sample_file_requires_featurecounts_version_evidence(
    tmp_path: Path,
) -> None:
    request_path = _featurecounts_per_sample(tmp_path)
    path = tmp_path / "s1.featureCounts.txt"
    text = path.read_text(encoding="utf-8").replace("featureCounts v2.0.6", "unknown")
    _write(path, text)

    with pytest.raises(InputEvidenceRequiredError):
        validate_request(request_path)


def test_per_sample_annotation_fields_must_match_exactly(tmp_path: Path) -> None:
    request_path = _featurecounts_per_sample(tmp_path)
    path = tmp_path / "s2.featureCounts.txt"
    text = path.read_text(encoding="utf-8").replace(
        "ENSG1\tchr1\t1\t10\t+\t10", "ENSG1\tchr1\t1\t10\t+\t999"
    )
    _write(path, text)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_full_length_salmon_accepts_fractional_estimated_counts(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)

    result = validate_request(request_path)

    salmon = result["data"]["salmon"]
    assert salmon["estimated_counts_may_be_fractional"] is True
    assert salmon["index_hash_identity"] == {"index_seq_hash": "a" * 64}
    assert salmon["inferential_replicates"]["state"] == "none"
    route = result["data"]["route"]
    assert route["tximport"] == {"countsFromAbundance": "no", "dropInfReps": False}
    assert route["edgeR_constructor"] == "edgeR::DGEListFromTximport"
    assert route["transcript_length_offset"] is True
    assert route["edgeR_options"] == {"divide": False}
    assert route["inferential_overdispersion"] == {
        "enabled": False,
        "source": "tximport.infReps",
        "relative_abundance_adjustment": False,
    }


def test_salmon_custom_cmd_info_aux_dir_is_used_consistently(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(2, 2))
    for sample_id in ("s1", "s2"):
        quant_dir = tmp_path / "salmon" / sample_id
        (quant_dir / "aux_info").rename(quant_dir / "custom_aux")
        cmd_path = quant_dir / "cmd_info.json"
        cmd = json.loads(cmd_path.read_text(encoding="utf-8"))
        cmd["auxDir"] = "custom_aux"
        _write_json(cmd_path, cmd)

    result = validate_request(request_path)

    assert {
        sample["aux_dir_declared"] for sample in result["data"]["salmon"]["samples"]
    } == {"custom_aux"}
    assert all(
        artifact["declared_path"].endswith("/custom_aux/meta_info.json")
        for artifact in result["artifacts"]
        if artifact["role"].startswith("salmon.meta_info.")
    )
    assert result["data"]["route"]["edgeR_options"] == {"divide": True}


def test_salmon_custom_aux_dir_must_exist(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)
    cmd_path = tmp_path / "salmon" / "s1" / "cmd_info.json"
    cmd = json.loads(cmd_path.read_text(encoding="utf-8"))
    cmd["auxDir"] = "missing_aux"
    _write_json(cmd_path, cmd)

    with pytest.raises(InputEvidenceRequiredError):
        validate_request(request_path)


def test_salmon_aux_dir_must_not_be_absolute(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)
    cmd_path = tmp_path / "salmon" / "s1" / "cmd_info.json"
    cmd = json.loads(cmd_path.read_text(encoding="utf-8"))
    cmd["auxDir"] = str((tmp_path / "salmon" / "s1" / "aux_info").resolve())
    _write_json(cmd_path, cmd)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_salmon_aux_dir_must_not_escape_quant_dir(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)
    cmd_path = tmp_path / "salmon" / "s1" / "cmd_info.json"
    cmd = json.loads(cmd_path.read_text(encoding="utf-8"))
    cmd["auxDir"] = "../s2/aux_info"
    _write_json(cmd_path, cmd)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_salmon_aux_dir_symlink_must_not_escape_quant_dir(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)
    quant_dir = tmp_path / "salmon" / "s1"
    outside = tmp_path / "outside_aux"
    outside.mkdir()
    (quant_dir / "linked_aux").symlink_to(outside, target_is_directory=True)
    cmd_path = quant_dir / "cmd_info.json"
    cmd = json.loads(cmd_path.read_text(encoding="utf-8"))
    cmd["auxDir"] = "linked_aux"
    _write_json(cmd_path, cmd)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


@pytest.mark.parametrize(
    "invalid_count",
    ["-0.1", "NA", "NaN", "Infinity", "", "1e400", "1_0", " 1"],
)
def test_salmon_rejects_nonfinite_missing_or_negative_estimated_counts(
    tmp_path: Path, invalid_count: str
) -> None:
    request_path, _ = _salmon_request(tmp_path, numreads=(invalid_count, "20.5"))

    with pytest.raises(CountValuesInvalidError):
        validate_request(request_path)


def test_full_length_all_inferential_replicates_enable_dgelist_divide(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))

    result = validate_request(request_path)

    summary = result["data"]["salmon"]["inferential_replicates"]
    assert summary["state"] == "all"
    assert summary["method"] == "bootstrap"
    assert summary["replicate_count"] == 10
    assert {record["archive_numeric_encoding"] for record in summary["per_sample"]} == {
        "float64_le"
    }
    assert result["data"]["route"]["edgeR_options"] == {"divide": True}
    overdispersion = result["data"]["route"]["inferential_overdispersion"]
    assert overdispersion == {
        "enabled": True,
        "source": "tximport.infReps",
        "relative_abundance_adjustment": True,
    }
    replicate_roles = {
        artifact["role"]
        for artifact in result["artifacts"]
        if "inferential_replicates" in artifact["role"]
    }
    assert len(replicate_roles) == 4


def test_full_length_single_inferential_replicate_is_inspectable_but_ineligible(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(1, 1))

    result = inspect_request(request_path)

    summary = result["data"]["salmon"]["inferential_replicates"]
    assert summary["state"] == "all"
    assert summary["consistent_method_and_count"] is True
    assert summary["replicate_count"] == 1
    assert result["data"]["input_certification_eligible"] is False
    assert result["data"]["input_certification_status"] == "ineligible"
    assert result["data"]["input_certification_reasons"] == [
        "INFERENTIAL_REPLICATES_INSUFFICIENT"
    ]
    assert result["data"]["route"]["edgeR_options"] == {"divide": False}
    assert result["data"]["route"]["inferential_overdispersion"] == {
        "enabled": False,
        "source": "tximport.infReps",
        "relative_abundance_adjustment": False,
    }
    warning = result["warnings"][-1]
    assert warning["code"] == "INFERENTIAL_REPLICATES_INSUFFICIENT"
    assert warning["severity"] == "high"
    assert warning["details"]["minimum_replicates_per_sample"] == 2


def test_full_length_single_inferential_replicate_fails_validation(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(1, 1))

    with pytest.raises(InputValidationError) as raised:
        validate_request(request_path)

    assert raised.value.code.value == "INPUT_VALIDATION_FAILED"
    assert raised.value.details["reason"] == "inferential_replicate_count_below_minimum"
    assert raised.value.details["observed_replicates_per_sample"] == 1
    assert raised.value.details["minimum_replicates_per_sample"] == 2


def test_three_prime_single_inferential_replicate_remains_valid(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(
        tmp_path,
        three_prime=True,
        replicate_counts=(1, 1),
    )

    result = validate_request(request_path)

    assert result["data"]["input_certification_eligible"] is True
    assert result["data"]["salmon"]["inferential_replicates"]["replicate_count"] == 1
    assert "edgeR_options" not in result["data"]["route"]
    assert "inferential_overdispersion" not in result["data"]["route"]


def test_gibbs_replicates_are_detected_without_conflating_bootstrap_zero(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    for sample_id in ("s1", "s2"):
        cmd_path = tmp_path / "salmon" / sample_id / "cmd_info.json"
        cmd = json.loads(cmd_path.read_text(encoding="utf-8"))
        cmd["numBootstraps"] = 0
        cmd["numGibbsSamples"] = 10
        _write_json(cmd_path, cmd)
        meta_path = tmp_path / "salmon" / sample_id / "aux_info" / "meta_info.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        meta["num_bootstraps"] = 10
        meta["samp_type"] = "gibbs"
        _write_json(meta_path, meta)

    result = validate_request(request_path)

    summary = result["data"]["salmon"]["inferential_replicates"]
    assert summary["state"] == "all"
    assert summary["method"] == "gibbs"
    assert summary["replicate_count"] == 10


@pytest.mark.parametrize(
    ("replicate_count", "samp_type"),
    [(0, "bootstrap"), (0, "unknown"), (10, "none"), (10, "unknown")],
)
def test_salmon_samp_type_must_be_known_and_consistent_with_replicate_count(
    tmp_path: Path, replicate_count: int, samp_type: str
) -> None:
    request_path, _ = _salmon_request(
        tmp_path, replicate_counts=(replicate_count, replicate_count)
    )
    meta_path = tmp_path / "salmon" / "s1" / "aux_info" / "meta_info.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["samp_type"] = samp_type
    _write_json(meta_path, meta)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_positive_salmon_replicates_require_relevant_command_count_evidence(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    cmd_path = tmp_path / "salmon" / "s1" / "cmd_info.json"
    cmd = json.loads(cmd_path.read_text(encoding="utf-8"))
    cmd.pop("numBootstraps")
    _write_json(cmd_path, cmd)

    with pytest.raises(InputEvidenceRequiredError) as raised:
        validate_request(request_path)

    assert raised.value.details["method"] == "bootstrap"


def test_salmon_samp_type_is_required_even_when_replicate_count_is_zero(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path)
    meta_path = tmp_path / "salmon" / "s1" / "aux_info" / "meta_info.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta.pop("samp_type")
    _write_json(meta_path, meta)

    with pytest.raises(InputEvidenceRequiredError):
        validate_request(request_path)


def test_full_length_zero_effective_length_has_offset_specific_error(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path)
    quant_path = tmp_path / "salmon" / "s1" / "quant.sf"
    text = quant_path.read_text(encoding="utf-8").replace(
        "tx1\t100\t80.5", "tx1\t100\t0"
    )
    _write(quant_path, text)

    with pytest.raises(SalmonOffsetRequiredError):
        validate_request(request_path)


def test_full_length_malformed_quant_header_has_offset_specific_error(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path)
    quant_path = tmp_path / "salmon" / "s1" / "quant.sf"
    text = quant_path.read_text(encoding="utf-8").replace(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads",
        "Name\tLength\teffective_length\tTPM\tNumReads",
    )
    _write(quant_path, text)

    with pytest.raises(SalmonOffsetRequiredError) as raised:
        validate_request(request_path)

    assert raised.value.code.value == "SALMON_OFFSET_REQUIRED"
    assert raised.value.details["observed_header"][2] == "effective_length"


def test_three_prime_zero_effective_length_does_not_create_an_offset_requirement(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, three_prime=True)
    for sample_id in ("s1", "s2"):
        quant_path = tmp_path / "salmon" / sample_id / "quant.sf"
        text = quant_path.read_text(encoding="utf-8").replace(
            "tx1\t100\t80.5", "tx1\t100\t0"
        )
        _write(quant_path, text)

    result = validate_request(request_path)

    assert result["data"]["route"]["transcript_length_offset"] is False


def test_mixed_inferential_replicates_are_visible_to_inspect_and_fail_validate(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 0))

    inspected = inspect_request(request_path)

    assert inspected["data"]["salmon"]["inferential_replicates"]["state"] == "mixed"
    assert inspected["data"]["input_certification_eligible"] is False
    assert inspected["data"]["input_certification_status"] == "ineligible"
    assert inspected["data"]["input_certification_reasons"] == [
        "INFERENTIAL_REPLICATES_INCONSISTENT"
    ]
    assert inspected["warnings"][-1]["code"] == "INFERENTIAL_REPLICATES_INCONSISTENT"
    with pytest.raises(InputValidationError):
        validate_request(request_path)


def test_declared_replicates_require_bootstrap_archive(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    archive = tmp_path / "salmon" / "s1" / "aux_info" / "bootstrap" / "bootstraps.gz"
    archive.unlink()

    with pytest.raises(InputEvidenceRequiredError):
        validate_request(request_path)


def test_corrupt_replicate_gzip_is_rejected_before_divide_is_enabled(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    archive = tmp_path / "salmon" / "s1" / "aux_info" / "bootstrap" / "bootstraps.gz"
    _write(archive, "not-a-gzip-stream\n")

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_zlib_decompression_errors_map_to_input_integrity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    original_read = gzip.GzipFile.read
    raised_once = False

    def raise_zlib_once(
        handle: gzip.GzipFile, *args: object, **kwargs: object
    ) -> bytes:
        nonlocal raised_once
        if not raised_once:
            raised_once = True
            raise zlib.error("synthetic corrupt deflate payload")
        return original_read(handle, *args, **kwargs)

    monkeypatch.setattr(gzip.GzipFile, "read", raise_zlib_once)

    with pytest.raises(InputIntegrityError) as raised:
        validate_request(request_path)

    assert raised.value.cause is not None
    assert isinstance(raised.value.cause, zlib.error)


def test_replicate_names_must_match_quant_transcript_order(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    names = tmp_path / "salmon" / "s1" / "aux_info" / "bootstrap" / "names.tsv.gz"
    with gzip.open(names, "wt", encoding="utf-8") as handle:
        handle.write("tx2\ttx1\n")

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_replicate_target_count_must_match_quant_rows(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    meta_path = tmp_path / "salmon" / "s1" / "aux_info" / "meta_info.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["num_valid_targets"] = 3
    _write_json(meta_path, meta)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_replicate_binary_dimensions_must_match_metadata(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    archive = tmp_path / "salmon" / "s1" / "aux_info" / "bootstrap" / "bootstraps.gz"
    with gzip.open(archive, "wb") as handle:
        handle.write(struct.pack("<d", 1.0))

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_legacy_int32_replicate_payload_supported_by_locked_tximport_is_valid(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, replicate_counts=(10, 10))
    for sample_id in ("s1", "s2"):
        archive = (
            tmp_path / "salmon" / sample_id / "aux_info" / "bootstrap" / "bootstraps.gz"
        )
        with gzip.open(archive, "wb") as handle:
            handle.write(struct.pack("<20i", *range(1, 21)))

    result = validate_request(request_path)

    summary = result["data"]["salmon"]["inferential_replicates"]
    assert {record["archive_numeric_encoding"] for record in summary["per_sample"]} == {
        "int32_le"
    }


def test_explicit_zero_and_nonzero_replicate_metadata_disagreement_fails(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path)
    meta_path = tmp_path / "salmon" / "s1" / "aux_info" / "meta_info.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["num_bootstraps"] = 10
    _write_json(meta_path, meta)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_salmon_oversized_string_integer_metadata_maps_to_integrity_error(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path)
    meta_path = tmp_path / "salmon" / "s1" / "aux_info" / "meta_info.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["num_bootstraps"] = "9" * 5000
    _write_json(meta_path, meta)

    with pytest.raises(InputIntegrityError) as raised:
        validate_request(request_path)

    assert raised.value.details["maximum_exact_integer"] == MAX_EXACT_INTEGER_COUNT


def test_salmon_requires_index_hash_identity_for_every_sample(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)
    for relative in ("cmd_info.json", "aux_info/meta_info.json"):
        path = tmp_path / "salmon" / "s1" / relative
        document = json.loads(path.read_text(encoding="utf-8"))
        document.pop("index_seq_hash")
        _write_json(path, document)

    with pytest.raises(InputEvidenceRequiredError):
        validate_request(request_path)


def test_salmon_rejects_malformed_index_hash_evidence(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)
    path = tmp_path / "salmon" / "s1" / "aux_info" / "meta_info.json"
    document = json.loads(path.read_text(encoding="utf-8"))
    document["index_seq_hash"] = "abc123"
    _write_json(path, document)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_salmon_version_is_required_in_both_metadata_documents(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)
    path = tmp_path / "salmon" / "s1" / "cmd_info.json"
    document = json.loads(path.read_text(encoding="utf-8"))
    document.pop("salmon_version")
    _write_json(path, document)

    with pytest.raises(InputEvidenceRequiredError):
        validate_request(request_path)


def test_salmon_index_hashes_must_match_across_samples(tmp_path: Path) -> None:
    request_path, _ = _salmon_request(tmp_path)
    for relative in ("cmd_info.json", "aux_info/meta_info.json"):
        path = tmp_path / "salmon" / "s2" / relative
        document = json.loads(path.read_text(encoding="utf-8"))
        document["index_seq_hash"] = "b" * 64
        _write_json(path, document)

    with pytest.raises(InputIntegrityError):
        validate_request(request_path)


def test_salmon_index_path_difference_is_provenance_warning_when_hashes_match(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path)
    path = tmp_path / "salmon" / "s2" / "cmd_info.json"
    document = json.loads(path.read_text(encoding="utf-8"))
    document["index"] = "/different/index"
    _write_json(path, document)

    result = validate_request(request_path)

    assert result["data"]["input_certification_eligible"] is True
    assert [warning["code"] for warning in result["warnings"]] == [
        "SALMON_INDEX_PATHS_DIFFER"
    ]


def test_tx2gene_must_exactly_cover_quantified_transcripts(tmp_path: Path) -> None:
    request_path, request = _salmon_request(tmp_path)
    tx2gene = _write(
        tmp_path / "tx2gene.tsv",
        "transcript_id\tgene_id\n" "tx1\tENSG1\n" "tx3\tENSG3\n",
    )
    request["salmon"]["tx2gene_sha256"] = _sha256(tx2gene)
    _rewrite_request(request_path, request)

    with pytest.raises(GeneIdentifierError) as raised:
        validate_request(request_path)

    assert raised.value.details["missing_count"] == 1
    assert raised.value.details["extra_count"] == 1


def test_tx2gene_rejects_duplicate_transcript_mapping(tmp_path: Path) -> None:
    request_path, request = _salmon_request(tmp_path)
    tx2gene = _write(
        tmp_path / "tx2gene.tsv",
        "transcript_id\tgene_id\n" "tx1\tENSG1\n" "tx1\tENSG2\n",
    )
    request["salmon"]["tx2gene_sha256"] = _sha256(tx2gene)
    _rewrite_request(request_path, request)

    with pytest.raises(GeneIdentifierError):
        validate_request(request_path)


def test_three_prime_protocol_must_be_literal(tmp_path: Path) -> None:
    request_path, request = _salmon_request(tmp_path, three_prime=True)
    request["assay_protocol"] = "3prime"
    _rewrite_request(request_path, request)

    with pytest.raises(AssayProtocolRequiredError):
        validate_request(request_path)


def test_full_length_request_rejects_three_prime_protocol_route_confusion(
    tmp_path: Path,
) -> None:
    request_path, request = _salmon_request(tmp_path)
    request["assay_protocol"] = "three_prime"
    _rewrite_request(request_path, request)

    with pytest.raises(AssayProtocolRequiredError):
        validate_request(request_path)


def test_three_prime_request_rejects_full_length_protocol_route_confusion(
    tmp_path: Path,
) -> None:
    request_path, request = _salmon_request(tmp_path, three_prime=True)
    request["assay_protocol"] = "full_length"
    _rewrite_request(request_path, request)

    with pytest.raises(AssayProtocolRequiredError):
        validate_request(request_path)


def test_salmon_request_rejects_featurecounts_route_declaration(
    tmp_path: Path,
) -> None:
    request_path, request = _salmon_request(tmp_path)
    request["featurecounts"] = {"layout": "per_sample_files"}
    _rewrite_request(request_path, request)

    with pytest.raises(InvalidRequestError):
        validate_request(request_path)


def test_featurecounts_request_rejects_salmon_route_declaration(
    tmp_path: Path,
) -> None:
    request_path, request = _featurecounts_combined(tmp_path)
    request["salmon"] = {"tx2gene": "tx2gene.tsv"}
    _rewrite_request(request_path, request)

    with pytest.raises(InvalidRequestError):
        validate_request(request_path)


def test_three_prime_route_uses_raw_counts_without_length_offset(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path, three_prime=True)

    result = validate_request(request_path)

    route = result["data"]["route"]
    assert route["count_source"] == "txi$counts"
    assert route["tximport"]["countsFromAbundance"] == "no"
    assert route["transcript_length_offset"] is False
    assert route["gene_length_correction"] is False
    assert route["route_interpretation"] == "certified_three_prime_default"
    assert route["certified_path_execution_permitted"] is True
    assert result["data"]["input_certification_eligible"] is True


def test_three_prime_override_requires_reason(tmp_path: Path) -> None:
    request_path, request = _salmon_request(tmp_path, three_prime=True)
    request["analysis_options"] = {
        "three_prime_length_correction_override": {"enabled": True}
    }
    _rewrite_request(request_path, request)

    with pytest.raises(InputEvidenceRequiredError):
        validate_request(request_path)


def test_three_prime_override_is_high_risk_and_noncertifying(tmp_path: Path) -> None:
    request_path, request = _salmon_request(tmp_path, three_prime=True)
    request["analysis_options"] = {
        "three_prime_length_correction_override": {
            "enabled": True,
            "reason": "External protocol evidence requires a sensitivity-only comparison.",
        }
    }
    _rewrite_request(request_path, request)

    result = validate_request(request_path)

    assert result["data"]["input_certification_eligible"] is False
    assert result["warnings"][-1]["severity"] == "high"
    assert result["warnings"][-1]["code"] == (
        "HIGH_RISK_THREE_PRIME_LENGTH_CORRECTION_OVERRIDE"
    )
    route = result["data"]["route"]
    assert route["certified_path_execution_permitted"] is False
    assert route["high_risk_override"]["execution_policy"] == (
        "not_executable_in_certified_path"
    )


def test_unknown_request_fields_are_rejected(tmp_path: Path) -> None:
    request_path, request = _featurecounts_combined(tmp_path)
    request["guess_input_type"] = True
    _rewrite_request(request_path, request)

    with pytest.raises(InvalidRequestError) as raised:
        validate_request(request_path)

    assert raised.value.details["unknown_keys"] == ["guess_input_type"]


def test_duplicate_json_keys_are_rejected_instead_of_last_value_winning(
    tmp_path: Path,
) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)
    text = request_path.read_text(encoding="utf-8").replace(
        '"schema_version": "1.0"',
        '"schema_version": "1.0",\n  "schema_version": "1.0"',
        1,
    )
    _write(request_path, text)

    with pytest.raises(InvalidRequestError) as raised:
        validate_request(request_path)

    assert raised.value.details["duplicate_key"] == "schema_version"


def test_nonstandard_json_numeric_constants_are_rejected(tmp_path: Path) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)
    text = request_path.read_text(encoding="utf-8").replace(
        '"schema_version": "1.0"', '"schema_version": NaN', 1
    )
    _write(request_path, text)

    with pytest.raises(InvalidRequestError) as raised:
        validate_request(request_path)

    assert raised.value.details["value"] == "NaN"


def test_oversized_json_integer_maps_to_stable_invalid_request_error(
    tmp_path: Path,
) -> None:
    request_path, _ = _featurecounts_combined(tmp_path)
    text = request_path.read_text(encoding="utf-8").replace(
        '"strip_version": false', f'"strip_version": {"9" * 5000}', 1
    )
    _write(request_path, text)

    with pytest.raises(InvalidRequestError) as raised:
        validate_request(request_path)

    assert raised.value.code.value == "INVALID_REQUEST"
    assert raised.value.details["maximum_absolute_integer"] == (MAX_EXACT_INTEGER_COUNT)


def test_oversized_raw_integer_in_salmon_metadata_maps_to_stable_json_error(
    tmp_path: Path,
) -> None:
    request_path, _ = _salmon_request(tmp_path)
    meta_path = tmp_path / "salmon" / "s1" / "aux_info" / "meta_info.json"
    text = meta_path.read_text(encoding="utf-8").replace(
        '"num_bootstraps": 0', f'"num_bootstraps": {"9" * 5000}', 1
    )
    _write(meta_path, text)

    with pytest.raises(InvalidRequestError) as raised:
        validate_request(request_path)

    assert raised.value.details["document_role"] == "salmon_meta_info"


@pytest.mark.parametrize("control", ["\x00", "\x01", "\x7f"])
def test_schema_identity_strings_reject_unicode_control_characters(
    tmp_path: Path, control: str
) -> None:
    request_path, request = _featurecounts_combined(tmp_path)
    request["producer"]["version"] = f"2.0.6{control}hidden"
    _rewrite_request(request_path, request)

    with pytest.raises(InvalidRequestError):
        validate_request(request_path)


def test_sequence_inventory_digest_uses_unambiguous_length_prefixes() -> None:
    assert input_semantics_module._sequence_digest(["a\x00b"]) != (
        input_semantics_module._sequence_digest(["a", "b"])
    )


def test_schema_strings_are_not_silently_trimmed(tmp_path: Path) -> None:
    request_path, request = _featurecounts_combined(tmp_path)
    request["producer"]["version"] = " 2.0.6"
    _rewrite_request(request_path, request)

    with pytest.raises(InvalidRequestError):
        validate_request(request_path)
