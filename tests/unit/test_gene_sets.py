"""Request-schema and source-evidence tests for optional frozen gene sets."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from rnaseq_downstream.errors import (
    InputIntegrityError,
    InputReadError,
    InvalidRequestError,
)
from rnaseq_downstream.gene_sets import parse_gene_sets_request


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sources(root: Path) -> tuple[Path, Path]:
    gmt = root / "sets.gmt"
    annotation = root / "annotation.tsv"
    gmt.write_text("SET_A\tdescription\tA\tB\n", encoding="utf-8")
    annotation.write_text("gene_id\tsymbol\nENSG1\tA\nENSG2\tB\n", encoding="utf-8")
    return gmt, annotation


def _request(root: Path) -> tuple[Path, dict[str, object]]:
    gmt, annotation = _sources(root)
    request_path = root / "analysis.json"
    request_path.write_text("{}\n", encoding="utf-8")
    return request_path, {
        "gmt": {
            "path": gmt.name,
            "sha256": _sha256(gmt),
            "collection": "synthetic_hallmarks",
            "version": "2026.1",
            "identifier_type": "symbol",
        },
        "annotation": {
            "path": annotation.name,
            "sha256": _sha256(annotation),
            "name": "synthetic_ensembl",
            "version": "1",
            "gene_id_column": "gene_id",
            "symbol_column": "symbol",
        },
        "minimum_tested_genes": 2,
    }


def test_v11_pathway_example_binds_the_bundled_frozen_sources() -> None:
    example_root = Path(__file__).parents[2] / "examples" / "analysis-requests"
    request_path = example_root / "v1.1-pathways.request.json"
    document = json.loads(request_path.read_text(encoding="utf-8"))

    assert document["schema_version"] == "1.1"
    assert set(document) == {
        "schema_version",
        "validated_input_bundle",
        "design",
        "contrasts",
        "display",
        "gene_sets",
    }
    gene_sets = document["gene_sets"]
    for role in ("gmt", "annotation"):
        source = example_root / gene_sets[role]["path"]
        assert source.is_file()
        assert gene_sets[role]["sha256"] == _sha256(source)


def test_gene_sets_sources_are_resolved_and_captured_from_declared_bytes(
    tmp_path: Path,
) -> None:
    request_path, document = _request(tmp_path)

    normalized = parse_gene_sets_request(document, request_path=request_path)

    for role, filename in (("gmt", "sets.gmt"), ("annotation", "annotation.tsv")):
        source = tmp_path / filename
        assert normalized[role]["declared_path"] == filename
        assert normalized[role]["path"] == str(source.resolve())
        assert normalized[role]["sha256"] == _sha256(source)
        assert normalized[role]["size_bytes"] == source.stat().st_size
    assert normalized["gmt"]["identifier_type"] == "symbol"
    assert normalized["annotation"]["gene_id_column"] == "gene_id"
    assert normalized["annotation"]["symbol_column"] == "symbol"
    assert normalized["minimum_tested_genes"] == 2


def test_gene_sets_contract_does_not_interpret_valid_utf8_source_contents(
    tmp_path: Path,
) -> None:
    request_path, document = _request(tmp_path)
    gmt = tmp_path / "sets.gmt"
    gmt.write_text("valid UTF-8 but not parsed at this boundary\n", encoding="utf-8")
    document["gmt"]["sha256"] = _sha256(gmt)  # type: ignore[index]

    normalized = parse_gene_sets_request(document, request_path=request_path)

    assert normalized["gmt"]["size_bytes"] == gmt.stat().st_size


@pytest.mark.parametrize(
    ("location", "mutation"),
    [
        ("root", {"unexpected": True}),
        ("gmt", {"unexpected": True}),
        ("annotation", {"unexpected": True}),
    ],
)
def test_gene_sets_exact_objects_reject_unknown_keys(
    tmp_path: Path,
    location: str,
    mutation: dict[str, object],
) -> None:
    request_path, document = _request(tmp_path)
    target = document if location == "root" else document[location]
    assert isinstance(target, dict)
    target.update(mutation)

    with pytest.raises(InvalidRequestError) as caught:
        parse_gene_sets_request(document, request_path=request_path)

    assert caught.value.details["unknown_keys"] == ["unexpected"]


@pytest.mark.parametrize(
    ("location", "key"),
    [
        ("root", "minimum_tested_genes"),
        ("gmt", "version"),
        ("annotation", "name"),
    ],
)
def test_gene_sets_exact_objects_reject_missing_keys(
    tmp_path: Path,
    location: str,
    key: str,
) -> None:
    request_path, document = _request(tmp_path)
    target = document if location == "root" else document[location]
    assert isinstance(target, dict)
    del target[key]

    with pytest.raises(InvalidRequestError) as caught:
        parse_gene_sets_request(document, request_path=request_path)

    assert caught.value.details["missing_keys"] == [key]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("identifier_type", "ensembl_gene_id"),
        ("gene_id_column", "Geneid"),
        ("symbol_column", "gene_name"),
    ],
)
def test_gene_sets_fixed_identity_fields_are_not_configurable(
    tmp_path: Path,
    field: str,
    value: str,
) -> None:
    request_path, document = _request(tmp_path)
    location = "gmt" if field == "identifier_type" else "annotation"
    nested = document[location]
    assert isinstance(nested, dict)
    nested[field] = value

    with pytest.raises(InvalidRequestError):
        parse_gene_sets_request(document, request_path=request_path)


@pytest.mark.parametrize("minimum", [True, False, 0, -1, 1.0, "10"])
def test_gene_sets_minimum_tested_genes_requires_positive_json_integer(
    tmp_path: Path,
    minimum: object,
) -> None:
    request_path, document = _request(tmp_path)
    document["minimum_tested_genes"] = minimum

    with pytest.raises(InvalidRequestError):
        parse_gene_sets_request(document, request_path=request_path)


@pytest.mark.parametrize("digest", ["A" * 64, "a" * 63, "g" * 64, 123])
def test_gene_sets_digest_must_be_canonical_lowercase_sha256(
    tmp_path: Path,
    digest: object,
) -> None:
    request_path, document = _request(tmp_path)
    gmt = document["gmt"]
    assert isinstance(gmt, dict)
    gmt["sha256"] = digest

    with pytest.raises(InvalidRequestError):
        parse_gene_sets_request(document, request_path=request_path)


def test_gene_sets_digest_mismatch_is_an_integrity_failure(tmp_path: Path) -> None:
    request_path, document = _request(tmp_path)
    gmt = document["gmt"]
    assert isinstance(gmt, dict)
    gmt["sha256"] = "0" * 64

    with pytest.raises(InputIntegrityError) as caught:
        parse_gene_sets_request(document, request_path=request_path)

    assert caught.value.details["role"] == "gmt"
    assert caught.value.details["expected_sha256"] == "0" * 64
    assert caught.value.details["observed_sha256"] == _sha256(tmp_path / "sets.gmt")


def test_gene_sets_source_change_during_capture_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request_path, document = _request(tmp_path)
    gmt = (tmp_path / "sets.gmt").resolve()
    original_open = Path.open

    class ChangingReader:
        def __init__(self, handle: object) -> None:
            self.handle = handle

        def __enter__(self) -> "ChangingReader":
            self.handle.__enter__()  # type: ignore[attr-defined]
            return self

        def read(self) -> bytes:
            content = self.handle.read()  # type: ignore[attr-defined]
            with original_open(gmt, "ab") as mutation:
                mutation.write(b"changed\n")
            return content

        def __exit__(self, *args: object) -> object:
            return self.handle.__exit__(*args)  # type: ignore[attr-defined]

    def changing_open(path: Path, mode: str = "r", *args: object, **kwargs: object):
        handle = original_open(path, mode, *args, **kwargs)
        if path == gmt and mode == "rb":
            return ChangingReader(handle)
        return handle

    monkeypatch.setattr(Path, "open", changing_open)

    with pytest.raises(InputIntegrityError) as caught:
        parse_gene_sets_request(document, request_path=request_path)

    assert caught.value.details["role"] == "gmt"
    assert caught.value.details["source_identity_stable"] is False


def test_gene_sets_source_must_be_utf8(tmp_path: Path) -> None:
    request_path, document = _request(tmp_path)
    gmt = tmp_path / "sets.gmt"
    gmt.write_bytes(b"SET_A\t\xff\n")
    nested = document["gmt"]
    assert isinstance(nested, dict)
    nested["sha256"] = _sha256(gmt)

    with pytest.raises(InputReadError) as caught:
        parse_gene_sets_request(document, request_path=request_path)

    assert caught.value.details["operation"] == "decode_gene_set_source"


@pytest.mark.parametrize(
    "path", ["https://example.invalid/sets.gmt", "file:///sets.gmt"]
)
def test_gene_sets_network_and_uri_paths_are_rejected(
    tmp_path: Path,
    path: str,
) -> None:
    request_path, document = _request(tmp_path)
    gmt = document["gmt"]
    assert isinstance(gmt, dict)
    gmt["path"] = path

    with pytest.raises(InvalidRequestError) as caught:
        parse_gene_sets_request(document, request_path=request_path)

    assert "local files" in caught.value.message


def test_gene_sets_control_character_in_path_is_rejected(tmp_path: Path) -> None:
    request_path, document = _request(tmp_path)
    gmt = document["gmt"]
    assert isinstance(gmt, dict)
    gmt["path"] = "sets\x00.gmt"

    with pytest.raises(InvalidRequestError):
        parse_gene_sets_request(document, request_path=request_path)


@pytest.mark.parametrize("kind", ["symlink", "directory", "missing"])
def test_gene_sets_source_must_be_an_existing_non_symlink_regular_file(
    tmp_path: Path,
    kind: str,
) -> None:
    request_path, document = _request(tmp_path)
    nested = document["gmt"]
    assert isinstance(nested, dict)
    if kind == "symlink":
        source = tmp_path / "sets-link.gmt"
        source.symlink_to("sets.gmt")
    elif kind == "directory":
        source = tmp_path / "sets-dir.gmt"
        source.mkdir()
    else:
        source = tmp_path / "missing.gmt"
    nested["path"] = source.name

    with pytest.raises(InputReadError) as caught:
        parse_gene_sets_request(document, request_path=request_path)

    assert caught.value.details["operation"] == "capture_gene_set_source"


@pytest.mark.parametrize("value", [None, [], "sets.gmt"])
def test_gene_sets_root_must_be_an_object(tmp_path: Path, value: object) -> None:
    request_path, _ = _request(tmp_path)

    with pytest.raises(InvalidRequestError):
        parse_gene_sets_request(value, request_path=request_path)
