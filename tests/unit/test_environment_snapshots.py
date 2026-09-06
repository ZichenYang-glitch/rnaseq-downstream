"""Fail-closed tests for append-only environment evidence snapshots."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil

import pytest

from scripts.benchmark.evidence_resolver import (
    BenchmarkError,
    load_environment_snapshots,
    manifest_evidence,
    resolve_archived_implementation_path,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT_ID = "p0-c1-c2-c6b6cd9"
SNAPSHOT_ROOT = PROJECT_ROOT / "environment" / "snapshots"
SNAPSHOT_MANIFEST_SHA256 = (
    "6630cee71660c44a0264a47c81fb2dcce7dbc5d3152536d148b7594a1b64488f"
)
SNAPSHOT_MANIFEST_SIZE = 854
SOURCE_REVISION = "c6b6cd970c1a0c6807145dd01033505efae60215"
EXPECTED_FILES = {
    "renv.lock": (
        "8d8ca4c7b9ab411150f207846d098cad799dcefe7219696a06e3ebf17cbe19ff",
        268638,
    ),
    "environment/r-sources.lock": (
        "0632d66bd4aee5cc00e7ee33f6061df0db3803843bd5bc98a91fb83c82e30f4d",
        1312,
    ),
    "environment/verify.R": (
        "0f13c98256143ae1febcb3750cd96b1c9d897d8a5cdfee1481018e12b1fc3a7e",
        8494,
    ),
}


@pytest.mark.unit
def test_p0_c1_c2_snapshot_is_immutable_and_manifest_bound() -> None:
    records = load_environment_snapshots()
    selected = [record for record in records if record.snapshot_id == SNAPSHOT_ID]

    assert len(selected) == len(EXPECTED_FILES)
    assert {record.source_revision for record in selected} == {SOURCE_REVISION}
    assert {
        record.source_path: (record.sha256, record.size_bytes) for record in selected
    } == EXPECTED_FILES
    assert all(record.snapshot_path == record.source_path for record in selected)
    assert manifest_evidence(SNAPSHOT_ID) == {
        "name": f"environment/snapshots/{SNAPSHOT_ID}/manifest.json",
        "sha256": SNAPSHOT_MANIFEST_SHA256,
        "size_bytes": SNAPSHOT_MANIFEST_SIZE,
    }


@pytest.mark.unit
@pytest.mark.parametrize("source_path", sorted(EXPECTED_FILES))
def test_archived_environment_evidence_resolves_to_snapshot(
    source_path: str,
) -> None:
    digest, size = EXPECTED_FILES[source_path]
    current = PROJECT_ROOT / source_path

    assert current.is_file()
    assert hashlib.sha256(current.read_bytes()).hexdigest() != digest
    resolved = resolve_archived_implementation_path(
        current,
        expected_sha256=digest,
        expected_size=size,
    )

    assert resolved == SNAPSHOT_ROOT / SNAPSHOT_ID / source_path
    assert hashlib.sha256(resolved.read_bytes()).hexdigest() == digest
    assert resolved.stat().st_size == size


@pytest.mark.unit
def test_resolver_prefers_an_exact_current_file() -> None:
    current = PROJECT_ROOT / "conda-lock.yml"
    payload = current.read_bytes()

    resolved = resolve_archived_implementation_path(
        current,
        expected_sha256=hashlib.sha256(payload).hexdigest(),
        expected_size=len(payload),
    )

    assert resolved == current.resolve()


@pytest.mark.unit
def test_resolver_rejects_unrecorded_evidence() -> None:
    with pytest.raises(BenchmarkError, match="No current or immutable snapshot"):
        resolve_archived_implementation_path(
            PROJECT_ROOT / "renv.lock",
            expected_sha256="0" * 64,
            expected_size=0,
        )


@pytest.mark.unit
@pytest.mark.parametrize("mutation", ["payload", "manifest", "extra"])
def test_snapshot_validation_fails_closed_on_mutation(
    tmp_path: Path,
    mutation: str,
) -> None:
    copied_root = tmp_path / "snapshots"
    shutil.copytree(SNAPSHOT_ROOT, copied_root)
    copied_snapshot = copied_root / SNAPSHOT_ID
    if mutation == "payload":
        (copied_snapshot / "environment" / "verify.R").write_bytes(b"tampered\n")
    elif mutation == "manifest":
        manifest_path = copied_snapshot / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["files"][0]["sha256"] = "0" * 64
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    else:
        (copied_snapshot / "unexpected.txt").write_text(
            "unexpected\n", encoding="utf-8"
        )

    with pytest.raises(BenchmarkError):
        load_environment_snapshots(copied_root)
