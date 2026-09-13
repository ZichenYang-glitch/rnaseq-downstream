"""Same-content source relocation must preserve historical evidence identity."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from scripts.benchmark import build_environment_compatibility_report as compatibility
from scripts.benchmark import evidence_resolver as resolver


pytestmark = pytest.mark.unit
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT_ID = "p1-source-urls-0dba865"
SNAPSHOT = resolver.SNAPSHOT_ROOT / SNAPSHOT_ID
SOURCE_PATH = "environment/r-sources.lock"
SOURCE_SHA256 = "cbe8cdc6dae73a629dadf6a6c2e743f0fb95c20c9b85bb9acecbe23c64bd8659"
SOURCE_SIZE = 1676
MANIFEST_SHA256 = "8da8a8f69ce5f588deee4ccc66fbc436ade48a9c574bd5c90a1f74fb46a5ae78"
REVISION = "0dba86589365fcb366c98a84c3bc812a6beab21c"
COMPCODER_SHA256 = "9890c63d8f6cb585ef9311fa888d162ebe1148809c9082d8710f39a07a013b07"
RELEASE_URL = (
    "https://bioconductor.org/packages/3.23/bioc/src/contrib/compcodeR_1.48.0.tar.gz"
)
ARCHIVE_URL = (
    "https://bioconductor.org/packages/3.23/bioc/src/contrib/Archive/"
    "compcodeR/compcodeR_1.48.0.tar.gz"
)


def _assert_snapshot_identity(snapshot: Path) -> None:
    manifest_bytes = (snapshot / "manifest.json").read_bytes()
    payload = (snapshot / SOURCE_PATH).read_bytes()
    assert len(manifest_bytes) == 452
    assert hashlib.sha256(manifest_bytes).hexdigest() == MANIFEST_SHA256
    assert len(payload) == SOURCE_SIZE
    assert hashlib.sha256(payload).hexdigest() == SOURCE_SHA256


def _write_records(path: Path, records: dict[str, dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["package", "version", "repository", "role", "url", "sha256"],
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records.values())


def test_pre_migration_snapshot_has_frozen_bytes_and_source_revision() -> None:
    _assert_snapshot_identity(SNAPSHOT)
    records = [
        record
        for record in resolver.load_environment_snapshots()
        if record.snapshot_id == SNAPSHOT_ID
    ]
    assert len(records) == 1
    record = records[0]
    assert record.source_revision == REVISION
    assert record.source_path == record.snapshot_path == SOURCE_PATH
    assert record.sha256 == SOURCE_SHA256
    assert record.size_bytes == SOURCE_SIZE


def test_source_lock_change_is_only_the_reviewed_same_hash_url_relocation() -> None:
    before = compatibility._source_manifest_records(SNAPSHOT / SOURCE_PATH)
    after = compatibility._source_manifest_records(PROJECT_ROOT / SOURCE_PATH)
    assert len(before) == len(after) == 9
    assert before["compcodeR"]["url"] == RELEASE_URL
    assert after["compcodeR"]["url"] == ARCHIVE_URL
    assert before["compcodeR"]["sha256"] == COMPCODER_SHA256
    before["compcodeR"]["url"] = ARCHIVE_URL
    assert before == after


def test_source_evidence_resolves_current_before_historical_snapshot() -> None:
    current = PROJECT_ROOT / SOURCE_PATH
    assert (
        resolver.resolve_archived_implementation_path(
            current,
            expected_sha256=hashlib.sha256(current.read_bytes()).hexdigest(),
            expected_size=current.stat().st_size,
            snapshot_root=PROJECT_ROOT / "does-not-exist",
        )
        == current
    )
    assert (
        resolver.resolve_archived_implementation_path(
            current, expected_sha256=SOURCE_SHA256, expected_size=SOURCE_SIZE
        )
        == SNAPSHOT / SOURCE_PATH
    )


@pytest.mark.parametrize("mutation", ["payload", "manifest", "unlisted_file"])
def test_corrupted_copied_snapshot_fails_closed(tmp_path: Path, mutation: str) -> None:
    copied = tmp_path / "snapshots" / SNAPSHOT_ID
    shutil.copytree(SNAPSHOT, copied)
    if mutation == "payload":
        with (copied / SOURCE_PATH).open("ab") as handle:
            handle.write(b"modified\n")
    elif mutation == "manifest":
        manifest = json.loads((copied / "manifest.json").read_text())
        manifest["files"][0]["sha256"] = "0" * 64
        compatibility.write_json(copied / "manifest.json", manifest)
    else:
        (copied / "unexpected.txt").write_text("unregistered payload\n")
    with pytest.raises(compatibility.BenchmarkError):
        resolver.load_environment_snapshots(copied.parent)
    with pytest.raises(compatibility.BenchmarkError):
        resolver.resolve_archived_implementation_path(
            PROJECT_ROOT / SOURCE_PATH,
            expected_sha256=SOURCE_SHA256,
            expected_size=SOURCE_SIZE,
            snapshot_root=copied.parent,
        )


def test_self_consistently_rewritten_snapshot_still_violates_frozen_identity(
    tmp_path: Path,
) -> None:
    copied = tmp_path / SNAPSHOT_ID
    shutil.copytree(SNAPSHOT, copied)
    payload = copied / SOURCE_PATH
    payload.write_bytes((PROJECT_ROOT / SOURCE_PATH).read_bytes())
    manifest = json.loads((copied / "manifest.json").read_text())
    manifest["files"][0]["sha256"] = hashlib.sha256(payload.read_bytes()).hexdigest()
    manifest["files"][0]["size_bytes"] = payload.stat().st_size
    compatibility.write_json(copied / "manifest.json", manifest)
    with pytest.raises(AssertionError):
        _assert_snapshot_identity(copied)


@pytest.mark.parametrize("url", [RELEASE_URL, ARCHIVE_URL])
def test_compatibility_accepts_unchanged_source_bytes_at_permitted_locations(
    tmp_path: Path, url: str
) -> None:
    records = compatibility._source_manifest_records(SNAPSHOT / SOURCE_PATH)
    records["compcodeR"]["url"] = url
    source_lock = tmp_path / "r-sources.lock"
    _write_records(source_lock, records)
    result = compatibility._unchanged_source_archives(expanded_source_lock=source_lock)
    assert result["compcodeR"] == records["compcodeR"]
    assert result["compcodeR"]["sha256"] == COMPCODER_SHA256
    assert result["airway"] == records["airway"]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("sha256", "0" * 64),
        ("sha256", ""),
        ("version", "1.48.1"),
        ("repository", "unapproved repository"),
        ("role", "unapproved role"),
        ("package", "renamed_compcodeR"),
        ("url", ARCHIVE_URL.replace("bioconductor.org", "unofficial.example")),
        ("url", ARCHIVE_URL.replace("bioconductor.org", "bioconductor.org.evil.test")),
        ("url", ARCHIVE_URL.replace("Archive/compcodeR/", "Archive/")),
        ("url", ARCHIVE_URL.replace("Archive/compcodeR/", "Archive/other/")),
        ("url", ARCHIVE_URL.replace("1.48.0.tar.gz", "1.48.1.tar.gz")),
        ("url", ARCHIVE_URL.replace("compcodeR_", "other_")),
        ("url", ARCHIVE_URL + "?download=1"),
        ("url", ARCHIVE_URL + "#fragment"),
        ("url", ARCHIVE_URL.replace("https://", "http://")),
        ("url", ARCHIVE_URL.replace("/3.23/", "/3.24/")),
        ("url", ARCHIVE_URL.replace("/Archive/", "/archive/")),
        ("url", ARCHIVE_URL.replace("/Archive/", "/./Archive/")),
    ],
)
def test_archive_relocation_cannot_authorize_other_source_changes(
    tmp_path: Path, field: str, value: str
) -> None:
    records = compatibility._source_manifest_records(SNAPSHOT / SOURCE_PATH)
    records["compcodeR"]["url"] = ARCHIVE_URL
    records["compcodeR"][field] = value
    source_lock = tmp_path / "r-sources.lock"
    _write_records(source_lock, records)
    with pytest.raises(compatibility.BenchmarkError):
        compatibility._unchanged_source_archives(expanded_source_lock=source_lock)


def test_reverse_archive_to_release_migration_is_not_implicitly_authorized(
    tmp_path: Path,
) -> None:
    copied_root = tmp_path / "snapshots"
    baseline = copied_root / compatibility.BASELINE_SNAPSHOT_ID
    shutil.copytree(
        resolver.SNAPSHOT_ROOT / compatibility.BASELINE_SNAPSHOT_ID, baseline
    )
    records = compatibility._source_manifest_records(baseline / SOURCE_PATH)
    records["compcodeR"]["url"] = ARCHIVE_URL
    _write_records(baseline / SOURCE_PATH, records)
    manifest = json.loads((baseline / "manifest.json").read_text())
    for item in manifest["files"]:
        if item["source_path"] == SOURCE_PATH:
            payload = (baseline / SOURCE_PATH).read_bytes()
            item["sha256"] = hashlib.sha256(payload).hexdigest()
            item["size_bytes"] = len(payload)
    compatibility.write_json(baseline / "manifest.json", manifest)
    records["compcodeR"]["url"] = RELEASE_URL
    current = tmp_path / "r-sources.lock"
    _write_records(current, records)
    with pytest.raises(
        compatibility.BenchmarkError, match="official Archive relocation"
    ):
        compatibility._unchanged_source_archives(
            expanded_source_lock=current, snapshot_root=copied_root
        )


def test_official_cran_relocation_uses_same_package_version_and_https_host() -> None:
    before = compatibility._source_manifest_records(SNAPSHOT / SOURCE_PATH)["renv"]
    after = dict(
        before,
        url="https://cran.r-project.org/src/contrib/Archive/renv/renv_1.2.4.tar.gz",
    )
    assert compatibility._is_official_archive_relocation(before, after)
    for url in (
        after["url"].replace("cran.r-project.org", "cloud.r-project.org"),
        after["url"].replace("https://", "http://"),
        after["url"] + "?download=1",
        after["url"].replace("1.2.4", "1.2.5"),
    ):
        assert not compatibility._is_official_archive_relocation(
            before, dict(after, url=url)
        )


def test_url_audit_covers_every_locked_source_and_proves_only_migrated_hash() -> None:
    report = json.loads(
        (
            PROJECT_ROOT / "environment/source-audits/2026-09-13-url-audit.json"
        ).read_text()
    )
    before = compatibility._source_manifest_records(SNAPSHOT / SOURCE_PATH)
    after = compatibility._source_manifest_records(PROJECT_ROOT / SOURCE_PATH)
    assert report["status"] == "pass"
    assert report["source_revision"] == REVISION
    assert report["source_lock_before_sha256"] == SOURCE_SHA256
    assert report["package_count"] == 9
    assert report["unchanged_url_count"] == 8
    assert report["same_hash_migration_count"] == 1
    assert len(report["packages"]) == len(before)
    assert {item["package"] for item in report["packages"]} == set(before)
    for item in report["packages"]:
        old = before[item["package"]]
        new = after[item["package"]]
        assert item["version"] == old["version"] == new["version"]
        assert item["repository"] == old["repository"] == new["repository"]
        assert item["old_url"] == old["url"]
        assert item["locked_sha256"] == old["sha256"] == new["sha256"]
        if item["package"] == "compcodeR":
            assert item["http_status"] == 404
            assert item["candidate_status"] == item["download_http_status"] == 200
            assert item["new_url"] == new["url"] == ARCHIVE_URL
            assert item["hash_check"] == "match"
            assert item["downloaded_sha256"] == COMPCODER_SHA256
            assert item["downloaded_size_bytes"] == 4973426
            assert item["action"] == "migrate_to_official_same_version_archive"
        else:
            assert item["http_status"] == 200
            assert item["new_url"] is None
            assert item["old_url"] == new["url"]
            assert item["hash_check"] == "not_performed_unchanged"
            assert item["downloaded_sha256"] is None
            assert item["action"] == "keep_current_url"
