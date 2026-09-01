"""Unit tests for non-overwriting input-validation evidence bundles."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import pytest

from rnaseq_downstream.errors import (
    InternalToolkitError,
    InvalidRequestError,
    OutputWriteError,
)
from rnaseq_downstream import validation_bundle


@pytest.mark.unit
def test_canonical_bundle_json_is_utf8_safe_for_surrogateescaped_text() -> None:
    payload = validation_bundle._canonical_bytes({"path": "bad-\udcff-name"})

    assert payload.decode("utf-8").isascii()
    assert json.loads(payload)["path"] == "bad-\udcff-name"


def _write_bundle(path: Path) -> dict:
    return validation_bundle._write_validation_bundle(
        path,
        input_data={
            "validation_level": "validate",
            "certification_scope": "input_semantics_only",
            "input_semantics": "featurecounts_integer",
            "input_certification_eligible": True,
            "input_certification_status": "passed",
            "request": {"path": "/declared/request.json", "sha256": "a" * 64},
            "route": {"backend_input": "edgeR::DGEList"},
        },
        warnings=[],
        consumed_artifacts=[
            {
                "kind": "consumed_input",
                "role": "request",
                "path": "/declared/request.json",
                "sha256": "a" * 64,
                "size_bytes": 1,
            }
        ],
    )


@pytest.mark.unit
def test_bundle_is_complete_strict_json_and_deterministic(tmp_path: Path) -> None:
    first = _write_bundle(tmp_path / "first")
    second = _write_bundle(tmp_path / "second")

    assert first["plan_id"] == second["plan_id"]
    assert first["validation"] == {
        "input_semantics": "passed",
        "design": "not_run",
        "backend": "not_run",
        "runnable": False,
        "analysis_path_certified": False,
    }
    assert first["publication_status"] == "durability_confirmed"
    assert first["warnings"] == []
    first_digests = {item["role"]: item["sha256"] for item in first["artifacts"]}
    second_digests = {item["role"]: item["sha256"] for item in second["artifacts"]}
    assert first_digests == second_digests

    names = sorted(path.name for path in (tmp_path / "first").iterdir())
    assert names == [
        "bundle_manifest.json",
        "input_plan.json",
        "provenance.json",
        "validated_request.json",
    ]
    observed_plan_ids = set()
    for path in (tmp_path / "first").iterdir():
        text = path.read_text(encoding="utf-8")
        assert text.endswith("\n")
        assert "NaN" not in text
        assert "Infinity" not in text
        document = json.loads(text)
        assert isinstance(document, dict)
        observed_plan_ids.add(document["plan_id"])
    assert observed_plan_ids == {first["plan_id"]}
    for artifact in first["artifacts"]:
        payload = Path(artifact["path"]).read_bytes()
        assert len(payload) == artifact["size_bytes"]
        assert hashlib.sha256(payload).hexdigest() == artifact["sha256"]

    manifest = json.loads(
        (tmp_path / "first" / "bundle_manifest.json").read_text(encoding="utf-8")
    )
    assert [member["path"] for member in manifest["members"]] == [
        "input_plan.json",
        "provenance.json",
        "validated_request.json",
    ]


@pytest.mark.unit
def test_ineligible_validation_bundle_never_claims_input_passed(
    tmp_path: Path,
) -> None:
    result = validation_bundle._write_validation_bundle(
        tmp_path / "ineligible",
        input_data={
            "validation_level": "validate",
            "certification_scope": "input_semantics_only",
            "input_semantics": "salmon_quant_dirs_three_prime",
            "input_certification_eligible": False,
            "input_certification_status": "ineligible",
            "request": {"path": "/declared/request.json", "sha256": "a" * 64},
            "route": {"certified_path_execution_permitted": False},
        },
        warnings=[{"code": "HIGH_RISK_OVERRIDE", "severity": "high"}],
        consumed_artifacts=[],
    )

    assert result["validation"]["input_semantics"] == "ineligible"
    plan = json.loads(
        (tmp_path / "ineligible" / "input_plan.json").read_text(encoding="utf-8")
    )
    assert plan["validation"]["input_semantics"] == "ineligible"


@pytest.mark.unit
def test_bundle_never_overwrites_existing_output(tmp_path: Path) -> None:
    output = tmp_path / "evidence"
    output.mkdir()
    marker = output / "archive.txt"
    marker.write_text("keep", encoding="utf-8")

    with pytest.raises(InvalidRequestError) as caught:
        _write_bundle(output)

    assert caught.value.details["reason"] == "output_exists"
    assert marker.read_text(encoding="utf-8") == "keep"


@pytest.mark.unit
def test_invalid_json_data_creates_no_output(tmp_path: Path) -> None:
    output = tmp_path / "evidence"

    with pytest.raises(InternalToolkitError):
        validation_bundle._write_validation_bundle(
            output,
            input_data={
                "validation_level": "validate",
                "certification_scope": "input_semantics_only",
                "input_semantics": "featurecounts_integer",
                "input_certification_eligible": True,
                "input_certification_status": "passed",
                "request": {},
                "route": {"invalid": float("nan")},
            },
            warnings=[],
            consumed_artifacts=[],
        )

    assert not output.exists()


@pytest.mark.unit
def test_inspection_data_cannot_mint_validated_evidence(tmp_path: Path) -> None:
    output = tmp_path / "evidence"

    with pytest.raises(InternalToolkitError) as caught:
        validation_bundle._write_validation_bundle(
            output,
            input_data={
                "validation_level": "inspect",
                "certification_scope": "input_semantics_only",
                "input_semantics": "featurecounts_integer",
                "input_certification_eligible": None,
                "input_certification_status": "not_evaluated",
                "request": {},
                "route": {},
            },
            warnings=[],
            consumed_artifacts=[],
        )

    assert "validation_level" in caught.value.details["mismatches"]
    assert not output.exists()


@pytest.mark.unit
@pytest.mark.parametrize("live_target", [False, True])
def test_output_symlink_is_rejected_before_resolution(
    tmp_path: Path,
    live_target: bool,
) -> None:
    referent = tmp_path / "referent"
    if live_target:
        referent.mkdir()
    link = tmp_path / "evidence"
    link.symlink_to(referent, target_is_directory=True)

    with pytest.raises(InvalidRequestError) as caught:
        _write_bundle(link)

    assert caught.value.details["reason"] == "output_symlink"
    assert link.is_symlink()
    if live_target:
        assert list(referent.iterdir()) == []
    else:
        assert not referent.exists()


@pytest.mark.unit
def test_atomic_publish_does_not_replace_racing_empty_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "evidence"
    real_publish = validation_bundle._publish_directory_noreplace
    racing_inode: int | None = None

    def create_racing_target(source: Path, target: Path) -> None:
        nonlocal racing_inode
        target.mkdir()
        racing_inode = target.stat().st_ino
        real_publish(source, target)

    monkeypatch.setattr(
        validation_bundle,
        "_publish_directory_noreplace",
        create_racing_target,
    )

    with pytest.raises(InvalidRequestError) as caught:
        _write_bundle(output)

    assert caught.value.details["reason"] == "output_race"
    assert output.is_dir()
    assert output.stat().st_ino == racing_inode
    assert list(output.iterdir()) == []
    assert sorted(path.name for path in tmp_path.iterdir()) == ["evidence"]


@pytest.mark.unit
def test_bundle_synchronizes_staging_and_parent_directories(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_fsync = validation_bundle._fsync_directory
    operations: list[str] = []

    def record_fsync(path: Path, *, operation: str) -> None:
        operations.append(operation)
        real_fsync(path, operation=operation)

    monkeypatch.setattr(validation_bundle, "_fsync_directory", record_fsync)

    _write_bundle(tmp_path / "evidence")

    assert operations == [
        "synchronize_staging_directory",
        "synchronize_output_parent",
    ]


@pytest.mark.unit
def test_parent_sync_failure_returns_published_recovery_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "evidence"
    real_fsync = validation_bundle._fsync_directory

    def fail_parent_sync(path: Path, *, operation: str) -> None:
        if operation == "synchronize_output_parent":
            raise OutputWriteError(
                "simulated post-publish sync failure",
                path=path,
                operation=operation,
                cause=OSError("simulated"),
            )
        real_fsync(path, operation=operation)

    monkeypatch.setattr(validation_bundle, "_fsync_directory", fail_parent_sync)

    result = _write_bundle(output)

    assert result["publication_status"] == "published_durability_unconfirmed"
    assert result["warnings"] == [
        {
            "code": "OUTPUT_DURABILITY_UNCONFIRMED",
            "severity": "high",
            "message": (
                "The complete bundle is visible, but synchronizing its parent "
                "directory failed; crash durability is not confirmed."
            ),
            "details": {
                "output_dir": str(output),
                "plan_id": result["plan_id"],
                "manifest": str(output / "bundle_manifest.json"),
                "operation": "synchronize_output_parent",
                "cause_type": "OSError",
            },
        }
    ]
    assert sorted(path.name for path in output.iterdir()) == [
        "bundle_manifest.json",
        "input_plan.json",
        "provenance.json",
        "validated_request.json",
    ]
    assert all(Path(artifact["path"]).is_file() for artifact in result["artifacts"])


@pytest.mark.unit
def test_failed_staging_write_is_cleaned_without_partial_bundle(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "evidence"
    real_write = validation_bundle._write_json
    calls = 0

    def fail_second_write(path: Path, document: dict) -> tuple[str, int]:
        nonlocal calls
        calls += 1
        if calls == 2:
            cause = OSError("simulated write failure")
            raise OutputWriteError(
                "simulated",
                path=path,
                operation="test",
                cause=cause,
            )
        return real_write(path, document)

    monkeypatch.setattr(validation_bundle, "_write_json", fail_second_write)

    with pytest.raises(OutputWriteError):
        _write_bundle(output)

    assert not output.exists()
    assert list(tmp_path.iterdir()) == []
