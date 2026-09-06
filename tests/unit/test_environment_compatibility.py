"""Tests for expanded-environment compatibility evidence generation."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from scripts.benchmark import build_environment_compatibility_report as compatibility


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = (
    PROJECT_ROOT / "scripts/benchmark/environment-compatibility-report-v1.schema.json"
)
ARCHIVED_COMPATIBILITY_REPORT = (
    PROJECT_ROOT / "tests/compatibility/p1-environment-compatibility-report.json"
)
CERTIFICATION_WORKFLOW = PROJECT_ROOT / ".github/workflows/p0-certification.yml"
ARCHIVED_COMPATIBILITY_SHA256 = (
    "7c46297c8337a1fe3ac66040928c5fbb5ed83b612d0f098fd5f7d2444a8661b3"
)
ARCHIVED_COMPATIBILITY_SIZE = 10824


def _expanded_report_directory(tmp_path: Path) -> tuple[Path, Path]:
    reports = tmp_path / "expanded-reports"
    reports.mkdir()
    rscript = tmp_path / "locked-prefix" / "bin" / "Rscript"
    rscript.parent.mkdir(parents=True)
    rscript.write_bytes(b"relocated locked Rscript test fixture\n")
    rscript_digest = hashlib.sha256(rscript.read_bytes()).hexdigest()
    for specification in compatibility.BENCHMARKS:
        archived_path = specification["archived_path"]
        assert isinstance(archived_path, Path)
        document = compatibility.read_json(archived_path)
        implementation = {item["name"]: item for item in document["implementation"]}
        for name, path in compatibility.ENVIRONMENT_IMPLEMENTATION_PATHS.items():
            implementation[name] = compatibility.file_evidence(path, name=name)
        document["implementation"] = [
            implementation[item["name"]] for item in document["implementation"]
        ]
        document["runtime"]["rscript_sha256"] = rscript_digest
        rebuild_metadata_input = specification["rebuild_metadata_input"]
        if rebuild_metadata_input is not None:
            assert isinstance(rebuild_metadata_input, str)
            document["inputs"] = [
                {
                    "name": item["name"],
                    "sha256": "1" * 64,
                    "size_bytes": item["size_bytes"] + 1,
                }
                if item["name"] == rebuild_metadata_input
                else item
                for item in document["inputs"]
            ]
        filename = specification["filename"]
        assert isinstance(filename, str)
        compatibility.write_json(reports / filename, document)
    return reports, rscript


@pytest.mark.unit
def test_archived_compatibility_report_is_immutable_passing_evidence() -> None:
    payload = ARCHIVED_COMPATIBILITY_REPORT.read_bytes()
    report = compatibility.read_json(ARCHIVED_COMPATIBILITY_REPORT)

    assert len(payload) == ARCHIVED_COMPATIBILITY_SIZE
    assert hashlib.sha256(payload).hexdigest() == ARCHIVED_COMPATIBILITY_SHA256
    compatibility.assert_compatibility_report_shape(report)
    assertions = {item["name"]: item for item in report["assertions"]}
    assert assertions["numeric_artifact_sha256_and_size_identity"]["observed"] == 98
    assert (
        assertions["installed_description_drift_is_source_archive_bound"]["observed"]
        == 2
    )
    assert (
        assertions["rscript_relocation_is_conda_lock_bound"]["observed"]
        == "d56eb9cb63c0c727c0b3bd025d6c5559cf1362273a6b913bf4498a0bd5e5636c"
    )
    runtime_changes = [
        record["runtime_metadata_changes"][0] for record in report["benchmarks"]
    ]
    assert {change["before"] for change in runtime_changes} == {
        "a09d33dc7b786eb86191e428b3f75406d3a25c5b8fdf8c8e4282d0e55451c331"
    }
    assert {change["after"] for change in runtime_changes} == {
        "d56eb9cb63c0c727c0b3bd025d6c5559cf1362273a6b913bf4498a0bd5e5636c"
    }
    assert {change["reason"] for change in runtime_changes} == {
        "locked_conda_prefix_relocation"
    }
    assert {change["conda_lock"]["sha256"] for change in runtime_changes} == {
        "7cddb25b11bdd7d31a5dee156cfd2f8c456398005f3b2acc6c23d5df99e5b702"
    }


@pytest.mark.unit
def test_archived_compatibility_environment_resolves_after_a_later_evolution(
    tmp_path: Path,
) -> None:
    report = compatibility.read_json(ARCHIVED_COMPATIBILITY_REPORT)
    fake_project = tmp_path / "future-project"
    fake_snapshot_root = fake_project / "environment" / "snapshots"
    shutil.copytree(compatibility.SNAPSHOT_ROOT, fake_snapshot_root)

    environment_paths: dict[str, Path] = {}
    source_paths = {
        "renv.lock": "renv.lock",
        "r-sources.lock": "environment/r-sources.lock",
        "verify.R": "environment/verify.R",
    }
    p1_snapshot = fake_snapshot_root / "p1-environment-expanded"
    manifest_files: list[dict[str, object]] = []
    for name, source_path in source_paths.items():
        current = compatibility.ENVIRONMENT_IMPLEMENTATION_PATHS[name]
        snapshot_payload = p1_snapshot / source_path
        snapshot_payload.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(current, snapshot_payload)
        evidence = compatibility.file_evidence(current, name=source_path)
        manifest_files.append(
            {
                "source_path": source_path,
                "snapshot_path": source_path,
                "sha256": evidence["sha256"],
                "size_bytes": evidence["size_bytes"],
            }
        )

        future_current = fake_project / source_path
        future_current.parent.mkdir(parents=True, exist_ok=True)
        future_current.write_text(
            f"future environment bytes for {name}\n", encoding="utf-8"
        )
        environment_paths[name] = future_current

    compatibility.write_json(
        p1_snapshot / "manifest.json",
        {
            "schema_version": "rnaseq-downstream-environment-snapshot-v1",
            "snapshot_id": p1_snapshot.name,
            "source_revision": "1" * 40,
            "immutability": "append_only",
            "files": manifest_files,
        },
    )

    compatibility.assert_compatibility_report_shape(
        report,
        environment_paths=environment_paths,
        project_root=fake_project,
        snapshot_root=fake_snapshot_root,
    )


@pytest.mark.unit
def test_compatibility_builder_binds_four_real_report_shapes(tmp_path: Path) -> None:
    reports, rscript = _expanded_report_directory(tmp_path)
    report = compatibility.build_compatibility_report(
        reports,
        rscript_path=rscript,
    )

    assert report["schema_version"] == compatibility.SCHEMA_VERSION
    assert report["compatibility_id"] == compatibility.COMPATIBILITY_ID
    assert report["status"] == "pass"
    assert report["permitted_implementation_changes"] == [
        "r-sources.lock",
        "renv.lock",
        "verify.R",
    ]
    assert len(report["benchmarks"]) == 4
    assert all(
        len(record["implementation_changes"]) == 3 for record in report["benchmarks"]
    )
    assert all(
        len(record["runtime_metadata_changes"]) == 1 for record in report["benchmarks"]
    )
    assert (
        sum(len(record["rebuild_metadata_changes"]) for record in report["benchmarks"])
        == 2
    )
    assert (
        sum(record["artifact_inventory"]["count"] for record in report["benchmarks"])
        > 0
    )
    assert all(assertion["passed"] for assertion in report["assertions"])


@pytest.mark.unit
def test_compatibility_builder_rejects_non_environment_drift(tmp_path: Path) -> None:
    reports, rscript = _expanded_report_directory(tmp_path)
    specification = compatibility.BENCHMARKS[0]
    filename = specification["filename"]
    assert isinstance(filename, str)
    path = reports / filename
    document = compatibility.read_json(path)
    document["metrics"] = copy.deepcopy(document["metrics"])
    document["metrics"]["scope"] = "tampered"
    compatibility.write_json(path, document)

    with pytest.raises(compatibility.BenchmarkError, match="changed outside"):
        compatibility.build_compatibility_report(reports, rscript_path=rscript)


@pytest.mark.unit
def test_compatibility_builder_rejects_numeric_artifact_drift(tmp_path: Path) -> None:
    reports, rscript = _expanded_report_directory(tmp_path)
    specification = compatibility.BENCHMARKS[0]
    filename = specification["filename"]
    assert isinstance(filename, str)
    path = reports / filename
    document = compatibility.read_json(path)
    document["artifacts"][0]["sha256"] = "0" * 64
    compatibility.write_json(path, document)

    with pytest.raises(
        compatibility.BenchmarkError,
        match="numeric artifact SHA-256/size inventory changed",
    ):
        compatibility.build_compatibility_report(reports, rscript_path=rscript)


@pytest.mark.unit
@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("other_runtime_field", "runtime changes are not the approved Rscript"),
        ("wrong_rscript_digest", "does not match the live locked runtime"),
    ],
)
def test_compatibility_builder_rejects_unapproved_runtime_drift(
    tmp_path: Path,
    mutation: str,
    message: str,
) -> None:
    reports, rscript = _expanded_report_directory(tmp_path)
    specification = compatibility.BENCHMARKS[0]
    filename = specification["filename"]
    assert isinstance(filename, str)
    path = reports / filename
    document = compatibility.read_json(path)
    if mutation == "other_runtime_field":
        document["runtime"]["python"] = "0.0.0"
    else:
        document["runtime"]["rscript_sha256"] = "2" * 64
    compatibility.write_json(path, document)

    with pytest.raises(compatibility.BenchmarkError, match=message):
        compatibility.build_compatibility_report(reports, rscript_path=rscript)


@pytest.mark.unit
def test_compatibility_builder_rejects_a_stale_environment_record(
    tmp_path: Path,
) -> None:
    reports, rscript = _expanded_report_directory(tmp_path)
    specification = compatibility.BENCHMARKS[0]
    filename = specification["filename"]
    archived_path = specification["archived_path"]
    assert isinstance(filename, str)
    assert isinstance(archived_path, Path)
    path = reports / filename
    document = compatibility.read_json(path)
    archived = compatibility.read_json(archived_path)
    archived_by_name = {item["name"]: item for item in archived["implementation"]}
    document["implementation"] = [
        archived_by_name["verify.R"] if item["name"] == "verify.R" else item
        for item in document["implementation"]
    ]
    compatibility.write_json(path, document)

    with pytest.raises(
        compatibility.BenchmarkError,
        match="implementation changes are not environment-only",
    ):
        compatibility.build_compatibility_report(reports, rscript_path=rscript)


@pytest.mark.unit
def test_compatibility_schema_is_strict_and_versioned() -> None:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))

    assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    assert schema["additionalProperties"] is False
    assert schema["properties"]["schema_version"]["const"] == (
        compatibility.SCHEMA_VERSION
    )
    assert schema["properties"]["compatibility_id"]["const"] == (
        compatibility.COMPATIBILITY_ID
    )
    assert schema["properties"]["benchmarks"]["minItems"] == 4
    assert schema["properties"]["benchmarks"]["maxItems"] == 4
    assert schema["properties"]["assertions"]["minItems"] == 6
    assert schema["properties"]["assertions"]["maxItems"] == 6
    runtime_changes = schema["properties"]["benchmarks"]["items"]["properties"][
        "runtime_metadata_changes"
    ]
    assert runtime_changes["minItems"] == 1
    assert runtime_changes["maxItems"] == 1


@pytest.mark.unit
def test_ci_rebuilds_and_archives_live_compatibility_evidence() -> None:
    workflow = CERTIFICATION_WORKFLOW.read_text(encoding="utf-8")

    assert "build_environment_compatibility_report.py" in workflow
    assert "p1-environment-compatibility-report.json" in workflow
    assert "path: benchmark-results/*.json" in workflow
    assert "cmp --silent" not in workflow
