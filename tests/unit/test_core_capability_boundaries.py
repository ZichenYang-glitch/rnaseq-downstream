"""Regression tests for explicitly excluded checkpoint-A core capabilities."""

from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.unit
def test_experimental_workflow_is_grouped_under_legacy() -> None:
    legacy_root = PROJECT_ROOT / "legacy"

    for old_root_source in (
        "config.py",
        "environment.yaml",
        "main.py",
        "modules/data.py",
        "workflow/rules/rnaseq.smk",
        "workflow_config.yaml",
    ):
        assert not (PROJECT_ROOT / old_root_source).exists()

    for required_entry in (
        "config.py",
        "environment.yaml",
        "main.py",
        "modules",
        "workflow",
        "workflow_config.yaml",
    ):
        assert (legacy_root / required_entry).exists()

    root_snakefile = (PROJECT_ROOT / "Snakefile").read_text(encoding="utf-8")
    assert 'include: "legacy/Snakefile"' in root_snakefile

    legacy_snakefile = (legacy_root / "Snakefile").read_text(encoding="utf-8")
    assert 'configfile: "legacy/workflow_config.yaml"' in legacy_snakefile
    assert 'include: "workflow/rules/rnaseq.smk"' in legacy_snakefile


@pytest.mark.unit
def test_homer_is_retained_only_as_a_noncore_extra() -> None:
    legacy_root = PROJECT_ROOT / "legacy"
    assert not (legacy_root / "modules" / "motif.py").exists()
    assert not (legacy_root / "scripts" / "run_motif.py").exists()
    assert (legacy_root / "extras" / "homer" / "motif.py").is_file()
    assert (legacy_root / "extras" / "homer" / "run_motif.py").is_file()

    main_source = (legacy_root / "main.py").read_text(encoding="utf-8")
    workflow_source = (
        legacy_root / "workflow" / "rules" / "rnaseq.smk"
    ).read_text(encoding="utf-8")
    config_source = (legacy_root / "workflow_config.yaml").read_text(encoding="utf-8")

    assert "motif" not in main_source.lower()
    assert "rule motif" not in workflow_source.lower()
    assert "run_motif" not in workflow_source.lower()
    assert "run_motif" not in config_source.lower()
    assert "homer" not in config_source.lower()
