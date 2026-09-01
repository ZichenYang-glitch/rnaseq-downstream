"""Regression tests for explicitly excluded checkpoint-A core capabilities."""

from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.unit
def test_homer_is_retained_only_as_a_noncore_extra() -> None:
    assert not (PROJECT_ROOT / "modules" / "motif.py").exists()
    assert not (PROJECT_ROOT / "scripts" / "run_motif.py").exists()
    assert (PROJECT_ROOT / "extras" / "homer" / "motif.py").is_file()
    assert (PROJECT_ROOT / "extras" / "homer" / "run_motif.py").is_file()

    main_source = (PROJECT_ROOT / "main.py").read_text(encoding="utf-8")
    workflow_source = (PROJECT_ROOT / "workflow" / "rules" / "rnaseq.smk").read_text(
        encoding="utf-8"
    )
    config_source = (PROJECT_ROOT / "workflow_config.yaml").read_text(encoding="utf-8")

    assert "motif" not in main_source.lower()
    assert "rule motif" not in workflow_source.lower()
    assert "run_motif" not in workflow_source.lower()
    assert "run_motif" not in config_source.lower()
    assert "homer" not in config_source.lower()
