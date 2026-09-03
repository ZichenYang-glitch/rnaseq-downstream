"""Focused regressions for stable-ID display annotations in legacy reports."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import types

import pytest

pd = pytest.importorskip(
    "pandas",
    reason="Legacy report integration tests require the optional pandas runtime.",
)
pytest.importorskip(
    "yaml",
    reason="Importing the legacy data module requires the optional PyYAML runtime.",
)


pytestmark = pytest.mark.unit
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _stub_module(name: str, *, package: bool = False) -> types.ModuleType:
    module = types.ModuleType(name)
    if package:
        module.__path__ = []
    return module


@pytest.fixture
def legacy_report(monkeypatch: pytest.MonkeyPatch):
    matplotlib = _stub_module("matplotlib", package=True)
    pyplot = _stub_module("matplotlib.pyplot")
    seaborn = _stub_module("seaborn")
    monkeypatch.setitem(sys.modules, "matplotlib", matplotlib)
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", pyplot)
    monkeypatch.setitem(sys.modules, "seaborn", seaborn)

    spec = importlib.util.spec_from_file_location(
        "checkpoint_a_legacy_report",
        PROJECT_ROOT / "legacy" / "modules" / "report.py",
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_master_table_preserves_count_gene_names_as_display_only(
    legacy_report,
    tmp_path: Path,
) -> None:
    metadata = pd.DataFrame(
        {"group": ["A", "B"]},
        index=pd.Index(["s1", "s2"], name="sample_id"),
    )
    counts = pd.DataFrame(
        {"ENSG2": [2, 4], "ENSG1": [1, 3]},
        index=metadata.index,
    )
    count_annotation = pd.DataFrame(
        {"gene_name": ["SHARED", "SHARED"]},
        index=pd.Index(["ENSG2", "ENSG1"], name="gene_id"),
    )
    out_dir = tmp_path / "report"

    legacy_report.create_master_table(
        metadata,
        "group",
        {},
        tmp_path / "missing-tpm.tsv",
        counts,
        out_dir,
        count_annotation_df=count_annotation,
    )

    table = pd.read_csv(out_dir / "Master_Expression_Table.csv", index_col="gene_id")
    assert table.index.tolist() == ["ENSG2", "ENSG1"]
    assert table["gene_name"].tolist() == ["SHARED", "SHARED"]
    assert table.loc["ENSG2", ["s1", "s2"]].tolist() == [2, 4]


def test_both_legacy_report_callers_forward_embedded_count_annotations() -> None:
    for path in (
        PROJECT_ROOT / "legacy" / "main.py",
        PROJECT_ROOT / "legacy" / "scripts" / "run_report.py",
    ):
        source = path.read_text(encoding="utf-8")
        assert "return_gene_annotations=True" in source
        assert "count_annotation_df=count_annotation_df" in source


def test_contrast_report_omits_retired_motif_card(
    legacy_report,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    report_dir = output_dir / "05_Summary"
    (output_dir / "06_Motif" / "treated_vs_control").mkdir(parents=True)
    (output_dir / "04_GSEA" / "treated_vs_control" / "Hallmark").mkdir(
        parents=True
    )
    results = {
        "treated_vs_control": pd.DataFrame(
            {
                "log2FoldChange": [2.0, -1.0],
                "padj": [0.01, 0.2],
                "pvalue": [0.001, 0.1],
                "stat": [4.0, -1.5],
            },
            index=["ENSG1", "ENSG2"],
        )
    }

    legacy_report.create_contrast_report_pages(output_dir, report_dir, results)

    html = (report_dir / "contrast_treated_vs_control.html").read_text(
        encoding="utf-8"
    )
    assert "GSEA Outputs" in html
    assert "Hallmark" in html
    assert "Motif Outputs" not in html
    assert "Open motif results" not in html
    assert "06_Motif" not in html
