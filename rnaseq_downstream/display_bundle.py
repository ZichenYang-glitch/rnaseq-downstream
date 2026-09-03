"""Deterministic, inference-free display bundle generation and verification.

The R backend owns the route-aware logCPM transformation.  This module only
selects variable genes, runs the already-tested centered/unscaled PCA helper,
and renders static SVG views from verified numeric artifacts.  Nothing created
here is consumed by differential-expression inference.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
from html import escape
import io
import json
import math
import os
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from .errors import (
    BackendFailedError,
    InputIntegrityError,
    OutputWriteError,
    QCValidationError,
    ToolkitError,
)

DISPLAY_SCHEMA_VERSION = "1.0"
DISPLAY_RENDERER_ID = "rnaseq-downstream-static-svg-v1"
SVG_WIDTH = 960
SVG_HEIGHT = 640
FDR_FLOOR = 1e-300

_CORE_FILES = (
    "analysis.json",
    "backend_manifest.json",
    "coefficients.tsv",
    "design.tsv",
    "results.tsv",
)
_RESULT_HEADER = (
    "gene_id",
    "contrast_id",
    "status",
    "logFC",
    "unshrunk_logFC",
    "logCPM",
    "statistic",
    "statistic_type",
    "statistic_status",
    "PValue",
    "FDR",
    "test_method",
    "lfc_threshold",
)
_DECIMAL = re.compile(r"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$")


@dataclass(frozen=True)
class _LogCpm:
    gene_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    values: Any


@dataclass(frozen=True)
class _ContrastPoints:
    contrast_id: str
    lfc_threshold: float
    logfc: Any
    logcpm: Any
    fdr: Any
    excluded_gene_count: int


def _qc_error(message: str, *, reason: str, **details: Any) -> None:
    raise QCValidationError(message, details={"reason": reason, **details})


def _capture(path: Path) -> tuple[bytes, str, int]:
    try:
        if path.is_symlink():
            raise OSError("symbolic links are not accepted")
        resolved = path.resolve(strict=True)
        if not resolved.is_file():
            raise IsADirectoryError(str(resolved))
        before = resolved.stat()
        content = resolved.read_bytes()
        after = resolved.stat()
    except OSError as error:
        raise BackendFailedError(
            "A display source artifact could not be captured.",
            details={"path": str(path), "cause_type": type(error).__name__},
            cause=error,
        ) from error
    before_identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    after_identity = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if before_identity != after_identity or len(content) != after.st_size:
        raise BackendFailedError(
            "A display source artifact changed while it was captured.",
            details={"path": str(path)},
        )
    return content, sha256(content).hexdigest(), len(content)


def _strict_rows(content: bytes, *, role: str) -> list[list[str]]:
    try:
        text = content.decode("utf-8")
        rows = list(
            csv.reader(
                io.StringIO(text, newline=""),
                delimiter="\t",
                quoting=csv.QUOTE_NONE,
                strict=True,
            )
        )
    except (UnicodeError, csv.Error) as error:
        raise BackendFailedError(
            "A display source is not strict UTF-8 TSV.",
            details={"role": role, "cause_type": type(error).__name__},
            cause=error,
        ) from error
    if not rows or not rows[0]:
        raise BackendFailedError(
            "A display source TSV is empty.", details={"role": role}
        )
    width = len(rows[0])
    for number, row in enumerate(rows[1:], start=2):
        if len(row) != width:
            raise BackendFailedError(
                "A display source TSV row has an incompatible width.",
                details={"role": role, "row": number},
            )
    return rows


def _number(raw: str, *, role: str, row: int, field: str) -> float:
    if _DECIMAL.fullmatch(raw) is None:
        raise BackendFailedError(
            "A display source contains an invalid numeric field.",
            details={"role": role, "row": row, "field": field},
        )
    value = float(raw)
    if not math.isfinite(value):
        raise BackendFailedError(
            "A display source contains a non-finite numeric field.",
            details={"role": role, "row": row, "field": field},
        )
    return value


def _parse_logcpm(path: Path) -> tuple[_LogCpm, dict[str, Any]]:
    try:
        import numpy as np
    except ImportError as error:
        _qc_error(
            "Display generation requires the locked NumPy runtime.",
            reason="numpy_unavailable",
            cause_type=type(error).__name__,
        )
    content, digest, size = _capture(path)
    rows = _strict_rows(content, role="display_logcpm")
    header = rows[0]
    if len(header) < 3 or header[0] != "gene_id":
        raise BackendFailedError(
            "The display logCPM header is incompatible.",
            details={"observed_header": header},
        )
    sample_ids = tuple(header[1:])
    if any(not item for item in sample_ids) or len(set(sample_ids)) != len(sample_ids):
        raise BackendFailedError("The display logCPM sample identifiers are invalid.")
    gene_ids: list[str] = []
    observed_gene_ids: set[str] = set()
    values: list[list[float]] = []
    for row_number, row in enumerate(rows[1:], start=2):
        gene_id = row[0]
        if not gene_id or gene_id in observed_gene_ids:
            raise BackendFailedError(
                "The display logCPM gene identifiers are invalid.",
                details={"row": row_number, "gene_id": gene_id},
            )
        gene_ids.append(gene_id)
        observed_gene_ids.add(gene_id)
        values.append(
            [
                _number(raw, role="display_logcpm", row=row_number, field=sample)
                for sample, raw in zip(sample_ids, row[1:], strict=True)
            ]
        )
    if not gene_ids:
        raise BackendFailedError("The display logCPM matrix contains no genes.")
    matrix = np.asarray(values, dtype=float)
    return (
        _LogCpm(tuple(gene_ids), sample_ids, matrix),
        {
            "path": "logcpm.tsv",
            "sha256": digest,
            "size_bytes": size,
        },
    )


def _parse_results(
    path: Path, contrast_order: Sequence[Mapping[str, Any]]
) -> tuple[list[_ContrastPoints], dict[str, Any], set[str]]:
    try:
        import numpy as np
    except ImportError as error:
        _qc_error(
            "Display generation requires the locked NumPy runtime.",
            reason="numpy_unavailable",
            cause_type=type(error).__name__,
        )
    content, digest, size = _capture(path)
    rows = _strict_rows(content, role="results.tsv")
    if tuple(rows[0]) != _RESULT_HEADER:
        raise BackendFailedError(
            "The results table cannot be rendered because its schema is incompatible."
        )
    expected = [str(item["contrast_id"]) for item in contrast_order]
    thresholds = {
        str(item["contrast_id"]): float(item["lfc_threshold"])
        for item in contrast_order
    }
    by_contrast: dict[str, list[tuple[float, float, float]]] = {
        identifier: [] for identifier in expected
    }
    excluded_counts = {identifier: 0 for identifier in expected}
    tested_genes: dict[str, set[str]] = {identifier: set() for identifier in expected}
    for row_number, row in enumerate(rows[1:], start=2):
        values = dict(zip(_RESULT_HEADER, row, strict=True))
        identifier = values["contrast_id"]
        if identifier not in by_contrast:
            raise BackendFailedError(
                "The results table contains an undeclared display contrast.",
                details={"contrast_id": identifier},
            )
        if values["status"] != "tested":
            excluded_counts[identifier] += 1
            continue
        gene_id = values["gene_id"]
        if not gene_id or gene_id in tested_genes[identifier]:
            raise BackendFailedError(
                "The tested display gene inventory is invalid.",
                details={"contrast_id": identifier, "gene_id": gene_id},
            )
        tested_genes[identifier].add(gene_id)
        logfc = _number(
            values["logFC"], role="results.tsv", row=row_number, field="logFC"
        )
        logcpm = _number(
            values["logCPM"], role="results.tsv", row=row_number, field="logCPM"
        )
        fdr = _number(values["FDR"], role="results.tsv", row=row_number, field="FDR")
        if not 0.0 <= fdr <= 1.0:
            raise BackendFailedError(
                "A display FDR is outside [0, 1].",
                details={"contrast_id": identifier, "row": row_number},
            )
        by_contrast[identifier].append((logfc, logcpm, fdr))
    inventories = [tested_genes[identifier] for identifier in expected]
    if (
        not inventories
        or not inventories[0]
        or any(inventory != inventories[0] for inventory in inventories[1:])
    ):
        raise BackendFailedError(
            "The tested gene inventory is not identical across display contrasts."
        )
    points = []
    for identifier in expected:
        records = by_contrast[identifier]
        points.append(
            _ContrastPoints(
                contrast_id=identifier,
                lfc_threshold=thresholds[identifier],
                logfc=np.asarray([item[0] for item in records], dtype=float),
                logcpm=np.asarray([item[1] for item in records], dtype=float),
                fdr=np.asarray([item[2] for item in records], dtype=float),
                excluded_gene_count=excluded_counts[identifier],
            )
        )
    return (
        points,
        {"path": "../results.tsv", "sha256": digest, "size_bytes": size},
        inventories[0],
    )


def _write_bytes(path: Path, payload: bytes) -> dict[str, Any]:
    try:
        with path.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        path.chmod(0o400)
    except OSError as error:
        raise OutputWriteError(
            "A display artifact could not be written.",
            path=path,
            operation="write_display_artifact",
            cause=error,
        ) from error
    return {
        "path": path.name,
        "sha256": sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


def _write_tsv(
    path: Path, header: Sequence[str], rows: Sequence[Sequence[Any]]
) -> dict[str, Any]:
    buffer = io.StringIO(newline="")
    writer = csv.writer(
        buffer,
        delimiter="\t",
        lineterminator="\n",
        quoting=csv.QUOTE_NONE,
        escapechar=None,
    )
    writer.writerow(header)
    writer.writerows(rows)
    return _write_bytes(path, buffer.getvalue().encode("utf-8"))


def _bounds(values: Any) -> tuple[float, float]:
    minimum = float(values.min())
    maximum = float(values.max())
    if not math.isfinite(minimum) or not math.isfinite(maximum):
        _qc_error(
            "Display coordinates must have finite bounds.",
            reason="nonfinite_display_bounds",
        )
    if minimum == maximum:
        delta = max(0.5, abs(minimum) * 0.05)
        lower, upper = minimum - delta, maximum + delta
    else:
        span = maximum - minimum
        if not math.isfinite(span):
            _qc_error(
                "The display coordinate range cannot be represented safely.",
                reason="unrepresentable_display_range",
            )
        padding = span * 0.06
        lower, upper = minimum - padding, maximum + padding
    if not math.isfinite(lower) or not math.isfinite(upper) or lower >= upper:
        _qc_error(
            "The padded display coordinate range cannot be represented safely.",
            reason="unrepresentable_display_range",
        )
    return lower, upper


def _format_tick(value: float) -> str:
    return f"{value:.4g}"


def _scatter_svg(
    *,
    title: str,
    x_label: str,
    y_label: str,
    x: Any,
    y: Any,
    highlighted: Any,
    vertical_lines: Sequence[float] = (),
    horizontal_lines: Sequence[float] = (),
) -> bytes:
    x_min, x_max = _bounds(x)
    y_min, y_max = _bounds(y)
    left, right, top, bottom = 92.0, 32.0, 66.0, 76.0
    plot_width = SVG_WIDTH - left - right
    plot_height = SVG_HEIGHT - top - bottom

    def map_x(value: float) -> float:
        mapped = left + (float(value) - x_min) * plot_width / (x_max - x_min)
        if not math.isfinite(mapped):
            _qc_error(
                "An SVG x coordinate cannot be represented safely.",
                reason="nonfinite_svg_coordinate",
            )
        return mapped

    def map_y(value: float) -> float:
        mapped = top + (y_max - float(value)) * plot_height / (y_max - y_min)
        if not math.isfinite(mapped):
            _qc_error(
                "An SVG y coordinate cannot be represented safely.",
                reason="nonfinite_svg_coordinate",
            )
        return mapped

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{SVG_WIDTH}" '
            f'height="{SVG_HEIGHT}" viewBox="0 0 {SVG_WIDTH} {SVG_HEIGHT}" '
            'role="img">'
        ),
        f"<title>{escape(title)}</title>",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{SVG_WIDTH / 2:.1f}" y="34" text-anchor="middle" '
        f'font-family="sans-serif" font-size="20" fill="#16202a">{escape(title)}</text>',
    ]
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        x_value = x_min + fraction * (x_max - x_min)
        x_pixel = map_x(x_value)
        y_value = y_min + fraction * (y_max - y_min)
        y_pixel = map_y(y_value)
        parts.extend(
            [
                f'<line x1="{x_pixel:.3f}" y1="{top:.1f}" x2="{x_pixel:.3f}" y2="{top + plot_height:.1f}" stroke="#e7ebef" stroke-width="1"/>',
                f'<text x="{x_pixel:.3f}" y="{top + plot_height + 24:.1f}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#52606d">{escape(_format_tick(x_value))}</text>',
                f'<line x1="{left:.1f}" y1="{y_pixel:.3f}" x2="{left + plot_width:.1f}" y2="{y_pixel:.3f}" stroke="#e7ebef" stroke-width="1"/>',
                f'<text x="{left - 12:.1f}" y="{y_pixel + 4:.3f}" text-anchor="end" font-family="sans-serif" font-size="12" fill="#52606d">{escape(_format_tick(y_value))}</text>',
            ]
        )
    for value in vertical_lines:
        if x_min <= value <= x_max:
            pixel = map_x(value)
            parts.append(
                f'<line x1="{pixel:.3f}" y1="{top:.1f}" x2="{pixel:.3f}" y2="{top + plot_height:.1f}" stroke="#d97706" stroke-width="1.5" stroke-dasharray="6 5"/>'
            )
    for value in horizontal_lines:
        if y_min <= value <= y_max:
            pixel = map_y(value)
            parts.append(
                f'<line x1="{left:.1f}" y1="{pixel:.3f}" x2="{left + plot_width:.1f}" y2="{pixel:.3f}" stroke="#d97706" stroke-width="1.5" stroke-dasharray="6 5"/>'
            )
    for x_value, y_value, marked in zip(x, y, highlighted, strict=True):
        color = "#c2415d" if bool(marked) else "#5c7184"
        opacity = "0.72" if bool(marked) else "0.36"
        parts.append(
            f'<circle cx="{map_x(float(x_value)):.3f}" cy="{map_y(float(y_value)):.3f}" r="2.1" fill="{color}" fill-opacity="{opacity}"/>'
        )
    parts.extend(
        [
            f'<rect x="{left:.1f}" y="{top:.1f}" width="{plot_width:.1f}" height="{plot_height:.1f}" fill="none" stroke="#263746" stroke-width="1.2"/>',
            f'<text x="{left + plot_width / 2:.1f}" y="{SVG_HEIGHT - 22:.1f}" text-anchor="middle" font-family="sans-serif" font-size="15" fill="#263746">{escape(x_label)}</text>',
            f'<text x="22" y="{top + plot_height / 2:.1f}" text-anchor="middle" transform="rotate(-90 22 {top + plot_height / 2:.1f})" font-family="sans-serif" font-size="15" fill="#263746">{escape(y_label)}</text>',
            "</svg>",
            "",
        ]
    )
    return "\n".join(parts).encode("utf-8")


def _member(record: Mapping[str, Any], *, role: str, media_type: str) -> dict[str, Any]:
    return {
        "path": record["path"],
        "sha256": record["sha256"],
        "size_bytes": record["size_bytes"],
        "role": role,
        "media_type": media_type,
    }


def _canonical_digest(value: Any) -> str:
    payload = json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def build_display_bundle(
    *,
    display_dir: Path,
    logcpm_path: Path,
    core_dir: Path,
    core_artifacts: Sequence[Mapping[str, Any]],
    backend_document: Mapping[str, Any],
    backend_data: Mapping[str, Any],
    configuration: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Build and internally verify one display sidecar before publication."""

    try:
        import numpy as np
        from .qc_math import (
            PCA_METHOD_ID,
            centered_unscaled_pca,
            select_top_variable_features,
        )
    except ImportError as error:
        _qc_error(
            "Display generation requires the locked NumPy runtime.",
            reason="numpy_unavailable",
            cause_type=type(error).__name__,
        )

    if display_dir.exists():
        raise BackendFailedError("The private display stage already exists.")
    try:
        display_dir.mkdir(mode=0o700)
    except OSError as error:
        raise OutputWriteError(
            "The private display stage could not be created.",
            path=display_dir,
            operation="create_display_stage",
            cause=error,
        ) from error

    logcpm, logcpm_record = _parse_logcpm(logcpm_path)
    contrasts, results_source, tested_genes = _parse_results(
        core_dir / "results.tsv", backend_document["contrasts"]
    )
    if set(logcpm.gene_ids) != tested_genes:
        raise BackendFailedError(
            "The R-exported logCPM genes do not equal the tested results universe."
        )
    if (
        backend_data.get("tested_gene_count") != len(logcpm.gene_ids)
        or backend_data.get("gene_count") is None
    ):
        raise BackendFailedError("The backend response and display matrix disagree.")
    display_export = backend_data.get("display_export")
    if not isinstance(display_export, Mapping) or (
        display_export.get("gene_count") != len(logcpm.gene_ids)
        or display_export.get("sample_count") != len(logcpm.sample_ids)
        or display_export.get("purpose") != "display_only_not_for_inference"
    ):
        raise BackendFailedError("The backend display-export evidence is incompatible.")

    matrix_by_sample = logcpm.values.T
    selection = select_top_variable_features(
        matrix_by_sample,
        logcpm.gene_ids,
        int(configuration["pca_top_n"]),
    )
    selected_matrix = matrix_by_sample[:, selection.indices]
    displayed_components = tuple(
        int(value) for value in configuration["pca_components"]
    )
    component_count = max(displayed_components)
    pca = centered_unscaled_pca(selected_matrix, n_components=component_count)
    x_component, y_component = (value - 1 for value in displayed_components)

    logcpm_payload, repeated_digest, repeated_size = _capture(logcpm_path)
    if (
        repeated_digest != logcpm_record["sha256"]
        or repeated_size != logcpm_record["size_bytes"]
    ):
        raise BackendFailedError("The display logCPM source changed before staging.")
    copied_logcpm = _write_bytes(display_dir / "logcpm.tsv", logcpm_payload)
    coordinate_header = [
        "sample_id",
        *(f"PC{index + 1}" for index in range(component_count)),
    ]
    coordinate_rows = [
        [sample, *(format(float(value), ".17g") for value in coordinate)]
        for sample, coordinate in zip(logcpm.sample_ids, pca.coordinates, strict=True)
    ]
    coordinate_record = _write_tsv(
        display_dir / "pca_coordinates.tsv", coordinate_header, coordinate_rows
    )
    selected_record = _write_tsv(
        display_dir / "pca_selected_genes.tsv",
        ["rank", "gene_id", "sample_variance"],
        [
            [rank, gene, format(float(variance), ".17g")]
            for rank, (gene, variance) in enumerate(
                zip(selection.names, selection.variances, strict=True), start=1
            )
        ],
    )
    pca_svg = _scatter_svg(
        title="PCA of post-filter logCPM",
        x_label=(
            f"PC{displayed_components[0]} "
            f"({100 * float(pca.explained_variance_ratio[x_component]):.1f}%)"
        ),
        y_label=(
            f"PC{displayed_components[1]} "
            f"({100 * float(pca.explained_variance_ratio[y_component]):.1f}%)"
        ),
        x=pca.coordinates[:, x_component],
        y=pca.coordinates[:, y_component],
        highlighted=np.zeros(len(logcpm.sample_ids), dtype=bool),
        vertical_lines=(0.0,),
        horizontal_lines=(0.0,),
    )
    pca_record = _write_bytes(display_dir / "pca.svg", pca_svg)

    members = [
        _member(
            copied_logcpm, role="display_logcpm", media_type="text/tab-separated-values"
        ),
        _member(
            coordinate_record,
            role="pca_coordinates",
            media_type="text/tab-separated-values",
        ),
        _member(
            selected_record,
            role="pca_selected_genes",
            media_type="text/tab-separated-values",
        ),
        _member(pca_record, role="pca_plot", media_type="image/svg+xml"),
    ]
    plots: list[dict[str, Any]] = [
        {
            "path_id": "edger_ql_p0_v1.display.pca",
            "path": pca_record["path"],
            "sha256": pca_record["sha256"],
            "size_bytes": pca_record["size_bytes"],
            "source_artifacts": [dict(logcpm_record)],
            "gene_count": len(selection.names),
            "excluded_gene_count": int(backend_data["gene_count"])
            - len(logcpm.gene_ids),
            "sample_count": len(logcpm.sample_ids),
            "threshold_lines": [
                {"axis": "x", "kind": "principal_component_origin", "value": 0.0},
                {"axis": "y", "kind": "principal_component_origin", "value": 0.0},
            ],
            "explained_variance_ratio": [
                float(value) for value in pca.explained_variance_ratio
            ],
        }
    ]
    fdr_threshold = float(configuration["fdr_threshold"])
    for points in contrasts:
        highlighted = points.fdr <= fdr_threshold
        volcano_y = -np.log10(np.maximum(points.fdr, FDR_FLOOR))
        volcano_name = f"volcano--{points.contrast_id}.svg"
        volcano_record = _write_bytes(
            display_dir / volcano_name,
            _scatter_svg(
                title=f"Volcano: {points.contrast_id}",
                x_label="log2 fold change",
                y_label="-log10 FDR",
                x=points.logfc,
                y=volcano_y,
                highlighted=highlighted,
                vertical_lines=(-points.lfc_threshold, points.lfc_threshold),
                horizontal_lines=(-math.log10(max(fdr_threshold, FDR_FLOOR)),),
            ),
        )
        ma_name = f"ma--{points.contrast_id}.svg"
        ma_record = _write_bytes(
            display_dir / ma_name,
            _scatter_svg(
                title=f"MA: {points.contrast_id}",
                x_label="average logCPM",
                y_label="log2 fold change",
                x=points.logcpm,
                y=points.logfc,
                highlighted=highlighted,
                horizontal_lines=(-points.lfc_threshold, points.lfc_threshold),
            ),
        )
        members.extend(
            [
                _member(
                    volcano_record, role="volcano_plot", media_type="image/svg+xml"
                ),
                _member(ma_record, role="ma_plot", media_type="image/svg+xml"),
            ]
        )
        common = {
            "source_artifacts": [dict(results_source)],
            "gene_count": len(points.logfc),
            "excluded_gene_count": points.excluded_gene_count,
            "sample_count": len(logcpm.sample_ids),
            "highlight_thresholds": {"FDR": fdr_threshold},
        }
        plots.extend(
            [
                {
                    "path_id": f"edger_ql_p0_v1.display.volcano.{points.contrast_id}",
                    "path": volcano_record["path"],
                    "sha256": volcano_record["sha256"],
                    "size_bytes": volcano_record["size_bytes"],
                    **common,
                    "threshold_lines": [
                        {"axis": "y", "kind": "FDR", "value": fdr_threshold},
                        {
                            "axis": "x",
                            "kind": "absolute_log2_fold_change_test_threshold",
                            "values": [-points.lfc_threshold, points.lfc_threshold],
                        },
                    ],
                },
                {
                    "path_id": f"edger_ql_p0_v1.display.ma.{points.contrast_id}",
                    "path": ma_record["path"],
                    "sha256": ma_record["sha256"],
                    "size_bytes": ma_record["size_bytes"],
                    **common,
                    "threshold_lines": [
                        {
                            "axis": "y",
                            "kind": "absolute_log2_fold_change_test_threshold",
                            "values": [-points.lfc_threshold, points.lfc_threshold],
                        }
                    ],
                },
            ]
        )

    artifact_by_path = {
        str(item["relative_path"]): {
            "path": str(item["relative_path"]),
            "sha256": str(item["sha256"]),
            "size_bytes": int(item["size_bytes"]),
        }
        for item in core_artifacts
    }
    if set(artifact_by_path) != set(_CORE_FILES):
        raise BackendFailedError("The display source bundle is incomplete.")
    source_members = [artifact_by_path[name] for name in _CORE_FILES]
    source_bundle = {
        "source_bundle_id": _canonical_digest(source_members),
        "backend": "edgeR_QL",
        "execution_scope": str(backend_document["execution_scope"]),
        "plan_id": str(backend_document["input_evidence"].get("plan_id", "")),
        "analysis_request_sha256": str(backend_document["analysis_request"]["sha256"]),
        "analysis_request_schema_version": "1.1",
        "members": source_members,
    }
    ql_parameters = backend_data.get("ql_fit_parameters")
    if not isinstance(ql_parameters, Mapping):
        raise BackendFailedError("The backend response lacks QL parameter provenance.")
    manifest = {
        "schema_version": DISPLAY_SCHEMA_VERSION,
        "kind": "rnaseq_downstream_display_manifest",
        "renderer": {
            "renderer_id": DISPLAY_RENDERER_ID,
            "format": "static_svg",
            "width": SVG_WIDTH,
            "height": SVG_HEIGHT,
        },
        "source_bundle": source_bundle,
        "configuration": dict(configuration),
        "methods": {
            "ql_fit_parameters": dict(ql_parameters),
            "logcpm": dict(display_export),
            "pca": {
                "method_id": PCA_METHOD_ID,
                "selection": "top_positive_sample_variance_then_gene_id",
                "center": True,
                "scale": False,
                "requested_top_n": int(configuration["pca_top_n"]),
                "selected_gene_count": len(selection.names),
                "requested_components": list(displayed_components),
                "computed_component_count": component_count,
            },
            "volcano": {
                "x": "results.tsv:logFC",
                "y": "-log10(results.tsv:FDR)",
                "fdr_floor_for_rendering": FDR_FLOOR,
                "statistical_recalculation": False,
            },
            "ma": {
                "x": "results.tsv:logCPM",
                "y": "results.tsv:logFC",
                "statistical_recalculation": False,
            },
        },
        "plots": plots,
        "members": sorted(members, key=lambda item: item["path"]),
    }
    manifest_payload = (
        json.dumps(
            manifest,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    manifest_record = _write_bytes(
        display_dir / "display_manifest.json", manifest_payload
    )
    try:
        descriptor = os.open(display_dir, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except OSError as error:
        raise OutputWriteError(
            "The display bundle could not be synchronized.",
            path=display_dir,
            operation="fsync_display_bundle",
            cause=error,
        ) from error
    artifacts = [
        {
            "kind": "generated_display_artifact",
            "role": item["role"],
            "relative_path": f"display/{item['path']}",
            "sha256": item["sha256"],
            "size_bytes": item["size_bytes"],
            "media_type": item["media_type"],
            "schema_version": DISPLAY_SCHEMA_VERSION,
        }
        for item in manifest["members"]
    ]
    artifacts.append(
        {
            "kind": "generated_display_artifact",
            "role": "display_manifest",
            "relative_path": "display/display_manifest.json",
            "sha256": manifest_record["sha256"],
            "size_bytes": manifest_record["size_bytes"],
            "media_type": "application/json",
            "schema_version": DISPLAY_SCHEMA_VERSION,
        }
    )
    return sorted(artifacts, key=lambda item: item["relative_path"])


def _strict_json(content: bytes) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError(f"duplicate key: {key}")
            value[key] = item
        return value

    try:
        document = json.loads(
            content.decode("utf-8"),
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-standard constant: {value}")
            ),
            object_pairs_hook=no_duplicates,
        )
    except (UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise InputIntegrityError(
            "The display manifest is not strict UTF-8 JSON.",
            details={"cause_type": type(error).__name__},
            cause=error,
        ) from error
    if not isinstance(document, dict):
        raise InputIntegrityError("The display manifest root must be an object.")
    return document


def _exact_keys(value: Mapping[str, Any], expected: set[str], *, context: str) -> None:
    if set(value) != expected:
        raise InputIntegrityError(
            f"The {context} schema is incompatible.",
            details={
                "missing_fields": sorted(expected - set(value)),
                "unexpected_fields": sorted(set(value) - expected),
            },
        )


def _json_exact_equal(observed: Any, expected: Any) -> bool:
    """Compare JSON values without Python's bool/int numeric coercion."""

    if isinstance(expected, Mapping):
        return (
            isinstance(observed, Mapping)
            and set(observed) == set(expected)
            and all(_json_exact_equal(observed[key], expected[key]) for key in expected)
        )
    if isinstance(expected, list):
        return (
            isinstance(observed, list)
            and len(observed) == len(expected)
            and all(
                _json_exact_equal(observed_item, expected_item)
                for observed_item, expected_item in zip(observed, expected, strict=True)
            )
        )
    if expected is None:
        return observed is None
    if isinstance(expected, bool):
        return type(observed) is bool and observed is expected
    if isinstance(expected, int):
        return type(observed) is int and observed == expected
    if isinstance(expected, float):
        return type(observed) is float and observed == expected
    if isinstance(expected, str):
        return type(observed) is str and observed == expected
    return type(observed) is type(expected) and observed == expected


def _expected_display_roles(contrast_ids: Sequence[str]) -> dict[str, tuple[str, str]]:
    expected = {
        "logcpm.tsv": ("display_logcpm", "text/tab-separated-values"),
        "pca.svg": ("pca_plot", "image/svg+xml"),
        "pca_coordinates.tsv": ("pca_coordinates", "text/tab-separated-values"),
        "pca_selected_genes.tsv": (
            "pca_selected_genes",
            "text/tab-separated-values",
        ),
    }
    for identifier in contrast_ids:
        expected[f"volcano--{identifier}.svg"] = ("volcano_plot", "image/svg+xml")
        expected[f"ma--{identifier}.svg"] = ("ma_plot", "image/svg+xml")
    return expected


def _verify_svg(payload: bytes, *, expected: bytes, path: str) -> None:
    if payload != expected:
        raise InputIntegrityError(
            "A display SVG does not reproduce from its bound source values.",
            details={"path": path},
        )
    if b"<script" in payload.lower() or b"<!doctype" in payload.lower():
        raise InputIntegrityError(
            "A display SVG contains active or external document content.",
            details={"path": path},
        )


def verify_display_bundle(
    *,
    display_dir: Path,
    core_dir: Path,
    core_captured: Mapping[str, tuple[bytes, str, int]],
    analysis: Mapping[str, Any],
    backend_manifest: Mapping[str, Any],
    expected_sample_ids: Sequence[str],
) -> dict[str, Any]:
    """Fail closed over a published C1 display extension.

    The verifier reproduces coordinates and SVG bytes from the bound core
    results and R-exported logCPM matrix.  This makes source-value fidelity an
    executable invariant rather than trusting a self-described manifest.
    """

    try:
        from .analysis_contract import _parse_display
        from .qc_math import (
            PCA_METHOD_ID,
            centered_unscaled_pca,
            select_top_variable_features,
        )
        import numpy as np

        contrasts_raw = analysis.get("contrasts")
        if not isinstance(contrasts_raw, list) or not contrasts_raw:
            raise InputIntegrityError("The core contrast inventory is unavailable.")
        contrast_ids = [str(item.get("contrast_id", "")) for item in contrasts_raw]
        if any(not identifier for identifier in contrast_ids) or len(
            set(contrast_ids)
        ) != len(contrast_ids):
            raise InputIntegrityError("The core contrast inventory is invalid.")
        expected_roles = _expected_display_roles(contrast_ids)
        expected_names = set(expected_roles) | {"display_manifest.json"}

        try:
            if display_dir.is_symlink():
                raise OSError("symbolic links are not accepted")
            before = display_dir.stat()
            entries = list(display_dir.iterdir())
        except OSError as error:
            raise InputIntegrityError(
                "The display directory cannot be captured safely.",
                details={"path": str(display_dir)},
                cause=error,
            ) from error
        names = {entry.name for entry in entries}
        unsafe = sorted(
            entry.name for entry in entries if entry.is_symlink() or not entry.is_file()
        )
        if names != expected_names or unsafe:
            raise InputIntegrityError(
                "The display directory inventory is incomplete or unsafe.",
                details={
                    "missing_files": sorted(expected_names - names),
                    "unexpected_files": sorted(names - expected_names),
                    "unsafe_files": unsafe,
                },
            )
        captured = {
            name: _capture(display_dir / name) for name in sorted(expected_names)
        }
        after = display_dir.stat()
        final_names = {entry.name for entry in display_dir.iterdir()}
        if (before.st_dev, before.st_ino) != (
            after.st_dev,
            after.st_ino,
        ) or final_names != expected_names:
            raise InputIntegrityError("The display directory changed during capture.")

        manifest = _strict_json(captured["display_manifest.json"][0])
        _exact_keys(
            manifest,
            {
                "schema_version",
                "kind",
                "renderer",
                "source_bundle",
                "configuration",
                "methods",
                "plots",
                "members",
            },
            context="display manifest",
        )
        if (
            type(manifest["schema_version"]) is not str
            or manifest["schema_version"] != DISPLAY_SCHEMA_VERSION
            or type(manifest["kind"]) is not str
            or manifest["kind"] != "rnaseq_downstream_display_manifest"
            or not _json_exact_equal(
                manifest["renderer"],
                {
                    "renderer_id": DISPLAY_RENDERER_ID,
                    "format": "static_svg",
                    "width": SVG_WIDTH,
                    "height": SVG_HEIGHT,
                },
            )
        ):
            raise InputIntegrityError("The display manifest identity is incompatible.")
        configuration = _parse_display(manifest["configuration"])
        if not _json_exact_equal(manifest["configuration"], configuration):
            raise InputIntegrityError(
                "The normalized display configuration has incompatible JSON types."
            )

        expected_source_members = [
            {
                "path": name,
                "sha256": core_captured[name][1],
                "size_bytes": core_captured[name][2],
            }
            for name in _CORE_FILES
        ]
        expected_source_bundle = {
            "source_bundle_id": _canonical_digest(expected_source_members),
            "backend": "edgeR_QL",
            "execution_scope": "validated_p0_input",
            "plan_id": backend_manifest["input_evidence"]["plan_id"],
            "analysis_request_sha256": analysis["analysis_request"]["sha256"],
            "analysis_request_schema_version": "1.1",
            "members": expected_source_members,
        }
        if not _json_exact_equal(manifest["source_bundle"], expected_source_bundle):
            raise InputIntegrityError(
                "The display manifest does not bind the verified five-file core bundle."
            )

        members = manifest["members"]
        if not isinstance(members, list):
            raise InputIntegrityError("The display member inventory is invalid.")
        by_path: dict[str, Mapping[str, Any]] = {}
        for index, member in enumerate(members):
            if not isinstance(member, Mapping):
                raise InputIntegrityError(
                    "A display member is not an object.",
                    details={"member_index": index},
                )
            _exact_keys(
                member,
                {"path", "sha256", "size_bytes", "role", "media_type"},
                context=f"display member {index}",
            )
            path = member["path"]
            if path not in expected_roles or path in by_path:
                raise InputIntegrityError(
                    "A display member path is invalid.", details={"path": path}
                )
            role, media_type = expected_roles[path]
            _, digest, size = captured[path]
            if not _json_exact_equal(
                member,
                {
                    "path": path,
                    "sha256": digest,
                    "size_bytes": size,
                    "role": role,
                    "media_type": media_type,
                },
            ):
                raise InputIntegrityError(
                    "A display member does not match its manifest.",
                    details={"path": path},
                )
            by_path[path] = member
        if set(by_path) != set(expected_roles):
            raise InputIntegrityError("The display manifest omits required members.")

        logcpm, logcpm_record = _parse_logcpm(display_dir / "logcpm.tsv")
        if not _json_exact_equal(
            logcpm_record,
            {
                "path": "logcpm.tsv",
                "sha256": captured["logcpm.tsv"][1],
                "size_bytes": captured["logcpm.tsv"][2],
            },
        ):
            raise InputIntegrityError(
                "The display logCPM source changed after directory capture."
            )
        core_contrasts = [
            {
                "contrast_id": item["contrast_id"],
                "lfc_threshold": item["lfc_threshold"],
            }
            for item in contrasts_raw
        ]
        contrast_points, results_source, tested_genes = _parse_results(
            core_dir / "results.tsv", core_contrasts
        )
        if not _json_exact_equal(
            results_source,
            {
                "path": "../results.tsv",
                "sha256": core_captured["results.tsv"][1],
                "size_bytes": core_captured["results.tsv"][2],
            },
        ):
            raise InputIntegrityError(
                "The core results source changed after bundle capture."
            )
        if set(logcpm.gene_ids) != tested_genes:
            raise InputIntegrityError(
                "The display logCPM genes differ from the tested core universe."
            )
        sample_count = analysis.get("design", {}).get("sample_count")
        if sample_count != len(logcpm.sample_ids):
            raise InputIntegrityError(
                "The display sample count differs from the verified core design."
            )
        if tuple(expected_sample_ids) != logcpm.sample_ids:
            raise InputIntegrityError(
                "The display sample order differs from the verified core design."
            )
        selection = select_top_variable_features(
            logcpm.values.T, logcpm.gene_ids, configuration["pca_top_n"]
        )
        displayed_components = tuple(configuration["pca_components"])
        component_count = max(displayed_components)
        pca = centered_unscaled_pca(
            logcpm.values.T[:, selection.indices], n_components=component_count
        )
        x_component, y_component = (value - 1 for value in displayed_components)

        coordinate_header = [
            "sample_id",
            *(f"PC{index + 1}" for index in range(component_count)),
        ]
        coordinate_rows = [
            [sample, *(format(float(value), ".17g") for value in coordinate)]
            for sample, coordinate in zip(
                logcpm.sample_ids, pca.coordinates, strict=True
            )
        ]
        buffer = io.StringIO(newline="")
        writer = csv.writer(
            buffer,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            escapechar=None,
        )
        writer.writerow(coordinate_header)
        writer.writerows(coordinate_rows)
        if captured["pca_coordinates.tsv"][0] != buffer.getvalue().encode("utf-8"):
            raise InputIntegrityError(
                "PCA coordinates do not reproduce from the display logCPM matrix."
            )
        selected_buffer = io.StringIO(newline="")
        selected_writer = csv.writer(
            selected_buffer,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            escapechar=None,
        )
        selected_writer.writerow(["rank", "gene_id", "sample_variance"])
        selected_writer.writerows(
            [
                [rank, gene, format(float(variance), ".17g")]
                for rank, (gene, variance) in enumerate(
                    zip(selection.names, selection.variances, strict=True), start=1
                )
            ]
        )
        if captured["pca_selected_genes.tsv"][0] != selected_buffer.getvalue().encode(
            "utf-8"
        ):
            raise InputIntegrityError(
                "The selected PCA genes do not reproduce from logCPM variance."
            )
        expected_pca_svg = _scatter_svg(
            title="PCA of post-filter logCPM",
            x_label=(
                f"PC{displayed_components[0]} "
                f"({100 * float(pca.explained_variance_ratio[x_component]):.1f}%)"
            ),
            y_label=(
                f"PC{displayed_components[1]} "
                f"({100 * float(pca.explained_variance_ratio[y_component]):.1f}%)"
            ),
            x=pca.coordinates[:, x_component],
            y=pca.coordinates[:, y_component],
            highlighted=np.zeros(sample_count, dtype=bool),
            vertical_lines=(0.0,),
            horizontal_lines=(0.0,),
        )
        _verify_svg(captured["pca.svg"][0], expected=expected_pca_svg, path="pca.svg")

        fdr_threshold = configuration["fdr_threshold"]
        for points in contrast_points:
            highlighted = points.fdr <= fdr_threshold
            volcano_y = -np.log10(np.maximum(points.fdr, FDR_FLOOR))
            volcano_path = f"volcano--{points.contrast_id}.svg"
            expected_volcano = _scatter_svg(
                title=f"Volcano: {points.contrast_id}",
                x_label="log2 fold change",
                y_label="-log10 FDR",
                x=points.logfc,
                y=volcano_y,
                highlighted=highlighted,
                vertical_lines=(-points.lfc_threshold, points.lfc_threshold),
                horizontal_lines=(-math.log10(max(fdr_threshold, FDR_FLOOR)),),
            )
            _verify_svg(
                captured[volcano_path][0], expected=expected_volcano, path=volcano_path
            )
            ma_path = f"ma--{points.contrast_id}.svg"
            expected_ma = _scatter_svg(
                title=f"MA: {points.contrast_id}",
                x_label="average logCPM",
                y_label="log2 fold change",
                x=points.logcpm,
                y=points.logfc,
                highlighted=highlighted,
                horizontal_lines=(-points.lfc_threshold, points.lfc_threshold),
            )
            _verify_svg(captured[ma_path][0], expected=expected_ma, path=ma_path)

        methods = manifest["methods"]
        if not isinstance(methods, Mapping):
            raise InputIntegrityError("The display method provenance is invalid.")
        _exact_keys(
            methods,
            {"ql_fit_parameters", "logcpm", "pca", "volcano", "ma"},
            context="display methods",
        )
        ql = methods["ql_fit_parameters"]
        if not isinstance(ql, Mapping) or set(ql) != {
            "abundance.trend",
            "robust",
            "winsor.tail.p",
            "legacy",
            "top.proportion",
            "resolved_top.proportion",
            "keep.unit.mat",
        }:
            raise InputIntegrityError("The locked QL parameter provenance is invalid.")
        resolved_top = ql.get("resolved_top.proportion")
        if (
            isinstance(resolved_top, bool)
            or not isinstance(resolved_top, (int, float))
            or not math.isfinite(float(resolved_top))
            or float(resolved_top) <= 0
            or float(resolved_top) > 1
        ):
            raise InputIntegrityError("The resolved QL top proportion is invalid.")
        expected_ql = {
            "abundance.trend": True,
            "robust": True,
            "winsor.tail.p": [0.05, 0.1],
            "legacy": False,
            "top.proportion": None,
            "resolved_top.proportion": float(resolved_top),
            "keep.unit.mat": False,
        }
        if not _json_exact_equal(ql, expected_ql):
            raise InputIntegrityError("The locked QL parameter provenance is invalid.")
        log_method = methods["logcpm"]
        if not isinstance(log_method, Mapping) or not _json_exact_equal(
            log_method,
            {
                "method": "edgeR::cpm.DGEList",
                "source": "post_filter_post_TMM_observed_DGEList",
                "purpose": "display_only_not_for_inference",
                "arguments": {
                    "normalized.lib.sizes": True,
                    "log": True,
                    "prior.count": 2,
                },
                "scale": "log2",
                "gene_count": len(logcpm.gene_ids),
                "sample_count": len(logcpm.sample_ids),
            },
        ):
            raise InputIntegrityError("The logCPM method provenance is invalid.")
        if not _json_exact_equal(
            methods["pca"],
            {
                "method_id": PCA_METHOD_ID,
                "selection": "top_positive_sample_variance_then_gene_id",
                "center": True,
                "scale": False,
                "requested_top_n": configuration["pca_top_n"],
                "selected_gene_count": len(selection.names),
                "requested_components": list(displayed_components),
                "computed_component_count": component_count,
            },
        ):
            raise InputIntegrityError("The PCA method provenance is invalid.")
        if not _json_exact_equal(
            methods["volcano"],
            {
                "x": "results.tsv:logFC",
                "y": "-log10(results.tsv:FDR)",
                "fdr_floor_for_rendering": FDR_FLOOR,
                "statistical_recalculation": False,
            },
        ) or not _json_exact_equal(
            methods["ma"],
            {
                "x": "results.tsv:logCPM",
                "y": "results.tsv:logFC",
                "statistical_recalculation": False,
            },
        ):
            raise InputIntegrityError(
                "The contrast display method provenance is invalid."
            )

        core_gene_evidence = analysis.get("genes")
        if not isinstance(core_gene_evidence, Mapping):
            raise InputIntegrityError("The core gene-count evidence is invalid.")
        excluded_gene_count = core_gene_evidence.get("filtered")
        if (
            isinstance(excluded_gene_count, bool)
            or not isinstance(excluded_gene_count, int)
            or excluded_gene_count < 0
        ):
            raise InputIntegrityError("The core filtered-gene count is invalid.")
        expected_plots: list[dict[str, Any]] = [
            {
                "path_id": "edger_ql_p0_v1.display.pca",
                "path": "pca.svg",
                "sha256": captured["pca.svg"][1],
                "size_bytes": captured["pca.svg"][2],
                "source_artifacts": [dict(logcpm_record)],
                "gene_count": len(selection.names),
                "excluded_gene_count": excluded_gene_count,
                "sample_count": sample_count,
                "threshold_lines": [
                    {
                        "axis": "x",
                        "kind": "principal_component_origin",
                        "value": 0.0,
                    },
                    {
                        "axis": "y",
                        "kind": "principal_component_origin",
                        "value": 0.0,
                    },
                ],
                "explained_variance_ratio": [
                    float(value) for value in pca.explained_variance_ratio
                ],
            }
        ]
        for points in contrast_points:
            common = {
                "source_artifacts": [dict(results_source)],
                "gene_count": len(points.logfc),
                "excluded_gene_count": points.excluded_gene_count,
                "sample_count": sample_count,
                "highlight_thresholds": {"FDR": fdr_threshold},
            }
            volcano_path = f"volcano--{points.contrast_id}.svg"
            ma_path = f"ma--{points.contrast_id}.svg"
            expected_plots.extend(
                [
                    {
                        "path_id": (
                            f"edger_ql_p0_v1.display.volcano.{points.contrast_id}"
                        ),
                        "path": volcano_path,
                        "sha256": captured[volcano_path][1],
                        "size_bytes": captured[volcano_path][2],
                        **common,
                        "threshold_lines": [
                            {"axis": "y", "kind": "FDR", "value": fdr_threshold},
                            {
                                "axis": "x",
                                "kind": "absolute_log2_fold_change_test_threshold",
                                "values": [
                                    -points.lfc_threshold,
                                    points.lfc_threshold,
                                ],
                            },
                        ],
                    },
                    {
                        "path_id": f"edger_ql_p0_v1.display.ma.{points.contrast_id}",
                        "path": ma_path,
                        "sha256": captured[ma_path][1],
                        "size_bytes": captured[ma_path][2],
                        **common,
                        "threshold_lines": [
                            {
                                "axis": "y",
                                "kind": "absolute_log2_fold_change_test_threshold",
                                "values": [
                                    -points.lfc_threshold,
                                    points.lfc_threshold,
                                ],
                            }
                        ],
                    },
                ]
            )
        if not _json_exact_equal(manifest["plots"], expected_plots):
            raise InputIntegrityError(
                "The display plot metadata does not reproduce from bound sources."
            )

        artifacts = [
            {
                "kind": "consumed_display_artifact",
                "role": by_path[name]["role"],
                "path": str(display_dir / name),
                "sha256": captured[name][1],
                "size_bytes": captured[name][2],
            }
            for name in sorted(by_path)
        ]
        artifacts.append(
            {
                "kind": "consumed_display_artifact",
                "role": "display_manifest",
                "path": str(display_dir / "display_manifest.json"),
                "sha256": captured["display_manifest.json"][1],
                "size_bytes": captured["display_manifest.json"][2],
            }
        )
        return {
            "summary": {
                "schema_version": DISPLAY_SCHEMA_VERSION,
                "status": "verified_complete",
                "source_bundle_id": expected_source_bundle["source_bundle_id"],
                "plot_count": len(expected_plots),
                "plot_types": {
                    "pca": 1,
                    "volcano": len(contrast_ids),
                    "ma": len(contrast_ids),
                },
                "tested_gene_count": len(tested_genes),
                "pca_gene_count": len(selection.names),
                "sample_count": sample_count,
                "configuration": configuration,
            },
            "artifacts": artifacts,
        }
    except InputIntegrityError:
        raise
    except ToolkitError as error:
        raise InputIntegrityError(
            "The display bundle failed independent verification.",
            details={
                "cause_code": error.code.value,
                "cause_details": dict(error.details),
            },
            cause=error,
        ) from error
    except (KeyError, TypeError, ValueError, OSError) as error:
        raise InputIntegrityError(
            "The display bundle is structurally incompatible.",
            details={"cause_type": type(error).__name__},
            cause=error,
        ) from error


__all__ = [
    "DISPLAY_RENDERER_ID",
    "DISPLAY_SCHEMA_VERSION",
    "build_display_bundle",
    "verify_display_bundle",
]
