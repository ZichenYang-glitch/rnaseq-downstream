#!/usr/bin/env python3
"""Fail-closed compatibility stub for an unsafe duplicated legacy workflow.

The former standalone implementation independently rounded Salmon-derived
values, intersected sample sets, aggregated by display symbols, standardized
features before PCA, and continued after contrast failures. Keeping that code
callable would bypass the checkpoint-A integrity corrections in
``legacy/modules/``.
"""

from __future__ import annotations

import sys


def main() -> int:
    """Refuse execution and direct users to the retained experimental path."""

    sys.stderr.write(
        "run_integrated_pydeseq2.py is disabled because it bypasses the "
        "checkpoint-A data-integrity and QC corrections. Use `python -m legacy` "
        "or the legacy Snakemake workflow after archiving prior results. "
        "Neither path is P0-certified.\n"
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
