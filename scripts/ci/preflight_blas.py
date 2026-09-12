#!/usr/bin/env python3
"""Fail before environment restoration unless SkylakeX ISA flags are exposed."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys

if __package__:
    from .blas_isa import CORE_REQUIREMENTS, cpu_flags
else:
    from blas_isa import CORE_REQUIREMENTS, cpu_flags


CERTIFICATION_CORE = "SkylakeX"


def require_certification_isa(flags: set[str]) -> None:
    """Check the fixed certification core, never an environment-selected core."""
    missing = CORE_REQUIREMENTS[CERTIFICATION_CORE] - flags
    if missing:
        raise ValueError(
            "runner does not expose the required AVX-512 ISA subset for SkylakeX; "
            f"missing CPU flags: {', '.join(sorted(missing))}; "
            "rerun to obtain a compatible machine."
        )


def append_summary(message: str) -> None:
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        with Path(summary_path).open("a", encoding="utf-8") as handle:
            handle.write("\n### Certification BLAS ISA preflight\n\n" + message + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.parse_args(argv)
    try:
        result = subprocess.run(
            ["lscpu"],
            env={**os.environ, "LC_ALL": "C"},
            check=True,
            capture_output=True,
            text=True,
            stdin=subprocess.DEVNULL,
            timeout=30,
        )
        flags = cpu_flags(result.stdout)
        require_certification_isa(flags)
        message = (
            "PASS: runner exposes the required SkylakeX ISA flags: "
            + ", ".join(sorted(CORE_REQUIREMENTS[CERTIFICATION_CORE]))
            + "."
        )
        exit_code = 0
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        message = f"ERROR: Certification BLAS ISA preflight failed: {exc}"
        exit_code = 1
    print(message, file=sys.stderr if exit_code else sys.stdout, flush=True)
    try:
        append_summary(message)
    except OSError as exc:
        print(
            f"ERROR: Cannot write the BLAS preflight step summary: {exc}",
            file=sys.stderr,
        )
        return 1
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
