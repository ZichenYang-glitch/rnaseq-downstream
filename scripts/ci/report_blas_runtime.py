#!/usr/bin/env python3
"""Log hardware and loaded BLAS identity without changing benchmark reports."""

from __future__ import annotations

import argparse
import ctypes
import os
from pathlib import Path
import subprocess
import sys


THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OPENBLAS_DEFAULT_NUM_THREADS",
    "GOTO_NUM_THREADS",
    "OMP_NUM_THREADS",
)


def cpu_identity(lscpu_output: str) -> list[str]:
    """Retain CPU identity and the ISA flags relevant to BLAS dispatch."""
    fields = {}
    for line in lscpu_output.splitlines():
        key, separator, value = line.partition(":")
        if separator:
            fields[key.strip()] = value.strip()
    flags = set(fields.get("Flags", "").split())
    selected = sorted(
        flag for flag in flags if flag in {"avx2", "fma"} or flag.startswith("avx512")
    )
    return [
        f"CPU Model name: {fields.get('Model name', 'unavailable')}",
        f"CPU BLAS ISA flags: {' '.join(selected) or '(none reported)'}",
    ]


def link_chain(path: Path) -> str:
    """Describe every symlink hop, including missing or cyclic destinations."""
    parts = []
    visited = set()
    while True:
        path = Path(os.path.abspath(path))
        parts.append(str(path))
        if path in visited:
            parts.append("[symlink cycle]")
            break
        visited.add(path)
        if not path.is_symlink():
            parts.append("[file]" if path.is_file() else "[missing]")
            break
        target = Path(os.readlink(path))
        path = target if target.is_absolute() else path.parent / target
    return " -> ".join(parts)


def openblas_identity(library_path: Path) -> tuple[str, str, int]:
    """Query the actual loaded library, not a guessed core from the CPU name."""
    library = ctypes.CDLL(str(library_path))
    config = library.openblas_get_config
    core = library.openblas_get_corename
    threads = library.openblas_get_num_threads
    for function in (config, core, threads):
        function.argtypes = []
    config.restype = core.restype = ctypes.c_char_p
    threads.restype = ctypes.c_int
    raw_config, raw_core = config(), core()
    if not raw_config or not raw_core:
        raise ValueError("OpenBLAS identity API returned an empty value")
    return raw_config.decode("ascii"), raw_core.decode("ascii"), threads()


def validate_identity(config: str, core: str, expected: str | None) -> None:
    if "DYNAMIC_ARCH" not in config.split():
        raise ValueError("Loaded OpenBLAS configuration does not declare DYNAMIC_ARCH")
    if expected is not None:
        requested = os.environ.get("OPENBLAS_CORETYPE", "")
        if requested.casefold() != expected.casefold():
            raise ValueError(
                f"OPENBLAS_CORETYPE={requested!r} does not match expected {expected!r}"
            )
        if core.casefold() != expected.casefold():
            raise ValueError(
                f"OpenBLAS resolved core {core!r} does not match expected {expected!r}"
            )


def command_log(command: list[str]) -> str:
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
        env={**os.environ, "LC_ALL": "C"},
    )
    if result.returncode:
        raise ValueError(
            f"{' '.join(command)} exited {result.returncode}: {result.stderr.strip()}"
        )
    return result.stdout


def append_summary(title: str, lines: list[str]) -> None:
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        with Path(summary_path).open("a", encoding="utf-8") as handle:
            handle.write(f"\n### {title}\n\n```text\n" + "\n".join(lines) + "\n```\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prefix", type=Path, help="Omit for a pre-restore CPU-only log"
    )
    parser.add_argument(
        "--expect-core", help="Require both requested and resolved core"
    )
    args = parser.parse_args(argv)
    if args.expect_core is not None and args.prefix is None:
        parser.error("--expect-core requires --prefix")
    prefix = args.prefix.absolute() if args.prefix is not None else None
    lines = [f"Locked prefix: {prefix or '(not restored; CPU-only probe)'}"]
    failure = False
    try:
        lines.extend(cpu_identity(command_log(["lscpu"])))
        lines.append(
            "OPENBLAS_CORETYPE requested: "
            + os.environ.get("OPENBLAS_CORETYPE", "(automatic)")
        )
        for name in THREAD_VARIABLES:
            lines.append(f"{name}: {os.environ.get(name, '(unset; unchanged)')}")
        if prefix is None:
            print("\n".join(lines))
            append_summary(
                "CPU runtime provenance before environment restoration", lines
            )
            return 0
        paths = [
            prefix / "lib/R/lib/libRblas.so",
            prefix / "lib/libblas.so.3",
            prefix / "lib/libopenblas.so",
            prefix / "lib/libopenblas.so.0",
        ]
        lines.extend("Library chain: " + link_chain(path) for path in paths)
        config, core, threads = openblas_identity(prefix / "lib/libblas.so.3")
        lines.extend(
            [
                f"OpenBLAS configuration: {config}",
                f"OpenBLAS resolved core: {core}",
                f"OpenBLAS effective threads: {threads} (unchanged)",
            ]
        )
        validate_identity(config, core, args.expect_core)
        libraries = [prefix / "lib/R/lib/libR.so"]
        r_library = os.environ.get("RNASEQ_P0_R_LIBRARY")
        if r_library:
            libraries.append(Path(r_library) / "edgeR/libs/edgeR.so")
        for library in libraries:
            if library.is_file():
                lines.append(f"Linked dependencies ({library}):")
                lines.extend(command_log(["ldd", str(library)]).splitlines())
            else:
                lines.append(f"Linked dependencies unavailable: {library}")
    except (OSError, ValueError, AttributeError, subprocess.SubprocessError) as exc:
        lines.append(f"ERROR: {exc}")
        failure = True
    print("\n".join(lines))
    try:
        append_summary("CPU and BLAS runtime provenance", lines)
    except OSError as exc:
        print(f"ERROR: Could not append GitHub step summary: {exc}", file=sys.stderr)
        failure = True
    return int(failure)


if __name__ == "__main__":
    raise SystemExit(main())
