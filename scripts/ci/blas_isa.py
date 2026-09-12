"""Shared ISA requirements for BLAS diagnostics and certification preflight."""

from __future__ import annotations

import re


CORE_REQUIREMENTS = {
    "automatic": frozenset(),
    "Prescott": frozenset({"sse3"}),
    "Barcelona": frozenset({"sse3"}),
    "Core2": frozenset({"ssse3"}),
    "Nehalem": frozenset({"sse4_1", "sse4_2"}),
    "Sandybridge": frozenset({"avx"}),
    "Haswell": frozenset({"avx2", "fma"}),
    "Zen": frozenset({"avx2", "fma"}),
    "SkylakeX": frozenset(
        {"avx2", "fma", "avx512f", "avx512dq", "avx512bw", "avx512vl", "avx512cd"}
    ),
    "Cooperlake": frozenset(
        {
            "avx2",
            "fma",
            "avx512f",
            "avx512dq",
            "avx512bw",
            "avx512vl",
            "avx512cd",
            "avx512_bf16",
        }
    ),
}


def cpu_flags(output: str) -> set[str]:
    """Parse one nonempty Linux lscpu Flags field before loading any BLAS."""
    fields = []
    for line in output.splitlines():
        key, separator, value = line.partition(":")
        if separator and key.strip() == "Flags":
            fields.append(value.split())
    if not fields:
        raise ValueError(
            "lscpu did not provide CPU flags; explicit core probes are unsafe"
        )
    if len(fields) != 1 or not fields[0]:
        raise ValueError("lscpu must provide exactly one nonempty CPU Flags field")
    if any(re.fullmatch(r"[a-z0-9_]+", flag) is None for flag in fields[0]):
        raise ValueError("lscpu provided malformed CPU flags")
    flags = set(fields[0])
    if "pni" in flags:
        flags.add("sse3")
    return flags
