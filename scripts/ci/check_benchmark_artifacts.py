#!/usr/bin/env python3
"""Compare full artifact inventories: exit 0 equal, 1 different, 2 invalid input."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys


def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON object key: {key}")
        result[key] = value
    return result


def read_inventory(path: Path) -> tuple[str, dict[str, tuple[str, int]]]:
    report = json.loads(
        path.read_text(encoding="utf-8"), object_pairs_hook=unique_object
    )
    if not isinstance(report, dict) or report.get("status") != "pass":
        raise ValueError(f"{path}: report must have status pass")
    benchmark_id = report.get("benchmark_id")
    if not isinstance(benchmark_id, str) or not benchmark_id:
        raise ValueError(f"{path}: missing or invalid benchmark_id")
    artifacts = report.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ValueError(f"{path}: artifacts must be a nonempty array")
    inventory = {}
    for item in artifacts:
        if not isinstance(item, dict) or set(item) != {"name", "sha256", "size_bytes"}:
            raise ValueError(f"{path}: invalid artifact record shape")
        name, digest, size = item["name"], item["sha256"], item["size_bytes"]
        if not isinstance(name, str) or not name or name in inventory:
            raise ValueError(f"{path}: invalid or duplicate artifact name: {name!r}")
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(f"{path}: invalid artifact SHA-256: {name}")
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise ValueError(f"{path}: invalid artifact size: {name}")
        inventory[name] = (digest, size)
    return benchmark_id, inventory


def compare_inventories(baseline: Path, current: Path) -> tuple[bool, list[str]]:
    baseline_id, expected = read_inventory(baseline)
    current_id, observed = read_inventory(current)
    if current_id != baseline_id:
        raise ValueError(f"benchmark_id mismatch: {baseline_id!r} != {current_id!r}")
    differences = []
    for name in sorted(expected.keys() - observed.keys()):
        differences.append(f"MISSING {name}: baseline={expected[name]}")
    for name in sorted(observed.keys() - expected.keys()):
        differences.append(f"EXTRA {name}: current={observed[name]}")
    for name in sorted(expected.keys() & observed.keys()):
        if expected[name] != observed[name]:
            differences.append(
                f"CHANGED {name}: baseline={expected[name]} current={observed[name]}"
            )
    passed = not differences
    lines = [
        f"Benchmark: {baseline_id}",
        "Comparison: complete artifact name/SHA-256/size inventory; exact equality",
        f"Artifacts: baseline={len(expected)}, current={len(observed)}",
        *differences,
        f"Artifact byte identity: {'PASS' if passed else 'FAIL'}",
    ]
    return passed, lines


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--current", required=True, type=Path)
    args = parser.parse_args(argv)
    invalid = False
    try:
        passed, lines = compare_inventories(args.baseline, args.current)
    except (OSError, ValueError) as exc:
        passed, lines = False, [f"ERROR: {exc}"]
        invalid = True
    print("\n".join(lines))
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        try:
            with Path(summary_path).open("a", encoding="utf-8") as handle:
                handle.write(
                    "\n### Frozen benchmark artifact comparison\n\n```text\n"
                    + "\n".join(lines)
                    + "\n```\n"
                )
        except OSError as exc:
            print(
                f"ERROR: Could not append GitHub step summary: {exc}", file=sys.stderr
            )
            invalid = True
    return 2 if invalid else (0 if passed else 1)


if __name__ == "__main__":
    raise SystemExit(main())
