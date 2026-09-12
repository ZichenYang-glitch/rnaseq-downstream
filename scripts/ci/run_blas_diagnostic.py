#!/usr/bin/env python3
"""Sample ISA-safe OpenBLAS cores against the frozen airway artifact inventory."""

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


PROJECT_ROOT = Path(__file__).resolve().parents[2]
R_PROBE = (
    'cat("R BLAS: ", extSoftVersion()[["BLAS"]], "\\n", sep=""); '
    "invisible(crossprod(matrix(seq_len(4096), nrow=64)))"
)


def candidate_environment(
    base: dict[str, str], core: str, report_directory: Path
) -> dict[str, str]:
    """Preserve all thread controls; only automatic selection removes CORETYPE."""
    environment = base.copy()
    if core == "automatic":
        environment.pop("OPENBLAS_CORETYPE", None)
    else:
        environment["OPENBLAS_CORETYPE"] = core
    environment.pop("OPENBLAS_VERBOSE", None)
    environment["RNASEQ_P0_REQUIRE_BENCHMARKS"] = "1"
    environment["RNASEQ_P0_BENCHMARK_REPORT_DIR"] = str(report_directory)
    return environment


def announce(message: str, environment: dict[str, str]) -> None:
    print(message, flush=True)
    summary = environment.get("GITHUB_STEP_SUMMARY")
    if summary:
        with Path(summary).open("a", encoding="utf-8") as handle:
            handle.write(message + "\n")


def run_command(command: list[str], environment: dict[str, str]) -> int:
    """Inherit stdout/stderr so runtime and gate failures are visible in CI."""
    print("Running: " + " ".join(command), flush=True)
    return subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        env=environment,
        check=False,
        stdin=subprocess.DEVNULL,
        timeout=1800,
    ).returncode


def run_r_probe(prefix: Path, environment: dict[str, str]) -> int:
    """Forward the time-bounded R probe to both CI logs and the step summary."""
    command = [str(prefix / "bin/Rscript"), "--vanilla", "-e", R_PROBE]
    announce("\nR OpenBLAS verbose probe:\n\n```text", environment)
    process = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        env=environment,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=30,
    )
    for line in process.stdout.splitlines():
        announce(line, environment)
    announce("```\n", environment)
    return process.returncode


def run_diagnostic(
    prefix: Path, output: Path, environment: dict[str, str], flags: set[str]
) -> int:
    python = str(prefix / "bin/python")
    runtime_probe = str(PROJECT_ROOT / "scripts/ci/report_blas_runtime.py")
    comparator = str(PROJECT_ROOT / "scripts/ci/check_benchmark_artifacts.py")
    baseline = str(PROJECT_ROOT / "tests/oracle/airway-benchmark-report.json")
    outcomes: dict[str, str] = {}
    failed = False
    for core, requirements in CORE_REQUIREMENTS.items():
        missing = requirements - flags
        if missing:
            outcomes[core] = "skipped: missing ISA flags " + ", ".join(sorted(missing))
            announce(f"{core}: {outcomes[core]}", environment)
            continue
        report_directory = output / core.lower()
        candidate = candidate_environment(environment, core, report_directory)
        probe_environment = {**candidate, "OPENBLAS_VERBOSE": "2"}
        announce(f"\n### BLAS diagnostic candidate: {core}\n", candidate)
        try:
            probe_command = [python, runtime_probe, "--prefix", str(prefix)]
            if core != "automatic":
                probe_command.extend(["--expect-core", core])
            if run_command(probe_command, probe_environment) or run_r_probe(
                prefix, probe_environment
            ):
                outcomes[core] = "probe_failed"
                failed = True
                continue
            gate_command = [
                python,
                "-m",
                "pytest",
                "tests/oracle/test_airway_oracle.py",
                "-v",
            ]
            if run_command(gate_command, candidate):
                outcomes[core] = "gate_failed"
                failed = True
                continue
            compare_command = [
                python,
                comparator,
                "--baseline",
                baseline,
                "--current",
                str(report_directory / "airway-benchmark-report.json"),
            ]
            result = run_command(compare_command, candidate)
            if result == 0:
                outcomes[core] = "match"
            elif result == 1:
                outcomes[core] = "nonmatch"
            else:
                outcomes[core] = "comparison_failed"
                failed = True
        except (OSError, subprocess.SubprocessError) as exc:
            outcomes[core] = f"execution_failed: {exc}"
            failed = True
        finally:
            announce(f"{core}: {outcomes.get(core, 'execution_failed')}", candidate)
    announce("\n### BLAS frozen-artifact sampling outcomes\n", environment)
    for core, outcome in outcomes.items():
        announce(f"- {core}: {outcome}", environment)
    matched = any(outcome == "match" for outcome in outcomes.values())
    if not matched:
        announce("No sampled core reproduced the frozen airway bytes.", environment)
    if failed:
        announce(
            "At least one probe, gate, or comparison execution failed.", environment
        )
    return int(failed or not matched)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    prefix, output = args.prefix.absolute(), args.output.absolute()
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    try:
        for executable in (prefix / "bin/python", prefix / "bin/Rscript"):
            if not executable.is_file():
                raise ValueError(
                    f"Locked runtime executable is unavailable: {executable}"
                )
        library = environment.get("RNASEQ_P0_R_LIBRARY")
        if not library or not Path(library).is_dir():
            raise ValueError("RNASEQ_P0_R_LIBRARY must identify the restored R library")
        result = subprocess.run(
            ["lscpu"],
            env=environment,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        flags = cpu_flags(result.stdout)
        output.mkdir(parents=True, exist_ok=False)
        return run_diagnostic(prefix, output, environment, flags)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
