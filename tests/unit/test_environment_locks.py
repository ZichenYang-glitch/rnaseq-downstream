"""Fail-closed contract tests for the approved two-layer analysis runtime."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import re
import tomllib

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONDA_SPEC = PROJECT_ROOT / "environment.p0.yml"
CONDA_LOCK = PROJECT_ROOT / "conda-lock.yml"
RENV_LOCK = PROJECT_ROOT / "renv.lock"
SOURCE_LOCK = PROJECT_ROOT / "environment" / "r-sources.lock"
BOOTSTRAP_SCRIPT = PROJECT_ROOT / "environment" / "bootstrap.R"
VERIFY_SCRIPT = PROJECT_ROOT / "environment" / "verify.R"
ENVIRONMENT_README = PROJECT_ROOT / "environment" / "README.md"
BASELINE_RENV_LOCK = PROJECT_ROOT / "environment/snapshots/p0-c1-c2-c6b6cd9/renv.lock"
PYPROJECT = PROJECT_ROOT / "pyproject.toml"


EXPECTED_CONDA_PINS = {
    "python": "3.12.12",
    "r-base": "4.6.1",
    "pytest": "8.4.2",
    "numpy": "2.4.4",
    "coreutils": "9.5",
    "zlib": "1.3.2",
    "libuv": "1.52.1",
}
EXPECTED_R_PINS = {
    "renv": "1.2.4",
    "BiocManager": "1.30.27",
    "BiocVersion": "3.23.1",
    "limma": "3.68.0",
    "edgeR": "4.10.0",
    "tximport": "1.40.0",
    "DESeq2": "1.52.0",
    "apeglm": "1.34.0",
    "compcodeR": "1.48.0",
    "airway": "1.32.0",
}
EXPECTED_DESEQ2_CLOSURE = {
    "BH": "1.90.0-1",
    "BiocParallel": "1.46.0",
    "DESeq2": "1.52.0",
    "RcppArmadillo": "15.4.2-1",
    "RcppEigen": "0.3.4.0.2",
    "RcppNumerical": "0.7-0",
    "apeglm": "1.34.0",
    "bbmle": "1.0.26",
    "bdsmatrix": "1.3-7",
    "coda": "0.19-4.1",
    "emdbook": "1.3.14",
    "formatR": "1.14",
    "futile.logger": "1.4.9",
    "futile.options": "1.0.1",
    "lambda.r": "1.2.4",
    "mvtnorm": "1.4-2",
    "numDeriv": "2016.8-1.1",
    "plyr": "1.8.9",
    "snow": "0.4-4",
}
EXPECTED_CONDA_CONTENT_HASH = (
    "93e7d3c507fd94ebf0fb4e0d761b03048064fbb97a9157fd21345781218f4194"
)
EXPECTED_SOURCE_RECORDS = {
    "renv": {
        "version": "1.2.4",
        "repository": "CRAN",
        "role": "bootstrap",
        "url": "https://cran.r-project.org/src/contrib/renv_1.2.4.tar.gz",
        "sha256": "e63c637dc785d55848d9dbc6c9599378103803efd47c1f3f1f82057c00575e8c",
    },
    "BiocManager": {
        "version": "1.30.27",
        "repository": "CRAN",
        "role": "bootstrap",
        "url": "https://cran.r-project.org/src/contrib/BiocManager_1.30.27.tar.gz",
        "sha256": "acbfdcf09602c8279085556ca54531c2ada0ac3c1775d54ce2c3b9a3dc03fbb5",
    },
    "edgeR": {
        "version": "4.10.0",
        "repository": "Bioconductor 3.23",
        "role": "primary",
        "url": "https://bioconductor.org/packages/3.23/bioc/src/contrib/Archive/edgeR/edgeR_4.10.0.tar.gz",
        "sha256": "b02f02639e08ba978ecc0420a770121922f683ff3ebd0645df44c1e67b42a061",
    },
    "limma": {
        "version": "3.68.0",
        "repository": "Bioconductor 3.23",
        "role": "primary",
        "url": "https://bioconductor.org/packages/3.23/bioc/src/contrib/Archive/limma/limma_3.68.0.tar.gz",
        "sha256": "e58fc11741b3c1a3aa311f26101ef6a67404ce92fd4d1cecc89b27aa6f70e095",
    },
    "tximport": {
        "version": "1.40.0",
        "repository": "Bioconductor 3.23",
        "role": "primary",
        "url": "https://bioconductor.org/packages/3.23/bioc/src/contrib/tximport_1.40.0.tar.gz",
        "sha256": "72841a6b2ac1f64ccb05af2d9fba4eea781596771798cefcb3a9327f9dc929f9",
    },
    "DESeq2": {
        "version": "1.52.0",
        "repository": "Bioconductor 3.23",
        "role": "primary",
        "url": "https://bioconductor.org/packages/3.23/bioc/src/contrib/DESeq2_1.52.0.tar.gz",
        "sha256": "8c91699286336350e66eec132ce6fdf5bb4af78e2a4d015a5a61224f62a95984",
    },
    "apeglm": {
        "version": "1.34.0",
        "repository": "Bioconductor 3.23",
        "role": "primary",
        "url": "https://bioconductor.org/packages/3.23/bioc/src/contrib/apeglm_1.34.0.tar.gz",
        "sha256": "78aae4deb507580a491204f5a9462dc154d48e96899c4f1a55a9205aea6d481a",
    },
    "compcodeR": {
        "version": "1.48.0",
        "repository": "Bioconductor 3.23",
        "role": "primary",
        "url": "https://bioconductor.org/packages/3.23/bioc/src/contrib/compcodeR_1.48.0.tar.gz",
        "sha256": "9890c63d8f6cb585ef9311fa888d162ebe1148809c9082d8710f39a07a013b07",
    },
    "airway": {
        "version": "1.32.0",
        "repository": "Bioconductor 3.23",
        "role": "primary",
        "url": "https://bioconductor.org/packages/3.23/data/experiment/src/contrib/airway_1.32.0.tar.gz",
        "sha256": "829f654fc8075ac9d810574e8a8d4672e69e847c45fa2ce598d1ffbc30c7fe19",
    },
}


def _parse_conda_package_records(text: str) -> dict[str, dict[str, str]]:
    """Parse only the simple scalar fields needed from conda-lock YAML."""

    records: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    in_hash = False
    for line in text.splitlines():
        if match := re.fullmatch(r"- name: ([A-Za-z0-9_.+-]+)", line):
            current = {"name": match.group(1)}
            records[current["name"]] = current
            in_hash = False
        elif current is not None and (
            match := re.fullmatch(r"  version: '?([^']+)'?", line)
        ):
            current["version"] = match.group(1)
        elif current is not None and (
            match := re.fullmatch(r"  platform: ([A-Za-z0-9_-]+)", line)
        ):
            current["platform"] = match.group(1)
        elif current is not None and (
            match := re.fullmatch(r"  url: (https://\S+)", line)
        ):
            current["url"] = match.group(1)
        elif current is not None and line == "  hash:":
            in_hash = True
        elif (
            current is not None
            and in_hash
            and (match := re.fullmatch(r"    sha256: ([0-9a-f]{64})", line))
        ):
            current["sha256"] = match.group(1)
        elif current is not None and line and not line.startswith(" "):
            current = None
            in_hash = False
    return records


@pytest.mark.unit
def test_conda_spec_is_exact_and_single_channel() -> None:
    text = CONDA_SPEC.read_text(encoding="utf-8")

    assert re.search(r"(?m)^channels:\n  - conda-forge\ndependencies:$", text)
    for package, version in EXPECTED_CONDA_PINS.items():
        assert f"  - {package}={version}\n" in text
    assert not re.search(r"(?m)^  - (defaults|bioconda|r):", text)


@pytest.mark.unit
def test_conda_lock_matches_source_and_contains_hashed_linux_64_records() -> None:
    text = CONDA_LOCK.read_text(encoding="utf-8")
    records = _parse_conda_package_records(text)
    spec_bytes = CONDA_SPEC.read_bytes()

    assert re.search(r"(?m)^  platforms:\n  - linux-64\n", text)
    assert f"    linux-64: {EXPECTED_CONDA_CONTENT_HASH}\n" in text
    assert "  - environment.p0.yml\n" in text
    assert f"      md5: {hashlib.md5(spec_bytes).hexdigest()}\n" in text
    assert f"      sha256: {hashlib.sha256(spec_bytes).hexdigest()}\n" in text
    assert len(records) >= 50
    assert all(record["version"] for record in records.values())
    assert all(record["platform"] == "linux-64" for record in records.values())
    assert all(
        record["url"].startswith("https://conda.anaconda.org/conda-forge/")
        for record in records.values()
    )
    assert all(
        re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) for record in records.values()
    )
    for package, version in EXPECTED_CONDA_PINS.items():
        assert records[package]["version"] == version
        assert records[package]["platform"] == "linux-64"
        assert records[package]["url"].startswith(
            "https://conda.anaconda.org/conda-forge/"
        )
        assert re.fullmatch(r"[0-9a-f]{64}", records[package]["sha256"])


@pytest.mark.unit
def test_python_build_backend_matches_locked_toolchain() -> None:
    project = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    records = _parse_conda_package_records(CONDA_LOCK.read_text(encoding="utf-8"))

    assert project["build-system"]["requires"] == [
        "setuptools==84.0.0",
        "wheel==0.48.0",
    ]
    assert records["setuptools"]["version"] == "84.0.0"
    assert records["wheel"]["version"] == "0.48.0"
    assert project["project"]["dependencies"] == ["numpy==2.4.4"]


@pytest.mark.unit
def test_source_manifest_is_closed_exact_and_hashed() -> None:
    with SOURCE_LOCK.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert list(rows[0]) == [
        "package",
        "version",
        "repository",
        "role",
        "url",
        "sha256",
    ]
    assert {row["package"] for row in rows} == set(EXPECTED_SOURCE_RECORDS)
    assert len(rows) == len(EXPECTED_SOURCE_RECORDS)
    for row in rows:
        assert {
            "version": row["version"],
            "repository": row["repository"],
            "role": row["role"],
            "url": row["url"],
            "sha256": row["sha256"],
        } == EXPECTED_SOURCE_RECORDS[row["package"]]


@pytest.mark.unit
def test_renv_lock_pins_complete_bioconductor_closure() -> None:
    lock = json.loads(RENV_LOCK.read_text(encoding="utf-8"))

    assert lock["R"]["Version"] == "4.6.1"
    assert lock["Bioconductor"]["Version"] == "3.23"
    assert len(lock["Packages"]) >= 100
    assert all(record.get("Version") for record in lock["Packages"].values())
    for package, version in EXPECTED_R_PINS.items():
        assert lock["Packages"][package]["Version"] == version
    for package, version in EXPECTED_DESEQ2_CLOSURE.items():
        assert lock["Packages"][package]["Version"] == version
    assert len(lock["Packages"]) == 131
    assert "ashr" not in lock["Packages"]
    assert lock["Packages"]["DESeq2"]["Source"] == "Bioconductor"
    assert lock["Packages"]["apeglm"]["Source"] == "Bioconductor"


@pytest.mark.unit
def test_renv_expansion_preserves_every_baseline_record() -> None:
    baseline = json.loads(BASELINE_RENV_LOCK.read_text(encoding="utf-8"))
    current = json.loads(RENV_LOCK.read_text(encoding="utf-8"))

    assert len(baseline["Packages"]) == 112
    assert baseline["R"] == current["R"]
    assert baseline["Bioconductor"] == current["Bioconductor"]
    assert set(current["Packages"]) - set(baseline["Packages"]) == set(
        EXPECTED_DESEQ2_CLOSURE
    )
    assert set(baseline["Packages"]) <= set(current["Packages"])
    for package, record in baseline["Packages"].items():
        assert current["Packages"][package] == record


@pytest.mark.unit
def test_bootstrap_is_noninteractive_checksum_gated_and_ordered() -> None:
    text = BOOTSTRAP_SCRIPT.read_text(encoding="utf-8")

    assert 'TARGET_LIBRARY_ENV <- "RNASEQ_P0_R_LIBRARY"' in text
    assert "assert_vanilla()" in text
    assert 'EXPECTED_SHA256SUM_VERSION <- "9.5"' in text
    assert 'Sys.which("sha256sum")' in text
    assert "verify_archive(archive, record$sha256, sha256sum)" in text
    assert "Sys.setenv(RENV_PATHS_ROOT = state_root)" in text
    assert "file_sha256(lockfile, sha256sum)" in text
    assert "renv::restore(" in text
    assert "renv::refresh()" in text
    assert "exclude = excluded" in text
    assert "clean = TRUE" in text
    assert "strict = TRUE" in text
    assert "transactional = TRUE" in text
    assert "retry = FALSE" in text
    assert "prompt = FALSE" in text
    assert 'file.path(R.home("bin"), "R")' in text
    assert '"CMD",' in text
    assert '"INSTALL",' in text
    assert "if (!identical(status, 0L))" in text
    assert "R CMD INSTALL failed" in text
    assert "install.packages(" not in text
    assert "requireNamespace(package, quietly = TRUE, lib.loc = target)" in text
    assert 'file.path(conda_prefix, "include", "zlib.h")' in text
    assert 'c("--exists", "libuv")' in text
    assert 'TARGET_LIBRARY_MARKER <- ".rnaseq-downstream-p0-library"' in text
    assert "refusing to clean it" in text
    assert (
        '"limma", "edgeR", "tximport", "DESeq2", "apeglm", "compcodeR", "airway"'
        in text
    )
    assert 'DESeq2 = "1.52.0"' in text
    assert 'apeglm = "1.34.0"' in text
    assert "ashr" not in text
    assert "run_verifier(root, target)" in text
    assert "BiocManager::install" not in text


@pytest.mark.unit
def test_verifier_has_primary_pins_and_machine_readable_contract() -> None:
    text = VERIFY_SCRIPT.read_text(encoding="utf-8")

    for package, version in EXPECTED_R_PINS.items():
        assert f'{package} = "{version}"' in text
    assert 'EXPECTED_R_VERSION <- "4.6.1"' in text
    assert "assert_vanilla()" in text
    assert 'EXPECTED_BIOCONDUCTOR_VERSION <- "3.23"' in text
    for package, version in {
        "coreutils": "9.5",
        "zlib": "1.3.2",
        "libuv": "1.52.1",
        "numpy": "2.4.4",
    }.items():
        assert f'{package} = "{version}"' in text
    assert 'Sys.getenv("CONDA_PREFIX"' in text
    assert "load_expected_namespaces(target)" in text
    assert "requireNamespace(package, quietly = TRUE, lib.loc = target)" in text
    assert '"status":"success"' in text
    assert '"status":"error"' in text
    assert 'quit(save = "no", status = status' in text
    assert "jsonlite" not in text


@pytest.mark.unit
def test_environment_documentation_uses_reviewed_generation_command() -> None:
    text = ENVIRONMENT_README.read_text(encoding="utf-8")

    assert "research\n+preview" not in text  # Guard against patch artifacts.
    assert "not\nevidence that the analysis path is certified" in text
    assert "--metadata input_md5" in text
    assert "--metadata input_sha" in text
    assert "--platform linux-64" in text
    assert "conda-lock.yml` header" in text
    assert "RNASEQ_P0_R_LIBRARY" in text
    assert "Rscript --vanilla environment/bootstrap.R" in text
    assert "Rscript --vanilla environment/verify.R" in text
    assert "focused runtime smoke check" in text
    assert "not a complete post-hoc audit" in text
    assert "DESeq2 1.52.0 and apeglm 1.34.0" in text
    assert "preserves every one\nof the previous 112 package records" in text
    assert "`ashr` is intentionally absent" in text
