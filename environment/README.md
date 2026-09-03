# P0 locked runtime

This directory defines the approved two-layer runtime for the P0 research
preview. It is reproducibility infrastructure, not evidence that the analysis
path is certified. Certification still requires the locked oracle and
simulation gates.

The lock is currently limited to `linux-64` (x86-64 Linux):

1. `environment.p0.yml` is the human-reviewed top-level Conda specification.
2. `conda-lock.yml` pins the complete Python, R 4.6.1, compiler, and system
dependency closure, including package URLs and hashes.
   NumPy 2.4.4 is part of this locked layer for the optional C1 PCA display;
   differential-expression inference remains in R.
3. `renv.lock` pins the complete R/Bioconductor 3.23 package closure.
4. `environment/r-sources.lock` pins source URLs and SHA-256 digests for the
   bootstrap packages and the five primary P0 R packages.

`renv` and `BiocManager` are bootstrapped from their verified source archives
because compatible Conda packages are not available for R 4.6. The primary
packages are also installed from verified archives after `renv` restores their
transitive dependencies. There is no fallback to a latest package version.

## Install

Install the exact Conda layer with `conda-lock` 4.0.1. Use a dedicated absolute
prefix:

```bash
conda-lock install --prefix /absolute/path/to/rnaseq-p0 conda-lock.yml
conda activate /absolute/path/to/rnaseq-p0
python -m pip install --no-deps --no-build-isolation -e .
```

Create an exact R library in another dedicated absolute path:

```bash
Rscript --vanilla environment/bootstrap.R --library /absolute/path/to/rnaseq-p0-r-library
```

The bootstrap requires network access to the exact source URLs in
`r-sources.lock` and to the repositories recorded in `renv.lock`. It downloads
every source-manifest archive, checks it with GNU `sha256sum` 9.5 from the
locked Conda environment, restores the transitive closure, and then installs
primary packages in this order:
`limma`, `edgeR`, `tximport`, `compcodeR`, `airway`.

The bootstrap isolates the renv cache and availability metadata under a
lock-fingerprinted state directory next to the target library. This prevents a
host or older project cache from supplying stale package metadata while keeping
the restored library's cache links valid after the bootstrap process exits.
It refreshes repository availability metadata inside that isolated state, then
performs a strict transactional restore with retries disabled.

The target R library is reconciled with `clean = TRUE`. Do not point the script
at a shared, system, or site library; packages not present in `renv.lock` can be
removed from the dedicated target. The first run requires an empty directory
and writes an ownership marker. Later runs require that marker. The script also
rejects the R system and site library paths.

Instead of `--library`, automation may set an absolute path in
`RNASEQ_P0_R_LIBRARY`:

```bash
RNASEQ_P0_R_LIBRARY=/absolute/path/to/rnaseq-p0-r-library \
  Rscript --vanilla environment/bootstrap.R
```

## Verify

Primary-runtime verification is noninteractive and emits exactly one compact
JSON document on stdout. A mismatch produces `status: "error"` and a nonzero
exit status.

```bash
Rscript --vanilla environment/verify.R --library /absolute/path/to/rnaseq-p0-r-library
```

The verifier checks R 4.6.1, the linux-64 architecture, the active Conda prefix,
the pinned `coreutils`, `zlib`, and `libuv` records, the zlib development
header, pkg-config resolution of libuv, Bioconductor 3.23,
BiocVersion 3.23.1, and exact versions of `renv`, `BiocManager`, `limma`,
`edgeR`, `tximport`, `compcodeR`, and `airway`. Each expected namespace is
loaded from the target library so broken compiled libraries fail verification.

This command is a focused runtime smoke check, not a complete post-hoc audit of
all Conda and R transitive records. Exact closure installation is enforced by
`conda-lock install` and the strict, clean `renv::restore`; reinstall from both
locks to reassert full closure parity after an environment has been modified.

## Regenerate the Conda lock

After an intentional edit to `environment.p0.yml`, generate both input hashes
with this command:

```bash
conda-lock lock \
  --file environment.p0.yml \
  --platform linux-64 \
  --lockfile conda-lock.yml \
  --metadata input_md5 \
  --metadata input_sha
```

Do not copy the regeneration example in the generated `conda-lock.yml` header:
conda-lock 4.0.1 renders the final metadata option and `-f` without a separating
space. The command above is the reviewed invocation.

Any intentional R-package change requires review of both `renv.lock` and
`r-sources.lock`. New source-manifest entries must use an exact HTTPS archive
URL and a SHA-256 digest computed from the downloaded bytes; do not replace an
unavailable archive with a newer release.
