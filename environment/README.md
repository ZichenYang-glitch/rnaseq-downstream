# Locked analysis runtime

This directory defines the approved two-layer runtime for the research preview's
evidence-gated analysis paths. It is reproducibility infrastructure, not
evidence that the analysis path is certified. Certification still requires the
locked oracle and simulation gates for each backend.

The lock is currently limited to `linux-64` (x86-64 Linux).

The historical `environment.p0.yml`, `RNASEQ_P0_R_LIBRARY`, and P0 library
marker names remain compatibility identifiers; the locked runtime now serves
both the edgeR and DESeq2 backends.

1. `environment.p0.yml` is the human-reviewed top-level Conda specification.
2. `conda-lock.yml` pins the complete Python, R 4.6.1, compiler, and system
dependency closure, including package URLs and hashes.
   NumPy 2.4.4 is part of this locked layer for the optional C1 PCA display;
   differential-expression inference remains in R.
3. `renv.lock` pins the complete R/Bioconductor 3.23 package closure.
4. `environment/r-sources.lock` pins source URLs and SHA-256 digests for the
   bootstrap packages and the seven primary R packages.

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
`limma`, `edgeR`, `tximport`, `DESeq2`, `apeglm`, `compcodeR`, `airway`.
A source download that fails transiently (for example, a gateway timeout) is
retried with backoff, at most three attempts per archive. Every attempt is
accepted only after its SHA-256 check passes, and exhausting the retries
aborts the bootstrap.

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
`edgeR`, `tximport`, `DESeq2`, `apeglm`, `compcodeR`, and `airway`. Each
expected namespace is loaded from the target library so broken compiled
libraries fail verification.

This command is a focused runtime smoke check, not a complete post-hoc audit of
all Conda and R transitive records. Exact closure installation is enforced by
`conda-lock install` and the strict, clean `renv::restore`; reinstall from both
locks to reassert full closure parity after an environment has been modified.

## BLAS kernel policy and diagnostics

The lock contains `libblas 3.11.0` (`10_h4a7cf45_openblas`) and
`libopenblas 0.3.34` (`pthreads_hcf972fe_1`). The latter is a `DYNAMIC_ARCH`
build: its selected kernel can change with the runner CPU even when all package
hashes are identical. This R build links `libR.so` directly to `libblas.so.3`,
which resolves to `libopenblasp-r0.3.34.so`; there is no separate `libRblas.so`.

The `locked-oracle-and-simulation` certification job pins
`OPENBLAS_CORETYPE=SkylakeX` at job level. Before installing the Conda bootstrap
tool, `scripts/ci/preflight_blas.py` checks the runner's `lscpu` flags for
`avx2`, `fma`, `avx512f`, `avx512dq`, `avx512bw`, `avx512vl`, and `avx512cd`.
These requirements and the flag parser are shared with the diagnostic sampler
through `scripts/ci/blas_isa.py`; there is no separate preflight copy.
Missing or unreadable flags cause a blocking failure before environment
installation, not a warning, skipped certification, or alternate-kernel run.
After restoration, the runtime probe requires `--expect-core SkylakeX`, checking
both the requested environment value and the core actually resolved by BLAS.

Hosted runners are heterogeneous. If preflight reports that the runner does not
expose the required AVX-512 ISA, rerun to obtain a compatible machine. For this
acceptance checkpoint, stop and report after at most five total attempts without
a compatible runner. An ISA pass is only permission to execute the locked gates;
it is not an oracle, simulation, or compatibility pass. A later gate failure
still requires investigation and must not be hidden by runner resampling.

Normal local analysis does not require AVX-512: it uses automatic BLAS dispatch
and retains all ordinary input, design, statistical, and bundle checks. The
SkylakeX pin and preflight are certification-only restrictions for byte-level
replay of the frozen oracle, not new requirements for the toolkit's analysis
CLI. A non-AVX-512 machine can run analyses normally, but cannot reproduce the
frozen oracle bytes with this approved profile. For local evidence replay,
activate the locked environment and pass preflight before applying the pin:

```bash
python scripts/ci/preflight_blas.py &&
  OPENBLAS_CORETYPE=SkylakeX \
  RNASEQ_P0_REQUIRE_BENCHMARKS=1 \
  RNASEQ_P0_RSCRIPT=/absolute/path/to/rnaseq-p0/bin/Rscript \
  RNASEQ_P0_R_LIBRARY=/absolute/path/to/rnaseq-p0-r-library \
  RNASEQ_P0_BENCHMARK_REPORT_DIR=/absolute/path/to/new-benchmark-reports \
  python -m pytest tests -q
```

Replace the paths with the verified locked runtime and a new report directory,
as described in the gate documentation. Keep non-skippable benchmark mode
enabled. Do not use the pin on an incompatible machine.

Certification logs the CPU model, AVX2/AVX-512/FMA flags, loaded BLAS identity,
and unchanged thread settings to the CI log and step summary. These diagnostic
fields are not added to benchmark reports. A manual diagnostic dispatch samples
three runners independently, restores the same locks, prints the automatically
selected core with `OPENBLAS_VERBOSE=2`, and tests eligible explicit cores:

```bash
gh workflow run p0-certification.yml --ref BRANCH -f blas_diagnostic=true
```

This dispatch-only job runs the airway live oracle and separately compares its
complete numeric-artifact SHA-256/size inventory with the frozen report. A
within-run oracle pass is not a frozen-byte pass. Unsupported ISA candidates
are skipped before loading BLAS; forcing a core bypasses OpenBLAS's automatic
CPU checks and must never be used as an ISA compatibility test. Diagnostic
reports are archived separately, without modifying the frozen reports.

A certification attempt requires both frozen-byte reproduction and a supported
ISA on its runner. The approved policy rejects incompatible machines rather
than claiming that every machine in the pool supports the pin. Neither a
tolerance-based compatibility check nor a replacement baseline is permitted.
Do not change `OPENBLAS_NUM_THREADS` while investigating kernel selection.
Keep `OPENBLAS_VERBOSE=2` confined to the standalone diagnostic probes, not
bootstrap, verification, or live gates: extra BLAS messages would contaminate
the verifier's exact NumPy-version output check. The verifier is not relaxed.

The SkylakeX-only certification policy follows the reviewed diagnostic study.
In the local locked-environment sample, `SkylakeX` and `Cooperlake` reproduce all five
frozen airway artifact digests and sizes, while `Prescott`, `Barcelona`, `Core2`,
`Nehalem`, `Sandybridge`, `Haswell`, and `Zen` do not. All ten local runs, including
automatic selection, pass the unchanged within-run oracle. The matching cores
require AVX-512, which is not exposed by every sampled hosted runner. Forcing
either core on an unsupported machine is not a valid remedy: the upstream
[forced-core dispatch](https://github.com/OpenMathLib/OpenBLAS/blob/v0.3.34/driver/others/dynamic.c)
bypasses automatic ISA checks, and
[Cooperlake inherits the SkylakeX kernels](https://github.com/OpenMathLib/OpenBLAS/blob/v0.3.34/kernel/x86_64/KERNEL.COOPERLAKE).

The [dispatch-only diagnostic run](https://github.com/ZichenYang-glitch/rnaseq-downstream/actions/runs/34560991836)
records per-runner hardware and candidate outcomes in logs and step summaries,
with unchanged-schema reports as separate artifacts. The
`blas-kernel-diagnostic` job remains dispatch-only and has no job-level core pin;
it continues to sample automatic dispatch and ISA-eligible explicit candidates
across the runner pool. Do not force an unsupported kernel, change thread
counts, replace frozen reports, or substitute a tolerance check. Frozen-byte
claims remain limited to compatible machines whose pinned runs pass the gates.

## Regenerate the R lock

Author an R-package change only in a fresh library under the locked Conda
runtime. Restore the previous `renv.lock` first, install checksum-verified local
source archives with dependency resolution limited to `Depends`, `Imports`, and
`LinkingTo`, and snapshot the complete isolated library with
`renv::snapshot(type = "all", update = FALSE)`. Primary package records in the
candidate lock must identify `Bioconductor 3.23`, never a local path.

Before replacing the reviewed lock, compare the candidate structurally with the
previous lock. The DESeq2 1.52.0 and apeglm 1.34.0 expansion preserves every one
of the previous 112 package records and adds exactly 19 records, for 131 total.
`ashr` is intentionally absent because the certified shrinkage surface is
apeglm-only. A change to any pre-existing record is a separate environment
upgrade and requires its own review.

After changing either R lock, rebuild into a new absolute Conda prefix and a new
empty R library. Run `bootstrap.R`, then `verify.R`, then every existing oracle
and simulation gate with benchmark skipping disabled. Reusing a populated
library is not accepted as evidence of a clean rebuild.

## Locked source migration policy

A source URL is a retrieval location; the locked SHA-256 digest is the content
identity. If a pinned version leaves the release tree, audit every source URL
with HTTP HEAD. A release URL returning 404 may move to the same package and
version under the official Bioconductor or CRAN `Archive/` tree only after a
fresh download matches the existing SHA-256 exactly. A transport error or 5xx
response is not evidence of a move. A digest mismatch is a hard stop: never
change the version or hash to make a replacement pass, or silently use a mirror.

Record the old URL/status, candidate URL/status, and download digest in a
machine-readable migration audit. Preserve the previous source lock in a new
append-only environment snapshot before editing its URL, and keep every
previously committed snapshot and benchmark report byte-identical. Historical
evidence continues to resolve current file first, then snapshot by hash.
Compatibility checks retain every source identity field and allow only an
official same-package, same-version release-to-Archive locator change; the
resolved source record, including its URL, remains in each compatibility report.
The full fresh rebuild, runtime verification, live gates, and frozen numeric
artifact checks remain mandatory after a migration.

The [2026-09-13 URL audit](source-audits/2026-09-13-url-audit.json) checks all nine
source records. Only compcodeR 1.48.0 moved to its official Archive URL; the
download retains SHA-256
`9890c63d8f6cb585ef9311fa888d162ebe1148809c9082d8710f39a07a013b07`.
The other eight URLs remain unchanged. The prior D1/D2 source-lock bytes are
preserved in the new `snapshots/p1-source-urls-0dba865/` sibling; the existing
P0/C1/C2 snapshot is untouched. The exact-URL test is updated only for this
verified locator migration; package versions and digest assertions do not change.

## Preserve environment evidence across lock changes

An archived benchmark report is immutable: never rewrite its implementation
hashes merely because a current lock or verifier changed. Before an approved
environment evolution, copy every about-to-change environment evidence file
into a new sibling below `environment/snapshots/`, preserving its repository-
relative path. Each snapshot directory has a strict `manifest.json` binding its
payloads to their original paths, byte sizes, SHA-256 digests, and source Git
revision.

Snapshot directories are append-only. Never edit, replace, rename, or delete a
committed snapshot manifest or payload. A later environment evolution creates a
new sibling snapshot. The snapshot is an evidence object, not a standalone
runtime restore; use its `source_revision` to recover the complete historical
tree when a rebuild of that runtime is required.

`scripts/benchmark/evidence_resolver.py` verifies every snapshot before use.
`resolve_archived_implementation_path` first accepts the current repository file
when its path, SHA-256 digest, and size match an archived implementation record.
If the current bytes changed, it resolves only a snapshot payload with the same
original path, digest, and size. `verify_archived_implementation` applies that
rule to the complete implementation inventory. Missing files, malformed
manifests, unlisted files, symlinks, digest drift, size drift, and unresolved
records are hard failures.

The frozen P0/C1/C2 environment evidence is stored in
`environment/snapshots/p0-c1-c2-c6b6cd9/`. Its three payloads correspond to
`renv.lock`, `environment/r-sources.lock`, and `environment/verify.R` at source
revision `c6b6cd970c1a0c6807145dd01033505efae60215`.

After rebuilding an expanded environment, rerun the four existing edgeR and
limma gates into a new, otherwise empty report directory. Do not overwrite the
four frozen reports under `tests/oracle/` and `tests/simulation/`. Once all four
real reruns exist, create the compatibility evidence with:

```bash
python scripts/benchmark/build_environment_compatibility_report.py \
  --reports /absolute/path/to/expanded-environment-reports \
  --output tests/compatibility/p1-environment-compatibility-report.json
```

The generator refuses an existing output path, a missing or failing gate, a
changed frozen report, a stale current environment hash, or numeric artifact
drift. It separately classifies two kinds of non-scientific rebuild metadata.
First, `runtime.rscript_sha256` may change when Conda rewrites an executable for
a different absolute prefix. It must be the only changed runtime field, must
have one identical new value across all four reports, must match the live
`Rscript` beside the locked Python interpreter, and both reports must bind the
same unchanged `conda-lock.yml`. The compatibility report records the old and
new digests and the lock evidence. Second, an installed `airway/DESCRIPTION` or
`compcodeR/DESCRIPTION` record may change through its install-time `Built:`
metadata. That exception is accepted only when the package version and
checksum-pinned source archive are identical across the baseline and expanded
`r-sources.lock`; both DESCRIPTION hashes and the unchanged source identity are
disclosed in the compatibility report.

All other report content must remain byte-equivalent after the three approved
environment records and those explicitly classified metadata fields are
normalized. Because the disclosed runtime and DESCRIPTION hashes also change
the live gate-report hashes, CI validates each new compatibility report
structurally and archives it; it never requires a fresh report to be
byte-identical to the committed historical report. The output shape is defined by
`scripts/benchmark/environment-compatibility-report-v1.schema.json`. A
compatibility report must only be committed after this command consumes the
actual locked-runtime gate outputs; placeholder or hand-authored pass evidence
is forbidden.

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
`r-sources.lock`. New source-manifest entries must use an exact official HTTPS
source-tarball URL and a SHA-256 digest computed from the downloaded bytes with
the locked GNU `sha256sum`. A current Bioconductor tarball may live directly
under `src/contrib`; do not invent a non-existent `/Archive/` URL, and do not
replace an unavailable version with a newer release.
