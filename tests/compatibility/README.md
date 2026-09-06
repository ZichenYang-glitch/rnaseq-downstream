# Environment compatibility evidence

This directory holds generated, machine-readable evidence that an approved
environment expansion preserves earlier certified paths. It is not a location
for hand-authored benchmark claims.

For the DESeq2 1.52.0/apeglm 1.34.0 expansion, first run all four pre-existing
P0/C2 oracle and simulation gates in the freshly rebuilt 131-package R library.
Keep their new reports in a separate directory, then run:

```bash
python scripts/benchmark/build_environment_compatibility_report.py \
  --reports /absolute/path/to/expanded-environment-reports \
  --output tests/compatibility/p1-environment-compatibility-report.json
```

The generator requires every rerun to pass and proves exact SHA-256/size
identity for all 98 numeric artifacts. It permits only the implementation
records for `renv.lock`, `r-sources.lock`, and `verify.R` to change. A relocated
Conda prefix may change only `runtime.rscript_sha256`: all four reports must
carry the same new digest, it must match the live locked `Rscript`, and the old
and new reports must bind the same unchanged `conda-lock.yml`. The transition
is recorded explicitly. The two pathway reports also disclose their rebuilt
installed-package `DESCRIPTION` hashes; these files contain install-time
`Built:` metadata, so the exception is accepted only after proving that the
versioned airway/compcodeR source archive pin is unchanged. All other report
content remains byte-equivalent. Both sides of the environment transition bind
to the immutable P0/C1/C2 snapshot and the current environment files.

The relocated executable and installed `DESCRIPTION` hashes make live rerun-
report and compatibility-report hashes build-specific. CI therefore requires
the generator's structural checks to pass and archives its live output; it does
not byte-compare that output with the committed historical report.

The expected schema is
`scripts/benchmark/environment-compatibility-report-v1.schema.json`. The
committed JSON was generated only after the four real reports passed under the
locked Python 3.12.12/R 4.6.1 runtime rebuilt from both locks in a new Conda
prefix and a new empty R library. Never replace it with a pending, synthetic,
or hand-authored result.
