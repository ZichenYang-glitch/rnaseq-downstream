# CLI JSON contract, version 1.0

The `rnaseq-downstream` command is non-interactive. Standard output contains
exactly one JSON document followed by one newline. Human-readable logs, when
present, are written to standard error. The equivalent module entry point is
`python -m rnaseq_downstream`.

This contract covers the checkpoint-A input gate and the single checkpoint-B
edgeR v4 QL path. The toolkit remains a research preview: the path is narrowly
evidence-gated, but the CLI does not claim that a study or conclusion is
publication-ready.

## Response envelope

Every response includes all seven top-level fields:

```json
{"schema_version":"1.0","command":"capabilities","status":"success","data":{},"warnings":[],"errors":[],"artifacts":[]}
```

- `schema_version`: the response schema version.
- `command`: the requested command, or `cli` when no command was resolved.
- `status`: one of `success`, `error`, or `partial`.
- `data`: command-specific data, or `null` when no data is available.
- `warnings`: structured warning objects.
- `errors`: structured error objects with `code`, `message`, and `details`.
- `artifacts`: structured consumed-input and generated-output records. Every
  record includes a logical `kind`, role, resolved path, SHA-256 digest, and
  byte size; generated JSON records also include media type and schema version.

Envelope construction enforces the status/error invariant: `success` responses
must have an empty `errors` array, while `error` and `partial` responses must
contain at least one structured error. Invalid combinations are rejected before
serialization.

The serializer rejects IEEE-754 non-finite values by using `allow_nan=False`.
Non-ASCII code points are emitted as JSON escapes, so stdout is always a valid
UTF-8 (ASCII-subset) byte stream even when an operating-system argument cannot
be decoded as ordinary Unicode.
If response data cannot be serialized, the CLI emits a minimal valid envelope
with `INTERNAL_ERROR` and exits with code 70.

## Exit codes

| Exit code | Meaning |
| ---: | --- |
| 0 | Success |
| 2 | Invalid/unsupported request or requested input/output boundary failure |
| 3 | Scientific validation failure |
| 4 | Statistical backend failure |
| 5 | Explicitly partial run |
| 70 | Unexpected internal failure |

Exit code 5 must be paired with top-level `status: "partial"`. A backend or
contrast failure must never be represented as a successful response.

## Stable error codes

Contract version 1.0 reserves these codes:

- `INVALID_REQUEST`
- `INPUT_READ_FAILED`
- `OUTPUT_WRITE_FAILED`
- `INPUT_VALIDATION_FAILED`
- `INPUT_EVIDENCE_REQUIRED`
- `INPUT_INTEGRITY_FAILED`
- `COUNT_VALUES_INVALID`
- `SAMPLE_SET_MISMATCH`
- `GENE_ID_INVALID`
- `ASSAY_PROTOCOL_REQUIRED`
- `QC_VALIDATION_FAILED`
- `SALMON_OFFSET_REQUIRED`
- `DESIGN_RANK_DEFICIENT`
- `COVARIATE_CONFOUNDED`
- `CONTRAST_NOT_ESTIMABLE`
- `BACKEND_FAILED`
- `PARTIAL_RUN`
- `FEATURE_NOT_IMPLEMENTED`
- `INTERNAL_ERROR`

Additional codes require a contract revision or a backward-compatible
extension documented here.

## Commands

- `capabilities` returns the implemented surfaces and the exact scope of the
  evidence-gated `edger_ql_p0_v1` path. `certified_analysis_paths` remains
  empty: the airway and compcodeR reports certify `backend_kernel_only`, while
  public input routes are covered separately by the validation contract and
  locked integration tests. Machine fields state that combined-manifest origin
  is self-attested rather than authenticated.

  `analysis_request_schema_versions` lists the accepted public analysis
  request versions. `non_statistical_display_capabilities` describes the
  optional version 1.1 static-SVG sidecar without adding another analysis path.
  Its record fixes the compatible path and request version, plot inventory,
  output location, PCA input/scaling, locked-runtime determinism scope, and
  `summarize` verification mode. `statistical_role:
  "display_only_no_inference"` and `publication_grade_claim: false` are
  normative boundaries, not marketing labels.
- `inspect --request REQUEST.json` resolves, fingerprints, and structurally
  inspects one declared input request without writing files. Because full count
  numeric-domain validation is intentionally not run, its scope says
  `input_semantics: "inspected"` and `full_numeric_validation: "not_run"`.
- `validate --request REQUEST.json --output-dir DIRECTORY` runs the complete
  input-semantic gate and publishes a validation bundle after validation
  completes. Its `input_semantics` scope is `passed` for eligible input or
  `ineligible` for an explicitly recorded, non-certifying high-risk override.
  Its scope still says `design: "not_run"`, `backend: "not_run"`,
  `runnable: false`, and `analysis_path_certified: false`.
- `run --request ANALYSIS_REQUEST.json --output-dir DIRECTORY [--rscript
  EXECUTABLE] [--r-library DIRECTORY]` verifies an eligible validation bundle,
  runs the locked edgeR QL backend in an independent non-interactive R process,
  and atomically publishes five core result files. A version 1.1 request also
  publishes the requested `display/` sidecar in that same atomic operation. It
  never overwrites an existing output directory. Backend diagnostics captured
  by the adapter are forwarded to stderr; the CLI still writes exactly one JSON
  document to stdout.
  The success data includes `ql_fit_parameters` with the explicitly supplied
  `glmQLFit` arguments and resolved `top.proportion`, plus `display_export`
  provenance (or JSON `null` for a version 1.0 run).
- `summarize --run-dir DIRECTORY` independently verifies a public result
  bundle and returns per-contrast status and FDR-at-0.05 counts. When a display
  sidecar is present, it also verifies its inventory, source binding, PCA
  coordinates, and reproducible SVG bytes. It rejects incomplete, modified,
  schema-incompatible, wrong-runtime, and private `backend_kernel_only` output.
- `--help` and subcommand help return their text inside a successful JSON
  envelope. Parser failures return `INVALID_REQUEST` with exit code 2; argparse
  never writes usage text directly to standard output.

Input inspection/validation uses the [input request
contract](input-request-v1.md). Analysis execution uses the [analysis request
contract](analysis-request-v1.md).

## End-to-end invocation

```bash
rnaseq-downstream inspect --request input-request.json
rnaseq-downstream validate \
  --request input-request.json \
  --output-dir evidence/sample-set-1
rnaseq-downstream run \
  --request analysis-request.json \
  --output-dir results/sample-set-1
rnaseq-downstream summarize --run-dir results/sample-set-1
```

Every command is non-interactive. An agent should require process exit code 0,
top-level `status: "success"`, and the command-specific success state. For
`summarize`, that state is `data.status: "verified_complete"`. File presence
alone is never evidence of successful analysis.

## Validation evidence publication

`validate` writes these strict JSON files as one bundle:

- `validated_request.json`: normalized input receipt, warnings, and consumed
  artifact inventory;
- `input_plan.json`: immutable input route plus explicit not-run downstream
  states;
- `provenance.json`: toolkit/runtime identity and consumed file digests; and
- `bundle_manifest.json`: relative member names, sizes, and SHA-256 digests for
  the preceding three documents.

Every document carries the same content-derived `plan_id`. Files are staged in
a sibling directory, individually synchronized, and published with an atomic
no-replace directory rename in the locked Linux environment. The staging and
parent directories are synchronized as well. An existing output, a target
symlink, or a path that appears concurrently is never replaced. Validation
failures create no bundle. Relative output paths are interpreted from the
process working directory and the returned path is canonical and absolute.

After the atomic publish, a failure to synchronize the parent directory cannot
be rolled back safely. In that state the command returns success with the
high-severity warning `OUTPUT_DURABILITY_UNCONFIRMED`, bundle
`publication_status: "published_durability_unconfirmed"`, the `plan_id`, and
the manifest path. All four files are already visible and are returned as
artifacts; consumers should verify the manifest and copy the bundle to durable
storage. Ordinary successful publication reports `durability_confirmed`.

Consumers must distinguish input eligibility from complete analysis. A
successful `validate` response still says `runnable: false` because no design
has been supplied at that stage. Only `run` verifies design rank, residual
degrees of freedom, contrast estimability, exact runtime identity, and backend
completion.

## Result evidence and summary

`run` always publishes the five core regular files `analysis.json`,
`backend_manifest.json`, `coefficients.tsv`, `design.tsv`, and `results.tsv`.
The backend manifest binds the other four core members by relative name,
SHA-256 digest, and byte size. A version 1.1 analysis request additionally
publishes `display/`, whose independent manifest binds its configuration, the
five-file source bundle, and every display member. The run response returns all
generated artifacts with canonical paths and digests. If the bundle is visible
but parent-directory synchronization fails, the command succeeds only with the
high-severity warning
`OUTPUT_DURABILITY_UNCONFIRMED` and
`publication_status: "published_durability_unconfirmed"`.

`summarize` captures and hashes the files again. It requires the locked public
execution identity, checks every manifest record, parses all TSVs with exact
headers and row widths, rejects duplicate gene/contrast rows, verifies the
complete gene inventory, and checks tested/non-tested numeric rules. A tested
`glmQLFTest` row must have a reported finite non-negative F statistic. A tested
`glmTreat` row must explicitly leave the statistic empty and mark it
`not_reported`. Probabilities must be finite and lie in `[0, 1]`.
Within each contrast, `summarize` also recomputes the Benjamini-Hochberg
adjustment from all tested `PValue` entries and rejects any mismatched `FDR`.

The summary returns all five core files as `consumed_analysis_artifact` records
and a compact contrast array containing `status_counts` and
`significance_counts` (`fdr_le_0_05`, `fdr_gt_0_05`, and `not_tested`). If a
display sidecar is present, verified members are additionally returned as
`consumed_display_artifact` records with a display summary. A version 1.0
five-file directory remains valid without `display/`.

For compatibility, the unchanged five-file core does not record whether a
display was requested. `summarize` therefore verifies a display when present
but cannot distinguish a genuine version 1.0 core from a version 1.1 bundle
after the whole `display/` directory has been removed. Consumers that must
detect that omission must retain and verify the version 1.1 request and the
successful `run` receipt, whose artifact inventory includes the display.

## Evidence boundary

The validator verifies declared bytes and relationships; it does not establish
producer truth cryptographically. In particular, the combined
`featurecounts_integer` manifest is self-attested. A forged but internally
consistent manifest can mislabel another integer matrix and pass the input
gate. This limitation is normative, covered by a documented-limitation test,
and is not repaired by later design or backend validation.

Consumers must inspect both the process exit code and the response `status`.
They must not infer scientific success from the presence of an output file or
from a successful input-only validation.
