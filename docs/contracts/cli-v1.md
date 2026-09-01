# CLI JSON contract, version 1.0

The `rnaseq-downstream` command is non-interactive. Standard output contains
exactly one JSON document followed by one newline. Human-readable logs, when
present, are written to standard error. The equivalent module entry point is
`python -m rnaseq_downstream`.

This contract describes the foundation plus checkpoint-A input gate. It does
not certify an analysis path.

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

## Commands at this checkpoint

- `capabilities` returns the implemented and reserved command surfaces. Its
  `certified_analysis_paths` array remains empty.
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
- `run` and `summarize` remain reserved. They return
  `FEATURE_NOT_IMPLEMENTED` with exit code 2.
- `--help` and subcommand help return their text inside a successful JSON
  envelope. Parser failures return `INVALID_REQUEST` with exit code 2; argparse
  never writes usage text directly to standard output.

The request schema is defined separately in the
[input request contract](input-request-v1.md).

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

Consumers must distinguish input eligibility from complete analysis
certification. Even a successful `validate` response is not runnable until the
design and backend gates are implemented and pass.

Consumers must inspect both the process exit code and the response `status`.
They must not infer scientific success from the presence of an output file.
