# CLI JSON contract, version 1.0

The `rnaseq-downstream` command is non-interactive. Standard output contains
exactly one JSON document followed by one newline. Human-readable logs, when
present, are written to standard error. The equivalent module entry point is
`python -m rnaseq_downstream`.

This contract describes the P0 foundation checkpoint. It does not certify an
analysis path.

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
- `artifacts`: structured artifact records. Paths and digests will be defined by
  the certified execution contract.

Envelope construction enforces the status/error invariant: `success` responses
must have an empty `errors` array, while `error` and `partial` responses must
contain at least one structured error. Invalid combinations are rejected before
serialization.

The serializer rejects IEEE-754 non-finite values by using `allow_nan=False`.
If response data cannot be serialized, the CLI emits a minimal valid envelope
with `INTERNAL_ERROR` and exits with code 70.

## Exit codes

| Exit code | Meaning |
| ---: | --- |
| 0 | Success |
| 2 | Invalid or currently unsupported request |
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
- `SALMON_OFFSET_REQUIRED`
- `DESIGN_RANK_DEFICIENT`
- `CONTRAST_NOT_ESTIMABLE`
- `BACKEND_FAILED`
- `PARTIAL_RUN`
- `FEATURE_NOT_IMPLEMENTED`
- `INTERNAL_ERROR`

Additional codes require a contract revision or a backward-compatible
extension documented here.

## Commands at this checkpoint

- `capabilities` returns the implemented and reserved command surfaces.
- `inspect`, `validate`, `run`, and `summarize` reserve the intended P0 command
  names but return `FEATURE_NOT_IMPLEMENTED` with exit code 2 until their
  respective implementation work is complete.
- `--help` and subcommand help return their text inside a successful JSON
  envelope. Parser failures return `INVALID_REQUEST` with exit code 2; argparse
  never writes usage text directly to standard output.

Consumers must inspect both the process exit code and the response `status`.
They must not infer scientific success from the presence of an output file.
