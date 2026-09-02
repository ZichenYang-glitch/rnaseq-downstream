# Input request contract v1

Status: P0 research preview

Schema version: `1.0`

This contract covers input identity, provenance, count semantics, and metadata
alignment only. It does not validate a statistical design, execute an R backend,
or certify an analysis path.

The normative Python entry points are:

```python
from rnaseq_downstream.input_semantics import inspect_request, validate_request
```

Both functions are read-only. They do not create, replace, or normalize input
files and do not write analysis artifacts.

## Return value and certification scope

On success, both functions return a JSON-serializable object with exactly these
top-level members:

```json
{
  "data": {},
  "warnings": [],
  "artifacts": []
}
```

`data.validation_level` is either `inspect` or `validate`.
`data.input_certification_eligible` applies only to this input-semantics gate,
and `data.certification_scope` is always `input_semantics_only`. A clean inspect
sets eligibility to JSON `null`, status to `not_evaluated`, and records
`NUMERIC_DOMAIN_NOT_EVALUATED`; it cannot claim that validation passed. Known
high-risk blockers found during inspection may set eligibility to false and
status to `ineligible`. A completed validation returns a boolean eligibility and
status `passed` or `ineligible`. `input_certification_reasons` makes every
non-passing state explicit. None of these fields means that the design,
contrast, backend, or complete analysis is certified.

`artifacts` is a deterministic inventory of every consumed file. Each entry
contains its logical role, declared path, resolved absolute path, byte size, and
observed SHA-256 digest. Relative paths are resolved against the JSON document
that declares them:

- paths in the request are relative to the request file;
- paths in a featureCounts manifest are relative to the manifest file;
- paths are never interpreted relative to the process working directory.

Files that are only provenance evidence, such as the declared reference FASTA,
are SHA-256 hashed as a stream and are never retained wholesale. Each file that
must be parsed is captured as one immutable in-memory byte snapshot. Both its
reported SHA-256 and its JSON, TSV, or gzip parser consume that exact snapshot.
The snapshot is transferred to the parser and released before the next parsed
file is captured, so memory is not cumulative across samples or reference and
bootstrap artifacts. Parsing never reopens a path and therefore cannot validate
bytes that differ from the reported digest, even if a path is changed and
restored during the operation (an ABA mutation). Peak memory still includes one
largest parsed artifact plus its decoded parser representation; v1 does not yet
provide a streaming tabular parser with equivalent same-byte evidence coupling.
After parsing, every unique consumed path is also fingerprinted again. A
changed, replaced, or truncated current path fails with
`INPUT_INTEGRITY_FAILED`; this final check is additional drift detection, not
the mechanism that couples evidence to parsing.

Expected errors are raised as typed `ToolkitError` subclasses. Inputs are never
silently rounded, filled, intersected, renamed, or dropped.

## `inspect` versus `validate`

`inspect_request` verifies the JSON schema, resolves and hashes all consumed
files, checks typed upstream evidence, reads metadata and identifiers, validates
sample and identifier identity, validates Salmon metadata and inferential
archive structure, and normalizes the future route. It intentionally does not
claim that assay count values have passed their full numeric-domain checks.
Mixed inferential-replicate evidence is
reported as a high-risk warning and makes the input ineligible.

`validate_request` performs every inspection check and additionally:

- validates all `featurecounts_integer` values as literal nonnegative integers;
- validates Salmon numeric fields against an explicit decimal grammar, finite
  downstream IEEE-754 range, and their documented domain;
- fails when inferential replicates are mixed across samples or their method and
  replicate count are inconsistent.

Neither operation validates a design matrix or invokes edgeR.

## Common request object

Unknown keys are rejected. All three semantics require these members:

```json
{
  "schema_version": "1.0",
  "input_semantics": "featurecounts_integer",
  "metadata": {
    "path": "metadata.tsv",
    "sample_id_column": "sample_id"
  },
  "producer": {
    "name": "featureCounts",
    "version": "2.0.6"
  },
  "reference": {
    "name": "GENCODE",
    "version": "M36",
    "source": "reference/transcripts.fa",
    "sha256": "OPTIONAL_EXPECTED_64_HEX_DIGEST"
  },
  "gene_id": {
    "strip_version": false
  },
  "samples": [
    {"sample_id": "sample_a"},
    {"sample_id": "sample_b"}
  ]
}
```

`reference.sha256` is optional, but the source is always hashed and reported.
Producer and reference names and versions must be explicit non-empty strings.
Identity-bearing strings and paths must contain no Unicode control characters;
values are never silently trimmed. JSON integers outside the exactly
interoperable range `[-(2^53 - 1), 2^53 - 1]` are rejected with
`INVALID_REQUEST` instead of leaking implementation-level integer conversion
errors.

Sample identifiers must be unique, non-empty, and contain no whitespace. The
metadata sample set must exactly equal the declared assay sample set. A metadata
row permutation is permitted only because it is explicitly reported with a
`METADATA_ORDER_NORMALIZED` warning; the declared `samples` array remains the
analysis order. Missing or unexpected samples are a hard failure.

The internal key is always `gene_id`. `gene_name` or a symbol is display-only.
Repeated symbols across different stable gene IDs are valid and are never
aggregated. If `strip_version` is true, only a terminal numeric suffix such as
`.12` is removed. Any collision caused by stripping is a hard failure.
Identifier-inventory SHA-256 values use a domain-separated sequence encoding
with an explicit item count and an eight-byte big-endian length prefix before
every UTF-8 value. They are not delimiter-based and therefore cannot be
ambiguous when one identifier contains bytes resembling a separator.

## `featurecounts_integer`

The producer name must be exactly `featureCounts`. Two layouts are accepted.

### Typed combined matrix

The request adds:

```json
{
  "featurecounts": {
    "layout": "combined_matrix",
    "matrix": "counts.tsv",
    "manifest": "counts.manifest.json"
  }
}
```

A combined matrix without a typed manifest is rejected. The manifest has this
shape:

```json
{
  "schema_version": "1.0",
  "artifact_type": "featurecounts_integer_matrix",
  "artifact": {
    "path": "counts.tsv",
    "sha256": "REQUIRED_64_HEX_DIGEST"
  },
  "gene_id_column": "gene_id",
  "display_columns": ["gene_name"],
  "sample_columns": ["sample_a", "sample_b"],
  "producer": {
    "name": "featureCounts",
    "version": "2.0.6"
  },
  "reference": {
    "name": "GENCODE",
    "version": "M36",
    "source": "reference/transcripts.fa",
    "sha256": "REQUIRED_64_HEX_DIGEST"
  }
}
```

The request and manifest must resolve to the same matrix and reference. Their
producer and reference identities must agree. The actual matrix header must be
exactly `gene_id_column`, then `display_columns`, then `sample_columns` in the
declared order. The matrix digest must match the manifest.

#### Trust boundary: manifest self-attestation

The typed combined-matrix manifest is self-attested evidence supplied alongside
the matrix. Validation proves that the request, manifest, matrix bytes, declared
producer/version, reference, schema, and integer count domain are internally
consistent. It does **not** cryptographically prove that featureCounts produced
those bytes. In particular, a user could round a Salmon-derived merged matrix,
declare it as `featurecounts_integer`, and create a matching manifest and
digest; an internally consistent construction of that kind passes this
input-only gate even though its producer claim is false. The validator cannot
infer or authenticate producer origin from integer matrix contents. Users who
need stronger origin evidence must retain original per-sample featureCounts
files and external workflow provenance or use the `per_sample_files` layout.

### Original per-sample files

The request adds `featurecounts.layout = "per_sample_files"`, and each sample
adds both an original file and its exact count column:

```json
{
  "samples": [
    {
      "sample_id": "sample_a",
      "counts_file": "featurecounts/sample_a.txt",
      "count_column": "sample_a.bam"
    }
  ],
  "featurecounts": {"layout": "per_sample_files"}
}
```

Each file must contain one unambiguous `Program:featureCounts v...` header and
the canonical columns `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`, plus
exactly the declared count column. All files must contain the same stable gene
IDs in the same order.

For both layouts, accepted count text is limited to literal decimal digits.
Missing values, negative values, fractions, `NaN`, infinities, and scientific
notation are rejected. Values greater than `2^53 - 1` are also rejected because
they cannot be represented exactly by the downstream R numeric type. The bound
is reported as `route.maximum_exact_integer`. No coercion or rounding occurs.

The normalized future route is integer counts to `edgeR::DGEList`, with no
transcript-length offset.

## `salmon_quant_dirs_full_length`

The producer name must be exactly `Salmon`. Each sample declares its original
quantification directory, and the request declares the exact tx2gene file:

```json
{
  "input_semantics": "salmon_quant_dirs_full_length",
  "samples": [
    {"sample_id": "sample_a", "quant_dir": "salmon/sample_a"}
  ],
  "salmon": {
    "tx2gene": "reference/tx2gene.tsv",
    "tx2gene_sha256": "OPTIONAL_EXPECTED_64_HEX_DIGEST"
  }
}
```

Every quantification directory must contain and expose provenance for:

- `quant.sf`;
- `cmd_info.json`;
- `<auxDir>/meta_info.json`;
- both `<auxDir>/bootstrap/names.tsv.gz` and `bootstraps.gz` when metadata
  declares replicates.

`<auxDir>` is the exact `cmd_info.json` member `auxDir`, defaulting to
`aux_info` only when that member is absent, matching the locked tximport route.
It must be one non-empty relative POSIX path, may not contain empty, `.`, or
`..` segments, and must resolve to a directory inside that sample's resolved
quantification directory. Absolute paths, Windows drive paths, backslashes,
missing directories, and symlink escapes are hard failures. The same resolved
directory is used for both `meta_info.json` and inferential-replicate archives;
the validator never falls back to a separately hard-coded `aux_info` path.

Both Salmon metadata documents must agree with the declared Salmon version and
identify the index. Every sample must provide canonical hexadecimal
`index_seq_hash` evidence. Other known 256-bit and 512-bit Salmon index hashes
are validated and retained when present, and the complete hash dictionaries
must match exactly across samples. The
`cmd_info.index` value is retained as provenance but is commonly a filesystem
path. Differing paths with equal required hashes emit `SALMON_INDEX_PATHS_DIFFER`
instead of failing, so runs from different nodes or mounts remain portable.
The ordinary file SHA-256 of `reference.source` and Salmon's sequence-derived
`index_seq_hash` use different hash definitions. v1 records them independently
and labels this relationship in `reference.salmon_index_linkage_status`; it
does not falsely claim that a raw file digest proves how the Salmon index was
built.

`quant.sf` must have the canonical columns `Name`, `Length`, `EffectiveLength`,
`TPM`, and `NumReads`. Unlike featureCounts, Salmon `NumReads` are estimated
counts and may be fractional. They are accepted without rounding when finite
and nonnegative. Length must be finite and positive. Effective length must be
finite and nonnegative in general, and strictly positive for the full-length
offset route. TPM must be finite and nonnegative.

For inferential replicates, `names.tsv.gz` must contain exactly the `quant.sf`
transcript order. The uncompressed little-endian f64 payload in
`bootstraps.gz` must have exactly `replicate_count × target_count` finite,
nonnegative values. Legacy little-endian i32 payloads accepted by locked
tximport 1.40 are also recognized explicitly. In either encoding, the metadata
target count must equal the quantification row count. The locked tximport import
remains a backend integration check in the next work item; this structural gate
does not replace it.

All samples must contain the same transcripts in the same order. tx2gene must
map that exact transcript set, with every transcript appearing exactly once.
Many transcripts may map to one stable gene ID. An optional `gene_name` column
is display-only.

Salmon stores the replicate total in `meta_info.num_bootstraps` for both
bootstrap and Gibbs sampling and identifies the method with `samp_type`.
`samp_type` is required and must be exactly `none`, `bootstrap`, or `gibbs`.
`none` requires a zero metadata total, while `bootstrap` and `gibbs` require a
positive total. A positive total also requires the corresponding command-level
bootstrap or Gibbs count to be present and equal; an irrelevant positive count
or any contradiction is rejected. Integer evidence is bounded at `2^53 - 1`,
and an explicit zero is not treated as a missing field. Inferential replicates are
classified as `all`, `none`, or `mixed`. When every sample has the same
replicate method and count, the full-length route records
future use of `tximport.infReps` for relative-abundance-adjusted inferential
overdispersion. `edgeR_options.divide` is true only in that case. Mixed or
inconsistent replicate evidence fails `validate`.

The normalized full-length route fixes:

- `tximport` `countsFromAbundance = "no"`;
- `tximport` `dropInfReps = false`;
- `edgeR::DGEListFromTximport`;
- `DGEListFromTximport(..., divide = true)` only when every sample has
  consistent inferential replicates;
- a required transcript-length offset.

## `salmon_quant_dirs_three_prime`

This semantics requires the literal declaration:

```json
{
  "input_semantics": "salmon_quant_dirs_three_prime",
  "assay_protocol": "three_prime"
}
```

It uses the same quantification-directory and tx2gene evidence checks as the
full-length route. Its normalized route is deliberately different:

- use raw `txi$counts`;
- use `edgeR::DGEList`;
- `countsFromAbundance = "no"`;
- no transcript-length offset;
- no gene-length correction.

A high-risk override can be declared only with an explicit reason:

```json
{
  "analysis_options": {
    "three_prime_length_correction_override": {
      "enabled": true,
      "reason": "Sensitivity analysis required by a documented external protocol."
    }
  }
}
```

The override emits
`HIGH_RISK_THREE_PRIME_LENGTH_CORRECTION_OVERRIDE`, records the reason in the
normalized route, sets `input_certification_eligible` to false, and fixes
`certified_path_execution_permitted` to false. The unchanged no-correction
fields describe the certified default, while the override carries the explicit
policy `not_executable_in_certified_path`; a later backend must not silently
ignore or apply it. The override does not turn this input-only operation into a
runnable or certified analysis.

## Rejected legacy input

A bare merged gene-count matrix, including a rounded
`salmon.merged.gene_counts.tsv`, does not satisfy any v1 input semantics. dtype
or filename heuristics are never used to guess its origin.
