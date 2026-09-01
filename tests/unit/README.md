# Unit tests

These tests cover the stable error/JSON contracts, strict input semantics,
provenance and validation-bundle behavior, and dependency-light QC math.

The core input and contract tests are standard-library only. Legacy
data-integrity and QC-math modules may import Pandas or NumPy behind explicit
`pytest.importorskip` boundaries, but no unit test imports or executes the
PyDESeq2, plotting, enrichment, or motif stacks.
