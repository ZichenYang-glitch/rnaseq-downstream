# Unit tests

These tests cover the stable error/JSON contracts, strict input semantics,
provenance and validation-bundle behavior, and dependency-light QC math.

The core input and contract tests are standard-library only. NumPy is locked for
the C1 display-only PCA path; legacy data-integrity modules may import Pandas
behind explicit `pytest.importorskip` boundaries. No unit test imports or
executes the PyDESeq2, plotting, enrichment, or motif stacks.
