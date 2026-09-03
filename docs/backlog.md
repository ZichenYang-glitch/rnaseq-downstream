# Research-preview backlog

The following review findings remain outside the current evidence-gated path:

- record the joint adjusted-PCA design-matrix condition number alongside its
  rank diagnostics, without turning a heuristic near-confounding threshold
  into an undocumented hard failure;
- evaluate `csv.QUOTE_NONE` for remaining strict TSV readers where quoted fields are not
  part of the declared format.

These items are not part of the P0 statistical certification claim and require
their own review before implementation.
