# Research-preview backlog

The following review findings are intentionally outside checkpoint B. They are
recorded here so that completing the P0 edgeR path does not silently erase
them:

- cap visualization PCA components at `min(n_samples, n_features) - 1`, so a
  second component is not reported when centering leaves only one independent
  direction;
- record the joint adjusted-PCA design-matrix condition number alongside its
  rank diagnostics, without turning a heuristic near-confounding threshold
  into an undocumented hard failure;
- remove the now-empty legacy `Motif Outputs` report card; and
- evaluate `csv.QUOTE_NONE` for strict TSV readers where quoted fields are not
  part of the declared format.

These items are not part of the P0 statistical certification claim and must not
be mixed into the work-item 5--6 correctness diff without a separate review.
