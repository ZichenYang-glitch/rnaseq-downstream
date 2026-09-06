# Research-preview backlog

The following review findings remain outside the current evidence-gated paths:

- record the joint adjusted-PCA design-matrix condition number alongside its
  rank diagnostics, without turning a heuristic near-confounding threshold
  into an undocumented hard failure;
- evaluate `csv.QUOTE_NONE` for remaining strict TSV readers where quoted
  fields are not part of the declared format;
- extend the pathway simulation gate with diluted-positive and mixed-direction
  gene sets, so power and direction checks do not rely only on pure,
  same-direction DE sets;
- add complete-null gene-set simulations with controlled within-set dependence
  (for example, a negative-binomial latent-factor design) to assess fry/mroast
  calibration under correlated genes;
- report benchmark evidence per gene-set method in the `capabilities` JSON,
  explicitly separating the mroast/fry same-engine-plus-simulation evidence
  from camera's same-engine-only evidence;
- clarify the robustness claim for `limma::mroast`, distinguish it from
  upstream `glmQLFit(robust=TRUE)`, and explain the frozen roast options and
  their sensitivity assumptions;
- document finite-rotation P-value granularity for mroast: with 9,999 rotations
  and `midp=false`, the smallest reported mixed or two-sided P value is
  `1e-4`, including the consequence for BH-adjusted results;
- make `summarize` recompute the rank and residual degrees of freedom of a
  DESeq2 `design.tsv`, rather than treating the backend-attested rank as an
  internally consistent provenance field;
- add a controlled R-side fault-injection test for a nonzero or missing apeglm
  convergence code; synthetic result-bundle tests currently lock the
  fail-closed publishing and verification behavior;
- rename the P0-oriented certification workflow after the P1-DESeq2 gates are
  added in D2, so its operational name reflects the expanded environment and
  backend inventory.

These items are outside the current P0/C2 or D1 implementation evidence scope
and require their own review before implementation.
