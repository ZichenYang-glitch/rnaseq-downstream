# Airway edgeR oracle scaffold — not certified

This directory is reserved for the same-engine oracle required by P0:

- the toolkit edgeR v4 quasi-likelihood path and the direct oracle script must
  run in the same locked environment;
- the airway input, direct R script, toolkit request, and expected artifact
  hashes must be frozen;
- coefficients and log fold changes must meet the declared `rtol=1e-6` gate;
- the comparison must write a machine-readable benchmark report.

The current test is explicitly skipped because the edgeR backend and locked
fixture are later P0 work. This scaffold is **not** evidence of parity and
provides **no scientific certification**.
