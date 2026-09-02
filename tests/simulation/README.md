# compcodeR negative-binomial FDR/TPR gate

This directory contains the P0 backend-kernel simulation gate. It freezes
compcodeR 1.48.0, seeds 5001 through 5010, ten balanced negative-binomial
replicates, and the hard FDR/TPR thresholds documented in
`scripts/benchmark/README.md`.

The test regenerates all counts and truth labels, invokes the private backend
adapter, evaluates calls at BH FDR 0.05, and archives a strict JSON report. It
does not treat the synthetic matrices as featureCounts output and does not
certify the public input route.
