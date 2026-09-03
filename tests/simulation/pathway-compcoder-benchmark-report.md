# C2 pathway compcodeR gate: scope report

The adjacent
[`pathway-compcoder-benchmark-report.json`](pathway-compcoder-benchmark-report.json)
is the authoritative machine-generated evidence artifact. This companion is a
human-readable interpretation of its frozen design and does not alter the gate
output or add an assertion after the run.

## Recorded outcome

The archived run passed all 17 assertions. For both mroast and fry, the
complete-null family rejected in 3 of 40 replicates, mean and worst mixed-case
TPR were 1.0, and no significant positive set had the wrong direction. The JSON
contains the exact per-replicate metrics, artifacts, runtime, implementation
hashes, thresholds, and assertions.

## Coverage boundary

The compcodeR fixture deliberately introduces no set-specific latent factor or
other within-set correlation structure. Counts are drawn gene by gene
conditional on the frozen gene-level negative-binomial parameters and sample
library-size factors. Consequently, this gate cannot assess fry/mroast
calibration for correlated pathways and cannot expose anti-conservative
behavior caused by within-set dependence.

The positive arm is deliberately strong: each positive set contains 40 pure,
same-direction truth-DE genes. Its current TPR primarily verifies membership,
test routing, and direction handling; it is not a discriminating sensitivity
benchmark for diluted or mixed-direction biological sets.

## Complete-null cutoff

Let `X` be the number of complete-null replicates in which at least one of the
100 BH-adjusted directional tests rejects for one method. Using independent
replicate-level family indicators as a nominal reference model,
`X ~ Binomial(40, 0.05)`. The gate fails at `X >= 5` because
`P(X >= 5 | p = 0.05) = 0.0480`, giving an approximately 5% false-alarm rate
at nominal calibration.

This cutoff is intentionally coarse. If the true family rejection rate is
`0.10`, its probability of failure is only
`P(X >= 5 | p = 0.10) = 0.3710`. The observed `3/40` result has
`P(X >= 3 | p = 0.05) = 0.3233`. These finite-binomial results support a gross
regression check, not a precise claim of pathway-level FDR calibration. The
shared mroast rotation seed can violate replicate independence, so the binomial
calculation is an operating-characteristic reference rather than an exact
sampling model for mroast.

Correlated-null simulations, diluted and mixed-direction positive sets, and
the other planned evidence improvements are tracked in
[`docs/backlog.md`](../../docs/backlog.md).
