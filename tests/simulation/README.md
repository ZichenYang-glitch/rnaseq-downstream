# compcodeR negative-binomial regression gates

## P0 gene-level gate

The P0 backend-kernel simulation gate freezes compcodeR 1.48.0, seeds 5001
through 5010, ten balanced negative-binomial replicates, and the hard FDR/TPR
thresholds documented in `scripts/benchmark/README.md`.

The test regenerates all counts and truth labels, invokes the private backend
adapter, evaluates calls at BH FDR 0.05, and archives a strict JSON report. It
does not treat the synthetic matrices as featureCounts output and does not
certify the public input route.

## C2 pathway-level gate

The C2 gate evaluates the self-contained mroast and fry directional families.
Its machine-readable result is accompanied by a
[human-readable scope report](pathway-compcoder-benchmark-report.md).
It runs 20 mixed fixtures with 500 true DE genes and 40 complete-null fixtures;
each fixture has 5,000 genes and six samples per condition. Every fixture has
100 null sets across four sizes. Mixed fixtures additionally contain five pure
Up and five pure Down sets of 40 truth-DE genes each.

Membership is generated from truth labels using an RNG stream separate from
count simulation. It never uses observed counts, fitted statistics, P values,
or ranks. Sampling is without replacement within each null set; null sets may
overlap one another. Positive sets are disjoint within direction. The gate also
recalculates BH in Python for every contrast/method/hypothesis family.

The count generator deliberately introduces no set-specific latent factor or
other within-set dependence. Counts are drawn gene by gene conditional on the
frozen gene-level negative-binomial parameters and sample library-size factors.
This gate therefore does not assess fry/mroast calibration when genes within a
set are correlated, including the anti-conservative behavior that correlation
can expose. A correlated-null simulation requires a separate gate.

The frozen limits are mean/worst mixed FDP at most 0.10/0.25, mean/worst TPR at
least 0.80/0.60, no wrong-direction significant positive sets, and at most four
family rejections among 40 complete-null replicates for each method. These
limits were frozen after a disclosed exploratory design pilot. Because the
mroast production seed is reset identically across replicates and the Monte
Carlo study is finite, this is a gross-regression gate, not proof of universal
pathway-test error control or power. Camera is not included in this FDR/TPR gate;
its numerical behavior is covered by the C2 airway oracle.

Using independent replicate-level family indicators as a nominal reference,
a complete-null rejection count `X ~ Binomial(40, 0.05)` fails at `X >= 5`
with probability `0.0480`, which motivates accepting at most 4 of 40. The same
rule has only `P(X >= 5) = 0.3710` power when the true family rejection rate is
`0.10`. The shared mroast rotation seed can weaken replicate independence, so
these probabilities are operating-characteristic references rather than an
exact sampling model for mroast. The gate is a coarse alarm for substantial
miscalibration, not a sensitive test for mild inflation.

As with the P0 simulation, all matrices use the private benchmark-kernel scope.
No featureCounts provenance or public-input-route claim is made.
