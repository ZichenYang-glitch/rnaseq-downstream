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
It runs 20 mixed fixtures with 500 true DE genes and 40 complete-null fixtures;
each fixture has 5,000 genes and six samples per condition. Every fixture has
100 null sets across four sizes. Mixed fixtures additionally contain five pure
Up and five pure Down sets of 40 truth-DE genes each.

Membership is generated from truth labels using an RNG stream separate from
count simulation. It never uses observed counts, fitted statistics, P values,
or ranks. Sampling is without replacement within each null set; null sets may
overlap one another. Positive sets are disjoint within direction. The gate also
recalculates BH in Python for every contrast/method/hypothesis family.

The frozen limits are mean/worst mixed FDP at most 0.10/0.25, mean/worst TPR at
least 0.80/0.60, no wrong-direction significant positive sets, and at most four
family rejections among 40 complete-null replicates for each method. These
limits were frozen after a disclosed exploratory design pilot. Because the
mroast production seed is reset identically across replicates and the Monte
Carlo study is finite, this is a gross-regression gate, not proof of universal
pathway-test error control or power. Camera is not included in this FDR/TPR gate;
its numerical behavior is covered by the C2 airway oracle.

As with the P0 simulation, all matrices use the private benchmark-kernel scope.
No featureCounts provenance or public-input-route claim is made.
