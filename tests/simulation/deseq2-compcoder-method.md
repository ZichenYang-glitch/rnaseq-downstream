# DESeq2 compcodeR simulation method and exploratory result

## Scope

This D2 study evaluates the locked DESeq2 1.52.0 Wald path on balanced
negative-binomial simulations. It is a finite regression study, not proof of
universal FDR control or biological power.

Unlike the earlier edgeR backend-kernel simulation, every replicate traverses
the public toolkit chain:

```text
inspect -> validate -> run (backend=deseq2) -> summarize
```

The simulated integer matrix is routed through the
`featurecounts_integer/combined_matrix` contract. Its typed manifest is a
self-attested simulation adapter. It exercises request validation, immutable
input bundling, design lint, the independent R process, result publication,
and independent summary verification; it does not prove that featureCounts
produced the synthetic counts and is not evidence of producer authenticity.

## Frozen scenario and independent grids

Both grids use compcodeR 1.48.0 with 5,000 genes, 500 truth-DE genes, six
samples per condition, balanced up/down regulation, effect size 2, no injected
outliers, and `filter.threshold.total=0`. Truth-null genes are defined solely
by compcodeR's pre-analysis `differential.expression` annotation. Membership
does not consult counts, fitted statistics, P values, ranks, or results.

- Exploratory seeds: 61001 through 61020.
- Held-out certification seeds: 62001 through 62020.

The candidate limits were written into the runner before inspecting the
exploratory grid:

- mean replicate FDP at most 0.065;
- worst replicate FDP at most 0.12;
- mean TPR at least 0.45;
- worst replicate TPR at least 0.35.

The FDP limits allow finite-replicate variation around a nominal BH FDR of
0.05. The TPR limits are deliberately conservative regression floors for this
fixed 6-vs-6 scenario. They are not inferred from one observed run and must not
be revised after looking at the held-out grid.

## Disclosed exploratory result

The authoritative machine-readable record is
[`deseq2-compcoder-exploratory-report.json`](deseq2-compcoder-exploratory-report.json).
All 20 public CLI chains, locked-runtime checks, fixture checks, and independent
`summarize` verifications completed.

At BH FDR 0.05, the observed distribution was:

| Metric | Exploratory result |
| --- | ---: |
| Mean replicate FDP | 0.11821 |
| Median replicate FDP | 0.11842 |
| Minimum / maximum replicate FDP | 0.09117 / 0.13731 |
| Discovery-weighted pooled FDP (descriptive only) | 0.11856 |
| Mean replicate TPR | 0.65050 |
| Median replicate TPR | 0.64300 |
| Minimum / maximum replicate TPR | 0.61800 / 0.68200 |
| Total discoveries | 7,380 |
| Total false positives | 875 |

Across the 20 independent simulation seeds, the replicate-FDP standard
deviation was 0.01389 (standard error 0.00311). A descriptive normal-reference
95% interval for the mean is 0.11212 to 0.12429, still well above the candidate
0.065 bound. This interval is a Monte Carlo diagnostic, not a confidence
interval for performance on arbitrary biological data.

The exploratory report intentionally applies no FDR/TPR gate. It records
threshold-selection evidence and cannot serve as certification evidence.

## Audit of the unexpected FDP

The result was more aggressive than expected and exceeds the predeclared mean
FDP limit. The held-out grid was therefore not run.

For representative seed 61001, a direct locked DESeq2 invocation and the
toolkit agreed on 350 calls, 36 false positives, and 314 true positives. The
maximum absolute differences were `3.91e-14` for log2 fold change, `7.22e-16`
for raw P value, and `1.55e-15` for adjusted P value. An independent Python BH
calculation over the exact 4,976 finite-P family agreed with toolkit FDR to
`2.89e-15`. Thus the excess is not explained by CLI routing, contrast
direction, status inclusion, or BH-family reconstruction.

Truth annotations were internally consistent: exactly 500 genes had nonzero
true log2 fold change, exactly 4,500 had zero true log2 fold change, and the
truth-DE genes were split 250 up and 250 down. Across all 20 exploratory
replicates, 89,495 finite truth-null P values had `Pr(P < 0.05) = 0.06095`
(replicate range 0.05614 to 0.06540). A pooled diagnostic KS statistic was
0.01152; its asymptotic reference P value was approximately `9.7e-11`.

Representative estimated size factors closely tracked compcodeR's true depth
factors (log-scale correlation 0.99908; maximum absolute log-ratio error
0.01443). Median estimated/true dispersion was 1.0796 overall, while false
positive genes had median 0.6680 (interquartile range 0.4969 to 0.7370), which
is consistent with local dispersion underestimation contributing to the
excess. DESeq2 performed no automatic count replacement: the replacement count
was zero and no `replaceCounts` assay existed. Cook's filtering excluded 24
genes in the representative replicate and 544 across the 20-replicate grid.

A separate diagnostic using the earlier edgeR simulation seeds 5001 through
5010 gave direct-DESeq2 mean FDP 0.12366 (range 0.10482 to 0.14058) and mean TPR
0.6496. The behavior is therefore not an accident of the new exploratory seed
range.

## Current decision

The exploratory execution evidence is preserved, the candidate thresholds and
held-out seeds remain unchanged, and no held-out certification report has been
created. Under this frozen scenario, D2 cannot be called passing without a new
reviewed decision about the simulation model or certification claim. Raising
the FDP limits merely to fit these observations would not be an acceptable
gate.
