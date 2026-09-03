# Airway same-engine oracles

This directory contains two locked backend-kernel oracles.

## P0 gene-level edgeR oracle

- the toolkit edgeR v4 quasi-likelihood path and the direct oracle script must
  run in the same locked environment;
- airway is extracted from the locked data package at runtime and both routes
  consume those exact bytes;
- coefficients, log fold changes, quasi-likelihood F statistics, raw P values,
  and BH FDR values must meet the declared `rtol=1e-6`, `atol=1e-10` gate;
- the comparison must write a machine-readable benchmark report.

There are no checked-in expected coefficient values. The direct R script is an
independent literal edgeR implementation and the toolkit route uses the private
benchmark-kernel adapter. This gate does not assert featureCounts origin and
does not certify the public input route.

## C2 edgeR/limma pathway oracle

The C2 oracle uses the same airway counts and paired `~ cell + dex` design. It
creates a frozen local GMT and one-to-one synthetic symbol annotation, including
below-minimum and unmapped negative controls. These sets are test fixtures, not
curated biological pathways; sets chosen from crude airway group means exist to
exercise numerical branches and do not provide independent biological truth.

The independent R script and the complete toolkit backend path each build the
same locked DGEGLM route, then run one multi-set mroast call, fry, and camera.
The gate requires the full set/method/hypothesis grid, mapping/status metadata,
direction, gene counts, mroast proportions, raw/effective camera correlations,
VIF, P values, and FDR values to agree. Numeric fields use `rtol=1e-6` and
`atol=1e-10`; Python also recalculates BH within each
contrast/method/hypothesis family to `1e-12`.

This is a same-engine regression oracle. It can catch routing, mapping, dispatch,
or numerical drift between the independently written paths, but it does not
validate limma's methods against an external statistical implementation. Like
the P0 oracle, it enters through the private benchmark adapter and does not
certify featureCounts origin or the public input route.
