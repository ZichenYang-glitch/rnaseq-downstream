# Airway edgeR same-engine oracle

This directory contains the same-engine backend-kernel oracle required by P0:

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
