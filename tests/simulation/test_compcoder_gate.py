"""Reserved compcodeR FDR/TPR gate; intentionally not active yet."""

import pytest


@pytest.mark.simulation
@pytest.mark.skip(
    reason=(
        "P0 work item 6 is incomplete: the locked compcodeR scenario, fixed "
        "seed, FDR/TPR thresholds, and benchmark JSON artifact do not exist yet"
    )
)
def test_compcoder_negative_binomial_fdr_tpr_gate() -> None:
    """Enforce the future fixed FDR/TPR acceptance thresholds."""

    raise AssertionError("Remove the skip only after the complete simulation gate exists")
