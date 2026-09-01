"""Reserved airway same-engine oracle gate; intentionally not active yet."""

import pytest


@pytest.mark.oracle
@pytest.mark.skip(
    reason=(
        "P0 work items 5-6 are incomplete: no locked edgeR backend, airway "
        "fixture, direct oracle script, or benchmark artifact exists yet"
    )
)
def test_airway_edger_ql_oracle_parity_rtol_1e_6() -> None:
    """Compare toolkit and direct edgeR outputs in one locked environment."""

    raise AssertionError("Remove the skip only after the complete oracle gate exists")
