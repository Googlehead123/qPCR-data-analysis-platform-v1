"""Pins the ΔΔCt/fold-change contract the Streamlit app actually runs.

Ownership as of 2026-07-31 (the dual-copy hazard this file was written for is
gone):

* Parsing, QC, graphing and the ΔΔCt maths live ONLY in `qpcr/`. The app's
  `AnalysisEngine` *subclasses* `qpcr.analysis.AnalysisEngine` and adds just
  `run_full_analysis`, so `calculate_ddct` here and in the package are the same
  function object. `test_monolith_matches_package_fc_errors` therefore no longer
  detects divergence — it is kept as a cheap assertion that the subclass has not
  started shadowing the inherited maths.
* Reporting and export live ONLY in `streamlit qpcr analysis v1.py`. Their
  `qpcr/report.py` / `qpcr/export.py` duplicates were deleted; nothing imported
  them, so fixes applied there never shipped while their tests still passed.

What still earns its keep here is the behavioural contract, not the parity: the
fold-change error columns must exist and carry the per-well fold-change spread
(symmetric, in fold units, since 2026-10-28 — they were asymmetric 2^-x
transforms of the Ct spread before), and the cache-fronted replicate
fold-change wrapper must be wired to the right function with the right argument
order.
"""

from importlib import import_module

import numpy as np


def test_monolith_calculate_ddct_emits_fold_unit_error_bars(
    mock_streamlit, sample_qpcr_raw_data, sample_mapping
):
    """Every error-bar column graph.py and the writers read must exist and be
    the per-well fold-change SD / SEM / CI half-width — symmetric, because the
    spread is now taken in fold units, where the hand sheet takes it."""
    spec = import_module("streamlit qpcr analysis v1")
    result = spec.AnalysisEngine.calculate_ddct(
        data=sample_qpcr_raw_data,
        hk_gene="GAPDH",
        ref_sample="Non-treated",
        excluded_wells=set(),
        excluded_samples=set(),
        sample_mapping=sample_mapping,
    )

    for col in ("FC_Error_Upper", "FC_Error_Lower", "FC_SEM_Upper",
                "FC_SEM_Lower", "FC_CI_Upper", "FC_CI_Lower"):
        assert col in result.columns, col

    spread = result[result["SD"] > 0]
    assert not spread.empty, "fixture should yield replicate spread"
    np.testing.assert_allclose(spread["FC_Error_Upper"], spread["SD"])
    np.testing.assert_allclose(spread["FC_Error_Lower"], spread["SD"])
    np.testing.assert_allclose(spread["FC_SEM_Upper"], spread["SEM"])
    np.testing.assert_allclose(spread["FC_SEM_Lower"], spread["SEM"])
    assert (spread["FC_CI_Upper"] > spread["FC_SEM_Upper"]).all()


def test_monolith_matches_package_fc_errors(
    mock_streamlit, sample_qpcr_raw_data, sample_mapping
):
    """Monolith and qpcr-package calculate_ddct must agree on the FC error bars."""
    spec = import_module("streamlit qpcr analysis v1")
    from qpcr.analysis import AnalysisEngine as PkgEngine

    kw = dict(
        data=sample_qpcr_raw_data,
        hk_gene="GAPDH",
        ref_sample="Non-treated",
        excluded_wells=set(),
        excluded_samples=set(),
        sample_mapping=sample_mapping,
    )
    mono = spec.AnalysisEngine.calculate_ddct(**kw).sort_values(["Target", "Condition"]).reset_index(drop=True)
    pkg = PkgEngine.calculate_ddct(**kw).sort_values(["Target", "Condition"]).reset_index(drop=True)

    for col in ("FC_Error_Upper", "FC_Error_Lower"):
        np.testing.assert_allclose(
            mono[col].to_numpy(dtype=float), pkg[col].to_numpy(dtype=float),
            rtol=1e-9, atol=1e-9, err_msg=f"{col} diverged between monolith and package",
        )


def test_cache_fronted_replicate_fold_changes_matches_direct(
    mock_streamlit, sample_qpcr_raw_data, sample_mapping
):
    """The data-point overlay reads per-replicate fold changes through the
    cache-fronted `get_replicate_fold_changes`; it must return exactly what the
    underlying `compute_replicate_fold_changes` produces. (Under test the
    st.cache_data mock is identity, so this pins the wrapper's wiring — correct
    target function and argument order — not the runtime memoization itself.)
    Regression guard for the CPU-throttle fix that introduced the cache."""
    from pandas.testing import assert_frame_equal

    spec = import_module("streamlit qpcr analysis v1")
    args = (sample_qpcr_raw_data, "GAPDH", "Non-treated", sample_mapping, {})

    direct = spec.AnalysisEngine.compute_replicate_fold_changes(*args)
    fronted = spec.get_replicate_fold_changes(*args)

    assert not fronted.empty, "fixture should yield replicate fold changes"
    assert "Replicate_FC" in fronted.columns
    assert_frame_equal(
        direct.reset_index(drop=True), fronted.reset_index(drop=True)
    )
