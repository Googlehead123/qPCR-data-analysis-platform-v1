"""Guards against monolith / qpcr-package divergence for computations that the
Streamlit app actually runs from the monolith copy.

The app executes the INLINE classes in `streamlit qpcr analysis v1.py`, not the
`qpcr/` package. Bug fixes applied only to the package silently never ship. This
suite pins the parity that matters at runtime.
"""

from importlib import import_module

import numpy as np


def test_monolith_calculate_ddct_emits_asymmetric_fc_error_bars(
    mock_streamlit, sample_qpcr_raw_data, sample_mapping
):
    """The monolith's calculate_ddct must produce fold-change-domain (Livak)
    error columns, so graphs use asymmetric bars instead of the symmetric
    SD/SEM fallback. This fix previously lived only in qpcr/analysis.py."""
    spec = import_module("streamlit qpcr analysis v1")
    result = spec.AnalysisEngine.calculate_ddct(
        data=sample_qpcr_raw_data,
        hk_gene="GAPDH",
        ref_sample="Non-treated",
        excluded_wells=set(),
        excluded_samples=set(),
        sample_mapping=sample_mapping,
    )

    assert "FC_Error_Upper" in result.columns
    assert "FC_Error_Lower" in result.columns

    # For rows with real replicate spread (sd > 0), the 2^-x transform makes the
    # upper and lower fold-change errors asymmetric, and both must be >= 0.
    spread = result[result["SD"] > 0]
    assert not spread.empty, "fixture should yield replicate spread"
    assert (spread["FC_Error_Upper"] >= -1e-9).all()
    assert (spread["FC_Error_Lower"] >= -1e-9).all()
    # At least one row is meaningfully asymmetric.
    assert (np.abs(spread["FC_Error_Upper"] - spread["FC_Error_Lower"]) > 1e-6).any()


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
