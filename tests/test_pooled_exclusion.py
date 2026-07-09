"""Regression test for the pooled-Condition well-exclusion bug.

When one Condition pools multiple raw samples, `calculate_ddct` must apply
per-gene-sample well exclusions to EVERY pooled sample, not just the first one
(`.iloc[0]`). The old code silently ignored an exclusion on any sample but the
first, so the ΔΔCt point estimate disagreed with the p-value/scatter (which are
computed row-wise). Tested against BOTH the monolith (what the app runs) and the
qpcr package copy.
"""

from importlib import import_module

import pandas as pd
import pytest


def _pooled_raw():
    """COL1A1 target + GAPDH HK. 'Treatment' pools raw samples D1 and D2.
    Well B6 (sample D2, COL1A1) is an outlier CT=30 that will be excluded."""
    rows = []

    def add(well, sample, target, ct):
        rows.append({"Well": well, "Sample": sample, "Target": target, "CT": ct})

    # Reference condition "Ctrl"
    for w, ct in [("C1", 20.0), ("C2", 20.0), ("C3", 20.0)]:
        add(w, "Ctrl", "GAPDH", ct)
    for w, ct in [("C4", 25.0), ("C5", 25.0), ("C6", 25.0)]:
        add(w, "Ctrl", "COL1A1", ct)

    # Treatment, sample D1
    for w, ct in [("A1", 20.0), ("A2", 20.0), ("A3", 20.0)]:
        add(w, "D1", "GAPDH", ct)
    for w, ct in [("A4", 24.0), ("A5", 24.0), ("A6", 24.0)]:
        add(w, "D1", "COL1A1", ct)

    # Treatment, sample D2 (B6 is a bad well: CT=30)
    for w, ct in [("B1", 20.0), ("B2", 20.0), ("B3", 20.0)]:
        add(w, "D2", "GAPDH", ct)
    for w, ct in [("B4", 24.0), ("B5", 24.0), ("B6", 30.0)]:
        add(w, "D2", "COL1A1", ct)

    return pd.DataFrame(rows)


_MAPPING = {
    "Ctrl": {"condition": "Ctrl", "group": "Negative Control", "include": True},
    "D1": {"condition": "Treatment", "group": "Treatment", "include": True},
    "D2": {"condition": "Treatment", "group": "Treatment", "include": True},
}


def _run(engine):
    return engine.calculate_ddct(
        data=_pooled_raw(),
        hk_gene="GAPDH",
        ref_sample="Ctrl",
        # Exclude the outlier well B6, keyed to the SECOND pooled sample (D2).
        excluded_wells={("COL1A1", "D2"): {"B6"}},
        excluded_samples=set(),
        sample_mapping=_MAPPING,
    )


def _treatment_col1a1_mean(result):
    row = result[(result["Target"] == "COL1A1") & (result["Condition"] == "Treatment")]
    assert len(row) == 1, f"expected one Treatment/COL1A1 row, got {len(row)}"
    return float(row["Target_Ct_Mean"].iloc[0])


def test_monolith_pooled_exclusion_applies_to_all_samples(mock_streamlit):
    spec = import_module("streamlit qpcr analysis v1")
    mean = _treatment_col1a1_mean(_run(spec.AnalysisEngine))
    # Correct: B6 (30.0) excluded -> mean of five 24.0 wells = 24.0.
    # Buggy .iloc[0]: B6 kept -> mean of [24,24,24,24,24,30] = 25.0.
    assert mean == pytest.approx(24.0), (
        f"Treatment/COL1A1 Ct mean = {mean}; excluded well B6 (D2) was not "
        f"applied (25.0 means the pooled-exclusion bug is present)."
    )


def test_package_pooled_exclusion_applies_to_all_samples(mock_streamlit):
    from qpcr.analysis import AnalysisEngine
    mean = _treatment_col1a1_mean(_run(AnalysisEngine))
    assert mean == pytest.approx(24.0)
