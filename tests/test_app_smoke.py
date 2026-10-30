"""End-to-end runtime smoke tests for the Streamlit app.

Unit tests import classes directly; these instead run the *actual script*
headless via Streamlit's AppTest, so they catch runtime wiring bugs (session
state, tab rendering, the fragment editor, the Overview verdict) that unit tests
miss. Self-contained — no external data files needed.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Absolute, because AppTest.from_file resolves a relative path against the
# calling file (tests/) since Streamlit 1.61 (streamlit#16122); before that it
# tried the working directory first, which is what the bare name relied on.
SCRIPT = str(Path(__file__).resolve().parents[1] / "streamlit qpcr analysis v1.py")


@pytest.fixture(autouse=True)
def mock_streamlit():
    """Override conftest's autouse Streamlit mock.

    AppTest needs the REAL streamlit module and fresh app/package modules, so we
    evict the mock and any mock-tainted caches before each test here.

    ``streamlit_sortables`` goes too: it binds its component to whichever
    ``streamlit`` was loaded when it was first imported. Kept across tests it
    still points at the copy evicted here. Streamlit 1.53 tolerated that, but
    1.66 reads each run's state from a per-copy context variable that the
    stale copy never initialises, so ``sort_items`` raises
    ``FragmentThreadState not initialized`` on the Mapping tab.
    """
    for name in [
        k for k in list(sys.modules)
        if k == "streamlit" or k.startswith("streamlit.")
        or k == "streamlit_sortables" or k.startswith("streamlit_sortables.")
        or k.startswith("qpcr") or k == "streamlit qpcr analysis v1"
    ]:
        del sys.modules[name]
    import streamlit  # noqa: F401  — real module, reloaded from disk
    yield


def _triplicate(sample, target, base):
    return [{"Well": f"{target[:3]}{sample[:2]}{r}", "Sample": sample,
             "Target": target, "CT": base + 0.05 * r} for r in range(3)]


def _synthetic_data():
    """Small photoaging-style dataset: ACTIN (HK) + COL1A1 + MMP1, 2 conditions."""
    rows = []
    for sample, shift in [("Non-treated", 0.0), ("Active", 0.0)]:
        rows += _triplicate(sample, "ACTIN", 18.0)
        rows += _triplicate(sample, "COL1A1", 22.0 + (-0.6 if sample == "Active" else 0))
        rows += _triplicate(sample, "MMP1", 24.0 + (0.7 if sample == "Active" else 0))
    return pd.DataFrame(rows)


def test_full_pipeline_runs_without_exceptions():
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(SCRIPT, default_timeout=90)
    at.session_state["data"] = _synthetic_data()
    at.session_state["hk_gene"] = "ACTIN"
    at.session_state["selected_efficacy"] = "광노화"
    at.run()
    assert not at.exception, at.exception

    run_btns = [b for b in at.button if "Run Full Analysis" in (b.label or "")]
    assert run_btns, "Run Full Analysis button not found"
    run_btns[0].click().run()
    assert not at.exception, at.exception
    assert at.session_state["processed_data"], "analysis produced no results"

    # Re-render with results: exercises Overview + Graphs fragment + Export.
    at.run()
    assert not at.exception, at.exception
    graphs = at.session_state["graphs"]
    assert set(at.session_state["processed_data"]) <= set(graphs)


def test_pvalues_default_to_two_tailed_student():
    """With the T-test radio untouched, the shipped p-value is the hand sheet's:
    two-tailed Student's t-test on each well's fold change.

    Per-well ΔCt: Non-treated 4.00/4.05/3.95, Active 1.9/2.5/1.6 — the
    control-vs-A ΔCt of tests/test_fold_change_method.py shifted by a constant.
    That scales every fold change by one factor, which leaves t and p
    unchanged, so the expected value is the one pinned there by hand
    calculation (Welch would give 0.0478; the old ΔCt-domain test 0.00168).
    """
    from streamlit.testing.v1 import AppTest
    rows = []
    for sample, col1a1 in [("Non-treated", [22.0, 22.05, 21.95]),
                           ("Active", [19.9, 20.5, 19.6])]:
        for r in range(3):
            rows.append({"Well": f"AC{sample[:2]}{r}", "Sample": sample,
                         "Target": "ACTIN", "CT": 18.0})
            rows.append({"Well": f"CO{sample[:2]}{r}", "Sample": sample,
                         "Target": "COL1A1", "CT": col1a1[r]})
    at = AppTest.from_file(SCRIPT, default_timeout=90)
    at.session_state["data"] = pd.DataFrame(rows)
    at.session_state["hk_gene"] = "ACTIN"
    at.session_state["selected_efficacy"] = "광노화"
    at.run()
    assert not at.exception, at.exception
    assert at.radio(key="ttest_type").value == "student"

    run_btns = [b for b in at.button if "Run Full Analysis" in (b.label or "")]
    run_btns[0].click().run()
    assert not at.exception, at.exception
    col1a1 = at.session_state["processed_data"]["COL1A1"]
    row = col1a1[col1a1["p_value"].notna()].iloc[0]
    assert row["p_value"] == pytest.approx(0.0117003895956, rel=1e-9)
    assert row["significance"] == "*"
    assert row["stat_test_used"].startswith("Student t-test")


def _gene_df(up):
    folds = {"UVB only": 1.0, "TGFβ": 1.8, "Active": 1.56} if up \
        else {"UVB only": 1.0, "TGFβ": 0.45, "Active": 0.62}
    return pd.DataFrame([
        {"Condition": c, "Fold_Change": f, "Relative_Expression": f,
         "p_value": (np.nan if c == "UVB only" else 0.002),
         "significance": ("" if c == "UVB only" else ("***" if c == "TGFβ" else "**")),
         "SEM": 0.05, "SD": 0.08, "Group": "Treatment", "n_replicates": 3,
         "Target_Ct_Mean": 25.0}
        for c, f in folds.items()
    ])


def test_overview_verdict_and_benchmark():
    """Overview grades expected direction and auto-selects the benchmark."""
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file(SCRIPT, default_timeout=60)
    at.session_state["processed_data"] = {"COL1A1": _gene_df(up=True), "MMP1": _gene_df(up=False)}
    at.session_state["selected_efficacy"] = "광노화"   # expected: COL1A1 up, MMP1 down
    at.session_state["analysis_ref_condition"] = "UVB only"
    at.session_state["hk_gene"] = "ACTIN"
    at.session_state["sample_order"] = ["s1", "s2", "s3"]
    at.session_state["sample_mapping"] = {
        "s1": {"condition": "UVB only", "group": "Negative Control", "include": True},
        "s2": {"condition": "TGFβ", "group": "Positive Control", "include": True},
        "s3": {"condition": "Active", "group": "Treatment", "include": True},
    }
    at.run()
    assert not at.exception, at.exception

    assert at.session_state["overview_benchmark"] == "TGFβ"   # auto-matched positive control
    metrics = {m.label: m.value for m in at.metric}
    assert metrics.get("Markers as expected") == "2/2"
    # "Avg % of benchmark" until caca43e, which summarises with the MEDIAN
    # instead: these are ratios, and a mean let a negative ratio from a marker
    # moving the other way cancel a genuine positive.
    assert metrics.get("Median % of benchmark", "—") != "—"


def test_overview_measures_from_the_reference_mean_not_from_one():
    """Since 2026-10-28 the reference bar is the simple mean of its per-well
    fold changes — ~1, not exactly 1 (wells at FC 2, 1, 0.5 average 1.1667).
    The benchmark % and the up/down call must be measured from that reported
    baseline: (1.75 - 1.1667) / (2.3333 - 1.1667) = 50%, where a baseline of 1
    reported 56%."""
    from streamlit.testing.v1 import AppTest

    folds = {"UVB only": 7 / 6, "TGFβ": 7 / 3, "Active": 1.75}
    df = pd.DataFrame([
        {"Condition": c, "Fold_Change": f, "Relative_Expression": f,
         "p_value": (np.nan if c == "UVB only" else 0.002),
         "significance": ("" if c == "UVB only" else "**"),
         "SEM": 0.05, "SD": 0.08, "Group": "Treatment", "n_replicates": 3,
         "Target_Ct_Mean": 25.0}
        for c, f in folds.items()
    ])
    at = AppTest.from_file(SCRIPT, default_timeout=60)
    at.session_state["processed_data"] = {"COL1A1": df}
    at.session_state["selected_efficacy"] = "광노화"   # expected: COL1A1 up
    at.session_state["analysis_ref_condition"] = "UVB only"
    at.session_state["hk_gene"] = "ACTIN"
    at.session_state["sample_order"] = ["s1", "s2", "s3"]
    at.session_state["sample_mapping"] = {
        "s1": {"condition": "UVB only", "group": "Negative Control", "include": True},
        "s2": {"condition": "TGFβ", "group": "Positive Control", "include": True},
        "s3": {"condition": "Active", "group": "Treatment", "include": True},
    }
    at.run()
    assert not at.exception, at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics.get("Median % of benchmark") == "50%"
    assert metrics.get("Markers as expected") == "1/1"


def _validation_messages(ct_target):
    """Warnings and successes the Upload tab's Data Validation block shows."""
    from streamlit.testing.v1 import AppTest

    rows = []
    for sample in ("Non-treated", "Active"):
        rows += _triplicate(sample, "ACTIN", 18.0)
        rows += _triplicate(sample, "COL1A1", ct_target)
    at = AppTest.from_file(SCRIPT, default_timeout=90)
    at.session_state["data"] = pd.DataFrame(rows)
    at.session_state["hk_gene"] = "ACTIN"
    at.session_state["selected_efficacy"] = "광노화"
    at.run()
    assert not at.exception, at.exception
    return [w.value for w in at.warning], [s.value for s in at.success]


def test_upload_validation_high_ct_uses_the_38_cutoff():
    """Ct 35-38 is no longer a validation warning; Ct above 38 is, naming 38.

    The block is inline tab code, so the smallest reliable check is the real
    script under AppTest (the golden harness never executes it).
    """
    warnings, successes = _validation_messages(37.0)
    assert not [w for w in warnings if "High CT Values" in w], warnings
    assert any("Data validation passed" in s for s in successes), successes

    warnings, successes = _validation_messages(38.5)
    high = [w for w in warnings if "High CT Values" in w]
    assert len(high) == 1 and "6 wells" in high[0] and "CT > 38 " in high[0], warnings
    assert not any("Data validation passed" in s for s in successes)
