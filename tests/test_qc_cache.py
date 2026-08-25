"""Guards for the cache-fronted QC helpers in the monolith.

`st.tabs` executes every tab body on every rerun, so the QC dashboard's summary
and replicate table used to recompute on each widget interaction anywhere in the
app. They are now memoized.

Two things have to hold:

1. **Wiring** — each wrapper must return exactly what the underlying
   `QualityControl` method returns. (Under test `st.cache_data` is mocked to the
   identity decorator, so these tests pin the wiring, not the memoization.)
2. **Key completeness** — the QC thresholds change the helpers' output, so they
   must be part of every QC cache key or moving one would keep serving the old
   counts. That is the failure this suite exists for.

   They used to be MUTABLE CLASS ATTRIBUTES that the QC tab's widgets
   reassigned, which is process state shared by every browser session on the
   container. They now live in the widgets' own session_state keys and are
   passed to the helpers as arguments; the class attributes are immutable
   defaults. Three tests below pin that.
"""

from importlib import import_module

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


# NOTE: always reach QualityControl through the app module (`app.QualityControl`),
# never via a module-level `from qpcr.quality_control import ...`. The autouse
# mock_streamlit fixture drops every `qpcr*` entry from sys.modules per test, so a
# top-level import here would bind a DIFFERENT class object than the one the
# monolith uses — and mutating its thresholds would silently affect nothing.


@pytest.fixture
def qc_data():
    rows = []
    well = 1
    for sample in ("Non-treated", "Treatment1"):
        for target in ("GAPDH", "COL1A1"):
            for ct in (20.0, 20.1, 20.2):
                rows.append({
                    "Well": f"A{well}", "Sample": sample,
                    "Target": target, "CT": ct + (0.5 if target == "COL1A1" else 0.0),
                })
                well += 1
    return pd.DataFrame(rows)


@pytest.fixture
def app(mock_streamlit):
    return import_module("streamlit qpcr analysis v1")


@pytest.fixture(autouse=True)
def _restore_thresholds(app):
    """The tests below mutate class-level thresholds; put them back."""
    qc = app.QualityControl
    saved = (qc.CT_HIGH_THRESHOLD, qc.CT_LOW_THRESHOLD,
             qc.CV_THRESHOLD, qc.HK_VARIATION_THRESHOLD)
    yield
    (qc.CT_HIGH_THRESHOLD, qc.CT_LOW_THRESHOLD,
     qc.CV_THRESHOLD, qc.HK_VARIATION_THRESHOLD) = saved


# ---- wiring --------------------------------------------------------------

def test_qc_summary_wrapper_matches_direct_call(app, qc_data):
    assert app.get_qc_summary_stats(qc_data, set()) == \
        app.QualityControl.get_qc_summary_stats(qc_data, set())


def test_qc_summary_wrapper_honours_exclusions(app, qc_data):
    excluded = {"A1", "A2"}
    fronted = app.get_qc_summary_stats(qc_data, excluded)
    assert fronted == app.QualityControl.get_qc_summary_stats(qc_data, excluded)
    # and the exclusions genuinely changed the answer
    assert fronted != app.get_qc_summary_stats(qc_data, set())


def test_qc_summary_treats_none_and_empty_set_alike(app, qc_data):
    """Call sites pass both forms; the helper coerces `None` to an empty set."""
    assert app.get_qc_summary_stats(qc_data, None) == \
        app.get_qc_summary_stats(qc_data, set())


def test_replicate_stats_wrapper_matches_direct_call(app, qc_data):
    assert_frame_equal(
        app.get_replicate_stats(qc_data).reset_index(drop=True),
        app.QualityControl.get_replicate_stats(qc_data).reset_index(drop=True),
    )


def test_plate_heatmap_is_gone(app):
    """The 96-well heatmap was removed: an 8x12 grid cannot represent the
    concatenated multi-plate frames every real export produces."""
    assert not hasattr(app, "get_plate_heatmap")
    assert not hasattr(app.QualityControl, "create_plate_heatmap")


# ---- cache-key completeness ----------------------------------------------

@pytest.mark.parametrize("attr,new_value", [
    ("CT_HIGH_THRESHOLD", 30.0),
    ("CT_LOW_THRESHOLD", 12.0),
    ("CV_THRESHOLD", 0.25),
    ("HK_VARIATION_THRESHOLD", 2.5),
])
def test_threshold_change_changes_the_cache_key(app, attr, new_value):
    """Each slider-backed threshold must move the key, or QC output goes stale."""
    before = app._qc_threshold_key()
    setattr(app.QualityControl, attr, new_value)
    assert app._qc_threshold_key() != before, (
        f"{attr} is missing from the QC cache key — changing that slider "
        "would keep serving the previous QC results"
    )


def test_threshold_key_is_stable_when_nothing_changes(app):
    assert app._qc_threshold_key() == app._qc_threshold_key()


def test_excluded_key_is_order_independent(app):
    """Sets iterate in unstable order; the key must not."""
    assert app._excluded_key({"A3", "A1", "A2"}) == app._excluded_key({"A1", "A2", "A3"})


def test_excluded_key_distinguishes_different_contents(app):
    assert app._excluded_key({"A1", "A2"}) != app._excluded_key({"A1", "A3"})


def test_excluded_key_normalises_empty_forms(app):
    assert app._excluded_key(None) == app._excluded_key(set()) == ()


def test_ct_threshold_actually_reaches_the_summary(app, qc_data):
    """End-to-end: a threshold move must change the numbers the wrapper returns.

    Pairs with the key test above — together they show the setting both changes
    the key and changes the result, so a stale hit would be observable.

    Drives the SESSION store, not the class attribute. The widgets used to
    assign onto QualityControl's class attributes, which is process state shared
    by every browser session on the container; the threshold now lives in the
    widget's own session_state key and is passed to the helpers as an argument.
    """
    app.st.session_state["qc_settings_ct_high"] = 100.0
    lenient = app.get_qc_summary_stats(qc_data, set())
    app.st.session_state["qc_settings_ct_high"] = 15.0
    strict = app.get_qc_summary_stats(qc_data, set())
    assert lenient["high_ct_count"] != strict["high_ct_count"], (
        "moving the session threshold must change the counts the wrapper "
        "returns; if these match, the value is not reaching the helper"
    )


def test_thresholds_are_not_process_global(app):
    """The class attributes must stay at their documented defaults.

    This is the actual defect: a class attribute is PROCESS state, so one
    operator's QC settings became every concurrent operator's and survived a
    refresh. Setting a session threshold must not touch the class.
    """
    qc = app.QualityControl
    before = (qc.CT_HIGH_THRESHOLD, qc.CT_LOW_THRESHOLD,
              qc.CV_THRESHOLD, qc.HK_VARIATION_THRESHOLD)

    app.st.session_state["qc_settings_ct_high"] = 26.0
    app.st.session_state["qc_settings_cv"] = 12.0  # widget is in percent

    resolved = app.qc_thresholds()
    assert resolved["ct_high"] == 26.0
    assert resolved["cv"] == pytest.approx(0.12), "percent must be scaled"

    after = (qc.CT_HIGH_THRESHOLD, qc.CT_LOW_THRESHOLD,
             qc.CV_THRESHOLD, qc.HK_VARIATION_THRESHOLD)
    assert before == after == (35.0, 10.0, 0.05, 1.0), (
        "the class attributes are immutable defaults; a session threshold must "
        "not write through to them"
    )


def test_unset_thresholds_fall_back_to_the_defaults(app):
    """Before the QC Settings widgets have ever rendered, keys do not exist."""
    for key in ("qc_settings_ct_high", "qc_settings_ct_low",
                "qc_settings_cv", "qc_settings_hk_var"):
        app.st.session_state.pop(key, None)
    assert app.qc_thresholds() == {
        "ct_high": 35.0, "ct_low": 10.0, "cv": 0.05, "hk_variation": 1.0,
    }
