"""Regression tests for the per-gene graph editor's widget/settings coupling.

These drive the REAL app file through ``streamlit.testing.v1.AppTest`` so they
exercise the shipped widgets rather than a copy that can drift.

They must run in a SUBPROCESS: the repo-root ``conftest.py`` replaces
``sys.modules["streamlit"]`` with a MagicMock at import time, which would make
an in-process AppTest meaningless (and is exactly the blind spot that let these
bugs ship — the mock returns a canned value for every widget, so no test could
observe a widget ignoring its ``value=`` argument).

Background — the defect class these guard against: a keyed Streamlit widget
derives its element id from the ``key`` alone, so a recomputed ``value=`` or
``index=`` is ignored once the widget has state. Any code that wrote
``graph_settings`` and expected the widget to follow was a silent no-op.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
APP = REPO / "streamlit qpcr analysis v1.py"

DRIVER = r'''
import json, sys
import pandas as pd
from streamlit.testing.v1 import AppTest

REPO = %(repo)r
APP = %(app)r
sys.path.insert(0, REPO)

GENES = ["COL1A1", "ELN"]
CONDITIONS = ["Control", "TGFb", "Test A"]
G = "COL1A1"


def _processed():
    out = {}
    for gi, g in enumerate(GENES):
        out[g] = pd.DataFrame({
            "Condition": CONDITIONS,
            "Relative_Expression": [1.0, 2.4 + gi, 1.7 + gi],
            "Fold_Change": [1.0, 2.4 + gi, 1.7 + gi],
            "SEM": [0.05, 0.12, 0.09],
            "SD": [0.07, 0.18, 0.13],
            "FC_Error_Upper": [0.08, 0.22, 0.15],
            "FC_Error_Lower": [0.08, 0.22, 0.15],
            "n_replicates": [3, 3, 3],
            "significance": ["", "**", "*"],
            "significance_2": ["", "", ""],
            "significance_3": ["", "", ""],
            "p_value": [1.0, 0.004, 0.03],
        })
    return out


def _raw():
    rows = []
    wells = iter(["%%s%%d" %% (r, c) for r in "ABCDEFGH" for c in range(1, 13)])
    for g in GENES + ["GAPDH"]:
        for s in CONDITIONS:
            for _ in range(3):
                rows.append({"Well": next(wells), "Target": g, "Sample": s, "CT": 22.0})
    return pd.DataFrame(rows)


def app():
    at = AppTest.from_file(APP, default_timeout=300)
    ss = at.session_state
    ss["data"] = _raw()
    ss["raw_data"] = _raw()
    ss["processed_data"] = _processed()
    ss["selected_efficacy"] = "패릭" if False else "탄력"
    ss["hk_gene"] = "GAPDH"
    ss["analysis_ref_condition"] = "Control"
    ss["analysis_cmp_condition"] = "Control"
    ss["sample_mapping"] = {c: {"condition": c, "include": True,
                                "group": "Treatment"} for c in CONDITIONS}
    ss["sample_order"] = list(CONDITIONS)
    ss["excluded_wells"] = {}
    ss["excluded_samples"] = []
    ss["gene_display_names"] = {}
    ss["selected_gene_idx"] = 0
    at.run()
    return at


def gs(at):
    return dict(at.session_state["graph_settings"])


out = {}

# --- baseline: the app renders at all with results present
at = app()
out["baseline_exceptions"] = [e.message for e in at.exception]

# --- size preset buttons must actually move the dimensions
at.button(key="sp_PPT Half_%%s" %% G).click().run()
g = gs(at)
out["preset_half"] = {
    "width": g.get("%%s_figure_width" %% G),
    "height": g.get("%%s_figure_height" %% G),
    "slider_w": at.slider(key="fig_w_%%s" %% G).value,
    "slider_h": at.slider(key="fig_h_%%s" %% G).value,
    "exceptions": [e.message for e in at.exception],
}
at.button(key="sp_Wide_%%s" %% G).click().run()
g = gs(at)
out["preset_wide"] = {"width": g.get("%%s_figure_width" %% G),
                      "height": g.get("%%s_figure_height" %% G)}

# --- a hand-picked bar colour must survive the next interaction
at = app()
at.color_picker(key="cp_%%s_TGFb" %% G).set_value("#FF0000").run()
g = gs(at)
out["color_immediate"] = {
    "preset": g.get("%%s_color_preset" %% G),
    "bar": (g.get("bar_colors_per_sample") or {}).get("%%s_TGFb" %% G),
}
# touch an unrelated control, then re-read
at.slider(key="gf_%%s" %% G).set_value(18).run()
g = gs(at)
out["color_after_other_interaction"] = {
    "preset": g.get("%%s_color_preset" %% G),
    "bar": (g.get("bar_colors_per_sample") or {}).get("%%s_TGFb" %% G),
    "picker": at.color_picker(key="cp_%%s_TGFb" %% G).value,
    "exceptions": [e.message for e in at.exception],
}

# --- Reset this gene must restore defaults
at = app()
at.slider(key="fig_w_%%s" %% G).set_value(19.5).run()
at.selectbox(key="preset_%%s" %% G).select("Coral").run()
at.slider(key="gf_%%s" %% G).set_value(22).run()
g = gs(at)
out["before_reset"] = {"width": g.get("%%s_figure_width" %% G),
                       "preset": g.get("%%s_color_preset" %% G),
                       "font": g.get("%%s_font_size" %% G)}
at.button(key="reset_all_%%s" %% G).click().run()
g = gs(at)
out["after_reset"] = {
    "width": g.get("%%s_figure_width" %% G),
    "height": g.get("%%s_figure_height" %% G),
    "preset": g.get("%%s_color_preset" %% G),
    "font": g.get("%%s_font_size" %% G),
    "slider_w": at.slider(key="fig_w_%%s" %% G).value,
    "widget_preset": at.selectbox(key="preset_%%s" %% G).value,
    "exceptions": [e.message for e in at.exception],
}

# --- changing size then forcing a FULL rerun must not raise
#     (the old code called st.rerun(scope="fragment") from a fragment during a
#      full script run, which Streamlit rejects outright)
at = app()
at.slider(key="fig_w_%%s" %% G).set_value(24.0).run()
at.button(key="gene_btn_ELN").click().run()          # gene switch = full rerun
out["full_rerun_after_size_change"] = {
    "exceptions": [e.message for e in at.exception],
    "selected_idx": at.session_state["selected_gene_idx"],
}

print("@@JSON@@" + json.dumps(out, default=str))
'''


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    script = tmp_path_factory.mktemp("editor") / "driver.py"
    script.write_text(DRIVER % {"repo": str(REPO), "app": str(APP)})
    proc = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True, text=True, timeout=900, cwd=str(REPO),
    )
    marker = "@@JSON@@"
    if marker not in proc.stdout:
        pytest.fail(
            "AppTest driver produced no result payload.\n"
            f"stdout tail:\n{proc.stdout[-4000:]}\n\nstderr tail:\n{proc.stderr[-4000:]}"
        )
    return json.loads(proc.stdout.split(marker, 1)[1].splitlines()[0])


def test_app_renders_with_results_present(results):
    assert results["baseline_exceptions"] == []


def test_size_preset_button_applies_its_dimensions(results):
    """A size preset used to be overwritten by the Width/Height sliders read
    immediately afterwards, so every preset silently resolved to 28x16."""
    half = results["preset_half"]
    assert half["exceptions"] == []
    assert (half["width"], half["height"]) == (14, 10), half
    assert (half["slider_w"], half["slider_h"]) == (14.0, 10.0), half
    wide = results["preset_wide"]
    assert (wide["width"], wide["height"]) == (32, 14), wide


def test_custom_bar_color_survives_the_next_interaction(results):
    """The colour-preset selectbox used to re-assert its stored preset on the
    next rerun, silently repainting a hand-picked bar back to the preset tone."""
    assert results["color_immediate"]["bar"] == "#FF0000"
    assert results["color_immediate"]["preset"] == "Custom"
    after = results["color_after_other_interaction"]
    assert after["exceptions"] == []
    assert after["preset"] == "Custom", after
    assert after["bar"] == "#FF0000", after
    assert after["picker"] == "#FF0000", after


def test_reset_this_gene_restores_defaults(results):
    """Reset used to pop graph_settings only; the keyed widgets kept their values
    and repopulated it on the next run, making the button a no-op."""
    before = results["before_reset"]
    assert before["preset"] == "Coral" and before["font"] == 22, before
    after = results["after_reset"]
    assert after["exceptions"] == []
    assert (after["width"], after["height"]) == (28.0, 16.0), after
    assert after["preset"] == "Classic", after
    assert after["widget_preset"] == "Classic", after
    assert after["font"] in (None, 14), after
    assert after["slider_w"] == 28.0, after


def test_full_rerun_after_a_size_change_does_not_raise(results):
    """Adjusting the width and then switching gene (a full rerun) used to raise
    StreamlitAPIException from st.rerun(scope="fragment")."""
    res = results["full_rerun_after_size_change"]
    assert res["exceptions"] == [], res
    assert res["selected_idx"] == 1, res
