"""Guards for the figure-rebuild gate in `build_all_figures`.

`build_all_figures` skips rebuilding a gene's Plotly figure when a fingerprint of
its inputs is unchanged. That saves most of the per-rerun CPU on Streamlit Cloud,
but it is only safe while the fingerprint is COMPLETE: if an input is missing from
it, the user changes that input and the chart silently keeps the old value.

So the bulk of this suite is one assertion repeated over every input that can
change a figure — including the ones `create_gene_graph` reads straight out of
`st.session_state` rather than taking as arguments (sample_mapping,
error_bar_mode, `f"{gene}_bar_settings"`, the analysis_cmp_condition names), which
are invisible in its signature and therefore easy to forget.
"""

from importlib import import_module

import pandas as pd
import pytest


GENES = ["COL1A1", "ELN"]
PRIMARY = "COL1A1"
OTHER = "ELN"


def _gene_df(gene, expr=(1.0, 2.5)):
    return pd.DataFrame({
        "Target": [gene, gene],
        "Condition": ["Non-treated", "Treatment1"],
        "Group": ["Negative Control", "Treatment"],
        "Relative_Expression": list(expr),
        "SEM": [0.05, 0.08],
        "n_replicates": [3, 3],
    })


def _raw_df(ct=25.0):
    return pd.DataFrame({
        "Well": ["A1", "A2", "A3", "A4"],
        "Sample": ["Non-treated", "Non-treated", "Treatment1", "Treatment1"],
        "Target": ["GAPDH", PRIMARY, "GAPDH", PRIMARY],
        "CT": [18.5, ct, 18.4, ct - 1.5],
    })


@pytest.fixture
def app(mock_streamlit):
    """The monolith with a minimal, fully-populated session_state."""
    spec = import_module("streamlit qpcr analysis v1")
    ss = spec.st.session_state
    ss.clear()
    ss["processed_data"] = {g: _gene_df(g) for g in GENES}
    ss["graphs"] = {}
    ss["_graph_signatures"] = {}
    ss["graph_settings"] = {
        "font_size": 14,
        "bar_gap": 0.45,
        "bar_opacity": 0.85,
        "figure_width": 28,
        "figure_height": 16,
        "plot_bgcolor": "#FFFFFF",
        "bar_colors_per_sample": {},
    }
    ss["gene_display_names"] = {}
    ss["sample_order"] = ["Non-treated", "Treatment1"]
    ss["sample_mapping"] = {
        "Non-treated": {"condition": "Non-treated", "include": True},
        "Treatment1": {"condition": "Treatment1", "include": True},
    }
    ss["error_bar_mode"] = "livak_sd"
    ss["analysis_ref_condition"] = "Non-treated"
    ss["analysis_cmp_condition"] = "Treatment1"
    ss["analysis_cmp_condition_2"] = ""
    ss["analysis_cmp_condition_3"] = ""
    ss["data"] = _raw_df()
    ss["hk_gene"] = "GAPDH"
    ss["excluded_wells"] = {}
    return spec


@pytest.fixture
def builds(app, monkeypatch):
    """Record which genes actually got rebuilt, without real Plotly work."""
    calls = []

    def _fake_build(gene, gene_data, efficacy_config):
        calls.append(gene)
        return object()

    monkeypatch.setattr(app, "build_gene_figure", _fake_build)
    return calls


# ---- the gate actually gates ---------------------------------------------

def test_first_pass_builds_every_gene(app, builds):
    app.build_all_figures(GENES, {})
    assert sorted(builds) == sorted(GENES)


def test_unchanged_inputs_skip_every_rebuild(app, builds):
    """The whole point: an idle rerun must cost zero figure builds."""
    app.build_all_figures(GENES, {})
    builds.clear()
    for _ in range(5):
        app.build_all_figures(GENES, {})
    assert builds == []


def test_figures_are_still_present_after_a_skipped_pass(app, builds):
    app.build_all_figures(GENES, {})
    first = dict(app.st.session_state["graphs"])
    app.build_all_figures(GENES, {})
    assert app.st.session_state["graphs"] == first
    assert set(first) == set(GENES)


def test_one_gene_tweak_does_not_rebuild_the_others(app, builds):
    """Per-gene isolation — editing COL1A1 must not re-cost ELN."""
    app.build_all_figures(GENES, {})
    builds.clear()
    app.st.session_state["graph_settings"][f"{PRIMARY}_bar_gap"] = 0.2
    app.build_all_figures(GENES, {})
    assert builds == [PRIMARY]


def test_cleared_graphs_are_rebuilt_even_if_signature_matches(app, builds):
    """A stale signature must never leave the export path without a figure."""
    app.build_all_figures(GENES, {})
    builds.clear()
    app.st.session_state["graphs"] = {}
    app.build_all_figures(GENES, {})
    assert sorted(builds) == sorted(GENES)


# ---- every input must invalidate -----------------------------------------

def _set_gs(ss, key, value):
    ss["graph_settings"][key] = value


MUTATIONS = {
    # per-gene chart data
    "processed data values":
        lambda ss: ss["processed_data"].__setitem__(PRIMARY, _gene_df(PRIMARY, (1.0, 9.9))),
    # global settings
    "global font size": lambda ss: _set_gs(ss, "font_size", 22),
    "global figure width": lambda ss: _set_gs(ss, "figure_width", 40),
    "global y log scale": lambda ss: _set_gs(ss, "y_log_scale", True),
    "global show legend": lambda ss: _set_gs(ss, "show_legend", True),
    "global show n": lambda ss: _set_gs(ss, "show_n", True),
    "global color scheme": lambda ss: _set_gs(ss, "color_scheme", "simple_white"),
    # gene-prefixed settings that graph.py reads straight out of `settings`
    "gene bar gap": lambda ss: _set_gs(ss, f"{PRIMARY}_bar_gap", 0.05),
    "gene bg color": lambda ss: _set_gs(ss, f"{PRIMARY}_bg_color", "#101010"),
    "gene tick size": lambda ss: _set_gs(ss, f"{PRIMARY}_tick_size", 20),
    "gene ylabel size": lambda ss: _set_gs(ss, f"{PRIMARY}_ylabel_size", 22),
    "gene font size": lambda ss: _set_gs(ss, f"{PRIMARY}_font_size", 9),
    "gene y min": lambda ss: _set_gs(ss, f"{PRIMARY}_y_min", 0.5),
    "gene y max": lambda ss: _set_gs(ss, f"{PRIMARY}_y_max", 12.0),
    "gene label mode": lambda ss: _set_gs(ss, f"{PRIMARY}_label_mode", "Angled 45°"),
    "gene show significance": lambda ss: _set_gs(ss, f"{PRIMARY}_show_sig", False),
    "gene show error": lambda ss: _set_gs(ss, f"{PRIMARY}_show_err", False),
    "gene reference line": lambda ss: _set_gs(ss, f"{PRIMARY}_ref_line", "Treatment1"),
    "gene colour preset": lambda ss: _set_gs(ss, f"{PRIMARY}_color_preset", "Steel"),
    "gene data-point toggle": lambda ss: _set_gs(ss, f"{PRIMARY}_show_data_points", True),
    # per-bar colour override for THIS gene
    "per-bar colour": lambda ss: ss["graph_settings"]["bar_colors_per_sample"].update(
        {f"{PRIMARY}_Treatment1": "#123456"}),
    # per-bar toggles live outside graph_settings entirely
    "per-bar settings": lambda ss: ss[f"{PRIMARY}_bar_settings"].update(
        {f"{PRIMARY}_Treatment1": {"color": "#ABCDEF", "show_sig": False,
                                   "show_sig_1": False, "show_sig_2": False,
                                   "show_sig_3": False, "show_err": False}}),
    # display naming
    "gene display name": lambda ss: ss["gene_display_names"].__setitem__(PRIMARY, "Collagen I"),
    # shared context — several of these are read from session_state INSIDE
    # create_gene_graph and appear nowhere in its argument list
    "sample order": lambda ss: ss.__setitem__("sample_order", ["Treatment1", "Non-treated"]),
    "sample mapping condition": lambda ss: ss["sample_mapping"]["Treatment1"].__setitem__(
        "condition", "Renamed"),
    "sample mapping include flag": lambda ss: ss["sample_mapping"]["Treatment1"].__setitem__(
        "include", False),
    # graph_settings["error_bar_mode"], not a bare session_state key. This
    # mutated ss["error_bar_mode"], which the APP never sets — the mode lives in
    # graph_settings. It passed only because the fingerprint read the same dead
    # key, so test and code agreed on a value neither the widget nor the graph
    # ever used. Both were corrected 2026-08-24.
    "error bar mode": lambda ss: ss["graph_settings"].__setitem__(
        "error_bar_mode", "ci95"),
    "analysis reference condition": lambda ss: ss.__setitem__(
        "analysis_ref_condition", "Treatment1"),
    "comparison condition name": lambda ss: ss.__setitem__(
        "analysis_cmp_condition", "Something Else"),
    "second comparison name": lambda ss: ss.__setitem__(
        "analysis_cmp_condition_2", "Third Arm"),
    "third comparison name": lambda ss: ss.__setitem__(
        "analysis_cmp_condition_3", "Fourth Arm"),
}


@pytest.mark.parametrize("label", sorted(MUTATIONS))
def test_mutation_invalidates_the_cached_figure(app, builds, label):
    app.build_all_figures(GENES, {})
    builds.clear()
    MUTATIONS[label](app.st.session_state)
    app.build_all_figures(GENES, {})
    assert PRIMARY in builds, (
        f"changing {label!r} left {PRIMARY}'s figure cached — the chart would go stale"
    )


def test_global_bar_gap_is_masked_by_resolve_gene_settings(app, builds):
    """Documents a PRE-EXISTING quirk, not a caching bug.

    `resolve_gene_settings` sets ``cs["bar_gap"] = gs.get(f"{gene}_bar_gap", 0.45)``
    — a hardcoded 0.45 default, unlike its sibling lines (font_size, figure_width,
    y_log_scale...) which all chain to the global value. The global "bar_gap" key
    therefore cannot reach a gene chart, so skipping the rebuild is correct.
    If that line is ever fixed to chain to the global, this test should flip to a
    normal invalidation case in MUTATIONS.
    """
    app.build_all_figures(GENES, {})
    builds.clear()
    app.st.session_state["graph_settings"]["bar_gap"] = 0.1
    app.build_all_figures(GENES, {})
    assert builds == []

    # ...whereas the per-gene key that DOES reach the chart must invalidate.
    app.st.session_state["graph_settings"][f"{PRIMARY}_bar_gap"] = 0.1
    app.build_all_figures(GENES, {})
    assert builds == [PRIMARY]


def test_efficacy_config_change_invalidates(app, builds):
    """The config is an argument, not session state, so it needs its own check."""
    app.build_all_figures(GENES, {"expected_direction": {PRIMARY: "up"}})
    builds.clear()
    app.build_all_figures(GENES, {"expected_direction": {PRIMARY: "down"}})
    assert sorted(builds) == sorted(GENES)


# ---- the replicate overlay's extra inputs --------------------------------

def test_raw_data_change_invalidates_when_data_points_shown(app, builds):
    """With the overlay on, the figure depends on the raw CT table too."""
    app.st.session_state["graph_settings"][f"{PRIMARY}_show_data_points"] = True
    app.build_all_figures(GENES, {})
    builds.clear()
    app.st.session_state["data"] = _raw_df(ct=31.0)
    app.build_all_figures(GENES, {})
    assert PRIMARY in builds


def test_qc_exclusion_change_invalidates_when_data_points_shown(app, builds):
    app.st.session_state["graph_settings"][f"{PRIMARY}_show_data_points"] = True
    app.build_all_figures(GENES, {})
    builds.clear()
    app.st.session_state["excluded_wells"] = {PRIMARY: {"A2"}}
    app.build_all_figures(GENES, {})
    assert PRIMARY in builds


def test_hk_gene_change_invalidates_when_data_points_shown(app, builds):
    app.st.session_state["graph_settings"][f"{PRIMARY}_show_data_points"] = True
    app.build_all_figures(GENES, {})
    builds.clear()
    app.st.session_state["hk_gene"] = "ACTB"
    app.build_all_figures(GENES, {})
    assert PRIMARY in builds


def test_exclusion_set_ordering_does_not_cause_spurious_rebuilds(app, builds):
    """Sets have unstable iteration order; an unstable key would miss every rerun."""
    app.st.session_state["graph_settings"][f"{PRIMARY}_show_data_points"] = True
    app.st.session_state["excluded_wells"] = {PRIMARY: {"A2", "A4", "A1"}}
    app.build_all_figures(GENES, {})
    builds.clear()
    # Same contents, different insertion order.
    app.st.session_state["excluded_wells"] = {PRIMARY: {"A4", "A1", "A2"}}
    app.build_all_figures(GENES, {})
    assert builds == []


# ---- safe degradation ----------------------------------------------------

def test_unfingerprintable_data_falls_back_to_rebuilding(app, builds, monkeypatch):
    """If the inputs can't be fingerprinted we must rebuild, never serve stale."""
    monkeypatch.setattr(app, "_df_fingerprint", lambda df: None)
    app.build_all_figures(GENES, {})
    builds.clear()
    app.build_all_figures(GENES, {})
    assert sorted(builds) == sorted(GENES)


def test_failed_build_is_not_cached_and_is_retried(app, monkeypatch):
    """A gene that raises must be retried next rerun, not pinned as 'done'."""
    attempts = []

    def _flaky(gene, gene_data, efficacy_config):
        attempts.append(gene)
        if gene == PRIMARY:
            raise ValueError("render exploded")
        return object()

    monkeypatch.setattr(app, "build_gene_figure", _flaky)
    app.build_all_figures(GENES, {})
    assert PRIMARY not in app.st.session_state["graphs"]
    attempts.clear()
    app.build_all_figures(GENES, {})
    assert PRIMARY in attempts, "failed gene must be retried"


def test_a_gene_with_empty_data_is_skipped(app, builds):
    app.st.session_state["processed_data"][PRIMARY] = pd.DataFrame()
    app.build_all_figures(GENES, {})
    assert builds == [OTHER]
