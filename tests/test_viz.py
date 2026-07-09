"""Tests for Phase 4 visualization + reporting rigor.

Plotly figures are introspectable, so error-bar modes, the n= annotations, and
the caption are verified by inspecting the returned go.Figure (no browser).
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from qpcr.analysis import AnalysisEngine
from qpcr.graph import GraphGenerator
from qpcr.constants import GRAPH_PRESETS
from qpcr.auto import build_miqe_checklist


# ---------------- 95% CI columns ----------------

def test_calculate_ddct_emits_ci_columns(mock_streamlit, sample_qpcr_raw_data, sample_mapping):
    res = AnalysisEngine.calculate_ddct(
        data=sample_qpcr_raw_data, hk_gene="GAPDH", ref_sample="Non-treated",
        excluded_wells=set(), excluded_samples=set(), sample_mapping=sample_mapping,
    )
    assert {"FC_CI_Upper", "FC_CI_Lower"} <= set(res.columns)
    spread = res[res["SD"] > 0]
    assert not spread.empty
    assert (spread["FC_CI_Upper"] >= -1e-9).all()
    assert (spread["FC_CI_Lower"] >= -1e-9).all()
    # For n=3, 95% CI half-width = t(.975,2)/sqrt(3) * SD ≈ 2.48*SD > 1*SD,
    # so the CI bounds must exceed the ±SD Livak bounds.
    n3 = spread[spread["n_replicates"] == 3]
    if not n3.empty:
        assert (n3["FC_CI_Upper"] >= n3["FC_Error_Upper"] - 1e-9).all()


# ---------------- graph error-bar modes / n= / caption ----------------

def _gene_df():
    return pd.DataFrame({
        "Condition": ["Ctrl", "Tx"],
        "Relative_Expression": [1.0, 2.0],
        "Fold_Change": [1.0, 2.0],
        "n_replicates": [3, 3],
        "SEM": [0.05, 0.10],
        "SD": [0.09, 0.17],
        "FC_Error_Upper": [0.0, 0.30], "FC_Error_Lower": [0.0, 0.25],
        "FC_CI_Upper": [0.0, 0.70], "FC_CI_Lower": [0.0, 0.55],
        "p_value": [np.nan, 0.01], "significance": ["", "*"],
    })


def _bar_trace(fig):
    return next(t for t in fig.data if isinstance(t, go.Bar))


def test_graph_ci_mode_uses_ci_columns(mock_streamlit):
    fig = GraphGenerator.create_gene_graph(
        _gene_df(), "COL1A1", {"show_error": True, "error_bar_mode": "ci95"},
        color_preset="Steel", ref_condition="Ctrl")
    up = list(_bar_trace(fig).error_y.array)
    assert up[1] == 0.70  # FC_CI_Upper for Tx
    assert any("95% CI" in (a.text or "") for a in fig.layout.annotations)


def test_graph_default_mode_uses_livak_sd(mock_streamlit):
    fig = GraphGenerator.create_gene_graph(
        _gene_df(), "COL1A1", {"show_error": True}, color_preset="Steel", ref_condition="Ctrl")
    up = list(_bar_trace(fig).error_y.array)
    assert up[1] == 0.30  # FC_Error_Upper
    assert any("Livak" in (a.text or "") for a in fig.layout.annotations)


def test_graph_sd_mode_uses_ct_domain(mock_streamlit):
    fig = GraphGenerator.create_gene_graph(
        _gene_df(), "COL1A1", {"show_error": True, "error_bar_mode": "sd"},
        color_preset="Steel", ref_condition="Ctrl")
    up = list(_bar_trace(fig).error_y.array)
    assert up[1] == 0.17  # SD column
    assert any("Ct domain" in (a.text or "") for a in fig.layout.annotations)


def test_graph_n_annotation_optin(mock_streamlit):
    without = GraphGenerator.create_gene_graph(
        _gene_df(), "COL1A1", {"show_error": True}, color_preset="Steel", ref_condition="Ctrl")
    assert not any("n=" in (a.text or "") for a in without.layout.annotations)
    with_n = GraphGenerator.create_gene_graph(
        _gene_df(), "COL1A1", {"show_error": True, "show_n": True},
        color_preset="Steel", ref_condition="Ctrl")
    assert any("n=3" in (a.text or "") for a in with_n.layout.annotations)


# ---------------- palette + MIQE ----------------

def test_colorblind_safe_palette_present():
    assert "Colorblind-Safe" in GRAPH_PRESETS
    assert GRAPH_PRESETS["Colorblind-Safe"]["color"] == "#0072B2"


def test_miqe_checklist_autofills_and_lists_todo(mock_streamlit):
    prov = {
        "generated": "2026-07-09", "software": "qPCR Analysis Suite v3.1",
        "method": "Livak 2^-ddCt, single reference gene", "fdr_correction": "Benjamini-Hochberg",
        "reference_gene": "GAPDH", "reference_condition": "Non-treated",
        "comparison_conditions": ["TGFb"], "statistical_test": "Welch t-test",
        "n_genes": 3, "n_samples": 6, "excluded_wells_count": 2,
    }
    md = build_miqe_checklist(prov)
    assert "MIQE-style reporting checklist" in md
    assert "GAPDH" in md and "Welch t-test" in md
    assert "- [ ]" in md  # unchecked wet-lab items
    assert "Amplification efficiency" in md
