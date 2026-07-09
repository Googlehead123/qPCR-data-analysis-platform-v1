"""Phase 1 guard: the monolith must USE the qpcr package classes (single source
of truth), not re-inline them, and the unified pipeline must run end-to-end.
"""

from importlib import import_module

import pandas as pd
import plotly.graph_objects as go


def test_monolith_reuses_package_classes(mock_streamlit):
    """Regression guard against re-drift: the app's classes must BE the package
    classes (or a subclass), not independent inline copies."""
    spec = import_module("streamlit qpcr analysis v1")
    import qpcr.parser, qpcr.quality_control, qpcr.graph, qpcr.analysis

    assert spec.QPCRParser is qpcr.parser.QPCRParser
    assert spec.QualityControl is qpcr.quality_control.QualityControl
    assert spec.GraphGenerator is qpcr.graph.GraphGenerator
    # AnalysisEngine is a thin subclass that adds run_full_analysis.
    assert issubclass(spec.AnalysisEngine, qpcr.analysis.AnalysisEngine)
    assert spec.AnalysisEngine.calculate_ddct is qpcr.analysis.AnalysisEngine.calculate_ddct
    assert "run_full_analysis" in spec.AnalysisEngine.__dict__


def test_unified_pipeline_runs_with_ref_condition(mock_streamlit):
    """DDCt -> graph with the ported ref_condition kwarg (the path PPTGenerator uses)."""
    spec = import_module("streamlit qpcr analysis v1")
    rows = []
    for s, cond in [("Ctrl", "Ctrl"), ("Tx", "Treatment")]:
        for i in range(3):
            rows.append({"Well": f"{s}G{i}", "Sample": s, "Target": "GAPDH", "CT": 20.0 + i * 0.1})
            rows.append({"Well": f"{s}C{i}", "Sample": s, "Target": "COL1A1",
                         "CT": (25.0 if s == "Ctrl" else 24.0) + i * 0.1})
    data = pd.DataFrame(rows)
    mapping = {"Ctrl": {"condition": "Ctrl", "group": "Negative Control", "include": True},
               "Tx": {"condition": "Treatment", "group": "Treatment", "include": True}}

    result = spec.AnalysisEngine.calculate_ddct(
        data=data, hk_gene="GAPDH", ref_sample="Ctrl",
        excluded_wells={}, excluded_samples=set(), sample_mapping=mapping,
    )
    assert "FC_Error_Upper" in result.columns  # inherited compute method works
    gene_df = result[result["Target"] == "COL1A1"]

    fig = spec.GraphGenerator.create_gene_graph(
        data=gene_df, gene="COL1A1", settings={"show_error": True},
        color_preset="Steel", ref_condition="Ctrl",
    )
    assert isinstance(fig, go.Figure)
    assert any(isinstance(t, go.Bar) for t in fig.data)
