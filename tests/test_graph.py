import pandas as pd
import pytest


class TestGraphGeneratorCreateGeneGraph:
    def test_create_gene_graph_returns_figure(self, mock_streamlit, processed_gene_data, graph_settings):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        go = spec.go
        
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data,
            gene='COL1A1',
            settings=graph_settings,
            sample_order=None
        )
        
        assert fig is not None
        assert isinstance(fig, go.Figure)

    def test_create_gene_graph_handles_empty_data(self, mock_streamlit, graph_settings):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        go = spec.go
        
        empty_df = pd.DataFrame()
        
        fig = GraphGenerator.create_gene_graph(
            data=empty_df,
            gene='COL1A1',
            settings=graph_settings,
            sample_order=None
        )
        
        assert fig is not None
        assert isinstance(fig, go.Figure)

    def test_create_gene_graph_handles_none_data(self, mock_streamlit, graph_settings):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        go = spec.go
        
        fig = GraphGenerator.create_gene_graph(
            data=None,
            gene='COL1A1',
            settings=graph_settings,
            sample_order=None
        )
        
        assert fig is not None
        assert isinstance(fig, go.Figure)

    def test_create_gene_graph_uses_fold_change_fallback(self, mock_streamlit, graph_settings):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        go = spec.go
        
        data_without_rel_expr = pd.DataFrame({
            'Target': ['COL1A1', 'COL1A1'],
            'Condition': ['Non-treated', 'Treatment1'],
            'Fold_Change': [1.0, 2.5],
            'SEM': [0.1, 0.2],
            'Group': ['Negative Control', 'Treatment']
        })
        
        fig = GraphGenerator.create_gene_graph(
            data=data_without_rel_expr,
            gene='COL1A1',
            settings=graph_settings,
            sample_order=None
        )
        
        assert fig is not None
        assert len(fig.data) > 0


class TestGraphGeneratorBugFixes:
    def test_significance_symbols_do_not_crash(self, mock_streamlit, graph_settings):
        from qpcr.graph import GraphGenerator
        import plotly.graph_objects as go

        data = pd.DataFrame({
            "Target": ["COL1A1", "COL1A1", "COL1A1"],
            "Condition": ["Non-treated", "Treatment1", "Treatment2"],
            "Group": ["Negative Control", "Treatment", "Treatment"],
            "Relative_Expression": [1.0, 2.5, 0.4],
            "SEM": [0.1, 0.2, 0.05],
            "FC_Error_Upper": [0.15, 0.3, 0.08],
            "FC_Error_Lower": [0.12, 0.25, 0.06],
            "significance": ["", "**", "*"],
            "significance_2": ["", "", ""],
        })
        graph_settings["show_significance"] = True
        graph_settings["show_error"] = True

        fig = GraphGenerator.create_gene_graph(
            data=data, gene="COL1A1", settings=graph_settings
        )
        assert isinstance(fig, go.Figure)
        assert len(fig.layout.annotations) >= 1

    def test_create_gene_graph_respects_sample_order(self, mock_streamlit, processed_gene_data, graph_settings, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        
        mock_streamlit.session_state['sample_mapping'] = sample_mapping
        
        custom_order = ['Treatment2', 'Treatment1', 'Non-treated']
        
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data,
            gene='COL1A1',
            settings=graph_settings,
            sample_order=custom_order
        )
        
        assert fig is not None

    def test_create_gene_graph_includes_error_bars(self, mock_streamlit, processed_gene_data, graph_settings):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        
        graph_settings['show_error'] = True
        
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data,
            gene='COL1A1',
            settings=graph_settings,
            sample_order=None
        )
        
        bar_trace = fig.data[0]
        assert bar_trace.error_y is not None

    def test_create_gene_graph_hides_error_bars_when_disabled(self, mock_streamlit, processed_gene_data, graph_settings):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        
        graph_settings['show_error'] = False
        
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data,
            gene='COL1A1',
            settings=graph_settings,
            sample_order=None
        )
        
        bar_trace = fig.data[0]
        error_values = bar_trace.error_y.array
        assert all(v == 0 for v in error_values)


class TestGraphGeneratorWrapText:
    def test_wrap_text_short_string(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        
        result = GraphGenerator._wrap_text("Short", width=15)
        assert result == "Short"

    def test_wrap_text_long_string(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator
        
        result = GraphGenerator._wrap_text("This is a very long condition name", width=15)
        assert "<br>" in result
