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
        # Expectation updated 2026-08-24: a switched-off bar is now None, not 0.
        # A 0-length error bar still draws its cap, so "disabled" looked
        # identical to "measured a spread of exactly zero". Plotly omits the
        # bar entirely for None, which is what this test is really asserting.
        assert all(v is None for v in error_values)

    def test_n1_undefined_spread_draws_no_error_bar(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """n=1 has no estimable variance, so it must get NO error bar.

        Regression test for the defect where calculate_ddct correctly stored
        undefined spread as NaN and the graph then coerced it back to 0 with
        .fillna(0) — drawing a zero-length bar, the visual signature of the
        most precise measurement in the panel, on the single well QC left.
        """
        import numpy as np
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        GraphGenerator = spec.GraphGenerator

        data = processed_gene_data.copy()
        # Middle condition is a lone well: NaN spread, as calculate_ddct stores.
        data.loc[1, "n_replicates"] = 1
        for col in ("FC_Error_Upper", "FC_Error_Lower"):
            data[col] = [0.3, np.nan, 0.2]

        graph_settings["show_error"] = True
        fig = GraphGenerator.create_gene_graph(
            data=data, gene='COL1A1', settings=graph_settings, sample_order=None
        )

        # Assert on the SERIALISED figure: that is what Plotly renders, and it
        # is where NaN becomes null (the in-memory array keeps NaN).
        import json
        rendered = json.loads(fig.to_json())["data"][0]["error_y"]["array"]
        assert rendered[1] is None, (
            "the n=1 bar must render no error bar; a 0 would read as a "
            f"measured spread of exactly zero (got {rendered[1]!r})"
        )
        assert rendered[0] is not None and rendered[2] is not None, (
            "suppressing the undefined bar must not suppress its neighbours"
        )
        # And the significance marker must still lay out (it does arithmetic
        # on this array, so a None here used to raise TypeError).
        assert fig.layout is not None


class TestStackedSignificanceSpacing:
    """Stacked significance symbols must clear each other in PIXELS.

    The step was y_max_auto * 0.05 — a fraction of the DATA range, which at the
    default 28x16cm geometry is ~13.5px, smaller than the 16px asterisk it has
    to clear. A dagger's stem was drawn straight through the middle asterisk of
    a "***" and corrupted it on client-facing charts.
    """

    SIG_CHARS = set("*#†")

    def _fig(self, **overrides):
        import pandas as pd
        from qpcr.graph import GraphGenerator

        data = pd.DataFrame({
            "Target": ["KI67"] * 3,
            "Condition": ["Non-treated", "EGF 25 ng/ml", "Test article"],
            "Relative_Expression": [1.0, 0.52, 0.46],
            "Fold_Change": [1.0, 0.52, 0.46],
            "FC_Error_Upper": [0.25, 0.04, 0.01],
            "FC_Error_Lower": [0.25, 0.04, 0.01],
            "n_replicates": [3, 3, 3],
            "significance": ["", "***", "***"],
            "significance_2": ["", "##", "#"],
            "significance_3": ["", "†", "††"],
        })
        settings = {"show_error": True, "show_significance": True,
                    "figure_height": 16, "figure_width": 28}
        settings.update(overrides)
        return GraphGenerator.create_gene_graph(
            data=data, gene="KI67", settings=settings, sample_order=None
        )

    def _stacks(self, fig):
        from qpcr.constants import CM_TO_PX
        m = fig.layout.margin
        plot_h = int(16 * CM_TO_PX) - (m.b or 0) - (m.t or 0)
        lo, hi = fig.layout.yaxis.range
        per_px = (hi - lo) / plot_h
        by_x = {}
        for a in fig.layout.annotations:
            if a.text and set(a.text) <= self.SIG_CHARS and a.yref == "y":
                by_x.setdefault(a.x, []).append((a.y, a.font.size, a.text))
        for items in by_x.values():
            items.sort()
        return by_x, per_px, hi

    def test_stacked_symbols_do_not_overlap(self, mock_streamlit):
        fig = self._fig()
        by_x, per_px, _ = self._stacks(fig)
        assert by_x, "expected stacked significance symbols on this fixture"
        for x, items in by_x.items():
            for i in range(1, len(items)):
                y, fs, txt = items[i]
                prev_y, prev_fs, prev_txt = items[i - 1]
                gap_px = (y - prev_y) / per_px
                assert gap_px >= max(prev_fs, fs), (
                    f"bar {x}: {prev_txt!r} -> {txt!r} step is {gap_px:.1f}px "
                    f"but must clear a {max(prev_fs, fs)}px glyph"
                )

    def test_stack_is_not_clipped_by_the_axis(self, mock_streamlit):
        """Widening the step must also widen the range that has to contain it."""
        fig = self._fig()
        by_x, per_px, y_hi = self._stacks(fig)
        for x, items in by_x.items():
            top_y, top_fs, _ = items[-1]
            assert top_y + top_fs * 1.2 * per_px <= y_hi + 1e-9, (
                f"bar {x}: the top symbol runs past the axis maximum {y_hi}"
            )

    def test_a_short_figure_still_clears_the_glyphs(self, mock_streamlit):
        """The step is geometry-dependent, so check a much shorter figure too."""
        fig = self._fig(figure_height=8)
        from qpcr.constants import CM_TO_PX
        m = fig.layout.margin
        plot_h = max(int(8 * CM_TO_PX) - (m.b or 0) - (m.t or 0), 100)
        lo, hi = fig.layout.yaxis.range
        per_px = (hi - lo) / plot_h
        by_x = {}
        for a in fig.layout.annotations:
            if a.text and set(a.text) <= self.SIG_CHARS and a.yref == "y":
                by_x.setdefault(a.x, []).append((a.y, a.font.size, a.text))
        for items in by_x.values():
            items.sort()
            for i in range(1, len(items)):
                gap_px = (items[i][0] - items[i - 1][0]) / per_px
                assert gap_px >= max(items[i - 1][1], items[i][1])


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


class TestColorPresets:
    def test_all_presets_have_color_and_ref(self):
        from qpcr.constants import GRAPH_PRESETS
        for preset_name, preset_def in GRAPH_PRESETS.items():
            assert "color" in preset_def, f"Preset '{preset_name}' missing 'color'"
            assert "ref" in preset_def, f"Preset '{preset_name}' missing 'ref'"

    def test_all_preset_colors_are_valid_hex(self):
        import re
        from qpcr.constants import GRAPH_PRESETS
        hex_re = re.compile(r'^#[0-9A-Fa-f]{6}$')
        for preset_name, preset_def in GRAPH_PRESETS.items():
            assert hex_re.match(preset_def["color"]), f"Invalid hex in {preset_name}/color: {preset_def['color']}"
            assert hex_re.match(preset_def["ref"]), f"Invalid hex in {preset_name}/ref: {preset_def['ref']}"

    def test_all_preset_refs_are_white(self):
        from qpcr.constants import GRAPH_PRESETS
        for preset_name, preset_def in GRAPH_PRESETS.items():
            assert preset_def["ref"] == "#FFFFFF", f"Preset '{preset_name}' ref should be white"

    def test_default_group_colors_all_white(self):
        from qpcr.constants import DEFAULT_GROUP_COLORS
        for group, color in DEFAULT_GROUP_COLORS.items():
            assert color == "#FFFFFF", f"Group '{group}' should be white, got {color}"

    def test_figure_size_presets_have_width_and_height(self):
        from qpcr.constants import FIGURE_SIZE_PRESETS
        for name, dims in FIGURE_SIZE_PRESETS.items():
            assert "width" in dims and "height" in dims, f"Preset '{name}' missing width/height"


class TestVisualPolish:
    def test_axis_color_is_dark(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings
        )
        assert fig.layout.yaxis.linecolor == "#2C3E50"

    def test_axis_line_width(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings
        )
        assert fig.layout.yaxis.linewidth == 1.5

    def test_default_bar_opacity(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        graph_settings.pop("bar_opacity", None)
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings
        )
        assert fig.data[0].marker.opacity == 0.85

    def test_error_bar_cap_width(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        graph_settings["show_error"] = True
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings
        )
        assert fig.data[0].error_y.width == 6


class TestDataPointOverlay:
    def test_no_scatter_trace_when_disabled(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings,
            show_data_points=False,
        )
        assert len(fig.data) == 1  # Only bar trace

    def test_scatter_trace_when_enabled(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        import pandas as pd
        replicate_data = pd.DataFrame({
            "Target": ["COL1A1"] * 9,
            "Condition": ["Non-treated"] * 3 + ["Treatment1"] * 3 + ["Treatment2"] * 3,
            "Well": [f"A{i}" for i in range(1, 10)],
            "Replicate_FC": [0.95, 1.02, 1.03, 2.4, 2.6, 2.55, 0.35, 0.38, 0.37],
        })
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings,
            show_data_points=True, replicate_data=replicate_data,
        )
        assert len(fig.data) == 2  # Bar + scatter
        assert fig.data[1].mode == "markers"

    def test_no_scatter_when_enabled_but_no_replicate_data(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings,
            show_data_points=True, replicate_data=None,
        )
        assert len(fig.data) == 1  # Only bar, no scatter


class TestSignificanceBrackets:
    def test_direct_mode_uses_annotations(self, mock_streamlit, graph_settings):
        from qpcr.graph import GraphGenerator
        import plotly.graph_objects as go
        data = pd.DataFrame({
            "Target": ["COL1A1", "COL1A1"],
            "Condition": ["Non-treated", "Treatment1"],
            "Group": ["Negative Control", "Treatment"],
            "Relative_Expression": [1.0, 2.5],
            "SEM": [0.1, 0.2],
            "FC_Error_Upper": [0.15, 0.3],
            "FC_Error_Lower": [0.12, 0.25],
            "significance": ["", "**"],
            "significance_2": ["", ""],
        })
        graph_settings["show_significance"] = True
        fig = GraphGenerator.create_gene_graph(data=data, gene="COL1A1", settings=graph_settings)
        sig_annotations = [a for a in fig.layout.annotations if a.text in ["*", "**", "***"]]
        assert len(sig_annotations) >= 1

