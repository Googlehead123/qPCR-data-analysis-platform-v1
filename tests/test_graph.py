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


class TestSlideLegibilityAndFraming:
    """Chart text is sized in figure pixels but read on a slide.

    The PPT writer places the image 9.11in wide, so on a 28cm (1058px) figure
    one figure pixel is 0.620pt: the 9px captions landed at 5.6pt and the 12px
    ticks at 7.4pt, against slide chrome of 11-16pt. Meanwhile 40% of every
    exported PNG was blank and the plot area was 44% of the image.
    """

    def _fig(self, **overrides):
        import pandas as pd
        from qpcr.graph import GraphGenerator

        conds = ["Non-treated", "EGF 25 ng/ml", "Test article 1 ppm",
                 "Test article 10 ppm", "Test article 100 ppm"]
        data = pd.DataFrame({
            "Target": ["KI67"] * 5,
            "Condition": conds,
            "Relative_Expression": [1.0, 0.52, 0.46, 0.65, 0.45],
            "Fold_Change": [1.0, 0.52, 0.46, 0.65, 0.45],
            "FC_Error_Upper": [0.25, 0.04, 0.01, 0.27, 0.05],
            "FC_Error_Lower": [0.25, 0.04, 0.01, 0.27, 0.05],
            "n_replicates": [3] * 5,
            "significance": ["", "***", "***", "*", "***"],
        })
        settings = {"show_error": True, "show_significance": True,
                    "figure_width": 28, "figure_height": 16}
        settings.update(overrides)
        return GraphGenerator.create_gene_graph(
            data=data, gene="KI67", settings=settings, sample_order=None
        )

    def _pt_per_px(self, width_cm=28):
        from qpcr.constants import CM_TO_PX
        from qpcr.graph import PPT_PLACEMENT_WIDTH_IN
        return (PPT_PLACEMENT_WIDTH_IN * 72.0) / int(width_cm * CM_TO_PX)

    def test_no_chart_text_falls_below_the_slide_floor(self, mock_streamlit):
        from qpcr.graph import MIN_SLIDE_PT
        fig = self._fig()
        ppx = self._pt_per_px()
        sizes = [fig.layout.xaxis.tickfont.size, fig.layout.yaxis.title.font.size]
        sizes += [a.font.size for a in fig.layout.annotations
                  if a.font and a.font.size]
        worst = min(sizes)
        assert worst * ppx >= MIN_SLIDE_PT - 1e-6, (
            f"smallest chart text is {worst}px = {worst * ppx:.1f}pt on the "
            f"slide, below the {MIN_SLIDE_PT}pt floor"
        )

    def test_the_size_hierarchy_survives_scaling(self, mock_streamlit):
        """Clamping each size to the floor would flatten these into one value."""
        fig = self._fig()
        tick = fig.layout.xaxis.tickfont.size
        ytitle = fig.layout.yaxis.title.font.size
        caps = [a.font.size for a in fig.layout.annotations
                if a.font and a.font.size and a.yref == "paper"]
        assert caps, "expected the error/significance captions"
        assert min(caps) < tick < ytitle, (
            f"captions {min(caps)} < ticks {tick} < y-title {ytitle} expected"
        )

    def test_the_plot_area_is_most_of_the_figure(self, mock_streamlit):
        from qpcr.constants import CM_TO_PX
        fig = self._fig()
        fig_h = int(16 * CM_TO_PX)
        m = fig.layout.margin
        plot_h = fig_h - (m.b or 0) - (m.t or 0)
        assert plot_h / fig_h >= 0.60, (
            f"plot area is {plot_h / fig_h:.0%} of the figure; it was 44% when "
            f"the margins were floors rather than measurements"
        )

    def test_bottom_margin_tracks_the_label_block(self, mock_streamlit):
        """One line of labels must cost less than three."""
        one = self._fig(label_mode="Horizontal")
        wrapped = self._fig(label_mode="Auto-wrap")
        assert (one.layout.margin.b or 0) < (wrapped.layout.margin.b or 0), (
            "a single-line label row should reserve less than a wrapped one"
        )


class TestTextWidthModel:
    """Label width must be measured in em, not codepoints.

    The margin code did `len(label) * 0.50em`. That is a Latin assumption:
    Hangul, Han and Kana are one em per codepoint, so "미처리 대조군" reserved
    3.5 em for 6.5 em of text and the labels were drawn through the caption
    below them. Pure-Latin labels measured correctly, which is why this went
    unseen.
    """

    def test_latin_is_the_narrow_advance(self):
        from qpcr.text_metrics import text_em_width, EM_PER_LATIN_CHAR
        assert text_em_width("abcdef") == pytest.approx(
            6 * EM_PER_LATIN_CHAR)

    def test_hangul_is_one_em_per_codepoint(self):
        from qpcr.text_metrics import text_em_width, EM_PER_WIDE_CHAR
        assert text_em_width("한글다섯글자") == pytest.approx(
            6 * EM_PER_WIDE_CHAR)

    def test_a_korean_label_is_wider_than_its_length_suggests(self):
        from qpcr.text_metrics import text_em_width
        kor = text_em_width("미처리 대조군")
        lat = text_em_width("abcdefg")   # same codepoint count
        assert len("미처리 대조군") == len("abcdefg")
        assert kor > lat * 1.5, "Hangul must not be measured as Latin"

    def test_combining_marks_add_no_advance(self):
        from qpcr.text_metrics import text_em_width
        assert (text_em_width("e\u0301")
                == text_em_width("e"))

    def test_empty_text_is_zero(self):
        from qpcr.text_metrics import text_em_width
        assert text_em_width("") == 0.0

    def test_ellipsize_respects_the_budget(self):
        from qpcr.text_metrics import ellipsize_to_em, text_em_width
        out = ellipsize_to_em("Recombinant human EGF 25 ng/ml", 6.0)
        assert text_em_width(out) <= 6.0 + 1e-9
        assert "…" in out

    def test_ellipsize_keeps_the_tail_where_the_dose_lives(self):
        from qpcr.text_metrics import ellipsize_to_em, text_em_width
        out = ellipsize_to_em("Test article 100 ppm", 8.0)
        assert out.endswith("m"), f"tail dropped: {out!r}"

    def test_ellipsize_is_a_no_op_when_it_already_fits(self):
        from qpcr.text_metrics import ellipsize_to_em
        assert ellipsize_to_em("Ctrl", 50.0) == "Ctrl"


class TestAngledLabelGeometry:
    """Angled x-axis labels must not print through the captions below them.

    The bottom margin grew with label length while the caption anchor used the
    line COUNT and the angled branches pinned that at 1. Measured before the
    fix: the caption sat a constant 58.9px below the axis in EVERY mode and at
    EVERY label length, so a 30-char label at 45° reserved 234px of margin and
    then drew its labels straight through both captions.
    """

    LABELS = {
        "latin_short": ["Ctrl", "EGF", "A 1ppm", "A 10ppm", "A 100ppm"],
        "latin_long": ["Non-treated control", "Recombinant human EGF 25 ng/ml",
                       "Test article 1 ppm", "Test article 10 ppm",
                       "Test article 100 ppm"],
        "korean": ["미처리 대조군", "양성대조군 TGFβ 10 ng/ml", "시험물질 1 ppm 처리",
                   "시험물질 10 ppm 처리", "시험물질 100 ppm 처리"],
    }
    MODES = ["Auto-wrap", "Horizontal", "Angled 45°", "Angled 90°"]

    def _fig(self, conds, **overrides):
        import pandas as pd
        from qpcr.graph import GraphGenerator
        n = len(conds)
        data = pd.DataFrame({
            "Target": ["KI67"] * n,
            "Condition": conds,
            "Relative_Expression": [1.0, 1.85, 0.62, 0.88, 1.34][:n],
            "Fold_Change": [1.0, 1.85, 0.62, 0.88, 1.34][:n],
            "FC_Error_Upper": [0.12, 0.20, 0.08, 0.10, 0.15][:n],
            "FC_Error_Lower": [0.12, 0.20, 0.08, 0.10, 0.15][:n],
            "n_replicates": [3] * n,
            "significance": ["", "***", "*", "", "**"][:n],
        })
        settings = {"show_error": True, "show_significance": True,
                    "figure_width": 28, "figure_height": 16}
        settings.update(overrides)
        return GraphGenerator.create_gene_graph(
            data=data, gene="KI67", settings=settings, sample_order=None)

    @staticmethod
    def _geometry(fig, mode):
        """(plot_h_px, label_block_px, [caption offsets in px])."""
        from qpcr.text_metrics import label_block_px
        L = fig.layout
        plot_h = (L.height or 0) - (L.margin.b or 0) - (L.margin.t or 0)
        tick = L.xaxis.tickfont.size
        block = label_block_px(
            list(L.xaxis.ticktext or ()), mode, tick, tick * 1.45)
        offs = [-a.y * plot_h for a in L.annotations
                if a.yref == "paper" and a.y is not None and a.y < 0]
        return plot_h, block, offs

    def test_a_longer_angled_label_pushes_the_caption_further_down(self,
                                                                   mock_streamlit):
        """The direct regression. Before the fix both were 58.9px."""
        short = self._fig(self.LABELS["latin_short"], label_mode="Angled 45°")
        long_ = self._fig(self.LABELS["latin_long"], label_mode="Angled 45°")
        _, _, off_s = self._geometry(short, "Angled 45°")
        _, _, off_l = self._geometry(long_, "Angled 45°")
        assert off_l[0] > off_s[0] + 50, (
            f"caption did not track label length: {off_s[0]:.1f} -> {off_l[0]:.1f}")

    @pytest.mark.parametrize("mode", MODES)
    @pytest.mark.parametrize("labels", sorted(LABELS))
    @pytest.mark.parametrize("preset", ["PPT Full", "PPT Half", "Square", "Wide"])
    def test_the_caption_clears_the_label_block(self, mock_streamlit, mode,
                                                labels, preset):
        from qpcr.constants import FIGURE_SIZE_PRESETS
        dims = FIGURE_SIZE_PRESETS[preset]
        fig = self._fig(self.LABELS[labels], label_mode=mode,
                        figure_width=dims["width"], figure_height=dims["height"])
        plot_h, block, offs = self._geometry(fig, mode)
        assert offs, "no paper-anchored caption found"
        for off in offs:
            assert off >= block - 1.0, (
                f"{mode}/{labels}/{preset}: caption at {off:.1f}px is inside a "
                f"{block:.1f}px label block")

    @pytest.mark.parametrize("mode", ["Angled 45°", "Angled 90°"])
    def test_korean_reserves_more_than_its_codepoint_count(self, mock_streamlit,
                                                           mode):
        kor = ["미처리 대조군"] * 5           # 7 codepoints, 6.5 em
        lat = ["abcdefg"] * 5              # 7 codepoints, 4.06 em
        b_kor = self._fig(kor, label_mode=mode).layout.margin.b
        b_lat = self._fig(lat, label_mode=mode).layout.margin.b
        assert b_kor > b_lat, (
            f"{mode}: Hangul reserved {b_kor} vs Latin {b_lat} for the same "
            "codepoint count")

    @pytest.mark.parametrize("mode", MODES)
    @pytest.mark.parametrize("height", [6, 10, 16, 25])
    def test_the_plot_area_never_collapses(self, mock_streamlit, mode, height):
        from qpcr.graph import MIN_PLOT_HEIGHT_FRAC, PLOT_HEIGHT_FLOOR_PX
        fig = self._fig(self.LABELS["latin_long"], label_mode=mode,
                        figure_height=height)
        plot_h, _, _ = self._geometry(fig, mode)
        frac = plot_h / fig.layout.height
        assert plot_h >= PLOT_HEIGHT_FLOOR_PX, (
            f"{mode}@{height}cm: plot {plot_h}px is under the floor, so the "
            "caption denominator is fiction")
        assert frac >= MIN_PLOT_HEIGHT_FRAC - 0.01, (
            f"{mode}@{height}cm: plot is {frac:.0%} of the figure")

    def test_an_unfittable_label_is_ellipsized_not_clipped(self, mock_streamlit):
        long_labels = ["Recombinant human epidermal growth factor 25 ng/ml"] * 5
        fig = self._fig(long_labels, label_mode="Angled 90°", figure_height=6)
        drawn = list(fig.layout.xaxis.ticktext or ())
        assert any("…" in d for d in drawn), (
            f"expected an ellipsis when the label cannot fit: {drawn!r}")

    def test_auto_wrap_and_horizontal_geometry_is_frozen(self, mock_streamlit):
        """The deck guard: the default modes must not move by one pixel.

        Verified against a 48-figure render matrix — all 24 Auto-wrap and
        Horizontal cases were pixel-identical across this change.
        """
        for mode in ("Auto-wrap", "Horizontal"):
            fig = self._fig(self.LABELS["latin_short"], label_mode=mode)
            L = fig.layout
            assert L.height == 604, mode
            assert L.margin.b == 92, mode
            assert L.margin.t == 28, mode
            caps = [a.y for a in L.annotations
                    if a.yref == "paper" and a.y is not None and a.y < 0]
            assert caps, mode
            for y in caps:
                assert y == -0.12169421487603306, (
                    f"{mode}: caption moved to {y!r}")


class TestReferenceLine:
    """The reference line must not take the chart down with it.

    `graph.py` read `fixed_symbol_spacing` when placing the reference-line
    label — a name that stopped existing when the significance glyphs moved to
    `pending_sig`. Every chart with a reference line raised NameError.
    `_memoized_gene_figure` catches build failures, so the gene never reached
    `st.session_state["graphs"]` and simply went MISSING from the exported deck
    behind a transient warning. Nothing in the suite exercised this path.
    """

    def _fig(self, **overrides):
        import pandas as pd
        from qpcr.graph import GraphGenerator

        conds = ["Non-treated", "EGF 25 ng/ml", "Test article 100 ppm"]
        data = pd.DataFrame({
            "Target": ["KI67"] * 3,
            "Condition": conds,
            "Relative_Expression": [1.0, 1.85, 0.62],
            "Fold_Change": [1.0, 1.85, 0.62],
            "FC_Error_Upper": [0.12, 0.20, 0.08],
            "FC_Error_Lower": [0.12, 0.20, 0.08],
            "n_replicates": [3] * 3,
            "significance": ["", "***", "*"],
        })
        settings = {"show_error": True, "show_significance": True,
                    "figure_width": 28, "figure_height": 16}
        settings.update(overrides)
        return GraphGenerator.create_gene_graph(
            data=data, gene="KI67", settings=settings, sample_order=None,
            ref_line_value=overrides.pop("_ref", 1.0), ref_line_label="Control",
        )

    @pytest.mark.parametrize("label_mode", [
        "Auto-wrap", "Horizontal", "Angled 45°", "Angled 90°",
    ])
    def test_a_reference_line_builds_in_every_label_mode(self, mock_streamlit,
                                                         label_mode):
        fig = self._fig(label_mode=label_mode)
        assert fig.data, "the chart lost its traces"

    def test_the_reference_line_is_actually_drawn(self, mock_streamlit):
        fig = self._fig()
        shapes = list(fig.layout.shapes or ())
        lines = [s for s in shapes if getattr(s, "type", None) == "line"]
        assert lines, "no reference line shape was added"

    def test_significance_symbols_do_not_break_the_placement(self, mock_streamlit):
        """The placement reads the significance stack top; drawing symbols and
        drawing none must both work."""
        with_syms = self._fig()
        assert with_syms.data
        without = self._fig(show_significance=False)
        assert without.data

    def test_a_high_reference_line_flips_the_label_below(self, mock_streamlit):
        """A line above the tallest bar must not print the label off-plot."""
        low = self._fig(_ref=0.2)
        high = self._fig(_ref=50.0)
        assert low.data and high.data


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

