"""Guards for the browser-geometry copy of a chart figure.

One figure per gene is memoized and it is the EXPORT figure: 1058x604px with
every font multiplied so the smallest text clears 10pt on a slide. Streamlit's
frontend replaces layout.width with the container width but keeps the height,
the fonts and the margins, so the editor drew a 1.75:1 chart at 1.16:1 and the
Overview panel at 0.58:1 — portrait — with the fixed 80px left margin taking
23% of the element.

The display path now scales a COPY uniformly. Uniform is the safety property:
both captions are anchored as a fraction of plot height and the significance
stack is spaced in data units derived from it, so one factor applied to height
and margins leaves those relationships describing the same place.

The test that matters most here is that the memoized export figure is never
mutated — every export path reads it directly.
"""
import ast
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
APP = REPO / "streamlit qpcr analysis v1.py"


@pytest.fixture
def export_fig(mock_streamlit):
    from qpcr.graph import GraphGenerator
    conds = ["Non-treated", "EGF 25 ng/ml", "시험물질 100 ppm 처리",
             "Test article 10 ppm", "Test article 100 ppm"]
    data = pd.DataFrame({
        "Target": ["KI67"] * 5,
        "Condition": conds,
        "Relative_Expression": [1.0, 1.85, 0.62, 0.88, 1.34],
        "Fold_Change": [1.0, 1.85, 0.62, 0.88, 1.34],
        "FC_Error_Upper": [0.12, 0.20, 0.08, 0.10, 0.15],
        "FC_Error_Lower": [0.12, 0.20, 0.08, 0.10, 0.15],
        "n_replicates": [3] * 5,
        "significance": ["", "***", "*", "", "**"],
    })
    return GraphGenerator.create_gene_graph(
        data=data, gene="KI67",
        settings={"show_error": True, "show_significance": True,
                  "figure_width": 28, "figure_height": 16},
        sample_order=None,
    )


def _app(mock_streamlit):
    from importlib import import_module
    return import_module("streamlit qpcr analysis v1")


class TestScreenFigureScaling:

    def test_the_export_figure_is_never_mutated(self, mock_streamlit, export_fig):
        """The one that matters: every export path reads this object."""
        spec = _app(mock_streamlit)
        before = export_fig.to_json()
        spec._scale_figure_for_screen(export_fig, 400)
        assert export_fig.to_json() == before, (
            "the memoized export figure was mutated; the deck would follow the "
            "browser's geometry"
        )

    def test_the_width_is_released_to_the_container(self, mock_streamlit,
                                                    export_fig):
        spec = _app(mock_streamlit)
        out = spec._scale_figure_for_screen(export_fig, 400)
        assert out.layout.width is None, (
            "a fixed width can overflow the element and be clipped by "
            "overflow:hidden"
        )

    def test_the_aspect_ratio_is_preserved(self, mock_streamlit, export_fig):
        spec = _app(mock_streamlit)
        src_w, src_h = export_fig.layout.width, export_fig.layout.height
        out = spec._scale_figure_for_screen(export_fig, 400)
        assert out.layout.height == pytest.approx(src_h * 400 / src_w, abs=1.5)

    def test_every_font_is_scaled_by_the_same_factor(self, mock_streamlit,
                                                     export_fig):
        spec = _app(mock_streamlit)
        target = 400
        s = target / export_fig.layout.width
        out = spec._scale_figure_for_screen(export_fig, target)

        pairs = [(export_fig.layout.font.size, out.layout.font.size),
                 (export_fig.layout.xaxis.tickfont.size,
                  out.layout.xaxis.tickfont.size)]
        for a, b in zip(export_fig.layout.annotations, out.layout.annotations):
            if a.font and a.font.size:
                pairs.append((a.font.size, b.font.size))
        for src, got in pairs:
            assert got == max(1, round(src * s)), (
                f"font {src} -> {got}, expected {max(1, round(src * s))}")

    def test_margins_scale_with_the_fonts(self, mock_streamlit, export_fig):
        spec = _app(mock_streamlit)
        target = 400
        s = target / export_fig.layout.width
        out = spec._scale_figure_for_screen(export_fig, target)
        for side in ("l", "r", "t", "b"):
            src = getattr(export_fig.layout.margin, side)
            got = getattr(out.layout.margin, side)
            assert got == max(1, round(src * s)), f"margin.{side}"

    def test_paper_anchored_captions_keep_their_fraction(self, mock_streamlit,
                                                         export_fig):
        """The invariant uniform scaling buys: no annotation is rewritten."""
        spec = _app(mock_streamlit)
        out = spec._scale_figure_for_screen(export_fig, 400)
        src = [a.y for a in export_fig.layout.annotations
               if a.yref == "paper" and a.y is not None]
        got = [a.y for a in out.layout.annotations
               if a.yref == "paper" and a.y is not None]
        assert src and src == got

    def test_the_left_margin_is_a_sane_share_of_a_thumbnail(self, mock_streamlit,
                                                            export_fig):
        """Was 80px of a ~350px element — 23% before a single bar was drawn."""
        spec = _app(mock_streamlit)
        out = spec._scale_figure_for_screen(export_fig, 400)
        assert out.layout.margin.l / 400 < 0.15

    def test_a_figure_narrower_than_the_target_is_returned_unchanged(
            self, mock_streamlit, export_fig):
        spec = _app(mock_streamlit)
        assert spec._scale_figure_for_screen(export_fig, 5000) is export_fig

    def test_a_figure_without_a_width_is_returned_unchanged(self, mock_streamlit):
        import plotly.graph_objects as go
        spec = _app(mock_streamlit)
        bare = go.Figure()
        assert spec._scale_figure_for_screen(bare, 400) is bare

    def test_the_editor_is_wider_than_the_overview_thumbnail(self, mock_streamlit):
        spec = _app(mock_streamlit)
        assert (spec.SCREEN_CHART_WIDTH_PX["editor"]
                > spec.SCREEN_CHART_WIDTH_PX["overview"])


class TestScreenFigureCache:

    def test_the_cache_is_keyed_on_the_export_signature(self, mock_streamlit,
                                                        export_fig):
        spec = _app(mock_streamlit)
        mock_streamlit.session_state["_graph_signatures"] = {"KI67": "sig-1"}
        mock_streamlit.session_state["_screen_graphs"] = {}
        first = spec._screen_figure("KI67", export_fig, "editor")
        again = spec._screen_figure("KI67", export_fig, "editor")
        assert again is first, "a signature hit should not rebuild"

        mock_streamlit.session_state["_graph_signatures"]["KI67"] = "sig-2"
        after = spec._screen_figure("KI67", export_fig, "editor")
        assert after is not first, "a changed signature must invalidate"

    def test_editor_and_overview_do_not_share_an_entry(self, mock_streamlit,
                                                       export_fig):
        spec = _app(mock_streamlit)
        mock_streamlit.session_state["_graph_signatures"] = {"KI67": "sig-1"}
        mock_streamlit.session_state["_screen_graphs"] = {}
        ed = spec._screen_figure("KI67", export_fig, "editor")
        ov = spec._screen_figure("KI67", export_fig, "overview")
        assert ed.layout.height != ov.layout.height

    def test_an_unfingerprintable_figure_is_always_copied(self, mock_streamlit,
                                                          export_fig):
        spec = _app(mock_streamlit)
        mock_streamlit.session_state["_graph_signatures"] = {}
        mock_streamlit.session_state["_screen_graphs"] = {}
        a = spec._screen_figure("KI67", export_fig, "editor")
        b = spec._screen_figure("KI67", export_fig, "editor")
        assert a is not b


class TestEveryChartGoesThroughTheScaler:
    """A third chart added later must not bypass the screen geometry.

    streamlit.testing.v1 has no plotly element accessor, so the shipped figure
    cannot be read back from an AppTest. Parsing the source is the practical
    substitute.
    """

    def test_every_plotly_chart_call_scales_first(self):
        tree = ast.parse(APP.read_text(encoding="utf-8"))
        bad = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            fn = node.func
            if not (isinstance(fn, ast.Attribute) and fn.attr == "plotly_chart"):
                continue
            if not node.args:
                bad.append((node.lineno, "no positional figure"))
                continue
            arg = node.args[0]
            ok = (isinstance(arg, ast.Call)
                  and isinstance(arg.func, ast.Name)
                  and arg.func.id == "_screen_figure")
            if not ok:
                bad.append((node.lineno, ast.dump(arg)[:60]))
        assert not bad, (
            "st.plotly_chart must be handed a _screen_figure(...) result; "
            f"offenders: {bad}"
        )
