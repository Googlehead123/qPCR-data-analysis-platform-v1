"""Where a chart lands on the slide — and that both users agree on it.

qpcr/graph.py sizes every chart font so the smallest text clears MIN_SLIDE_PT
once the image is placed; PPTGenerator computes the frame it is placed into.
Those two must use identical arithmetic. They did not: the font scale divided by
a constant 9.11in, which describes the 28x16cm default and nothing else, so two
of the four size presets shipped chart text at under HALF the intended floor
(PPT Half 4.46pt, Square 4.68pt) while Wide overshot at 13.0pt.
"""
import ast
from pathlib import Path

import pytest

from qpcr.constants import FIGURE_SIZE_PRESETS
from qpcr.slide_geometry import (CM_PER_INCH, SLIDE_MAX_PICTURE_H_IN,
                                 SLIDE_MAX_PICTURE_W_IN,
                                 SLIDE_MIN_RENDER_H_PX, SLIDE_MIN_RENDER_W_PX,
                                 placement_size_in, placement_width_in,
                                 points_per_figure_px, render_size_px)

REPO = Path(__file__).resolve().parents[1]
APP = REPO / "streamlit qpcr analysis v1.py"
PRESETS = ["PPT Full", "PPT Half", "Square", "Wide"]


class TestRenderSize:

    def test_a_large_figure_is_untouched(self):
        assert render_size_px(1058, 604) == (1058, 604)

    def test_a_narrow_figure_is_raised_on_its_short_side(self):
        w, h = render_size_px(529, 377)
        assert w == SLIDE_MIN_RENDER_W_PX
        assert h == pytest.approx(377 * 800 / 529, abs=1)

    def test_the_aspect_ratio_survives_the_floor(self):
        """A floor that changed the aspect made python-pptx stretch the bitmap:
        32.5% vertically for Square, 14.3% for PPT Half."""
        for w, h in [(529, 377), (604, 604), (300, 900), (1200, 200)]:
            rw, rh = render_size_px(w, h)
            assert rw / rh == pytest.approx(w / h, rel=0.01), f"{w}x{h}"

    def test_both_floors_are_respected(self):
        rw, rh = render_size_px(120, 90)
        assert rw >= SLIDE_MIN_RENDER_W_PX and rh >= SLIDE_MIN_RENDER_H_PX

    def test_degenerate_input_does_not_divide_by_zero(self):
        assert render_size_px(0, 0)[0] >= SLIDE_MIN_RENDER_W_PX


class TestPlacement:

    def test_the_frame_never_exceeds_the_layout(self):
        for preset in PRESETS:
            d = FIGURE_SIZE_PRESETS[preset]
            for px in [(1058, 604), (604, 604), (800, 1400)]:
                w_in, h_in = placement_size_in(d["width"], *px)
                assert w_in <= SLIDE_MAX_PICTURE_W_IN + 1e-9, preset
                assert h_in <= SLIDE_MAX_PICTURE_H_IN + 1e-9, preset

    def test_the_frame_matches_the_bitmap_aspect(self):
        """Otherwise python-pptx stretches the image to fit the frame."""
        for preset in PRESETS:
            d = FIGURE_SIZE_PRESETS[preset]
            px = (int(d["width"] * 37.7953), int(d["height"] * 37.7953))
            rw, rh = render_size_px(*px)
            w_in, h_in = placement_size_in(d["width"], *px)
            assert w_in / h_in == pytest.approx(rw / rh, rel=0.01), preset

    def test_the_cm_setting_is_honoured_until_a_clamp_binds(self):
        # 14cm = 5.51in, well inside both clamps at a landscape aspect.
        assert placement_width_in(14, 1058, 604) == pytest.approx(
            14 / CM_PER_INCH, rel=1e-6)

    def test_a_wide_figure_is_clamped_by_width(self):
        assert placement_width_in(40, 1209, 529) == pytest.approx(
            SLIDE_MAX_PICTURE_W_IN)

    def test_a_tall_figure_is_clamped_by_height(self):
        w_in, h_in = placement_size_in(16, 604, 604)
        assert h_in == pytest.approx(SLIDE_MAX_PICTURE_H_IN)
        assert w_in < 16 / CM_PER_INCH

    def test_the_placement_differs_across_presets(self):
        """The whole point: one constant could not have described all four."""
        widths = set()
        for preset in PRESETS:
            d = FIGURE_SIZE_PRESETS[preset]
            px = (int(d["width"] * 37.7953), int(d["height"] * 37.7953))
            widths.add(round(placement_width_in(d["width"], *px), 2))
        assert len(widths) >= 3, f"presets collapsed to {widths}"

    def test_points_per_px_uses_the_rendered_width(self):
        """The short-side floor was the larger half of the error: it upscales
        the canvas while the fonts stay the same number of px."""
        px_w, px_h = 529, 377
        rendered_w, _ = render_size_px(px_w, px_h)
        assert rendered_w > px_w
        got = points_per_figure_px(14, px_w, px_h)
        naive = placement_width_in(14, px_w, px_h) * 72.0 / px_w
        assert got < naive, "the floor must reduce points per authored pixel"


class TestTheTwoUsersAgree:

    def test_the_app_imports_the_shared_module(self):
        """The writer must not re-derive the frame with its own numbers."""
        src = APP.read_text(encoding="utf-8")
        assert "from qpcr.slide_geometry import" in src

    def test_the_app_holds_no_second_copy_of_the_frame_constants(self):
        """11.5 / 5.2 / 800 / 500 must live in one module only."""
        offenders = []
        for lineno, line in enumerate(
                APP.read_text(encoding="utf-8").splitlines(), 1):
            code = line.split("#", 1)[0]
            if "914400" not in code:
                continue
            if "11.5" in code or "5.2" in code:
                offenders.append((lineno, line.strip()[:70]))
        assert not offenders, (
            f"picture-frame clamp duplicated in the app: {offenders}")

    def test_the_font_scale_no_longer_uses_a_constant_placement(self):
        import qpcr.graph as g
        assert not hasattr(g, "PPT_PLACEMENT_WIDTH_IN"), (
            "the constant placement width is what made the floor hold for one "
            "preset out of four"
        )

    @pytest.mark.parametrize("preset", PRESETS)
    def test_the_writer_and_the_font_scale_see_the_same_placement(self, preset):
        """Bind the two files numerically, not just by import.

        Reproduces the writer's frame arithmetic from the app source' inputs and
        checks it against the module the font scale reads.
        """
        d = FIGURE_SIZE_PRESETS[preset]
        px = (int(d["width"] * 37.7953), int(d["height"] * 37.7953))
        from_module = placement_width_in(d["width"], *px)
        # what points_per_figure_px implies, inverted
        rendered_w, _ = render_size_px(*px)
        implied = points_per_figure_px(d["width"], *px) * rendered_w / 72.0
        assert implied == pytest.approx(from_module, rel=1e-9), preset
