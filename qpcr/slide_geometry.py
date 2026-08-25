"""Where a chart image lands on a slide — the one definition.

Two places have to agree on this and previously did not:

* ``PPTGenerator.create_gene_slide`` computes the picture frame from the
  figure's own aspect ratio, clamps it to the space available on the layout,
  and raises the bitmap's short side before rasterising.
* ``qpcr.graph`` scales every chart font so the smallest text clears
  ``MIN_SLIDE_PT`` **once placed at that width**.

The second used a hardcoded ``PPT_PLACEMENT_WIDTH_IN = 9.11``. That number is
correct for exactly one of the four size presets — the 28x16cm default it was
measured from. Everywhere else it was wrong, and badly:

    preset      placed width   smallest text   vs the 10pt floor
    PPT Full        9.11"          10.5pt      ok
    PPT Half        5.51"           4.5pt      55% under
    Square          5.20"           4.7pt      53% under
    Wide           11.50"          13.0pt      30% over

So the guarantee asserted by the tests held for PPT Full alone, and two presets
shipped chart text at under half the intended minimum.

Nothing here imports plotly or streamlit; it is arithmetic about rectangles.
"""

# The picture frame available on the gene slide, in inches. Measured from the
# template layout: roughly 1.0" to 6.6" vertically, ~12" wide.
SLIDE_MAX_PICTURE_W_IN = 11.5
SLIDE_MAX_PICTURE_H_IN = 5.2

# Short-side floors applied before rasterising. These UPSCALE the canvas while
# every font stays the same number of pixels, so the text gets relatively
# smaller on the bitmap that reaches the slide — which is why the font scale has
# to know about them. They exist so a small figure is not rendered at a size
# where Plotly's own minimums distort the layout.
SLIDE_MIN_RENDER_W_PX = 800
SLIDE_MIN_RENDER_H_PX = 500

EMU_PER_INCH = 914400
CM_PER_INCH = 2.54


def render_size_px(px_w: int, px_h: int):
    """The pixel size the figure is actually rasterised at, aspect preserved.

    Raising the short side rather than the whole frame is deliberate: when a
    floor bound the raw frame instead, the rendered aspect no longer matched the
    frame computed from the cm setting and python-pptx stretched the image to
    fit (32.5% vertically for Square, 14.3% for PPT Half).
    """
    px_w = max(int(px_w or 0), 1)
    px_h = max(int(px_h or 0), 1)
    if px_w < SLIDE_MIN_RENDER_W_PX:
        px_h = int(round(px_h * SLIDE_MIN_RENDER_W_PX / px_w))
        px_w = SLIDE_MIN_RENDER_W_PX
    if px_h < SLIDE_MIN_RENDER_H_PX:
        px_w = int(round(px_w * SLIDE_MIN_RENDER_H_PX / px_h))
        px_h = SLIDE_MIN_RENDER_H_PX
    return px_w, px_h


def placement_size_in(width_cm: float, px_w: int, px_h: int):
    """``(width_in, height_in)`` the picture occupies on the slide.

    The frame follows the RENDERED aspect ratio so the bitmap is never
    stretched; the width still honours the operator's cm setting until one of
    the two clamps binds.
    """
    px_w, px_h = render_size_px(px_w, px_h)
    w_in = max(float(width_cm), 0.01) / CM_PER_INCH
    h_in = w_in * px_h / px_w
    if w_in > SLIDE_MAX_PICTURE_W_IN:
        h_in *= SLIDE_MAX_PICTURE_W_IN / w_in
        w_in = SLIDE_MAX_PICTURE_W_IN
    if h_in > SLIDE_MAX_PICTURE_H_IN:
        w_in *= SLIDE_MAX_PICTURE_H_IN / h_in
        h_in = SLIDE_MAX_PICTURE_H_IN
    return w_in, h_in


def placement_width_in(width_cm: float, px_w: int, px_h: int) -> float:
    """How wide the chart image ends up on the slide, in inches."""
    return placement_size_in(width_cm, px_w, px_h)[0]


def points_per_figure_px(width_cm: float, px_w: int, px_h: int) -> float:
    """Points on the slide per figure pixel.

    This is the conversion the font scale needs: a font set to N figure pixels
    is read at ``N * points_per_figure_px`` on the slide. Note the denominator
    is the RENDERED width, not the authored one — the short-side floors change
    the ratio and were the larger half of the error for PPT Half and Square.
    """
    rendered_w, _ = render_size_px(px_w, px_h)
    return placement_width_in(width_cm, px_w, px_h) * 72.0 / rendered_w
