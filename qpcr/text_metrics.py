"""Text measurement for chart labels — pure, no Plotly and no Streamlit.

Split out of ``qpcr/graph.py`` because it is arithmetic about text, not about
charts, and because ``graph.py`` was past the ~800-line ceiling. Everything here
is a pure function of its arguments, so it is testable without building a figure.

The reason this module exists at all: the label geometry used ``len(label)`` and
a fixed 0.50 em advance. That is a Latin assumption. Hangul, Han and Kana are one
em per codepoint, so a Korean condition name reserved roughly half the space it
needed and the x-axis labels were drawn straight through the caption beneath
them.
"""

import math
import unicodedata

# Mean horizontal advance per character, in em, for the chart font stack. 0.58
# is the value the auto-wrap width has always used; the bottom-margin code used
# 0.50 for the same quantity, so two estimates of one thing disagreed by 16%.
EM_PER_LATIN_CHAR = 0.58
# Hangul, Han, Kana and full-width punctuation are square by definition: one
# codepoint is one em, nearly twice the Latin advance.
EM_PER_WIDE_CHAR = 1.0

_SIN45 = math.sin(math.radians(45.0))
_ANGLED_MODES = ("Angled 45°", "Angled 90°")


def text_em_width(text: str) -> float:
    """Advance width of ``text`` in em, honouring East-Asian full-width glyphs.

    ``len()`` is the wrong unit for a mixed Korean/English/digit label:
    "미처리 대조군" is 7 codepoints but 6.5 em wide, while "abcdefg" is also 7
    codepoints and 4.06 em. Combining marks add no advance of their own.
    """
    total = 0.0
    for ch in str(text):
        if unicodedata.combining(ch):
            continue
        total += (EM_PER_WIDE_CHAR
                  if unicodedata.east_asian_width(ch) in ("W", "F")
                  else EM_PER_LATIN_CHAR)
    return total


def has_wide_chars(text: str) -> bool:
    """True if ``text`` contains any full-width (East-Asian) glyph."""
    return any(unicodedata.east_asian_width(ch) in ("W", "F")
               for ch in str(text))


def effective_wrap_chars(text: str, latin_chars: int) -> int:
    """A character budget for ``text`` costing the same ADVANCE as
    ``latin_chars`` Latin characters.

    ``textwrap`` counts codepoints, which is the wrong unit the moment a label
    is not Latin: a Hangul codepoint is ~1.0 em against Latin's 0.58, so a
    Korean label was allowed ~1.7x the visual width the wrap budget intended and
    ran past its own bar into the neighbouring one. Measured on a default 28cm
    figure with five conditions: the widest Korean line was 247px against 188px
    of bar, 1.32x, on four of the five labels.

    Text with no wide glyph returns ``latin_chars`` UNCHANGED — not merely
    equal by arithmetic, but by an early return — so Latin labels wrap exactly
    where they always did and decks using them cannot move.
    """
    text = str(text)
    if not text or not has_wide_chars(text):
        return latin_chars
    mean_em = text_em_width(text) / len(text)
    if mean_em <= 0:
        return latin_chars
    # Keep at least a few characters: a budget of one produces a column of
    # single glyphs, which is less readable than a slightly wide label.
    return max(int(round(latin_chars * EM_PER_LATIN_CHAR / mean_em)), 4)


def ellipsize_to_em(text: str, max_em: float) -> str:
    """Shorten ``text`` to ``max_em`` em of advance with a MIDDLE ellipsis.

    Middle, not tail: what distinguishes two condition names is usually the
    concentration at the END ("시험물질 1 ppm" vs "시험물질 100 ppm"), so a tail
    ellipsis would render several ticks that read identically.
    """
    text = str(text)
    if max_em <= 0:
        return ""
    if text_em_width(text) <= max_em:
        return text
    ell = "…"
    budget = max_em - text_em_width(ell)
    if budget <= 0 or len(text) <= 2:
        return ell
    head, tail = [], []
    head_em = tail_em = 0.0
    i, j = 0, len(text) - 1
    # Grow from both ends, taking the tail first so the dose survives.
    while i <= j:
        w = text_em_width(text[j])
        if head_em + tail_em + w > budget:
            break
        tail.append(text[j])
        tail_em += w
        j -= 1
        if i > j:
            break
        w = text_em_width(text[i])
        if head_em + tail_em + w > budget:
            break
        head.append(text[i])
        head_em += w
        i += 1
    return "".join(head) + ell + "".join(reversed(tail))


def label_block_px(labels, label_mode: str, tick_px: float,
                   line_px: float) -> float:
    """Height of the drawn tick-label block below the axis, in figure pixels.

    THE single definition. The bottom margin and the caption anchor are both
    derived from one call to this, because computing them separately is exactly
    how a 30-char label at 45° came to reserve 234px of margin while the captions
    stayed anchored a constant 58.9px down and printed through the labels.
    """
    labels = [str(l) for l in labels]
    n_lines = max((l.count("<br>") + 1 for l in labels), default=1)
    if label_mode in _ANGLED_MODES:
        widest_em = max((text_em_width(l) for l in labels),
                        default=5 * EM_PER_LATIN_CHAR)
        text_px = widest_em * tick_px
        if label_mode == "Angled 45°":
            # A rotated line projects BOTH its advance and its line height onto
            # the vertical. The previous formula omitted the second term.
            return (text_px + n_lines * line_px) * _SIN45
        # At 90° extra lines stack sideways, so they add width, not height.
        return text_px
    return n_lines * line_px


def _limit_lines(text: str, max_lines: int) -> str:
    """Keep at most ``max_lines`` wrapped lines, marking anything dropped."""
    lines = str(text).split("<br>")
    if len(lines) <= max_lines:
        return str(text)
    kept = lines[:max_lines]
    kept[-1] = kept[-1].rstrip() + "…"
    return "<br>".join(kept)


def fit_label_block(labels, label_mode: str, tick_px: float, line_px: float,
                    budget_px: float):
    """``(labels, block_px)`` with ``block_px <= budget_px``, by shortening.

    Ellipsising rather than letting it run: an over-long label used to eat the
    plot area until Plotly clipped the text at the image edge. A visible "…"
    tells the operator to widen the figure; silent truncation does not.

    The two modes shorten differently because they grow differently. An angled
    label grows by ADVANCE — one long line projecting further down — so it is
    ellipsized. A wrapped label grows by LINE COUNT, so its lines are capped
    instead. Auto-wrap was left unbounded until the wrap budget started using
    the plot width: narrower lines mean more of them, and a 6cm figure went to
    four lines, 327px of margin on a 452px canvas and 19% of the image left for
    the data.

    Neither path engages until the block exceeds the budget, so ordinary figures
    are untouched.
    """
    labels = [str(l) for l in labels]
    block = label_block_px(labels, label_mode, tick_px, line_px)
    if block <= budget_px or tick_px <= 0:
        return labels, block

    if label_mode not in _ANGLED_MODES:
        if line_px <= 0:
            return labels, block
        max_lines = max(int(budget_px // line_px), 1)
        fitted = [_limit_lines(l, max_lines) for l in labels]
        return fitted, label_block_px(fitted, label_mode, tick_px, line_px)
    # Invert the projection to get the advance budget the labels may use.
    if label_mode == "Angled 45°":
        n_lines = max((l.count("<br>") + 1 for l in labels), default=1)
        max_text_px = (budget_px / _SIN45) - n_lines * line_px
    else:
        max_text_px = budget_px
    max_em = max(max_text_px / tick_px, 1.0)
    fitted = [ellipsize_to_em(l, max_em) for l in labels]
    # Recompute FROM THE SHORTENED LABELS so the margin and the caption anchor
    # can never describe text that is not the text being drawn.
    return fitted, label_block_px(fitted, label_mode, tick_px, line_px)
