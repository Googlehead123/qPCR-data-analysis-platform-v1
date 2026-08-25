"""Utility functions for qPCR analysis.

Contains sorting helpers, well exclusion key generation, QC grid state
management, and the RdYlGn table gradient.
"""

import re
from typing import Optional


# ==================== TABLE GRADIENT ====================
# pandas' Styler.background_gradient(cmap=...) imports matplotlib at CALL time.
# Nothing in this project writes `import matplotlib`, so an audit that greps for
# imports concludes it is unused — and dropping the pin then breaks every styled
# results table with "background_gradient requires matplotlib", but only on a
# clean install. A dev machine that still has matplotlib from an older
# requirements.txt stays green, which is exactly how this reached CI.
#
# Reproducing the ramp here keeps that dependency gone (7 packages, ~15MB off
# every cold boot) and keeps the "Plotly only" rule true for charts.

# ColorBrewer RdYlGn, 11-class — the anchors matplotlib's "RdYlGn" interpolates.
_RDYLGN = (
    (165, 0, 38), (215, 48, 39), (244, 109, 67), (253, 174, 97),
    (254, 224, 139), (255, 255, 191), (217, 239, 139), (166, 217, 106),
    (102, 189, 99), (26, 152, 80), (0, 104, 55),
)

# pandas' own default: below this relative luminance it flips the text to light.
_TEXT_COLOR_THRESHOLD = 0.408


def rdylgn_hex(t: float) -> str:
    """Hex colour at position ``t`` (clamped to [0, 1]) on the RdYlGn ramp."""
    t = 0.0 if t < 0.0 else (1.0 if t > 1.0 else float(t))
    pos = t * (len(_RDYLGN) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(_RDYLGN) - 1)
    frac = pos - lo
    r, g, b = (
        round(_RDYLGN[lo][i] + (_RDYLGN[hi][i] - _RDYLGN[lo][i]) * frac)
        for i in range(3)
    )
    return f"#{r:02x}{g:02x}{b:02x}"


def _relative_luminance(hex_color: str) -> float:
    """WCAG relative luminance of a ``#rrggbb`` string."""
    r, g, b = (int(hex_color[i:i + 2], 16) / 255.0 for i in (1, 3, 5))

    def _linear(c: float) -> float:
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    return 0.2126 * _linear(r) + 0.7152 * _linear(g) + 0.0722 * _linear(b)


def gradient_styles(values, vmin: float = 0.0, vmax: float = 3.0) -> list:
    """CSS declarations painting ``values`` onto the RdYlGn ramp.

    Drop-in for pandas' RdYlGn ``Styler.background_gradient``, applied through
    ``Styler.apply``. Values outside [vmin, vmax] clamp to the ramp's ends;
    missing values get no styling at all, so an unmeasurable cell is never
    painted as though it sat somewhere on the scale.
    """
    span = (vmax - vmin) or 1.0
    out = []
    for value in values:
        # NaN is the only value not equal to itself — avoids importing numpy
        # just to classify one cell.
        if value is None or value != value:
            out.append("")
            continue
        try:
            t = (float(value) - vmin) / span
        except (TypeError, ValueError):
            out.append("")
            continue
        bg = rdylgn_hex(t)
        fg = "#f1f1f1" if _relative_luminance(bg) < _TEXT_COLOR_THRESHOLD else "#000000"
        out.append(f"background-color: {bg}; color: {fg}")
    return out


def natural_sort_key(sample_name):
    """Extract numbers from sample name for natural sorting (e.g., Sample2 < Sample10)"""
    parts = re.split(r"(\d+)", str(sample_name))
    # `part.isascii()` guards int(): "²".isdigit() is True but int("²")
    # raises, so a sample named e.g. "20 mJ/cm²" used to abort the whole
    # script run from the upload sort, with no way to load the file.
    return [int(part) if part.isascii() and part.isdigit() else part.lower()
            for part in parts]


def get_well_exclusion_key(gene: str, sample: str) -> tuple:
    """Generate key for per-gene-sample well exclusions."""
    return (gene, sample)


# ==================== QC GRID STATE MANAGEMENT ====================
def get_grid_cell_key(gene: str, sample: str) -> str:
    """Generate unique key for grid cell identification.

    Creates a consistent, unique string key from gene and sample names.
    Used internally to identify grid cells and ensure consistent state management.

    Args:
        gene: Target gene name (e.g., "GENE1")
        sample: Sample name (e.g., "Sample_A")

    Returns:
        Unique string key (e.g., "GENE1::Sample_A")
    """
    return f"{gene}::{sample}"


def get_selected_cell(session_state) -> Optional[tuple[str, str]]:
    """Retrieve currently selected (gene, sample) or None.

    Args:
        session_state: Streamlit session state object (or dict-like)

    Returns:
        Tuple of (gene, sample) if a cell is selected, None otherwise
    """
    if "qc_grid_selected_cell" not in session_state:
        return None

    cell = session_state["qc_grid_selected_cell"]
    if cell is None:
        return None

    return (cell["gene"], cell["sample"])


def set_selected_cell(session_state, gene: str, sample: str) -> None:
    """Store the selected cell in session state.

    Args:
        session_state: Streamlit session state object (or dict-like)
        gene: Target gene name
        sample: Sample name
    """
    session_state["qc_grid_selected_cell"] = {"gene": gene, "sample": sample}


def clear_selected_cell(session_state) -> None:
    """Remove the current selection, returning the grid to an unselected state.

    Args:
        session_state: Streamlit session state object (or dict-like)
    """
    session_state["qc_grid_selected_cell"] = None


def is_cell_selected(session_state, gene: str, sample: str) -> bool:
    """Check if the given gene and sample match the currently selected cell.

    Args:
        session_state: Streamlit session state object (or dict-like)
        gene: Target gene name to check
        sample: Sample name to check

    Returns:
        True if this cell is selected, False otherwise
    """
    selected = get_selected_cell(session_state)
    if selected is None:
        return False

    return selected[0] == gene and selected[1] == sample
