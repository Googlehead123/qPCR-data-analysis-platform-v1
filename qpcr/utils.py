"""Utility functions for qPCR analysis.

Contains sorting helpers, well exclusion key generation, and QC grid state management.
"""

import re
from typing import Optional


def natural_sort_key(sample_name):
    """Extract numbers from sample name for natural sorting (e.g., Sample2 < Sample10)"""
    parts = re.split(r"(\d+)", str(sample_name))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


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
