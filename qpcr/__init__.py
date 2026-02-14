"""qPCR Data Analysis Package.

Modular extraction of the qPCR Analysis Suite. Provides:
- QPCRParser: Raw CSV parsing (QuantStudio formats)
- QualityControl: QC checks, outlier detection, triplicate stats
- AnalysisEngine: DDCt calculations and statistical tests
- GraphGenerator: Plotly bar chart visualizations
- ReportGenerator: Legacy PowerPoint generation
- PPTGenerator: Modern Korean-centric PowerPoint generation
- export_to_excel: Multi-sheet Excel export
"""

from qpcr.constants import (
    EFFICACY_CONFIG,
    AnalysisConstants,
    DEFAULT_GROUP_COLORS,
    COSMAX_RED,
    COSMAX_BLACK,
    COSMAX_WHITE,
    COSMAX_LAB_WHITE,
    COSMAX_FROST_GREY,
    COSMAX_CREAM,
    PLOTLY_FONT_FAMILY,
    CM_TO_PX,
    CM_TO_EMU,
)
from qpcr.utils import (
    natural_sort_key,
    get_well_exclusion_key,
    get_grid_cell_key,
    get_selected_cell,
    set_selected_cell,
    clear_selected_cell,
    is_cell_selected,
)
from qpcr.parser import QPCRParser
from qpcr.quality_control import QualityControl
from qpcr.analysis import AnalysisEngine
from qpcr.graph import GraphGenerator
from qpcr.report import ReportGenerator, PPTGenerator
from qpcr.export import export_to_excel

__all__ = [
    "EFFICACY_CONFIG",
    "AnalysisConstants",
    "DEFAULT_GROUP_COLORS",
    "COSMAX_RED",
    "COSMAX_BLACK",
    "COSMAX_WHITE",
    "COSMAX_LAB_WHITE",
    "COSMAX_FROST_GREY",
    "COSMAX_CREAM",
    "PLOTLY_FONT_FAMILY",
    "CM_TO_PX",
    "CM_TO_EMU",
    "natural_sort_key",
    "get_well_exclusion_key",
    "get_grid_cell_key",
    "get_selected_cell",
    "set_selected_cell",
    "clear_selected_cell",
    "is_cell_selected",
    "QPCRParser",
    "QualityControl",
    "AnalysisEngine",
    "GraphGenerator",
    "ReportGenerator",
    "PPTGenerator",
    "export_to_excel",
]
