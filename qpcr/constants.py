"""Constants and configuration for qPCR analysis.

Contains efficacy database, analysis thresholds, color schemes, and font configuration.
"""

import re

# ==================== COLOR CONSTANTS ====================
DEFAULT_GROUP_COLORS = {
    "Baseline": "#FFFFFF",
    "Non-treated": "#FFFFFF",
    "Control": "#FFFFFF",
    "Negative Control": "#FFFFFF",
    "Inducer": "#909090",
    "Positive Control": "#909090",
    "Treatment": "#D3D3D3",
}

# COSMAX brand colors for PPT redesign
COSMAX_RED = "#EA1D22"
COSMAX_BLACK = "#000000"
COSMAX_WHITE = "#FFFFFF"
COSMAX_LAB_WHITE = "#F3F0ED"
COSMAX_FROST_GREY = "#C1C6C7"
COSMAX_CREAM = "#D4CEC1"

# ==================== FONT CONFIGURATION ====================
_DEFAULT_CJK_FONTS = [
    "Noto Sans CJK KR", "NanumGothic", "Malgun Gothic",
    "Apple SD Gothic Neo", "AppleGothic",
]
_FALLBACK_FONTS = ["Arial", "sans-serif"]


def _detect_available_fonts():
    """Check which CJK fonts are actually installed on this system."""
    try:
        from matplotlib import font_manager
        system_fonts = {f.name for f in font_manager.fontManager.ttflist}
        available_cjk = [f for f in _DEFAULT_CJK_FONTS if f in system_fonts]
        return available_cjk + _FALLBACK_FONTS
    except ImportError:
        return _DEFAULT_CJK_FONTS + _FALLBACK_FONTS


PLOTLY_FONT_FAMILY = ", ".join(_detect_available_fonts())

CM_TO_PX = 37.7953
CM_TO_EMU = 360000

# ==================== EFFICACY DATABASE ====================
EFFICACY_CONFIG = {
    "탄력": {
        "genes": ["COL1A1", "ELN", "FBN-1", "FBN1"],
        "cell": "HS68 fibroblast",
        "controls": {
            "negative": "Non-treated",
            "positive": "TGFb",
            "compare_to": "negative",
        },
        "description": "Elasticity - Non-treated vs TGFb (positive) vs Treatments",
    },
    "항노화": {
        "genes": ["COL1A1", "COL1", "MMP-1", "MMP1"],
        "cell": "HS68 fibroblast",
        "controls": {
            "baseline": "Non-treated (No UV)",
            "negative": "UVB only",
            "positive": "UVB+TGFb",
            "compare_to": "negative",
        },
        "description": "Anti-aging - COL1↑ (recovery), MMP1↓ (inhibition) after UVB damage",
        "expected_direction": {
            "COL1A1": "up",
            "COL1": "up",
            "MMP-1": "down",
            "MMP1": "down",
        },
    },
    "보습": {
        "genes": ["AQP3", "HAS3"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "negative": "Non-treated",
            "positive": "Retinoic acid",
            "compare_to": "negative",
        },
        "description": "Hydration - Non-treated vs Retinoic acid (positive) vs Treatments",
    },
    "장벽": {
        "genes": ["FLG", "CLDN", "IVL"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "negative": "Non-treated",
            "positive": "Retinoic acid",
            "compare_to": "negative",
        },
        "description": "Barrier function - Non-treated vs Retinoic acid (positive) vs Treatments",
    },
    "표피증식": {
        "genes": ["KI67", "PCNA"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "negative": "Non-treated",
            "positive": "TGFb or FBS",
            "compare_to": "negative",
        },
        "description": "Proliferation - Non-treated vs TGFb/FBS (positive) vs Treatments",
    },
    "멜라닌억제": {
        "genes": ["MITF", "TYR", "Melanin"],
        "cell": "B16F10 melanocyte",
        "controls": {
            "baseline": "Non-treated",
            "negative": "α-MSH only",
            "positive": "α-MSH+Arbutin",
            "compare_to": "negative",
        },
        "description": "Melanin inhibition - α-MSH induced vs α-MSH+Arbutin (positive) vs α-MSH+Treatments",
        "expected_direction": {"MITF": "down", "TYR": "down", "Melanin": "down"},
    },
    "진정": {
        "genes": ["IL1B", "IL-1β", "IL6", "TNFA", "TNFα"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "baseline": "Non-treated",
            "negative": "IL4+PolyIC (Inflammation)",
            "positive": "Inflammation+Dexamethasone",
            "compare_to": "negative",
        },
        "description": "Anti-inflammation - Reduce IL1β/IL6/TNFα (all should decrease)",
        "expected_direction": {
            "IL1B": "down",
            "IL-1β": "down",
            "IL6": "down",
            "TNFA": "down",
            "TNFα": "down",
        },
    },
    "지질억제": {
        "genes": ["SREBPA", "SREBPa", "SREBPC", "SREBPc", "PPARY", "PPARy"],
        "cell": "SZ95 sebocyte",
        "controls": {
            "baseline": "Non-treated",
            "negative": "IGF only",
            "positive": "IGF+Reference inhibitor",
            "compare_to": "negative",
        },
        "description": "Sebum inhibition - IGF induced vs IGF+Treatments",
        "expected_direction": {
            "SREBPA": "down",
            "SREBPa": "down",
            "SREBPC": "down",
            "SREBPc": "down",
            "PPARY": "down",
            "PPARy": "down",
        },
    },
    "냉감": {
        "genes": ["TRPM8", "CIRBP"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "negative": "Non-treated",
            "positive": "Menthol",
            "compare_to": "negative",
        },
        "description": "Cooling effect - Non-treated vs Menthol (positive) vs Treatments",
    },
    "모근 강화": {
        "genes": ["VEGF", "COL17A1", "HGF", "FGF7", "FLG"],
        "cell": "HFDPC / HaCaT",
        "controls": {
            "negative": "Non-treated",
            "positive": "Minoxidil",
            "compare_to": "negative",
        },
        "description": "Hair root strengthening - VEGF, COL17A1, HGF, FGF7, FLG expression",
    },
}


# ==================== ANALYSIS CONSTANTS ====================
class AnalysisConstants:
    MIN_REPLICATES_FOR_STATS = 2
    RECOMMENDED_REPLICATES = 3
    CT_UNDETERMINED_THRESHOLD = 40.0
    CT_HIGH_WARNING = 35.0
    CT_LOW_WARNING = 10.0
    CV_WARNING_THRESHOLD = 0.05
    CV_ERROR_THRESHOLD = 0.10
    HK_DEVIATION_WARNING = 1.0
    HK_DEVIATION_ERROR = 2.0
    P_VALUE_THRESHOLDS = {"*": 0.05, "**": 0.01, "***": 0.001}
