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
    "Inducer": "#FFFFFF",
    "Positive Control": "#FFFFFF",
    "Treatment": "#FFFFFF",
}

GRAPH_PRESETS = {
    "Classic": {"color": "#D3D3D3", "ref": "#FFFFFF"},
    "Steel": {"color": "#4A7A9F", "ref": "#FFFFFF"},
    "Warm Neutral": {"color": "#D4B896", "ref": "#FFFFFF"},
    "Sage": {"color": "#8BAF9A", "ref": "#FFFFFF"},
    "Slate": {"color": "#8E8EA0", "ref": "#FFFFFF"},
    "Coral": {"color": "#D4826A", "ref": "#FFFFFF"},
    "Plum": {"color": "#9B7EAF", "ref": "#FFFFFF"},
    # Okabe-Ito blue — colorblind-safe (deuteranopia/protanopia/tritanopia)
    "Colorblind-Safe": {"color": "#0072B2", "ref": "#FFFFFF"},
}

FIGURE_SIZE_PRESETS = {
    "PPT Full": {"width": 28, "height": 16},
    "PPT Half": {"width": 14, "height": 10},
    "Square": {"width": 16, "height": 16},
    "Wide": {"width": 32, "height": 14},
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
# Source of truth: 효능평가항목_update.xlsx (Cosmax efficacy assay catalog).
# Single definition — the Streamlit app imports this (do not re-define inline).
# Schema per item:
#   genes: markers for the assay (canonical symbol + common spelling variants,
#          combined across the file's representative + advanced columns) — used
#          for the gene->efficacy auto-suggest and analysis.
#   cell:  cell model.
#   controls: negative = the inducer / induced-untreated baseline (shown as the
#             PPT "Inducer:" field); positive = benchmark active (omitted when
#             the file lists none); compare_to = "negative".
#   expected_direction: per-marker up/down for a positive test-article result,
#             keyed by every spelling variant so the verdict matches the data's
#             gene name. Omitted for markers/items whose direction is ambiguous
#             (lip color; cBD103; NID1) pending confirmation.
# Dict order = auto-suggest priority (mainstream assays first, niche last).
EFFICACY_CONFIG = {
    "탄력": {
        "genes": ["COL1A1", "COL1", "ELN", "FBN1", "FBN"],
        "cell": "HS68 fibroblast",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "TGFβ", "compare_to": "negative"},
        "expected_direction": {"COL1A1": "up", "COL1": "up", "ELN": "up", "FBN1": "up", "FBN": "up"},
    },
    "광노화": {
        "genes": ["MMP1", "MMP-1", "COL1A1", "COL1"],
        "cell": "HS68 fibroblast",
        "treatment_time": "24 h",
        "controls": {"negative": "UVB only", "positive": "TGFβ", "compare_to": "negative"},
        "expected_direction": {"MMP1": "down", "MMP-1": "down", "COL1A1": "up", "COL1": "up"},
    },
    "보습/수분": {
        "genes": ["HAS3", "AQP3"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "Retinoic acid", "compare_to": "negative"},
        "expected_direction": {"HAS3": "up", "AQP3": "up"},
    },
    "장벽": {
        "genes": ["FLG", "CLDN1", "CLDN", "IVL"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        # Positive control kept as Retinoic acid per Min (file listed Calcium 1.2mM).
        "controls": {"negative": "Non-treated", "positive": "Retinoic acid", "compare_to": "negative"},
        "expected_direction": {"FLG": "up", "CLDN1": "up", "CLDN": "up", "IVL": "up"},
    },
    "속보습": {
        "genes": ["HAS2", "AQP1"],
        "cell": "HS68 fibroblast",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "compare_to": "negative"},
        "expected_direction": {"HAS2": "up", "AQP1": "up"},
    },
    "멜라닌 생성": {
        "genes": ["MITF", "TYR"],
        "cell": "B16F10 melanocyte",
        "treatment_time": "24 / 48 h",
        "controls": {"negative": "α-MSH only", "positive": "Arbutin", "compare_to": "negative"},
        "expected_direction": {"MITF": "down", "TYR": "down"},
    },
    "진정": {
        "genes": ["IL1B", "IL-1β", "IL6", "IL-6", "TNF", "TNFα", "TNFA"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "4 h",
        "controls": {"negative": "Poly(I:C)+IL-4", "positive": "Dexamethasone", "compare_to": "negative"},
        "expected_direction": {
            "IL1B": "down", "IL-1β": "down", "IL6": "down", "IL-6": "down",
            "TNF": "down", "TNFα": "down", "TNFA": "down",
        },
    },
    "가려움 개선": {
        "genes": ["TSLP"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "4 h",
        "controls": {"negative": "Poly(I:C)+IL-4", "positive": "Dexamethasone", "compare_to": "negative"},
        "expected_direction": {"TSLP": "down"},
    },
    "냉감": {
        "genes": ["TRPM8", "CIRBP", "CIRP"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "Menthol", "compare_to": "negative"},
        "expected_direction": {"TRPM8": "up", "CIRBP": "up", "CIRP": "up"},
    },
    "열감": {
        "genes": ["TRPV1"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "compare_to": "negative"},
        "expected_direction": {"TRPV1": "up"},
    },
    "여드름": {
        "genes": ["PPARG", "PPARY", "PPARy", "SREBF1", "SREBP1a", "SREBP1c", "SREBPA", "SREBPC"],
        "cell": "SZ95 sebocyte",
        "treatment_time": "48 h",
        "controls": {"negative": "IGF only", "compare_to": "negative"},
        "expected_direction": {
            "PPARG": "down", "PPARY": "down", "PPARy": "down", "SREBF1": "down",
            "SREBP1a": "down", "SREBP1c": "down", "SREBPA": "down", "SREBPC": "down",
        },
    },
    "과각화": {
        "genes": ["MKI67", "Ki67", "KI67"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        "controls": {"negative": "SZ95 supernatant", "compare_to": "negative"},
        "expected_direction": {"MKI67": "down", "Ki67": "down", "KI67": "down"},
    },
    "활력": {
        "genes": ["PCNA"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "EGF", "compare_to": "negative"},
        "expected_direction": {"PCNA": "up"},
    },
    "탈모 개선": {
        # New file markers (VEGFA, HGF) + the original hair panel (COL17A1, FGF7, FLG).
        "genes": ["VEGFA", "VEGF", "HGF", "COL17A1", "FGF7", "FLG"],
        "cell": "HFDPC",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "Minoxidil", "compare_to": "negative"},
        # FLG direction in the hair/scalp context is unclear — left unannotated.
        "expected_direction": {"VEGFA": "up", "VEGF": "up", "HGF": "up", "COL17A1": "up", "FGF7": "up"},
    },
    "모공 탄력": {
        "genes": ["EMILIN1", "MFAP2", "MAGP1"],
        "cell": "HDF (neonatal / aged)",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "Niacinamide", "compare_to": "negative"},
        "expected_direction": {"EMILIN1": "up", "MFAP2": "up", "MAGP1": "up"},
    },
    "열 노화": {
        "genes": ["TRPV1", "COL6A1", "NID1"],
        "cell": "HaCaT / HS68 fibroblast",
        "treatment_time": "24 h",
        "controls": {"negative": "Heat (41°C)", "positive": "Ascorbic acid", "compare_to": "negative"},
        # NID1 direction unresolved (low confidence) — omitted pending confirmation.
        "expected_direction": {"TRPV1": "down", "COL6A1": "up"},
    },
    "민감성 기전": {
        "genes": ["TRPV1", "FLG"],
        "cell": "Sensi-TV1",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "Panthenol", "compare_to": "negative"},
        "expected_direction": {"TRPV1": "down", "FLG": "up"},
    },
    "외이도염": {
        "genes": ["CXCL8", "IL8", "IL-8", "TNF", "TNF-α", "TNFα", "CBD103", "cBD103", "TSLP"],
        "cell": "CPEK",
        "treatment_time": "6 h",
        "controls": {"negative": "Bacterial LPS", "positive": "Salicylic acid", "compare_to": "negative"},
        # cBD103 direction ambiguous (anti-inflammatory vs antimicrobial-boost) — omitted.
        "expected_direction": {
            "CXCL8": "down", "IL8": "down", "IL-8": "down",
            "TNF": "down", "TNF-α": "down", "TNFα": "down", "TSLP": "down",
        },
    },
    "구강 개선": {
        "genes": ["CLDN1", "KRT10", "K10"],
        "cell": "CGEP",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "compare_to": "negative"},
        "expected_direction": {"CLDN1": "up", "KRT10": "up", "K10": "up"},
    },
    "선번 완화": {
        "genes": ["TRPV1", "VEGFA", "VEGF", "IL6", "IL-6", "AQP3"],
        "cell": "CPEK / HaCaT",
        "treatment_time": "6 h",
        "controls": {"negative": "UVB", "compare_to": "negative"},
        "expected_direction": {
            "TRPV1": "down", "VEGFA": "down", "VEGF": "down",
            "IL6": "down", "IL-6": "down", "AQP3": "up",
        },
    },
    "립 색상": {
        "genes": ["VEGFA", "VEGF", "NOS3", "EDN1", "MC1R", "TYR"],
        "cell": "Lip explant",
        "treatment_time": "24 h",
        # Whole-item direction ambiguous (UV-protection vs vascularity claim) —
        # no expected_direction until the intended claim is confirmed.
        "controls": {"negative": "UVB", "compare_to": "negative"},
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
