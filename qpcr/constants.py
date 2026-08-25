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

# The Korean face for PowerPoint and Excel output. Matches the lab reference
# workbook +qPCR_진정_20260205_1424.xlsx, whose default font is 맑은 고딕
# (family 2, charset 129) — decision (Min, 2026-08-24), item 8a. "Malgun Gothic"
# is the same face under its Latin name, which is what OOXML wants.
#
# Every PPT text run was Arial, which has ZERO Hangul coverage, so 진정 효능 평가
# and Results: 효능 有 rendered by host substitution — visibly a different
# typeface from the chart's Korean inside the same slide. The Excel workbook
# defaulted to Calibri, same problem.
#
# Note this face is absent on Linux and Streamlit Cloud, so a deck generated
# there still falls back — but packages.txt installs fonts-noto-cjk, so it falls
# back to Noto Sans CJK KR and Korean renders either way. The reference match
# holds where it matters, which is Windows PowerPoint.
KOREAN_FONT_NAME = "Malgun Gothic"

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

# What the Graphs tab OFFERS. The dict above keeps every preset so a stored
# per-gene choice still resolves — dropping an entry would make Streamlit reset
# that widget to options[0] and silently repaint a gene (see the reference
# selectbox fix). The offered set is the two house palettes plus one neutral.
#
# Steel leads, and is the default. "Classic" was the default while being the
# WORST preset for the job it has to do: #D3D3D3 at 0.85 opacity composites to
# #dadada, 1.40:1 against the paper, and its white reference bar sits at dE 15.4
# from the treated bars — the smallest control-vs-treated separation of all
# eight (Steel is 56.7). The default was the one where the control bar is
# hardest to pick out.
#
# Colour is not a data channel here: every preset is a single tone plus a white
# reference bar, so all of them are colourblind-safe by construction and the
# "Colorblind-Safe" name implied the others were not. It stays resolvable but
# is no longer offered.
GRAPH_PRESET_CHOICES = ("Steel", "Warm Neutral", "Classic")

DEFAULT_GRAPH_PRESET = "Steel"

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

# This is a CSS font-family list handed to Plotly, so a name that is not
# installed is simply skipped by the browser — there is nothing to gain by
# filtering the list first, and the filter actively hurt. It asked matplotlib
# which families exist, and matplotlib does not index .ttc collections, which is
# exactly how the fonts-noto-cjk package in packages.txt ships Noto Sans CJK. So
# on both a dev box and Streamlit Cloud the CJK entries were always dropped and
# every chart asked for "Arial, sans-serif": Korean still rendered, but only
# because Chrome falls back per-character, never because a Korean font was
# chosen. Emitting the full preference list lets a real Korean face (Malgun
# Gothic on Windows, Apple SD Gothic Neo on macOS, Noto Sans CJK KR on Linux)
# win when it is present, and costs a matplotlib import at startup less.
PLOTLY_FONT_FAMILY = ", ".join(_DEFAULT_CJK_FONTS + _FALLBACK_FONTS)

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
#             These strings are BOTH report text and matching keys: the mapping
#             tab's `_suggest_group` and the summary's `_match` compare them
#             against real condition names from the uploaded data, exact first
#             then substring. So a benchmark carries its dose ("Arbutin
#             100 ppm") — a condition named "Arbutin" still matches on the
#             substring pass — but a name is never changed to match the xlsx,
#             because that is the half of the string the data is keyed on.
#             Dose-less entries (Minoxidil, Ascorbic acid) have none on file.
#             `장벽` positive and `탄력` positive are deliberate overrides that
#             disagree with the ledger on purpose; see their comments.
#   controls.negative_display: OPTIONAL, display-only. What the PPT "Inducer:"
#             field prints when the ledger's wording cannot be used as the
#             matching key. `negative` stays the key ("UVB only"), this carries
#             the reportable text ("UVB 40 mj/cm2"). Absent or empty falls back
#             to `negative`, so entries without it behave exactly as before.
#             Nothing but the PPT writer may read this — never match on it.
#             Why the asymmetry with `positive`, which carries its dose inline:
#             a dose appended to a positive label still matches its bare
#             condition name on the substring pass, so no separate field was
#             needed. Inducers are not so lucky — "Poly(I:C)+IL-4" is NOT a
#             substring of "Poly(I:C) 10 μg/ml+IL-4 10 ng/ml" (the doses
#             interrupt it), and "Bacterial LPS" -> "녹농균 LPS" is a rename.
#             Deliberately absent: 과각화 (the ledger cell is a medium, not a
#             concentration), 열 노화 ("Heat (41C)" already equals the key), and
#             여드름 — its cell reads "IGF 50 ng/ml or Dexamethasone 10uM", an
#             either/or resolved per experiment, so no single static string is
#             correct and printing one would be the same class of quiet-wrong
#             value as the old "1 ppm" concentration default. Left blank
#             pending a per-run input.
#   expected_direction: per-marker up/down for a positive test-article result,
#             keyed by every spelling variant so the verdict matches the data's
#             gene name. Omitted for markers/items whose direction is ambiguous
#             (lip color; cBD103; NID1) pending confirmation.
# Dict order = auto-suggest priority (mainstream assays first, niche last).
EFFICACY_CONFIG = {
    "탄력": {
        # "FBN-1" is the hyphenated spelling that appears in real plate
        # exports; without it the marker is silently ungraded in the Overview.
        "genes": ["COL1A1", "COL1", "ELN", "FBN1", "FBN-1", "FBN"],
        "cell": "HS68 fibroblast",
        "treatment_time": "24 h",
        # Positive control confirmed by Min (2026-07-31): TGF-beta at 10 ng/ml.
        # The ledger row disagrees twice and is wrong on both counts — do NOT
        # "correct" this back from the xlsx:
        #   1. it lists 10 μg/ml, a 1000x error. Min confirmed 10 ng/ml is the
        #      actual final concentration in every experiment (and it matches
        #      광노화's TGFβ 10 ng/ml).
        #   2. it lists two further actives (Retinoic acid 5 μM, Niacinamide
        #      10 μg/ml) that are not used as the 탄력 benchmark.
        "controls": {"negative": "Non-treated", "positive": "TGFβ 10 ng/ml", "compare_to": "negative"},
        "expected_direction": {"COL1A1": "up", "COL1": "up", "ELN": "up",
                               "FBN1": "up", "FBN-1": "up", "FBN": "up"},
    },
    "광노화": {
        "genes": ["MMP1", "MMP-1", "COL1A1", "COL1"],
        "cell": "HS68 fibroblast",
        "treatment_time": "24 h",
        "controls": {"negative": "UVB only", "negative_display": "UVB 40 mj/cm2",
                     "positive": "TGFβ 10 ng/ml", "compare_to": "negative"},
        "expected_direction": {"MMP1": "down", "MMP-1": "down", "COL1A1": "up", "COL1": "up"},
    },
    "보습/수분": {
        "genes": ["HAS3", "AQP3"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "Retinoic acid 1μM", "compare_to": "negative"},
        "expected_direction": {"HAS3": "up", "AQP3": "up"},
    },
    "장벽": {
        "genes": ["FLG", "CLDN1", "CLDN", "IVL", "ABCA12"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        # Positive control kept as Retinoic acid per Min (file listed Calcium 1.2mM).
        "controls": {"negative": "Non-treated", "positive": "Retinoic acid", "compare_to": "negative"},
        "expected_direction": {"FLG": "up", "CLDN1": "up", "CLDN": "up",
                               "IVL": "up", "ABCA12": "up"},
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
        "controls": {"negative": "α-MSH only", "negative_display": "α-MSH 0.1 μM",
                     "positive": "Arbutin 100 ppm", "compare_to": "negative"},
        "expected_direction": {"MITF": "down", "TYR": "down"},
    },
    "진정": {
        "genes": ["IL1B", "IL-1β", "IL6", "IL-6", "TNF", "TNFα", "TNFA"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "4 h",
        "controls": {"negative": "Poly(I:C)+IL-4",
                     "negative_display": "Poly(I:C) 10 μg/ml + IL-4 10 ng/ml",
                     "positive": "Dexamethasone 1μM", "compare_to": "negative"},
        "expected_direction": {
            "IL1B": "down", "IL-1β": "down", "IL6": "down", "IL-6": "down",
            "TNF": "down", "TNFα": "down", "TNFA": "down",
        },
    },
    "가려움 개선": {
        "genes": ["TSLP"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "4 h",
        "controls": {"negative": "Poly(I:C)+IL-4",
                     "negative_display": "Poly(I:C) 10 μg/ml + IL-4 10 ng/ml",
                     "positive": "Dexamethasone 1μM", "compare_to": "negative"},
        "expected_direction": {"TSLP": "down"},
    },
    "냉감": {
        "genes": ["TRPM8", "CIRBP", "CIRP"],
        "cell": "HaCaT keratinocyte",
        "treatment_time": "24 h",
        "controls": {"negative": "Non-treated", "positive": "Menthol 100 ppm", "compare_to": "negative"},
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
        "controls": {"negative": "Non-treated", "positive": "EGF 25 ng/ml", "compare_to": "negative"},
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
        "controls": {"negative": "Non-treated", "positive": "Niacinamide 0.1%", "compare_to": "negative"},
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
        "controls": {"negative": "Non-treated", "positive": "Panthenol 500 ppm", "compare_to": "negative"},
        "expected_direction": {"TRPV1": "down", "FLG": "up"},
    },
    "외이도염": {
        "genes": ["CXCL8", "IL8", "IL-8", "TNF", "TNF-α", "TNFα", "CBD103", "cBD103", "TSLP"],
        "cell": "CPEK",
        "treatment_time": "6 h",
        "controls": {"negative": "Bacterial LPS", "negative_display": "녹농균 LPS 1 μg/ml",
                     "positive": "Salicylic acid 1 mM", "compare_to": "negative"},
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
        # Same "UVB" matching key as 립 색상 but a different dose — the display
        # values are deliberately per-efficacy. Do not hoist to a shared constant.
        "controls": {"negative": "UVB", "negative_display": "UVB 15 mJ/cm2",
                     "compare_to": "negative"},
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
        # 100 mJ/cm2 here vs 15 mJ/cm2 for 선번 완화, same "UVB" matching key.
        "controls": {"negative": "UVB", "negative_display": "UVB 100 mJ/cm2",
                     "compare_to": "negative"},
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
