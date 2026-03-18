# qPCR Graph Improvements — Design Spec

**Date:** 2026-03-18
**Author:** Min + Claude
**Status:** Approved
**Scope:** Graph visual upgrades + UI improvements for the qPCR Data Analysis Platform v1

---

## Context

The platform generates Plotly bar charts for qPCR relative expression data. Current graphs are functional but visually generic (default Plotly styling) and require excessive per-gene manual customization. The primary audience is internal lab reports and client presentations (cosmetics efficacy evaluation), not journal submissions.

### Goals
- Polished, presentation-ready graphs with minimal manual tweaking
- Better defaults that reduce per-gene customization friction
- Preserve full customizability for power users

### Non-Goals
- Journal-spec compliance (TIFF export, MIQE checklist, exact dpi requirements)
- New graph types (heatmaps, volcano plots) — deferred to future round
- Custom preset saving/import — may add later if needed

---

## Deliverables

| # | Feature | Files | Effort |
|---|---------|-------|--------|
| 1 | Color preset system (5 palettes) | constants.py, main UI | Medium |
| 2 | Significance bracket mode | graph.py | Medium |
| 3 | Individual data point overlay | graph.py, analysis.py | Medium |
| 4 | Visual polish uplift | graph.py, constants.py | Simple |
| 5 | Bug fix: error_visible_array | graph.py:207, :400 | Simple |
| 6 | "Apply to All Genes" button | main UI | Simple |
| 7 | Preset dropdown in toolbar | main UI, constants.py | Simple |
| 8 | Cleaner per-bar settings table | main UI | Medium |
| 9 | Quick-size presets | main UI, constants.py | Simple |
| 10 | Unified export button | main UI | Medium |

---

## 1. Color Preset System

### Preset Definitions (constants.py)

```python
GRAPH_PRESETS = {
    "Steel": {
        "Baseline": "#C8D6E0",
        "Treatment": "#4A7A9F",
        "Positive Control": "#2C5F7F",
        "Negative Control": "#C8D6E0",
    },
    "Warm Neutral": {
        "Baseline": "#FFFFFF",
        "Treatment": "#D4B896",
        "Positive Control": "#B89A70",
        "Negative Control": "#FFFFFF",
    },
    "Classic": {
        "Baseline": "#FFFFFF",
        "Treatment": "#D3D3D3",
        "Positive Control": "#909090",
        "Negative Control": "#FFFFFF",
    },
    "Sage": {
        "Baseline": "#E8E8E8",
        "Treatment": "#8BAF9A",
        "Positive Control": "#5C8A6E",
        "Negative Control": "#E8E8E8",
    },
    "Slate": {
        "Baseline": "#EDEDED",
        "Treatment": "#8E8EA0",
        "Positive Control": "#5B5B78",
        "Negative Control": "#EDEDED",
    },
}
```

### Behavior
- Selecting a preset overwrites all bar colors for the current gene based on each bar's Group assignment
- Per-bar color pickers still work — changing any individual color silently switches the dropdown to "Custom"
- "Apply to All Genes" respects the active preset
- Data point colors derive from bar color — darkened 30%, opacity 0.65 (calculated, not stored)

### Session State
- `graph_settings["color_preset"]` — string, one of preset names or "Custom"
- Per-gene override: `graph_settings["{gene}_color_preset"]`

---

## 2. Significance Bracket Mode

### Two Display Modes

**"Direct" (current behavior, default):**
- Asterisks/hashtags float directly above each bar
- No connecting lines
- Compact, good for many conditions

**"Bracketed" (new option):**
- Horizontal line connects the reference bar to the compared bar
- Small vertical ticks (4px) at each end of the bracket
- Significance symbol (*, **, ***, #, ##, ###) centered above the bracket
- Multiple brackets stack vertically with 12px spacing
- Y-position calculated from: max(bar_height + error_bar) + stacking_offset

### Implementation

New method in GraphGenerator:
```python
def _add_bracket_annotation(self, fig, x1_idx, x2_idx, y_base, symbol, offset_level, font_size=12):
    """Add a significance bracket between two bars."""
    # Draws: horizontal line + two vertical ticks + centered text annotation
    # Uses fig.add_shape() for lines, fig.add_annotation() for symbol
```

Bracket pairs are derived from existing significance data:
- First comparison: reference condition → each bar with significance_1
- Second comparison: second reference → each bar with significance_2

### UI
- `st.radio("Significance Style", ["Direct ✱", "Bracketed ┬"], horizontal=True)`
- Visible only when significance toggle is ON
- Session state: `graph_settings["sig_style"]` — "direct" or "bracketed"
- Per-gene override: `graph_settings["{gene}_sig_style"]`

---

## 3. Individual Data Point Overlay

### Display
- Jittered scatter dots overlaid on each bar
- Dot color: bar color darkened 30%, opacity 0.65
- Dot size: 5px (Plotly marker size)
- Jitter: random horizontal offset ±15% of bar width
- Jitter seeded by condition name for reproducibility across re-renders
- One dot per biological replicate

### Data Source
- Replicate-level fold change values already computed by AnalysisEngine
- GraphGenerator receives replicate data as additional parameter

### API Change
```python
def create_gene_graph(
    self,
    data: pd.DataFrame,
    gene: str,
    ...,
    show_data_points: bool = False,        # NEW
    replicate_data: Optional[pd.DataFrame] = None,  # NEW
) -> go.Figure:
```

### UI
- `st.toggle("Data Points", value=False)` in toolbar row 1
- Session state: `graph_settings["show_data_points"]` — bool, default False
- Per-gene override: `graph_settings["{gene}_show_data_points"]`

---

## 4. Visual Polish Uplift

### New Defaults (graph.py)

| Property | Current | New |
|----------|---------|-----|
| Axis line width | 1.2 | 1.5 |
| Axis color | Plotly default | #2C3E50 |
| Y-axis title font weight | normal | 600 (semi-bold) |
| Error bar cap width | default | 6px |
| Bar opacity | 0.95 | 0.85 |
| Gridline opacity | hidden or default | 0.3 (subtle) |

### Template
- Stays `plotly_white`
- Gridlines enabled at 0.3 opacity for subtle reference

### Font Stack
- No change — keeps existing CJK-aware fallback (Noto Sans CJK KR → Malgun Gothic → Arial)

---

## 5. Bug Fix

### Issue
`graph.py:207` and `graph.py:400` reference `error_visible_array` which is undefined. Variable was renamed to `error_visible_upper` and `error_visible_lower` during a prior refactor.

### Fix
- Line 207: `error_visible_array[idx]` → `error_visible_upper[idx]`
- Line 400: same pattern → `error_visible_upper[idx]`

---

## 6. "Apply to All Genes" Button

### Behavior
- Button in toolbar area, below gene pill selector
- Copies the active gene's full settings to all other genes:
  - Color preset (or custom colors)
  - Figure dimensions
  - Font sizes
  - Bar opacity, outline width
  - Label mode, bar gap
  - Significance style
  - Data points toggle
  - Background color
- Does NOT copy: Y-axis min/max (gene-specific), gene display names
- Confirmation: `st.warning("Apply current settings to all N genes?")` + confirm button

### Session State
- Iterates all gene keys in `graph_settings` and overwrites per-gene overrides

---

## 7. Preset Dropdown in Toolbar

### Placement
- Toolbar row 1, first position (before sig/error/bar-gap toggles)
- `st.selectbox("Color Preset", ["Steel", "Warm Neutral", "Classic", "Sage", "Slate", "Custom"])`

### Behavior
- Selecting a preset: looks up each bar's Group from the mapping data, assigns the preset color for that group
- Modifying any individual bar color via color picker: dropdown auto-switches to "Custom"
- Default: "Classic" (preserves backward compatibility)

---

## 8. Cleaner Per-Bar Settings Table

### Current Layout (6 columns)
```
Sample | Color | * | # | † | ±
```
Individual checkboxes for each significance symbol and error bar toggle.

### New Layout (3 columns)
```
Condition | Color | Options
```

- **Condition**: text label (read-only)
- **Color**: color picker (same as current)
- **Options**: `st.multiselect` with compact choices: `["✱ sig1", "# sig2", "† sig3", "± error"]`
  - Default: all enabled (matching current behavior)
  - Removing an item disables that annotation for the bar

### Benefits
- 50% fewer columns — less cramped
- Scales better with many conditions
- Easier to scan

---

## 9. Quick-Size Presets

### Size Preset Definitions (constants.py)

```python
FIGURE_SIZE_PRESETS = {
    "PPT Full": {"width": 28, "height": 16},
    "PPT Half": {"width": 14, "height": 10},
    "Square": {"width": 16, "height": 16},
    "Wide": {"width": 32, "height": 14},
    "Custom": None,  # enables manual sliders
}
```

### UI
- `st.selectbox("Figure Size", ["PPT Full", "PPT Half", "Square", "Wide", "Custom"])`
- Placed above the width/height sliders in the Settings expander
- Selecting a named preset fills both sliders with the preset values
- Modifying either slider switches dropdown to "Custom"
- Default: "PPT Full" (28×16, matches current default)

---

## 10. Unified Export Button

### Placement
- Top of Export tab, full width, primary button style
- Label: "Export All (ZIP)"

### Contents of ZIP
```
qPCR_Export_{efficacy}_{timestamp}.zip
├── qPCR_Report.xlsx           (full Excel report)
├── qPCR_Report.pptx           (PowerPoint with all genes)
├── figures/
│   ├── {gene1}.png            (300 DPI)
│   ├── {gene2}.png
│   └── ...
└── figures_html/
    └── all_graphs.html        (interactive Plotly HTML)
```

### Behavior
- Uses existing generation logic for each component
- PPT uses current experiment description defaults (editable below in the tab)
- PNG export at 300 DPI (scale=2, 1200px width)
- Progress bar during generation
- Individual export buttons remain below for selective downloads

---

## Testing Strategy

### Existing Tests to Update
- `tests/test_graph.py` — add cases for:
  - Data point overlay rendering
  - Bracket annotation positioning
  - Preset color application
  - New default styling values

### New Test Cases
- Preset application: verify bar colors match preset for each group
- Bracket mode: verify annotations are shapes (not just text) when bracketed
- Data points: verify scatter trace present when enabled, absent when disabled
- Apply-to-all: verify settings propagation across gene keys
- Size presets: verify width/height values set correctly

---

## Migration / Backward Compatibility

- "Classic" preset matches current defaults exactly — no visual change for existing users
- All new features default to OFF or to current behavior
- No session state schema changes — new keys are additive
- PPT export auto-includes new graph features (data points, brackets) when enabled
