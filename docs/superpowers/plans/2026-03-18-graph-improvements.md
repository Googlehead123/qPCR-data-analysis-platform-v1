# Graph Improvements Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Upgrade qPCR graph visuals with color presets, significance brackets, data point overlays, and UI improvements for presentation-ready output.

**Architecture:** All changes stay within the existing Streamlit + Plotly architecture. New constants go in `qpcr/constants.py`, graph rendering changes in `qpcr/graph.py`, new analysis method in `qpcr/analysis.py`, UI changes in the main Streamlit file. No new dependencies.

**Tech Stack:** Python 3.12, Streamlit, Plotly, pandas, numpy, scipy

**Spec:** `docs/superpowers/specs/2026-03-18-graph-improvements-design.md`

---

### Task 1: Bug Fix — `error_visible_array` undefined

**Files:**
- Modify: `qpcr/graph.py:207` and `qpcr/graph.py:400`
- Test: `tests/test_graph.py`

- [ ] **Step 1: Write a failing test that triggers the bug**

In `tests/test_graph.py`, add a test that provides significance data — this hits the code path at line 207:

```python
class TestGraphGeneratorBugFixes:
    def test_significance_symbols_do_not_crash(self, mock_streamlit, graph_settings):
        from qpcr.graph import GraphGenerator
        import plotly.graph_objects as go

        data = pd.DataFrame({
            "Target": ["COL1A1", "COL1A1", "COL1A1"],
            "Condition": ["Non-treated", "Treatment1", "Treatment2"],
            "Group": ["Negative Control", "Treatment", "Treatment"],
            "Relative_Expression": [1.0, 2.5, 0.4],
            "SEM": [0.1, 0.2, 0.05],
            "FC_Error_Upper": [0.15, 0.3, 0.08],
            "FC_Error_Lower": [0.12, 0.25, 0.06],
            "significance": ["", "**", "*"],
            "significance_2": ["", "", ""],
        })
        graph_settings["show_significance"] = True
        graph_settings["show_error"] = True

        fig = GraphGenerator.create_gene_graph(
            data=data, gene="COL1A1", settings=graph_settings
        )
        assert isinstance(fig, go.Figure)
        # Should have significance annotations
        assert len(fig.layout.annotations) >= 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m pytest tests/test_graph.py::TestGraphGeneratorBugFixes::test_significance_symbols_do_not_crash -v`
Expected: FAIL with `NameError: name 'error_visible_array' is not defined`

- [ ] **Step 3: Fix the bug**

In `qpcr/graph.py`, line 207:
```python
# BEFORE:
error_bar_height = error_visible_array[idx]
# AFTER:
error_bar_height = error_visible_upper[idx]
```

In `qpcr/graph.py`, line 400:
```python
# BEFORE:
_err_h = error_visible_array[_idx] if _idx < len(error_visible_array) else 0
# AFTER:
_err_h = error_visible_upper[_idx] if _idx < len(error_visible_upper) else 0
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m pytest tests/test_graph.py -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add qpcr/graph.py tests/test_graph.py
git commit -m "fix: resolve NameError for error_visible_array in graph.py"
```

---

### Task 2: Color Preset Definitions

**Files:**
- Modify: `qpcr/constants.py` (add after `DEFAULT_GROUP_COLORS` block, ~line 17)
- Test: `tests/test_graph.py`

- [ ] **Step 1: Write a test for preset structure**

```python
class TestColorPresets:
    def test_all_presets_cover_all_group_keys(self):
        from qpcr.constants import GRAPH_PRESETS, DEFAULT_GROUP_COLORS
        required_groups = set(DEFAULT_GROUP_COLORS.keys())
        for preset_name, colors in GRAPH_PRESETS.items():
            assert required_groups.issubset(set(colors.keys())), \
                f"Preset '{preset_name}' missing groups: {required_groups - set(colors.keys())}"

    def test_classic_preset_matches_defaults(self):
        from qpcr.constants import GRAPH_PRESETS, DEFAULT_GROUP_COLORS
        for group, color in DEFAULT_GROUP_COLORS.items():
            assert GRAPH_PRESETS["Classic"][group] == color, \
                f"Classic preset mismatch for '{group}': {GRAPH_PRESETS['Classic'][group]} != {color}"

    def test_all_preset_colors_are_valid_hex(self):
        import re
        from qpcr.constants import GRAPH_PRESETS
        hex_re = re.compile(r'^#[0-9A-Fa-f]{6}$')
        for preset_name, colors in GRAPH_PRESETS.items():
            for group, color in colors.items():
                assert hex_re.match(color), f"Invalid hex in {preset_name}/{group}: {color}"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_graph.py::TestColorPresets -v`
Expected: FAIL with `ImportError: cannot import name 'GRAPH_PRESETS'`

- [ ] **Step 3: Add preset definitions to constants.py**

Add after line 17 in `qpcr/constants.py` (after `DEFAULT_GROUP_COLORS`):

```python
GRAPH_PRESETS = {
    "Steel": {
        "Baseline": "#C8D6E0",
        "Non-treated": "#C8D6E0",
        "Control": "#C8D6E0",
        "Negative Control": "#C8D6E0",
        "Treatment": "#4A7A9F",
        "Inducer": "#2C5F7F",
        "Positive Control": "#2C5F7F",
    },
    "Warm Neutral": {
        "Baseline": "#FFFFFF",
        "Non-treated": "#FFFFFF",
        "Control": "#FFFFFF",
        "Negative Control": "#FFFFFF",
        "Treatment": "#D4B896",
        "Inducer": "#B89A70",
        "Positive Control": "#B89A70",
    },
    "Classic": {
        "Baseline": "#FFFFFF",
        "Non-treated": "#FFFFFF",
        "Control": "#FFFFFF",
        "Negative Control": "#FFFFFF",
        "Treatment": "#D3D3D3",
        "Inducer": "#909090",
        "Positive Control": "#909090",
    },
    "Sage": {
        "Baseline": "#E8E8E8",
        "Non-treated": "#E8E8E8",
        "Control": "#E8E8E8",
        "Negative Control": "#E8E8E8",
        "Treatment": "#8BAF9A",
        "Inducer": "#5C8A6E",
        "Positive Control": "#5C8A6E",
    },
    "Slate": {
        "Baseline": "#EDEDED",
        "Non-treated": "#EDEDED",
        "Control": "#EDEDED",
        "Negative Control": "#EDEDED",
        "Treatment": "#8E8EA0",
        "Inducer": "#5B5B78",
        "Positive Control": "#5B5B78",
    },
}

FIGURE_SIZE_PRESETS = {
    "PPT Full": {"width": 28, "height": 16},
    "PPT Half": {"width": 14, "height": 10},
    "Square": {"width": 16, "height": 16},
    "Wide": {"width": 32, "height": 14},
}
```

Also update `qpcr/__init__.py` — add to the imports from constants:
```python
from qpcr.constants import GRAPH_PRESETS, FIGURE_SIZE_PRESETS
```
And add both names to the `__all__` list if one exists.

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_graph.py::TestColorPresets -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add qpcr/constants.py qpcr/__init__.py tests/test_graph.py
git commit -m "feat: add color and figure size presets to constants"
```

---

### Task 3: Visual Polish Uplift

**Files:**
- Modify: `qpcr/graph.py:155-177` (bar trace), `qpcr/graph.py:268-284` (y-axis), `qpcr/graph.py:343-377` (layout)
- Test: `tests/test_graph.py`

- [ ] **Step 1: Write tests for new default styling**

```python
class TestVisualPolish:
    def test_axis_color_is_dark(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings
        )
        assert fig.layout.yaxis.linecolor == "#2C3E50"

    def test_axis_line_width(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings
        )
        assert fig.layout.yaxis.linewidth == 1.5

    def test_default_bar_opacity(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        # Don't set bar_opacity in settings — test the default
        graph_settings.pop("bar_opacity", None)
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings
        )
        assert fig.data[0].marker.opacity == 0.85

    def test_error_bar_cap_width(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        graph_settings["show_error"] = True
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings
        )
        assert fig.data[0].error_y.width == 6
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_graph.py::TestVisualPolish -v`
Expected: FAIL (current defaults don't match)

- [ ] **Step 3: Update graph.py defaults**

In `qpcr/graph.py`, update the `create_gene_graph` method:

1. **Error bar width** (~line 165): change `width=4` → `width=6`
2. **Bar opacity default** (~line 174): change `settings.get("bar_opacity", 0.95)` → `settings.get("bar_opacity", 0.85)`
3. **Y-axis config** (~line 278-281): change:
   ```python
   showline=True,
   linewidth=1,
   linecolor="black",
   ```
   to:
   ```python
   showline=True,
   linewidth=1.5,
   linecolor="#2C3E50",
   ```
4. **Y-axis title font** (~line 270): Plotly's `Font` object does not support `weight`. Instead, wrap the Y-axis title in `<b>` tags. The existing `y_label_html` already uses HTML; change to: `y_label_html = f"<b>Relative <span style='color:red;'>{gene_label}</span> Expression Level</b>"`
5. **Gridline** (~line 274): change `showgrid=False` → `showgrid=True, gridcolor="rgba(0,0,0,0.08)", gridwidth=0.5`

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_graph.py::TestVisualPolish -v`
Expected: ALL PASS

- [ ] **Step 5: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS (no regressions)

- [ ] **Step 6: Commit**

```bash
git add qpcr/graph.py tests/test_graph.py
git commit -m "feat: visual polish — axis styling, bar opacity, gridlines, error bar caps"
```

---

### Task 4: Replicate Fold Change Computation

**Files:**
- Modify: `qpcr/analysis.py` (add new method after `calculate_ddct`)
- Test: `tests/test_analysis.py`

- [ ] **Step 1: Write test for replicate fold changes**

Read `tests/test_analysis.py` first to understand existing test patterns, then add:

```python
class TestReplicateFoldChanges:
    def test_returns_one_row_per_replicate(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from qpcr.analysis import AnalysisEngine
        result = AnalysisEngine.compute_replicate_fold_changes(
            raw_data=sample_qpcr_raw_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            sample_mapping=sample_mapping,
            excluded_wells=set(),
        )
        # 3 samples × 3 replicates each for COL1A1 = 9 rows
        assert len(result) == 9
        assert "Condition" in result.columns
        assert "Replicate_FC" in result.columns
        assert "Target" in result.columns

    def test_reference_condition_fcs_average_near_one(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from qpcr.analysis import AnalysisEngine
        result = AnalysisEngine.compute_replicate_fold_changes(
            raw_data=sample_qpcr_raw_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            sample_mapping=sample_mapping,
            excluded_wells=set(),
        )
        ref_fcs = result[result["Condition"] == "Non-treated"]["Replicate_FC"]
        assert abs(ref_fcs.mean() - 1.0) < 0.3  # Mean should be ~1.0

    def test_handles_excluded_wells(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from qpcr.analysis import AnalysisEngine
        # Exclude a COL1A1 well (index 3 = first COL1A1 replicate, after 3 GAPDH wells)
        excluded = {sample_qpcr_raw_data["Well"].iloc[3]}
        result = AnalysisEngine.compute_replicate_fold_changes(
            raw_data=sample_qpcr_raw_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            sample_mapping=sample_mapping,
            excluded_wells=excluded,
        )
        assert len(result) == 8  # 9 - 1 excluded COL1A1 well
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_analysis.py::TestReplicateFoldChanges -v`
Expected: FAIL with `AttributeError: type object 'AnalysisEngine' has no attribute 'compute_replicate_fold_changes'`

- [ ] **Step 3: Implement the method**

Add to `qpcr/analysis.py` after `calculate_ddct` method:

```python
@staticmethod
def compute_replicate_fold_changes(
    raw_data: pd.DataFrame,
    hk_gene: str,
    ref_sample: str,
    sample_mapping: dict,
    excluded_wells=None,
) -> pd.DataFrame:
    """Compute per-replicate fold change values for data point overlay.

    For each replicate i in condition c:
      dCt_i = Ct_target_i - Ct_hk_mean_of_condition
      ddCt_i = dCt_i - dCt_mean_of_reference
      FC_i = 2^(-ddCt_i)

    Uses per-condition HK mean (not per-replicate HK) to isolate
    target gene variability from housekeeping variability.
    """
    data = raw_data.copy()

    # Apply exclusions
    if isinstance(excluded_wells, dict):
        # Per-gene-sample exclusion dict
        mask = data.apply(
            lambda r: r["Well"] not in excluded_wells.get(
                (r["Target"], r["Sample"]), set()
            ),
            axis=1,
        )
        data = data[mask]
    elif excluded_wells:
        data = data[~data["Well"].isin(excluded_wells)]

    # Map conditions
    data["Condition"] = data["Sample"].map(
        lambda x: sample_mapping.get(x, {}).get("condition", x)
    )

    # Compute per-condition HK means
    hk_data = data[data["Target"].str.upper() == hk_gene.upper()]
    hk_means = hk_data.groupby("Condition")["CT"].mean().to_dict()

    # Compute reference dCt (mean across reference condition)
    results = []
    for target in data["Target"].unique():
        if target.upper() == hk_gene.upper():
            continue

        target_data = data[data["Target"] == target]

        # Reference dCt
        ref_rows = target_data[target_data["Condition"] == ref_sample]
        ref_hk_mean = hk_means.get(ref_sample, np.nan)
        if np.isnan(ref_hk_mean) or ref_rows.empty:
            continue
        ref_dct_mean = ref_rows["CT"].mean() - ref_hk_mean

        # Per-replicate fold changes
        for _, row in target_data.iterrows():
            condition = row["Condition"]
            cond_hk_mean = hk_means.get(condition, np.nan)
            if np.isnan(cond_hk_mean):
                continue
            dct_i = row["CT"] - cond_hk_mean
            ddct_i = dct_i - ref_dct_mean
            ddct_clamped = np.clip(ddct_i, -50, 50)
            fc_i = 2 ** (-ddct_clamped)
            results.append({
                "Target": target,
                "Condition": condition,
                "Well": row["Well"],
                "Replicate_FC": fc_i,
            })

    return pd.DataFrame(results)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_analysis.py::TestReplicateFoldChanges -v`
Expected: ALL PASS

- [ ] **Step 5: Commit**

```bash
git add qpcr/analysis.py tests/test_analysis.py
git commit -m "feat: add compute_replicate_fold_changes for data point overlay"
```

---

### Task 5: Data Point Overlay in GraphGenerator

**Files:**
- Modify: `qpcr/graph.py` (add scatter trace in `create_gene_graph`)
- Test: `tests/test_graph.py`

- [ ] **Step 1: Write tests for data point overlay**

```python
class TestDataPointOverlay:
    def test_no_scatter_trace_when_disabled(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings,
            show_data_points=False,
        )
        # Only bar trace, no scatter
        assert len(fig.data) == 1

    def test_scatter_trace_when_enabled(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        replicate_data = pd.DataFrame({
            "Target": ["COL1A1"] * 9,
            "Condition": ["Non-treated"] * 3 + ["Treatment1"] * 3 + ["Treatment2"] * 3,
            "Well": [f"A{i}" for i in range(1, 10)],
            "Replicate_FC": [0.95, 1.02, 1.03, 2.4, 2.6, 2.55, 0.35, 0.38, 0.37],
        })
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings,
            show_data_points=True, replicate_data=replicate_data,
        )
        # Bar trace + scatter trace
        assert len(fig.data) == 2
        assert fig.data[1].mode == "markers"

    def test_no_scatter_when_enabled_but_no_replicate_data(self, mock_streamlit, processed_gene_data, graph_settings):
        from qpcr.graph import GraphGenerator
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data, gene="COL1A1", settings=graph_settings,
            show_data_points=True, replicate_data=None,
        )
        assert len(fig.data) == 1
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_graph.py::TestDataPointOverlay -v`
Expected: FAIL with `TypeError: create_gene_graph() got an unexpected keyword argument 'show_data_points'`

- [ ] **Step 3: Add data point parameters and scatter trace to graph.py**

In `qpcr/graph.py`, modify `create_gene_graph` signature (add two params):
```python
def create_gene_graph(
    data, gene, settings, efficacy_config=None, sample_order=None,
    per_sample_overrides=None, condition_colors=None, display_gene_name=None,
    ref_line_value=None, ref_line_label=None,
    show_data_points: bool = False,          # NEW
    replicate_data: pd.DataFrame = None,     # NEW
) -> go.Figure:
```

After the bar trace (`fig.add_trace(go.Bar(...))`, ~line 178), add:

```python
# Data point overlay (jittered scatter on top of bars)
if show_data_points and replicate_data is not None and not replicate_data.empty:
    import hashlib
    scatter_x = []
    scatter_y = []
    scatter_colors = []
    for idx, condition in enumerate(condition_names):
        cond_replicates = replicate_data[replicate_data["Condition"] == condition]
        if cond_replicates.empty:
            continue
        # Deterministic jitter seeded by gene + condition
        seed = int(hashlib.md5(f"{gene}_{condition}".encode()).hexdigest()[:8], 16)
        rng = np.random.RandomState(seed)
        n_pts = len(cond_replicates)
        jitter = rng.uniform(-0.15, 0.15, size=n_pts)
        scatter_x.extend([idx + j for j in jitter])
        scatter_y.extend(cond_replicates["Replicate_FC"].tolist())
        # Darken bar color by 30%
        base_color = bar_colors[idx] if idx < len(bar_colors) else "#666666"
        scatter_colors.extend([_darken_hex(base_color, 0.3)] * n_pts)

    if scatter_x:
        fig.add_trace(go.Scatter(
            x=scatter_x,
            y=scatter_y,
            mode="markers",
            marker=dict(
                size=5,
                color=scatter_colors,
                opacity=0.65,
                line=dict(width=0),
            ),
            showlegend=False,
            hoverinfo="y",
        ))
```

Also add a helper function at the top of the file (after imports):

```python
def _darken_hex(hex_color: str, factor: float = 0.3) -> str:
    """Darken a hex color by the given factor (0-1)."""
    hex_color = hex_color.lstrip("#")
    r, g, b = int(hex_color[:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    r = int(r * (1 - factor))
    g = int(g * (1 - factor))
    b = int(b * (1 - factor))
    return f"#{r:02x}{g:02x}{b:02x}"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_graph.py::TestDataPointOverlay -v`
Expected: ALL PASS

- [ ] **Step 5: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add qpcr/graph.py tests/test_graph.py
git commit -m "feat: add data point overlay with jittered scatter dots"
```

---

### Task 6: Significance Bracket Mode

**Files:**
- Modify: `qpcr/graph.py` (add `_add_bracket_annotation` method, modify significance rendering)
- Test: `tests/test_graph.py`

- [ ] **Step 1: Write tests for bracket mode**

```python
class TestSignificanceBrackets:
    def test_direct_mode_uses_annotations(self, mock_streamlit, graph_settings):
        from qpcr.graph import GraphGenerator
        import plotly.graph_objects as go
        data = pd.DataFrame({
            "Target": ["COL1A1", "COL1A1"],
            "Condition": ["Non-treated", "Treatment1"],
            "Group": ["Negative Control", "Treatment"],
            "Relative_Expression": [1.0, 2.5],
            "SEM": [0.1, 0.2],
            "FC_Error_Upper": [0.15, 0.3],
            "FC_Error_Lower": [0.12, 0.25],
            "significance": ["", "**"],
            "significance_2": ["", ""],
        })
        graph_settings["show_significance"] = True
        graph_settings["sig_style"] = "direct"
        fig = GraphGenerator.create_gene_graph(data=data, gene="COL1A1", settings=graph_settings)
        # Direct mode: annotations only, no shapes for brackets
        sig_annotations = [a for a in fig.layout.annotations if a.text in ["*", "**", "***"]]
        assert len(sig_annotations) >= 1

    def test_bracketed_mode_uses_shapes(self, mock_streamlit, graph_settings):
        from qpcr.graph import GraphGenerator
        data = pd.DataFrame({
            "Target": ["COL1A1", "COL1A1"],
            "Condition": ["Non-treated", "Treatment1"],
            "Group": ["Negative Control", "Treatment"],
            "Relative_Expression": [1.0, 2.5],
            "SEM": [0.1, 0.2],
            "FC_Error_Upper": [0.15, 0.3],
            "FC_Error_Lower": [0.12, 0.25],
            "significance": ["", "**"],
            "significance_2": ["", ""],
        })
        graph_settings["show_significance"] = True
        graph_settings["sig_style"] = "bracketed"
        mock_streamlit.session_state["analysis_cmp_condition"] = "Non-treated"
        fig = GraphGenerator.create_gene_graph(data=data, gene="COL1A1", settings=graph_settings)
        # Bracketed mode: should have shapes (lines) for brackets
        bracket_shapes = [s for s in (fig.layout.shapes or []) if s.type == "line"]
        assert len(bracket_shapes) >= 3  # 1 horizontal + 2 vertical ticks per bracket

    def test_bracket_fallback_when_too_many_comparisons(self, mock_streamlit, graph_settings):
        from qpcr.graph import GraphGenerator
        # 8 conditions = 7 comparisons > 6 max brackets → should fall back to direct
        conditions = ["Non-treated"] + [f"T{i}" for i in range(1, 8)]
        data = pd.DataFrame({
            "Target": ["COL1A1"] * 8,
            "Condition": conditions,
            "Group": ["Negative Control"] + ["Treatment"] * 7,
            "Relative_Expression": [1.0] + [float(i) for i in range(2, 9)],
            "SEM": [0.1] * 8,
            "FC_Error_Upper": [0.15] * 8,
            "FC_Error_Lower": [0.12] * 8,
            "significance": [""] + ["**"] * 7,
            "significance_2": [""] * 8,
        })
        graph_settings["show_significance"] = True
        graph_settings["sig_style"] = "bracketed"
        mock_streamlit.session_state["analysis_cmp_condition"] = "Non-treated"
        fig = GraphGenerator.create_gene_graph(data=data, gene="COL1A1", settings=graph_settings)
        # Should fall back to direct (annotations, no bracket shapes)
        bracket_shapes = [s for s in (fig.layout.shapes or []) if s.type == "line"]
        assert len(bracket_shapes) == 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_graph.py::TestSignificanceBrackets -v`
Expected: FAIL

- [ ] **Step 3: Implement bracket mode**

In `qpcr/graph.py`, add a new method to `GraphGenerator`:

```python
@staticmethod
def _add_bracket_annotation(fig, x1, x2, y_base, symbol, offset_level,
                            font_size=12, spacing=0.05):
    """Add a significance bracket between two x positions."""
    y = y_base + (offset_level * spacing * y_base)
    tick_h = spacing * y_base * 0.3

    # Horizontal line
    fig.add_shape(type="line", x0=x1, x1=x2, y0=y, y1=y,
                  line=dict(color="#2C3E50", width=1), xref="x", yref="y")
    # Left tick
    fig.add_shape(type="line", x0=x1, x1=x1, y0=y - tick_h, y1=y,
                  line=dict(color="#2C3E50", width=1), xref="x", yref="y")
    # Right tick
    fig.add_shape(type="line", x0=x2, x1=x2, y0=y - tick_h, y1=y,
                  line=dict(color="#2C3E50", width=1), xref="x", yref="y")
    # Symbol centered above
    fig.add_annotation(
        x=(x1 + x2) / 2, y=y + (tick_h * 0.5),
        text=symbol, showarrow=False,
        font=dict(size=font_size, color="#2C3E50", family=PLOTLY_FONT_FAMILY),
        xref="x", yref="y", xanchor="center", yanchor="bottom",
    )
```

Then modify the significance annotation block (~lines 194-263). Wrap the existing direct-mode logic in a conditional:

```python
sig_style = settings.get("sig_style", "direct")

# Count significant comparisons to check bracket limit
if sig_style == "bracketed":
    ref_condition = st.session_state.get("analysis_cmp_condition", "")
    sig_count = sum(1 for idx in range(n_bars)
                    if gene_data_indexed.iloc[idx].get("significance", "") in ["*", "**", "***"])
    if sig_count > 6:
        sig_style = "direct"  # fallback

if sig_style == "bracketed" and show_sig_global:
    # Find reference bar index
    ref_condition = st.session_state.get("analysis_cmp_condition", "")
    ref_idx = None
    for idx in range(n_bars):
        if gene_data_indexed.iloc[idx]["Condition"] == ref_condition:
            ref_idx = idx
            break

    if ref_idx is not None:
        offset = 0
        for idx in range(n_bars):
            row = gene_data_indexed.iloc[idx]
            bar_key = f"{gene}_{row['Condition']}"
            bar_config = gene_bar_settings.get(bar_key, {"show_sig": True})
            sig = row.get("significance", "")
            if sig in ["*", "**", "***"] and bar_config.get("show_sig", True):
                max_bar_y = max(
                    gene_data_indexed.iloc[ref_idx]["Relative_Expression"] + error_visible_upper[ref_idx],
                    row["Relative_Expression"] + error_visible_upper[idx],
                )
                GraphGenerator._add_bracket_annotation(
                    fig, ref_idx, idx, max_bar_y * 1.08, sig, offset
                )
                offset += 1

        # Second comparison brackets (if present)
        ref_condition_2 = st.session_state.get("analysis_cmp_condition_2", "")
        if ref_condition_2 and "significance_2" in gene_data_indexed.columns:
            ref_idx_2 = None
            for idx in range(n_bars):
                if gene_data_indexed.iloc[idx]["Condition"] == ref_condition_2:
                    ref_idx_2 = idx
                    break
            if ref_idx_2 is not None:
                for idx in range(n_bars):
                    row = gene_data_indexed.iloc[idx]
                    sig2 = row.get("significance_2", "")
                    if sig2 in ["#", "##", "###"]:
                        bar_key = f"{gene}_{row['Condition']}"
                        bar_config = gene_bar_settings.get(bar_key, {"show_sig": True})
                        if bar_config.get("show_sig", True):
                            max_bar_y = max(
                                gene_data_indexed.iloc[ref_idx_2]["Relative_Expression"] + error_visible_upper[ref_idx_2],
                                row["Relative_Expression"] + error_visible_upper[idx],
                            )
                            GraphGenerator._add_bracket_annotation(
                                fig, ref_idx_2, idx, max_bar_y * 1.08, sig2, offset, font_size=10
                            )
                            offset += 1
else:
    # Original direct-mode significance annotation code (existing lines 195-263)
    ...
```

Also update `y_max_auto` to account for bracket stacking when in bracket mode.

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_graph.py::TestSignificanceBrackets -v`
Expected: ALL PASS

- [ ] **Step 5: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add qpcr/graph.py tests/test_graph.py
git commit -m "feat: add significance bracket mode with auto-fallback"
```

---

### Task 7: Preset Dropdown + Apply-to-All Button in UI

**Files:**
- Modify: `streamlit qpcr analysis v1.py:5489-5529` (toolbar section)

- [ ] **Step 1: Add preset dropdown to toolbar row 1**

In the main file, modify toolbar row 1 (~line 5490). Change column layout from `[1.2, 1.2, 1.5, 0.8]` to `[1.5, 1.0, 1.0, 1.5, 0.7]` to fit the preset dropdown:

```python
from qpcr.constants import GRAPH_PRESETS, FIGURE_SIZE_PRESETS

# Initialize preset key
color_preset_key = f"{current_gene}_color_preset"
if color_preset_key not in st.session_state.graph_settings:
    st.session_state.graph_settings[color_preset_key] = "Classic"

tb_row1 = st.columns([1.5, 1.0, 1.0, 1.5, 0.7])
with tb_row1[0]:
    preset_names = list(GRAPH_PRESETS.keys()) + ["Custom"]
    current_preset = st.session_state.graph_settings.get(color_preset_key, "Classic")
    if current_preset not in preset_names:
        current_preset = "Custom"
    selected_preset = st.selectbox(
        "Color Preset",
        preset_names,
        index=preset_names.index(current_preset),
        key=f"preset_{current_gene}",
    )
    if selected_preset != "Custom" and selected_preset != st.session_state.graph_settings.get(color_preset_key):
        # Apply preset colors to all bars
        preset_colors = GRAPH_PRESETS[selected_preset]
        for _, row in gene_data.iterrows():
            condition = row["Condition"]
            group = row.get("Group", "Treatment")
            bar_key = f"{current_gene}_{condition}"
            color = preset_colors.get(group, preset_colors.get("Treatment", "#D3D3D3"))
            if bar_key in st.session_state.get(f"{current_gene}_bar_settings", {}):
                st.session_state[f"{current_gene}_bar_settings"][bar_key]["color"] = color
            st.session_state.graph_settings["bar_colors_per_sample"][bar_key] = color
    st.session_state.graph_settings[color_preset_key] = selected_preset
# ... sig toggle, error toggle, bar gap, reset in remaining columns
```

- [ ] **Step 2: Add sig style radio**

After toolbar row 1 toggles, add significance style radio (only visible when significance is ON):

```python
if st.session_state.graph_settings.get(show_sig_key, True):
    sig_style_key = f"{current_gene}_sig_style"
    if sig_style_key not in st.session_state.graph_settings:
        st.session_state.graph_settings[sig_style_key] = "direct"
    sig_style = st.radio(
        "Significance Style",
        ["Direct ✱", "Bracketed ┬"],
        index=0 if st.session_state.graph_settings.get(sig_style_key, "direct") == "direct" else 1,
        key=f"sig_style_{current_gene}",
        horizontal=True,
    )
    st.session_state.graph_settings[sig_style_key] = "direct" if "Direct" in sig_style else "bracketed"
```

- [ ] **Step 3: Add data points toggle**

Add to toolbar row 1 (or row 2):

```python
show_dp_key = f"{current_gene}_show_data_points"
if show_dp_key not in st.session_state.graph_settings:
    st.session_state.graph_settings[show_dp_key] = False
dp_on = st.toggle(
    "Data Points",
    st.session_state.graph_settings[show_dp_key],
    key=f"tgl_dp_{current_gene}",
)
st.session_state.graph_settings[show_dp_key] = dp_on
```

- [ ] **Step 4: Add "Apply to All Genes" button**

After the gene pill selector (~line 5465), add:

```python
PER_GENE_SUFFIXES = [
    "_color_preset", "_figure_width", "_figure_height", "_font_size",
    "_tick_size", "_ylabel_size", "_bar_opacity", "_marker_line_width",
    "_bg_color", "_bar_gap", "_sig_style", "_show_data_points",
    "_label_mode", "_ref_line", "_show_sig", "_show_err",
]

if len(gene_list) > 1:
    apply_col1, apply_col2 = st.columns([3, 1])
    with apply_col2:
        if st.button(f"Apply to All {len(gene_list)} Genes", key="apply_all_genes", use_container_width=True):
            source_gene = current_gene
            source_preset = st.session_state.graph_settings.get(f"{source_gene}_color_preset", "Classic")
            for target_gene in gene_list:
                if target_gene == source_gene:
                    continue
                # Copy per-gene settings
                for suffix in PER_GENE_SUFFIXES:
                    src_key = f"{source_gene}{suffix}"
                    tgt_key = f"{target_gene}{suffix}"
                    if src_key in st.session_state.graph_settings:
                        st.session_state.graph_settings[tgt_key] = st.session_state.graph_settings[src_key]
                # Apply preset colors to target gene
                if source_preset != "Custom":
                    preset_colors = GRAPH_PRESETS.get(source_preset, {})
                    target_data = st.session_state.processed_data.get(target_gene)
                    if target_data is not None:
                        if f"{target_gene}_bar_settings" not in st.session_state:
                            st.session_state[f"{target_gene}_bar_settings"] = {}
                        for _, row in target_data.iterrows():
                            condition = row["Condition"]
                            group = row.get("Group", "Treatment")
                            bar_key = f"{target_gene}_{condition}"
                            color = preset_colors.get(group, preset_colors.get("Treatment", "#D3D3D3"))
                            if bar_key not in st.session_state[f"{target_gene}_bar_settings"]:
                                st.session_state[f"{target_gene}_bar_settings"][bar_key] = {
                                    "color": color, "show_sig": True, "show_err": True,
                                    "show_sig_1": True, "show_sig_2": True, "show_sig_3": True,
                                }
                            else:
                                st.session_state[f"{target_gene}_bar_settings"][bar_key]["color"] = color
                            st.session_state.graph_settings["bar_colors_per_sample"][bar_key] = color
            st.success(f"Settings applied to all {len(gene_list)} genes.")
            st.rerun()
```

- [ ] **Step 5: Wire new settings into graph call**

Update the `create_gene_graph` call (~line 5798) to pass new params:

```python
# Compute replicate data if data points enabled
replicate_df = None
show_dp = st.session_state.graph_settings.get(f"{current_gene}_show_data_points", False)
if show_dp:
    raw = st.session_state.get("data")
    hk = st.session_state.get("hk_gene")
    ref = st.session_state.get("analysis_ref_condition")
    mapping = st.session_state.get("sample_mapping", {})
    excl = st.session_state.get("excluded_wells", set())
    if raw is not None and hk and ref:
        from qpcr.analysis import AnalysisEngine
        all_replicates = AnalysisEngine.compute_replicate_fold_changes(
            raw, hk, ref, mapping, excl,
        )
        replicate_df = all_replicates[all_replicates["Target"] == current_gene]

current_settings["sig_style"] = st.session_state.graph_settings.get(
    f"{current_gene}_sig_style", "direct"
)

fig = GraphGenerator.create_gene_graph(
    gene_data, current_gene, current_settings, efficacy_config,
    sample_order=st.session_state.get("sample_order"),
    condition_colors=st.session_state.get("condition_colors", {}),
    display_gene_name=display_gene_name,
    ref_line_value=ref_line_val, ref_line_label=ref_line_lbl,
    show_data_points=show_dp, replicate_data=replicate_df,
)
```

- [ ] **Step 6: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 7: Commit**

```bash
git add "streamlit qpcr analysis v1.py"
git commit -m "feat: add preset dropdown, sig style radio, data points toggle, apply-to-all button"
```

---

### Task 8: Quick-Size Presets + Cleaner Per-Bar Table

**Files:**
- Modify: `streamlit qpcr analysis v1.py:5584-5758` (Settings expander)

- [ ] **Step 1: Add figure size preset above sliders**

In the Settings expander, in `s_col2` (~line 5610), add before the width/height sliders:

```python
with s_col2:
    st.caption("**Figure Size**")
    size_preset_key = f"{current_gene}_size_preset"
    size_preset_names = list(FIGURE_SIZE_PRESETS.keys()) + ["Custom"]
    if size_preset_key not in st.session_state.graph_settings:
        st.session_state.graph_settings[size_preset_key] = "PPT Full"
    current_size_preset = st.session_state.graph_settings.get(size_preset_key, "PPT Full")
    if current_size_preset not in size_preset_names:
        current_size_preset = "Custom"
    selected_size = st.selectbox(
        "Size Preset",
        size_preset_names,
        index=size_preset_names.index(current_size_preset),
        key=f"size_preset_{current_gene}",
    )
    if selected_size != "Custom" and selected_size in FIGURE_SIZE_PRESETS:
        preset = FIGURE_SIZE_PRESETS[selected_size]
        st.session_state.graph_settings[f"{current_gene}_figure_width"] = preset["width"]
        st.session_state.graph_settings[f"{current_gene}_figure_height"] = preset["height"]
    st.session_state.graph_settings[size_preset_key] = selected_size

    fig_width_cm = st.slider(
        "Width (cm)", 10.0, 40.0,
        value=float(st.session_state.graph_settings.get(f"{current_gene}_figure_width", 28)),
        step=0.5, key=f"fig_w_{current_gene}",
    )
    fig_height_cm = st.slider(
        "Height (cm)", 6.0, 25.0,
        value=float(st.session_state.graph_settings.get(f"{current_gene}_figure_height", 16)),
        step=0.5, key=f"fig_h_{current_gene}",
    )
    # Check if manual slider change → switch to Custom
    if selected_size != "Custom" and selected_size in FIGURE_SIZE_PRESETS:
        preset = FIGURE_SIZE_PRESETS[selected_size]
        if fig_width_cm != preset["width"] or fig_height_cm != preset["height"]:
            st.session_state.graph_settings[size_preset_key] = "Custom"
    st.session_state.graph_settings[f"{current_gene}_figure_width"] = fig_width_cm
    st.session_state.graph_settings[f"{current_gene}_figure_height"] = fig_height_cm
```

- [ ] **Step 2: Refactor per-bar settings table to 3 columns**

Replace the 6-column per-bar table (lines 5718-5760) with a cleaner 3-column layout:

```python
st.markdown("**Per-Bar Settings**")

hdr = st.columns([3, 0.8, 2.5])
hdr[0].markdown("<small>**Condition**</small>", unsafe_allow_html=True)
hdr[1].markdown("<small>**Color**</small>", unsafe_allow_html=True)
hdr[2].markdown("<small>**Options**</small>", unsafe_allow_html=True)

option_labels = ["✱", "#", "†", "±"]
option_keys = ["show_sig_1", "show_sig_2", "show_sig_3", "show_err"]

for idx, (_, row) in enumerate(_bar_display_data.iterrows()):
    condition = row["Condition"]
    group = row.get("Group", "Treatment")
    bar_key = f"{current_gene}_{condition}"
    bs = st.session_state[f"{current_gene}_bar_settings"][bar_key]

    rc = st.columns([3, 0.8, 2.5])
    lbl = condition if len(condition) <= 22 else condition[:19] + "..."
    rc[0].markdown(
        f"<small>{lbl} <span style='color:#888;'>({group})</span></small>",
        unsafe_allow_html=True,
    )
    new_color = rc[1].color_picker(
        "c", bs["color"],
        key=f"cp_{current_gene}_{idx}",
        label_visibility="collapsed",
    )
    bs["color"] = new_color
    st.session_state.graph_settings["bar_colors_per_sample"][bar_key] = new_color

    # Detect if color was manually changed → switch preset to Custom
    color_preset_key = f"{current_gene}_color_preset"
    current_preset = st.session_state.graph_settings.get(color_preset_key, "Classic")
    if current_preset != "Custom":
        expected_color = GRAPH_PRESETS.get(current_preset, {}).get(group, "#D3D3D3")
        if new_color != expected_color:
            st.session_state.graph_settings[color_preset_key] = "Custom"

    # Options as inline checkboxes in a sub-row
    with rc[2]:
        opt_cols = st.columns(4)
        for oi, (ol, ok) in enumerate(zip(option_labels, option_keys)):
            bs[ok] = opt_cols[oi].checkbox(
                ol, bs.get(ok, True),
                key=f"{ok}_{current_gene}_{idx}",
                label_visibility="collapsed",
            )
    bs["show_sig"] = bs["show_sig_1"] or bs["show_sig_2"] or bs["show_sig_3"]
```

- [ ] **Step 3: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 4: Commit**

```bash
git add "streamlit qpcr analysis v1.py"
git commit -m "feat: add size presets and cleaner per-bar settings table"
```

---

### Task 9: Unified Export Button

**Files:**
- Modify: `streamlit qpcr analysis v1.py:5884+` (Export tab)

- [ ] **Step 1: Add "Export All (ZIP)" button at top of Export tab**

After `st.header("Export Results")` (~line 5886), add:

```python
if st.session_state.processed_data and st.session_state.get("graphs"):
    if st.button("📦 Export All (ZIP)", key="export_all_zip", type="primary", use_container_width=True):
        import zipfile
        import io
        from datetime import datetime

        zip_buffer = io.BytesIO()
        errors = []
        efficacy = st.session_state.selected_efficacy

        with st.spinner("Generating complete export..."):
            progress = st.progress(0, text="Starting...")

            with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
                # 1. Excel report
                progress.progress(10, text="Generating Excel...")
                try:
                    from qpcr.export import export_to_excel
                    excel_buf = export_to_excel(
                        st.session_state.processed_data,
                        st.session_state.get("analysis_params", {}),
                        st.session_state.get("sample_mapping", {}),
                        st.session_state.get("data"),
                    )
                    zf.writestr("qPCR_Report.xlsx", excel_buf.getvalue())
                except Exception as e:
                    errors.append(f"Excel: {e}")

                # 2. PNG figures
                for i, (gene, fig) in enumerate(st.session_state.graphs.items()):
                    progress.progress(
                        10 + int(60 * (i + 1) / len(st.session_state.graphs)),
                        text=f"Exporting {gene}...",
                    )
                    try:
                        img_bytes = fig.to_image(format="png", scale=2, width=1200, height=800)
                        zf.writestr(f"figures/{gene}.png", img_bytes)
                    except Exception as e:
                        errors.append(f"PNG {gene}: {e}")

                # 3. PPT report
                progress.progress(80, text="Generating PowerPoint...")
                try:
                    from qpcr.report import PPTGenerator
                    ppt_gen = PPTGenerator()
                    ppt_buf = ppt_gen.generate_presentation(
                        processed_data=st.session_state.processed_data,
                        graphs=st.session_state.graphs,
                        efficacy_type=efficacy,
                        hk_gene=st.session_state.get("hk_gene", ""),
                        ref_sample=st.session_state.get("analysis_ref_condition", ""),
                        compare_condition=st.session_state.get("analysis_cmp_condition", ""),
                        graph_settings=st.session_state.get("graph_settings", {}),
                    )
                    zf.writestr("qPCR_Report.pptx", ppt_buf.getvalue())
                except Exception as e:
                    errors.append(f"PPT: {e}")

                # 4. Interactive HTML
                progress.progress(90, text="Generating HTML...")
                try:
                    html_parts = ["<html><head><title>qPCR Graphs</title></head><body>"]
                    for gene, fig in st.session_state.graphs.items():
                        html_parts.append(f"<h2>{gene}</h2>")
                        html_parts.append(fig.to_html(full_html=False, include_plotlyjs="cdn"))
                    html_parts.append("</body></html>")
                    zf.writestr("figures_html/all_graphs.html", "\n".join(html_parts))
                except Exception as e:
                    errors.append(f"HTML: {e}")

                # 5. Error manifest (if any)
                if errors:
                    zf.writestr("_errors.txt", "\n".join(errors))

            progress.progress(100, text="Done!")

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        st.download_button(
            label=f"⬇️ Download ZIP ({len(st.session_state.graphs)} genes)",
            data=zip_buffer.getvalue(),
            file_name=f"qPCR_Export_{efficacy}_{timestamp}.zip",
            mime="application/zip",
            key="dl_all_zip",
        )
        if errors:
            st.warning(f"Some exports failed: {', '.join(errors)}")

    st.markdown("---")
```

- [ ] **Step 2: Test manually** (Streamlit UI — no automated test for this)

Run: `streamlit run "streamlit qpcr analysis v1.py"` and verify the Export tab shows the ZIP button.

- [ ] **Step 3: Commit**

```bash
git add "streamlit qpcr analysis v1.py"
git commit -m "feat: add unified Export All (ZIP) button"
```

---

### Task 10: Final Integration + Full Test Run

**Files:**
- All modified files
- Test: full suite

- [ ] **Step 1: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 2: Run the Streamlit app and verify visually**

Run: `streamlit run "streamlit qpcr analysis v1.py"`

Verify:
1. Color preset dropdown works and applies colors
2. Significance brackets display correctly in bracketed mode
3. Data points toggle shows/hides jittered dots
4. Apply to All Genes copies settings
5. Size presets populate width/height sliders
6. Per-bar table is cleaner (3 columns)
7. Export All ZIP generates a valid ZIP file
8. PPT export includes new graph features

- [ ] **Step 3: Final commit**

```bash
git add qpcr/ tests/ "streamlit qpcr analysis v1.py"
git commit -m "feat: graph improvements — presets, brackets, data points, UI polish"
```

**Note:** Do NOT use `git add -A` — the repo has untracked files (screenshots, .claude/) that should not be committed.
