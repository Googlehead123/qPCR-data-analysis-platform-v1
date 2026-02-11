# Five Improvements to qPCR Analysis App

## TL;DR

> **Quick Summary**: Five targeted improvements to the Streamlit qPCR app: triplicate health status grid + UI polish, PPT export image fix, graph "undefined" title fix, axis spacing + color fix, and redundant tab removal.
> 
> **Deliverables**:
> - REQ1: Interactive triplicate health grid with real-time status, improved well browser UI
> - REQ2: PPT export images no longer cut off
> - REQ3: No "undefined" text on graphs
> - REQ4: Proper axis label spacing + correct positive control color
> - REQ5: Redundant PPT Report tab removed
> 
> **Estimated Effort**: Medium
> **Parallel Execution**: YES - 3 waves
> **Critical Path**: Tasks 2-5 are independent quick fixes; Task 1 is the main effort

---

## Context

### Original Request
Five improvements to `streamlit qpcr analysis v1.py` (~5600 lines). All changes are in a single file. Root causes for REQ2-5 are already identified with exact line numbers. REQ1 requires integrating the existing but unused `render_triplicate_grid()` function and enhancing the QC tab UI.

### Key Constraints
- NO test file modifications
- NO new dependencies
- NO module extraction (single-file architecture)
- All existing QC functionality must remain intact
- Tab5 PPT export must stay fully functional after removing Tab6
- File has spaces in name: `"streamlit qpcr analysis v1.py"`

---

## Work Objectives

### Core Objective
Apply five targeted improvements to improve UX and fix bugs in the qPCR analysis app.

### Definition of Done
- [ ] `pytest tests/ --ignore=tests/test_ppt_report.py -q` → 79 tests pass
- [ ] No "undefined" text visible on graphs
- [ ] PPT export images include full significance labels (not cut off)
- [ ] Axis labels have spacing from axes
- [ ] Positive Control bars are `#909090`
- [ ] Only 6 tabs remain (no "📑 PPT Report")
- [ ] Triplicate health grid renders in QC tab with color-coded status

### Must NOT Have
- New pip dependencies
- Changes to test files
- Multi-file extraction
- Breaking changes to existing QC functions
- Removal of `ReportGenerator._fig_to_image()` (used by PPTGenerator)

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (pytest)
- **User wants tests**: NO modifications to test files allowed
- **QA approach**: Automated test suite + manual verification

### Verification
```bash
pytest tests/ --ignore=tests/test_ppt_report.py -q
# Expected: 79 passed
```

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (All independent - start immediately):
├── Task 1: REQ3 - Fix "undefined" graph title
├── Task 2: REQ4 - Axis spacing + positive control color
├── Task 3: REQ5 - Remove redundant PPT Report tab
└── Task 4: REQ2 - Fix PPT export image cutoff

Wave 2 (After Wave 1 - depends on nothing but is largest):
└── Task 5: REQ1 - Triplicate health grid + UI improvements

Wave 3 (Final):
└── Task 6: Run full test suite, verify all changes
```

### Dependency Matrix

| Task | Depends On | Blocks | Can Parallelize With |
|------|------------|--------|---------------------|
| 1 | None | 6 | 2, 3, 4 |
| 2 | None | 6 | 1, 3, 4 |
| 3 | None | 6 | 1, 2, 4 |
| 4 | None | 6 | 1, 2, 3 |
| 5 | None | 6 | 1, 2, 3, 4 |
| 6 | 1-5 | None | None (final) |

---

## TODOs

- [ ] 1. Fix "undefined" graph title text (REQ3)

  **What to do**:
  - At line 2132, change `title=None` to `title=""` in `fig.update_layout()`
  - This prevents Plotly.js from serializing `None` → `null` → "undefined"

  **Must NOT do**:
  - Don't change any other layout properties in this edit

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 2, 3, 4)
  - **Blocks**: Task 6
  - **Blocked By**: None

  **References**:
  - `streamlit qpcr analysis v1.py:2131-2132` - The `fig.update_layout(title=None, ...)` call in `GraphGenerator.create_gene_graph()`

  **Acceptance Criteria**:
  - [ ] Line 2132: `title=""` instead of `title=None`
  - [ ] `pytest tests/ --ignore=tests/test_ppt_report.py -q` → 79 passed

  **Commit**: YES
  - Message: `fix(graphs): replace title=None with empty string to prevent "undefined" text`
  - Files: `streamlit qpcr analysis v1.py`

---

- [ ] 2. Fix axis label spacing and positive control color (REQ4)

  **What to do**:
  - In `DEFAULT_GROUP_COLORS` (line 30-38): Change `"Positive Control": "#333333"` to `"#909090"` and `"Inducer": "#333333"` to `"#909090"`
  - In y-axis config (line 2059-2074): Add `standoff=10` to the `title=dict(...)` inside `y_axis_config`
  - In x-axis config (line 2133-2146): Add `ticklabelstandoff=10` or equivalent padding parameter

  **Must NOT do**:
  - Don't change graph dimensions or other color mappings

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 3, 4)
  - **Blocks**: Task 6
  - **Blocked By**: None

  **References**:
  - `streamlit qpcr analysis v1.py:30-38` - `DEFAULT_GROUP_COLORS` dict with `"Positive Control": "#333333"` and `"Inducer": "#333333"`
  - `streamlit qpcr analysis v1.py:2059-2074` - `y_axis_config = dict(title=dict(text=..., font=...), ...)` — add `standoff=10` inside title dict
  - `streamlit qpcr analysis v1.py:2133-2146` - x-axis config in `fig.update_layout()` — add tick padding

  **Acceptance Criteria**:
  - [ ] `DEFAULT_GROUP_COLORS["Positive Control"]` == `"#909090"`
  - [ ] `DEFAULT_GROUP_COLORS["Inducer"]` == `"#909090"`
  - [ ] y-axis title dict contains `standoff=10`
  - [ ] x-axis config has tick padding (e.g., `ticklabelstandoff=10` or `tickson`/`dtick` spacing)
  - [ ] `pytest tests/ --ignore=tests/test_ppt_report.py -q` → 79 passed

  **Commit**: YES
  - Message: `fix(graphs): add axis label spacing and update positive control color to #909090`
  - Files: `streamlit qpcr analysis v1.py`

---

- [ ] 3. Remove redundant PPT Report tab (REQ5)

  **What to do**:
  - Line 2962: Change `tab1, tab_qc, tab2, tab3, tab4, tab5, tab6 = st.tabs([...])` to remove `tab6` and `"📑 PPT Report"` from the list
  - Delete the entire `with tab6:` block (lines 5430-5592)
  - Keep `ReportGenerator` class intact (its `_fig_to_image()` is used by `PPTGenerator`)

  **Must NOT do**:
  - Don't delete the `ReportGenerator` class or any of its methods (tests reference them)
  - Don't modify Tab 5 export functionality

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 2, 4)
  - **Blocks**: Task 6
  - **Blocked By**: None

  **References**:
  - `streamlit qpcr analysis v1.py:2962-2971` - Tab definition with 7 tabs → change to 6
  - `streamlit qpcr analysis v1.py:5430-5592` - Entire `with tab6:` block to delete
  - `streamlit qpcr analysis v1.py:5056-5092` - Tab 5 PPT export (must remain untouched)

  **Acceptance Criteria**:
  - [ ] Tab definition has 6 tabs, no `tab6` variable
  - [ ] No `"📑 PPT Report"` in tab labels
  - [ ] Lines 5430-5592 (`with tab6:` block) removed
  - [ ] `ReportGenerator` class still exists with all methods
  - [ ] `pytest tests/ --ignore=tests/test_ppt_report.py -q` → 79 passed

  **Commit**: YES
  - Message: `refactor(ui): remove redundant PPT Report tab (tab6), export already in tab5`
  - Files: `streamlit qpcr analysis v1.py`

---

- [ ] 4. Fix PPT export image cutoff (REQ2)

  **What to do**:
  - `PPTGenerator.generate_presentation()` (line 2832-2835): Change `margin=dict(l=80, r=80, t=60, b=80)` to `margin=dict(l=80, r=80, t=60, b=180)` and `height=600` to `height=700`
  - `ReportGenerator.create_presentation()` (line 2238-2241): Same fix — `b=80` → `b=180`, `height=600` → `height=700`
  - `PPTGenerator.create_gene_slide()` (line 2721): Change `fig.to_image(format="png", scale=2, width=1200, height=800)` to `height=900` to accommodate bottom margin

  **Must NOT do**:
  - Don't change PNG export margins (those already work with `b=180`)
  - Don't change slide dimensions or image positioning ratios significantly

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 2, 3)
  - **Blocks**: Task 6
  - **Blocked By**: None

  **References**:
  - `streamlit qpcr analysis v1.py:2831-2838` - `PPTGenerator.generate_presentation()` layout with `b=80`
  - `streamlit qpcr analysis v1.py:2237-2243` - `ReportGenerator.create_presentation()` layout with `b=80`
  - `streamlit qpcr analysis v1.py:2721` - `PPTGenerator.create_gene_slide()` `fig.to_image()` with `height=800`

  **Acceptance Criteria**:
  - [ ] All three PPT export paths use `b=180` bottom margin
  - [ ] Heights increased proportionally (600→700 for layout, 800→900 for to_image)
  - [ ] `pytest tests/ --ignore=tests/test_ppt_report.py -q` → 79 passed

  **Commit**: YES
  - Message: `fix(export): increase PPT export bottom margin to prevent significance label cutoff`
  - Files: `streamlit qpcr analysis v1.py`

---

- [ ] 5. Add triplicate health status grid + UI improvements (REQ1)

  **What to do**:
  - **Activate `render_triplicate_grid()`**: Call the existing function (line 1168) from `qc_tab1` to display a color-coded gene×sample grid showing triplicate health status at a glance
  - **Add health summary banner**: At the top of qc_tab1, show an aggregated health summary (e.g., "✅ 45 OK | ⚠️ 3 Warning | ❌ 1 Error") that updates when wells are excluded/included
  - **Integrate grid into well browser**: Place the triplicate grid ABOVE the existing per-gene expanders in qc_tab1 so users see overall status first, then drill down
  - **Make grid reactive**: The grid should use the current `excluded_wells` state so it reflects real-time exclusion changes
  - **UI polish for well browser**: 
    - Add colored status badges to per-sample statistics (green/yellow/red based on health)
    - Add subtle section dividers between genes
    - Use `st.container()` with border for each gene section instead of plain expanders

  **Must NOT do**:
  - Don't remove or break existing well browser functionality (data_editor, Include checkboxes, per-sample stats)
  - Don't modify `QualityControl.get_triplicate_data()` logic
  - Don't add new dependencies
  - Don't create new files

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: [`frontend-ui-ux`]
    - `frontend-ui-ux`: UI layout and visual design decisions for the grid and status badges

  **Parallelization**:
  - **Can Run In Parallel**: YES (independent of other tasks)
  - **Parallel Group**: Wave 1-2 (can start immediately, but is largest task)
  - **Blocks**: Task 6
  - **Blocked By**: None

  **References**:

  **Pattern References**:
  - `streamlit qpcr analysis v1.py:1168-1261` - `render_triplicate_grid()` — existing grid implementation with button-per-cell, emoji status, color coding. Currently NOT called anywhere. This is the base to activate and enhance.
  - `streamlit qpcr analysis v1.py:978-1165` - Grid state management functions: `build_grid_matrix()`, `get_cell_status_color()`, `get_cell_display_text()`, `is_cell_selected()`, `set_selected_cell()`, `get_grid_cell_key()` — all exist and are ready to use
  - `streamlit qpcr analysis v1.py:3201-3358` - Current qc_tab1 well browser: filter controls, per-gene expanders with data_editor, per-sample statistics — this is what gets enhanced

  **API/Type References**:
  - `streamlit qpcr analysis v1.py:443-530` - `QualityControl.get_triplicate_data()` — returns DataFrame with columns: Sample, Target, Mean_CT, SD, n, Min_CT, Max_CT, CT_Values, Wells, CV_pct, Range, Status, Severity
  - `streamlit qpcr analysis v1.py:488-520` - `get_health_status()` inner function — health checks: Low n, High CV, High CT, Low CT, High range, Grubbs outlier; severity: ok/warning/error

  **Color References**:
  - `streamlit qpcr analysis v1.py:1111` - `get_cell_status_color()` returns: green `#d4edda` / yellow `#fff3cd` / red `#f8d7da`

  **Acceptance Criteria**:
  - [ ] Triplicate health grid is visible in qc_tab1 above the well browser
  - [ ] Grid shows gene×sample cells with color coding (green/yellow/red)
  - [ ] Health summary banner shows counts of OK/Warning/Error triplicates
  - [ ] Excluding a well and refreshing updates the grid status
  - [ ] All existing well browser functionality works (data_editor, Include checkboxes, per-sample stats, filter controls, Include/Exclude All, Undo)
  - [ ] `pytest tests/ --ignore=tests/test_ppt_report.py -q` → 79 passed

  **Commit**: YES
  - Message: `feat(qc): activate triplicate health grid with real-time status and UI improvements`
  - Files: `streamlit qpcr analysis v1.py`

---

- [ ] 6. Final verification

  **What to do**:
  - Run full test suite: `pytest tests/ --ignore=tests/test_ppt_report.py -q`
  - Verify 79 tests pass
  - Spot-check the file for any syntax errors introduced

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 3 (final)
  - **Blocks**: None
  - **Blocked By**: Tasks 1-5

  **Acceptance Criteria**:
  - [ ] `pytest tests/ --ignore=tests/test_ppt_report.py -q` → 79 passed
  - [ ] No Python syntax errors in file

  **Commit**: NO

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `fix(graphs): replace title=None with empty string to prevent "undefined" text` | `streamlit qpcr analysis v1.py` | pytest |
| 2 | `fix(graphs): add axis label spacing and update positive control color to #909090` | `streamlit qpcr analysis v1.py` | pytest |
| 3 | `refactor(ui): remove redundant PPT Report tab` | `streamlit qpcr analysis v1.py` | pytest |
| 4 | `fix(export): increase PPT export bottom margin to prevent significance label cutoff` | `streamlit qpcr analysis v1.py` | pytest |
| 5 | `feat(qc): activate triplicate health grid with real-time status and UI improvements` | `streamlit qpcr analysis v1.py` | pytest |

---

## Success Criteria

### Verification Commands
```bash
pytest tests/ --ignore=tests/test_ppt_report.py -q  # Expected: 79 passed
```

### Final Checklist
- [ ] No "undefined" text on graphs (REQ3)
- [ ] Axis labels have spacing, positive control is #909090 (REQ4)
- [ ] Only 6 tabs, no "📑 PPT Report" (REQ5)
- [ ] PPT export images show full significance labels (REQ2)
- [ ] Triplicate health grid visible and reactive in QC tab (REQ1)
- [ ] 79 tests pass
