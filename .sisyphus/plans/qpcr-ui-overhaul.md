# qPCR Platform: Apple UI + Triplicate Form + Per-Gene Heatmaps

## TL;DR

> **Quick Summary**: Overhaul the entire visual design to Apple-inspired monochrome, replace the Triplicate Browser's data_editor with form-based checkboxes, and add per-gene/per-sample heatmaps to QC Overview.
> 
> **Deliverables**:
> - Modern monochrome CSS theme replacing all purple gradients and emoji clutter
> - Form-based well exclusion UX with explicit submit button
> - Per-gene and per-sample plate heatmaps in QC Overview
> - Dead code cleanup (render_triplicate_grid, build_grid_matrix)
> 
> **Estimated Effort**: Large
> **Parallel Execution**: NO - sequential (single file, overlapping line ranges)
> **Critical Path**: Task 1 (CSS) → Task 2 (emoji cleanup) → Task 3 (dead code) → Task 4 (triplicate form) → Task 5 (heatmaps) → Task 6 (verify)

---

## Context

### Original Request
Three changes to the qPCR analysis platform:
1. Modern Apple-like UI (monochrome, clean, minimal emojis)
2. Triplicate Browser: delete Health Status Grid, replace data_editor with st.form + checkboxes
3. Per-gene/per-sample heatmaps in QC Overview

### Research Findings
- `render_triplicate_grid` (line 1168) and `build_grid_matrix` (line 1065) are ONLY used at line 3385 — safe to delete
- `well_editor_{gene}` session state keys are ONLY used as data_editor keys, never read elsewhere — safe to remove
- `create_plate_heatmap(data, value_col="CT", excluded_wells=set())` works with filtered DataFrames
- Additional emoji/CSS locations found beyond user's initial list (see Task 2 references)

---

## Work Objectives

### Core Objective
Transform the visual identity from emoji-heavy purple-gradient to Apple-inspired monochrome, fix Triplicate Browser UX with form submission, and add comprehensive heatmap views.

### Concrete Deliverables
- Updated `streamlit qpcr analysis v1.py` with all 3 changes

### Definition of Done
- [ ] `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"` passes
- [ ] `pytest tests/ -v` all tests pass
- [ ] No purple gradients (#667eea, #764ba2) remain in file
- [ ] No `render_triplicate_grid` or `build_grid_matrix` functions in file
- [ ] No `st.data_editor` in the Triplicate Browser section
- [ ] `st.form` present in Triplicate Browser section

### Must Have
- All existing business logic unchanged (QPCRParser, AnalysisEngine, GraphGenerator, QualityControl classes)
- All session state variables work identically
- All existing functionality preserved
- Korean text labels preserved

### Must NOT Have (Guardrails)
- NO changes to class method signatures or business logic
- NO changes to data flow (parse → map → analyze → graph → export)
- NO new dependencies added
- NO splitting into multiple files
- NO removal of functional UI elements (only visual/UX changes)
- NO leftover purple gradient strings (#667eea, #764ba2) anywhere

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (pytest, 8 test files)
- **User wants tests**: Run existing tests as regression check
- **Framework**: pytest
- **QA approach**: `ast.parse` for syntax + `pytest tests/ -v` for regression

---

## Execution Strategy

### Sequential Execution (Single File)

All tasks modify the same file and have overlapping concerns. Must be sequential.

```
Task 1: Global CSS theme replacement
  ↓
Task 2: Emoji cleanup across all headers/labels
  ↓
Task 3: Delete dead code (render_triplicate_grid, build_grid_matrix)
  ↓
Task 4: Replace Triplicate Browser data_editor with st.form
  ↓
Task 5: Add per-gene/per-sample heatmaps to QC Overview
  ↓
Task 6: Final verification (ast.parse + pytest)
```

### Dependency Matrix

| Task | Depends On | Blocks |
|------|------------|--------|
| 1 | None | 2 |
| 2 | 1 | 3 |
| 3 | 2 | 4 |
| 4 | 3 | 5 |
| 5 | 4 | 6 |
| 6 | 5 | None (final) |

---

## TODOs

- [ ] 1. Replace all CSS injection blocks with unified Apple-inspired theme

  **What to do**:
  - Insert a single comprehensive `<style>` block near the top of the file (after line 27, before CONSTANTS) that defines the entire theme
  - Apple design principles: lots of whitespace, `font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Segoe UI', sans-serif`, monochrome palette (white/near-white backgrounds, #1d1d1f text, #86868b secondary text), one subtle accent (#0071e3 blue or keep a muted indigo), subtle `box-shadow` instead of borders, generous padding/border-radius
  - Replace the QC banner (lines 3124-3132): remove purple gradient div, replace with clean minimal header (no background gradient, just clean typography)
  - Replace graph CSS (lines 4236-4273): update chart shadow, expander bg, toolbar gradient, stat-highlight to monochrome
  - Replace gene pill CSS (lines 4309-4321): change `.gene-pill-active` from purple gradient to monochrome accent, `.gene-pill-inactive` to light grey, `.compact-control` to clean white
  - Replace condition color card CSS (lines 4412-4431): update `.condition-color-card` from #f8f9fa/border to cleaner white/subtle-shadow
  - Replace color editor CSS (lines 4524-4552): update `.color-editor-card`, `.color-editor-title`, `.color-editor-subtitle` to match theme
  - Replace analysis summary card (lines 4093-4102): remove #f0f2f6 background, use clean white card with subtle shadow
  - Replace sample mapping header (lines 3758-3774): update table header styling from #f8f9fa to clean design
  - Replace footer (lines 5220-5227): simplify, remove emojis, clean typography
  - Update `st.set_page_config` title (line 24): change to "qPCR Analysis Suite" (drop "Pro")
  - Address inline HTML styling at lines 3801-3804 (mapping order), 3884 (hr), 4464-4471 (condition card HTML), 4572-4580 (color editor HTML), 4623-4626 (spacer div)

  **Must NOT do**:
  - Change any Plotly graph generation code in GraphGenerator class
  - Change any computation logic
  - Remove functional elements — only restyle them
  - Add external CSS files or CDN links

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: [`frontend-ui-ux`]
    - `frontend-ui-ux`: CSS theming, Apple design system knowledge, visual polish

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential — Task 1
  - **Blocks**: Task 2
  - **Blocked By**: None

  **References**:

  **Pattern References** (current CSS to replace):
  - `streamlit qpcr analysis v1.py:3124-3132` - QC banner with purple gradient — replace with minimal header
  - `streamlit qpcr analysis v1.py:4236-4273` - Graph tab CSS (chart shadow, expander, toolbar, stat-highlight) — retheme
  - `streamlit qpcr analysis v1.py:4309-4321` - Gene pill CSS with purple gradient — change to monochrome accent
  - `streamlit qpcr analysis v1.py:4412-4431` - Condition color card CSS — update to white/shadow
  - `streamlit qpcr analysis v1.py:4524-4552` - Color editor card CSS — update to match theme
  - `streamlit qpcr analysis v1.py:4093-4102` - Analysis summary card HTML — simplify
  - `streamlit qpcr analysis v1.py:3758-3774` - Sample mapping header table — restyle
  - `streamlit qpcr analysis v1.py:5220-5227` - Footer HTML — simplify
  - `streamlit qpcr analysis v1.py:23-27` - Page config

  **Inline HTML to update**:
  - `streamlit qpcr analysis v1.py:3801-3804` - Mapping order number div
  - `streamlit qpcr analysis v1.py:3884` - Styled hr element
  - `streamlit qpcr analysis v1.py:4464-4471` - Condition color card HTML instance
  - `streamlit qpcr analysis v1.py:4572-4580` - Color editor subtitle HTML
  - `streamlit qpcr analysis v1.py:4623-4626` - Spacer div

  **Acceptance Criteria**:
  - [ ] `grep -c "667eea\|764ba2" "streamlit qpcr analysis v1.py"` returns 0
  - [ ] `grep -c "linear-gradient(135deg, #667eea" "streamlit qpcr analysis v1.py"` returns 0
  - [ ] New unified `<style>` block exists near top of file with Apple-inspired CSS
  - [ ] `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"` exits 0

  **Commit**: YES
  - Message: `style(ui): replace purple gradient theme with Apple-inspired monochrome design`
  - Files: `streamlit qpcr analysis v1.py`
  - Pre-commit: `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"`

---

- [ ] 2. Remove emojis from all headers, labels, sidebar, and tab names

  **What to do**:
  - Clean up `st.title` (line 2916): remove 🧬, just "qPCR Analysis Suite"
  - Clean up sidebar (lines 2920-2963):
    - Line 2921: `st.header("💬 Quick Guide")` → `st.header("Quick Guide")`
    - Lines 2922-2929: Remove all 📁🔍🗺️🔬📊📤 from guide steps, keep text
    - Line 2932: Remove ⚡ from "Quick Actions"
    - Line 2935: Remove ✅ from success message
    - Line 2940: Remove ⚠️ from warning
    - Line 2943: Remove ✅ from success
    - Line 2945: Remove 📁 from info
    - Line 2948: Remove ⌨️ from "Navigation Tips"
  - Clean up tab labels (lines 2966-2974): Remove all emojis from tab names — use plain text: "Upload", "QC Check", "Mapping", "Analysis", "Graphs", "Export"
  - Clean up section headers throughout:
    - Line 3070: Remove 📊 from "Data Preview"
    - Line 3073: Remove ⚠️ from "Data Validation"  
    - Line 3384: Remove 🔍 from "Health Status Grid" (this block gets deleted in Task 3, but clean label if remnants)
    - Line 3482: Remove 🧪 from "Plate Heatmap"
    - Line 3513: Remove ⚠️ from "Flagged Wells"
    - Line 3524: Remove 🚫 from "Exclude All Flagged" button
    - Line 3547: Remove ✅ from success message
    - Line 3552: Remove 📋 from "Pre-Analysis Summary"
    - Line 3633: Remove ⚠️ from warning
    - Line 3733: Remove 🗺️ from "Sample Condition Mapping"
    - Line 3898: Remove 📊 from "Mapping Summary"
    - Line 3950: Remove 🔬 from "Run Full Analysis"
    - Line 3966: Remove 📊⚙️ from "Analysis Configuration"
    - Line 4031: Remove ⚙️ from "Statistical Options"
    - Line 4152: Remove 📊 from "Analysis Summary"
    - Line 4166: Remove 🧬 from "Gene-by-Gene Results"
    - Line 4234: Remove emojis from "Individual Gene Graphs" header
    - Line 4323: Remove 🧬 from "Select Gene"
    - Line 4411: Remove 🎨 from "Condition Colors" expander
    - Line 4524: Remove 🎨 from "Bar Color & Visibility Editor" expander
    - Lines 4772-5028: Remove emojis from all Export tab headers (📦📑📋📸)
  - Keep functional emojis that serve as status indicators in `st.success`, `st.warning`, `st.error` calls ONLY where Streamlit itself doesn't already provide icons (Streamlit adds its own icons to these widgets)

  **Must NOT do**:
  - Remove Korean text labels
  - Change any text content beyond emoji removal
  - Remove emojis from data values or EFFICACY_CONFIG

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential — Task 2
  - **Blocks**: Task 3
  - **Blocked By**: Task 1

  **References**:

  **Pattern References** (all emoji locations):
  - `streamlit qpcr analysis v1.py:2916` - st.title with 🧬
  - `streamlit qpcr analysis v1.py:2920-2963` - Sidebar with 💬⚡✅⚠️📁⌨️
  - `streamlit qpcr analysis v1.py:2966-2974` - Tab labels with 📁🔍🗺️🔬📊📤
  - `streamlit qpcr analysis v1.py:3070,3073` - Upload tab headers
  - `streamlit qpcr analysis v1.py:3482,3513,3547,3552` - QC Overview headers
  - `streamlit qpcr analysis v1.py:3733,3898,3950,3966` - Mapping/Analysis headers
  - `streamlit qpcr analysis v1.py:4031,4152,4166,4234,4323` - Analysis/Graphs headers
  - `streamlit qpcr analysis v1.py:4411,4524` - Color editor expanders
  - `streamlit qpcr analysis v1.py:4772-5028` - Export tab headers

  **Acceptance Criteria**:
  - [ ] Tab labels contain no emoji characters (verify: `grep "st.tabs" "streamlit qpcr analysis v1.py"` shows plain text)
  - [ ] `st.title` contains no emoji
  - [ ] Sidebar `st.header` and `st.markdown` contain no emoji
  - [ ] `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"` exits 0

  **Commit**: YES
  - Message: `style(ui): remove emojis from headers, tabs, sidebar for clean Apple aesthetic`
  - Files: `streamlit qpcr analysis v1.py`
  - Pre-commit: `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"`

---

- [ ] 3. Delete dead code: Health Status Grid, render_triplicate_grid, build_grid_matrix

  **What to do**:
  - Delete `build_grid_matrix` function (lines 1065-1166)
  - Delete `render_triplicate_grid` function (lines 1168-1262)
  - Delete the entire Health Status Grid block (lines 3342-3389):
    - This includes: `_tri_data` computation, severity counting, banner, `_gene_severity` lookup, the expander with `render_triplicate_grid`, and the `st.markdown("---")`
  - After deletion, the `_gene_severity` dict no longer exists. Fix the gene label (line 3403-3404):
    - Remove `gene_status = _gene_severity.get(gene, "")` (line 3403)
    - Simplify `gene_label` (line 3404) to remove `{gene_status}` prefix:
      ```python
      gene_label = f"{gene}  ({len(gene_wells)} wells, {gene_excluded} excluded)" if gene_excluded > 0 else f"{gene}  ({len(gene_wells)} wells)"
      ```
  - Set `_gene_severity = {}` default (line 3389) is also deleted since the whole block goes away

  **Must NOT do**:
  - Delete any code that IS referenced elsewhere
  - Change the per-gene expander sections (lines 3391-3474) other than the gene_label fix
  - Remove the per-sample statistics display (lines 3456-3474)

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential — Task 3
  - **Blocks**: Task 4
  - **Blocked By**: Task 2

  **References**:

  **Code to delete**:
  - `streamlit qpcr analysis v1.py:1065-1166` - `build_grid_matrix` function (only used by render_triplicate_grid)
  - `streamlit qpcr analysis v1.py:1168-1262` - `render_triplicate_grid` function (only used at line 3385)
  - `streamlit qpcr analysis v1.py:3342-3389` - Health Status Grid block (calls render_triplicate_grid)

  **Code to modify**:
  - `streamlit qpcr analysis v1.py:3403-3404` - gene_label construction references `_gene_severity` — must remove reference

  **Acceptance Criteria**:
  - [ ] `grep -c "render_triplicate_grid\|build_grid_matrix" "streamlit qpcr analysis v1.py"` returns 0
  - [ ] `grep -c "_gene_severity" "streamlit qpcr analysis v1.py"` returns 0
  - [ ] `grep -c "Health Status Grid" "streamlit qpcr analysis v1.py"` returns 0
  - [ ] `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"` exits 0

  **Commit**: YES
  - Message: `refactor(qc): remove Health Status Grid and dead code (render_triplicate_grid, build_grid_matrix)`
  - Files: `streamlit qpcr analysis v1.py`
  - Pre-commit: `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"`

---

- [ ] 4. Replace Triplicate Browser data_editor with st.form + checkboxes

  **What to do**:
  - In the per-gene expander loop (previously lines ~3406-3454, line numbers shifted after Task 3 deletions), replace the `st.data_editor` block with an `st.form` approach:
  
  **New implementation pattern** (inside each gene expander):
  ```python
  with st.form(key=f"well_form_{gene}"):
      # Display wells as rows with checkboxes
      form_cols_header = st.columns([1, 1, 2, 1, 1])
      form_cols_header[0].markdown("**Include**")
      form_cols_header[1].markdown("**Well**")
      form_cols_header[2].markdown("**Sample**")
      form_cols_header[3].markdown("**CT**")
      form_cols_header[4].markdown("**Dev**")
      
      checkbox_states = {}
      for idx, row in gene_editor_df.iterrows():
          cols = st.columns([1, 1, 2, 1, 1])
          is_included = not is_well_excluded(row["Well"], gene, row["Sample"])
          checkbox_states[idx] = cols[0].checkbox(
              "incl", value=is_included, key=f"cb_{gene}_{row['Well']}", label_visibility="collapsed"
          )
          cols[1].text(row["Well"])
          cols[2].text(row["Sample"])
          cols[3].text(f"{row['CT']:.2f}")
          cols[4].text(f"{row['Deviation']:.3f}")
      
      submitted = st.form_submit_button("Apply Changes", use_container_width=True)
      
      if submitted:
          for idx, row in gene_editor_df.iterrows():
              well = row["Well"]
              sample = row["Sample"]
              include = checkbox_states[idx]
              if not include and not is_well_excluded(well, gene, sample):
                  st.session_state.excluded_wells_history.append(
                      {k: v.copy() for k, v in st.session_state.excluded_wells.items()}
                  )
                  exclude_well(well, gene, sample)
              elif include and is_well_excluded(well, gene, sample):
                  st.session_state.excluded_wells_history.append(
                      {k: v.copy() for k, v in st.session_state.excluded_wells.items()}
                  )
                  include_well(well, gene, sample)
          st.rerun()
  ```
  
  - Keep the gene_editor_df construction (computing Deviation, CT rounding) — it's still needed for display data
  - Keep the per-sample statistics display below the form
  - The per-sample statistics should use the CURRENT session state (not form state) so they update after submit+rerun

  **Must NOT do**:
  - Change `exclude_well()` or `include_well()` function signatures
  - Change `excluded_wells_history` append pattern
  - Change `is_well_excluded()` logic
  - Remove per-sample statistics (lines ~3456-3474)

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: [`frontend-ui-ux`]
    - `frontend-ui-ux`: Form layout, UX flow, checkbox alignment

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential — Task 4
  - **Blocks**: Task 5
  - **Blocked By**: Task 3

  **References**:

  **Pattern References** (current code to replace):
  - `streamlit qpcr analysis v1.py:3407-3454` (post-Task-3 shifted lines) - Current data_editor + change processing block
  - `streamlit qpcr analysis v1.py:3409-3411` - Include column computation pattern (reuse for checkbox default)
  - `streamlit qpcr analysis v1.py:3413-3416` - Deviation computation (keep as-is for display)
  - `streamlit qpcr analysis v1.py:3440-3454` - Exclude/include_well logic (preserve exactly)

  **API References**:
  - `streamlit qpcr analysis v1.py:3445-3453` - `exclude_well` / `include_well` call pattern with history
  - `streamlit qpcr analysis v1.py:3398-3400` - `is_well_excluded(well, gene, sample)` usage pattern

  **Acceptance Criteria**:
  - [ ] `grep -c "st.data_editor" "streamlit qpcr analysis v1.py"` returns 0 (no data_editor in triplicate section — note: check there's no other data_editor elsewhere that should remain)
  - [ ] `grep -c "st.form_submit_button" "streamlit qpcr analysis v1.py"` returns at least 1
  - [ ] `grep -c "well_form_" "streamlit qpcr analysis v1.py"` returns at least 1
  - [ ] `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"` exits 0

  **Commit**: YES
  - Message: `feat(qc): replace data_editor with form-based checkboxes for well exclusion`
  - Files: `streamlit qpcr analysis v1.py`
  - Pre-commit: `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"`

---

- [ ] 5. Add per-gene and per-sample heatmaps to QC Overview

  **What to do**:
  - After the existing overall plate heatmap section (around line 3492, shifted after prior tasks), add two new subsections:
  
  **Per-Gene Heatmaps**:
  ```python
  st.markdown("### Per-Gene Plate Heatmaps")
  all_genes_heatmap = sorted(data["Target"].unique())
  selected_heatmap_gene = st.selectbox(
      "Select gene", all_genes_heatmap, key="heatmap_gene_select"
  )
  gene_filtered = data[data["Target"] == selected_heatmap_gene]
  if not gene_filtered.empty:
      gene_hm_fig = QualityControl.create_plate_heatmap(
          gene_filtered, value_col="CT", excluded_wells=get_all_excluded_wells()
      )
      st.plotly_chart(gene_hm_fig, use_container_width=True)
      st.caption(f"Showing CT values for {selected_heatmap_gene} only")
  ```
  
  **Per-Sample Heatmaps**:
  ```python
  st.markdown("### Per-Sample Plate Heatmaps")
  all_samples_heatmap = sorted(data["Sample"].unique(), key=natural_sort_key)
  selected_heatmap_sample = st.selectbox(
      "Select sample", all_samples_heatmap, key="heatmap_sample_select"
  )
  sample_filtered = data[data["Sample"] == selected_heatmap_sample]
  if not sample_filtered.empty:
      sample_hm_fig = QualityControl.create_plate_heatmap(
          sample_filtered, value_col="CT", excluded_wells=get_all_excluded_wells()
      )
      st.plotly_chart(sample_hm_fig, use_container_width=True)
      st.caption(f"Showing CT values for sample {selected_heatmap_sample} only")
  ```
  
  - Insert these AFTER the existing heatmap (after the caption at line ~3492) and BEFORE the `st.markdown("---")` separator (line ~3510)
  - Place both new subsections inside `heatmap_col1` (the wider column) so layout is consistent
  - Use `st.selectbox` (not tabs/expanders) for gene/sample selection — clean, minimal

  **Must NOT do**:
  - Remove the existing overall plate heatmap
  - Modify `QualityControl.create_plate_heatmap` method
  - Change the Replicate Statistics panel in `heatmap_col2`

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: [`frontend-ui-ux`]
    - `frontend-ui-ux`: Layout decisions for heatmap placement

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential — Task 5
  - **Blocks**: Task 6
  - **Blocked By**: Task 4

  **References**:

  **Pattern References** (existing heatmap code to follow):
  - `streamlit qpcr analysis v1.py:3485-3492` - Existing `create_plate_heatmap` call pattern — follow exactly
  - `streamlit qpcr analysis v1.py:899-930` - `create_plate_heatmap` method signature and how it handles filtered data (iterates wells, maps to 8x12 grid)

  **API References**:
  - `streamlit qpcr analysis v1.py:900` - `create_plate_heatmap(data, value_col="CT", excluded_wells=set())` — exact signature
  - `streamlit qpcr analysis v1.py:3487` - `get_all_excluded_wells()` — helper to get excluded wells set

  **Data shape reference**:
  - DataFrame columns: Well, Sample, Target, CT (minimum needed for heatmap)
  - `natural_sort_key` function available for sample sorting

  **Acceptance Criteria**:
  - [ ] `grep -c "heatmap_gene_select\|heatmap_sample_select" "streamlit qpcr analysis v1.py"` returns 2
  - [ ] `grep -c "Per-Gene Plate Heatmaps\|Per-Sample Plate Heatmaps" "streamlit qpcr analysis v1.py"` returns 2
  - [ ] `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"` exits 0

  **Commit**: YES
  - Message: `feat(qc): add per-gene and per-sample plate heatmaps in QC Overview`
  - Files: `streamlit qpcr analysis v1.py`
  - Pre-commit: `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"`

---

- [ ] 6. Final verification: syntax check + full test suite

  **What to do**:
  - Run syntax validation: `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"`
  - Run full test suite: `pytest tests/ -v`
  - If any test fails, investigate and fix (likely import-level issues if function signatures changed, which they shouldn't have)
  - Verify no purple gradient strings remain: `grep "667eea\|764ba2" "streamlit qpcr analysis v1.py"`
  - Verify no dead code remains: `grep "render_triplicate_grid\|build_grid_matrix\|_gene_severity" "streamlit qpcr analysis v1.py"`

  **Must NOT do**:
  - Skip any failing test — all must pass
  - Modify test files to make tests pass (fix the source instead)

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential — Task 6 (final)
  - **Blocks**: None
  - **Blocked By**: Task 5

  **References**:
  - `tests/test_parser.py` - Parser tests
  - `tests/test_analysis.py` - Analysis engine tests
  - `tests/test_quality_control.py` - QC tests (most likely to be affected)
  - `tests/test_qc_grid.py` - Grid tests (may reference deleted functions!)
  - `tests/test_utils.py` - Utility tests
  - `tests/test_graph.py` - Graph tests
  - `tests/test_ppt_report.py` - Report tests

  **IMPORTANT**: `tests/test_qc_grid.py` may import `render_triplicate_grid` or `build_grid_matrix`. If so, those test imports/tests must be removed or updated since the functions are deleted in Task 3.

  **Acceptance Criteria**:
  - [ ] `python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"` exits 0
  - [ ] `pytest tests/ -v` exits 0 with all tests passing
  - [ ] `grep -c "667eea\|764ba2" "streamlit qpcr analysis v1.py"` returns 0
  - [ ] `grep -c "render_triplicate_grid\|build_grid_matrix\|_gene_severity" "streamlit qpcr analysis v1.py"` returns 0

  **Commit**: NO (verification only — if fixes needed, amend prior commits or create fix commit)

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `style(ui): replace purple gradient theme with Apple-inspired monochrome design` | streamlit qpcr analysis v1.py | ast.parse |
| 2 | `style(ui): remove emojis from headers, tabs, sidebar for clean Apple aesthetic` | streamlit qpcr analysis v1.py | ast.parse |
| 3 | `refactor(qc): remove Health Status Grid and dead code` | streamlit qpcr analysis v1.py | ast.parse |
| 4 | `feat(qc): replace data_editor with form-based checkboxes for well exclusion` | streamlit qpcr analysis v1.py | ast.parse |
| 5 | `feat(qc): add per-gene and per-sample plate heatmaps in QC Overview` | streamlit qpcr analysis v1.py | ast.parse |
| 6 | — | — | pytest tests/ -v |

---

## Success Criteria

### Verification Commands
```bash
python3 -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"  # Expected: no output (success)
pytest tests/ -v  # Expected: all tests pass
grep -c "667eea\|764ba2" "streamlit qpcr analysis v1.py"  # Expected: 0
grep -c "render_triplicate_grid\|build_grid_matrix" "streamlit qpcr analysis v1.py"  # Expected: 0
grep -c "st.data_editor" "streamlit qpcr analysis v1.py"  # Expected: 0 (or count of non-triplicate usages if any)
```

### Final Checklist
- [ ] All purple gradients removed — monochrome theme applied
- [ ] All emojis removed from headers, tabs, sidebar
- [ ] Health Status Grid deleted, dead functions removed
- [ ] data_editor replaced with st.form + checkboxes + submit button
- [ ] Per-gene and per-sample heatmaps added to QC Overview
- [ ] All existing tests pass
- [ ] No business logic changed
- [ ] Korean text labels preserved
