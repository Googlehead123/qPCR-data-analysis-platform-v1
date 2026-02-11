# QC Tab Consolidation & Triplicate Browser UI Improvement

## TL;DR

> **Quick Summary**: Consolidate the 5 QC sub-tabs into 2 tabs (enhanced Triplicate Browser + compact QC Overview), improving usability while preserving all business logic and session state behavior.
> 
> **Deliverables**:
> - Restructured QC section: 2 tabs instead of 5
> - Enhanced Triplicate Browser with better UX
> - Consolidated "QC Overview" tab merging useful parts of tabs 2-5
> 
> **Estimated Effort**: Medium
> **Parallel Execution**: NO - sequential (single file, overlapping code regions)
> **Critical Path**: Task 1 → Task 2 → Task 3 → Task 4

---

## Context

### Original Request
Consolidate the QC Check section from 5 sub-tabs to 2 tabs. The user only uses the Triplicate Browser (tab 1). Tabs 2-5 (Plate Heatmap, Flagged Wells, Settings, Summary Check) should be merged into one simple "QC Overview" tab. The Triplicate Browser UI should be improved. All existing functions, calculations, and logic must remain unbroken.

### Current Architecture
- Single file: `streamlit qpcr analysis v1.py` (5480 lines)
- QC section: lines 3110-3953
- 5 sub-tabs created at line 3193
- Tab 2 (Plate Heatmap) contains a DUPLICATE well selection editor (lines 3434-3598) that duplicates Triplicate Browser functionality
- Tab 3 (Flagged Wells) has auto-detection + exclude checkboxes
- Tab 4 (Settings) has 4 QC threshold number_inputs
- Tab 5 (Summary Check) has per-gene-sample inclusion summary with metrics

### Key Functions/Classes (MUST NOT MODIFY)
- `QualityControl.get_qc_summary_stats()`, `.get_triplicate_data()`, `.create_plate_heatmap()`, `.get_replicate_stats()`, `.detect_outliers()`
- `render_triplicate_grid()` (line 1168)
- `is_well_excluded()`, `exclude_well()`, `include_well()`, `get_all_excluded_wells()`
- Session state: `excluded_wells` (dict), `excluded_wells_history` (list)

---

## Work Objectives

### Core Objective
Restructure 5 QC tabs into 2 while improving the primary Triplicate Browser UX.

### Concrete Deliverables
- Modified `streamlit qpcr analysis v1.py` with 2-tab QC structure
- Tab 1: Enhanced "🔬 Triplicate Browser" 
- Tab 2: Consolidated "📋 QC Overview" (heatmap + flagged wells + settings + summary)

### Definition of Done
- [ ] App launches without errors: `streamlit run "streamlit qpcr analysis v1.py"`
- [ ] QC section shows exactly 2 tabs
- [ ] All well exclusion/inclusion functionality works identically
- [ ] All QC threshold settings are accessible and functional
- [ ] Plate heatmap renders correctly
- [ ] Flagged wells detection works
- [ ] Summary check data is visible
- [ ] Session state `excluded_wells` behavior unchanged

### Must Have
- All `QualityControl` class methods called in the same way
- All session state variables work identically
- Undo functionality preserved
- Bottom status bar preserved

### Must NOT Have (Guardrails)
- NO changes to any class definitions or helper functions
- NO changes to session state variable structure
- NO removal of QC algorithm functionality (just UI reorganization)
- NO new dependencies or imports
- NO changes to code outside the QC tab section (lines ~3193-3953)
- DO NOT duplicate the well selection editor (remove the duplicate from Tab 2)

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (pytest)
- **User wants tests**: Manual verification (UI refactoring)
- **Framework**: pytest (existing, but this is UI work)

### Automated Verification

Each task includes verification via Playwright browser automation to confirm the Streamlit app renders correctly.

---

## Task Dependency Graph

| Task | Depends On | Reason |
|------|------------|--------|
| Task 1 | None | Restructure tab declarations - foundation for all other tasks |
| Task 2 | Task 1 | Enhanced Triplicate Browser content goes into qc_tab1 |
| Task 3 | Task 1 | Consolidated QC Overview content goes into qc_tab2 |
| Task 4 | Task 2, Task 3 | Final verification requires both tabs complete |

## Parallel Execution Graph

```
Wave 1 (Start immediately):
└── Task 1: Tab structure refactoring (foundation)

Wave 2 (After Wave 1):
├── Task 2: Enhanced Triplicate Browser (independent of Task 3)
└── Task 3: Consolidated QC Overview tab (independent of Task 2)

Wave 3 (After Wave 2):
└── Task 4: Integration verification & cleanup
```

**NOTE**: Despite the theoretical parallelism of Tasks 2 and 3, since they modify the SAME FILE in overlapping regions, they MUST be executed sequentially to avoid merge conflicts. Recommended order: Task 1 → Task 2 → Task 3 → Task 4.

Critical Path: Task 1 → Task 2 → Task 3 → Task 4

---

## TODOs

- [ ] 1. Restructure QC Tab Declarations (5 tabs → 2 tabs)

  **What to do**:
  - Change line 3193 from `qc_tab1, qc_tab2, qc_tab3, qc_tab4, qc_tab5 = st.tabs([...5 items...])` to `qc_tab1, qc_tab2 = st.tabs(["🔬 Triplicate Browser", "📋 QC Overview"])`
  - Remove references to `qc_tab3`, `qc_tab4`, `qc_tab5` variables
  - Keep the `with qc_tab1:` block placeholder for Task 2
  - Create `with qc_tab2:` block placeholder for Task 3

  **Must NOT do**:
  - Do not modify any content inside the tab blocks yet
  - Do not change the bottom status bar (lines 3914-3953)

  **Recommended Agent Profile**:
  - **Category**: `quick` - Simple variable rename and line change
  - **Skills**: [`python-programmer`]
    - `python-programmer`: Python syntax and Streamlit patterns

  **Skills Evaluated but Omitted**:
  - `frontend-ui-ux`: Not needed for structural change
  - `git-master`: Will commit in Task 4

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 1 (solo)
  - **Blocks**: Tasks 2, 3, 4
  - **Blocked By**: None

  **References**:

  **Pattern References**:
  - `streamlit qpcr analysis v1.py:3193-3201` - Current 5-tab declaration to modify

  **Acceptance Criteria**:

  ```bash
  # Verify the tab declaration was changed
  grep -n "qc_tab1, qc_tab2 = st.tabs" "streamlit qpcr analysis v1.py"
  # Assert: Returns exactly 1 match
  
  # Verify old 5-tab pattern is gone
  grep -c "qc_tab3\|qc_tab4\|qc_tab5" "streamlit qpcr analysis v1.py"
  # Assert: Returns 0
  ```

  **Commit**: NO (groups with Task 4)

---

- [ ] 2. Enhance Triplicate Browser Tab (qc_tab1)

  **What to do**:
  - **Keep all existing functionality** from current Tab 1 (lines 3204-3398)
  - **Integrate flagged wells alert banner** at the top: Run `QualityControl.detect_outliers(data, hk_gene)` and show a compact warning banner if flagged wells exist (e.g., `st.warning(f"⚠️ {n} wells flagged by QC algorithms")`) with an "Auto-exclude flagged" button
  - **Move QC Settings inline** as a collapsed expander at the top: Wrap the 4 threshold inputs (CT high, CT low, CV%, HK variation) in `st.expander("⚙️ QC Settings", expanded=False)` with a "Re-run QC" button, placed BEFORE the filter controls
  - **Improve the Health Status Grid**: Change from `expanded=False` to `expanded=True` by default so users see triplicate health immediately
  - **Improve per-gene expander**: Add color-coded status indicator to expander label (✅/⚠️/❌ based on triplicate health data for that gene)
  - **Add per-sample inline stats improvement**: Show CV% alongside Mean CT and SD in the per-sample statistics section
  - **Remove the duplicate well selection editor** that existed in old Tab 2 (lines 3434-3598) - this code is simply deleted since Tab 2 no longer exists

  **Must NOT do**:
  - Do not change `render_triplicate_grid()` function definition (line 1168)
  - Do not change `QualityControl` class methods
  - Do not change session state variable names or structure
  - Do not change `is_well_excluded()`, `exclude_well()`, `include_well()` logic
  - Do not change the data_editor column config structure (just add CV% to stats display)

  **Recommended Agent Profile**:
  - **Category**: `unspecified-low` - Moderate UI refactoring with careful code preservation
  - **Skills**: [`python-programmer`, `frontend-ui-ux`]
    - `python-programmer`: Python/Streamlit patterns, data manipulation
    - `frontend-ui-ux`: UI layout decisions, user experience improvements

  **Skills Evaluated but Omitted**:
  - `data-scientist`: Not modifying calculations
  - `git-master`: Commit in Task 4

  **Parallelization**:
  - **Can Run In Parallel**: NO (same file as Task 3)
  - **Parallel Group**: Wave 2 (sequential with Task 3)
  - **Blocks**: Task 3, Task 4
  - **Blocked By**: Task 1

  **References**:

  **Pattern References**:
  - `streamlit qpcr analysis v1.py:3204-3398` - Current Triplicate Browser implementation (KEEP and enhance)
  - `streamlit qpcr analysis v1.py:3721-3795` - QC Settings code to MOVE into expander inside Tab 1
  - `streamlit qpcr analysis v1.py:3601-3672` - Flagged wells detection logic to reuse for alert banner
  - `streamlit qpcr analysis v1.py:3674-3690` - "Exclude All Flagged" button pattern to reuse
  - `streamlit qpcr analysis v1.py:1168-1250` - `render_triplicate_grid()` function (DO NOT MODIFY, just change expanded parameter when calling)
  - `streamlit qpcr analysis v1.py:3283` - `QualityControl.get_triplicate_data()` call pattern
  - `streamlit qpcr analysis v1.py:3607` - `QualityControl.detect_outliers()` call pattern

  **API/Type References**:
  - `QualityControl.CT_HIGH_THRESHOLD`, `.CT_LOW_THRESHOLD`, `.CV_THRESHOLD`, `.HK_VARIATION_THRESHOLD` - Class attributes for settings inputs
  - `QualityControl.GRUBBS_ALPHA` - Read-only display in settings info table

  **Acceptance Criteria**:

  ```
  # Agent executes via playwright browser automation:
  1. Run: streamlit run "streamlit qpcr analysis v1.py" --server.port 8501
  2. Navigate to: http://localhost:8501
  3. Upload a test CSV file
  4. Navigate to the QC Check tab
  5. Assert: "Triplicate Browser" tab is visible
  6. Click on Triplicate Browser tab
  7. Assert: QC Settings expander exists (collapsed by default)
  8. Assert: Health Status Grid section is visible (expanded)
  9. Assert: Gene filter and Sample filter dropdowns exist
  10. Assert: Per-gene expanders with well editors are visible
  11. Screenshot: .sisyphus/evidence/task-2-triplicate-browser.png
  ```

  ```bash
  # Verify no syntax errors
  python -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"
  # Assert: Exit code 0
  
  # Verify old duplicate editor is removed
  grep -c "qc_summary_well_editor" "streamlit qpcr analysis v1.py"
  # Assert: Returns 0
  
  # Verify settings are now in Tab 1 area
  grep -n "QC Threshold Settings\|CT_HIGH_THRESHOLD" "streamlit qpcr analysis v1.py" | head -5
  # Assert: Line numbers fall within Tab 1 region
  ```

  **Commit**: NO (groups with Task 4)

---

- [ ] 3. Build Consolidated QC Overview Tab (qc_tab2)

  **What to do**:
  - Build new `with qc_tab2:` block containing consolidated content from old tabs 2-5
  - **Section 1: Plate Heatmap** (from old Tab 2, lines 3400-3431 ONLY - the heatmap + replicate stats, NOT the duplicate well editor)
    - Keep: `QualityControl.create_plate_heatmap()` call with plotly chart
    - Keep: `QualityControl.get_replicate_stats()` table with conditional styling
    - OMIT: The entire "Individual Well Selection" section (lines 3434-3598) - this was a duplicate of Triplicate Browser
  - **Section 2: Flagged Wells Summary** (from old Tab 3, simplified)
    - Show the flagged wells table as READ-ONLY (no exclude checkboxes - exclusion is done in Triplicate Browser now)
    - Keep: `QualityControl.detect_outliers(data, hk_gene)` call
    - Show: `st.dataframe()` instead of `st.data_editor()` for flagged wells (read-only view)
    - Keep: The "Exclude All Flagged" button (useful quick action)
    - OMIT: Individual per-well exclude checkboxes (done in Triplicate Browser)
    - OMIT: "Clear All Exclusions" and "Undo" buttons (available in Triplicate Browser and bottom bar)
  - **Section 3: Pre-Analysis Summary** (from old Tab 5, lines 3798-3912)
    - Keep: The full summary check with metrics and per-gene expanders
    - This is READ-ONLY informational, no changes needed to logic

  **Must NOT do**:
  - Do not include Settings here (moved to Triplicate Browser in Task 2)
  - Do not duplicate any well exclusion editor functionality
  - Do not modify `QualityControl` class methods
  - Do not change session state handling

  **Recommended Agent Profile**:
  - **Category**: `unspecified-low` - Moderate code reorganization
  - **Skills**: [`python-programmer`, `frontend-ui-ux`]
    - `python-programmer`: Python/Streamlit code patterns
    - `frontend-ui-ux`: Clean layout for consolidated view

  **Skills Evaluated but Omitted**:
  - `data-scientist`: Not modifying data logic

  **Parallelization**:
  - **Can Run In Parallel**: NO (same file)
  - **Parallel Group**: Wave 2 (after Task 2)
  - **Blocks**: Task 4
  - **Blocked By**: Task 1, Task 2

  **References**:

  **Pattern References**:
  - `streamlit qpcr analysis v1.py:3400-3431` - Plate heatmap + replicate stats (COPY this section, omit lines 3434-3598)
  - `streamlit qpcr analysis v1.py:3601-3718` - Flagged wells display (simplify to read-only)
  - `streamlit qpcr analysis v1.py:3798-3912` - Summary check (COPY as-is)
  - `streamlit qpcr analysis v1.py:3674-3690` - "Exclude All Flagged" button pattern to keep

  **API/Type References**:
  - `QualityControl.create_plate_heatmap(data, value_col="CT", excluded_wells=get_all_excluded_wells())` - Exact call signature
  - `QualityControl.get_replicate_stats(data)` - Returns DataFrame with Status column
  - `QualityControl.detect_outliers(data, hk_gene)` - Returns DataFrame with Flagged column

  **Acceptance Criteria**:

  ```
  # Agent executes via playwright browser automation:
  1. Navigate to: http://localhost:8501
  2. Navigate to QC Check tab
  3. Assert: "QC Overview" tab is visible
  4. Click on QC Overview tab
  5. Assert: Plate heatmap chart is visible
  6. Assert: Replicate statistics table is visible
  7. Assert: Flagged wells section exists
  8. Assert: Pre-Analysis Summary section exists with metrics
  9. Screenshot: .sisyphus/evidence/task-3-qc-overview.png
  ```

  ```bash
  # Verify syntax
  python -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"
  # Assert: Exit code 0
  
  # Verify all QC methods are still called
  grep -c "create_plate_heatmap\|get_replicate_stats\|detect_outliers\|get_triplicate_data" "streamlit qpcr analysis v1.py"
  # Assert: Returns >= 4 (all methods still referenced)
  ```

  **Commit**: NO (groups with Task 4)

---

- [ ] 4. Integration Verification & Commit

  **What to do**:
  - Run existing tests: `pytest tests/ -v`
  - Verify the app starts without errors
  - Verify the bottom status bar (lines 3914-3953) still works correctly (it references `excluded_wells` session state)
  - Fix any key conflicts (ensure all `st.button`, `st.data_editor`, etc. have unique keys across both tabs)
  - Ensure no orphaned references to `qc_tab3`, `qc_tab4`, `qc_tab5`
  - Commit all changes

  **Must NOT do**:
  - Do not add new features beyond the consolidation scope
  - Do not refactor code outside the QC section

  **Recommended Agent Profile**:
  - **Category**: `quick` - Verification and commit
  - **Skills**: [`python-programmer`, `git-master`]
    - `python-programmer`: Syntax verification
    - `git-master`: Clean atomic commit

  **Skills Evaluated but Omitted**:
  - `frontend-ui-ux`: Not needed for verification
  - `python-debugger`: Only if errors found

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 3 (final)
  - **Blocks**: None
  - **Blocked By**: Tasks 1, 2, 3

  **References**:

  **Pattern References**:
  - `streamlit qpcr analysis v1.py:3914-3953` - Bottom status bar (verify still works)
  - `tests/` directory - Existing test suite

  **Acceptance Criteria**:

  ```bash
  # Run existing tests
  pytest tests/ -v
  # Assert: All tests pass
  
  # Verify app starts
  timeout 15 streamlit run "streamlit qpcr analysis v1.py" --server.port 8502 &
  sleep 8
  curl -s http://localhost:8502 | head -5
  # Assert: Returns HTML content
  kill %1
  
  # Verify no references to old tabs
  grep -c "qc_tab3\|qc_tab4\|qc_tab5" "streamlit qpcr analysis v1.py"
  # Assert: Returns 0
  
  # Verify syntax
  python -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"
  # Assert: Exit code 0
  ```

  **Commit**: YES
  - Message: `refactor(qc): consolidate 5 QC tabs into 2 and improve Triplicate Browser UX`
  - Files: `streamlit qpcr analysis v1.py`
  - Pre-commit: `pytest tests/ -v`

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 4 (all tasks) | `refactor(qc): consolidate 5 QC tabs into 2 and improve Triplicate Browser UX` | `streamlit qpcr analysis v1.py` | `pytest tests/ -v` |

---

## Success Criteria

### Verification Commands
```bash
pytest tests/ -v                    # All existing tests pass
python -c "import ast; ast.parse(open('streamlit qpcr analysis v1.py').read())"  # No syntax errors
grep -c "qc_tab3\|qc_tab4\|qc_tab5" "streamlit qpcr analysis v1.py"            # Returns 0
```

### Final Checklist
- [ ] QC section shows exactly 2 tabs (Triplicate Browser + QC Overview)
- [ ] Triplicate Browser has inline QC Settings expander
- [ ] Triplicate Browser has flagged wells alert banner with auto-exclude button
- [ ] Health Status Grid expanded by default
- [ ] Per-gene expanders show status indicators
- [ ] QC Overview tab shows plate heatmap, flagged wells (read-only), and summary
- [ ] Duplicate well selection editor removed
- [ ] All `QualityControl` methods still called
- [ ] Session state `excluded_wells` works identically
- [ ] Bottom status bar functions correctly
- [ ] All existing tests pass
