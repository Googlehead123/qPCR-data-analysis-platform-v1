# QC Tab Grid/Matrix Redesign - Selection Independence Fix

## TL;DR

> **Quick Summary**: Redesign the QC Triplicate Browser from a dropdown-based selection to an interactive Grid/Matrix UI where genes are rows, samples are columns, and each cell represents an independent (gene, sample) triplicate with its own selection state.
> 
> **Deliverables**:
> - New `TriplicateGrid` component with visual matrix layout
> - Independent state management per (gene, sample) combination  
> - Full TDD test suite for state logic and UI behavior
> - Backward compatibility with existing `excluded_wells` set
> 
> **Estimated Effort**: Medium (8-12 hours)
> **Parallel Execution**: YES - 3 waves
> **Critical Path**: Task 1 (State Tests) -> Task 2 (State Logic) -> Task 4 (UI Tests) -> Task 5 (Grid UI)

---

## Context

### Original Request
Fix QC tab bug where well selections are incorrectly linked between genes/samples. Root cause: static selectbox key `"selected_triplicate"` while data_editor uses dynamic key `f"editor_{selected_sample}_{selected_target}"`. User wants a complete Grid/Matrix UI redesign for the best possible user experience.

### Interview Summary
**Key Discussions**:
- **UI Direction**: User chose Option C - Grid/Matrix UI (best possible redesign, not quick fix)
- **Filter Behavior**: Reset selection to clean slate when filters change
- **Other Issues**: None - focus solely on selection independence bug + UI redesign
- **Testing**: Full TDD approach with unit tests, integration tests

**Research Findings**:
- `QualityControl.get_triplicate_data()` returns DataFrame grouped by (Sample, Target) with stats
- `QualityControl.get_wells_for_triplicate(data, sample, target)` returns individual wells for editing
- Existing test infrastructure in `tests/test_quality_control.py` uses pytest with Streamlit mocking via `conftest.py`
- `excluded_wells` is a global set that must be maintained for backward compatibility
- `excluded_wells_history` provides undo functionality

### Current Implementation Analysis
**File**: `streamlit qpcr analysis v1.py`  
**QC Tab Location**: Lines 2676-3356

**Bug Flow**:
```
User selects Gene A, Sample 1 -> selected_triplicate = 0 (index in filtered list)
User changes filter to Gene B -> filtered list changes  
Index 0 now points to different (Gene B, Sample X) -> WRONG!
```

**Key Methods**:
- `QualityControl.get_triplicate_data(data, excluded_wells)` - Lines 388-477
- `QualityControl.get_wells_for_triplicate(data, sample, target)` - Lines 480-530

---

## Work Objectives

### Core Objective
Replace the dropdown-based triplicate selector with a visual Grid/Matrix interface where each (gene, sample) cell is independently selectable, eliminating the index-based selection bug and improving UX.

### Concrete Deliverables
- New helper functions for grid state management
- Interactive Grid/Matrix UI component in Triplicate Browser sub-tab
- Complete TDD test suite (state logic + UI integration)
- Updated session_state schema with per-cell selection tracking

### Definition of Done
- [x] `pytest tests/test_qc_grid.py -v` passes (all new tests GREEN)
- [x] Clicking Cell(Gene_A, Sample_1) does NOT affect Cell(Gene_B, Sample_2)
- [x] `excluded_wells` set continues to work globally across all views
- [x] Filters (gene/sample/status) filter the grid but reset selection state
- [x] Streamlit app runs without errors: `streamlit run "streamlit qpcr analysis v1.py"`

### Must Have
- Grid with genes as rows, samples as columns
- Visual status indicators per cell (OK=green, Warning=yellow, Error=red)
- Click cell to expand detail editor for that triplicate
- Independent state: editing Cell A never affects Cell B
- Backward compatible with `excluded_wells` global set

### Must NOT Have (Guardrails)
- **NO dropdown selectbox** for triplicate selection (replace entirely)
- **NO shared index state** between different (gene, sample) combinations
- **NO changes to QualityControl class methods** (get_triplicate_data, get_wells_for_triplicate)
- **NO modifications to other tabs** (only qc_tab1: Triplicate Browser)
- **NO changes to excluded_wells data structure** (must remain a set)
- **NO over-engineering**: Avoid creating new classes when simple functions suffice
- **NO premature optimization**: Handle 20x20 grid first, optimize later if needed

---

## Verification Strategy (MANDATORY)

### Test Decision
- **Infrastructure exists**: YES (pytest with conftest.py)
- **User wants tests**: TDD (RED-GREEN-REFACTOR)
- **Framework**: pytest with Streamlit mocking

### TDD Workflow for Each Task

**Task Structure:**
1. **RED**: Write failing test first
   - Test file: `tests/test_qc_grid.py`
   - Test command: `pytest tests/test_qc_grid.py -v`
   - Expected: FAIL (test exists, implementation doesn't)
2. **GREEN**: Implement minimum code to pass
   - Command: `pytest tests/test_qc_grid.py -v`
   - Expected: PASS
3. **REFACTOR**: Clean up while keeping green
   - Command: `pytest tests/ -v` (full suite)
   - Expected: PASS (no regressions)

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Start Immediately):
+-- Task 1: Write state management tests (RED)
+-- Task 3: Write UI helper tests (RED)

Wave 2 (After Wave 1):
+-- Task 2: Implement state management (GREEN for Task 1)
+-- Task 4: Implement UI helpers (GREEN for Task 3)

Wave 3 (After Wave 2):
+-- Task 5: Write grid UI integration tests (RED)
+-- Task 6: Implement grid UI component (GREEN for Task 5)

Wave 4 (Final):
+-- Task 7: Integration testing and manual QA

Critical Path: Task 1 -> Task 2 -> Task 5 -> Task 6 -> Task 7
Parallel Speedup: ~35% faster than sequential
```

### Dependency Matrix

| Task | Depends On | Blocks | Can Parallelize With |
|------|------------|--------|---------------------|
| 1 | None | 2 | 3 |
| 2 | 1 | 5, 6 | 4 |
| 3 | None | 4 | 1 |
| 4 | 3 | 5, 6 | 2 |
| 5 | 2, 4 | 6 | None |
| 6 | 5 | 7 | None |
| 7 | 6 | None | None |

---

## TODOs

### Task 1: Write State Management Unit Tests (RED Phase)

**What to do**:
- Create `tests/test_qc_grid.py` with comprehensive test cases
- Test `get_grid_cell_key(gene, sample)` - generates unique key
- Test `get_selected_cell(session_state)` - retrieves current selection
- Test `set_selected_cell(session_state, gene, sample)` - sets selection
- Test `clear_selected_cell(session_state)` - clears selection
- Test `is_cell_selected(session_state, gene, sample)` - checks if specific cell selected
- Test independence: setting cell A does NOT affect cell B

**Must NOT do**:
- Don't implement the actual functions yet (TDD: tests first)
- Don't import functions that don't exist (use import guards)

**Recommended Agent Profile**:
- **Category**: `quick`
  - Reason: Test file creation is straightforward, well-defined scope
- **Skills**: [`git-master`]
  - `git-master`: Atomic commit after test file creation

**Parallelization**:
- **Can Run In Parallel**: YES
- **Parallel Group**: Wave 1 (with Task 3)
- **Blocks**: Task 2
- **Blocked By**: None

**References**:

**Pattern References** (existing code to follow):
- `tests/test_quality_control.py:16-82` - Test class structure with mock_streamlit fixture
- `tests/test_quality_control.py:156-193` - TestQualityControlGetTriplicateData class pattern
- `conftest.py:18-35` - MockSessionState class for testing session_state

**API/Type References** (contracts to implement against):
- `conftest.py:45-86` - `_create_mock_streamlit()` shows available mock attributes
- `conftest.py:92-102` - `mock_streamlit` fixture pattern

**Test References** (testing patterns to follow):
- `tests/test_quality_control.py:16-28` - Import pattern with `import_module`
- `tests/test_quality_control.py:176-193` - Testing with excluded_wells pattern

**WHY Each Reference Matters**:
- `test_quality_control.py` shows exact import pattern needed to test functions in main file
- `conftest.py` shows how to mock session_state for state management tests
- Test structure shows class-based organization by feature

**Acceptance Criteria**:

**RED Phase Verification:**
- [x] Test file created: `tests/test_qc_grid.py`
- [x] Test imports use `import_module("streamlit qpcr analysis v1")` pattern
- [x] At least 6 test methods covering: key generation, get/set/clear/check selection, independence
- [x] `pytest tests/test_qc_grid.py -v` -> FAILS with `AttributeError: module has no attribute` (functions don't exist yet)

**Commit**: YES
- Message: `test(qc): add RED phase tests for grid state management`
- Files: `tests/test_qc_grid.py`
- Pre-commit: None (tests expected to fail)

---

### Task 2: Implement State Management Functions (GREEN Phase)

**What to do**:
- Add state management helper functions to main file after QualityControl class (around line 940)
- Implement `get_grid_cell_key(gene: str, sample: str) -> str`
- Implement `get_selected_cell(session_state) -> tuple[str, str] | None`
- Implement `set_selected_cell(session_state, gene: str, sample: str) -> None`
- Implement `clear_selected_cell(session_state) -> None`
- Implement `is_cell_selected(session_state, gene: str, sample: str) -> bool`
- Use session_state key: `qc_grid_selected_cell` storing `{"gene": str, "sample": str}` or `None`

**Must NOT do**:
- Don't modify QualityControl class methods
- Don't create a new class (simple functions are sufficient)
- Don't add UI code yet (pure state logic only)

**Recommended Agent Profile**:
- **Category**: `quick`
  - Reason: Simple function implementations, well-defined contracts from tests
- **Skills**: [`git-master`]
  - `git-master`: Atomic commit after making tests pass

**Parallelization**:
- **Can Run In Parallel**: NO (depends on Task 1)
- **Parallel Group**: Wave 2 (with Task 4)
- **Blocks**: Task 5, Task 6
- **Blocked By**: Task 1

**References**:

**Pattern References** (existing code to follow):
- `streamlit qpcr analysis v1.py:55-62` - Session state initialization pattern
- `streamlit qpcr analysis v1.py:2684-2687` - Pattern for checking/initializing session_state keys

**API/Type References** (contracts to implement against):
- Tests from Task 1 define the exact function signatures
- `st.session_state` behaves like a dict with attribute access

**Implementation Location**:
- Insert after line ~940 (after QualityControl class, before AnalysisEngine)
- Add section comment: `# ==================== QC GRID STATE MANAGEMENT ====================`

**WHY Each Reference Matters**:
- Session state patterns show how to safely check/initialize keys
- Tests define exact expected behavior - implement to satisfy tests

**Acceptance Criteria**:

**GREEN Phase Verification:**
- [x] Functions added to `streamlit qpcr analysis v1.py` around line 940
- [x] `pytest tests/test_qc_grid.py -v` -> PASS (all tests GREEN)
- [x] `pytest tests/ -v` -> PASS (no regressions in existing tests)

**Manual Verification:**
- [x] Using Python REPL:
  ```python
  > from importlib import import_module
  > spec = import_module("streamlit qpcr analysis v1")
  > spec.get_grid_cell_key("GAPDH", "Sample1")
  Expected: "GAPDH::Sample1" (or similar unique key)
  ```

**Commit**: YES
- Message: `feat(qc): implement grid state management functions (GREEN)`
- Files: `streamlit qpcr analysis v1.py`
- Pre-commit: `pytest tests/test_qc_grid.py -v`

---

### Task 3: Write UI Helper Tests (RED Phase)

**What to do**:
- Add tests to `tests/test_qc_grid.py` for UI helper functions
- Test `build_grid_matrix(triplicate_data) -> dict` - builds {gene: {sample: cell_data}}
- Test `get_cell_status_color(status) -> str` - returns CSS color for status
- Test `get_cell_display_text(cell_data) -> str` - returns compact display text
- Test grid handles empty data gracefully
- Test grid handles single gene/sample edge cases

**Must NOT do**:
- Don't implement the helpers yet (TDD: tests first)
- Don't test actual Streamlit rendering (only data transformation)

**Recommended Agent Profile**:
- **Category**: `quick`
  - Reason: Test additions are straightforward extensions
- **Skills**: [`git-master`]
  - `git-master`: Atomic commit for test additions

**Parallelization**:
- **Can Run In Parallel**: YES
- **Parallel Group**: Wave 1 (with Task 1)
- **Blocks**: Task 4
- **Blocked By**: None

**References**:

**Pattern References** (existing code to follow):
- `tests/test_qc_grid.py` (from Task 1) - Add to same file
- `streamlit qpcr analysis v1.py:433-465` - `get_health_status()` returns status strings to handle

**Test Data References**:
- `conftest.py:107-147` - `sample_qpcr_raw_data` fixture structure
- `streamlit qpcr analysis v1.py:412-421` - Column names from get_triplicate_data output

**WHY Each Reference Matters**:
- Status strings from `get_health_status()` determine color mapping
- Triplicate data structure defines what `build_grid_matrix` receives

**Acceptance Criteria**:

**RED Phase Verification:**
- [x] Tests added to `tests/test_qc_grid.py`
- [x] Test class: `TestQCGridUIHelpers`
- [x] At least 5 test methods covering: matrix building, color mapping, display text, edge cases
- [x] `pytest tests/test_qc_grid.py::TestQCGridUIHelpers -v` -> FAILS

**Commit**: YES
- Message: `test(qc): add RED phase tests for grid UI helpers`
- Files: `tests/test_qc_grid.py`
- Pre-commit: None (tests expected to fail)

---

### Task 4: Implement UI Helper Functions (GREEN Phase)

**What to do**:
- Add UI helper functions to main file (same section as Task 2)
- Implement `build_grid_matrix(triplicate_data: pd.DataFrame) -> dict`
  - Returns `{gene: {sample: {"mean_ct": float, "cv": float, "status": str, "n": int}}}`
- Implement `get_cell_status_color(status: str) -> str`
  - "OK" -> "#d4edda" (green), Warning -> "#fff3cd" (yellow), Error -> "#f8d7da" (red)
- Implement `get_cell_display_text(cell_data: dict) -> str`
  - Returns compact text like "n=3, CV=2.1%"

**Must NOT do**:
- Don't render any Streamlit components yet
- Don't modify the existing styled_overview dataframe code

**Recommended Agent Profile**:
- **Category**: `quick`
  - Reason: Simple helper implementations with clear test contracts
- **Skills**: [`git-master`]
  - `git-master`: Atomic commit after making tests pass

**Parallelization**:
- **Can Run In Parallel**: YES
- **Parallel Group**: Wave 2 (with Task 2)
- **Blocks**: Task 5, Task 6
- **Blocked By**: Task 3

**References**:

**Pattern References** (existing code to follow):
- `streamlit qpcr analysis v1.py:2864-2872` - `style_triplicate_row()` color mapping pattern
- `streamlit qpcr analysis v1.py:412-421` - Triplicate data column structure

**API/Type References** (contracts to implement against):
- Tests from Task 3 define exact function signatures
- `triplicate_data` DataFrame has columns: Sample, Target, n, Mean_CT, SD, CV_pct, Range, Status, Severity

**WHY Each Reference Matters**:
- Color mapping pattern in `style_triplicate_row()` ensures visual consistency
- Column names must match what `get_triplicate_data()` returns

**Acceptance Criteria**:

**GREEN Phase Verification:**
- [x] Functions added to `streamlit qpcr analysis v1.py`
- [x] `pytest tests/test_qc_grid.py::TestQCGridUIHelpers -v` -> PASS
- [x] `pytest tests/test_qc_grid.py -v` -> PASS (all grid tests)
- [x] `pytest tests/ -v` -> PASS (no regressions)

**Commit**: YES
- Message: `feat(qc): implement grid UI helper functions (GREEN)`
- Files: `streamlit qpcr analysis v1.py`
- Pre-commit: `pytest tests/test_qc_grid.py -v`

---

### Task 5: Write Grid UI Integration Tests (RED Phase)

**What to do**:
- Add integration tests to `tests/test_qc_grid.py`
- Test `render_triplicate_grid()` function behavior (not actual rendering)
- Test that selecting a cell updates session_state correctly
- Test that filter changes clear selection (reset behavior)
- Test that excluding wells from grid detail updates `excluded_wells` globally
- Test independence: verify Cell A selection is isolated from Cell B

**Must NOT do**:
- Don't test actual Streamlit UI rendering (mock it)
- Don't implement the render function yet (TDD)

**Recommended Agent Profile**:
- **Category**: `quick`
  - Reason: Integration tests follow established patterns
- **Skills**: [`git-master`]
  - `git-master`: Atomic commit for test additions

**Parallelization**:
- **Can Run In Parallel**: NO (depends on Task 2, Task 4)
- **Parallel Group**: Wave 3
- **Blocks**: Task 6
- **Blocked By**: Task 2, Task 4

**References**:

**Pattern References** (existing code to follow):
- `tests/test_analysis.py:62-88` - Testing with excluded_wells integration
- `conftest.py:45-86` - Mocking Streamlit components

**Test Data References**:
- `conftest.py:107-147` - `sample_qpcr_raw_data` fixture
- Functions from Tasks 2 and 4 are now available

**WHY Each Reference Matters**:
- `test_analysis.py` shows how to test excluded_wells interactions
- Mock patterns ensure tests don't require Streamlit runtime

**Acceptance Criteria**:

**RED Phase Verification:**
- [x] Tests added to `tests/test_qc_grid.py`
- [x] Test class: `TestQCGridIntegration`
- [x] At least 5 test methods covering: selection update, filter reset, excluded_wells sync, independence
- [x] `pytest tests/test_qc_grid.py::TestQCGridIntegration -v` -> FAILS

**Commit**: YES
- Message: `test(qc): add RED phase integration tests for grid UI`
- Files: `tests/test_qc_grid.py`
- Pre-commit: None (tests expected to fail)

---

### Task 6: Implement Grid UI Component (GREEN Phase + UI Replacement)

**What to do**:
- Replace existing Triplicate Browser content (lines 2770-3106) with Grid/Matrix UI
- Implement `render_triplicate_grid(data, excluded_wells, session_state)` function
- Create visual grid using `st.columns()` for matrix layout:
  - First column: Gene names (row headers)
  - Remaining columns: Samples (column headers + cells)
- Each cell shows: status color background, n value, CV% on hover
- Clicking cell calls `set_selected_cell()` and shows detail editor below
- Detail editor uses existing `st.data_editor` pattern with dynamic key
- Add filter reset logic: when filter changes, call `clear_selected_cell()`

**Must NOT do**:
- Don't modify other sub-tabs (qc_tab2, qc_tab3, qc_tab4)
- Don't change QualityControl class methods
- Don't remove the Quick Action buttons (keep them in detail view)

**Recommended Agent Profile**:
- **Category**: `visual-engineering`
  - Reason: UI component implementation with layout considerations
- **Skills**: [`frontend-ui-ux`, `git-master`]
  - `frontend-ui-ux`: Visual grid layout, color-coded cells, interactive elements
  - `git-master`: Atomic commit after passing tests

**Parallelization**:
- **Can Run In Parallel**: NO (depends on Task 5)
- **Parallel Group**: Wave 3 (after Task 5 tests written)
- **Blocks**: Task 7
- **Blocked By**: Task 5

**References**:

**Pattern References** (existing code to follow):
- `streamlit qpcr analysis v1.py:2778-2800` - Filter controls pattern (keep these)
- `streamlit qpcr analysis v1.py:2969-3006` - data_editor pattern for well editing
- `streamlit qpcr analysis v1.py:3032-3105` - Quick Action buttons (keep in detail view)
- `streamlit qpcr analysis v1.py:2926-2942` - Summary metrics display pattern

**State Management References** (use functions from Task 2):
- `get_selected_cell()`, `set_selected_cell()`, `clear_selected_cell()`

**UI Helper References** (use functions from Task 4):
- `build_grid_matrix()`, `get_cell_status_color()`, `get_cell_display_text()`

**Code to Replace** (lines 2894-3106):
- Remove: Selectbox triplicate selector (line 2907-2912)
- Remove: Static key `"selected_triplicate"` usage
- Keep: Filter controls at top
- Keep: data_editor for well editing (use in detail panel)
- Keep: Quick Action buttons (move to detail panel)

**Layout Structure**:
```
+-- Filter Row (gene, sample, status filters) [KEEP]
|
+-- Grid Matrix
|   +-- Header Row: [Gene] [Sample1] [Sample2] [Sample3] ...
|   +-- Row 1:      [GAPDH] [cell]   [cell]    [cell]   ...
|   +-- Row 2:      [COL1A1][cell]   [cell]    [cell]   ...
|
+-- Detail Panel (shown when cell selected)
|   +-- Summary metrics (gene, sample, status, CV%)
|   +-- data_editor for individual wells [KEEP existing pattern]
|   +-- Quick Action buttons [KEEP existing]
```

**WHY Each Reference Matters**:
- Filter controls pattern must be preserved for consistency
- data_editor pattern is proven to work with excluded_wells
- Quick Action buttons provide important QC workflows

**Acceptance Criteria**:

**GREEN Phase Verification:**
- [x] `pytest tests/test_qc_grid.py::TestQCGridIntegration -v` -> PASS
- [x] `pytest tests/ -v` -> PASS (no regressions)

**Manual Verification (Playwright browser):**
- [x] Run: `streamlit run "streamlit qpcr analysis v1.py"` (no errors)
- [x] Navigate to QC tab -> Triplicate Browser sub-tab
- [x] Verify: Grid matrix visible with genes as rows, samples as columns
- [x] Click Cell(GAPDH, Sample1) -> Detail panel appears
- [x] Click Cell(COL1A1, Sample2) -> Selection moves, detail updates
- [x] Verify: Previous selection does NOT affect new selection (independence)
- [x] In detail panel: Uncheck a well's "Include" checkbox
- [x] Verify: Well added to excluded_wells (check bottom status bar count)
- [x] Change gene filter -> Verify selection resets (clean slate)

**Commit**: YES
- Message: `feat(qc): implement Grid/Matrix UI for Triplicate Browser`
- Files: `streamlit qpcr analysis v1.py`
- Pre-commit: `pytest tests/test_qc_grid.py -v`

---

### Task 7: Final Integration Testing and Manual QA

**What to do**:
- Run full test suite and verify no regressions
- Perform comprehensive manual QA with Playwright browser automation
- Test edge cases: empty data, single gene, single sample, 100+ wells
- Verify excluded_wells persists correctly to Analysis tab
- Test undo functionality still works
- Document any issues found and fix them

**Must NOT do**:
- Don't skip any QA steps
- Don't ignore test failures (fix them)
- Don't leave console errors

**Recommended Agent Profile**:
- **Category**: `visual-engineering`
  - Reason: Browser-based verification with visual confirmation
- **Skills**: [`playwright`, `git-master`]
  - `playwright`: Browser automation for interactive testing
  - `git-master`: Final commit with all fixes

**Parallelization**:
- **Can Run In Parallel**: NO (final integration)
- **Parallel Group**: Wave 4 (final)
- **Blocks**: None (completion)
- **Blocked By**: Task 6

**References**:

**Test Commands**:
- `pytest tests/ -v` - Full test suite
- `pytest tests/ --cov="streamlit qpcr analysis v1" --cov-report=term-missing` - Coverage

**Manual QA Checklist**:
1. Upload sample CSV data
2. Navigate to QC tab
3. Test all filter combinations
4. Test cell selection independence (critical bug fix)
5. Test excluded_wells propagation
6. Test undo button
7. Proceed to Mapping tab and verify excluded wells carried over
8. Run analysis and verify excluded wells respected

**Acceptance Criteria**:

**Automated Verification:**
- [x] `pytest tests/ -v` -> ALL PASS
- [x] `pytest tests/ --cov="streamlit qpcr analysis v1"` -> Coverage maintained or improved

**Manual Verification (Playwright):**
- [x] QC Grid displays correctly for sample data
- [x] Cell selection is independent (CRITICAL)
- [x] Excluded wells appear in Analysis results
- [x] No console errors in browser
- [x] No Streamlit exceptions

**Evidence Required:**
- [x] Screenshot of Grid UI with selections
- [x] Screenshot of independence test (two cells, different states)
- [x] Test output log showing all tests pass

**Commit**: YES
- Message: `test(qc): verify Grid/Matrix UI integration complete`
- Files: Any final fixes
- Pre-commit: `pytest tests/ -v`

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `test(qc): add RED phase tests for grid state management` | tests/test_qc_grid.py | None (expected fail) |
| 2 | `feat(qc): implement grid state management functions (GREEN)` | streamlit qpcr analysis v1.py | pytest tests/test_qc_grid.py |
| 3 | `test(qc): add RED phase tests for grid UI helpers` | tests/test_qc_grid.py | None (expected fail) |
| 4 | `feat(qc): implement grid UI helper functions (GREEN)` | streamlit qpcr analysis v1.py | pytest tests/test_qc_grid.py |
| 5 | `test(qc): add RED phase integration tests for grid UI` | tests/test_qc_grid.py | None (expected fail) |
| 6 | `feat(qc): implement Grid/Matrix UI for Triplicate Browser` | streamlit qpcr analysis v1.py | pytest tests/test_qc_grid.py |
| 7 | `test(qc): verify Grid/Matrix UI integration complete` | any fixes | pytest tests/ |

---

## Success Criteria

### Verification Commands
```bash
# All tests pass
pytest tests/ -v

# Coverage check
pytest tests/ --cov="streamlit qpcr analysis v1" --cov-report=term-missing

# App runs without errors
streamlit run "streamlit qpcr analysis v1.py" --server.headless true &
sleep 5 && curl -s http://localhost:8501 | grep -q "qPCR"
# Expected: Grep finds "qPCR" (app loaded)
```

### Final Checklist
- [x] **Bug Fixed**: Cell selection is completely independent per (gene, sample)
- [x] **Grid UI**: Genes as rows, samples as columns, status color-coded
- [x] **TDD Complete**: All tests pass (RED->GREEN->REFACTOR cycle completed)
- [x] **Backward Compatible**: excluded_wells set works as before
- [x] **No Regressions**: All existing tests still pass
- [x] **Manual QA**: Browser testing confirms correct behavior
