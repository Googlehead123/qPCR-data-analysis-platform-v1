# Six UI & Analysis Fixes

## TL;DR

> **Quick Summary**: Apply 6 independent fixes to the qPCR Streamlit app — auto-rerun on exclusion change, graph export margin, NameError crash fix, triplicate stats display, graph styling, and color picker visibility.
> 
> **Deliverables**: All changes in `streamlit qpcr analysis v1.py` only
> - FIX 1: Auto-rerun analysis when QC exclusions change
> - FIX 2: Graph export legend cutoff fix
> - FIX 3: gen_col2 NameError crash fix
> - FIX 4: Mean CT & CT SD in triplicate browser
> - FIX 5: Graph styling (remove title, y-axis line, default colors)
> - FIX 6: Color picker expander visible by default
> 
> **Estimated Effort**: Medium
> **Parallel Execution**: NO — sequential (single file, overlapping line ranges)
> **Critical Path**: FIX 3 → FIX 2 → FIX 5 → FIX 6 → FIX 4 → FIX 1 → pytest

---

## Context

### Original Request
User requested 6 specific improvements identified during QC workflow testing. All requirements were fully defined in the previous planning session with exact line numbers and implementation details.

### Metis Review
**Identified Gaps** (addressed):
- Snapshot comparison for FIX 1 needs to handle both `excluded_wells` dict AND `excluded_samples` set → included both
- Test count is 94 not 79 → corrected
- FIX 1 edge case: don't auto-rerun if no `processed_data` exists → added guard
- FIX 3: ensure download button renders in correct container after removing `gen_col2` → noted

---

## Work Objectives

### Core Objective
Fix 6 UI/analysis issues in the qPCR Streamlit app to improve usability and correctness.

### Definition of Done
- [x] All 6 fixes applied to `streamlit qpcr analysis v1.py`
- [x] `pytest tests/ --ignore=tests/test_ppt_report.py -q` → all pass, 0 fail
- [x] `grep "gen_col2" "streamlit qpcr analysis v1.py"` → 0 matches

### Must NOT Have (Guardrails)
- NO test file modifications
- NO new dependencies in requirements.txt
- NO module extraction (single-file architecture)
- NO debug prints or excessive comments
- NO refactoring of surrounding code
- NO UI additions beyond what's specified (no "stale analysis" banners, no editable stats)

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (pytest)
- **User wants tests**: NO new tests for these fixes (UI/styling changes)
- **Framework**: pytest
- **QA approach**: Automated grep verification per fix + final pytest gate

---

## Execution Strategy

### Sequential Order (safest-first)
```
FIX 3 (crash fix) → FIX 2 (margin) → FIX 5 (styling) → FIX 6 (expander) → FIX 4 (stats display) → FIX 1 (auto-rerun) → pytest
```
Rationale: FIX 3 is a crash bug (highest priority). FIX 1 is most complex (last, so earlier fixes don't shift line numbers unpredictably).

---

## TODOs

- [x] 1. FIX 3: Remove gen_col2 NameError crash

  **What to do**:
  - Find lines ~5495-5513 where `with gen_col2:` is referenced
  - The column layout was removed in a previous refactor but `gen_col2` reference was left behind
  - Replace the `with gen_col2:` block with a simple inline download button (no column wrapper)
  - Keep the `if "ppt_bytes" in st.session_state` guard logic
  - Remove ALL references to `gen_col2` in the file

  **Must NOT do**:
  - Don't define `gen_col2` to "fix" it — the column layout is gone by design
  - Don't restructure the PPT export section

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]
    - `git-master`: Will commit after this fix
  - **Skills Evaluated but Omitted**:
    - `frontend-ui-ux`: Not needed — just removing a broken reference

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (Task 1 of 6)
  - **Blocks**: None
  - **Blocked By**: None

  **References**:
  - `streamlit qpcr analysis v1.py:5495-5513` — The broken `gen_col2` block to replace
  - `streamlit qpcr analysis v1.py:5470-5494` — Surrounding PPT export context (the `gen_col1` was also removed previously)

  **Acceptance Criteria**:
  ```bash
  grep "gen_col2" "streamlit qpcr analysis v1.py"
  # Assert: 0 matches (exit code 1)
  ```

  **Commit**: YES
  - Message: `fix(export): remove gen_col2 NameError in PPT export section`
  - Files: `streamlit qpcr analysis v1.py`

---

- [x] 2. FIX 2: Graph export legend cutoff

  **What to do**:
  - There are TWO bottom margin values in `create_gene_graph`:
    1. **Line 2080**: `gene_margins` default dict has `"b": 100` — change to `"b": 160`
    2. **Line 2158**: `gene_margins.get("b", 120)` fallback — change to `gene_margins.get("b", 160)`
  - Update ALL THREE export paths to add `margin=dict(b=180)`:
    1. **Per-gene download** (line 5145): `fig_copy.update_layout(...)` — add `margin=dict(b=180)`
    2. **Batch ZIP export** (line 5219): `fig_copy.update_layout(...)` — add `margin=dict(b=180)`
    3. **Full report ZIP export** (line 5310): `fig_copy.update_layout(...)` — add `margin=dict(b=180)`

  **Must NOT do**:
  - Don't change other margin values (left, right, top)
  - Don't modify the significance legend position (y=-0.15 stays)

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (Task 2 of 6)
  - **Blocks**: None
  - **Blocked By**: Task 1

  **References**:
  - `streamlit qpcr analysis v1.py:2080` — `gene_margins` default dict: `{"l": 80, "r": 80, "t": 100, "b": 100}` — change `"b": 100` → `"b": 160`
  - `streamlit qpcr analysis v1.py:2158` — Fallback: `gene_margins.get("b", 120)` — change `120` → `160`
  - `streamlit qpcr analysis v1.py:2109` — Significance legend at y=-0.15 (context only, don't change)
  - `streamlit qpcr analysis v1.py:5145-5150` — Per-gene export `fig_copy.update_layout` — add `margin=dict(b=180)`
  - `streamlit qpcr analysis v1.py:5219-5221` — Batch ZIP export `fig_copy.update_layout` — add `margin=dict(b=180)`
  - `streamlit qpcr analysis v1.py:5310-5312` — Full report ZIP export `fig_copy.update_layout` — add `margin=dict(b=180)`

  **Acceptance Criteria**:
  ```bash
  grep -n '"b": 160\|"b":160' "streamlit qpcr analysis v1.py"
  # Assert: finds updated default in gene_margins (line ~2080)
  grep -n 'get("b", 160)' "streamlit qpcr analysis v1.py"
  # Assert: finds updated fallback (line ~2158)
  grep -c "b=180\|b=.*180" "streamlit qpcr analysis v1.py"
  # Assert: ≥3 (three export paths)
  ```

  **Commit**: YES
  - Message: `fix(graph): increase bottom margin to prevent legend cutoff in export`
  - Files: `streamlit qpcr analysis v1.py`

---

- [x] 3. FIX 5: Graph styling changes

  **What to do**:
  - **5a: Remove graph title** — In `create_gene_graph` around line 2121-2130, set `title=None` or remove the title configuration. Gene name is already in y-axis label.
  - **5b: Add y-axis black line** — In `y_axis_config` dict around line 2052-2063, add: `showline=True, linewidth=1, linecolor="black", mirror=False`
  - **5c: Update DEFAULT_GROUP_COLORS** — Around line 30-38, change:
    - Negative Control → `#FFFFFF` (white)
    - Positive Control / Inducer → `#333333` (dark gray)

  **Must NOT do**:
  - Don't restructure the graph generation pipeline
  - Don't add new color themes or theme system
  - Don't change non-default colors (user-customized colors are separate)

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (Task 3 of 6)
  - **Blocks**: None
  - **Blocked By**: Task 2

  **References**:
  - `streamlit qpcr analysis v1.py:30-38` — `DEFAULT_GROUP_COLORS` dict
  - `streamlit qpcr analysis v1.py:2052-2063` — `y_axis_config` dict
  - `streamlit qpcr analysis v1.py:2121-2130` — Graph title configuration

  **Acceptance Criteria**:
  ```bash
  grep -n "FFFFFF" "streamlit qpcr analysis v1.py"
  # Assert: found in DEFAULT_GROUP_COLORS
  grep -n "showline.*True" "streamlit qpcr analysis v1.py"
  # Assert: found in y_axis_config
  grep -n "linecolor.*black" "streamlit qpcr analysis v1.py"
  # Assert: found in y_axis_config
  ```

  **Commit**: YES
  - Message: `style(graph): remove title, add y-axis line, update default control colors`
  - Files: `streamlit qpcr analysis v1.py`

---

- [x] 4. FIX 6: Color picker expander visible by default

  **What to do**:
  - Find the color picker expander around line 4579
  - Change `expanded=False` → `expanded=True`
  - Update the expander label to `"🎨 Condition Colors (applies to all genes)"`
  - Add `help="Check this to customize colors for individual bars within each gene graph"` to the per-bar edit checkbox around line 4562

  **Must NOT do**:
  - Don't redesign the graph settings UI
  - Don't move the expander to a different location

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (Task 4 of 6)
  - **Blocks**: None
  - **Blocked By**: Task 3

  **References**:
  - `streamlit qpcr analysis v1.py:4579` — Color picker expander definition
  - `streamlit qpcr analysis v1.py:4562` — Per-bar edit checkbox

  **Acceptance Criteria**:
  ```bash
  grep -n "expanded=True" "streamlit qpcr analysis v1.py" | grep -i "color"
  # Assert: finds the color picker expander with expanded=True
  ```

  **Commit**: YES
  - Message: `ui(graph): expand color picker by default, add help text`
  - Files: `streamlit qpcr analysis v1.py`

---

- [x] 5. FIX 4: Add Mean CT and CT SD to triplicate browser

  **What to do**:
  - In the QC triplicate browser, after line ~3336 (end of data_editor change processing)
  - Inside each gene expander, for each sample, add a summary showing:
    - Number of included wells (n)
    - Mean CT of included wells (rounded to 2 decimal places)
    - CT SD of included wells (sample std, `ddof=1`, rounded to 3 decimal places)
  - Use `st.caption()` with format: `f"n={n_included}, Mean CT={mean_ct:.2f}, SD={ct_sd:.3f}"`
  - If `n_included < 2`, show SD as "N/A" (can't compute sample SD with 1 value)
  - If `n_included == 0`, show "No wells included"
  - Use `lambda w, g=gene, s=sample:` pattern to avoid closure variable bugs
  - Only show stats for wells that are currently included (not excluded)
  - Compute from CT values in the data_editor DataFrame, filtering by the `Include` column

  **Must NOT do**:
  - Don't add editable stats or interactive elements
  - Don't add filtering or export capability
  - Don't modify the data_editor itself

  **Recommended Agent Profile**:
  - **Category**: `unspecified-low`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (Task 5 of 6)
  - **Blocks**: None
  - **Blocked By**: Task 4

  **References**:
  - `streamlit qpcr analysis v1.py:3336` — End of data_editor change processing in triplicate browser
  - `streamlit qpcr analysis v1.py:3200-3350` — Full triplicate browser section for context
  - `streamlit qpcr analysis v1.py:100-117` — Session state initialization pattern (for reference)

  **Acceptance Criteria**:
  ```bash
  grep -n "Mean CT\|mean_ct\|CT SD\|ct_sd" "streamlit qpcr analysis v1.py"
  # Assert: finds stat display code in triplicate browser section
  ```

  **Commit**: YES
  - Message: `feat(qc): show Mean CT and CT SD per sample in triplicate browser`
  - Files: `streamlit qpcr analysis v1.py`

---

- [x] 6. FIX 1: Auto-rerun analysis when QC exclusions change

  **What to do**:
  - **Step A: Store snapshot after analysis** — After `run_full_analysis` succeeds (~line 1770), store a snapshot of the current exclusion state:
    ```python
    st.session_state['_exclusion_snapshot'] = {
        'excluded_wells': {str(k): sorted(v) for k, v in st.session_state.excluded_wells.items()},
        'excluded_samples': sorted(st.session_state.get('excluded_samples', set())),
        'ttest_type': st.session_state.get('ttest_type', 'welch')
    }
    ```
    (Convert to JSON-serializable form for reliable comparison)
  - **Step B: Store ref/cmp sample keys** — Before the "Run Full Analysis Now" button (~line 4299), store `ref_sample_key`, `cmp_sample_key`, `cmp_sample_key_2`, and `use_second_comparison` in session state so auto-rerun can access them. The app supports a secondary comparison (lines 4195-4220) and auto-rerun must preserve it.
  - **Step C: Compare at tab3 entry** — At the top of `with tab3:` (~line 4320), compare current state vs snapshot. If different AND `processed_data` exists:
    1. Build current snapshot same way as Step A
    2. Compare with stored snapshot
    3. If different, call `run_full_analysis` with the stored ref/cmp keys
    4. Update the snapshot
  - **Edge case guards**:
    - Don't auto-rerun if `processed_data` is None/empty (analysis never run)
    - Don't auto-rerun if `_exclusion_snapshot` doesn't exist yet (first run)
    - Handle gracefully if session state was reset (data re-upload)

  **Must NOT do**:
  - Don't add "analysis stale" UI banners or notifications
  - Don't add debounce or delay logic
  - Don't modify the manual "Run Full Analysis Now" button behavior

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
  - **Skills**: [`git-master`]
    - `git-master`: Final commit with all verification

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (Task 6 of 6, most complex)
  - **Blocks**: None
  - **Blocked By**: Tasks 1-5

  **References**:
  - `streamlit qpcr analysis v1.py:1770` — After `run_full_analysis` succeeds (snapshot storage point)
  - `streamlit qpcr analysis v1.py:4299` — "Run Full Analysis Now" button (ref/cmp key storage point)
  - `streamlit qpcr analysis v1.py:4320` — Top of `with tab3:` (comparison/auto-rerun point)
  - `streamlit qpcr analysis v1.py:100-117` — Session state initialization pattern
  - `streamlit qpcr analysis v1.py:1700-1780` — `run_full_analysis` method for understanding parameters
  - `streamlit qpcr analysis v1.py:4195-4220` — Secondary comparison UI (`use_second_comparison`, `cmp_sample_key_2`)
  - `streamlit qpcr analysis v1.py:4291-4307` — How secondary comparison is passed to `run_full_analysis`

  **Acceptance Criteria**:
  ```bash
  grep -c "_exclusion_snapshot\|_snapshot" "streamlit qpcr analysis v1.py"
  # Assert: ≥2 (store + compare)
  grep -n "auto.*rerun\|auto.*analysis\|snapshot" "streamlit qpcr analysis v1.py"
  # Assert: finds snapshot logic
  ```

  **Commit**: YES
  - Message: `feat(analysis): auto-rerun when QC exclusions or t-test type change`
  - Files: `streamlit qpcr analysis v1.py`

---

- [x] 7. Final verification

  **What to do**:
  - Run full test suite: `pytest tests/ --ignore=tests/test_ppt_report.py -q`
  - Verify all 94 tests pass
  - Run all grep verification commands from tasks 1-6
  - Push to origin main

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (final gate)
  - **Blocks**: None
  - **Blocked By**: Tasks 1-6

  **Acceptance Criteria**:
  ```bash
  pytest tests/ --ignore=tests/test_ppt_report.py -q
  # Assert: all passed, 0 failed
  grep "gen_col2" "streamlit qpcr analysis v1.py"
  # Assert: 0 matches
  ```

  **Commit**: NO (already committed per-task)

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `fix(export): remove gen_col2 NameError in PPT export section` | `streamlit qpcr analysis v1.py` | grep gen_col2 → 0 |
| 2 | `fix(graph): increase bottom margin to prevent legend cutoff in export` | `streamlit qpcr analysis v1.py` | grep margin values |
| 3 | `style(graph): remove title, add y-axis line, update default control colors` | `streamlit qpcr analysis v1.py` | grep styling changes |
| 4 | `ui(graph): expand color picker by default, add help text` | `streamlit qpcr analysis v1.py` | grep expanded=True |
| 5 | `feat(qc): show Mean CT and CT SD per sample in triplicate browser` | `streamlit qpcr analysis v1.py` | grep Mean CT |
| 6 | `feat(analysis): auto-rerun when QC exclusions or t-test type change` | `streamlit qpcr analysis v1.py` | grep snapshot |
| 7 | — | — | pytest → all pass |

---

## Success Criteria

### Verification Commands
```bash
pytest tests/ --ignore=tests/test_ppt_report.py -q  # Expected: 94 passed
grep "gen_col2" "streamlit qpcr analysis v1.py"      # Expected: no matches
```

### Final Checklist
- [x] All 6 fixes applied
- [x] No gen_col2 references remain
- [x] All tests pass
- [x] Each fix committed separately with descriptive message
- [x] Pushed to origin main
