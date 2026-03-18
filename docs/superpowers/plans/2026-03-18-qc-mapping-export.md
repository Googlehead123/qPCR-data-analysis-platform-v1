# QC, Mapping & Export Improvements — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add QC auto-exclusion for high SD triplicates, make mapping finalize reversible with stale analysis badge, fix Excel QC_Report sheet, and add replicate fold change sheet to Excel.

**Architecture:** All changes within existing Streamlit + Plotly architecture. New QC method in `qpcr/quality_control.py`, Excel changes in `qpcr/export.py`, UI changes in the monolithic main file. No new dependencies.

**Tech Stack:** Python 3.12, Streamlit, Plotly, pandas, numpy, scipy

**Spec:** `docs/superpowers/specs/2026-03-18-qc-mapping-export-design.md`

---

### Task 1: `find_high_sd_outliers()` method in QualityControl

**Files:**
- Modify: `qpcr/quality_control.py` (add new static method after `suggest_exclusions`, ~line 289)
- Test: `tests/test_quality_control.py`

- [ ] **Step 1: Write failing tests**

Add to `tests/test_quality_control.py`:

```python
class TestFindHighSdOutliers:
    def test_finds_outlier_in_high_sd_group(self, mock_streamlit):
        from qpcr.quality_control import QualityControl
        data = pd.DataFrame({
            "Well": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "Sample": ["S1", "S1", "S1", "S2", "S2", "S2"],
            "Target": ["COL1A1", "COL1A1", "COL1A1", "COL1A1", "COL1A1", "COL1A1"],
            "CT": [20.0, 20.1, 22.0, 18.0, 18.1, 18.2],  # S1 has SD ~1.1, S2 is clean
        })
        results = QualityControl.find_high_sd_outliers(data, excluded_wells={}, sd_threshold=0.5)
        assert len(results) == 1
        assert results[0]["Sample"] == "S1"
        assert results[0]["Well"] == "A3"  # Furthest from mean (~20.7)
        assert results[0]["group_sd"] > 0.5

    def test_skips_clean_groups(self, mock_streamlit):
        from qpcr.quality_control import QualityControl
        data = pd.DataFrame({
            "Well": ["A1", "A2", "A3"],
            "Sample": ["S1", "S1", "S1"],
            "Target": ["COL1A1", "COL1A1", "COL1A1"],
            "CT": [20.0, 20.05, 20.1],  # SD ~0.05, well under threshold
        })
        results = QualityControl.find_high_sd_outliers(data, excluded_wells={}, sd_threshold=0.5)
        assert len(results) == 0

    def test_skips_groups_with_fewer_than_3_replicates(self, mock_streamlit):
        from qpcr.quality_control import QualityControl
        data = pd.DataFrame({
            "Well": ["A1", "A2"],
            "Sample": ["S1", "S1"],
            "Target": ["COL1A1", "COL1A1"],
            "CT": [20.0, 25.0],  # High SD but only 2 replicates
        })
        results = QualityControl.find_high_sd_outliers(data, excluded_wells={}, sd_threshold=0.5)
        assert len(results) == 0

    def test_respects_excluded_wells_dict(self, mock_streamlit):
        from qpcr.quality_control import QualityControl
        data = pd.DataFrame({
            "Well": ["A1", "A2", "A3"],
            "Sample": ["S1", "S1", "S1"],
            "Target": ["COL1A1", "COL1A1", "COL1A1"],
            "CT": [20.0, 20.1, 25.0],  # A3 is outlier
        })
        # Exclude A3 already — only 2 replicates remain, should skip
        excluded = {("COL1A1", "S1"): {"A3"}}
        results = QualityControl.find_high_sd_outliers(data, excluded_wells=excluded, sd_threshold=0.5)
        assert len(results) == 0

    def test_respects_flat_set_excluded_wells(self, mock_streamlit):
        from qpcr.quality_control import QualityControl
        data = pd.DataFrame({
            "Well": ["A1", "A2", "A3"],
            "Sample": ["S1", "S1", "S1"],
            "Target": ["COL1A1", "COL1A1", "COL1A1"],
            "CT": [20.0, 20.1, 25.0],
        })
        excluded = {"A3"}  # flat set format
        results = QualityControl.find_high_sd_outliers(data, excluded_wells=excluded, sd_threshold=0.5)
        assert len(results) == 0

    def test_gene_filter_limits_scope(self, mock_streamlit):
        from qpcr.quality_control import QualityControl
        data = pd.DataFrame({
            "Well": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "Sample": ["S1", "S1", "S1", "S1", "S1", "S1"],
            "Target": ["COL1A1", "COL1A1", "COL1A1", "ELN", "ELN", "ELN"],
            "CT": [20.0, 20.1, 25.0, 18.0, 18.1, 23.0],  # Both genes have outliers
        })
        results = QualityControl.find_high_sd_outliers(
            data, excluded_wells={}, sd_threshold=0.5, gene_filter="COL1A1"
        )
        assert len(results) == 1
        assert results[0]["Target"] == "COL1A1"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_quality_control.py::TestFindHighSdOutliers -v`
Expected: FAIL with `AttributeError`

- [ ] **Step 3: Implement the method**

Add to `qpcr/quality_control.py` after `suggest_exclusions` method (~line 289):

```python
@staticmethod
def find_high_sd_outliers(
    data: pd.DataFrame,
    excluded_wells,
    sd_threshold: float = 0.5,
    gene_filter: str = None,
) -> list:
    """Find the worst replicate in each gene-sample group exceeding SD threshold.

    Args:
        data: Raw CT data with Well, Sample, Target, CT columns.
        excluded_wells: Either dict {(gene, sample): set(well_ids)} or flat set.
        sd_threshold: CT SD threshold (default 0.5).
        gene_filter: If set, only analyze this gene.

    Returns:
        List of dicts with Target, Sample, Well, CT, deviation, group_sd, group_mean, n_replicates.
    """
    if data is None or data.empty:
        return []

    df = data.copy()

    # Normalize excluded_wells and filter
    if isinstance(excluded_wells, dict):
        # Dict format: {(gene, sample): set(well_ids)}
        excl_flat = set()
        for well_set in excluded_wells.values():
            excl_flat.update(well_set)
        df = df[~df["Well"].isin(excl_flat)]
    elif excluded_wells:
        df = df[~df["Well"].isin(excluded_wells)]

    if gene_filter:
        df = df[df["Target"] == gene_filter]

    suggestions = []
    for (target, sample), group in df.groupby(["Target", "Sample"]):
        if len(group) < 3:
            continue

        ct_values = group["CT"].values
        sd = np.std(ct_values, ddof=1)

        if sd <= sd_threshold:
            continue

        mean_ct = np.mean(ct_values)
        deviations = np.abs(ct_values - mean_ct)
        worst_idx = np.argmax(deviations)
        worst_row = group.iloc[worst_idx]

        suggestions.append({
            "Target": target,
            "Sample": sample,
            "Well": worst_row["Well"],
            "CT": round(float(worst_row["CT"]), 2),
            "deviation": round(float(worst_row["CT"] - mean_ct), 3),
            "group_sd": round(float(sd), 3),
            "group_mean": round(float(mean_ct), 2),
            "n_replicates": len(group),
        })

    return suggestions
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_quality_control.py::TestFindHighSdOutliers -v`
Expected: ALL PASS

- [ ] **Step 5: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add qpcr/quality_control.py tests/test_quality_control.py
git commit -m "feat: add find_high_sd_outliers for QC auto-exclusion"
```

---

### Task 2: QC Auto-Exclude UI (Global + Per-Gene Buttons)

**Files:**
- Modify: `streamlit qpcr analysis v1.py` (QC tab section)

This is a Streamlit UI task. No automated tests.

- [ ] **Step 1: Find the QC Sub-Tab 1 section**

Search for `qc_tab1` or the Triplicate Browser section in the main file. The QC tab has nested sub-tabs. Find where the per-gene/sample filter dropdowns are.

- [ ] **Step 2: Add SD threshold slider and global button**

At the top of QC Sub-Tab 1 (before the gene/sample filters), add:

```python
# Auto-exclude high SD outliers
st.markdown("#### Auto-Exclude High SD Outliers")
auto_excl_cols = st.columns([2, 1, 1])
with auto_excl_cols[0]:
    sd_thresh = st.slider(
        "SD Threshold (Ct)",
        min_value=0.3, max_value=1.0, value=0.5, step=0.1,
        key="sd_threshold_slider",
        help="Groups with CT SD above this threshold will be flagged",
    )
with auto_excl_cols[1]:
    find_all_btn = st.button("🔍 Find All Outliers", key="find_all_sd_outliers", use_container_width=True)
```

- [ ] **Step 3: Add per-gene "Clean Triplicates" button**

Next to wherever the gene filter dropdown is, add:

```python
with auto_excl_cols[2]:
    # Only show if a specific gene is selected (not "All")
    if selected_gene_filter and selected_gene_filter != "All":
        clean_gene_btn = st.button(
            f"🧹 Clean {selected_gene_filter}", key="clean_gene_btn", use_container_width=True
        )
    else:
        clean_gene_btn = False
```

Note: `selected_gene_filter` is whatever variable holds the current gene filter selection. Read the QC tab code to find the exact variable name.

- [ ] **Step 4: Add preview table and apply logic**

After the buttons, add the preview/apply logic:

```python
# Handle global find
if find_all_btn:
    from qpcr.quality_control import QualityControl
    suggestions = QualityControl.find_high_sd_outliers(
        st.session_state.data,
        st.session_state.get("excluded_wells", {}),
        sd_threshold=sd_thresh,
    )
    st.session_state["_sd_outlier_suggestions"] = suggestions

# Handle per-gene find
if clean_gene_btn:
    from qpcr.quality_control import QualityControl
    suggestions = QualityControl.find_high_sd_outliers(
        st.session_state.data,
        st.session_state.get("excluded_wells", {}),
        sd_threshold=sd_thresh,
        gene_filter=selected_gene_filter,
    )
    st.session_state["_sd_outlier_suggestions"] = suggestions

# Display preview table if suggestions exist
suggestions = st.session_state.get("_sd_outlier_suggestions", [])
if suggestions:
    st.markdown(f"**Found {len(suggestions)} high-SD outlier(s):**")
    preview_df = pd.DataFrame(suggestions)
    preview_df["Include"] = True
    edited = st.data_editor(
        preview_df[["Include", "Target", "Sample", "Well", "CT", "deviation", "group_sd", "n_replicates"]],
        disabled=["Target", "Sample", "Well", "CT", "deviation", "group_sd", "n_replicates"],
        key="sd_outlier_preview",
        use_container_width=True,
    )

    apply_cols = st.columns([3, 1])
    with apply_cols[1]:
        if st.button("Apply Selected Exclusions", key="apply_sd_exclusions", type="primary", use_container_width=True):
            selected = edited[edited["Include"]]
            excl_dict = st.session_state.get("excluded_wells", {})
            if not isinstance(excl_dict, dict):
                excl_dict = {}
            for _, row in selected.iterrows():
                key = (row["Target"], row["Sample"])
                if key not in excl_dict:
                    excl_dict[key] = set()
                excl_dict[key].add(row["Well"])
            st.session_state.excluded_wells = excl_dict
            st.session_state.pop("_sd_outlier_suggestions", None)
            st.success(f"Excluded {len(selected)} well(s).")
            st.rerun()
    with apply_cols[0]:
        if st.button("Dismiss", key="dismiss_sd_suggestions"):
            st.session_state.pop("_sd_outlier_suggestions", None)
            st.rerun()
elif find_all_btn or clean_gene_btn:
    st.info("No triplicates exceed the SD threshold.")
```

- [ ] **Step 5: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add "streamlit qpcr analysis v1.py"
git commit -m "feat: add QC auto-exclude UI with preview table"
```

---

### Task 3: Reversible Finalize Mapping + Stale Badge

**Files:**
- Modify: `streamlit qpcr analysis v1.py` (Mapping tab ~line 4922, Analysis tab ~line 5303, Graphs tab ~line 5411)

- [ ] **Step 1: Initialize `analysis_stale` in session state**

Near the top of the file where other session state keys are initialized (~line 400), add:

```python
if "analysis_stale" not in st.session_state:
    st.session_state.analysis_stale = False
```

- [ ] **Step 2: Modify "Finalize Mapping" to "Finalize & Run Analysis" (rename only)**

At ~line 4922, change:
```python
# BEFORE:
if st.button("Finalize Mapping", type="primary", use_container_width=True, key="finalize_mapping"):
    st.session_state.mapping_finalized = True
    st.rerun()
st.caption("Lock in condition names and groups, then reorder samples below.")
```
to:
```python
# AFTER:
if st.button("Finalize & Run Analysis", type="primary", use_container_width=True, key="finalize_mapping"):
    st.session_state.mapping_finalized = True
    st.session_state.analysis_stale = False
    st.rerun()
st.caption("Lock in conditions and groups. You can edit later with the 'Edit Mapping' button.")
```

- [ ] **Step 3: Add "Edit Mapping" button after finalize**

At ~line 4927, where `if st.session_state.get("mapping_finalized", False):` starts the reorder section, add an "Edit Mapping" button BEFORE the reorder UI:

```python
if st.session_state.get("mapping_finalized", False):
    # Edit Mapping button
    edit_col1, edit_col2 = st.columns([3, 1])
    with edit_col2:
        if st.button("✏️ Edit Mapping", key="edit_mapping", use_container_width=True):
            st.session_state.mapping_finalized = False
            st.session_state.analysis_stale = True
            st.rerun()

    st.markdown("### 🔀 Drag to Reorder Samples")
    # ... rest of existing reorder code
```

- [ ] **Step 4: Add stale badge to Analysis tab**

At the top of the Analysis tab (~line 5303, `with tab3:`), after the auto-rerun check, add:

```python
if st.session_state.get("analysis_stale", False):
    st.warning("⚠️ Mapping changed since last analysis. Results may be outdated.")
    if st.button("Re-run Analysis", key="rerun_stale_analysis", type="primary"):
        ref_key = st.session_state.get("_last_ref_sample_key")
        cmp_key = st.session_state.get("_last_cmp_sample_key")
        cmp_key_2 = st.session_state.get("_last_cmp_sample_key_2")
        cmp_key_3 = st.session_state.get("_last_cmp_sample_key_3")
        if ref_key and cmp_key:
            ok = AnalysisEngine.run_full_analysis(ref_key, cmp_key, cmp_key_2, cmp_key_3)
            if ok:
                st.session_state.analysis_stale = False
                st.success("Analysis re-run complete.")
                st.rerun()
            else:
                st.error("Re-run failed. Check mapping configuration.")
        else:
            st.error("Cannot re-run — no previous analysis configuration found. Go to Mapping tab and run analysis.")
```

- [ ] **Step 5: Add stale badge to Graphs tab**

At the top of the Graphs tab (~line 5411), add the same stale warning (simpler, no re-run button — just a notice):

```python
if st.session_state.get("analysis_stale", False):
    st.warning("⚠️ Mapping changed since last analysis. Go to Analysis tab to re-run.")
```

- [ ] **Step 6: Clear stale flag on successful analysis**

Find where `run_full_analysis` is called and succeeds (~line 5290). After `st.success(success_msg)`, add:

```python
st.session_state.analysis_stale = False
```

- [ ] **Step 7: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 8: Commit**

```bash
git add "streamlit qpcr analysis v1.py"
git commit -m "feat: reversible Finalize Mapping with analysis stale badge"
```

---

### Task 4: Excel — Fix QC_Report + Add Replicate_FC Sheet

**Files:**
- Modify: `qpcr/export.py` (add `excluded_wells` param, add Replicate_FC sheet)
- Modify: `streamlit qpcr analysis v1.py` (3 `export_to_excel` call sites)
- Test: `tests/test_package.py` or new test

- [ ] **Step 1: Write tests**

Add to `tests/test_quality_control.py` (or a new test file):

```python
class TestExcelExportFixes:
    def test_export_with_qc_stats_populates_qc_sheet(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from qpcr.export import export_to_excel
        from qpcr.analysis import AnalysisEngine
        from qpcr.quality_control import QualityControl
        import openpyxl
        import io

        processed = {"COL1A1": AnalysisEngine.calculate_ddct(
            sample_qpcr_raw_data, "GAPDH", "Non-treated", set(), set(), sample_mapping
        )}
        params = {"Housekeeping_Gene": "GAPDH", "Reference_Sample": "Non-treated"}
        replicate_stats = QualityControl.get_replicate_stats(sample_qpcr_raw_data)
        qc_stats = {"total_wells": 18, "excluded_wells": 0, "active_wells": 18,
                     "ct_mean": 20.0, "ct_std": 3.0, "ct_min": 18.0, "ct_max": 26.0,
                     "high_ct_count": 0, "low_ct_count": 0, "total_triplicates": 6,
                     "healthy_triplicates": 6, "warning_triplicates": 0, "error_triplicates": 0,
                     "avg_cv_pct": 0.5, "max_cv_pct": 1.0, "health_score": 100.0}

        result = export_to_excel(
            sample_qpcr_raw_data, processed, params, sample_mapping,
            qc_stats=qc_stats, replicate_stats=replicate_stats,
        )
        wb = openpyxl.load_workbook(io.BytesIO(result))
        assert "QC_Report" in wb.sheetnames
        ws = wb["QC_Report"]
        assert ws.cell(1, 1).value == "Metric"  # Header row populated
        assert ws.cell(2, 2).value == 18  # total_wells

    def test_export_with_excluded_wells_adds_replicate_fc_sheet(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from qpcr.export import export_to_excel
        from qpcr.analysis import AnalysisEngine
        import openpyxl
        import io

        processed = {"COL1A1": AnalysisEngine.calculate_ddct(
            sample_qpcr_raw_data, "GAPDH", "Non-treated", set(), set(), sample_mapping
        )}
        params = {"Housekeeping_Gene": "GAPDH", "Reference_Sample": "Non-treated"}

        result = export_to_excel(
            sample_qpcr_raw_data, processed, params, sample_mapping,
            excluded_wells=set(),
        )
        wb = openpyxl.load_workbook(io.BytesIO(result))
        assert "Replicate_FC" in wb.sheetnames
        ws = wb["Replicate_FC"]
        assert ws.cell(1, 1).value == "Target"
        assert ws.cell(1, 4).value == "Replicate_FC"
        # Should have 9 data rows (3 conditions × 3 replicates for COL1A1)
        assert ws.max_row >= 10  # header + 9 data rows
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m pytest tests/test_quality_control.py::TestExcelExportFixes -v`
Expected: FAIL (no Replicate_FC sheet, or QC_Report empty)

- [ ] **Step 3: Add `excluded_wells` param and Replicate_FC sheet to export.py**

In `qpcr/export.py`, modify `export_to_excel` signature (add `excluded_wells=None`):

```python
def export_to_excel(
    raw_data: pd.DataFrame,
    processed_data: Dict[str, pd.DataFrame],
    params: dict,
    mapping: dict,
    qc_stats: dict = None,
    replicate_stats: pd.DataFrame = None,
    excluded_wells=None,
) -> bytes:
```

After the FC_Matrix block (~line 112) and before the QC_Report call (~line 115), add:

```python
        # Replicate-level fold changes
        if raw_data is not None:
            hk_gene = params.get("Housekeeping_Gene")
            ref_sample = params.get("Reference_Sample")
            if hk_gene and ref_sample:
                try:
                    from qpcr.analysis import AnalysisEngine
                    replicate_fc = AnalysisEngine.compute_replicate_fold_changes(
                        raw_data=raw_data,
                        hk_gene=hk_gene,
                        ref_sample=ref_sample,
                        sample_mapping=mapping,
                        excluded_wells=excluded_wells,
                    )
                    if not replicate_fc.empty:
                        replicate_fc.to_excel(writer, sheet_name="Replicate_FC", index=False)
                except Exception:
                    pass  # Best-effort; don't break export if replicate FC fails
```

- [ ] **Step 4: Update all 3 call sites in main file**

Create a helper snippet to compute QC stats (avoid repeating at each call site). Add a helper function near the export section:

```python
def _build_export_extras():
    """Build qc_stats, replicate_stats, and excluded_wells for export."""
    qc_stats = None
    replicate_stats_df = None
    excl = st.session_state.get("excluded_wells", {})

    if st.session_state.get("data") is not None:
        from qpcr.quality_control import QualityControl

        # Flatten excluded_wells dict to flat set for QC functions that expect set
        excl_flat = set()
        if isinstance(excl, dict):
            for well_set in excl.values():
                excl_flat.update(well_set)
        elif excl:
            excl_flat = set(excl)

        # Pre-filter data for replicate stats (get_replicate_stats takes no exclusion param)
        qc_data = st.session_state.data.copy()
        if excl_flat:
            qc_data = qc_data[~qc_data["Well"].isin(excl_flat)]

        replicate_stats_df = QualityControl.get_replicate_stats(qc_data)

        # get_qc_summary_stats accepts a flat set
        qc_summary = QualityControl.get_qc_summary_stats(
            st.session_state.data, excl_flat or None
        )
        qc_stats = qc_summary if qc_summary else None

    return qc_stats, replicate_stats_df, excl
```

Then at each of the 3 call sites (~lines 6070, 6151, 6429), change from:
```python
export_to_excel(
    st.session_state.data,
    st.session_state.processed_data,
    analysis_params,
    st.session_state.sample_mapping,
)
```
to:
```python
_qc, _rep, _excl = _build_export_extras()
export_to_excel(
    st.session_state.data,
    st.session_state.processed_data,
    analysis_params,
    st.session_state.sample_mapping,
    qc_stats=_qc,
    replicate_stats=_rep,
    excluded_wells=_excl,
)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `python3 -m pytest tests/test_quality_control.py::TestExcelExportFixes -v`
Expected: ALL PASS

- [ ] **Step 6: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 7: Commit**

```bash
git add qpcr/export.py "streamlit qpcr analysis v1.py" tests/test_quality_control.py
git commit -m "feat: fix Excel QC_Report and add Replicate_FC sheet"
```

---

### Task 5: Final Integration + Full Test Run

**Files:** All modified files

- [ ] **Step 1: Run full test suite**

Run: `python3 -m pytest tests/ -v`
Expected: ALL PASS

- [ ] **Step 2: Verify Streamlit app manually**

Run: `streamlit run "streamlit qpcr analysis v1.py"`

Verify:
1. QC tab: SD threshold slider + "Find All Outliers" button shows preview table
2. QC tab: Per-gene "Clean" button works when a gene is selected
3. QC tab: "Apply Selected Exclusions" adds wells to excluded_wells
4. Mapping tab: "Finalize & Run Analysis" button works
5. Mapping tab: "Edit Mapping" button appears after finalize, re-enables editing
6. Analysis tab: stale badge appears when mapping is edited after analysis
7. Analysis tab: "Re-run Analysis" button works from stale badge
8. Graphs tab: stale warning appears
9. Export tab: Excel download has populated QC_Report sheet
10. Export tab: Excel download has Replicate_FC sheet with per-well fold changes

- [ ] **Step 3: Final commit**

```bash
git add qpcr/ tests/ "streamlit qpcr analysis v1.py"
git commit -m "feat: QC auto-exclude, reversible mapping, Excel export fixes"
```

**Note:** Do NOT use `git add -A` — the repo has untracked files that should not be committed.
