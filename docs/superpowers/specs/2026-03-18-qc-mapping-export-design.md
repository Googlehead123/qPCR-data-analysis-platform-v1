# QC, Mapping & Export Improvements — Design Spec

**Date:** 2026-03-18
**Author:** Min + Claude
**Status:** Approved
**Scope:** QC auto-exclusion, reversible mapping, Excel export fixes

---

## Context

Following the graph improvements round, three other platform areas need attention:
1. **QC tab** — well exclusion is tedious for high SD triplicates; no targeted auto-exclusion
2. **Mapping tab** — "Finalize Mapping" is irreversible, making iteration painful
3. **Excel export** — QC_Report sheet is always empty; no replicate-level data

### Goals
- Reduce QC manual work with smart, preview-based auto-exclusion
- Make mapping iteration painless without losing analysis safety
- Fix Excel exports to include QC stats and replicate fold changes

### Non-Goals
- Full QC redesign or visual overhaul
- Mapping presets/templates (deferred)
- Graph settings metadata in Excel (keep it lean)

---

## Deliverables

| # | Feature | Files | Effort |
|---|---------|-------|--------|
| 1 | QC auto-exclude high SD outliers (global + per-gene) | main UI, qpcr/quality_control.py | Medium |
| 2 | Reversible Finalize Mapping + "Analysis stale" badge | main UI | Medium |
| 3 | Excel: fix empty QC_Report sheet | main UI (call sites only) | Simple |
| 4 | Excel: add Replicate_FC sheet | qpcr/export.py | Simple |

---

## 1. QC Auto-Exclude High SD Outliers

### 1a. New Method: `QualityControl.find_high_sd_outliers()`

```python
@staticmethod
def find_high_sd_outliers(
    data: pd.DataFrame,
    excluded_wells: dict,
    sd_threshold: float = 0.5,
    gene_filter: str = None,
) -> list[dict]:
    """Find the worst replicate in each gene-sample group exceeding SD threshold.

    For each gene-sample group with SD > sd_threshold:
    - Identifies the replicate furthest from the group mean
    - Returns suggestion (not applied) for user review

    Args:
        data: Raw CT data with Well, Sample, Target, CT columns
        excluded_wells: Dict of {(gene, sample): set(well_ids)} already excluded
        sd_threshold: CT SD threshold (default 0.5)
        gene_filter: If set, only analyze this gene (for per-gene button)

    Returns:
        List of dicts: [{
            "Target": str,
            "Sample": str,
            "Well": str,
            "CT": float,
            "deviation": float,  # CT - group_mean (signed)
            "group_sd": float,
            "group_mean": float,
            "n_replicates": int,
        }, ...]
    """
```

**Logic:**
1. Group raw data by (Target, Sample)
2. Apply exclusions: `excluded_wells` is a dict `{(gene, sample): set(well_ids)}`. For each group, filter out wells in the corresponding exclusion set. If `excluded_wells` is a flat `set` instead, filter globally (backward compat).
3. For each group with n >= 3 remaining replicates: compute SD
4. If SD > sd_threshold: find the well with max |CT - mean|. On ties, take first encountered.
5. Return that well as a suggestion
6. If gene_filter is set, only process groups where Target == gene_filter

**Edge cases:**
- Groups with n < 3 (after exclusion): skip (can't meaningfully compute SD with n=2; users should handle those manually)
- Groups with all wells already excluded: skip
- HK gene groups: include (HK SD matters too)
- `excluded_wells` format: accept both dict `{(gene, sample): set}` and flat `set` (normalize internally, matching `AnalysisEngine.calculate_ddct` pattern)

### 1b. Global Button: "Auto-exclude High SD Outliers"

**Placement:** Top of QC Sub-Tab 1 (Triplicate Browser), before the per-gene/sample filters.

**UI Layout:**
```
[SD Threshold: [0.5] slider 0.3-1.0] [🔍 Find High SD Outliers] button
```

**When clicked:**
1. Calls `QualityControl.find_high_sd_outliers(data, excluded_wells, sd_threshold)`
2. If no outliers found: `st.info("No triplicates exceed SD threshold.")`
3. If outliers found: shows preview table with checkboxes

**Preview table columns:**
```
☑ | Gene | Sample | Well | CT | Deviation | Group SD | Replicates
```

- All rows checked by default
- User unchecks rows they want to keep
- "Apply Selected Exclusions" button executes checked rows
- Applied exclusions go through existing `excluded_wells` dict and undo history

### 1c. Per-Gene Button: "Clean Triplicates"

**Placement:** Next to the gene filter dropdown in QC Sub-Tab 1.

**Behavior:** Same as global button but with `gene_filter=current_gene`. Same preview table, same confirm flow. Only shows when a specific gene is selected (not "All Genes").

### Session State
- `st.session_state.sd_threshold` — float, default 0.5 (persists slider value)
- No new persistent state — suggestions are ephemeral, applied through existing exclusion mechanism

---

## 2. Reversible Finalize Mapping

### State Changes

**New session state keys:**
- `st.session_state.analysis_stale` — bool, default False
- Set True when: mapping edited after finalize, efficacy type changed, comparison conditions changed
- Set False when: analysis completes successfully

### Mapping Tab UI Changes

**Before finalize (current behavior, unchanged):**
- All mapping controls editable
- "Finalize & Run Analysis" button at bottom (renamed from "Finalize Mapping" — combines finalize + analysis trigger)

**After finalize (new behavior):**
- Currently: all controls disabled, no way back
- New: "Edit Mapping" button appears prominently
- Clicking "Edit Mapping":
  - Sets `mapping_finalized = False`
  - Sets `analysis_stale = True`
  - Re-enables all mapping controls with current values preserved
  - Does NOT clear analysis results or graphs

### "Analysis Stale" Badge

**Appears on:** Analysis tab and Graphs tab headers, when `analysis_stale == True`

**Display:**
```python
if st.session_state.get("analysis_stale", False):
    st.warning("⚠️ Mapping changed since last analysis. Results may be outdated.")
    if st.button("Re-run Analysis", key="rerun_stale"):
        # Trigger run_full_analysis() with current mapping
        # Set analysis_stale = False on success
```

**Triggers for stale flag:**
- User clicks "Edit Mapping" after finalize
- User changes efficacy type dropdown
- User changes any comparison condition selector
- User modifies condition names or groups after prior analysis

**Clears stale flag:**
- Successful `run_full_analysis()` completion

### Behavior Guarantees
- Old analysis results and graphs remain visible (with stale badge) until re-run
- Re-running analysis uses the CURRENT mapping state (whatever user edited)
- Undo history for QC exclusions is unaffected by mapping changes
- Changing efficacy type preserves custom condition names (only resets group suggestions)

---

## 3. Excel: Fix Empty QC_Report Sheet

### Problem
`export_to_excel()` accepts `qc_stats=None` parameter but all 3 call sites in the main file never pass it. The `QC_Report` sheet exists in the code but is always empty.

### Fix

**At each call site** (3 locations in main file: standalone Excel download, ZIP export, Complete Report ZIP), compute and pass QC stats:

**IMPORTANT:** `QualityControl.get_replicate_stats(data)` accepts ONLY a DataFrame — no `excluded_wells` parameter. Pre-filter the data before passing:

```python
qc_stats = None
replicate_stats_df = None
if st.session_state.get("data") is not None:
    from qpcr.quality_control import QualityControl

    # Pre-filter data to exclude wells before computing stats
    qc_data = st.session_state.data.copy()
    excl = st.session_state.get("excluded_wells", {})
    if isinstance(excl, dict):
        # Dict format: {(gene, sample): set(well_ids)}
        excl_flat = set()
        for well_set in excl.values():
            excl_flat.update(well_set)
        qc_data = qc_data[~qc_data["Well"].isin(excl_flat)]
    elif excl:
        qc_data = qc_data[~qc_data["Well"].isin(excl)]

    replicate_stats_df = QualityControl.get_replicate_stats(qc_data)

    qc_stats = {
        "replicate_stats": replicate_stats_df,
        "excluded_wells": excl,
        "excluded_count": sum(
            len(ws) for ws in excl.values()
        ) if isinstance(excl, dict) else len(excl or set()),
    }
```

Then pass `qc_stats=qc_stats` and `replicate_stats=replicate_stats_df` to `export_to_excel()`.

**No changes to `get_replicate_stats()` needed** — we pre-filter the data instead.

---

## 4. Excel: Add Replicate_FC Sheet

### New Sheet: "Replicate_FC"

**Contents:** Per-replicate fold change values — one row per well per gene.

| Target | Condition | Well | Replicate_FC |
|--------|-----------|------|-------------|
| COL1A1 | Non-treated | A4 | 0.97 |
| COL1A1 | Non-treated | A5 | 1.02 |
| COL1A1 | Non-treated | A6 | 1.01 |
| COL1A1 | Treatment1 | B4 | 2.41 |
| ... | ... | ... | ... |

### Implementation

In `qpcr/export.py`, after the FC_Matrix sheet block, add:

```python
# Replicate-level fold changes
if raw_data is not None:
    hk_gene = params.get("Housekeeping_Gene")
    ref_sample = params.get("Reference_Sample")
    if hk_gene and ref_sample:
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
```

### Parameter Change

`export_to_excel()` needs access to `excluded_wells` for the replicate computation. Add it as a new parameter:

```python
def export_to_excel(
    raw_data, processed_data, params, mapping,
    qc_stats=None, replicate_stats=None,
    excluded_wells=None,  # NEW
) -> bytes:  # NOTE: returns bytes via output.getvalue(), not BytesIO
```

All call sites already have `excluded_wells` available in session state — pass it through.

---

## Testing Strategy

### New Tests
- `find_high_sd_outliers`: returns correct suggestions for known high-SD data, respects gene_filter, skips groups with n<3
- `Replicate_FC` sheet: verify it appears in exported Excel, contains expected row count
- QC_Report sheet: verify it's populated (not empty) when qc_stats passed

### Existing Tests
- All 109 existing tests must continue passing
- No changes to graph tests

---

## Migration / Backward Compatibility

- "Finalize Mapping" renamed to "Finalize & Run Analysis" (cosmetic)
- "Edit Mapping" is additive — no existing flow changes
- Excel export gains 1-2 new sheets — existing sheets unchanged
- `export_to_excel()` new parameter has default `None` — existing callers unaffected
- QC auto-exclude is opt-in (button click) — no behavior change without user action
