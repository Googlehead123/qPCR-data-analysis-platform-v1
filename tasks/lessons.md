# Lessons Learned

## 2026-03-18: Monolith + Package Dual-Class Architecture

**What went wrong:** Applied changes only to `qpcr/` package modules (graph.py, quality_control.py) but the Streamlit app uses INLINE copies of those classes defined in the monolith file (`streamlit qpcr analysis v1.py`). Tests passed because they import from `qpcr/` directly, masking the runtime failure.

**Rule:** ANY change to a class in `qpcr/*.py` MUST also be applied to the corresponding inline class in the monolith file. The classes exist in both locations:
- `qpcr/graph.py` ↔ monolith line ~2028 `class GraphGenerator`
- `qpcr/quality_control.py` ↔ monolith line ~750 `class QualityControl`
- `qpcr/analysis.py` ↔ monolith line ~1383 `class AnalysisEngine`
- `qpcr/report.py` ↔ monolith line ~2497/2938 `class ReportGenerator`/`PPTGenerator`
- `qpcr/export.py` ↔ monolith line ~3429 `def export_to_excel`

**Prevention:** After modifying any `qpcr/*.py` file, always verify the Streamlit app uses the package version OR sync the monolith copy. Long-term fix: refactor the monolith to import from `qpcr/` instead of defining its own copies.

**Hit again on 2026-03-18:** `AnalysisEngine.compute_replicate_fold_changes` was added to `qpcr/analysis.py` but not the monolith. Data points toggle crashed at runtime. Fix: added method to monolith's AnalysisEngine at line ~1909.

## 2026-04-01: Streamlit Cloud vs Local Environment Drift

**What went wrong:** Three separate bugs caused Streamlit Cloud failures that didn't appear locally:

1. **Wrong env var for kaleido Chrome path** — `_fig_to_image` in both monolith (~line 2700) and `qpcr/report.py` set `os.environ["CHROME_PATH"]` but `choreographer` (kaleido's browser driver) reads `os.environ["BROWSER_PATH"]`. Image/PPT export would fail on Cloud because the Chromium path override was silently ignored.

2. **numpy constraint too loose + conflicting with scipy** — requirements.txt had `numpy>=1.24.0,<2.0` but scipy 1.17.0 requires `numpy>=1.26.4`. Pip's resolver handles this but it produced warnings and was fragile. Combined with all other packages being unpinned, Cloud got unpredictable versions.

3. **No runtime.txt** — Streamlit Cloud didn't know to use Python 3.12.

**Rule:** Any time kaleido/plotly image export is used, the env var is `BROWSER_PATH`, not `CHROME_PATH`. Always guard with `if not os.environ.get("BROWSER_PATH")` to avoid overriding a user-set value.

**Rule:** requirements.txt must pin exact versions for all packages. Check that numpy lower bound satisfies scipy's requirement (`>=1.26.4`). Always include `runtime.txt` for Streamlit Cloud Python version pinning.
