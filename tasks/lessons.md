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
