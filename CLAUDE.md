# qPCR Data Analysis Platform v1

## Overview
Comprehensive qPCR (quantitative PCR) data analysis platform for cosmetics/dermatology efficacy evaluation. Supports 10 Korean efficacy categories with automated DDCt calculations, quality control, and report generation.

## Tech Stack
- **Framework:** Streamlit
- **Language:** Python 3.12
- **Visualization:** Plotly, Matplotlib, Seaborn
- **Reports:** python-pptx (PowerPoint), kaleido (image export)
- **Data:** Pandas, NumPy, SciPy

## How to Run
```bash
streamlit run "streamlit qpcr analysis v1.py"    # http://localhost:8501
pytest tests/                                      # Run 80+ tests
```

## Key Files
- `streamlit qpcr analysis v1.py` — the Streamlit app (~5.7k lines). Holds the UI
  for all seven tabs plus the report/export writers that live nowhere else:
  - `ReportGenerator` / `PPTGenerator` — PowerPoint generation
  - `export_to_excel` — workbook + native Excel charts
  - `AnalysisEngine` — a thin subclass of `qpcr.analysis.AnalysisEngine` adding
    only `run_full_analysis` (the Streamlit orchestration); the ΔΔCt maths is in
    the package, not here.
- `qpcr/` — the computational core, and the single definition of each of:
  `parser.py` (QPCRParser), `quality_control.py` (QualityControl), `analysis.py`
  (ΔΔCt + statistics), `graph.py` (GraphGenerator), `constants.py`
  (EFFICACY_CONFIG, presets, fonts), `export_utils.py` (Kaleido image export),
  `auto/` (screening, test advice, MIQE).
  **Do not re-declare any of these in the monolith** — that dual-copy hazard is
  documented in `tasks/lessons.md` and has bitten this project repeatedly.
- `tests/` — 19 modules, 256 tests. Note `conftest.py` replaces `streamlit` with
  a MagicMock, so widget behaviour is invisible to most of the suite; the tests
  that need real Streamlit drive the app through `AppTest` **in a subprocess**
  (see `tests/test_gene_editor_state.py`).
- `requirements.txt` — exact pins; see `tasks/lessons.md` on Cloud Python drift.

## Korean Efficacy Categories (21)
Generated from `qpcr/constants.py::EFFICACY_CONFIG` — that dict is the source of
truth; regenerate this list rather than editing it by hand:

탄력, 광노화, 보습/수분, 장벽, 속보습, 멜라닌 생성, 진정, 가려움 개선, 냉감, 열감,
여드름, 과각화, 활력, 탈모 개선, 모공 탄력, 열 노화, 민감성 기전, 외이도염,
구강 개선, 선번 완화, 립 색상

Each entry carries `genes` (with spelling variants), `cell`, `treatment_time`,
`controls` (`negative` / optional `positive` / `compare_to`, plus display-only
`negative_display`) and `expected_direction` per marker.

## GitHub
Repository: Googlehead123/qPCR-data-analysis-platform-v1

## Notes
- **Architecture:** partially modularized. The computational layer lives in
  `qpcr/`; the app file keeps the UI plus the report/export writers. Remaining
  refactor target is those writers.
- Validation shows <0.001% DDCt error vs reference calculations. Independently
  re-verified 2026-08-24: ΔΔCt, the asymmetric fold-change error transform and
  the Benjamini-Hochberg implementation all match hand/reference computation.
- Significance markers shown on charts, slides and summary sheets are
  **uncorrected** p-values. BH q-values are computed and reach the per-gene
  Excel sheets, but nothing reads them for display — the provenance record says
  so explicitly.
- Error-bar spread is computed from the **target** replicates only; housekeeping
  variability is not propagated into it. All error-bar modes are plotted in the
  fold-change domain.
- The t-test's n is **technical wells**, not biological samples — see the open
  question in `tasks/lessons.md`.
- Uses streamlit-sortables for drag-and-drop sample ordering
- python-pptx for automated PowerPoint report generation
