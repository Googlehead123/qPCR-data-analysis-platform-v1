# qPCR Data Analysis Platform v1

## Overview
Comprehensive qPCR (quantitative PCR) data analysis platform for cosmetics/dermatology efficacy evaluation. Supports **21 efficacy assay items** (효능평가항목 catalog: anti-aging, whitening, hydration/barrier, soothing, sebum/acne, hair, pet, lip, pores, …) with automated ΔΔCt calculations, automatic best-2-of-3 QC, a Results Overview, and Excel/PowerPoint reports.

## Tech Stack
- **Framework:** Streamlit (Python 3.12)
- **Visualization:** Plotly (UI font: Pretendard via CDN)
- **Reports:** python-pptx (PowerPoint), kaleido + Chromium (image export)
- **Data:** Pandas, NumPy, SciPy

## How to Run
```bash
streamlit run "streamlit qpcr analysis v1.py"    # http://localhost:8501
python3 -m pytest                                 # test suite (~170 tests); note: no `python` on PATH
```

## Architecture
The app is a Streamlit UI shell (`streamlit qpcr analysis v1.py`) that **imports its logic from the `qpcr/` package** — the package is the single source of truth for computation and export:
- `qpcr/parser.py` — `QPCRParser` (raw instrument CSV parsing)
- `qpcr/quality_control.py` — `QualityControl` (QC, auto best-2-of-3 replicate selection)
- `qpcr/analysis.py` — `AnalysisEngine` core (ΔΔCt / stats); the monolith subclasses it for Streamlit orchestration
- `qpcr/graph.py` — `GraphGenerator.create_gene_graph` (Plotly charts)
- `qpcr/report.py` — `ReportGenerator` (chart images) + `PPTGenerator` (decks)
- `qpcr/export.py` — `export_to_excel` (+ native Excel chart post-processing)
- `qpcr/export_utils.py` — `export_figure_to_bytes` (route ALL image export here: Kaleido + `BROWSER_PATH`), `build_zip`
- `qpcr/constants.py` — `EFFICACY_CONFIG` (the 21-item catalog — **single definition**, imported by the app), presets, thresholds
- `qpcr/auto/` — deterministic screening / interpretation helpers

The monolith holds the UI (tabs, widgets, session-state orchestration) and a thin `AnalysisEngine` subclass. Prefer editing the `qpcr/` package; the monolith mostly wires it into the UI.

## Development
- **Branch from `origin/main`, not local `main`.** Local `main` can lag origin (it was once ~21 commits behind); always `git fetch && git checkout -b <branch> origin/main` so you build on shipped code.
- **CI** (`.github/workflows/ci.yml`) runs `pytest` on every push/PR — keep it green before merging.
- **Runtime check, not just units:** `tests/test_app_smoke.py` drives the whole Streamlit pipeline headless via `AppTest` (boot → QC → mapping → analysis → Overview → graphs → export). Run it after any UI/session-state change.
- Deploy: Streamlit Cloud auto-redeploys from `main`. Deploy config: `runtime.txt` (python-3.12), `packages.txt` (chromium + fonts-noto-cjk + lxml), pinned `requirements.txt`.

## GitHub
Repository: Googlehead123/qPCR-data-analysis-platform-v1

## Notes
- Validation shows <0.001% ΔΔCt error vs reference calculations.
- Uses `streamlit-sortables` for drag-and-drop sample ordering.
- Korean localization throughout (efficacy names, report labels).
