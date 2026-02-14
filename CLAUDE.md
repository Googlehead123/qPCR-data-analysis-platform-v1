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
- `streamlit qpcr analysis v1.py` — **Single monolithic file (5,472 lines)** containing:
  - `QPCRParser` — Raw data parsing
  - `QualityControl` — QC checks
  - `AnalysisEngine` — DDCt calculations
  - `GraphGenerator` — Plotly visualizations
  - `ReportGenerator` — Summary reports
  - `PPTGenerator` — PowerPoint generation
- `tests/` — 7 test modules (80+ tests)
- `requirements.txt` — Dependencies

## Korean Efficacy Categories (10)
탄력 (Elasticity), 항노화 (Anti-aging), 보습 (Moisturizing), 미백 (Whitening), 진정 (Soothing), 장벽 (Barrier), 모공 (Pores), 주름 (Wrinkles), 각질 (Keratin), 피지 (Sebum)

## GitHub
Repository: Googlehead123/qPCR-data-analysis-platform-v1

## Notes
- **Architecture:** Monolithic single file — candidate for refactoring into modules
- Validation shows <0.001% DDCt error vs reference calculations
- Uses streamlit-sortables for drag-and-drop group ordering
- python-pptx for automated PowerPoint report generation
