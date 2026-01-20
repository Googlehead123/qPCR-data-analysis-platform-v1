# AGENTS.md - qPCR Data Analysis Platform

This document provides essential information for AI coding agents working on this codebase.

## Project Overview

A Streamlit-based qPCR (quantitative PCR) data analysis application for gene expression analysis.
Primary use case: Cosmetics efficacy testing with gene-by-gene delta-delta Ct (DDCt) calculations.

**Tech Stack:**
- Python 3.x with Streamlit
- pandas, numpy, scipy for data processing
- plotly for interactive visualizations
- Single-file architecture: `streamlit qpcr analysis v1.py`

---

## Build/Run Commands

### Development Server
```bash
# Install dependencies
pip install -r requirements.txt

# Run the Streamlit application
streamlit run "streamlit qpcr analysis v1.py"

# Run on specific port
streamlit run "streamlit qpcr analysis v1.py" --server.port 8501
```

### Testing
```bash
# Run full test suite
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov="streamlit qpcr analysis v1" --cov-report=term-missing

# Run single test file
pytest tests/test_parser.py -v

# Run single test
pytest tests/test_parser.py::TestQPCRParserDetectFormat::test_detect_format1_with_well_position -v
```

### Linting (Recommended)
```bash
# No linter configured. Recommended setup:
# pip install ruff
# ruff check .
# ruff format .
```

---

## Code Style Guidelines

### File Structure
```
/
  streamlit qpcr analysis v1.py  # Main application (all code)
  requirements.txt               # Dependencies
  README.md                      # Project docs
  AGENTS.md                      # This file
```

### Import Organization
Standard library first, then third-party, then local. Group by category:
```python
# Standard library
import io
import json
import re
from datetime import datetime
from typing import Dict, List, Tuple

# Third-party - Data processing
import pandas as pd
import numpy as np
from scipy import stats

# Third-party - Visualization
import plotly.graph_objects as go
import plotly.express as px

# Third-party - Framework
import streamlit as st
```

### Naming Conventions

| Type | Convention | Example |
|------|------------|---------|
| Classes | PascalCase | `QPCRParser`, `AnalysisEngine`, `GraphGenerator` |
| Functions/Methods | snake_case | `calculate_ddct`, `parse_format1` |
| Constants | UPPER_SNAKE | `EFFICACY_CONFIG` |
| Variables | snake_case | `gene_data`, `sample_mapping` |
| Session state keys | snake_case strings | `'processed_data'`, `'sample_order'` |

### Section Comments
Use banner-style section delimiters:
```python
# ==================== SECTION NAME ====================
```

### Class Structure
Classes are organized by responsibility:
- `QPCRParser` - Data parsing and format detection
- `AnalysisEngine` - DDCt calculations and statistics
- `GraphGenerator` - Plotly visualization creation

Methods are typically `@staticmethod` when they don't need instance state.

### Type Annotations
Use typing module for function signatures:
```python
def calculate_ddct(data: pd.DataFrame, hk_gene: str, ref_sample: str, 
                   excluded_wells: set, sample_mapping: dict) -> pd.DataFrame:
```

### Error Handling
- Use try/except with specific error display via `st.error()`
- Return None on parse failures instead of raising
- Provide user-friendly error messages
```python
try:
    # operation
except Exception as e:
    st.error(f"Parse error: {e}")
    return None
```

### Streamlit Patterns

**Session State Initialization:**
```python
for key in ['data', 'processed_data', 'sample_mapping']:
    if key not in st.session_state:
        st.session_state[key] = {} if key in ['sample_mapping'] else None
```

**Tab-based Layout:**
```python
tab1, tab2, tab3 = st.tabs(["Tab 1", "Tab 2", "Tab 3"])
with tab1:
    # content
```

**Metrics Display:**
```python
col1, col2, col3 = st.columns(3)
col1.metric("Label", value)
```

**User Feedback:**
```python
st.success("Operation complete")
st.warning("Warning message")
st.error("Error message")
st.info("Information")
```

### DataFrame Operations
- Use `.copy()` when modifying DataFrames to avoid SettingWithCopyWarning
- Prefer `.query()` for filtering: `df.query('Sample.notna() & Target.notna()')`
- Use `pd.to_numeric(col, errors='coerce')` for safe numeric conversion

### Configuration Patterns
Use dictionaries for configuration with nested structure:
```python
EFFICACY_CONFIG = {
    'test_type': {
        'genes': ['GENE1', 'GENE2'],
        'cell': 'cell_type',
        'controls': {
            'negative': 'control_name',
            'positive': 'positive_name',
        },
        'description': 'Description text'
    }
}
```

---

## Architecture Notes

### Data Flow
1. **Upload** - CSV files parsed by `QPCRParser`
2. **Mapping** - Samples mapped to conditions via UI
3. **Analysis** - `AnalysisEngine.calculate_ddct()` + `calculate_statistics()`
4. **Visualization** - `GraphGenerator.create_gene_graph()` per gene
5. **Export** - Excel/HTML/JSON output

### Key Session State Variables
```python
st.session_state.data              # Raw parsed DataFrame
st.session_state.processed_data    # Dict[gene: DataFrame] with results
st.session_state.sample_mapping    # Dict[sample: {condition, group, include}]
st.session_state.sample_order      # List[str] sample display order
st.session_state.hk_gene           # Housekeeping gene name
st.session_state.excluded_wells    # Set of excluded well IDs
st.session_state.excluded_samples  # Set of excluded sample names
st.session_state.graphs            # Dict[gene: plotly.Figure]
st.session_state.graph_settings    # Graph customization settings
```

### Statistical Methods
- **DDCt**: `2^(-(target_ct - hk_ct) - (ref_target_ct - ref_hk_ct))`
- **P-values**: Welch's t-test (unequal variance) via `scipy.stats.ttest_ind`
- **Significance**: * p<0.05, ** p<0.01, *** p<0.001

---

## Common Tasks

### Adding a New Efficacy Type
Add to `EFFICACY_CONFIG` dictionary with required keys:
- `genes`: List of target gene names
- `cell`: Cell line description
- `controls`: Dict with control conditions
- `description`: User-facing description

### Modifying Graph Appearance
Edit `GraphGenerator.create_gene_graph()` method. Key parameters:
- Colors via `bar_colors` dict
- Error bars via `error_y` config
- Significance annotations via `fig.add_annotation()`

### Adding New Export Format
Add to TAB 5 section. Use `io.BytesIO()` or `io.StringIO()` for buffers.

---

## Known Patterns to Follow

1. **Filename has space** - Always quote: `"streamlit qpcr analysis v1.py"`
2. **Korean text** - Used in efficacy names (e.g., ''), handle with UTF-8
3. **Plotly graphs** - Use `go.Figure()` with `go.Bar()` traces
4. **Natural sorting** - Use regex-based key for sample ordering
5. **Dual p-value support** - Primary (*) and secondary (#) comparisons

---

## Do NOT

- Add type suppressions (`# type: ignore`)
- Use `st.experimental_*` deprecated APIs
- Modify session state directly in callbacks without `st.rerun()`
- Break the single-file architecture without discussion
- Remove Korean text labels (internationalization in progress)
