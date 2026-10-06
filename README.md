# qPCR Data Analysis Platform v1

A Streamlit app for analysing qPCR data from cosmetics / dermatology efficacy
studies. It parses instrument exports, runs quality control, computes ΔΔCt
fold changes with statistics, draws Plotly charts, and exports Excel workbooks
and PowerPoint reports. Efficacy presets cover 21 Korean efficacy categories
(`qpcr/constants.py::EFFICACY_CONFIG`).

## Run

```bash
pip install -r requirements.txt
streamlit run "streamlit qpcr analysis v1.py"    # http://localhost:8501
```

The filename contains spaces, so keep the quotes. `packages.txt` lists the
system packages Streamlit Cloud installs (Chromium and CJK/emoji fonts, used for
Kaleido image export and Korean text).

## Test

```bash
pip install -r requirements.txt -r requirements-dev.txt
pytest tests/
```

CI (`.github/workflows/ci.yml`) runs the suite on Python 3.12 and 3.13.

## Layout

- `streamlit qpcr analysis v1.py` — the app: UI plus the PowerPoint and Excel
  report writers.
- `qpcr/` — the computational core (parser, quality control, ΔΔCt and
  statistics, graphs, constants, image export, auto-analysis).
- `tests/` — pytest suite.
- `tasks/lessons.md` — lessons learned and design decisions.

## Dependencies

Edit `requirements.in` (direct runtime dependencies) only. `requirements.txt`
is generated from it and pins every transitive package; do not hand-edit it.
Regenerate with:

```bash
uv pip compile requirements.in --universal --python-version 3.12 \
    --output-file requirements.txt
```

then verify the install and tests on both Python 3.12 and 3.13 (recipe in the
header of `requirements.in`). Test-only pins live in `requirements-dev.txt`.

## Deploy

Deployed on Streamlit Community Cloud, which installs from `requirements.txt`.
The Python version (3.13 in production) is chosen in the app's deploy dialog,
Advanced settings — no file in the repo controls it.

## More

- [CLAUDE.md](CLAUDE.md) — authoritative project notes: architecture, key
  files, statistical conventions.
- [tasks/lessons.md](tasks/lessons.md) — lessons learned.
- [AGENTS.md](AGENTS.md) — conventions for AI coding agents.
