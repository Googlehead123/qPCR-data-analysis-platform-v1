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

**Rule:** requirements.txt must pin exact versions for all packages. Check that numpy lower bound satisfies scipy's requirement (`>=1.26.4`).

**Superseded 2026-08-07:** the `runtime.txt` claim below is wrong — see the next entry. `runtime.txt` does NOT control the Python version on Streamlit Community Cloud.

## 2026-08-07: `numpy<2.0` is unsatisfiable on Python 3.13 (Cloud boot failure)

**What went wrong:** The Cloud app never reached startup. Logs ended at `Resolved 63 packages in 1.61s` and never printed an install completion or app boot.

Root cause: `requirements.txt` pinned `numpy>=1.26.4,<2.0`, but Cloud provisioned **Python 3.13.14**. numpy jumped `1.26.4` → `2.0.0` — there is no 1.27/1.28/1.29 — and **no numpy 1.26.x ships a cp313 wheel** (tags are cp39/cp310/cp311/cp312/pp39 only). So the `<2.0` ceiling left exactly one candidate, 1.26.4, with no usable wheel, forcing a from-source numpy build (meson + C/Fortran toolchain) that the container could not complete.

Reproduce the failure without deploying:
```bash
uv pip install --dry-run --only-binary :all: --python-version 3.13 \
  --python-platform x86_64-unknown-linux-gnu "numpy>=1.26.4,<2.0"
# error: numpy==1.26.4 has no usable wheels ... requirements are unsatisfiable
```

**Rule — `runtime.txt` is NOT honored by Streamlit Community Cloud.** It is a Heroku convention. Cloud takes the Python version from the app's *Advanced settings* at creation time, and it cannot be changed in-place afterward. `runtime.txt` said `python-3.12` while production ran 3.13.14 — that mismatch is what let the bad pin ship. Treat `runtime.txt` as documentation only; verify the real version in the deploy log line `Using Python X.Y.Z environment at /home/adminuser/venv`.

**Rule — every dependency pin must have a wheel for the Python version Cloud *actually* runs, not the one you intended.** Before changing any pin, run the `--only-binary :all: --python-version <cloud-version>` dry-run above. A pin that resolves is not the same as a pin that installs: resolution succeeded here in 1.61s and the build still killed the container.

**Rule — prefer pins that carry wheels for both the local and Cloud interpreter.** `numpy==2.4.1` has cp312 and cp313 wheels, so local (3.12) and Cloud (3.13) run identical bytes. Local had silently drifted to numpy 2.4.1 while requirements.txt still said `<2.0`, so the pin was untested in dev — that drift hid the problem.

**Rule — `chromium-driver` is not needed for kaleido image export.** kaleido v1 / `choreographer` drives Chrome directly over CDP; `chromium-driver` is Selenium's chromedriver. `packages.txt` was installing 299 apt packages / 323 MB download / 1119 MB disk to get a browser; dropping `chromium-driver` trims that with no export impact.
