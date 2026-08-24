# Lessons Learned

## 2026-03-18: Monolith + Package Dual-Class Architecture

> **Partly superseded 2026-08-24.** The line numbers and the "inline copies
> exist" claim below are no longer true, but the RULE still is — keep reading it
> as a hazard, not as a map. Verified current state: `QPCRParser`,
> `QualityControl` and `GraphGenerator` are imported from `qpcr/` and defined
> nowhere else; the monolith's `AnalysisEngine` *subclasses*
> `qpcr.analysis.AnalysisEngine` and adds only `run_full_analysis`, so the ΔΔCt
> maths has one home; `qpcr/report.py` and `qpcr/export.py` were deleted, so
> `ReportGenerator` / `PPTGenerator` / `export_to_excel` live only in the
> monolith. The last real duplicates — `PLOTLY_FONT_FAMILY`, `CM_TO_PX`,
> `CM_TO_EMU`, declared byte-identically in both places with `qpcr/graph.py`
> reading one copy and the PPT/report writers the other — were removed the same
> day. `natural_sort_key` is still duplicated (monolith + `qpcr/utils.py`),
> byte-identical, and `tests/test_package.py` compares them.

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

## 2026-08-21: A keyed Streamlit widget ignores a recomputed `value=` / `index=`

**What went wrong:** Four separate defects in the per-gene graph editor, one root cause. "Size preset" never applied (every preset resolved to 28x16 and reported itself as "Custom"), a hand-picked bar colour was silently repainted to the preset tone on the next interaction, "Reset this gene" did nothing, and adjusting Width then switching gene crashed the Graphs tab with `StreamlitAPIException`.

Root cause: the editor was written assuming that writing `graph_settings` would make the widget follow. It never does. Verified in the installed streamlit 1.53 source — `compute_and_register_element_id` in `streamlit/elements/lib/utils.py` takes a `key_as_main_identity` argument, and when a `key=` is supplied the element id is derived from the key alone (`st.slider` additionally from min/max/step, `st.select_slider` from options, `st.selectbox` from accept_new_options). Because the id does not change when `value=`/`index=` changes, the widget keeps its stored state and the recomputed argument is dead after the first render. Confirmed empirically for slider, number_input, selectbox, select_slider, color_picker, radio, checkbox and toggle.

**Rule — the widget key is the source of truth; `value=`/`index=` is a first-render default only.** Any code that writes a settings dict and expects a widget to display or honor the new value is a silent no-op. This includes: applying a preset, a "Reset"/"Clear"/"Select all" button, and clamping a stored value back into a valid range.

**Rule — to move a keyed widget from code, assign `st.session_state[widget_key]` BEFORE the widget is instantiated.** Streamlit raises `StreamlitAPIException` if you assign it after (`session_state.py` guards on `ctx.widget_ids_this_run`). If the trigger sits below the widget in the script, set a flag and do the work at the top of the next run — that is why "Reset this gene" raises a flag consumed at the top of `render_gene_editor` instead of clearing inline.

**Rule — a two-way binding needs a sentinel.** Forcing the widget from settings unconditionally destroys the user's own selection. `_settings_backed_selectbox` remembers the last value it pushed and only pushes when the settings value has changed since then, which is exactly the case where something other than the widget changed it. Use it for any selectbox whose value code needs to set, and for any selectbox whose OPTIONS can change under it (a stored value that vanishes from the option list will not re-seed itself).

**Rule — a "Reset" must clear the widget keys, not just the settings.** Popping `graph_settings` alone leaves every widget holding its old value, and they repopulate the settings on the next run.

**Rule — never call `st.rerun(scope="fragment")` where the fragment body can run during a full script run.** Streamlit rejects it outright (`_new_fragment_id_queue` in `commands/execution_control.py` raises when `ctx.fragment_ids_this_run` is empty), and a fragment body executes as part of every full rerun. Use a plain `st.rerun()` unless you are certain the code path is only ever reached from a widget interaction inside the fragment.

**Rule — the repo-root `conftest.py` makes widget behaviour invisible to the test suite.** It replaces `sys.modules["streamlit"]` with a `MagicMock` whose widgets return canned values, so no test could observe a widget ignoring its argument — all four defects above shipped with 250 tests passing. To test real widget behaviour, drive the app with `streamlit.testing.v1.AppTest` **in a subprocess** (the mock is installed at conftest import time and would otherwise defeat it). `AppTest.from_file("streamlit qpcr analysis v1.py")` works on the whole app: inject `data`, `raw_data`, `processed_data`, `selected_efficacy`, `hk_gene`, `analysis_ref_condition`, `analysis_cmp_condition`, `sample_mapping`, `sample_order`, `excluded_wells`, `excluded_samples`, `gene_display_names`, `selected_gene_idx` into `at.session_state`, then `.run()` lands you on a fully populated Graphs tab. See `tests/test_gene_editor_state.py`. Note `at.session_state` has no `.get()` — index it.

**Rule — figure-memo fingerprints must be stable AND minimal.** Two entries moved on their own: `bar_colors_per_sample` was absent from the fingerprint until some *other* gene created the dict and then appeared as `{}`, and the whole `<gene>_bar_settings` dict was hashed including each bar's `color`, which `graph.py` never reads from there. Always emit a key rather than omitting it conditionally, and hash only the fields the renderer actually consumes. Counting `create_gene_graph` calls under AppTest is the way to check this: steady state must be 0 rebuilds, and a one-gene edit exactly 1.

---

## 2026-08-24: Open questions for Min (found in the platform-wide audit)

Both change published numbers, so neither was changed unilaterally.

### 1. The t-test's n is technical wells, not biological samples

`qpcr/analysis.py::calculate_statistics` flattens every well of every Sample
mapped to a condition and feeds them to the t-test. Nothing aggregates to
biological-sample level first, so the degrees of freedom are inflated by the
technical replication.

Measured on a standard 3 biological x 3 technical design:

| n used | p-value | marker |
|---|---|---|
| 9 wells (current) | 2.6965e-06 | *** |
| 3 biological samples | 2.1229e-03 | ** |

Three orders of magnitude. A real effect survives either way, but a marginal one
is manufactured by the current choice. Note the pooling is deliberate elsewhere:
"Samples sharing a condition name are treated as biological replicates" is the
Mapping tab's stated contract, and the ΔΔCt bar means are computed the same way.

Suggested shape if Min wants it changed: a `replicate_unit` argument defaulting
to `"technical"` (so no published number moves without an explicit choice) plus a
radio where the dead SEM/SD one used to be. Until then `build_provenance`'s
`statistical_test` string should say which unit was used.

### 2. Which significance source should be displayed

BH q-values (`p_value_fdr` / `significance_fdr`) are computed correctly but have
no display consumer — every asterisk on every chart, slide, summary sheet,
Overview matrix and results table is an uncorrected p-value. The provenance
record and MIQE checklist now state this plainly rather than claiming BH.

The proper fix is a single `graph_settings["sig_source"]` toggle next to the
error-bar selector, read by `graph.py`, the PPT writer, the Overview matrix and
the results table. **Which way it should default is Min's call** — corrected is
the defensible default for a multi-gene panel, uncorrected is what every number
shipped so far has used, so flipping it silently would change conclusions in
past-comparable reports.

### 3. Not investigated

The whole-platform health agent (dependency pinning, Cloud-vs-local drift,
per-rerun cost) was stopped before reporting. `requirements.txt` pinning and the
Cloud Python-version trap are covered by the 2026-04-01 and 2026-08-07 entries
above; the rest of that sweep is still open.
