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

## 2026-08-24: Two statistical choices, now DECIDED (Min)

Both affect published numbers. Both were surfaced by the platform-wide audit and
decided deliberately rather than drifted into — record them here so a future
reader does not "fix" either one.

### 1. DECIDED: the t-test's n is technical wells — keep it that way

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

**Decision (Min, 2026-08-24): keep technical wells.** No code change to the
maths. Rationale for keeping it: the pooling is already the app's stated
contract — "Samples sharing a condition name are treated as biological
replicates" in the Mapping tab — and the ΔΔCt bar means are computed the same
way, so aggregating only inside the t-test would have made the bars and the
p-values describe different units.

What DID change: the choice is now declared instead of implied.
`build_provenance` reports `statistical_test` as "<Welch|Student> t-test on ΔCt,
n = technical replicate wells pooled per condition" and adds a `replicate_unit`
field, and the MIQE checklist prints "Replicate unit for the test" as its own
item — MIQE requires it, because the replicate unit is what sets the degrees of
freedom. A reviewer can now see the n's meaning without reading the source.

If this is ever revisited, the shape to use is a `replicate_unit` argument
defaulting to `"technical"` so no existing number moves silently.

### 2. DECIDED: display uncorrected p-values

BH q-values (`p_value_fdr` / `significance_fdr`) are computed correctly but have
no display consumer — every asterisk on every chart, slide, summary sheet,
Overview matrix and results table is an uncorrected p-value. The provenance
record and MIQE checklist now state this plainly rather than claiming BH.

**Decision (Min, 2026-08-24): the displayed markers stay uncorrected.** This is
what every number shipped so far has used, so flipping it would have changed
conclusions relative to past-comparable reports. No display code changed.

What DID change: the record no longer claims otherwise. `fdr_correction` now
reads "none applied to the reported markers ... Benjamini-Hochberg q-values are
computed and included per gene in the Excel export (p_value_fdr /
significance_fdr)", and the MIQE checklist prints the same. So the corrected view
remains available in the workbook for anyone who asks, and nothing asserts a
correction it does not apply.

Do NOT delete `_apply_fdr_correction` or the `p_value_fdr` /
`significance_fdr` columns as "unused" — they are the corrected view the
provenance promises. If a reviewer ever wants corrected markers on the charts,
the shape is one `graph_settings["sig_source"]` toggle read by `graph.py`, the
PPT writer, the Overview matrix and the results table.

### 3. Not investigated

The whole-platform health agent (dependency pinning, Cloud-vs-local drift,
per-rerun cost) was stopped before reporting. `requirements.txt` pinning and the
Cloud Python-version trap are covered by the 2026-04-01 and 2026-08-07 entries
above; the rest of that sweep is still open.

---

## 2026-08-25: "No `import X` anywhere" does NOT mean X is unused

**What went wrong:** The audit dropped the `matplotlib` and `seaborn` pins as
dead weight (7 packages, ~15MB per cold boot) after establishing that nothing in
the project writes `import matplotlib`. True — and irrelevant. The per-gene
results table styled its Fold_Change column with

```python
display_df.style.background_gradient(subset=["Fold_Change"], cmap="RdYlGn", ...)
```

and `Styler.background_gradient` does `import matplotlib` **inside the function,
at call time**. Removing the pin made every styled results table raise
`ImportError: background_gradient requires matplotlib`, taking the Analysis tab
down with it.

**Why it survived local testing:** the dev machine still had matplotlib 3.10.8
installed from the older requirements.txt. The full suite was green locally —
290 tests — while a clean install could not render the tab. The 2026-08-07 entry
warns about pins that resolve but do not install; this is the mirror image, a
package that installs *only* because it was already there.

**Rule — a dependency audit must consider transitive, call-time imports, not
just `import` statements.** Grep is necessary and not sufficient. Libraries that
lazily import an optional backend inside a method (pandas Styler, plotting
backends, `pandas.read_*` engines) are invisible to an import grep. Before
dropping a pin, run the suite in an environment that does NOT have it.

**Rule — verify dependency changes in a clean environment, never the dev one.**
```bash
uv venv /tmp/clean --python 3.12
VIRTUAL_ENV=/tmp/clean uv pip install -r requirements.txt pytest pytest-mock pytest-timeout
/tmp/clean/bin/python -m pytest tests/
```
This is what CI does, and it is the only run that means anything for a pin
change. A green local suite proves nothing when the question is "is this package
still needed".

**Fix:** `qpcr/utils.py::gradient_styles` reproduces the RdYlGn ramp directly
(ColorBrewer 11-class anchors, linear interpolation, pandas' own
`text_color_threshold=0.408` for the light/dark text flip) and is applied via
`Styler.apply`. Verified against matplotlib's own output across the table's
0..3 range: max per-channel difference 3/255, endpoints exact. Two guards in
`tests/test_utils.py` — one renders the styled table with `matplotlib` forced
unimportable, one fails if any source file calls `.background_gradient(` again.

**Caught by CI**, which existed for one day at that point (PR #17). It was
merged ahead of the audit PR precisely so the audit would be gated, and this is
the defect it caught.

---

## 2026-08-26: DECLINED — the significance stack stays in DATA units, not pixel `yshift`

This trailed the #21–#26 chart-geometry series as "still deferred, never
started" and kept resurfacing. It is now **considered and declined**, not
pending. Do not re-raise it without the trigger at the bottom.

**The proposal.** `qpcr/graph.py` (~line 670) places each stacked significance
glyph in data units, converting the pixel step through the plot height:

```python
_data_per_px = (y_max_auto / plot_h_px) if plot_h_px else 0.0
```

Moving those annotations to pixel `yshift` would remove the last consumer of the
hand-computed `plot_h_px` and make `xaxis automargin` usable. Both are real
gains. Neither is worth the cost.

**Why declined — 1. The blast radius is larger than the entire series that
preceded it.** The mechanism does not just position glyphs; when a stack would
overflow it **rewrites the axis range** (~lines 697-700):

```python
if (settings.get("y_min") is None and settings.get("y_max") is None
        and not is_log and _stack_top > y_max_auto):
    y_max_auto = _stack_top
```

That extension is why bars read as compressed today, and dropping it makes them
visibly taller. It fires on every chart carrying a significance marker — which
is most real reports. #26 moved 4 of 48 matrix cases and could name them; this
moves nearly every deck, to fix nothing anyone reported.

**Why declined — 2. The 48-case pixel gate is blind to it.** That matrix is
4 label modes x 4 `FIGURE_SIZE_PRESETS` x {Latin short, Latin long, Korean}.
There is **no significance dimension in it at all.** The one instrument that made
#21–#26 safe cannot observe the thing this change alters, so "the gate is green"
would be meaningless here. Extending the matrix is a prerequisite, not a
follow-up.

**Why declined — 3. It trades a guarantee for a hazard.** Data units *guarantee*
containment: the range is extended until the stack fits. Pixel-shifted
annotations are not considered by autorange, so a tall stack on a short figure
can render into the top margin or clip. That swaps a cosmetic, self-correcting
issue (slightly compressed bars) for a correctness defect in an exported deck —
the wrong direction for a report generator.

**Rule — "removes a coupling" is not on its own a reason to change rendering
code that ships images.** Internal tidiness does not justify reprinting every
deck. Weigh it against blast radius, and against whether the existing gate can
even see the change.

**What would reopen this:** a real chart where the axis extension visibly ruins
the output — a tall three-line stack squashing bars until the effect is
unreadable. If that appears, the first step is **adding a significance dimension
to the render matrix**, not editing `graph.py`.
