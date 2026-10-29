"""Drive the real app pipeline for one recipe and normalize every output.

The app is imported fresh under the mocked Streamlit from ``conftest.py``, its
clock frozen at ``FIXED_NOW``, and state is seeded the way the tabs seed it:
upload -> QC -> mapping -> analysis -> graphs -> export. Logic that is still
inline in a tab body is reached through TEMPORARY verbatim copies (only
``st.session_state`` -> ``state``); Task 3 replaces each with the app function
it extracts, and the goldens must not move when it does.
"""

import csv
import hashlib
import importlib
import io
import json
import re
import sys
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Callable
from unittest.mock import MagicMock

import pytest

from tests.golden import paths, stub_render
from tests.golden.normalize import (
    _image_summary,
    mask_text,
    normalize_figures,
    normalize_frames,
    normalize_pptx,
    normalize_xlsx,
)
from tests.golden.recipes import SYNTHETIC_CSVS, Recipe, synthetic_csv

APP_MODULE = stub_render.APP_MODULE
REPO = Path(__file__).resolve().parents[2]
FIXED_NOW = datetime(2026, 1, 1, 9, 0)
MESSAGE_KINDS = ("error", "warning", "info", "success")
_FLOAT_CELL = re.compile(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+[eE][-+]?\d+")


class _FrozenDatetime(datetime):
    @classmethod
    def now(cls, tz=None):
        return FIXED_NOW


class _Upload(io.BytesIO):
    """What ``st.file_uploader`` hands the app: a named BytesIO."""

    def __init__(self, name: str, data: bytes):
        super().__init__(data)
        self.name = name


def _first_render_selectbox(label, options, index=0, *args, **kwargs):
    """A selectbox with no stored widget state returns ``options[index]``."""
    opts = list(options)
    return opts[index] if opts and index is not None else None


def _load_app(mp):
    """Fresh app module under a fresh mock Streamlit; returns ``(app, message log)``."""
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    for name in list(sys.modules):
        if name.startswith(("streamlit", "qpcr")) or name == APP_MODULE:
            del sys.modules[name]
    import conftest

    sys.modules["streamlit"] = conftest._create_mock_streamlit()
    app = importlib.import_module(APP_MODULE)
    mp.setattr(app, "datetime", _FrozenDatetime)
    st = app.st
    st.selectbox = MagicMock(side_effect=_first_render_selectbox)
    # Recording starts after import: the import-time render has no data, so
    # what it says is tab chrome, not pipeline output.
    log: list[tuple[str, str]] = []
    for kind in MESSAGE_KINDS:
        setattr(st, kind, lambda body="", *a, _k=kind, **kw: log.append((_k, str(body))))
    return app, log


def _uploads(recipe: Recipe) -> list[_Upload]:
    files = []
    for entry in recipe.csv_files:
        if entry in SYNTHETIC_CSVS:
            files.append(_Upload(entry, synthetic_csv(entry)))
            continue
        root = paths.test_data_root()
        if root is None or not root.is_dir():
            raise FileNotFoundError(f"test_data root not found for local CSV {entry!r}")
        files.append(_Upload(Path(entry).name, (root / entry).read_bytes()))
    return files


# ---- Upload tab ---------------------------------------------------------------


def _ingest(app, uploaded_files) -> None:
    st, state, pd = app.st, app.st.session_state, app.pd
    QPCRParser, make_well_ids_unique = app.QPCRParser, app.make_well_ids_unique
    natural_sort_key, apply_auto_qc = app.natural_sort_key, app.apply_auto_qc
    _reset_session_for_new_data = app._reset_session_for_new_data
    _clear_all_gene_style_state = app._clear_all_gene_style_state
    # As at app L4160–4165, which compute the hash list the copy records.
    current_file_hashes = sorted(
        [hashlib.md5(f.getvalue()).hexdigest() for f in uploaded_files]
    )
    for f in uploaded_files:
        f.seek(0)

    # TEMPORARY verbatim copy of app L4170–4260; Task 3 extracts
    all_data = []
    for file in uploaded_files:
        parsed = QPCRParser.parse(file)
        # FIX-03: Reject empty DataFrames (parsed but no valid rows)
        if parsed is not None and not parsed.empty:
            parsed["Source_File"] = file.name
            all_data.append(parsed)
            st.success(f"{file.name}: {len(parsed)} wells parsed")
        elif parsed is not None and parsed.empty:
            st.warning(f"{file.name}: parsed successfully but contained no valid data rows.")
        else:
            # parse() returned None. It reports the reason itself, but
            # this branch did not exist, so the file's failure was
            # invisible in the upload list.
            st.error(f"{file.name}: could NOT be loaded — see the error above.")

    if not all_data:
        # Nothing usable arrived. This used to fall through leaving the
        # PREVIOUS file's data in session_state while the uploader showed
        # the new filename, so every metric, the preview, QC and the
        # analysis kept describing the old experiment. Clearing also lets
        # the hash record advance so the failed file is not re-parsed on
        # every rerun.
        if state.get("data") is not None:
            st.error(
                "None of the selected files could be loaded. The previously "
                "loaded data has been cleared so nothing below describes the "
                "wrong experiment."
            )
        state.data = None
        state.processed_data = {}
        state.graphs = {}
        _clear_all_gene_style_state()
        state._uploaded_file_hashes = current_file_hashes

    if all_data:
        _combined_raw = pd.concat(all_data, ignore_index=True)
        # Measure overlap BEFORE the Well IDs are made unique (below).
        _dup_check = ["Target", "Well", "Sample"]
        _dup_mask = _combined_raw.duplicated(subset=_dup_check, keep=False)
        _cross_file_dups = _within_file_dups = 0
        if _dup_mask.any():
            _dup_sources = _combined_raw[_dup_mask].groupby(_dup_check)[
                "Source_File"
            ].apply(lambda x: list(x.unique()))
            _cross_file_dups = int((_dup_sources.apply(len) > 1).sum())
            _within_file_dups = int((_dup_sources.apply(len) == 1).sum())
        # A repeated Target+Well+Sample made exclusions and auto-QC
        # ambiguous (excluding one copy removed all of them).
        state.data, _n_wells_renamed = make_well_ids_unique(
            _combined_raw
        )

        # FIX-02: Complete state reset on new upload to prevent stale data
        had_previous_analysis = bool(state.get("processed_data"))
        _reset_session_for_new_data()

        state._uploaded_file_hashes = current_file_hashes

        # Automatic replicate QC (best-2-of-3 at SD 0.3) so the user can
        # go straight to mapping; fully reviewable/reversible in the QC tab.
        apply_auto_qc(0.3)

        if had_previous_analysis:
            st.info("Previous analysis results cleared due to new data upload.")

        # FIX-05: Detect overlapping data across files AND within one
        # (counted above, before Well IDs were made unique).
        if _cross_file_dups:
            st.warning(
                f"⚠️ Found {_cross_file_dups} overlapping data point(s) "
                f"across files (same Target+Well+Sample in different "
                f"files). Each is kept as its own well (repeats get a "
                f"' (2)' suffix) and counted as an extra replicate. "
                f"Consider deduplicating your input files."
            )
        if _within_file_dups:
            st.warning(
                f"⚠️ Found {_within_file_dups} repeated Target+Well+Sample "
                f"row(s) within a single file. They are counted as extra "
                f"replicates (repeats get a ' (2)' suffix), which inflates "
                f"n and narrows the error bars."
            )

        unique_samples = sorted(
            state.data["Sample"].unique(), key=natural_sort_key
        )
        state.sample_order = unique_samples


def _detect_hk(app) -> None:
    st, state = app.st, app.st.session_state
    # TEMPORARY verbatim copy of app L4274–4305; Task 3 extracts
    # (through L4310, where the selectbox call that L4305 opens closes)
    KNOWN_HK_GENES = [
        "ACTIN",
        "B-ACTIN",
        "BACTIN",
        "GAPDH",
        "ACTB",
        "BETA-ACTIN",
        "BETAACTIN",
        "18S",
        "18S RRNA",
        "HPRT",
        "HPRT1",
        "B2M",
        "RPLP0",
        "TBP",
        "PPIA",
        "RPL13A",
        "YWHAZ",
        "SDHA",
        "HMBS",
        "UBC",
        "GUSB",
        "PGK1",
    ]
    all_genes = list(state.data["Target"].unique())
    hk_genes = [g for g in all_genes if g.upper() in KNOWN_HK_GENES]

    if hk_genes:
        default_idx = 0
        if state.get("hk_gene") in hk_genes:
            default_idx = hk_genes.index(state.hk_gene)
        state.hk_gene = st.selectbox(
            "🔬 Housekeeping Gene (auto-detected)",
            hk_genes,
            index=default_idx,
            key="hk_select",
        )


# ---- QC tab -------------------------------------------------------------------


def _override_every_auto_qc_decision(state) -> None:
    """QC tab: untick every "Keep auto-selection", then "Apply overrides" (app L4647–4669)."""
    audit = state.get("auto_qc_audit", []) or []
    for a in audit:
        dropped = a["dropped_wells"]
        if not dropped:
            continue
        cur = state.excluded_wells.get((a["Target"], a["Sample"]))
        if cur:
            for w in dropped:
                cur.discard(w)
            if not cur:
                del state.excluded_wells[(a["Target"], a["Sample"])]
    _owned = state.get("_auto_qc_owned", {}) or {}
    for a in audit:
        if a["dropped_wells"]:
            _owned.pop((a["Target"], a["Sample"]), None)
    state._auto_qc_owned = _owned


def _apply_qc(app, recipe: Recipe) -> None:
    state = app.st.session_state
    if recipe.auto_qc_threshold is None:
        _override_every_auto_qc_decision(state)
    elif abs(float(recipe.auto_qc_threshold) - float(state.get("qc_sd_threshold", 0.3))) > 1e-9:
        app.apply_auto_qc(float(recipe.auto_qc_threshold))  # the threshold box (app L4582–4584)
    for gene, sample, well in recipe.manual_exclusions:
        app.exclude_well(well, gene, sample)


# ---- Mapping tab --------------------------------------------------------------


def _suggest_efficacy(app):
    state, EFFICACY_CONFIG = app.st.session_state, app.EFFICACY_CONFIG
    # TEMPORARY verbatim copy of app L5147–5150; Task 3 extracts
    # (plus L5143–5146, which define the names the loop reads)
    detected_genes = set(state.data["Target"].unique())
    # FIX-07: Case-insensitive efficacy gene matching
    detected_genes_upper = {g.upper() for g in detected_genes}
    suggested = None
    for eff, cfg in EFFICACY_CONFIG.items():
        if any(g.upper() in detected_genes_upper for g in cfg["genes"]):
            suggested = eff
            break
    return suggested


def _group_suggester(config: dict) -> Callable[[str], str]:
    # TEMPORARY verbatim copy of app L5211–5256; Task 3 extracts
    # Group type options
    group_types = ["Negative Control", "Positive Control", "Treatment"]
    if "baseline" in config["controls"]:
        group_types.insert(0, "Baseline")

    # Smart default group: match a new sample's name against the efficacy
    # category's known control labels so e.g. "Non-treated" pre-selects
    # "Negative Control" instead of "Treatment". Always overridable below.
    _role_to_group = {
        "baseline": "Baseline",
        "negative": "Negative Control",
        "positive": "Positive Control",
    }

    def _suggest_group(sample_name):
        sl = str(sample_name).strip().lower()
        if not sl:
            return "Treatment"
        controls = config.get("controls", {})
        # Exact label match first (highest confidence), then substring.
        for role in ("baseline", "negative", "positive"):
            label = controls.get(role)
            if label and sl == str(label).strip().lower():
                grp = _role_to_group[role]
                return grp if grp in group_types else "Treatment"
        # Substring pass. This used to match in BOTH directions on any
        # length, so a short sample name inside a dosed control label won:
        # every real export here names samples "1".."21", and under 탄력 the
        # sample named "1" became a Positive Control because "1" is inside
        # "tgfβ 10 ng/ml" (under 진정, "4" inside "poly(i:c)+il-4"). The
        # derived group reaches processed_data["Group"], the results table and
        # the Excel Summary grouping, and analysis.py takes a pooled
        # condition's group from its first sample — so one mis-grouped sample
        # relabels a whole condition. Require at least 3 characters and a
        # prefix/containment relation, which keeps "Non-treated 24h" and
        # "Dexamethasone" matching while a bare number no longer does.
        for role in ("baseline", "negative", "positive"):
            label = controls.get(role)
            ll = str(label).strip().lower() if label else ""
            if not ll or len(sl) < 3:
                continue
            if ll.startswith(sl) or sl.startswith(ll) or ll in sl:
                grp = _role_to_group[role]
                return grp if grp in group_types else "Treatment"
        return "Treatment"

    return _suggest_group


def _map_samples(app, recipe: Recipe) -> None:
    state = app.st.session_state
    # The efficacy selectbox defaults to the suggestion, else the first entry (app L5152–5157).
    efficacy = recipe.efficacy or _suggest_efficacy(app) or next(iter(app.EFFICACY_CONFIG))
    state.selected_efficacy = efficacy
    _suggest_group = _group_suggester(app.EFFICACY_CONFIG[efficacy])

    # One row per sample in the upload's natural order (app L5267–5278), with
    # the recipe typed into the condition box (stripped; blank -> sample, app
    # L5338–5344) and the include checkbox, then the group re-derived.
    for sample in state.sample_order:
        cond = (recipe.conditions.get(sample, sample) or "").strip() or sample
        state.sample_mapping[sample] = {
            "condition": cond,
            "group": _suggest_group(cond or sample),
            "include": bool(recipe.include.get(sample, True)),
        }

    if recipe.sample_order is not None:
        if sorted(recipe.sample_order) != sorted(state.sample_order):
            raise ValueError(f"{recipe.name}: sample_order must list every sample exactly once")
        # The drag list holds included samples; excluded ones follow (app L5396–5408).
        inc = {s: state.sample_mapping[s]["include"] for s in recipe.sample_order}
        state.sample_order = [s for s in recipe.sample_order if inc[s]] + [
            s for s in recipe.sample_order if not inc[s]
        ]

    # TEMPORARY verbatim copy of app L5411–5417; Task 3 extracts
    # (through L5418, which closes the call)
    # Update excluded_samples from include flags
    state.excluded_samples = set(
        [
            s
            for s, v in state.sample_mapping.items()
            if not v.get("include", True)
        ]
    )


def _run_analysis(app, recipe: Recipe, log) -> None:
    state = app.st.session_state
    if recipe.ttest_type not in app.TTEST_TYPES:
        raise ValueError(f"{recipe.name}: ttest_type must be one of {app.TTEST_TYPES}")
    state.ttest_type = recipe.ttest_type
    ref, cmp1, cmp2, cmp3 = recipe.ref_sample, *recipe.cmp_samples
    # As the Mapping tab records them for auto-rerun (app L5743–5746).
    state["_last_ref_sample_key"] = ref
    state["_last_cmp_sample_key"] = cmp1
    state["_last_cmp_sample_key_2"] = cmp2
    state["_last_cmp_sample_key_3"] = cmp3
    ok = app.AnalysisEngine.run_full_analysis(ref, cmp1, cmp2, cmp3)
    assert ok, f"{recipe.name}: run_full_analysis failed: {log}"


# ---- Export tab ---------------------------------------------------------------


def _export_params(app, concentration: str) -> dict:
    state = app.st.session_state
    datetime, DEFAULT_TTEST_TYPE = app.datetime, app.DEFAULT_TTEST_TYPE
    EFFICACY_CONFIG = app.EFFICACY_CONFIG
    _conc = concentration  # the "Sample concentration" box (app L6439–6445)
    # TEMPORARY verbatim copy of app L6381–6398 and L6447; Task 3 extracts
    # (plus L6446, which strips the box's value)
    analysis_params = {
        "Date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "Efficacy_Type": state.selected_efficacy,
        "Housekeeping_Gene": state.hk_gene,
        "Reference_Sample": state.get("analysis_ref_condition", "N/A"),
        "Compare_To": state.get("analysis_cmp_condition", "N/A"),
        "Compare_To_2": state.get("analysis_cmp_condition_2"),
        "Compare_To_3": state.get("analysis_cmp_condition_3"),
        "Excluded_Wells": sum(len(ws) for ws in state.excluded_wells.values()) if isinstance(state.excluded_wells, dict) else len(state.excluded_wells),
        "Excluded_Samples": len(state.excluded_samples),
        "Genes_Analyzed": len(state.processed_data),
        "ttest_type": state.get("ttest_type", DEFAULT_TTEST_TYPE),
        # Treatment time comes from the efficacy catalog (per-assay) so the
        # PPT "Treatment time:" field is accurate (4 h / 6 h / 24 h / 48 h).
        "treatment_time": EFFICACY_CONFIG.get(
            state.selected_efficacy, {}
        ).get("treatment_time", "24 h"),
    }
    _conc = (_conc or "").strip()
    analysis_params["concentration"] = _conc
    return analysis_params


def _export_extras(app):
    state = app.st.session_state
    QualityControl = app.QualityControl
    get_replicate_stats, get_qc_summary_stats = app.get_replicate_stats, app.get_qc_summary_stats

    # TEMPORARY verbatim copy of app L6400–6421; Task 3 extracts
    def _build_export_extras():
        """Build qc_stats, replicate_stats, and excluded_wells for export."""
        qc_stats = None
        replicate_stats_df = None
        excl = state.get("excluded_wells", {})

        if state.get("data") is not None:
            # Dict-aware, matching the QC tab: flattening the per-(gene,
            # sample) exclusions to bare well IDs dropped that coordinate for
            # every gene sharing it, so the workbook's QC sheet disagreed
            # with both the screen and the DDCt the report is built on.
            qc_data = state.data.copy()
            _exp_mask = QualityControl.excluded_mask(qc_data, excl)
            if _exp_mask.any():
                qc_data = qc_data[~_exp_mask]

            replicate_stats_df = get_replicate_stats(qc_data)
            qc_summary = get_qc_summary_stats(state.data, excl or None)
            qc_stats = qc_summary if qc_summary else None

        return qc_stats, replicate_stats_df, excl

    return _build_export_extras()


def _publication_images(app, image_export: dict) -> tuple[str, dict]:
    state, go = app.st.session_state, app.go
    PLOTLY_FONT_FAMILY, ReportGenerator = app.PLOTLY_FONT_FAMILY, app.ReportGenerator
    img_format = image_export["fmt"]  # the Gene Images widgets (app L6604–6613)
    img_width, img_height = int(image_export["width"]), int(image_export["height"])
    # TEMPORARY verbatim copy of app L6616–6622; Task 3 extracts
    fmt = "png" if "PNG" in img_format else "svg" if "SVG" in img_format else "pdf"
    mime = {"png": "image/png", "svg": "image/svg+xml", "pdf": "application/pdf"}[fmt]
    scale = 3 if fmt == "png" else 1
    # A 3000x2000 canvas at 3x is 54 MP PER GENE — enough to exhaust a
    # Streamlit Cloud worker. Step the scale down to stay under ~16 MP.
    while scale > 1 and img_width * img_height * scale * scale > 16_000_000:
        scale -= 1

    rendered = {}
    for gene, fig in state.graphs.items():  # the render loop (app L6634–6647)
        # TEMPORARY verbatim copy of app L6637–6644; Task 3 extracts
        fig_copy = go.Figure(fig)
        _orig_m = fig.layout.margin
        _pub_b = max(180, _orig_m.b if _orig_m and _orig_m.b else 180)
        _adj_h = img_height + max(0, _pub_b - 180)
        fig_copy.update_layout(width=img_width, height=_adj_h, margin=dict(b=_pub_b),
            font=dict(size=14, family=PLOTLY_FONT_FAMILY, color="black"))
        rendered[gene] = ReportGenerator._fig_to_image(
            fig_copy, format=fmt, scale=scale, width=img_width, height=_adj_h)
    return fmt, rendered


def _images_zip(app, fmt: str, rendered: dict) -> bytes:
    _safe_image_stem, build_zip = app._safe_image_stem, app.build_zip
    cached, _cfmt = {"fmt": fmt, "images": rendered}, fmt  # the cached render (app L6659–6660)
    # TEMPORARY verbatim copy of app L6673–6680; Task 3 extracts
    _zip_entries = {}
    _used_stems: dict = {}
    for _g, _b in cached["images"].items():
        _stem = _safe_image_stem(_g)
        _n = _used_stems.get(_stem, 0) + 1
        _used_stems[_stem] = _n
        _name = _stem if _n == 1 else f"{_stem}_{_n}"
        _zip_entries[f"{_name}.{_cfmt}"] = _b
    return build_zip(_zip_entries)


# ---- Analysis tab: per-gene results CSV ----------------------------------------


def _gene_csvs(app) -> dict[str, str]:
    st, state, datetime = app.st, app.st.session_state, app.datetime
    EFFICACY_CONFIG, expected_direction_for = app.EFFICACY_CONFIG, app.expected_direction_for
    gradient_styles = app.gradient_styles
    csvs = {}
    # TEMPORARY verbatim copy of app L6146–6207; Task 3 extracts
    # (through L6211, which closes the download call; the csvs[...] line is the driver's)
    st.subheader("Gene-by-Gene Results")

    for gene, gene_df in state.processed_data.items():
        with st.expander(f"📍 {gene}", expanded=False):
            # Show expected direction if available
            efficacy_config = EFFICACY_CONFIG.get(
                state.selected_efficacy, {}
            )
            if "expected_direction" in efficacy_config:
                direction = expected_direction_for(gene, efficacy_config["expected_direction"])
                if direction:
                    st.caption(
                        f"Expected: {'↑ Increase' if direction == 'up' else '↓ Decrease'}"
                    )

            # Display columns
            display_cols = [
                "Condition",
                "Group",
                "Fold_Change",
                "p_value",
                "significance",
                "n_replicates",
                "Target_Ct_Mean",
                "Target_Ct_SD",
                "Target_Ct_CV%",
                "HK_Ct_Mean",
                "HK_Ct_SD",
                "HK_Ct_CV%",
                "Delta_Ct",
                "SEM",
                "SD",
            ]

            # Filter to existing columns
            display_df = gene_df[[c for c in display_cols if c in gene_df.columns]]

            # Style the dataframe. Not Styler.background_gradient: it
            # imports matplotlib at call time, which this project does not
            # ship (see qpcr/utils.gradient_styles).
            styled = display_df.style.apply(
                gradient_styles, subset=["Fold_Change"], vmin=0, vmax=3
            ).format(
                {
                    "Fold_Change": "{:.3f}",
                    "p_value": "{:.4f}",
                    "Target_Ct_Mean": "{:.2f}",
                    "Target_Ct_SD": "{:.3f}",
                    "Target_Ct_CV%": "{:.1f}%",
                    "HK_Ct_Mean": "{:.2f}",
                    "HK_Ct_SD": "{:.3f}",
                    "HK_Ct_CV%": "{:.1f}%",
                    "Delta_Ct": "{:.2f}",
                    "SEM": "{:.3f}",
                    "SD": "{:.3f}",
                },
                na_rep="—",
            )

            st.dataframe(styled, width='stretch')
            st.download_button(
                "⬇️ Download CSV",
                data=display_df.to_csv(index=False).encode("utf-8-sig"),
                file_name=f"{gene}_results_{datetime.now().strftime('%Y%m%d')}.csv",
                mime="text/csv", key=f"csv_{gene}", width='stretch',
            )
            csvs[gene] = display_df.to_csv(index=False)
    return csvs


# ---- normalization -------------------------------------------------------------


def _round_csv(text: str) -> str:
    """Floats to 12 significant digits, as every other part; all else verbatim."""
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    for row in csv.reader(io.StringIO(text)):
        writer.writerow(["%.12g" % float(c) if _FLOAT_CELL.fullmatch(c) else c for c in row])
    return mask_text(out.getvalue())


def _pretty_json(text: str) -> str:
    """One JSON value per line, so a diff names the attribute that moved."""
    return json.dumps(json.loads(text), indent=1, sort_keys=True, ensure_ascii=False) + "\n"


def _zip_members(data: bytes) -> dict[str, str]:
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        return {name: _image_summary(z.read(name)) for name in z.namelist()}


def _message_lines(log) -> str:
    """One ``kind: text`` line per message, in call order; newlines escaped."""
    lines = []
    for kind, text in log:
        flat = mask_text(text).replace("\r", "\\r").replace("\n", "\\n")
        lines.append(f"{kind}: {flat}\n")
    return "".join(lines)


# ---- entry points --------------------------------------------------------------


def produce_outputs(
    recipe: Recipe, *, real_render: bool, mutate: Callable[[object], None] | None = None
) -> dict[str, dict[str, str]]:
    """Run ``recipe`` through the app and return every output, normalized.

    ``mutate(app)`` runs right after ``run_full_analysis`` and before figures
    and exports; the harness self-test uses it to plant a change.
    """
    files = _uploads(recipe)
    with pytest.MonkeyPatch.context() as mp:
        app, log = _load_app(mp)
        if not real_render:
            stub_render.install_stub(mp)
        state = app.st.session_state

        _ingest(app, files)
        assert state.get("data") is not None, f"{recipe.name}: no data parsed: {log}"
        _detect_hk(app)
        if recipe.hk_gene is not None:
            state.hk_gene = recipe.hk_gene
        _apply_qc(app, recipe)
        _map_samples(app, recipe)
        _run_analysis(app, recipe, log)
        if mutate is not None:
            mutate(app)

        # Graphs tab: the operator's renames and chart settings, then every figure.
        state.gene_display_names.update(recipe.gene_display_names)
        state.graph_settings.update(recipe.graph_settings)
        genes = list(state.processed_data.keys())
        _figs, errors = app.figures_for_genes(genes, app.EFFICACY_CONFIG[state.selected_efficacy])
        assert not errors, f"{recipe.name}: figures failed: {errors}"

        # Export tab, called exactly as at app L6542 and L6574.
        analysis_params = _export_params(app, recipe.concentration)
        _qc, _rep, _excl = _export_extras(app)
        xlsx = app.export_to_excel(
            state.data, state.processed_data,
            analysis_params, state.sample_mapping,
            qc_stats=_qc, replicate_stats=_rep, excluded_wells=_excl,
            gene_display_names=state.get("gene_display_names", {}),
            graph_settings=state.get("graph_settings", {}),
            provenance=app._current_provenance(),
        )
        ppt = app.PPTGenerator.generate_presentation(
            state.graphs, state.processed_data,
            analysis_params, graph_settings=state.get("graph_settings"),
            gene_display_names=state.get("gene_display_names", {}),
        )
        assert ppt, f"{recipe.name}: PPT generation returned nothing: {log}"
        ppt = ppt.getvalue() if hasattr(ppt, "getvalue") else ppt
        fmt, rendered = _publication_images(app, recipe.image_export)
        zipped = _images_zip(app, fmt, rendered)
        gene_csvs = _gene_csvs(app)

        return {
            "processed": normalize_frames(state.processed_data),
            "figures": {g: _pretty_json(t) for g, t in normalize_figures(state.graphs).items()},
            "xlsx": normalize_xlsx(xlsx),
            "pptx": normalize_pptx(ppt),
            "images": {g: _image_summary(b) + "\n" for g, b in rendered.items()},
            "images_zip": {n: s + "\n" for n, s in _zip_members(zipped).items()},
            "gene_csv": {g: _round_csv(t) for g, t in gene_csvs.items()},
            "messages": {"log": _message_lines(log)},
        }


def detect_defaults(csv_files: list[str]) -> dict:
    """What the app proposes for a plate set: HK gene, efficacy, sample order, t-test."""
    probe = Recipe(name="probe", csv_files=csv_files, ref_sample="", cmp_samples=[None] * 3)
    files = _uploads(probe)
    with pytest.MonkeyPatch.context() as mp:
        app, log = _load_app(mp)
        state = app.st.session_state
        _ingest(app, files)
        if state.get("data") is None:
            raise ValueError(f"no data parsed from {len(csv_files)} file(s)")
        _detect_hk(app)
        return {
            "hk_gene": state.get("hk_gene"),
            "efficacy": _suggest_efficacy(app) or next(iter(app.EFFICACY_CONFIG)),
            "samples": list(state.sample_order),
            "ttest_type": app.DEFAULT_TTEST_TYPE,
        }
