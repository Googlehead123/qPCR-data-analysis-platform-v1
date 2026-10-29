"""Drive the real app pipeline for one recipe and normalize every output.

The app is imported fresh under the mocked Streamlit from ``conftest.py``, its
clock frozen at ``FIXED_NOW``, and state is seeded the way the tabs seed it:
upload -> QC -> mapping -> analysis -> graphs -> export. Every decision a tab
makes is the app's own function (``parse_uploads``, ``ingest_frames``,
``suggest_group``, ``build_analysis_params``, ...), the same one the tab calls;
what is left here is the widgets: the values a recipe types, ticks and picks,
and the loops a tab runs over genes. ``tests/test_ui_golden_parity.py`` checks
that emulation against the real UI under AppTest.
"""

import codecs
import copy
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


class PipelineError(RuntimeError):
    """A pipeline step failed. ``detail`` holds the evidence (the message log).

    The text names the recipe by ``label`` and carries the detail only for a
    synthetic recipe: a local recipe's log quotes client file names.
    """

    def __init__(self, recipe: Recipe, what: str, detail: str = ""):
        self.detail = detail
        shown = f": {detail}" if detail and not recipe.is_local else ""
        super().__init__(f"{recipe.label}: {what}{shown}")


class _FrozenDatetime(datetime):
    @classmethod
    def now(cls, tz=None):
        return FIXED_NOW


class _Upload(io.BytesIO):
    """What ``st.file_uploader`` hands the app: a named BytesIO."""

    def __init__(self, name: str, data: bytes):
        super().__init__(data)
        self.name = name


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
            raise PipelineError(recipe, "test_data root not found for its CSVs")
        try:
            data = (root / entry).read_bytes()
        except OSError as exc:
            raise PipelineError(recipe, "a CSV could not be read", str(exc)) from None
        files.append(_Upload(Path(entry).name, data))
    return files


# ---- Upload tab ---------------------------------------------------------------


def _ingest(app, uploaded_files) -> None:
    # As at app L4532–4537, which compute the hash list the tab records.
    current_file_hashes = sorted(
        [hashlib.md5(f.getvalue()).hexdigest() for f in uploaded_files]
    )
    for f in uploaded_files:
        f.seek(0)
    app.ingest_frames(
        app.parse_uploads(uploaded_files), app.st.session_state,
        file_hashes=current_file_hashes,
    )


def _detect_hk(app) -> None:
    state = app.st.session_state
    hk_genes = app.detect_hk_genes(list(state.data["Target"].unique()))
    if hk_genes:
        # Widget first-render emulation of the "hk_select" selectbox (app
        # L4565–4574): the stored gene if it is still offered, else the first
        # one detected. test_ui_golden_parity pins it against the real widget.
        current = state.get("hk_gene")
        state.hk_gene = current if current in hk_genes else hk_genes[0]


# ---- QC tab -------------------------------------------------------------------


def _override_every_auto_qc_decision(state) -> None:
    """QC tab: untick every "Keep auto-selection", then "Apply overrides" (app L4911–4933)."""
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
        app.apply_auto_qc(float(recipe.auto_qc_threshold))  # the threshold box (app L4846–4848)
    for gene, sample, well in recipe.manual_exclusions:
        app.exclude_well(well, gene, sample)


# ---- Mapping tab --------------------------------------------------------------


def _map_samples(app, recipe: Recipe) -> None:
    state = app.st.session_state
    # The efficacy selectbox defaults to the suggestion, else the first entry (app L5409–5414).
    targets = list(state.data["Target"].unique())
    efficacy = recipe.efficacy or app.suggest_efficacy(targets) or next(iter(app.EFFICACY_CONFIG))
    state.selected_efficacy = efficacy
    config = app.EFFICACY_CONFIG[efficacy]
    group_types = app.group_types_for(config)

    # One row per sample in the upload's natural order (app L5480–5493), with
    # the recipe typed into the condition box (stripped; blank -> sample, app
    # L5553–5559) and the include checkbox, then the group re-derived.
    for sample in state.sample_order:
        cond = (recipe.conditions.get(sample, sample) or "").strip() or sample
        state.sample_mapping[sample] = {
            "condition": cond,
            "group": app.suggest_group(cond or sample, config, group_types),
            "include": bool(recipe.include.get(sample, True)),
        }

    if recipe.sample_order is not None:
        if sorted(recipe.sample_order) != sorted(state.sample_order):
            raise PipelineError(recipe, "sample_order must list every sample exactly once")
        # The drag list holds included samples; excluded ones follow (app L5611–5623).
        inc = {s: state.sample_mapping[s]["include"] for s in recipe.sample_order}
        state.sample_order = [s for s in recipe.sample_order if inc[s]] + [
            s for s in recipe.sample_order if not inc[s]
        ]

    state.excluded_samples = app.excluded_samples_from_mapping(state.sample_mapping)


def _run_analysis(app, recipe: Recipe, log) -> None:
    state = app.st.session_state
    if recipe.ttest_type not in app.TTEST_TYPES:
        raise PipelineError(recipe, f"ttest_type must be one of {app.TTEST_TYPES}")
    state.ttest_type = recipe.ttest_type
    ref, cmp1, cmp2, cmp3 = recipe.ref_sample, *recipe.cmp_samples
    # As the Mapping tab records them for auto-rerun (app L5954–5957).
    state["_last_ref_sample_key"] = ref
    state["_last_cmp_sample_key"] = cmp1
    state["_last_cmp_sample_key_2"] = cmp2
    state["_last_cmp_sample_key_3"] = cmp3
    if not app.AnalysisEngine.run_full_analysis(ref, cmp1, cmp2, cmp3):
        raise PipelineError(recipe, "run_full_analysis failed", _message_lines(log))


# ---- Export tab ---------------------------------------------------------------


def _export_params(app, concentration: str) -> dict:
    state = app.st.session_state
    state["exp_concentration"] = concentration  # typed into the "Sample concentration" box
    return app.build_analysis_params(state, now=app.datetime.now())


def _publication_images(app, image_export: dict) -> tuple[str, dict]:
    state = app.st.session_state
    img_format = image_export["fmt"]  # the Gene Images widgets (app L6758–6767)
    img_width, img_height = int(image_export["width"]), int(image_export["height"])
    fmt, _mime, scale = app.publication_image_format(img_format, img_width, img_height)
    rendered = {}
    for gene, fig in state.graphs.items():  # the render loop (app L6782–6790)
        fig_copy, adj_h = app.publication_figure(fig, img_width, img_height)
        rendered[gene] = app.ReportGenerator._fig_to_image(
            fig_copy, format=fmt, scale=scale, width=img_width, height=adj_h)
    return fmt, rendered


def _images_zip(app, fmt: str, rendered: dict) -> bytes:
    # The cached render's images and format (app L6802–6803), zipped as the tab does.
    return app.build_zip(app.image_zip_entries(rendered, fmt))


# ---- Analysis tab: per-gene results CSV ----------------------------------------


def _gene_csvs(app) -> dict[str, tuple[str, bytes]]:
    """``gene -> (file name, payload)`` exactly as the tab hands them to its download button."""
    state = app.st.session_state
    served = {}
    for gene, gene_df in state.processed_data.items():  # one expander per gene (app L6359)
        display_df = app.gene_results_frame(gene, gene_df)
        served[gene] = app.gene_results_csv(gene, display_df, now=app.datetime.now())
    return served


# ---- normalization -------------------------------------------------------------


def _round_csv(text: str) -> str:
    """Floats to 12 significant digits (-0 as 0), as every other part; all else verbatim."""
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    for row in csv.reader(io.StringIO(text)):
        writer.writerow(["%.12g" % (float(c) + 0.0) if _FLOAT_CELL.fullmatch(c) else c for c in row])
    return mask_text(out.getvalue())


def _served_csv(file_name: str, payload: bytes) -> str:
    """The served file: its name, whether it opens with a UTF-8 BOM, then its rows."""
    bom = payload.startswith(codecs.BOM_UTF8)
    body = payload[len(codecs.BOM_UTF8):] if bom else payload
    return f"file: {file_name}\nbom: {'utf-8' if bom else 'none'}\n" + _round_csv(body.decode("utf-8"))


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
        if state.get("data") is None:
            raise PipelineError(recipe, "no data parsed", _message_lines(log))
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
        if errors:
            raise PipelineError(recipe, f"{len(errors)} figure(s) failed", repr(errors))

        # Export tab, called exactly as at app L6696 and L6728.
        analysis_params = _export_params(app, recipe.concentration)
        _qc, _rep, _excl = app.build_export_extras(state)
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
        if not ppt:
            raise PipelineError(recipe, "PPT generation returned nothing", _message_lines(log))
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
            "gene_csv": {g: _served_csv(*served) for g, served in gene_csvs.items()},
            "messages": {"log": _message_lines(log)},
        }


def ingested_state(recipe: Recipe, keys) -> dict:
    """``keys`` of the session state right after upload ingestion, for seeding a real-Streamlit AppTest."""
    files = _uploads(recipe)
    with pytest.MonkeyPatch.context() as mp:
        app, log = _load_app(mp)
        state = app.st.session_state
        _ingest(app, files)
        if state.get("data") is None:
            raise PipelineError(recipe, "no data parsed", _message_lines(log))
        return {k: copy.deepcopy(state[k]) for k in keys}


def detect_defaults(csv_files: list[str]) -> dict:
    """What the app proposes for a plate set: HK gene, efficacy, sample order, t-test."""
    probe = Recipe(name="probe", csv_files=csv_files, ref_sample="", cmp_samples=[None] * 3)
    files = _uploads(probe)
    with pytest.MonkeyPatch.context() as mp:
        app, log = _load_app(mp)
        state = app.st.session_state
        _ingest(app, files)
        if state.get("data") is None:
            raise PipelineError(probe, "no data parsed", _message_lines(log))
        _detect_hk(app)
        return {
            "hk_gene": state.get("hk_gene"),
            "efficacy": app.suggest_efficacy(list(state.data["Target"].unique()))
            or next(iter(app.EFFICACY_CONFIG)),
            "samples": list(state.sample_order),
            "ttest_type": app.DEFAULT_TTEST_TYPE,
        }
