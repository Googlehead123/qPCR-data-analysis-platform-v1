"""The golden harness against the real UI: same recipe, same workbook, same charts.

``tests/golden/driver.py`` runs the app's own functions under the repo's
Streamlit mock and emulates the widgets. These tests check that emulation
against the actual script under real Streamlit (``AppTest``), driven through
its widgets from the state right after upload ingestion:

- the Upload tab, fed ``synthetic_a``'s files, ingests them into that state,
  and its housekeeping-gene selectbox opens on the gene the driver assumes;
- ``synthetic_a`` through the UI gives the driver's workbook, and the driver's
  canonical figures (the real Graphs editor seeds per-gene settings the driver
  never renders, so the charts are compared too);
- what the tabs hand out gives the driver's outputs: the per-gene CSV, image
  and zip download buttons (recorded as the script calls them) and the deck;
- a run that is then edited (a well excluded, the SD threshold moved, a
  condition renamed) and re-run gives the workbook and charts that a fresh run
  of the final state gives. The fresh run is its own process (``st.cache_data``
  is process-global) and starts from the final QC state, renaming before its
  first analysis, so no cache entry it reads was filled from an earlier state.

AppTest runs in SUBPROCESSES (this file, run as a module): the repo-root
``conftest.py`` replaces ``streamlit`` with a MagicMock at import time. Figures
render through the driver's stub (``real_render=False``), patched where the
script imports it. AppTest cannot freeze the clock, so the analysis date, the
provenance "Generated" stamps, the slide "Date:" lines and the dates in
download file names are masked, and ``docProps`` dropped; nothing else is.
"""

import dataclasses
import json
import pickle
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from tests.golden import capture, driver, recipes
from tests.golden.normalize import _image_summary, normalize_pptx, normalize_xlsx

REPO = Path(__file__).resolve().parents[1]
APP = REPO / "streamlit qpcr analysis v1.py"

# Everything ingest_frames (with the reset and auto-QC it runs) leaves in the
# session, except _uploaded_file_hashes: with it set and the uploader empty, the
# Upload tab takes its "uploader emptied" branch and clears the data again.
INGESTED_KEYS = (
    "data", "sample_order", "processed_data", "graphs", "sample_mapping",
    "excluded_wells", "excluded_wells_history", "_auto_qc_owned", "excluded_samples",
    "hk_gene", "selected_efficacy", "gene_display_names", "qc_reviewed",
    "graph_settings", "auto_qc_audit", "qc_sd_threshold",
)
# What the Upload tab's own ingestion must leave exactly as the driver's does.
# The other keys INGESTED_KEYS lists are rewritten by the first render of the
# later tabs (the HK and efficacy selectboxes, the Mapping rows).
UPLOAD_KEYS = (
    "data", "sample_order", "excluded_wells", "_auto_qc_owned", "auto_qc_audit",
    "qc_sd_threshold", "_uploaded_file_hashes",
)
# The "hk_select" first render, per stored hk_gene, on a plate with two
# housekeeping candidates (synthetic_a has one, so it cannot tell a default
# index from any other). None is what ingestion leaves.
HK_STORED = {"none_stored": None, "second_stored": "GAPDH"}
# Edits made through the UI after the first run in "mutated"; "fresh_final" starts from their result.
MUTATION = {"exclude": ("COL1A1", "S1", "A2"), "sd_threshold": 0.25, "rename": ("S2", "Active A2")}
_EB_LABELS = {"livak_sd": "±SD", "sem": "±SEM", "ci95": "95% CI"}
_STAMP = re.compile(r"(Generated:\s*)\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}")
_SLIDE_DATE = re.compile(r"(<a:t>Date: )\d{4}-\d{2}-\d{2}(</a:t>)")
_FILE_DATE = re.compile(r"_\d{8}(\.\w+)$")


def _synthetic_a() -> recipes.Recipe:
    return next(r for r in recipes.synthetic_recipes() if r.name == "synthetic_a")


def _final_qc(app) -> None:
    """The QC decisions "mutated" makes through the UI, in its order, via the
    functions its handlers call: the recipe's exclusions (QC form), the
    mutation's exclusion (QC form), then the SD threshold box."""
    gene, sample, well = MUTATION["exclude"]
    for g, s, w in [*_synthetic_a().manual_exclusions, [gene, sample, well]]:
        app.exclude_well(w, g, s)
    app.apply_auto_qc(MUTATION["sd_threshold"])


def _two_hk(data: pd.DataFrame) -> pd.DataFrame:
    """``data`` plus a GAPDH copy of its ACTIN wells: HK candidates ACTIN, GAPDH."""
    gapdh = data[data["Target"] == "ACTIN"].assign(Target="GAPDH")
    return pd.concat([data, gapdh], ignore_index=True)


# ---- the subprocess: real Streamlit ------------------------------------------------


def _drive(state_path: str, out_dir: str, scenario: str) -> None:
    import qpcr.export_utils
    import streamlit
    from streamlit.testing.v1 import AppTest

    from tests.golden import stub_render
    from tests.golden.normalize import normalize_figures

    seed = pickle.loads(Path(state_path).read_bytes())
    recipe = seed["recipe"]
    # "fresh" starts from the final QC state; its rename is typed before its first analysis.
    state = seed["final_state"] if scenario == "fresh" else seed["state"]
    conditions = dict(recipe["conditions"])
    if scenario == "fresh":
        conditions[MUTATION["rename"][0]] = MUTATION["rename"][1]
    out = Path(out_dir)

    # The script re-imports export_figure_to_bytes from qpcr.export_utils on
    # every run, so this routes its renders through the driver's stub.
    qpcr.export_utils.export_figure_to_bytes = stub_render._export_figure_to_bytes
    # Every st.download_button call, as the script makes it.
    served = []
    real_download = streamlit.download_button

    def recording_download(*args, **kwargs):
        named = dict(zip(("label", "data", "file_name"), args), **kwargs)
        served.append({k: named.get(k) for k in ("key", "label", "data", "file_name", "mime")})
        return real_download(*args, **kwargs)

    streamlit.download_button = recording_download

    def run(target):
        """Run the script (``target`` is the AppTest or a widget just set on it); fail on any exception."""
        at = target.run()
        if at.exception:
            raise RuntimeError("; ".join(e.message for e in at.exception))
        return at

    def uploaded():
        """The real Upload tab on synthetic_a's files; AppTest has no file_uploader to drive."""
        import io

        import streamlit

        class Upload(io.BytesIO):
            def __init__(self, name, data):
                super().__init__(data)
                self.name = name

        real_uploader = streamlit.file_uploader
        streamlit.file_uploader = lambda *a, **k: [Upload(n, b) for n, b in seed["uploads"].items()]
        try:
            at = run(AppTest.from_file(str(APP), default_timeout=300))
        finally:
            streamlit.file_uploader = real_uploader
        after = {k: at.session_state[k] for k in UPLOAD_KEYS}
        after["success"] = [e.value for e in at.success]
        (out / "uploaded.pkl").write_bytes(pickle.dumps(after))

    def hk_first_render():
        picked = {}
        for case, stored in HK_STORED.items():
            at = AppTest.from_file(str(APP), default_timeout=300)
            at.session_state["data"] = seed["hk_data"]
            at.session_state["hk_gene"] = stored
            picked[case] = run(at).session_state["hk_gene"]
        (out / "hk.json").write_text(json.dumps(picked), encoding="utf-8")

    def seeded():
        at = AppTest.from_file(str(APP), default_timeout=300)
        for key, value in state.items():
            at.session_state[key] = value
        return run(at)

    def click(at, label):
        buttons = [b for b in at.button if b.label == label]
        assert len(buttons) == 1, f"{len(buttons)} buttons labelled {label!r}"
        return run(buttons[0].click())

    def exclude_well(at, gene, sample, well):
        """Untick the well in the QC tab's Triplicate Browser, then Apply Changes."""
        prefix = f"cb_{gene}_{well}_{sample}_"
        boxes = [c for c in at.checkbox if (c.key or "").startswith(prefix)]
        assert len(boxes) == 1, f"{len(boxes)} checkboxes for {prefix}"
        boxes[0].uncheck()
        run(at.button(key=f"FormSubmitter:well_form_{gene}-Apply Changes").click())

    def apply_recipe_qc(at):
        if recipe["auto_qc_threshold"] != 0.3:
            raise NotImplementedError("only the upload's own SD 0.3 is driven here")
        for gene, sample, well in recipe["manual_exclusions"]:
            exclude_well(at, gene, sample, well)

    def map_samples(at):
        efficacy = [s for s in at.selectbox if s.label == "🎯 Efficacy Test Type"]
        run(efficacy[0].select(recipe["efficacy"]))
        for sample, cond in conditions.items():
            at.text_input(key=f"cond_{sample}").input(cond)
        for sample, inc in recipe["include"].items():
            box = at.checkbox(key=f"include_{sample}")
            box.check() if inc else box.uncheck()
        run(at)
        # The drag list is a custom component AppTest cannot drag; it writes
        # sample_order, included samples first (app L5611-5623).
        inc = recipe["include"]
        order = recipe["sample_order"]
        at.session_state["sample_order"] = [s for s in order if inc.get(s, True)] + [
            s for s in order if not inc.get(s, True)
        ]
        run(at)
        choose_analysis(at)

    def choose_analysis(at):
        """The efficacy, reference, comparisons and t-test the recipe picks."""
        efficacy = [s for s in at.selectbox if s.label == "🎯 Efficacy Test Type"]
        run(efficacy[0].select(recipe["efficacy"]))
        ref, cmp1, cmp2, cmp3 = recipe["ref_sample"], *recipe["cmp_samples"]
        at.selectbox(key="ref_choice_ddct").select(ref)
        at.selectbox(key="cmp_choice_pval").select(cmp1)
        at.radio(key="ttest_type").set_value(recipe["ttest_type"])
        if cmp2:
            run(at.checkbox(key="use_second_pval").check())
            at.selectbox(key="cmp_choice_pval_2").select(cmp2)
        if cmp3:
            run(at.checkbox(key="use_third_pval").check())
            at.selectbox(key="cmp_choice_pval_3").select(cmp3)
        run(at)

    def chosen(at) -> dict:
        ss = at.session_state

        def get(key):
            return ss[key] if key in ss else None

        return {
            "efficacy": get("selected_efficacy"),
            "ref_choice_ddct": get("ref_choice_ddct"),
            "cmp_choice_pval": get("cmp_choice_pval"),
            "use_second_pval": get("use_second_pval"),
            "cmp_choice_pval_2": get("cmp_choice_pval_2"),
            "use_third_pval": get("use_third_pval"),
            "cmp_choice_pval_3": get("cmp_choice_pval_3"),
            "ttest_type": get("ttest_type"),
        }

    def style_graphs(at):
        # The editor shows one gene at a time; its pill button selects it.
        for gene, name in recipe["gene_display_names"].items():
            run(at.button(key=f"gene_btn_{gene}").click())
            run(at.text_input(key=f"gene_display_{gene}").input(name))
        for key, value in recipe["graph_settings"].items():
            if key == "error_bar_mode":
                run(at.selectbox(key="eb_mode_global").select(_EB_LABELS[value]))
            elif key.endswith("_y_max"):
                gene = key[: -len("_y_max")]
                run(at.button(key=f"gene_btn_{gene}").click())
                run(at.number_input(key=f"ymax_{gene}").set_value(value))
            else:
                raise NotImplementedError(f"graph setting {key!r} is not driven here")

    def mutate(at):
        exclude_well(at, *MUTATION["exclude"])
        run(at.number_input(key="qc_sd_threshold_input").set_value(MUTATION["sd_threshold"]))
        sample, name = MUTATION["rename"]
        run(at.text_input(key=f"cond_{sample}").input(name))

    def export(at, name):
        run(at.text_input(key="exp_concentration").input(recipe["concentration"]))
        run(at.button(key="gen_excel").click())
        (out / f"{name}.xlsx").write_bytes(at.session_state["_excel_export"]["data"])
        figures = normalize_figures(dict(at.session_state["graphs"]))
        (out / f"{name}.figures.json").write_text(json.dumps(figures), encoding="utf-8")

    def deck_and_downloads(at, name):
        """Generate the deck and the images; keep what the tabs hand out on that run."""
        run(at.button(key="gen_ppt").click())
        deck = at.session_state["_ppt_export"]["data"]
        (out / f"{name}.pptx").write_bytes(deck.getvalue() if hasattr(deck, "getvalue") else deck)
        served.clear()
        run(at.button(key="gen_images").click())
        (out / f"{name}.served.pkl").write_bytes(pickle.dumps(list(served)))

    def prepared():
        at = seeded()
        if scenario == "main":
            apply_recipe_qc(at)
        map_samples(at)
        return at

    if scenario == "main":
        uploaded()
        hk_first_render()
        at = prepared()
        click(at, "Run Full Analysis")
        style_graphs(at)
        export(at, "parity")
        deck_and_downloads(at, "parity")
        picked = chosen(at)
        mutate(at)
        # A QC-tab st.rerun() drops the state of widgets the run had not yet
        # rendered, so the Mapping selections can revert to their defaults
        # (a pre-existing app hazard, Ruling 16). Record what survived, then
        # re-pick them as the operator would.
        after = chosen(at)
        (out / "survival.json").write_text(
            json.dumps({k: after[k] == v for k, v in picked.items()}), encoding="utf-8"
        )
        choose_analysis(at)
        assert chosen(at) == picked, (picked, chosen(at))
        click(at, "Run Full Analysis")
        export(at, "mutated")
    elif scenario == "fresh":
        at = prepared()
        click(at, "Run Full Analysis")
        style_graphs(at)
        export(at, "fresh_final")
    else:
        raise ValueError(f"unknown scenario {scenario!r}")


# ---- the tests ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def ui_out(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("ui_parity")
    recipe = _synthetic_a()
    seed = {
        "recipe": dataclasses.asdict(recipe),
        "state": driver.ingested_state(recipe, INGESTED_KEYS),
        "uploads": {name: recipes.synthetic_csv(name) for name in recipe.csv_files},
    }
    seed["hk_data"] = _two_hk(seed["state"]["data"])
    seed["final_state"] = driver.ingested_state(recipe, INGESTED_KEYS, edit=_final_qc)
    (tmp / "seed.pkl").write_bytes(pickle.dumps(seed))
    # One process per scenario: "fresh" must not inherit "main"'s st.cache_data.
    for scenario in ("main", "fresh"):
        proc = subprocess.run(
            [sys.executable, "-m", "tests.test_ui_golden_parity", str(tmp / "seed.pkl"), str(tmp),
             scenario],
            capture_output=True, text=True, timeout=1800, cwd=str(REPO),
        )
        if proc.returncode != 0:
            pytest.fail(
                f"AppTest driver failed ({scenario}).\n"
                f"stdout tail:\n{proc.stdout[-4000:]}\n\nstderr tail:\n{proc.stderr[-4000:]}"
            )
    return tmp


@pytest.fixture(scope="module")
def driver_out():
    return driver.produce_outputs(_synthetic_a(), real_render=False)


def _masked(parts: dict[str, str]) -> dict[str, str]:
    """Drop docProps; mask the Analysis_Parameters date and the Provenance stamps."""
    out = {}
    for name, text in parts.items():
        if name.startswith("docProps/"):
            continue
        if name == "sheet:Analysis_Parameters":
            col = re.search(r"^([A-Z]+)1\tDate$", text, flags=re.M)
            assert col, "Analysis_Parameters has no Date column"
            text = re.sub(rf"^({col.group(1)}2\t).*$", r"\1<DATE>", text, flags=re.M)
        if name == "sheet:Provenance":
            assert _STAMP.search(text), "Provenance has no Generated stamp"
            text = _STAMP.sub(r"\1<TS>", text)
        out[name] = text
    return out


def _masked_deck(parts: dict[str, str]) -> dict[str, str]:
    """Drop docProps; mask each gene slide's "Date:" line."""
    return {n: _SLIDE_DATE.sub(r"\1<DATE>\2", t) for n, t in parts.items() if not n.startswith("docProps/")}


def _masked_file(text: str) -> str:
    """A served CSV part with the date in its file name masked."""
    head, _, rest = text.partition("\n")
    return _FILE_DATE.sub(r"_<DATE>\1", head) + "\n" + rest


def _cell(text: str, label: str) -> str:
    """The value beside ``label`` in a two-column sheet part."""
    row = re.search(rf"^A(\d+)\t{re.escape(label)}$", text, flags=re.M)
    assert row, f"no {label!r} row"
    value = re.search(rf"^B{row.group(1)}\t(.*)$", text, flags=re.M)
    return value.group(1) if value else ""


def _ui_xlsx(ui_out, name):
    return _masked(normalize_xlsx((ui_out / f"{name}.xlsx").read_bytes()))


def _ui_figures(ui_out, name):
    figures = json.loads((ui_out / f"{name}.figures.json").read_text(encoding="utf-8"))
    return {g: driver._pretty_json(t) for g, t in figures.items()}


def test_ui_upload_ingests_like_the_driver(ui_out, driver_out):
    ui = pickle.loads((ui_out / "uploaded.pkl").read_bytes())
    expected = driver.ingested_state(_synthetic_a(), UPLOAD_KEYS)
    pd.testing.assert_frame_equal(ui.pop("data"), expected.pop("data"))
    assert ui.pop("_uploaded_file_hashes") == expected.pop("_uploaded_file_hashes") != []
    shown = ui.pop("success")
    assert ui == expected
    # parse_uploads' per-file report, in file order, as the driver records it.
    parsed = [ln[len("success: "):] for ln in driver_out["messages"]["log"].splitlines()
              if ln.startswith("success: ") and ln.endswith(" wells parsed")]
    assert len(parsed) == len(_synthetic_a().csv_files)
    assert [m for m in shown if m in parsed] == parsed


def _driver_hk(data: pd.DataFrame, stored) -> str | None:
    """The driver's emulation of the "hk_select" first render on ``data``."""
    with pytest.MonkeyPatch.context() as mp:
        app, _log = driver._load_app(mp)
        app.st.session_state.data = data
        app.st.session_state.hk_gene = stored
        driver._detect_hk(app)
        return app.st.session_state.hk_gene


def test_ui_hk_selectbox_default_matches_the_driver(ui_out):
    ui = json.loads((ui_out / "hk.json").read_text(encoding="utf-8"))
    data = pickle.loads((ui_out / "seed.pkl").read_bytes())["hk_data"]
    assert ui == {case: _driver_hk(data, stored) for case, stored in HK_STORED.items()}
    # Not vacuous: the two cases pick different genes, so a changed default shows.
    assert ui == {"none_stored": "ACTIN", "second_stored": "GAPDH"}


def test_ui_workbook_equals_driver_workbook(ui_out, driver_out):
    actual = {"xlsx": _ui_xlsx(ui_out, "parity")}
    expected = {"xlsx": _masked(driver_out["xlsx"])}
    assert capture.compare(actual, expected, allow=[]) == []


def test_ui_figures_equal_driver_figures(ui_out, driver_out):
    actual = {"figures": _ui_figures(ui_out, "parity")}
    assert actual["figures"], "the Graphs tab built no figures"
    assert capture.compare(actual, {"figures": driver_out["figures"]}, allow=[]) == []


def _served(ui_out, name) -> dict[str, dict]:
    """Recorded download calls by key; the keyless report buttons by label."""
    calls = pickle.loads((ui_out / f"{name}.served.pkl").read_bytes())
    served = {c["key"] or c["label"]: c for c in calls}
    assert len(served) == len(calls), "two download buttons share a key or label"
    return served


def test_ui_served_csvs_equal_driver_csvs(ui_out, driver_out):
    """What the Analysis tab hands each per-gene CSV button: name and bytes."""
    calls = _served(ui_out, "parity")
    served = {
        k[len("csv_"):]: _masked_file(driver._served_csv(c["file_name"], c["data"]))
        for k, c in calls.items() if k.startswith("csv_")
    }
    expected = {g: _masked_file(t) for g, t in driver_out["gene_csv"].items()}
    assert served, "no CSV download buttons recorded"
    assert capture.compare({"gene_csv": served}, {"gene_csv": expected}, allow=[]) == []


def test_ui_served_images_equal_driver_images(ui_out, driver_out):
    """The images zip and per-gene image buttons the Export tab renders."""
    calls = _served(ui_out, "parity")
    zipped = {n: s + "\n" for n, s in driver._zip_members(calls["dl_images_zip"]["data"]).items()}
    images = {
        k[len("img_"):]: _image_summary(c["data"]) + "\n"
        for k, c in calls.items() if k.startswith("img_")
    }
    actual = {"images_zip": zipped, "images": images}
    expected = {"images_zip": driver_out["images_zip"], "images": driver_out["images"]}
    assert capture.compare(actual, expected, allow=[]) == []
    # The MIME type comes from publication_image_format, for the cached format.
    assert {c["mime"] for k, c in calls.items() if k.startswith("img_")} == {"image/png"}
    assert calls["dl_images_zip"]["mime"] == "application/zip"


def test_ui_deck_equals_driver_deck(ui_out, driver_out):
    """The deck "Generate PowerPoint" builds from the tab's own arguments."""
    actual = {"pptx": _masked_deck(normalize_pptx((ui_out / "parity.pptx").read_bytes()))}
    expected = {"pptx": _masked_deck(driver_out["pptx"])}
    assert any(_SLIDE_DATE.search(t) for t in driver_out["pptx"].values()), "no slide Date: line"
    assert capture.compare(actual, expected, allow=[]) == []
    # The report buttons serve exactly what Generate built.
    calls = _served(ui_out, "parity")
    assert calls["Download PowerPoint"]["data"] == (ui_out / "parity.pptx").read_bytes()
    assert calls["Download Excel report"]["data"] == (ui_out / "parity.xlsx").read_bytes()


def test_mutation_sequence_equals_fresh_run(ui_out):
    mutated = _ui_xlsx(ui_out, "mutated")
    fresh = _ui_xlsx(ui_out, "fresh_final")
    assert capture.compare({"xlsx": mutated}, {"xlsx": fresh}, allow=[]) == []
    # The chart path's memo, checked the same way.
    assert capture.compare(
        {"figures": _ui_figures(ui_out, "mutated")}, {"figures": _ui_figures(ui_out, "fresh_final")},
        allow=[],
    ) == []
    # Not vacuous: each edit reached the workbook.
    parity = _ui_xlsx(ui_out, "parity")
    gene, sample, well = MUTATION["exclude"]
    assert int(_cell(mutated["sheet:QC_Report"], "Excluded Wells")) == (
        int(_cell(parity["sheet:QC_Report"], "Excluded Wells")) + 1
    )
    assert f"{gene} / {sample} / well {well}" in mutated["sheet:Provenance"]
    assert f"{gene} / {sample} / well {well}" not in parity["sheet:Provenance"]
    assert f"SD <= {MUTATION['sd_threshold']:.2f}" in mutated["sheet:Provenance"]
    assert "SD <= 0.30" in parity["sheet:Provenance"]
    assert MUTATION["rename"][1] in mutated["sheet:Sample_Mapping"]
    assert MUTATION["rename"][1] not in parity["sheet:Sample_Mapping"]


if __name__ == "__main__":
    _drive(*sys.argv[1:])
