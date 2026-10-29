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
- a run that is then edited (a well excluded, the SD threshold moved, a
  condition renamed) and re-run gives the workbook a fresh run of that final
  state gives, so no cache or leftover state from the first run leaks into it.

AppTest runs in a SUBPROCESS (this file, run as a module): the repo-root
``conftest.py`` replaces ``streamlit`` with a MagicMock at import time.
AppTest cannot freeze the clock, so the analysis date and the provenance
"Generated" stamps are masked, and ``docProps`` dropped; nothing else is.
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
from tests.golden.normalize import normalize_xlsx

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
# Edits applied after the first run in "mutated", before the only run in "fresh_final".
MUTATION = {"exclude": ("COL1A1", "S1", "A2"), "sd_threshold": 0.25, "rename": ("S2", "Active A2")}
_EB_LABELS = {"livak_sd": "±SD", "sem": "±SEM", "ci95": "95% CI"}
_STAMP = re.compile(r"(Generated:\s*)\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}")


def _synthetic_a() -> recipes.Recipe:
    return next(r for r in recipes.synthetic_recipes() if r.name == "synthetic_a")


def _two_hk(data: pd.DataFrame) -> pd.DataFrame:
    """``data`` plus a GAPDH copy of its ACTIN wells: HK candidates ACTIN, GAPDH."""
    gapdh = data[data["Target"] == "ACTIN"].assign(Target="GAPDH")
    return pd.concat([data, gapdh], ignore_index=True)


# ---- the subprocess: real Streamlit ------------------------------------------------


def _drive(state_path: str, out_dir: str) -> None:
    from streamlit.testing.v1 import AppTest

    from tests.golden.normalize import normalize_figures

    seed = pickle.loads(Path(state_path).read_bytes())
    recipe, state = seed["recipe"], seed["state"]
    out = Path(out_dir)

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
        for sample, cond in recipe["conditions"].items():
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

    def prepared():
        at = seeded()
        apply_recipe_qc(at)
        map_samples(at)
        return at

    uploaded()
    hk_first_render()

    at = prepared()
    click(at, "Run Full Analysis")
    style_graphs(at)
    export(at, "parity")
    mutate(at)
    click(at, "Run Full Analysis")
    export(at, "mutated")

    at = prepared()
    mutate(at)
    click(at, "Run Full Analysis")
    style_graphs(at)
    export(at, "fresh_final")


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
    (tmp / "seed.pkl").write_bytes(pickle.dumps(seed))
    proc = subprocess.run(
        [sys.executable, "-m", "tests.test_ui_golden_parity", str(tmp / "seed.pkl"), str(tmp)],
        capture_output=True, text=True, timeout=1800, cwd=str(REPO),
    )
    if proc.returncode != 0:
        pytest.fail(
            "AppTest driver failed.\n"
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


def test_mutation_sequence_equals_fresh_run(ui_out):
    mutated = _ui_xlsx(ui_out, "mutated")
    fresh = _ui_xlsx(ui_out, "fresh_final")
    assert capture.compare({"xlsx": mutated}, {"xlsx": fresh}, allow=[]) == []
    # Not vacuous: the edits reached the workbook.
    assert capture.compare({"xlsx": mutated}, {"xlsx": _ui_xlsx(ui_out, "parity")}, allow=[])
    assert MUTATION["rename"][1] in mutated["sheet:Sample_Mapping"]


if __name__ == "__main__":
    _drive(*sys.argv[1:])
