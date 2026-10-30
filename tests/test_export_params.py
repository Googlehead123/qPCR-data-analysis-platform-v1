"""The tab decisions the app and the golden harness share, pinned as functions.

``build_analysis_params`` feeds the Analysis_Parameters sheet and every slide's
metadata; ``suggest_group`` reaches ``processed_data["Group"]`` and the Excel
Summary grouping. Both used to be inline in tab bodies, so nothing could call
them without rendering the tab.
"""

import io
from importlib import import_module

import pandas as pd
import pytest

from tests.golden.driver import FIXED_NOW


def _app():
    return import_module("streamlit qpcr analysis v1")


def test_build_analysis_params_fields(mock_streamlit):
    app = _app()
    state = {
        "selected_efficacy": "진정",
        "hk_gene": "ACTIN",
        "analysis_ref_condition": "Non-treated",
        "analysis_cmp_condition": "Non-treated",
        "analysis_cmp_condition_2": "Active A",
        "excluded_wells": {("MMP1", "S2"): {"A6"}, ("MMP1", "S3"): {"B1", "B2"}},
        "excluded_samples": {"S6"},
        "processed_data": {"COL1A1": pd.DataFrame(), "MMP1": pd.DataFrame()},
        "ttest_type": "welch",
        # What the "Sample concentration" box holds; build_analysis_params strips it.
        "exp_concentration": "  0.1%  ",
    }
    params = app.build_analysis_params(state, now=FIXED_NOW)
    assert params == {
        "Date": "2026-01-01 09:00",
        "Efficacy_Type": "진정",
        "Housekeeping_Gene": "ACTIN",
        "Reference_Sample": "Non-treated",
        "Compare_To": "Non-treated",
        "Compare_To_2": "Active A",
        "Compare_To_3": None,
        "Excluded_Wells": 3,
        "Excluded_Samples": 1,
        "Genes_Analyzed": 2,
        "ttest_type": "welch",
        # Per-assay, from the catalog: 진정 is a 4 h assay, not the 24 h fallback.
        "treatment_time": app.EFFICACY_CONFIG["진정"]["treatment_time"],
        "concentration": "0.1%",
    }
    assert params["treatment_time"] == "4 h"
    # The writers read the dict in this order; concentration is appended last.
    assert list(params)[-1] == "concentration"


def test_build_analysis_params_fallbacks(mock_streamlit):
    """Unset comparisons, an unknown efficacy, a flat exclusion set and an
    untouched concentration box each fall back exactly as the tab did."""
    app = _app()
    state = {
        "selected_efficacy": None,
        "hk_gene": None,
        "excluded_wells": {"A1", "A2"},
        "excluded_samples": set(),
        "processed_data": {"COL1A1": pd.DataFrame()},
    }
    params = app.build_analysis_params(state, now=FIXED_NOW)
    assert params["Reference_Sample"] == "N/A"
    assert params["Compare_To"] == "N/A"
    assert params["Compare_To_2"] is None and params["Compare_To_3"] is None
    assert params["Excluded_Wells"] == 2
    assert params["ttest_type"] == app.DEFAULT_TTEST_TYPE
    assert params["treatment_time"] == "24 h"
    assert params["concentration"] == ""


def test_suggest_group_bare_number_is_treatment(mock_streamlit):
    """Under 탄력 the sample named "1" sits inside "TGFβ 10 ng/ml" and used to
    become a Positive Control (the regression suggest_group's substring-pass
    comment describes). Real exports name samples "1".."21"."""
    app = _app()
    config = app.EFFICACY_CONFIG["탄력"]
    group_types = app.group_types_for(config)
    assert group_types == ["Negative Control", "Positive Control", "Treatment"]
    assert app.suggest_group("1", config, group_types) == "Treatment"
    # Not vacuous: the real labels still match, exactly and by prefix.
    assert app.suggest_group("TGFβ 10 ng/ml", config, group_types) == "Positive Control"
    assert app.suggest_group("Non-treated 24h", config, group_types) == "Negative Control"
    assert app.suggest_group("  ", config, group_types) == "Treatment"


def test_suggest_group_role_outside_group_types_is_treatment(mock_streamlit):
    app = _app()
    config = app.EFFICACY_CONFIG["탄력"]
    assert app.suggest_group("Non-treated", config, ["Treatment"]) == "Treatment"


def test_excluded_samples_from_mapping(mock_streamlit):
    app = _app()
    mapping = {
        "S1": {"condition": "A", "include": True},
        "S2": {"condition": "B", "include": False},
        "S3": {"condition": "C"},  # no flag: included, as the tab defaults it
    }
    assert app.excluded_samples_from_mapping(mapping) == {"S2"}
    assert app.excluded_samples_from_mapping({}) == set()


def test_hk_and_efficacy_detection(mock_streamlit):
    app = _app()
    assert app.detect_hk_genes(["COL1A1", "gapdh", "MMP1", "ACTIN"]) == ["gapdh", "ACTIN"]
    assert app.detect_hk_genes([]) == []
    # The first catalog entry naming any detected gene, case-insensitively.
    assert app.suggest_efficacy(["actin", "col1a1"]) == "탄력"
    assert app.suggest_efficacy(["ACTIN"]) is None


def test_publication_image_format_steps_scale_down(mock_streamlit):
    app = _app()
    assert app.publication_image_format("PNG (3× pixels)", 1200, 800) == ("png", "image/png", 3)
    # 3000x2000 is 54 MP at 3x and 24 MP at 2x; 2000x1500 is 27 MP at 3x and
    # 12 MP at 2x. The scale steps down until the image is under ~16 MP.
    assert app.publication_image_format("PNG (3× pixels)", 3000, 2000)[2] == 1
    assert app.publication_image_format("PNG (3× pixels)", 2000, 1500)[2] == 2
    assert app.publication_image_format("SVG (Vector)", 3000, 2000) == ("svg", "image/svg+xml", 1)
    assert app.publication_image_format("PDF (Vector)", 400, 300) == ("pdf", "application/pdf", 1)


def _quiet(mock_st) -> None:
    """Forget the messages the app's import-time render emitted."""
    for kind in ("success", "warning", "error"):
        getattr(mock_st, kind).reset_mock()


class _Upload(io.BytesIO):
    def __init__(self, name: str):
        super().__init__(b"")
        self.name = name


def test_parse_uploads_reports_each_outcome(mock_streamlit, monkeypatch):
    app = _app()
    rows = pd.DataFrame({"Well": ["A1"], "Sample": ["S1"], "Target": ["ACTIN"], "CT": [18.0]})
    outcome = {"good.csv": rows, "empty.csv": rows.iloc[0:0], "broken.csv": None}
    monkeypatch.setattr(app.QPCRParser, "parse", staticmethod(lambda f: outcome[f.name]))
    _quiet(mock_streamlit)
    frames = app.parse_uploads([_Upload("good.csv"), _Upload("empty.csv"), _Upload("broken.csv")])
    assert len(frames) == 1 and frames[0]["Source_File"].tolist() == ["good.csv"]
    mock_streamlit.success.assert_called_once_with("good.csv: 1 wells parsed")
    mock_streamlit.warning.assert_called_once_with(
        "empty.csv: parsed successfully but contained no valid data rows."
    )
    mock_streamlit.error.assert_called_once_with(
        "broken.csv: could NOT be loaded — see the error above."
    )


def test_ingest_frames_with_nothing_usable_clears_and_advances(mock_streamlit):
    app = _app()
    state = app.st.session_state
    state["data"] = pd.DataFrame({"Sample": ["S1"]})
    state["processed_data"] = {"COL1A1": pd.DataFrame()}
    state["graphs"] = {"COL1A1": object()}
    state["_uploaded_file_hashes"] = ["old"]
    _quiet(mock_streamlit)
    app.ingest_frames([], state, file_hashes=["h"])
    assert state["data"] is None
    assert state["processed_data"] == {} and state["graphs"] == {}
    assert state["_uploaded_file_hashes"] == ["h"]
    mock_streamlit.error.assert_called_once()
    assert "previously loaded data has been cleared" in mock_streamlit.error.call_args.args[0]


def test_ingest_frames_refuses_any_state_but_the_session(mock_streamlit):
    """The reset and auto-QC it runs act on st.session_state; a stand-in would
    be left half-updated (data written to it, auto-QC run on the session)."""
    app = _app()
    with pytest.raises(TypeError, match="st.session_state"):
        app.ingest_frames([], {}, file_hashes=["h"])
    assert "_uploaded_file_hashes" not in app.st.session_state


def test_image_zip_entries_keep_colliding_names(mock_streamlit):
    app = _app()
    images = {"A/B": b"1", "A_B": b"2", "C": b"3"}
    entries = app.image_zip_entries(images, "png")
    stem = app._safe_image_stem("A/B")
    assert "/" not in stem
    assert entries == {f"{stem}.png": b"1", f"{stem}_2.png": b"2", "C.png": b"3"}
    assert app.image_zip_entries({}, "svg") == {}


@pytest.mark.parametrize("keys, expected_plain", [
    (["A", "A ", "A_2"], ["A", "A_2"]),
    (["A", "A ", "A  ", "A_2", "A_3"], ["A", "A_2", "A_3"]),
    (["A_2", "A", "A "], ["A_2", "A"]),
])
def test_image_zip_entries_names_are_unique_when_a_suffix_collides(
        mock_streamlit, keys, expected_plain):
    """A stem that already looks like ``X_2`` must not reuse a name allocated
    to a duplicate of ``X``: the later member overwrote the earlier one.
    Keys differ only by trailing spaces, which the stem strips."""
    app = _app()
    images = {k: str(i).encode() for i, k in enumerate(keys)}
    entries = app.image_zip_entries(images, "png")
    assert len(entries) == len(images)
    assert sorted(entries.values()) == sorted(images.values())
    # a stem's first occurrence keeps its plain name
    for stem in expected_plain:
        assert f"{stem}.png" in entries


def test_image_zip_members_survive_a_nul_in_the_gene_name(mock_streamlit):
    """Python's zip writer truncates a member name at NUL, so two stems that
    differed only after a NUL became duplicate members."""
    import io
    import zipfile

    app = _app()
    images = {"A\x00x": b"1", "A\x00y": b"2"}
    blob = app.build_zip(app.image_zip_entries(images, "png"))
    names = zipfile.ZipFile(io.BytesIO(blob)).namelist()
    assert len(names) == len(set(names)) == len(images)


@pytest.mark.parametrize("keys", [["A", "a"], ["A", "a", "A_2"], ["a", "A", "a_2", "A_2"]])
def test_image_zip_names_are_unique_case_insensitively(mock_streamlit, keys):
    """A.png and a.png overwrite each other when extracted on Windows/macOS."""
    app = _app()
    images = {k: str(i).encode() for i, k in enumerate(keys)}
    entries = app.image_zip_entries(images, "png")
    assert len(entries) == len(images)
    assert len({n.casefold() for n in entries}) == len(images)
    # the first occurrence keeps its name unchanged
    assert f"{keys[0]}.png" in entries


def test_image_zip_names_unique_case_insensitively_stay_unchanged(mock_streamlit):
    app = _app()
    entries = app.image_zip_entries({"Abc": b"1", "abd": b"2", "ABE": b"3"}, "png")
    assert list(entries) == ["Abc.png", "abd.png", "ABE.png"]


@pytest.mark.parametrize("name, expected", [
    ("CON", "_CON"), ("con", "_con"), ("Nul", "_Nul"), ("COM1", "_COM1"),
    ("lpt9", "_lpt9"), ("CON.txt", "_CON.txt"), ("AUX.x.y", "_AUX.x.y"),
    ("A\x01B", "A_B"), ("A\x1fB", "A_B"), ("A\tB", "A_B"),
    # unaffected: byte-identical
    ("CONSOLE", "CONSOLE"), ("COM10", "COM10"), ("COM0", "COM0"),
    ("ACTB", "ACTB"), ("IL-1B", "IL-1B"), ("col1a1", "col1a1"),
])
def test_safe_image_stem_windows_safe(mock_streamlit, name, expected):
    app = _app()
    assert app._safe_image_stem(name) == expected


def test_safe_image_stem_fallback_is_sanitised(mock_streamlit):
    """A blank display name falls back to the gene key, which must be cleaned
    too: it used to come back raw (with "/" or NUL)."""
    app = _app()
    app.st.session_state["gene_display_names"] = {"A/B\x00": "   ", "  ": " "}
    assert app._safe_image_stem("A/B\x00") == "A_B_"
    assert app._safe_image_stem("  ") == "image"


def test_image_zip_exact_names_for_one_collision_ordering(mock_streamlit):
    app = _app()
    entries = app.image_zip_entries({"A_2": b"0", "A": b"1", "A ": b"2"}, "png")
    assert list(entries) == ["A_2.png", "A.png", "A_3.png"]
