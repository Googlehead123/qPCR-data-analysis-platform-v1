"""Tests for the golden-output normalizer, lock check and test-data root."""

import importlib.metadata
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import openpyxl
from openpyxl.chart import BarChart, Reference
import pandas as pd
import plotly.graph_objects as go
from PIL import Image

from tests.golden import env, normalize, paths


REPO = Path(__file__).resolve().parents[1]


def test_lock_mismatches_reports_difference(monkeypatch):
    pins = (REPO / "requirements.txt").read_text().splitlines()
    locked = next(ln.split("==", 1)[1].split()[0] for ln in pins if ln.startswith("numpy=="))
    real_version = importlib.metadata.version
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: "0.0.1" if name == "numpy" else real_version(name),
    )
    assert f"numpy: installed 0.0.1 != locked {locked}" in env.lock_mismatches(REPO / "requirements.txt")


def test_lock_mismatches_reports_missing_and_skips_foreign_markers(tmp_path):
    req = tmp_path / "requirements.txt"
    req.write_text(
        "# comment\n"
        "    # via something\n"
        "not-a-real-pkg-xyz==1.0\n"
        "win-only-pkg==1.0 ; sys_platform == 'win32' and python_version < '2'\n"
    )
    assert env.lock_mismatches(req) == [
        "not-a-real-pkg-xyz: installed <missing> != locked 1.0"
    ]


def test_test_data_root_prefers_env(monkeypatch, tmp_path):
    monkeypatch.setenv("QPCR_TEST_DATA", str(tmp_path))
    assert paths.test_data_root() == tmp_path
    assert paths.recipes_dir() == tmp_path / "golden" / "recipes"
    assert paths.local_out("abc123") == tmp_path / "golden" / "out" / "abc123"


def test_test_data_root_none_outside_repo(monkeypatch, tmp_path):
    # tmp_path may itself sit inside some repository; stop git's search at it.
    fresh = tmp_path / "fresh"
    fresh.mkdir()
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))
    monkeypatch.delenv("QPCR_TEST_DATA", raising=False)
    monkeypatch.chdir(fresh)
    assert paths.test_data_root() is None
    assert paths.recipes_dir() is None
    assert paths.local_out("abc123") is None


def _xlsx_bytes() -> bytes:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Data"
    ws["A1"] = "2026-10-29 14:32"
    ws["B2"] = 0.1234567890123456
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def test_xlsx_cells_one_line_each_and_dates_preserved():
    parts = normalize.normalize_xlsx(_xlsx_bytes())
    lines = parts["sheet:Data"].splitlines()
    assert "A1\t2026-10-29 14:32" in lines
    assert "B2\t0.123456789012" in lines
    core = parts["docProps/core.xml"]
    assert "<dcterms:created" in core
    assert core.count("<TS>") == 2
    assert "<TS></dcterms:created>" in core
    assert "<TS></dcterms:modified>" in core


def _styled_xlsx(*, as_text=False, fmt=None, merge=False, width=None) -> bytes:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Data"
    ws["A1"] = "1.5" if as_text else 1.5
    if fmt:
        ws["A1"].number_format = fmt
    if merge:
        ws.merge_cells("B1:C1")
    if width:
        ws.column_dimensions["A"].width = width
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def test_xlsx_layout_changes_show_without_cell_values():
    base = normalize.normalize_xlsx(_styled_xlsx())
    sheet_xml = base["xl/worksheets/sheet1.xml"]
    assert "<v>" not in sheet_xml  # values live in the sheet:<name> part only
    for variant in (
        _styled_xlsx(as_text=True),  # same text, now a string cell
        _styled_xlsx(fmt="0.000"),
        _styled_xlsx(merge=True),
        _styled_xlsx(width=30),
    ):
        changed = normalize.normalize_xlsx(variant)
        assert changed["sheet:Data"] == base["sheet:Data"]
        assert (changed["xl/worksheets/sheet1.xml"], changed["xl/styles.xml"]) != (
            sheet_xml, base["xl/styles.xml"]
        )


def test_pptx_keeps_presentation_and_relationship_parts():
    def deck(order: str, target: str) -> bytes:
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr("ppt/presentation.xml", f"<p:presentation xmlns:p='urn:x'>{order}</p:presentation>")
            z.writestr("ppt/slides/_rels/slide1.xml.rels", f"<Relationships><R Target='{target}'/></Relationships>")
            z.writestr("ppt/slides/slide1.xml", "<p:sld xmlns:p='urn:x'/>")
            z.writestr("docProps/thumbnail.jpeg", b"\xff\xd8 not decodable")
        return buf.getvalue()

    a = normalize.normalize_pptx(deck("<p:s id='1'/><p:s id='2'/>", "../media/image1.png"))
    b = normalize.normalize_pptx(deck("<p:s id='2'/><p:s id='1'/>", "../media/image2.png"))
    assert a["ppt/presentation.xml"] != b["ppt/presentation.xml"]
    assert a["ppt/slides/_rels/slide1.xml.rels"] != b["ppt/slides/_rels/slide1.xml.rels"]
    assert a["docProps/thumbnail.jpeg"].startswith("<binary ")


def test_version_token_masked():
    assert normalize.mask_text("Software: qPCR Analysis Suite v3.1") == "Software: <APP_VERSION>"
    assert (
        normalize.mask_text("qPCR Analysis Suite 2026.10 (commit 1a2b3c4-dirty)")
        == "<APP_VERSION>"
    )


def _png(pixel: tuple) -> bytes:
    img = Image.new("RGBA", (10, 20), (255, 255, 255, 255))
    img.putpixel((3, 4), pixel)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def _pptx(png: bytes) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("ppt/slides/slide1.xml", "<p:sld xmlns:p='urn:x'><p:a>hi</p:a></p:sld>")
        z.writestr("ppt/media/image1.png", png)
        z.writestr("[Content_Types].xml", "<Types/>")
    return buf.getvalue()


def test_pptx_media_hash_changes_with_pixels():
    a = normalize.normalize_pptx(_pptx(_png((0, 0, 0, 255))))
    b = normalize.normalize_pptx(_pptx(_png((1, 0, 0, 255))))
    ka = a["ppt/media/image1.png"]
    kb = b["ppt/media/image1.png"]
    assert ka.startswith("<image 10x20 sha256:")
    assert kb.startswith("<image 10x20 sha256:")
    assert ka != kb
    assert a["ppt/slides/slide1.xml"] == b["ppt/slides/slide1.xml"]
    assert "\n" in a["ppt/slides/slide1.xml"]  # pretty-printed


def test_figure_json_is_canonical():
    f1 = go.Figure(go.Bar(x=["a", "b"], y=[1.0, 2.0], name="s"))
    f2 = go.Figure(go.Bar(name="s", y=[1.0, 2.0], x=["a", "b"]))
    f3 = go.Figure(go.Bar(x=["a", "b"], y=[1.0, 2.5], name="s"))
    out1 = normalize.normalize_figures({"f": f1})
    assert out1 == normalize.normalize_figures({"f": f2})
    assert out1 != normalize.normalize_figures({"f": f3})


def test_figure_floats_rounded_to_12_digits():
    f1 = go.Figure(go.Bar(x=["a"], y=[0.1234567890123456]))
    f2 = go.Figure(go.Bar(x=["a"], y=[0.1234567890123459]))
    assert normalize.normalize_figures({"f": f1}) == normalize.normalize_figures({"f": f2})


def test_frames_csv_rounded():
    df = pd.DataFrame({"a": [0.1234567890123456], "b": ["x"]})
    assert normalize.normalize_frames({"t": df})["t"] == ",a,b\n0,0.123456789012,x\n"


def _np_fig(y: list[float]) -> go.Figure:
    return go.Figure(go.Bar(x=["a", "b"], y=np.array(y), error_y={"array": np.array(y)}))


def test_numpy_backed_figure_floats_rounded_and_decoded():
    base = normalize.normalize_figures({"f": _np_fig([0.1234567890123456, 2.0])})["f"]
    ulp = normalize.normalize_figures({"f": _np_fig([0.1234567890123459, 2.0])})["f"]
    real = normalize.normalize_figures({"f": _np_fig([0.1234567891, 2.0])})["f"]
    assert "bdata" not in base
    assert base == ulp
    assert base != real


def _sheet_lines(values: dict) -> list[str]:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Data"
    for coord, value in values.items():
        ws[coord] = value
    buf = io.BytesIO()
    wb.save(buf)
    return normalize.normalize_xlsx(buf.getvalue())["sheet:Data"].splitlines()


def _raw_number_lines(stored: list[str]) -> list[str]:
    """Sheet lines for cells whose ``<v>`` text is exactly ``stored``.

    openpyxl writes 1.0 as "1", so the stored text is patched into the XML: a
    writer that emits "1.0" loads back as a float, "1" as an int.
    """
    wb = openpyxl.Workbook()
    wb.active.title = "Data"
    for i in range(len(stored)):
        wb.active[f"A{i + 1}"] = 1000.5 + i
    buf = io.BytesIO()
    wb.save(buf)
    out = io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(buf.getvalue())) as src, zipfile.ZipFile(out, "w") as dst:
        for item in src.namelist():
            data = src.read(item)
            if item == "xl/worksheets/sheet1.xml":
                text = data.decode()
                for i, value in enumerate(stored):
                    text = text.replace(f"<v>{1000.5 + i}</v>", f"<v>{value}</v>")
                data = text.encode()
            dst.writestr(item, data)
    return normalize.normalize_xlsx(out.getvalue())["sheet:Data"].splitlines()


def test_integral_float_cell_differs_from_int_cell():
    lines = _raw_number_lines(["1", "1.0", "2.5", "1E+20"])
    assert lines == ["A1\t1", "A2\t1.0", "A3\t2.5", "A4\t1e+20"]


def test_negative_zero_normalizes_to_zero():
    assert _raw_number_lines(["-0.0", "-0"]) == ["A1\t0.0", "A2\t0"]
    assert normalize.normalize_frames({"t": pd.DataFrame({"a": [-0.0]})})["t"] == ",a\n0,0\n"
    neg = normalize.normalize_figures({"f": go.Figure(go.Bar(y=np.array([-0.0, 1.0])))})
    pos = normalize.normalize_figures({"f": go.Figure(go.Bar(y=np.array([0.0, 1.0])))})
    assert neg == pos


def test_datetime64_figure_values_render_as_iso():
    x = np.array(["2026-01-01T09:00:00"], dtype="datetime64[ns]")
    canon = json.loads(normalize.normalize_figures({"f": go.Figure(go.Scatter(x=x, y=[1.0]))})["f"])
    assert canon["data"][0]["x"] == ["2026-01-01T09:00:00.000000000"]


def test_version_mask_reaches_cells_frames_and_figures():
    stamp = "qPCR Analysis Suite v3.1 (commit abc1234)"
    assert _sheet_lines({"A1": f"Software: {stamp}"}) == ["A1\tSoftware: <APP_VERSION>"]
    frame = normalize.normalize_frames({"t": pd.DataFrame({"v": [stamp]})})["t"]
    assert "<APP_VERSION>" in frame and "v3.1" not in frame
    fig = normalize.normalize_figures({"f": go.Figure(layout={"title": {"text": stamp}})})["f"]
    assert "<APP_VERSION>" in fig and "v3.1" not in fig


def test_xlsx_chart_part_is_kept_pretty_printed():
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Data"
    for row in ([1], [2], [3]):
        ws.append(row)
    chart = BarChart()
    chart.add_data(Reference(ws, min_col=1, min_row=1, max_row=3))
    ws.add_chart(chart, "C1")
    buf = io.BytesIO()
    wb.save(buf)
    chart_part = normalize.normalize_xlsx(buf.getvalue())["chart:xl/charts/chart1.xml"]
    assert chart_part.startswith("<chartSpace") and "\n  <chart>" in chart_part


def _identity_deck() -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("docProps/core.xml", (
            "<cp:coreProperties xmlns:cp='http://schemas.openxmlformats.org/package/2006/metadata/core-properties'"
            " xmlns:dc='http://purl.org/dc/elements/1.1/'>"
            "<dc:title>Deck</dc:title><dc:creator>Jane Roe</dc:creator>"
            "<cp:lastModifiedBy>John Doe</cp:lastModifiedBy></cp:coreProperties>"))
        z.writestr("docProps/app.xml", (
            "<Properties xmlns='http://schemas.openxmlformats.org/officeDocument/2006/extended-properties'>"
            "<Manager>Jane Roe</Manager><Company>Acme</Company></Properties>"))
        z.writestr("ppt/changesInfos/changesInfo1.xml", (
            "<pc:chgInfo xmlns:pc='urn:pc'><pc:chgData name='John Doe' userId='1234-abcd-userid'"
            " providerId='ADAL'/><pc:sldChg name='Results slide'/></pc:chgInfo>"))
        z.writestr("ppt/commentAuthors.xml", (
            "<p:cmAuthorLst xmlns:p='urn:p'><p:cmAuthor id='1' name='Jane Roe' initials='JR'/>"
            "</p:cmAuthorLst>"))
        z.writestr("ppt/authors.xml", (
            "<p188:authorLst xmlns:p188='urn:a'><p188:author id='{A}' name='Jane Roe' initials='JR'"
            " userId='jane@example.invalid' providerId='AD'/></p188:authorLst>"))
        z.writestr("ppt/slides/slide1.xml", "<p:sld xmlns:p='urn:p'><p:cNvPr id='2' name='Title 1'/></p:sld>")
    return buf.getvalue()


def test_personal_identity_fields_masked():
    parts = normalize.normalize_pptx(_identity_deck())
    everything = "".join(parts.values())
    people = ("Jane Roe", "John Doe", "JR", "1234-abcd-userid", "jane@example.invalid")
    assert not [p for p in people if p in everything]
    assert everything.count("<PERSON>") == 10
    # Shape names, titles and the company are not personal and stay, and so
    # does a "name" on any changesInfo element other than chgData.
    assert 'name="Title 1"' in parts["ppt/slides/slide1.xml"]
    assert 'name="Results slide"' in parts["ppt/changesInfos/changesInfo1.xml"]
    assert "<dc:title>Deck</dc:title>" in parts["docProps/core.xml"]
    assert "<Company>Acme</Company>" in parts["docProps/app.xml"]

    wb = openpyxl.Workbook()
    wb.properties.creator, wb.properties.lastModifiedBy = "Jane Roe", "John Doe"
    buf = io.BytesIO()
    wb.save(buf)
    core = normalize.normalize_xlsx(buf.getvalue())["docProps/core.xml"]
    assert "Jane Roe" not in core and "John Doe" not in core and core.count("<PERSON>") == 2
