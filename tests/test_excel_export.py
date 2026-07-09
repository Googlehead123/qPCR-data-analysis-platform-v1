"""Regression tests for the monolith's Excel export (native openpyxl charts).

Focus: gene -> chart Y-axis-title pairing must stay correct when there are 10+
genes. A lexicographic sort of chart file names (chart1, chart10, chart11, ...,
chart2) previously mispaired gene names onto the wrong chart's axis title; the
fix sorts the chart files numerically (creation order).
"""

import io
import re
import zipfile
from importlib import import_module

import pandas as pd


def _make_processed(genes):
    """Minimal per-gene processed frame with the columns the chart builder reads."""
    out = {}
    for i, g in enumerate(genes):
        out[g] = pd.DataFrame({
            "Target": [g, g],
            "Condition": ["Non-treated", "Treatment"],
            "Group": ["Negative Control", "Treatment"],
            "Fold_Change": [1.0, 2.0 + i],
            "Relative_Expression": [1.0, 2.0 + i],
            "SEM": [0.05, 0.1],
        })
    return out


def _make_raw(genes):
    rows = []
    for g in genes:
        for s in ("Non-treated", "Treatment"):
            rows.append({"Well": "A1", "Sample": s, "Target": g, "CT": 20.0})
    return pd.DataFrame(rows)


def _yaxis_gene_of_chart(chart_xml, gene_names):
    """Return which gene name appears in the chart's valAx (Y-axis) title, or None."""
    C = "http://schemas.openxmlformats.org/drawingml/2006/chart"
    A = "http://schemas.openxmlformats.org/drawingml/2006/main"
    from lxml import etree

    root = etree.fromstring(chart_xml)
    val_ax = root.find(".//c:valAx", {"c": C})
    if val_ax is None:
        return None
    title = val_ax.find(".//c:title", {"c": C})
    if title is None:
        return None
    texts = "".join(t.text or "" for t in title.findall(f".//{{{A}}}t"))
    for g in gene_names:
        if g in texts:
            return g
    return None


def test_chart_error_bars_are_asymmetric_fold_change(mock_streamlit):
    """Excel native chart error bars must use asymmetric FC bounds (plus != minus
    columns), matching the Plotly graph — not symmetric SEM in both directions."""
    spec = import_module("streamlit qpcr analysis v1")

    genes = ["COL1A1", "MMP1"]
    processed = _make_processed(genes)
    # Give distinct upper/lower fold-change error bounds.
    for g in genes:
        processed[g]["FC_Error_Upper"] = [0.30, 0.55]
        processed[g]["FC_Error_Lower"] = [0.20, 0.40]
    raw = _make_raw(genes)
    params = {"Housekeeping_Gene": "GAPDH", "Efficacy_Type": "Anti-Aging"}
    mapping = {"Non-treated": {"condition": "Non-treated", "group": "Negative Control"},
               "Treatment": {"condition": "Treatment", "group": "Treatment"}}

    xlsx = spec.export_to_excel(raw, processed, params, mapping)
    xlsx_bytes = xlsx.getvalue() if hasattr(xlsx, "getvalue") else xlsx

    C = "http://schemas.openxmlformats.org/drawingml/2006/chart"
    from lxml import etree
    zf = zipfile.ZipFile(io.BytesIO(xlsx_bytes))
    chart_files = [n for n in zf.namelist() if re.match(r"xl/charts/chart\d+\.xml$", n)]
    assert chart_files, "no chart XML found"

    checked = 0
    for cf in chart_files:
        root = etree.fromstring(zf.read(cf))
        errbars = root.find(".//c:errBars", {"c": C})
        if errbars is None:
            continue
        plus_f = errbars.find(".//c:plus//c:f", {"c": C})
        minus_f = errbars.find(".//c:minus//c:f", {"c": C})
        assert plus_f is not None and minus_f is not None, "errBars missing plus/minus refs"
        # Asymmetric: the two directions must reference DIFFERENT columns
        # (the old bug pointed both at the single SEM column).
        assert plus_f.text != minus_f.text, (
            f"error bars are symmetric (plus==minus: {plus_f.text}); "
            f"expected distinct FC_Err_Upper/Lower columns"
        )
        checked += 1
    assert checked > 0, "no gene chart with error bars was verified"


def test_chart_axis_titles_match_genes_beyond_ten(mock_streamlit):
    """With 12 genes, each chart (in numeric file order) must carry the correct
    gene name — i.e. creation order is preserved, not lexicographic order."""
    spec = import_module("streamlit qpcr analysis v1")

    genes = [f"GENE{i:02d}" for i in range(1, 13)]  # GENE01 .. GENE12
    processed = _make_processed(genes)
    raw = _make_raw(genes)
    params = {"Housekeeping_Gene": "GAPDH", "Efficacy_Type": "Anti-Aging"}
    mapping = {"Non-treated": {"condition": "Non-treated", "group": "Negative Control"},
               "Treatment": {"condition": "Treatment", "group": "Treatment"}}

    xlsx = spec.export_to_excel(raw, processed, params, mapping)
    xlsx_bytes = xlsx.getvalue() if hasattr(xlsx, "getvalue") else xlsx

    zf = zipfile.ZipFile(io.BytesIO(xlsx_bytes))
    chart_files = sorted(
        (n for n in zf.namelist() if re.match(r"xl/charts/chart\d+\.xml$", n)),
        key=lambda n: int(re.search(r"(\d+)", n.rsplit("/", 1)[-1]).group(1)),
    )
    assert len(chart_files) == 12, f"expected 12 gene charts, got {len(chart_files)}"

    titles_in_numeric_order = [_yaxis_gene_of_chart(zf.read(f), genes) for f in chart_files]
    # Charts are created one-per-gene in processed_data order, so reading charts
    # in numeric file order must yield the gene names in that same order.
    assert titles_in_numeric_order == genes, (
        f"chart/gene pairing scrambled:\n  got:      {titles_in_numeric_order}\n"
        f"  expected: {genes}"
    )
