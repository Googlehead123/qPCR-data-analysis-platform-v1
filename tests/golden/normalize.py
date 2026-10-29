"""Semantic normalizers: turn writer output into stable, diffable text parts.

Only the software-version token is masked (and the OOXML created/modified
stamps in ``docProps/core.xml``). Time is frozen by the driver, so a user
string that looks like a date must survive normalization untouched.

Masking is applied to string leaves (cell values, XML text and attributes,
frame cells, figure strings) rather than to serialized text, because the
version regex stops only at whitespace or parentheses and would otherwise eat
a closing quote, comma or tag.
"""

import datetime as dt
import hashlib
import io
import json
import re
import zipfile
from typing import Any

import numpy as np
import openpyxl
import pandas as pd
import plotly.graph_objects as go
from lxml import etree
from PIL import Image

APP_VERSION_RE = re.compile(r"qPCR Analysis Suite [^\s()]+( \(commit [^)]+\))?")
APP_VERSION_TOKEN = "<APP_VERSION>"
TS_TOKEN = "<TS>"
_FLOAT_FMT = "%.12g"
_CORE_XML = "docProps/core.xml"
_TS_TAGS = (
    "{http://purl.org/dc/terms/}created",
    "{http://purl.org/dc/terms/}modified",
)
_IMAGE_SUFFIXES = (".png", ".jpg", ".jpeg", ".gif", ".bmp", ".tif", ".tiff", ".webp")
_XML_SUFFIXES = (".xml", ".rels", ".vml")


def mask_text(text: str) -> str:
    """Replace every software-version token in ``text`` with ``<APP_VERSION>``."""
    return APP_VERSION_RE.sub(APP_VERSION_TOKEN, text)


def _round_float(value: float) -> float:
    return float(_FLOAT_FMT % value)


def _fmt_cell(value: Any) -> str:
    if isinstance(value, float):
        text = _FLOAT_FMT % value
    elif isinstance(value, dt.datetime):
        text = value.isoformat(sep=" ")
    elif isinstance(value, (dt.date, dt.time)):
        text = value.isoformat()
    else:
        text = str(value)
    return mask_text(text).replace("\r", "\\r").replace("\n", "\\n").replace("\t", "\\t")


def _pretty_xml(data: bytes, *, mask_timestamps: bool = False) -> str:
    parser = etree.XMLParser(remove_blank_text=True, resolve_entities=False)
    root = etree.fromstring(data, parser)
    for el in root.iter():
        if not isinstance(el.tag, str):
            continue
        if mask_timestamps and el.tag in _TS_TAGS:
            el.text = "@@TS@@"
            continue
        if el.text:
            el.text = mask_text(el.text)
        if el.tail:
            el.tail = mask_text(el.tail)
        for key, val in list(el.attrib.items()):
            el.set(key, mask_text(val))
    out = etree.tostring(root, pretty_print=True, encoding="unicode")
    # Tokens contain "<", which XML serialization escapes; restore them literally.
    return out.replace("@@TS@@", TS_TOKEN).replace("&lt;APP_VERSION&gt;", APP_VERSION_TOKEN)


def _binary_summary(data: bytes) -> str:
    return f"<binary {len(data)} bytes sha256:{hashlib.sha256(data).hexdigest()}>"


def _image_summary(data: bytes) -> str:
    """``<image WxH sha256:...>`` over decoded RGBA pixels; raw-bytes hash if undecodable."""
    try:
        with Image.open(io.BytesIO(data)) as img:
            rgba = img.convert("RGBA")
            digest = hashlib.sha256(rgba.tobytes()).hexdigest()
            return f"<image {rgba.width}x{rgba.height} sha256:{digest}>"
    except Exception:  # noqa: BLE001 - non-raster media (emf, svg) falls back to raw bytes
        return _binary_summary(data)


def _sheet_text(ws) -> str:
    lines = [
        f"{cell.coordinate}\t{_fmt_cell(cell.value)}"
        for row in ws.iter_rows()
        for cell in row
        if cell.value is not None
    ]
    return "\n".join(lines) + ("\n" if lines else "")


def normalize_xlsx(data: bytes) -> dict[str, str]:
    """One part per sheet (``sheet:<name>``) and chart (``chart:<part name>``), plus other parts.

    Worksheet XML and sharedStrings are replaced by the per-cell lines; every
    other XML part is kept, pretty-printed; media is pixel-hashed.
    """
    parts: dict[str, str] = {}
    wb = openpyxl.load_workbook(io.BytesIO(data), data_only=False)
    for ws in wb.worksheets:
        parts[f"sheet:{ws.title}"] = _sheet_text(ws)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for name in sorted(z.namelist()):
            if name.startswith("xl/worksheets/") and name.endswith(".xml"):
                continue
            if name == "xl/sharedStrings.xml":
                continue
            raw = z.read(name)
            if name.startswith("xl/charts/") and name.endswith(".xml"):
                parts[f"chart:{name}"] = _pretty_xml(raw)
            elif name.lower().endswith(_IMAGE_SUFFIXES):
                parts[name] = _image_summary(raw)
            elif name.endswith(_XML_SUFFIXES):
                parts[name] = _pretty_xml(raw, mask_timestamps=(name == _CORE_XML))
            else:
                parts[name] = _binary_summary(raw)
    return parts


def normalize_pptx(data: bytes) -> dict[str, str]:
    """One part per slide XML (pretty-printed) and one per media file (pixel hash)."""
    parts: dict[str, str] = {}
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for name in sorted(z.namelist()):
            if name.startswith("ppt/slides/") and name.endswith(".xml"):
                parts[name] = _pretty_xml(z.read(name))
            elif name.startswith("ppt/media/"):
                parts[name] = _image_summary(z.read(name))
    return parts


def normalize_frames(frames: dict[str, pd.DataFrame]) -> dict[str, str]:
    """CSV per frame with 12-significant-digit floats and masked version strings."""
    out: dict[str, str] = {}
    for name, df in frames.items():
        masked = df.copy()
        for col in masked.columns:
            if masked[col].dtype == object:
                masked[col] = masked[col].map(
                    lambda v: mask_text(v) if isinstance(v, str) else v
                )
        out[name] = masked.to_csv(float_format=_FLOAT_FMT)
    return out


def _canonical(obj: Any) -> Any:
    """Plain-JSON form of a plotly structure: floats to 12 digits, arrays to lists."""
    if isinstance(obj, dict):
        return {str(k): _canonical(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_canonical(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_canonical(v) for v in obj.tolist()]
    if isinstance(obj, np.generic):
        return _canonical(obj.item())
    if isinstance(obj, bool) or obj is None or isinstance(obj, int):
        return obj
    if isinstance(obj, float):
        return _round_float(obj) if np.isfinite(obj) else obj
    if isinstance(obj, str):
        return mask_text(obj)
    if isinstance(obj, (dt.datetime, dt.date)):
        return obj.isoformat()
    return mask_text(str(obj))


def normalize_figures(figs: dict[str, go.Figure]) -> dict[str, str]:
    """Canonical sorted-key JSON of each figure, floats at 12 significant digits."""
    return {
        name: json.dumps(_canonical(fig.to_plotly_json()), sort_keys=True)
        for name, fig in figs.items()
    }
