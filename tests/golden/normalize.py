"""Semantic normalizers: turn writer output into stable, diffable text parts.

Only the software-version token is masked, plus the OOXML created/modified
stamps in ``docProps/core.xml`` and the personal identities a template carries
(author, last editor, reviewer names and account ids) as ``<PERSON>``. Time is
frozen by the driver, so a user string that looks like a date must survive
normalization untouched.

Masking is applied to string leaves (cell values, XML text and attributes,
frame cells, figure strings) rather than to serialized text, because the
version regex stops only at whitespace or parentheses and would otherwise eat
a closing quote, comma or tag.
"""

import base64
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
_SHEET_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
PERSON_TOKEN = "<PERSON>"
# Elements whose text is a person, wherever they occur (core.xml, app.xml).
_PERSON_TAGS = (
    "{http://purl.org/dc/elements/1.1/}creator",
    "{http://schemas.openxmlformats.org/package/2006/metadata/core-properties}lastModifiedBy",
    "{http://schemas.openxmlformats.org/officeDocument/2006/extended-properties}Manager",
)
# Parts that list people: their name/initials/account-id attributes are masked.
# Elsewhere "name" is a shape or sheet name and is kept.
_PERSON_PARTS = re.compile(
    r"^(ppt/changesInfos/.*|ppt/authors\.xml|ppt/commentAuthors\.xml|xl/persons/.*)$"
)
# The identity attributes of the elements that carry a person in those parts.
# Any other attribute, a "name" on another element included, is kept.
_PERSON_ATTRS = {
    "chgData": ("name", "userId"),               # ppt/changesInfos/*
    "author": ("name", "initials", "userId"),     # ppt/authors.xml
    "cmAuthor": ("name", "initials"),             # ppt/commentAuthors.xml
    "person": ("displayName", "userId"),          # xl/persons/*
}


def mask_text(text: str) -> str:
    """Replace every software-version token in ``text`` with ``<APP_VERSION>``."""
    return APP_VERSION_RE.sub(APP_VERSION_TOKEN, text)


def _round_float(value: float) -> float:
    # + 0.0 turns -0.0 into 0.0. A value straddling a 12th-digit rounding
    # boundary can still round two ways across platforms; that is accepted
    # (Ruling 10).
    return float(_FLOAT_FMT % value) + 0.0


def _fmt_cell(value: Any) -> str:
    if isinstance(value, float):
        text = _FLOAT_FMT % (value + 0.0)
        if text.lstrip("-").isdigit():
            text += ".0"  # a float cell must not read like an int cell
    elif isinstance(value, dt.datetime):
        text = value.isoformat(sep=" ")
    elif isinstance(value, (dt.date, dt.time)):
        text = value.isoformat()
    else:
        text = str(value)
    return mask_text(text).replace("\r", "\\r").replace("\n", "\\n").replace("\t", "\\t")


def _pretty_xml(data: bytes, *, mask_timestamps: bool = False, mask_person_attrs: bool = False) -> str:
    parser = etree.XMLParser(remove_blank_text=True, resolve_entities=False)
    root = etree.fromstring(data, parser)
    for el in root.iter():
        if not isinstance(el.tag, str):
            continue
        if mask_timestamps and el.tag in _TS_TAGS:
            el.text = "@@TS@@"
            continue
        if el.tag in _PERSON_TAGS and el.text:
            el.text = "@@PERSON@@"
        if el.text:
            el.text = mask_text(el.text)
        if el.tail:
            el.tail = mask_text(el.tail)
        identity = _PERSON_ATTRS.get(etree.QName(el).localname, ()) if mask_person_attrs else ()
        for key, val in list(el.attrib.items()):
            person = etree.QName(key).localname in identity
            el.set(key, "@@PERSON@@" if person else mask_text(val))
    out = etree.tostring(root, pretty_print=True, encoding="unicode")
    # Tokens contain "<", which XML serialization escapes; restore them literally.
    return (
        out.replace("@@TS@@", TS_TOKEN)
        .replace("@@PERSON@@", PERSON_TOKEN)
        .replace("&lt;APP_VERSION&gt;", APP_VERSION_TOKEN)
    )


def _sheet_layout_xml(data: bytes) -> str:
    """Worksheet XML without cell values: types, styles, merges, widths, panes stay.

    Values are already one line per cell in ``sheet:<name>``; dropping them here
    keeps a value change out of this part while a format change still shows.
    """
    parser = etree.XMLParser(remove_blank_text=True, resolve_entities=False)
    root = etree.fromstring(data, parser)
    for cell in root.iter(f"{{{_SHEET_NS}}}c"):
        for child in list(cell):
            if child.tag in (f"{{{_SHEET_NS}}}v", f"{{{_SHEET_NS}}}is"):
                cell.remove(child)
    return _pretty_xml(etree.tostring(root))


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

    Cell values become the per-cell ``sheet:<name>`` lines, so sharedStrings is
    dropped and each worksheet's XML is kept without its values (cell types,
    style indices, merges, column widths, panes). Every other XML part is kept,
    pretty-printed; media is pixel-hashed.
    """
    parts: dict[str, str] = {}
    wb = openpyxl.load_workbook(io.BytesIO(data), data_only=False)
    for ws in wb.worksheets:
        parts[f"sheet:{ws.title}"] = _sheet_text(ws)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for name in sorted(z.namelist()):
            if name == "xl/sharedStrings.xml":
                continue
            raw = z.read(name)
            if name.startswith("xl/worksheets/") and name.endswith(".xml"):
                parts[name] = _sheet_layout_xml(raw)
            elif name.startswith("xl/charts/") and name.endswith(".xml"):
                parts[f"chart:{name}"] = _pretty_xml(raw)
            else:
                parts[name] = _package_part(name, raw)
    return parts


def _package_part(name: str, raw: bytes) -> str:
    if name.lower().endswith(_IMAGE_SUFFIXES):
        return _image_summary(raw)
    if name.endswith(_XML_SUFFIXES):
        return _pretty_xml(
            raw,
            mask_timestamps=(name == _CORE_XML),
            mask_person_attrs=bool(_PERSON_PARTS.match(name)),
        )
    return _binary_summary(raw)


def normalize_pptx(data: bytes) -> dict[str, str]:
    """Every package part: XML (slides, presentation, rels, layouts) pretty-printed, media pixel-hashed.

    Slide order and size live in ``ppt/presentation.xml`` and which picture a
    slide shows lives in its ``.rels``, so neither is safe to drop.
    """
    parts: dict[str, str] = {}
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        for name in sorted(z.namelist()):
            parts[name] = _package_part(name, z.read(name))
    return parts


def normalize_frames(frames: dict[str, pd.DataFrame]) -> dict[str, str]:
    """CSV per frame with 12-significant-digit floats (-0 as 0) and masked version strings."""
    out: dict[str, str] = {}
    for name, df in frames.items():
        masked = df.copy()
        for col in masked.columns:
            if masked[col].dtype == object:
                masked[col] = masked[col].map(
                    lambda v: mask_text(v) if isinstance(v, str) else v
                )
            elif masked[col].dtype.kind == "f":
                masked[col] = masked[col] + 0.0
        out[name] = masked.to_csv(float_format=_FLOAT_FMT)
    return out


def _decode_typed_array(obj: dict) -> np.ndarray:
    """Plotly 6 typed array ``{"dtype", "bdata"[, "shape"]}`` back to an ndarray."""
    raw = obj["bdata"]
    buf = base64.b64decode(raw) if isinstance(raw, str) else bytes(raw)
    arr = np.frombuffer(buf, dtype=np.dtype(obj["dtype"]))
    shape = obj.get("shape")
    if shape:
        dims = [int(s) for s in shape.split(",")] if isinstance(shape, str) else list(shape)
        arr = arr.reshape(dims)
    return arr


def _canonical(obj: Any) -> Any:
    """Plain-JSON form of a plotly structure: floats to 12 digits, arrays to lists."""
    if isinstance(obj, dict) and "bdata" in obj and "dtype" in obj:
        return _canonical(_decode_typed_array(obj))
    if isinstance(obj, dict):
        return {str(k): _canonical(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_canonical(v) for v in obj]
    if isinstance(obj, np.ndarray):
        if obj.dtype.kind == "M":  # .tolist() of datetime64[ns] gives epoch ints
            return np.datetime_as_string(obj).tolist()
        return [_canonical(v) for v in obj.tolist()]
    if isinstance(obj, np.datetime64):
        return str(np.datetime_as_string(obj))
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
