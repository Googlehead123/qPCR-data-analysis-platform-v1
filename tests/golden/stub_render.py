"""A Chrome-free stand-in for Kaleido, for the CI goldens.

The stub returns a solid PNG whose colour is the first three bytes of
sha256(canonical figure JSON), so any change to a chart changes the slide and
image media even though nothing is really drawn. It is sized as Kaleido sizes
its output, (``width`` argument, else ``layout.width``, else 64) x ``scale``, so
a change to the requested size or scale changes the media too. A non-PNG
format gets a small tagged byte string instead, so a format switch shows.
"""

import hashlib
import importlib
import io

from PIL import Image

from tests.golden.normalize import normalize_figures

APP_MODULE = "streamlit qpcr analysis v1"
DEFAULT_SIZE = 64


def stub_image(fig, fmt="png", scale=1, width=None, height=None) -> bytes:
    canonical = normalize_figures({"fig": fig})["fig"]
    colour = hashlib.sha256(canonical.encode("utf-8")).digest()[:3]
    w = int(width or fig.layout.width or DEFAULT_SIZE) * int(scale or 1)
    h = int(height or fig.layout.height or DEFAULT_SIZE) * int(scale or 1)
    if str(fmt).lower() != "png":
        return f"<stub {fmt} {w}x{h} #{colour.hex()}>".encode("ascii")
    buf = io.BytesIO()
    Image.new("RGB", (w, h), tuple(colour)).save(buf, format="PNG")
    return buf.getvalue()


def _export_figure_to_bytes(fig, fmt="png", scale=2, width=None, height=None):
    return stub_image(fig, fmt, scale, width, height)


def _fig_to_image(fig, format="png", scale=2, width=None, height=None):  # noqa: A002 - app's signature
    return stub_image(fig, format, scale, width, height)


def install_stub(monkeypatch_like) -> None:
    """Route every figure render in ``qpcr.export_utils`` and the loaded app through the stub."""
    export_utils = importlib.import_module("qpcr.export_utils")
    app = importlib.import_module(APP_MODULE)
    monkeypatch_like.setattr(export_utils, "export_figure_to_bytes", _export_figure_to_bytes)
    # The app imported the name directly, so its module attribute needs patching too.
    monkeypatch_like.setattr(app, "export_figure_to_bytes", _export_figure_to_bytes)
    monkeypatch_like.setattr(app.ReportGenerator, "_fig_to_image", staticmethod(_fig_to_image))
