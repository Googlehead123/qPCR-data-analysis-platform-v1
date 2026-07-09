"""Robust Plotly image export for Kaleido v1 / headless Chromium.

Plotly 6.x renders static images (PNG/SVG/PDF) through Kaleido v1, which drives a
real headless Chrome/Chromium over the DevTools protocol via ``choreographer``.
That is fragile across environments (local WSL, Streamlit Cloud, containers): the
browser may be absent, or may "close immediately after starting" (the classic
symptom of a broken distro Chromium or a stale ``BROWSER_PATH``).

This module centralizes ALL image export so every call site gets identical,
hardened behaviour:

1. Deterministic browser selection — prefer real Google Chrome over distro
   Chromium (Chrome launches headless far more reliably), and only honour a
   pre-existing ``BROWSER_PATH`` if it actually exists on disk.
2. A one-time, cached fallback that downloads a known-good Chrome-for-Testing
   (``plotly.io.get_chrome`` / ``plotly_get_chrome``) when the system browser is
   missing or refuses to launch. choreographer auto-discovers that download on
   subsequent calls.
3. A clear, actionable ``RuntimeError`` instead of an opaque choreographer
   traceback when everything fails.

``choreographer`` reads the browser location from the ``BROWSER_PATH`` env var
(NOT ``CHROME_PATH``); see the choreographer source ``utils/_which.py``.
"""

from __future__ import annotations

import io
import os
import shutil
import threading
import zipfile

# Preference order: real Chrome first (most reliable headless), then Chromium.
_CHROME_EXE_NAMES = (
    "google-chrome-stable",
    "google-chrome",
    "chrome",
    "chromium",
    "chromium-browser",
)

# Common absolute locations, probed only if the names above aren't on PATH.
_CHROME_COMMON_PATHS = (
    "/usr/bin/google-chrome-stable",
    "/usr/bin/google-chrome",
    "/opt/google/chrome/chrome",
    "/usr/bin/chromium",
    "/usr/bin/chromium-browser",
    "/snap/bin/chromium",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
)

_export_lock = threading.Lock()
_downloaded_chrome_path: str | None = None  # cache for plotly_get_chrome result


def _find_system_browser() -> str | None:
    """Return the path to a usable Chrome/Chromium executable, or ``None``."""
    for name in _CHROME_EXE_NAMES:
        found = shutil.which(name)
        if found:
            return found
    for candidate in _CHROME_COMMON_PATHS:
        if os.path.exists(candidate):
            return candidate
    return None


def _ensure_browser_path() -> str | None:
    """Point ``BROWSER_PATH`` at an existing browser for choreographer.

    A pre-existing ``BROWSER_PATH`` is respected only if it still exists on disk;
    a stale/non-existent override is replaced (or dropped so choreographer can run
    its own search). Returns the resolved path, or ``None`` if none was found.
    """
    current = os.environ.get("BROWSER_PATH")
    if current and os.path.exists(current):
        return current
    browser = _find_system_browser()
    if browser:
        os.environ["BROWSER_PATH"] = browser
    elif current:
        # Stale override that no longer exists — drop it so choreographer's own
        # discovery (incl. any previously downloaded Chrome) can take over.
        os.environ.pop("BROWSER_PATH", None)
    return browser


def _download_fallback_chrome() -> str:
    """Download (once, cached) a known-good Chrome-for-Testing; return its path.

    Uses ``plotly.io.get_chrome`` which fetches into choreographer's browser_exe
    cache. Subsequent renders find it automatically. Raises on failure.
    """
    global _downloaded_chrome_path
    if _downloaded_chrome_path and os.path.exists(_downloaded_chrome_path):
        return _downloaded_chrome_path
    import plotly.io as pio

    path = str(pio.get_chrome())
    _downloaded_chrome_path = path
    return path


def build_zip(files: dict) -> bytes:
    """Bundle ``{filename: content}`` into a ZIP archive, returned as bytes.

    Accepts ``bytes``, ``str`` (encoded UTF-8), or objects with ``.getvalue()``
    (e.g. ``io.BytesIO``). Entries whose content is ``None`` are skipped, so a
    partially-failed report bundle still yields a usable archive.
    """
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, data in files.items():
            if data is None:
                continue
            if hasattr(data, "getvalue"):
                data = data.getvalue()
            if isinstance(data, str):
                data = data.encode("utf-8")
            zf.writestr(name, data)
    return buf.getvalue()


def _looks_like_browser_error(err: Exception) -> bool:
    """Heuristic: does this exception look like a browser/Kaleido launch failure?

    Guards the (expensive, ~150 MB) download fallback so it only triggers for
    browser problems, not unrelated rendering errors.
    """
    msg = str(err).lower()
    return any(
        keyword in msg
        for keyword in (
            "chrome",
            "chromium",
            "browser",
            "kaleido",
            "choreographer",
            "close",
            "not found",
            "browser_path",
            "executable",
            "timeout",
        )
    )


def export_figure_to_bytes(fig, fmt: str = "png", scale: int = 2,
                           width: int | None = None, height: int | None = None) -> bytes:
    """Render a Plotly figure to image bytes, robust across environments.

    Tries the system browser first; on a browser-launch failure, downloads a
    known-good Chrome once and retries. Raises ``RuntimeError`` with actionable
    guidance if both attempts fail.

    Args:
        fig: a ``plotly.graph_objects.Figure``.
        fmt: ``"png"``, ``"svg"``, or ``"pdf"``.
        scale: pixel scale factor (higher = higher DPI); ignored for vector fmts.
        width / height: optional pixel dimensions.

    Returns:
        The encoded image as ``bytes``.
    """
    kwargs: dict = {"format": fmt, "scale": scale}
    if width is not None:
        kwargs["width"] = width
    if height is not None:
        kwargs["height"] = height

    with _export_lock:
        _ensure_browser_path()

    try:
        return fig.to_image(**kwargs)
    except Exception as first_err:  # noqa: BLE001 — need broad catch to classify
        if not _looks_like_browser_error(first_err):
            raise
        # Fallback: fetch a known-good Chrome and retry exactly once.
        try:
            with _export_lock:
                chrome = _download_fallback_chrome()
                os.environ["BROWSER_PATH"] = chrome
            return fig.to_image(**kwargs)
        except Exception as second_err:  # noqa: BLE001
            resolved = os.environ.get("BROWSER_PATH", "(none found)")
            raise RuntimeError(
                "Image export requires a headless Chrome/Chromium that can launch.\n"
                f"• System browser tried: {resolved}\n"
                f"• Automatic Chrome download also failed: {second_err}\n"
                "Fixes: run `plotly_get_chrome` in the app environment, or on "
                "Streamlit Cloud keep 'chromium' in packages.txt. You can still use "
                "the 'Download Interactive HTML' export as an alternative."
            ) from second_err
