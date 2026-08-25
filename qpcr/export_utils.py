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

import hashlib
import io
import os
import shutil
import threading
import zipfile
from collections import OrderedDict

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

# ---------------------------------------------------------------------------
# Rendered-image memo cache
# ---------------------------------------------------------------------------
# Every render launches headless Chrome (~1-3 s of pinned CPU). Exports re-render
# the SAME figures repeatedly: a PPT run, an Excel run and a ZIP bundle each walk
# st.session_state.graphs, and Streamlit re-executes the script on every rerun.
# Rendering is deterministic — identical figure JSON + identical parameters always
# yield identical bytes — so memoizing on that pair is exact, never stale.
#
# Bounded twice (entry count AND total bytes) because Streamlit Cloud caps RAM at
# ~2.7 GB; eviction is least-recently-used.
_RENDER_CACHE_MAX_ENTRIES = 32
_RENDER_CACHE_MAX_BYTES = 64 * 1024 * 1024  # 64 MB
_render_cache: "OrderedDict[str, bytes]" = OrderedDict()
_render_cache_bytes = 0
_render_cache_lock = threading.Lock()


def _render_cache_key(fig, kwargs: dict) -> str | None:
    """Content key for a render, or ``None`` if the figure can't be serialized.

    ``None`` disables caching for that call rather than risking a wrong hit.
    """
    try:
        payload = fig.to_json()
    except Exception:  # noqa: BLE001 — exotic figure: skip the cache, still render
        return None
    if payload is None:
        return None
    h = hashlib.md5(payload.encode("utf-8"))
    for k in sorted(kwargs):
        h.update(f"|{k}={kwargs[k]}".encode("utf-8"))
    return h.hexdigest()


def _render_cache_get(key: str | None) -> bytes | None:
    if key is None:
        return None
    with _render_cache_lock:
        hit = _render_cache.get(key)
        if hit is not None:
            _render_cache.move_to_end(key)  # mark most-recently-used
        return hit


def _render_cache_put(key: str | None, value: bytes) -> None:
    if key is None or not isinstance(value, (bytes, bytearray)):
        return
    global _render_cache_bytes
    with _render_cache_lock:
        if key in _render_cache:
            _render_cache_bytes -= len(_render_cache.pop(key))
        _render_cache[key] = bytes(value)
        _render_cache_bytes += len(value)
        while _render_cache and (
            len(_render_cache) > _RENDER_CACHE_MAX_ENTRIES
            or _render_cache_bytes > _RENDER_CACHE_MAX_BYTES
        ):
            _, evicted = _render_cache.popitem(last=False)
            _render_cache_bytes -= len(evicted)


def clear_render_cache() -> None:
    """Drop every memoized render (used by tests; safe to call at any time)."""
    global _render_cache_bytes
    with _render_cache_lock:
        _render_cache.clear()
        _render_cache_bytes = 0


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
    # Phrases, not bare words. "close" alone also matched "closest", so a plain
    # Plotly property error — hovermode: 'closest' is not a valid value, or an
    # invalid 'closest' property on a Bar — triggered the ~150 MB Chrome
    # download and then re-raised the same unrelated error. "not found" and
    # "timeout" stay as-is: both only appear here in browser-launch messages.
    return any(
        keyword in msg
        for keyword in (
            "chrome",
            "chromium",
            "browser",
            "kaleido",
            "choreographer",
            "closed unexpectedly",
            "closed the connection",
            "connection closed",
            "not found",
            "browser_path",
            "executable",
            "timeout",
        )
    )


_kaleido_server_started = False


def _ensure_kaleido_server() -> None:
    """Start Kaleido's sync server once per process, lazily.

    Without it every ``fig.to_image`` launches its own headless Chrome, so an
    N-gene panel pays N browser launches. Measured back to back on the same box,
    7 figures: 21.58s total (mean 3.08s, first 5.54s) versus 3.33s with the
    server held open (mean 0.48s, steady state 0.17s) — a 6.5x speedup, with
    byte-identical output and identical dimensions. A 14-gene panel was ~45s of
    pinned CPU, on a Cloud container with far less CPU than the box that
    measured it.

    Held open for the process rather than started and stopped per batch
    (decision: Min, 2026-08-24, item 12). One Chrome for the app's lifetime is
    the cost; Cloud caps RAM around 2.7 GB, so if that ever bites, the
    per-batch alternative recovers most of the win.

    Best-effort by design: this is a speedup, not a requirement. Any failure
    leaves the per-figure launch path working exactly as before, so the render
    still succeeds — which is why the flag is set even when the call raises.
    """
    global _kaleido_server_started
    if _kaleido_server_started:
        return
    _kaleido_server_started = True
    try:
        import atexit
        import kaleido

        start = getattr(kaleido, "start_sync_server", None)
        if start is None:
            return  # older kaleido: nothing to hold open
        # silence_warnings: kaleido's server is a process singleton while this
        # flag is module state, so a module reload (which the test suite does)
        # would re-enter here and warn "Server already open". An open server is
        # exactly the state we want, so that warning is noise.
        try:
            start(silence_warnings=True)
        except TypeError:
            start()  # older signature without the keyword
        stop = getattr(kaleido, "stop_sync_server", None)
        if stop is not None:
            atexit.register(_stop_kaleido_server_quietly, stop)
    except Exception:  # noqa: BLE001
        # Falls back to one browser launch per figure. Slow, still correct.
        pass


def _stop_kaleido_server_quietly(stop) -> None:
    """Shut the sync server down at exit without noise on the way out."""
    try:
        try:
            stop(silence_warnings=True)
        except TypeError:
            stop()
    except Exception:  # noqa: BLE001
        pass


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

    # Identical figure + identical parameters => identical bytes. Serving that
    # from memory skips a headless-Chrome launch entirely.
    cache_key = _render_cache_key(fig, kwargs)
    cached = _render_cache_get(cache_key)
    if cached is not None:
        return cached

    with _export_lock:
        _ensure_browser_path()
        _ensure_kaleido_server()

    try:
        out = fig.to_image(**kwargs)
        _render_cache_put(cache_key, out)
        return out
    except Exception as first_err:  # noqa: BLE001 — need broad catch to classify
        if not _looks_like_browser_error(first_err):
            raise
        # Fallback: fetch a known-good Chrome and retry exactly once.
        try:
            with _export_lock:
                chrome = _download_fallback_chrome()
                os.environ["BROWSER_PATH"] = chrome
            out = fig.to_image(**kwargs)
            _render_cache_put(cache_key, out)
            return out
        except Exception as second_err:  # noqa: BLE001
            resolved = os.environ.get("BROWSER_PATH", "(none found)")
            raise RuntimeError(
                "Image export requires a headless Chrome/Chromium that can launch.\n"
                f"• System browser tried: {resolved}\n"
                f"• Automatic Chrome download also failed: {second_err}\n"
                # No "Download Interactive HTML" alternative is offered: there
                # is no such export anywhere in the app, so pointing at it sent
                # the operator looking for a button that does not exist.
                "Fixes: run `plotly_get_chrome` in the app environment, or on "
                "Streamlit Cloud keep 'chromium' in packages.txt. The Excel "
                "report does not need Chrome, so it remains available."
            ) from second_err
