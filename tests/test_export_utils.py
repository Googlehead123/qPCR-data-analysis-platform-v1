"""Tests for qpcr.export_utils — robust Plotly image export.

Covers browser resolution, the happy render path, and the Chrome-download
fallback control flow (without performing the real ~150 MB download).
"""

import os

import plotly.graph_objects as go
import pytest

from qpcr import export_utils


@pytest.fixture(autouse=True)
def _clean_browser_env(monkeypatch):
    """Isolate BROWSER_PATH, the download cache, and the render memo per test.

    Every test here renders the same ``_fig()``, so without clearing the render
    cache one test's successful render would be served to the next (a correct
    cache hit, but it defeats the error-path assertions).
    """
    monkeypatch.delenv("BROWSER_PATH", raising=False)
    monkeypatch.setattr(export_utils, "_downloaded_chrome_path", None)
    export_utils.clear_render_cache()
    yield
    export_utils.clear_render_cache()


def _fig():
    return go.Figure(data=[go.Bar(x=["A", "B", "C"], y=[1, 2, 3])])


# ---- browser resolution -------------------------------------------------

def test_ensure_browser_path_prefers_chrome_over_chromium(monkeypatch):
    """When both exist, real Chrome must win (more reliable headless)."""
    def fake_which(name):
        return {
            "chromium": "/usr/bin/chromium",
            "google-chrome": "/usr/bin/google-chrome",
        }.get(name)

    monkeypatch.setattr(export_utils.shutil, "which", fake_which)
    resolved = export_utils._ensure_browser_path()
    assert resolved == "/usr/bin/google-chrome"
    assert os.environ["BROWSER_PATH"] == "/usr/bin/google-chrome"


def test_ensure_browser_path_drops_stale_override(monkeypatch):
    """A BROWSER_PATH pointing at a non-existent binary is dropped."""
    monkeypatch.setenv("BROWSER_PATH", "/does/not/exist/chromium")
    monkeypatch.setattr(export_utils.shutil, "which", lambda name: None)
    monkeypatch.setattr(export_utils.os.path, "exists", lambda p: False)
    resolved = export_utils._ensure_browser_path()
    assert resolved is None
    assert "BROWSER_PATH" not in os.environ


def test_ensure_browser_path_keeps_valid_override(monkeypatch):
    """An existing BROWSER_PATH override is respected as-is."""
    monkeypatch.setenv("BROWSER_PATH", "/custom/chrome")
    monkeypatch.setattr(export_utils.os.path, "exists",
                        lambda p: p == "/custom/chrome")
    assert export_utils._ensure_browser_path() == "/custom/chrome"


# ---- error classification ----------------------------------------------

def test_looks_like_browser_error():
    assert export_utils._looks_like_browser_error(
        RuntimeError("The browser seemed to close immediately after starting"))
    assert export_utils._looks_like_browser_error(
        ValueError("Chrome not found at /usr/bin/chromium"))
    assert not export_utils._looks_like_browser_error(
        ValueError("some unrelated numerical error"))


# ---- fallback control flow ----------------------------------------------

def test_fallback_downloads_chrome_then_retries(monkeypatch):
    """First render raises a browser error -> download fallback -> retry succeeds."""
    calls = {"to_image": 0, "download": 0}

    def fake_to_image(**kwargs):
        calls["to_image"] += 1
        if calls["to_image"] == 1:
            raise RuntimeError("The browser seemed to close immediately after starting")
        return b"PNGDATA"

    def fake_download():
        calls["download"] += 1
        return "/downloaded/chrome"

    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)
    monkeypatch.setattr(export_utils, "_download_fallback_chrome", fake_download)

    fig = _fig()
    monkeypatch.setattr(fig, "to_image", fake_to_image)

    out = export_utils.export_figure_to_bytes(fig, fmt="png")
    assert out == b"PNGDATA"
    assert calls["to_image"] == 2   # failed once, retried once
    assert calls["download"] == 1
    assert os.environ["BROWSER_PATH"] == "/downloaded/chrome"


def test_non_browser_error_not_retried(monkeypatch):
    """Unrelated errors must NOT trigger the expensive download fallback."""
    calls = {"download": 0}
    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)
    monkeypatch.setattr(export_utils, "_download_fallback_chrome",
                        lambda: calls.__setitem__("download", calls["download"] + 1))

    fig = _fig()
    monkeypatch.setattr(fig, "to_image",
                        lambda **k: (_ for _ in ()).throw(ValueError("bad data shape")))

    with pytest.raises(ValueError, match="bad data shape"):
        export_utils.export_figure_to_bytes(fig, fmt="png")
    assert calls["download"] == 0


def test_both_attempts_fail_raises_actionable_error(monkeypatch):
    """When download also fails, a clear RuntimeError is raised."""
    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)

    def boom():
        raise RuntimeError("network down")

    monkeypatch.setattr(export_utils, "_download_fallback_chrome", boom)
    fig = _fig()
    monkeypatch.setattr(fig, "to_image",
                        lambda **k: (_ for _ in ()).throw(RuntimeError("chromium closed")))

    with pytest.raises(RuntimeError, match="headless Chrome/Chromium"):
        export_utils.export_figure_to_bytes(fig, fmt="png")


# ---- render memo cache ---------------------------------------------------

def test_identical_figure_renders_once(monkeypatch):
    """A repeat export of the same figure must not relaunch Chrome."""
    calls = {"n": 0}

    def fake_to_image(**kwargs):
        calls["n"] += 1
        return b"PNG-" + str(calls["n"]).encode()

    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)
    fig = _fig()
    monkeypatch.setattr(fig, "to_image", fake_to_image)

    first = export_utils.export_figure_to_bytes(fig, fmt="png", scale=2)
    second = export_utils.export_figure_to_bytes(fig, fmt="png", scale=2)
    assert first == second == b"PNG-1"
    assert calls["n"] == 1  # second call served from memo


def test_render_params_are_part_of_the_key(monkeypatch):
    """Same figure at a different scale/size must render fresh, not reuse bytes."""
    calls = {"n": 0}

    def fake_to_image(**kwargs):
        calls["n"] += 1
        return b"PNG-" + str(calls["n"]).encode()

    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)
    fig = _fig()
    monkeypatch.setattr(fig, "to_image", fake_to_image)

    export_utils.export_figure_to_bytes(fig, fmt="png", scale=2)
    export_utils.export_figure_to_bytes(fig, fmt="png", scale=3)      # scale differs
    export_utils.export_figure_to_bytes(fig, fmt="svg", scale=2)      # format differs
    export_utils.export_figure_to_bytes(fig, fmt="png", scale=2, width=800)
    assert calls["n"] == 4


def test_changed_figure_content_invalidates(monkeypatch):
    """Editing the figure changes its JSON, so the memo must miss."""
    calls = {"n": 0}

    def fake_to_image(**kwargs):
        calls["n"] += 1
        return b"PNG-" + str(calls["n"]).encode()

    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)
    fig = _fig()
    monkeypatch.setattr(fig, "to_image", fake_to_image)

    export_utils.export_figure_to_bytes(fig, fmt="png")
    fig.update_layout(title="now different")
    export_utils.export_figure_to_bytes(fig, fmt="png")
    assert calls["n"] == 2


def test_failed_render_is_not_cached(monkeypatch):
    """An error must never be memoized as a result."""
    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)
    fig = _fig()
    monkeypatch.setattr(fig, "to_image",
                        lambda **k: (_ for _ in ()).throw(ValueError("bad data shape")))
    with pytest.raises(ValueError):
        export_utils.export_figure_to_bytes(fig, fmt="png")

    monkeypatch.setattr(fig, "to_image", lambda **k: b"OK")
    assert export_utils.export_figure_to_bytes(fig, fmt="png") == b"OK"


def test_cache_evicts_beyond_entry_cap(monkeypatch):
    """The memo is bounded — it must not grow without limit."""
    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)
    monkeypatch.setattr(export_utils, "_RENDER_CACHE_MAX_ENTRIES", 4)
    for i in range(10):
        fig = _fig()
        fig.update_layout(title=f"fig-{i}")
        monkeypatch.setattr(fig, "to_image", lambda **k: b"X" * 10)
        export_utils.export_figure_to_bytes(fig, fmt="png")
    assert len(export_utils._render_cache) <= 4


def test_unserializable_figure_still_renders(monkeypatch):
    """If to_json fails the render must proceed uncached, not blow up."""
    monkeypatch.setattr(export_utils, "_ensure_browser_path", lambda: None)

    class _Weird:
        def to_json(self):
            raise TypeError("not serializable")

        def to_image(self, **kwargs):
            return b"RAW"

    assert export_utils.export_figure_to_bytes(_Weird(), fmt="png") == b"RAW"


# ---- build_zip ----------------------------------------------------------

def test_build_zip_mixed_content_types():
    import io
    import zipfile

    z = export_utils.build_zip({
        "report.txt": "hello",           # str
        "data.bin": b"\x00\x01\x02",     # bytes
        "stream.dat": io.BytesIO(b"xyz"),  # file-like
        "dropped.pptx": None,             # skipped
    })
    zf = zipfile.ZipFile(io.BytesIO(z))
    assert zf.namelist() == ["report.txt", "data.bin", "stream.dat"]
    assert zf.read("report.txt") == b"hello"
    assert zf.read("stream.dat") == b"xyz"


def test_build_zip_empty():
    import io
    import zipfile

    zf = zipfile.ZipFile(io.BytesIO(export_utils.build_zip({})))
    assert zf.namelist() == []


# ---- real end-to-end render (skips if no working browser present) --------

def test_real_export_succeeds_if_browser_available():
    """Smoke test: on a machine with a working Chrome, PNG bytes are produced."""
    if export_utils._find_system_browser() is None:
        pytest.skip("no system browser available")
    try:
        out = export_utils.export_figure_to_bytes(_fig(), fmt="png", scale=1,
                                                  width=400, height=300)
    except RuntimeError as e:
        pytest.skip(f"browser present but could not launch headless: {e}")
    assert isinstance(out, bytes) and out[:8] == b"\x89PNG\r\n\x1a\n"
