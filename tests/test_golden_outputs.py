"""Golden-output gate: the app's real pipeline and every export, frozen.

Synthetic recipes run in CI against goldens committed under
``tests/golden/synthetic``. Local recipes run over the client plates in
``test_data/`` with real Kaleido rendering; they live, and are compared, only
under the test-data root (marker ``local_golden``).
"""

import io
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from tests.golden import capture, driver, recipes, stub_render

REPO = Path(__file__).resolve().parents[1]
SYNTH_DIR = REPO / "tests" / "golden" / "synthetic"


@pytest.fixture(scope="module")
def synthetic_outputs():
    return {r.name: driver.produce_outputs(r, real_render=False) for r in recipes.synthetic_recipes()}


def _recipe(name):
    return next(r for r in recipes.synthetic_recipes() if r.name == name)


def test_recipe_round_trips_through_json():
    for r in recipes.synthetic_recipes():
        assert recipes.Recipe.from_json(r.to_json()) == r


def test_synthetic_goldens_match(synthetic_outputs):
    for r in recipes.synthetic_recipes():
        golden = capture.load_golden(SYNTH_DIR / r.name)
        assert golden, f"no committed goldens for {r.name}"
        assert capture.compare(synthetic_outputs[r.name], golden, allow=[]) == []


def _bump(columns):
    def mutate(app):
        df = app.st.session_state.processed_data["COL1A1"]
        row = df.index[df["Condition"] == "Active A"][0]
        df.loc[row, columns] *= 1.01

    return mutate


def test_harness_detects_changed_bar():
    r = _recipe("synthetic_a")
    golden = capture.load_golden(SYNTH_DIR / r.name)
    # The chart draws Relative_Expression (Fold_Change is only graph.py's
    # fallback); the analysis writes both from the same value, so a moved bar
    # moves both.
    bar = capture.compare(
        driver.produce_outputs(r, real_render=False, mutate=_bump(["Fold_Change", "Relative_Expression"])),
        golden, allow=[],
    )
    for key in ("processed", "figures", "pptx"):
        assert any(v.startswith(f"{key}/") for v in bar), key
    # Same-size slide picture, different pixels: the stub colour carries the change.
    assert any(v.startswith("pptx/ppt/media/") for v in bar)

    # Fold_Change alone never reaches the chart, but the tables and workbook show it.
    fc_only = capture.compare(
        driver.produce_outputs(r, real_render=False, mutate=_bump(["Fold_Change"])),
        golden, allow=[],
    )
    for key in ("processed", "xlsx", "gene_csv"):
        assert any(v.startswith(f"{key}/") for v in fc_only), key


def test_stub_render_sizes_like_kaleido_and_keeps_format():
    import plotly.graph_objects as go
    from PIL import Image

    def size(data):
        return Image.open(io.BytesIO(data)).size

    fig = go.Figure(go.Bar(y=[1, 2])).update_layout(width=300, height=200)
    assert size(stub_render._fig_to_image(fig, format="png", scale=2)) == (600, 400)
    assert size(stub_render._fig_to_image(fig, format="png", scale=3, width=100, height=50)) == (300, 150)
    assert size(stub_render._fig_to_image(go.Figure(go.Bar(y=[1])), scale=1)) == (64, 64)
    svg = stub_render._fig_to_image(fig, format="svg", scale=1)
    pdf = stub_render._fig_to_image(fig, format="pdf", scale=1)
    assert not svg.startswith(b"\x89PNG") and svg != pdf


def test_compare_allowlist_permits_only_matching_lines():
    golden = {"xlsx": {"sheet:QC": "A1\tHigh CT Wells (>35)\nB1\t3\n"}}
    actual = {"xlsx": {"sheet:QC": "A1\tHigh CT Wells (>38)\nB1\t4\n"}}
    allow = [{"part": "xlsx/sheet:QC", "line": r"^A1\tHigh CT Wells \(>\d+\)$"}]
    assert capture.compare(actual, golden, allow) == ["xlsx/sheet:QC: -B1\t3", "xlsx/sheet:QC: +B1\t4"]
    # The line regex alone is not enough: the part glob must match too.
    wrong_part = [{"part": "pptx/*", "line": ".*"}]
    assert len(capture.compare(actual, golden, wrong_part)) == 4
    # A part that appears or disappears is never allowlistable.
    assert capture.compare({"xlsx": {}}, golden, [{"part": "*", "line": ".*"}]) == [
        "xlsx/sheet:QC: part missing"
    ]


def test_exclusion_hits_only_its_gene():
    r = replace(_recipe("synthetic_a"), auto_qc_threshold=None)
    seen = {}

    def grab(app):
        state = app.st.session_state
        seen["excl"] = {k: set(v) for k, v in state.excluded_wells.items()}
        seen["n"] = {
            g: dict(zip(df["Condition"], df["n_replicates"]))
            for g, df in state.processed_data.items()
        }

    driver.produce_outputs(r, real_render=False, mutate=grab)
    assert seen["excl"][("MMP1", "S3")] == {"B1"}
    assert "B1" not in seen["excl"].get(("COL1A1", "S3"), set())
    # ...and the analysis agrees: only MMP1 lost a well in Active B (= S3).
    assert seen["n"]["MMP1"]["Active B"] == 2
    assert seen["n"]["COL1A1"]["Active B"] == 3


def _capture_cli(*args, env_extra=None):
    env = {**os.environ, **(env_extra or {})}
    return subprocess.run(
        [sys.executable, "-m", "tests.golden.capture", *args],
        capture_output=True, text=True, timeout=600, cwd=str(REPO), env=env,
    )


def test_stability_across_processes():
    for seed in ("1", "2"):
        proc = _capture_cli("--recipe", "synthetic_a", env_extra={"PYTHONHASHSEED": seed})
        assert proc.returncode == 0, f"PYTHONHASHSEED={seed}\n{proc.stdout}\n{proc.stderr}"


def test_no_error_messages(synthetic_outputs):
    for name, out in synthetic_outputs.items():
        lines = [ln for text in out["messages"].values() for ln in text.splitlines()]
        assert lines, f"{name}: no messages recorded"
        assert not [ln for ln in lines if ln.startswith("error:")], name


def test_only_synthetic_goldens_are_tracked():
    names = {r.name for r in recipes.synthetic_recipes()}
    entries = sorted(SYNTH_DIR.iterdir())
    assert entries and all(p.is_dir() and p.name in names for p in entries), entries
    tracked = subprocess.run(
        ["git", "ls-files", str(SYNTH_DIR.relative_to(REPO))],
        capture_output=True, text=True, check=True, cwd=str(REPO),
    ).stdout.split()
    assert {Path(t).parts[3] for t in tracked} <= names

    local = replace(_recipe("synthetic_b"), name="local_leak", csv_files=["plate 1.csv"])
    dest = SYNTH_DIR / "local_leak"
    with pytest.raises(ValueError):
        capture.write_golden(local, {"messages": {"log": "x\n"}}, dest)
    assert not dest.exists()


def test_local_goldens_skip_reason_and_require_flag(monkeypatch, tmp_path):
    missing = tmp_path / "nowhere"
    monkeypatch.setenv("QPCR_TEST_DATA", str(missing))
    reason = capture.local_skip_reason()
    assert reason and "test_data" in reason

    node = "tests/test_golden_outputs.py::test_local_goldens_match"
    base = [sys.executable, "-m", "pytest", node, "-m", "local_golden", "-q", "-p", "no:cacheprovider"]
    env = {**os.environ, "QPCR_TEST_DATA": str(missing)}
    env.pop("QPCR_REQUIRE_LOCAL_GOLDENS", None)
    skipped = subprocess.run(base + ["-rs"], capture_output=True, text=True, timeout=300,
                             cwd=str(REPO), env=env)
    assert skipped.returncode == 0 and "skipped" in skipped.stdout, skipped.stdout
    required = subprocess.run(base, capture_output=True, text=True, timeout=300, cwd=str(REPO),
                              env={**env, "QPCR_REQUIRE_LOCAL_GOLDENS": "1"})
    assert required.returncode == 1 and "test_data" in required.stdout, required.stdout


@pytest.mark.local_golden
def test_local_goldens_match():
    reason = capture.local_skip_reason()
    if reason:
        if os.environ.get("QPCR_REQUIRE_LOCAL_GOLDENS") == "1":
            pytest.fail(reason)
        pytest.skip(reason)
    out_dir = capture.local_baseline_dir()
    failed = []
    for i, r in enumerate(recipes.load_local_recipes(), 1):
        violations = capture.compare(
            driver.produce_outputs(r, real_render=True), capture.load_golden(out_dir / r.name), allow=[]
        )
        if violations:
            capture.write_local_diff(r, violations)
            failed.append(f"local recipe #{i}: {len(violations)} changed line(s)")
    # Counts only: client values stay in <root>/golden/diffs/, never in test output.
    assert not failed, "; ".join(failed) + " (diffs under <test-data root>/golden/diffs/)"
