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


SECRET = "ZZSECRETSTUDY"  # stands in for a client study code


def _fake_local_root(tmp_path, csv_names, csv_bytes):
    root = tmp_path / "td"
    (root / "golden" / "recipes").mkdir(parents=True)
    for name in csv_names:
        (root / name).write_bytes(csv_bytes)
    return root


def _pytest_local(env_extra):
    node = "tests/test_golden_outputs.py::test_local_goldens_match"
    return subprocess.run(
        [sys.executable, "-m", "pytest", node, "-m", "local_golden", "-q", "-p", "no:cacheprovider"],
        capture_output=True, text=True, timeout=300, cwd=str(REPO), env={**os.environ, **env_extra},
    )


def test_local_failures_never_name_the_recipe_or_csv(tmp_path):
    csv = f"{SECRET} 1.csv"
    root = _fake_local_root(tmp_path, [csv], b"not,a,qpcr,export\n1,2,3,4\n")
    rec = replace(_recipe("synthetic_b"), name=SECRET, csv_files=[csv])
    (root / "golden" / "recipes" / f"{SECRET}.json").write_text(rec.to_json(), encoding="utf-8")
    (root / "golden" / "BASELINE").write_text("feedface\n")
    (root / "golden" / "out" / "feedface").mkdir(parents=True)
    env = {"QPCR_TEST_DATA": str(root), "QPCR_REQUIRE_LOCAL_GOLDENS": "1"}

    runs = {
        "compare": _capture_cli("--local", env_extra=env),
        "init": _capture_cli("--init-local", env_extra=env),
        "pytest": _pytest_local(env),
    }
    for label, proc in runs.items():
        shown = proc.stdout + proc.stderr
        assert proc.returncode != 0, (label, shown)
        assert SECRET not in shown, label
        # A runner without Chrome fails the pytest run on the skip reason instead.
        chrome_missing = label == "pytest" and "Chrome not found" in shown
        assert recipes.opaque_id(SECRET) in shown or chrome_missing, (label, shown)
    # The full story is kept, but only under the test-data root.
    details = "".join(p.read_text(encoding="utf-8") for p in (root / "golden" / "diffs").iterdir())
    assert SECRET in details and "could NOT be loaded" in details


def test_plate_sets_group_only_plate_indices():
    stem = capture.plate_set_stem
    assert stem("990101 KER10 1") == stem("990101 KER10 2") == "990101 KER10"
    assert stem("990101 KER11 1") == "990101 KER11"
    assert stem("X-ALPHA1") == stem("X-ALPHA2_data") == "X-ALPHA"
    assert stem("Y BETA1-1") == stem("Y BETA2_data-1") == "Y BETA"
    assert stem("Z GAMMA1") == stem("Z GAMMA 2") == "Z GAMMA"
    # A multi-digit tail is a study code, not a plate index.
    assert stem("STUDY12") != stem("STUDY13")


def test_init_local_prints_plate_sets_as_ids_and_counts(tmp_path):
    names = ["990101 KER10 1.csv", "990101 KER10 2.csv", "990101 KER11 1.csv"]
    root = _fake_local_root(tmp_path, names, recipes.synthetic_csv("synthetic_b.csv"))
    proc = _capture_cli("--init-local", env_extra={"QPCR_TEST_DATA": str(root)})
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "KER" not in proc.stdout
    assert f"{recipes.opaque_id('990101 KER10')}: 2 CSV(s)" in proc.stdout
    assert f"{recipes.opaque_id('990101 KER11')}: 1 CSV(s)" in proc.stdout
    written = recipes.Recipe.from_json(
        (root / "golden" / "recipes" / "990101 KER10.json").read_text(encoding="utf-8")
    )
    assert written.csv_files == names[:2]


def test_chrome_missing_is_a_named_skip_reason(monkeypatch, tmp_path):
    import qpcr.export_utils

    root = _fake_local_root(tmp_path, [], b"")
    rec = replace(_recipe("synthetic_b"), name="local_x", csv_files=["x 1.csv"])
    (root / "golden" / "recipes" / "local_x.json").write_text(rec.to_json(), encoding="utf-8")
    (root / "golden" / "BASELINE").write_text("feedface\n")
    (root / "golden" / "out" / "feedface").mkdir(parents=True)
    monkeypatch.setenv("QPCR_TEST_DATA", str(root))
    assert capture.local_skip_reason() is None or "Chrome" in capture.local_skip_reason()
    monkeypatch.setattr(qpcr.export_utils, "_find_system_browser", lambda: None)
    assert "Chrome" in capture.local_skip_reason()


def test_update_local_refuses_a_dirty_tree(monkeypatch, tmp_path):
    root = _fake_local_root(tmp_path, [], b"")
    monkeypatch.setenv("QPCR_TEST_DATA", str(root))
    monkeypatch.setattr(capture, "is_dirty", lambda repo=None: True)
    assert capture.main(["--update", "--local"]) == 2
    assert not (root / "golden" / "BASELINE").exists()


def test_is_dirty_tracks_uncommitted_changes(tmp_path):
    def git(*args):
        subprocess.run(["git", "-C", str(tmp_path), *args], check=True, capture_output=True)

    git("init", "-q")
    git("config", "user.email", "t@example.invalid")
    git("config", "user.name", "t")
    (tmp_path / "f.txt").write_text("a")
    assert capture.is_dirty(tmp_path)
    git("add", "f.txt")
    git("commit", "-qm", "init")
    assert not capture.is_dirty(tmp_path)
    (tmp_path / "f.txt").write_text("b")
    assert capture.is_dirty(tmp_path)


def test_gene_csv_is_the_served_payload(synthetic_outputs):
    text = synthetic_outputs["synthetic_a"]["gene_csv"]["COL1A1"]
    head = text.splitlines()[:3]
    assert head[0] == "file: COL1A1_results_20260101.csv"
    assert head[1] == "bom: utf-8"
    assert head[2].startswith("Condition,Group,Fold_Change,")


def test_recipe_ttest_default_is_the_apps():
    from qpcr.analysis import DEFAULT_TTEST_TYPE

    r = recipes.Recipe(name="x", csv_files=[], ref_sample="a", cmp_samples=["a", None, None])
    assert r.ttest_type == DEFAULT_TTEST_TYPE


@pytest.mark.local_golden
def test_local_goldens_match():
    reason = capture.local_skip_reason()
    if reason:
        if os.environ.get("QPCR_REQUIRE_LOCAL_GOLDENS") == "1":
            pytest.fail(reason)
        pytest.skip(reason)
    # Opaque ids and counts only: client values stay in <root>/golden/diffs/.
    failed = [msg for r in recipes.load_local_recipes() if (msg := capture.check_local(r))]
    assert not failed, "; ".join(failed) + " (details under <test-data root>/golden/diffs/)"
