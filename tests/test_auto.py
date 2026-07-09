"""Tests for the deterministic Auto-Analyze engine (qpcr.auto)."""

import numpy as np
import pandas as pd
import pytest

from qpcr.auto import screen_data, recommend_test, interpret_results, interpret_gene


# ---------------- screening ----------------

def _raw(n_rep=3, hk="GAPDH", with_nondetect=False):
    rows = []
    for s in ("Ctrl", "Tx"):
        for i in range(n_rep):
            rows.append({"Well": f"{s}g{i}", "Sample": s, "Target": hk, "CT": 20.0 + i * 0.1})
            ct = 24.0 + i * 0.1
            if with_nondetect and s == "Tx" and i == 0:
                ct = np.nan
            rows.append({"Well": f"{s}c{i}", "Sample": s, "Target": "COL1A1", "CT": ct})
    return pd.DataFrame(rows)


def test_screen_clean_data_ok():
    rep = screen_data(_raw(), "GAPDH")
    assert rep["ok"] is True
    assert rep["summary"]["hk_present"] is True
    assert rep["summary"]["n_targets"] == 2


def test_screen_missing_hk_is_error():
    rep = screen_data(_raw(), "ACTB")  # not present
    assert rep["ok"] is False
    assert any(i["level"] == "error" and "ACTB" in i["message"] for i in rep["issues"])


def test_screen_missing_column_is_error():
    bad = pd.DataFrame({"Sample": ["a"], "Target": ["x"], "CT": [1.0]})  # no Well
    rep = screen_data(bad, "x")
    assert rep["ok"] is False
    assert any("Well" in i["message"] for i in rep["issues"])


def test_screen_flags_nondetect_and_thin_replicates():
    rep = screen_data(_raw(n_rep=2, with_nondetect=True), "GAPDH")
    msgs = " ".join(i["message"] for i in rep["issues"])
    assert "non-detect" in msgs
    assert rep["summary"]["n_nondetect"] == 1
    # Tx/COL1A1 now has a single detected replicate -> single-replicate warning
    assert "single replicate" in msgs


# ---------------- stats advisor ----------------

def test_recommend_needs_two_groups():
    r = recommend_test([[1, 2, 3]])
    assert r["test"] is None and r["n_groups"] == 1


def test_recommend_nonnormal_two_groups_mann_whitney():
    rng = np.random.default_rng(0)
    normal = rng.normal(10, 1, 12)
    skewed = np.concatenate([np.ones(11), [500.0]])  # gross outlier -> non-normal
    r = recommend_test([skewed, normal])
    assert r["n_groups"] == 2
    assert r["test"] == "Mann-Whitney U test"
    assert r["normal"] is False


def test_recommend_normal_equal_variance_student():
    rng = np.random.default_rng(1)
    a = rng.normal(10, 1.0, 40)
    b = rng.normal(12, 1.0, 40)
    r = recommend_test([a, b])
    assert r["test"] == "Student's t-test"
    assert r["equal_var"] is True


def test_recommend_normal_unequal_variance_welch():
    rng = np.random.default_rng(2)
    a = rng.normal(10, 1.0, 40)
    b = rng.normal(12, 6.0, 40)  # much larger spread -> Levene significant
    r = recommend_test([a, b])
    assert r["test"] == "Welch's t-test"
    assert r["equal_var"] is False


def test_recommend_three_groups_omnibus_with_followup():
    rng = np.random.default_rng(3)
    groups = [rng.normal(m, 1.0, 30) for m in (10, 11, 12)]
    r = recommend_test(groups)
    assert r["n_groups"] == 3
    assert r["test"] in ("One-way ANOVA", "Kruskal-Wallis test")
    assert r["followup"] is not None


def test_recommend_small_n_flags_caution():
    r = recommend_test([[5.0, 5.1], [4.0, 4.2]])  # n=2 each
    assert r["normal"] is None
    assert "caution" in r["reason"].lower() or r["test"] is not None


# ---------------- interpretation ----------------

def _gene_df(fold_tx, sig="**", p=0.003):
    return pd.DataFrame({
        "Condition": ["Ctrl", "Tx"],
        "Fold_Change": [1.0, fold_tx],
        "p_value": [np.nan, p],
        "significance": ["", sig],
        "n_replicates": [3, 3],
    })


def test_interpret_matches_expected_up():
    res = interpret_gene("COL1A1", _gene_df(2.4), expected="up", ref_condition="Ctrl")
    tx = [r for r in res["rows"] if r["condition"] == "Tx"][0]
    assert tx["direction"] == "up"
    assert tx["verdict"] == "as expected"
    assert "COL1A1" in res["narrative_en"]
    assert "증가" in res["narrative_ko"]


def test_interpret_opposite_to_expected():
    res = interpret_gene("MMP1", _gene_df(3.0), expected="down", ref_condition="Ctrl")
    tx = [r for r in res["rows"] if r["condition"] == "Tx"][0]
    assert tx["verdict"] == "opposite to expected"
    assert "OPPOSITE" in res["narrative_en"]


def test_interpret_nonsignificant_is_ns():
    res = interpret_gene("COL1A1", _gene_df(2.4, sig="", p=0.4), expected="up", ref_condition="Ctrl")
    tx = [r for r in res["rows"] if r["condition"] == "Tx"][0]
    assert tx["verdict"] == "n.s."


def test_interpret_results_maps_expected_direction_by_gene():
    processed = {"COL1A1": _gene_df(2.0), "MMP1": _gene_df(0.4)}
    out = interpret_results(processed, expected_direction={"COL1A1": "up", "MMP1": "down"},
                            ref_condition="Ctrl")
    verdicts = {o["gene"]: o["rows"][0]["verdict"] for o in out}
    assert verdicts["COL1A1"] == "as expected"
    assert verdicts["MMP1"] == "as expected"  # 0.4-fold down, expected down
