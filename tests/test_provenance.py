"""Tests for the analysis provenance/reproducibility record (Phase 2b)."""

import json
from importlib import import_module


def _spec():
    return import_module("streamlit qpcr analysis v1")


def test_build_provenance_captures_run_parameters(mock_streamlit):
    prov = _spec().build_provenance(
        efficacy="탄력",
        hk_gene="GAPDH",
        ref_condition="Non-treated",
        cmp_conditions=["Non-treated", "TGFb", None],  # None must be dropped
        ttest_type="welch",
        excluded_wells={("COL1A1", "D2"): {"B6", "B5"}, ("GAPDH", "D1"): {"A3"}},
        excluded_samples={"BadSample"},
        n_genes=3,
        n_samples=6,
        timestamp="2026-07-09 12:00:00",
    )
    assert prov["reference_gene"] == "GAPDH"
    assert prov["reference_condition"] == "Non-treated"
    assert prov["comparison_conditions"] == ["Non-treated", "TGFb"]  # None dropped
    # The test description now carries the replicate unit, which MIQE requires
    # be declared because it sets the degrees of freedom.
    assert prov["statistical_test"].startswith("Welch t-test")
    assert "technical replicate wells" in prov["statistical_test"]
    assert prov["replicate_unit"].startswith("technical wells")
    # The record must describe what the numbers are: per-well fold change,
    # simple mean, spread over the same wells (hand-sheet method, 2026-10-28).
    assert "per-well fold change" in prov["method"]
    assert "arithmetic mean" in prov["method"]
    assert "per-well fold change" in prov["statistical_test"]
    assert "Livak" not in prov["method"]
    assert prov["excluded_wells_count"] == 3
    assert prov["excluded_samples"] == ["BadSample"]
    # Excluded wells are flattened, sorted, and carry gene+sample+well.
    assert prov["excluded_wells"][0] == {"gene": "COL1A1", "sample": "D2", "well": "B5"}
    assert all({"gene", "sample", "well"} == set(e) for e in prov["excluded_wells"])
    # Must be JSON-serialisable (it's offered as a download).
    json.dumps(prov)


def test_build_provenance_student_and_empty(mock_streamlit):
    prov = _spec().build_provenance(
        efficacy=None, hk_gene=None, ref_condition=None, cmp_conditions=None,
        ttest_type="student", excluded_wells={}, excluded_samples=set(),
        n_genes=0, n_samples=0, timestamp="t",
    )
    assert prov["statistical_test"].startswith("Student t-test")
    assert prov["comparison_conditions"] == []
    assert prov["excluded_wells_count"] == 0
    assert prov["excluded_samples"] == []


def test_format_provenance_text_lists_excluded_wells(mock_streamlit):
    spec = _spec()
    prov = spec.build_provenance(
        efficacy="Anti-Aging", hk_gene="GAPDH", ref_condition="Ctrl",
        cmp_conditions=["Ctrl"], ttest_type="welch",
        excluded_wells={("COL1A1", "S1"): {"A1"}}, excluded_samples=set(),
        n_genes=1, n_samples=2, timestamp="2026-07-09 12:00:00",
    )
    text = spec.format_provenance_text(prov)
    assert "PROVENANCE" in text
    assert "GAPDH" in text
    assert "Welch t-test" in text
    assert "COL1A1 / S1 / well A1" in text


def _clear_version_cache(spec):
    spec.app_version_string.cache_clear()


def test_app_version_string_format(mock_streamlit):
    import re
    spec = _spec()
    _clear_version_cache(spec)
    assert re.match(
        r"^qPCR Analysis Suite 2026\.10( \(commit [0-9a-g-]+(-dirty)?\))?$",
        spec.app_version_string(),
    )


def test_app_version_string_without_git(mock_streamlit, monkeypatch):
    import subprocess
    spec = _spec()

    def boom(*a, **k):
        raise OSError("no git")

    monkeypatch.setattr(subprocess, "run", boom)
    _clear_version_cache(spec)
    try:
        assert spec.app_version_string() == "qPCR Analysis Suite 2026.10"
    finally:
        _clear_version_cache(spec)


def test_build_provenance_stays_pure(mock_streamlit, monkeypatch):
    import subprocess
    spec = _spec()

    def boom(*a, **k):
        raise AssertionError("build_provenance must not spawn a process")

    monkeypatch.setattr(subprocess, "run", boom)
    prov = spec.build_provenance(
        efficacy=None, hk_gene=None, ref_condition=None, cmp_conditions=None,
        ttest_type="student", excluded_wells={}, excluded_samples=set(),
        n_genes=0, n_samples=0, timestamp="t",
    )
    assert prov["software"] == "qPCR Analysis Suite 2026.10"


def test_provenance_gene_count_rules(mock_streamlit):
    import pandas as pd
    spec = _spec()
    count = spec.provenance_gene_count
    raw = pd.DataFrame({"Target": ["actin", "PPARy", "SREBP1a", "actin"]})

    assert count({"PPARy": 1, "SREBP1a": 2}, raw, "ACTIN") == 2  # results win
    assert count({}, raw, "ACTIN") == 2          # reference gene, any case, excluded
    assert count(None, raw, "ACTIN") == 2
    assert count({}, raw, None) == 3
    assert count({}, raw, "") == 3
    assert count({}, raw, "ABSENT") == 3
    assert count({}, None, "ACTIN") == 0
    assert count({}, pd.DataFrame(), "ACTIN") == 0
    assert count({}, pd.DataFrame({"Sample": ["s"]}), "ACTIN") == 0  # no Target column
