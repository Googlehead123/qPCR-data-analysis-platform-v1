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
    assert prov["statistical_test"] == "Welch t-test"
    assert prov["method"].startswith("Livak")
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
    assert prov["statistical_test"] == "Student t-test"
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
