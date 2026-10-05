"""Regression tests for the 2026-10 audit (Claude subagents + Codex).

Each test pins one defect that was reproduced against the code before the fix.
"""
import io

import numpy as np
import pandas as pd

from qpcr.analysis import AnalysisEngine
from qpcr.graph import GraphGenerator
from qpcr.parser import QPCRParser
from qpcr.quality_control import QualityControl


def _csv(rows, header="Well Position,Sample Name,Target Name,CT"):
    return io.BytesIO((header + "\n" + "\n".join(rows) + "\n").encode("utf-8"))


class TestGrubbsKeyedPerGroup:
    def test_outlier_flag_does_not_spill_to_other_samples_sharing_a_well(self):
        data = pd.DataFrame(
            {
                "Well": ["A1", "A2", "A3", "A1", "A2", "A3"],
                "Sample": ["S1"] * 3 + ["S2"] * 3,
                "Target": ["G"] * 6,
                "CT": [20.0, 20.1, 28.0, 20.0, 20.0, 20.0],
            }
        )
        qc = QualityControl.detect_outliers(data)
        flagged = qc[qc["Issues"].str.contains("Grubbs")]
        assert set(zip(flagged["Sample"], flagged["Well"])) == {("S1", "A3")}


class TestThresholdsForwarded:
    def test_summary_stats_health_uses_custom_ct_high(self):
        wells = [f"A{i}" for i in range(1, 4)]
        data = pd.DataFrame(
            {"Well": wells, "Sample": ["S"] * 3, "Target": ["G"] * 3,
             "CT": [36.0, 36.1, 35.9]}
        )
        stats = QualityControl.get_qc_summary_stats(
            data, thresholds={"ct_high": 40}
        )
        assert stats["high_ct_count"] == 0
        assert stats["healthy_triplicates"] == stats["total_triplicates"]


class TestHousekeepingCaseConsistency:
    def test_lowercase_hk_name_still_normalises(self):
        rows = []
        for cond, ct in (("Ctl", 25.0), ("Trt", 24.0)):
            rows += [
                {"Well": f"A{i}", "Sample": cond, "Condition": cond,
                 "Target": "GAPDH", "CT": 18.0}
                for i in range(1, 4)
            ]
            rows += [
                {"Well": f"B{i}", "Sample": cond, "Condition": cond,
                 "Target": "G1", "CT": ct}
                for i in range(1, 4)
            ]
        data = pd.DataFrame(rows)
        mapping = {c: {"condition": c, "group": "g", "include": True}
                   for c in ("Ctl", "Trt")}
        res = AnalysisEngine.calculate_ddct(
            data, "gapdh", "Ctl", set(), set(), mapping
        )
        assert not res.empty
        trt = res[res["Condition"] == "Trt"].iloc[0]
        assert np.isclose(trt["Relative_Expression"], 2.0)


class TestParserRobustness:
    def test_sample_named_NA_is_kept(self):
        f = _csv(["A1,NA,G,20.0", "A2,NA,G,20.1", "A3,S2,G,21.0"])
        out = QPCRParser.parse(f)
        assert out is not None
        assert "NA" in set(out["Sample"])

    def test_undetermined_ct_is_still_dropped(self):
        f = _csv(["A1,S1,G,Undetermined", "A2,S1,G,20.1"])
        out = QPCRParser.parse(f)
        assert len(out) == 1

    def test_infinite_ct_is_dropped_not_turned_into_a_fold_change(self):
        f = _csv(["A1,S1,G,inf", "A2,S1,G,20.1", "A3,S1,G,20.2"])
        out = QPCRParser.parse(f)
        assert np.isfinite(out["CT"]).all()
        assert len(out) == 2

    def test_lowercase_headers_parse(self):
        f = _csv(["A1,S1,G,20.0"], header="well,sample,target,ct")
        out = QPCRParser.parse(f)
        assert out is not None and len(out) == 1


class TestRefLineWithoutErrorBars:
    def test_reference_line_with_error_bars_off_builds(self, processed_gene_data):
        fig = GraphGenerator.create_gene_graph(
            data=processed_gene_data,
            gene="COL1A1",
            settings={"show_error": False},
            ref_line_value=1.0,
        )
        assert len(fig.data) >= 1


class TestExcelWritesNamesVerbatim:
    def test_export_does_not_turn_uploaded_names_into_formulas(self, mock_streamlit):
        import openpyxl
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        processed = {
            "G1": pd.DataFrame({
                "Target": ["G1", "G1"],
                "Condition": ["Non-treated", "=1+1"],
                "Group": ["Negative Control", "Treatment"],
                "Fold_Change": [1.0, 2.0],
                "Relative_Expression": [1.0, 2.0],
                "SEM": [0.05, 0.1],
            })
        }
        raw = pd.DataFrame(
            [{"Well": "A1", "Sample": s_, "Target": "G1", "CT": 20.0}
             for s_ in ("Non-treated", "=1+1")]
        )
        mapping = {"Non-treated": {"condition": "Non-treated", "group": "Negative Control"},
                   "=1+1": {"condition": "=1+1", "group": "Treatment"}}
        xlsx = spec.export_to_excel(
            raw, processed,
            {"Housekeeping_Gene": "GAPDH", "Efficacy_Type": "Anti-Aging"}, mapping,
        )
        wb = openpyxl.load_workbook(io.BytesIO(xlsx.getvalue() if hasattr(xlsx, "getvalue") else xlsx))
        formulas = [
            (ws.title, c.coordinate, c.value)
            for ws in wb.worksheets for row in ws.iter_rows() for c in row
            if c.data_type == "f"
        ]
        assert formulas == []


class TestCodexReviewFollowups:
    def test_well_header_matched_case_insensitively_when_not_first(self):
        f = _csv(["S1,G,20.0,A1", "S1,G,20.1,A2"], header="sample,target,ct,well")
        out = QPCRParser.parse(f)
        assert list(out["Well"]) == ["A1", "A2"]

    def test_hk_exclusion_applies_to_case_variant_rows(self):
        rows = []
        for cond in ("Ctl", "Trt"):
            for i, ct in enumerate((18.0, 18.0, 30.0), start=1):
                rows.append({"Well": f"A{i}", "Sample": cond, "Condition": cond,
                             "Target": "gapdh", "CT": ct})
            for i in range(1, 4):
                rows.append({"Well": f"B{i}", "Sample": cond, "Condition": cond,
                             "Target": "G1", "CT": 25.0})
        data = pd.DataFrame(rows)
        mapping = {c: {"condition": c, "group": "g", "include": True}
                   for c in ("Ctl", "Trt")}
        excl = {("gapdh", c): {"A3"} for c in ("Ctl", "Trt")}
        res = AnalysisEngine.calculate_ddct(
            data, "GAPDH", "Ctl", excl, set(), mapping
        )
        trt = res[res["Condition"] == "Trt"].iloc[0]
        assert np.isclose(trt["Relative_Expression"], 1.0)
