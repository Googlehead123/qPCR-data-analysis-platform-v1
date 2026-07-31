"""Tests for the qpcr modular package.

Verifies that all public classes and functions are importable from the
qpcr package and produce identical results to the monolithic source.
"""

import numpy as np
import pandas as pd


class TestPackageImports:
    """Verify all expected symbols are importable from the qpcr package."""

    def test_import_constants(self):
        from qpcr import (
            EFFICACY_CONFIG, AnalysisConstants, DEFAULT_GROUP_COLORS,
            COSMAX_RED, COSMAX_BLACK, COSMAX_WHITE, COSMAX_LAB_WHITE,
            COSMAX_FROST_GREY, COSMAX_CREAM, PLOTLY_FONT_FAMILY,
            CM_TO_PX, CM_TO_EMU,
        )
        assert isinstance(EFFICACY_CONFIG, dict)
        assert len(EFFICACY_CONFIG) == 21
        # Spot-check the updated catalog (효능평가항목_update.xlsx).
        assert "MMP1" in EFFICACY_CONFIG["광노화"]["genes"]
        assert EFFICACY_CONFIG["광노화"]["expected_direction"]["MMP1"] == "down"
        assert EFFICACY_CONFIG["멜라닌 생성"]["expected_direction"]["MITF"] == "down"
        assert EFFICACY_CONFIG["장벽"]["controls"]["positive"] == "Retinoic acid"
        assert "expected_direction" not in EFFICACY_CONFIG["립 색상"]  # ambiguous → omitted
        assert AnalysisConstants.CT_UNDETERMINED_THRESHOLD == 40.0
        assert CM_TO_PX > 37

    def test_import_utils(self):
        from qpcr import (
            natural_sort_key, get_well_exclusion_key, get_grid_cell_key,
            get_selected_cell, set_selected_cell, clear_selected_cell, is_cell_selected,
        )
        assert callable(natural_sort_key)
        assert callable(get_grid_cell_key)

    def test_import_classes(self):
        from qpcr import (
            QPCRParser, QualityControl, AnalysisEngine, GraphGenerator,
        )
        assert hasattr(QPCRParser, 'parse')
        assert hasattr(QualityControl, 'grubbs_test')
        assert hasattr(AnalysisEngine, 'calculate_ddct')
        assert hasattr(GraphGenerator, 'create_gene_graph')

    def test_report_and_export_live_in_the_monolith(self):
        """Reporting/export are owned by the app module, not the package.

        `qpcr/report.py` and `qpcr/export.py` were deleted: the app never
        imported them, so any fix applied there silently never shipped. This
        pins the surface at its single real location.
        """
        from importlib import import_module
        spec = import_module("streamlit qpcr analysis v1")

        assert hasattr(spec.ReportGenerator, 'create_presentation')
        assert hasattr(spec.PPTGenerator, 'generate_presentation')
        assert callable(spec.export_to_excel)

        # And they must NOT come back as package modules.
        import pytest
        for dead in ("qpcr.report", "qpcr.export"):
            with pytest.raises(ImportError):
                import_module(dead)


class TestControlLabelsCarryDoses:
    """Control labels must gain doses without losing their matching key.

    `controls.positive` / `controls.negative` are dual-purpose: they are printed
    on the PPT ("Positive control:" / "Inducer:") AND compared against real
    condition names from the uploaded data by the mapping tab's `_suggest_group`
    and the summary's `_match` (exact, then substring). A dose may therefore be
    *appended* — "Arbutin" -> "Arbutin 100 ppm" still matches a condition named
    "Arbutin" on the substring pass — but the substance name itself must never be
    rewritten to match the xlsx ledger, or the key the data is matched on moves.
    """

    # The bare substance names the app has always matched on. Values here are
    # prefixes, not the full labels: a dose may follow, nothing may replace.
    HISTORICAL_NAMES = {
        "탄력": ("Non-treated", "TGFβ"),
        "광노화": ("UVB only", "TGFβ"),
        "보습/수분": ("Non-treated", "Retinoic acid"),
        "장벽": ("Non-treated", "Retinoic acid"),
        "속보습": ("Non-treated", None),
        "멜라닌 생성": ("α-MSH only", "Arbutin"),
        "진정": ("Poly(I:C)+IL-4", "Dexamethasone"),
        "가려움 개선": ("Poly(I:C)+IL-4", "Dexamethasone"),
        "냉감": ("Non-treated", "Menthol"),
        "열감": ("Non-treated", None),
        "여드름": ("IGF only", None),
        "과각화": ("SZ95 supernatant", None),
        "활력": ("Non-treated", "EGF"),
        "탈모 개선": ("Non-treated", "Minoxidil"),
        "모공 탄력": ("Non-treated", "Niacinamide"),
        "열 노화": ("Heat (41°C)", "Ascorbic acid"),
        "민감성 기전": ("Non-treated", "Panthenol"),
        "외이도염": ("Bacterial LPS", "Salicylic acid"),
    }

    def test_names_are_prefixes_never_replaced(self):
        from qpcr import EFFICACY_CONFIG

        for efficacy, (neg, pos) in self.HISTORICAL_NAMES.items():
            controls = EFFICACY_CONFIG[efficacy]["controls"]
            for role, expected in (("negative", neg), ("positive", pos)):
                if expected is None:
                    assert role not in controls, (
                        f"{efficacy}.{role} gained a control; add it to this map"
                    )
                    continue
                actual = controls[role]
                assert actual.startswith(expected), (
                    f"{efficacy}.{role} was renamed {expected!r} -> {actual!r}. "
                    "Doses append; names do not change (they are matching keys)."
                )

    def test_a_bare_condition_name_still_matches_a_dosed_label(self):
        """Mirrors the exact-then-substring rule the app's matchers use.

        Not a call into those closures (they are defined inside the tab bodies),
        but the same predicate applied to the catalog: if this fails, a user
        whose sample is named after the substance stops being recognised as the
        positive control and silently lands in the Treatment group.
        """
        from qpcr import EFFICACY_CONFIG

        def matches(label, condition):
            ll, cl = label.strip().lower(), condition.strip().lower()
            return cl == ll or ll in cl or cl in ll

        for efficacy, (neg, pos) in self.HISTORICAL_NAMES.items():
            controls = EFFICACY_CONFIG[efficacy]["controls"]
            for role, bare in (("negative", neg), ("positive", pos)):
                if bare is None:
                    continue
                assert matches(controls[role], bare), (
                    f"{efficacy}.{role}={controls[role]!r} no longer matches a "
                    f"condition named {bare!r}"
                )

    def test_doses_present_where_the_ledger_has_them(self):
        """The point of the change: benchmarks reach the slide with a dose."""
        from qpcr import EFFICACY_CONFIG

        expected = {
            "광노화": "TGFβ 10 ng/ml",
            "보습/수분": "Retinoic acid 1μM",
            "멜라닌 생성": "Arbutin 100 ppm",
            "진정": "Dexamethasone 1μM",
            "가려움 개선": "Dexamethasone 1μM",
            "냉감": "Menthol 100 ppm",
            "모공 탄력": "Niacinamide 0.1%",
            "외이도염": "Salicylic acid 1 mM",
            "민감성 기전": "Panthenol 500 ppm",
            "활력": "EGF 25 ng/ml",
        }
        for efficacy, label in expected.items():
            assert EFFICACY_CONFIG[efficacy]["controls"]["positive"] == label

        # Min's override survives, and it is deliberately dose-less: the ledger
        # row for 장벽 is a different substance entirely (Calcium 1.2mM), so
        # there is no dose on file for Retinoic acid here to borrow.
        assert EFFICACY_CONFIG["장벽"]["controls"]["positive"] == "Retinoic acid"
        # No dose on file for these two.
        assert EFFICACY_CONFIG["탈모 개선"]["controls"]["positive"] == "Minoxidil"
        assert EFFICACY_CONFIG["열 노화"]["controls"]["positive"] == "Ascorbic acid"


class TestPackageParity:
    """Verify package classes produce identical results to monolith."""

    def _get_monolith(self):
        from importlib import import_module
        return import_module("streamlit qpcr analysis v1")

    def test_grubbs_test_parity(self):
        from qpcr import QualityControl as QC_new
        QC_old = self._get_monolith().QualityControl

        values = np.array([18.5, 18.3, 18.6, 25.0])
        assert QC_old.grubbs_test(values) == QC_new.grubbs_test(values)

        values_no_outlier = np.array([18.5, 18.3, 18.6])
        assert QC_old.grubbs_test(values_no_outlier) == QC_new.grubbs_test(values_no_outlier)

    def test_natural_sort_key_parity(self):
        from qpcr import natural_sort_key as nsk_new
        nsk_old = self._get_monolith().natural_sort_key

        names = ["Sample10", "Sample2", "Sample1", "Sample20", "abc", "ABC"]
        assert sorted(names, key=nsk_new) == sorted(names, key=nsk_old)

    def test_ddct_parity(self, sample_qpcr_raw_data, sample_mapping):
        from qpcr import AnalysisEngine as AE_new
        AE_old = self._get_monolith().AnalysisEngine

        old = AE_old.calculate_ddct(
            sample_qpcr_raw_data.copy(), "GAPDH", "Non-treated", set(), set(), sample_mapping
        )
        new = AE_new.calculate_ddct(
            sample_qpcr_raw_data.copy(), "GAPDH", "Non-treated", set(), set(), sample_mapping
        )

        assert len(old) == len(new)
        for col in ["Delta_Ct", "Delta_Delta_Ct", "Relative_Expression"]:
            if col in old.columns and col in new.columns:
                np.testing.assert_allclose(
                    old[col].values, new[col].values, atol=1e-10, err_msg=f"{col} mismatch"
                )

    def test_replicate_stats_parity(self, sample_qpcr_raw_data):
        from qpcr import QualityControl as QC_new
        QC_old = self._get_monolith().QualityControl

        old_stats = QC_old.get_replicate_stats(sample_qpcr_raw_data)
        new_stats = QC_new.get_replicate_stats(sample_qpcr_raw_data)

        assert len(old_stats) == len(new_stats)
        assert list(old_stats.columns) == list(new_stats.columns)

    def test_get_grid_cell_key_parity(self):
        from qpcr import get_grid_cell_key as new_fn
        old_fn = self._get_monolith().get_grid_cell_key

        assert old_fn("GENE1", "Sample_A") == new_fn("GENE1", "Sample_A")
        assert old_fn("", "") == new_fn("", "")

    def test_detect_outliers_parity(self, sample_qpcr_raw_data):
        from qpcr import QualityControl as QC_new
        QC_old = self._get_monolith().QualityControl

        old_qc = QC_old.detect_outliers(sample_qpcr_raw_data, hk_gene="GAPDH")
        new_qc = QC_new.detect_outliers(sample_qpcr_raw_data, hk_gene="GAPDH")

        assert len(old_qc) == len(new_qc)
        assert set(old_qc.columns) == set(new_qc.columns)
