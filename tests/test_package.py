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
        assert len(EFFICACY_CONFIG) == 10
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
            QPCRParser, QualityControl, AnalysisEngine,
            GraphGenerator, ReportGenerator, PPTGenerator,
        )
        assert hasattr(QPCRParser, 'parse')
        assert hasattr(QualityControl, 'grubbs_test')
        assert hasattr(AnalysisEngine, 'calculate_ddct')
        assert hasattr(GraphGenerator, 'create_gene_graph')
        assert hasattr(ReportGenerator, 'create_presentation')
        assert hasattr(PPTGenerator, 'generate_presentation')

    def test_import_export(self):
        from qpcr import export_to_excel
        assert callable(export_to_excel)


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
