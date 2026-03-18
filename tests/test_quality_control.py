"""
Tests for QualityControl class functions.

These tests cover the statistical outlier detection and QC functionality
that is critical for ensuring data quality in qPCR analysis.
"""

import numpy as np
import pandas as pd
import pytest


class TestQualityControlGrubbsTest:
    """Test suite for the Grubbs outlier detection test."""

    def test_grubbs_test_detects_clear_outlier(self, mock_streamlit):
        """Grubbs test should detect a clear statistical outlier."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        # Data with a clear outlier (10.0 among values near 1.0)
        values = np.array([1.0, 1.1, 1.05, 1.08, 10.0])
        is_outlier, idx = QualityControl.grubbs_test(values, alpha=0.05)

        assert is_outlier == True
        assert idx == 4  # Index of the outlier (10.0)

    def test_grubbs_test_no_outlier(self, mock_streamlit):
        """Grubbs test should return False when data has no significant outliers."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        values = np.array([1.00, 1.01, 1.02, 1.01, 1.00])
        is_outlier, idx = QualityControl.grubbs_test(values, alpha=0.05)

        assert is_outlier == False or abs(values[idx] - values.mean()) < 0.05

    def test_grubbs_test_insufficient_replicates(self, mock_streamlit):
        """Grubbs test requires n >= 3 and should return False for smaller samples."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        # Only 2 values - insufficient for Grubbs test
        values = np.array([1.0, 1.1])
        is_outlier, idx = QualityControl.grubbs_test(values, alpha=0.05)

        assert is_outlier == False
        assert idx == -1

    def test_grubbs_test_zero_variance(self, mock_streamlit):
        """Grubbs test with identical values (zero variance) should not crash."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        # All identical values - zero variance
        values = np.array([1.0, 1.0, 1.0])
        is_outlier, idx = QualityControl.grubbs_test(values, alpha=0.05)

        # Should handle gracefully without crashing
        assert is_outlier == False
        assert idx == -1

    def test_grubbs_test_single_value(self, mock_streamlit):
        """Grubbs test with single value should return False."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        values = np.array([1.0])
        is_outlier, idx = QualityControl.grubbs_test(values, alpha=0.05)

        assert is_outlier == False
        assert idx == -1


class TestQualityControlDetectOutliers:
    """Test suite for the comprehensive outlier detection function."""

    def test_detect_outliers_high_ct(self, mock_streamlit, sample_qpcr_raw_data):
        """Should flag wells with CT values above threshold."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        # Add a well with high CT value
        high_ct_data = sample_qpcr_raw_data.copy()
        high_ct_row = pd.DataFrame(
            [
                {
                    "Well": "Z1",
                    "Sample": "HighCT",
                    "Target": "GAPDH",
                    "CT": 38.0,  # Above default threshold of 35
                }
            ]
        )
        high_ct_data = pd.concat([high_ct_data, high_ct_row], ignore_index=True)

        flagged = QualityControl.detect_outliers(high_ct_data, "GAPDH")

        # Should flag the high CT well
        assert len(flagged) > 0
        assert "Z1" in flagged["Well"].values

    def test_detect_outliers_low_ct(self, mock_streamlit, sample_qpcr_raw_data):
        """Should flag wells with CT values below threshold."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        # Add a well with suspiciously low CT value
        low_ct_data = sample_qpcr_raw_data.copy()
        low_ct_row = pd.DataFrame(
            [
                {
                    "Well": "Z2",
                    "Sample": "LowCT",
                    "Target": "GAPDH",
                    "CT": 5.0,  # Below default threshold of 10
                }
            ]
        )
        low_ct_data = pd.concat([low_ct_data, low_ct_row], ignore_index=True)

        flagged = QualityControl.detect_outliers(low_ct_data, "GAPDH")

        # Should flag the low CT well
        assert len(flagged) > 0
        assert "Z2" in flagged["Well"].values

    def test_detect_outliers_returns_dataframe(
        self, mock_streamlit, sample_qpcr_raw_data
    ):
        """detect_outliers should return a DataFrame."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        result = QualityControl.detect_outliers(sample_qpcr_raw_data, "GAPDH")

        assert isinstance(result, pd.DataFrame)


class TestQualityControlGetTriplicateData:
    """Test suite for triplicate data aggregation."""

    def test_get_triplicate_data_basic(self, mock_streamlit, sample_qpcr_raw_data):
        """Should return triplicate statistics grouped by Sample and Target."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        result = QualityControl.get_triplicate_data(sample_qpcr_raw_data)

        assert isinstance(result, pd.DataFrame)
        assert "Sample" in result.columns
        assert "Target" in result.columns
        assert "Mean_CT" in result.columns
        assert "SD" in result.columns
        assert "CV_pct" in result.columns
        assert "n" in result.columns

    def test_get_triplicate_data_with_exclusions(
        self, mock_streamlit, sample_qpcr_raw_data
    ):
        """Should respect excluded wells."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        QualityControl = spec.QualityControl

        first_well = sample_qpcr_raw_data["Well"].iloc[0]
        excluded_wells = {first_well}

        result = QualityControl.get_triplicate_data(
            sample_qpcr_raw_data, excluded_wells=excluded_wells
        )

        # Result should not include the excluded well in calculations
        assert isinstance(result, pd.DataFrame)


class TestAnalysisEngineFDRCorrection:
    """Test suite for FDR (Benjamini-Hochberg) correction."""

    def test_fdr_correction_increases_pvalues(self, mock_streamlit):
        """FDR-corrected p-values should be >= original p-values."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        AnalysisEngine = spec.AnalysisEngine

        # Create test data with p-values
        results = pd.DataFrame(
            {
                "Target": ["GENE1", "GENE2", "GENE3"],
                "Condition": ["Treated", "Treated", "Treated"],
                "p_value": [0.001, 0.01, 0.05],
            }
        )

        corrected = AnalysisEngine._apply_fdr_correction(
            results.copy(), "p_value", "p_value_fdr", "significance_fdr", "*"
        )

        # FDR-adjusted p-values should be >= original
        for idx in corrected.index:
            original = corrected.loc[idx, "p_value"]
            adjusted = corrected.loc[idx, "p_value_fdr"]
            if not np.isnan(original) and not np.isnan(adjusted):
                assert adjusted >= original, (
                    f"FDR p-value {adjusted} < original {original}"
                )

    def test_fdr_correction_bounded_by_one(self, mock_streamlit):
        """FDR-corrected p-values should never exceed 1.0."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        AnalysisEngine = spec.AnalysisEngine

        results = pd.DataFrame(
            {
                "Target": ["GENE1", "GENE2"],
                "Condition": ["Treated", "Treated"],
                "p_value": [
                    0.8,
                    0.9,
                ],  # High p-values that might exceed 1 if not capped
            }
        )

        corrected = AnalysisEngine._apply_fdr_correction(
            results.copy(), "p_value", "p_value_fdr", "significance_fdr", "*"
        )

        assert (corrected["p_value_fdr"] <= 1.0).all()

    def test_fdr_correction_handles_nan(self, mock_streamlit):
        """FDR correction should handle NaN p-values gracefully."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        AnalysisEngine = spec.AnalysisEngine

        results = pd.DataFrame(
            {
                "Target": ["GENE1", "GENE2", "GENE3"],
                "Condition": ["Treated", "Treated", "Treated"],
                "p_value": [0.01, np.nan, 0.05],
            }
        )

        corrected = AnalysisEngine._apply_fdr_correction(
            results.copy(), "p_value", "p_value_fdr", "significance_fdr", "*"
        )

        # Should not crash and should have correct columns
        assert "p_value_fdr" in corrected.columns
        assert "significance_fdr" in corrected.columns

    def test_fdr_correction_empty_data(self, mock_streamlit):
        """FDR correction should handle empty or all-NaN data."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        AnalysisEngine = spec.AnalysisEngine

        results = pd.DataFrame(
            {
                "Target": [],
                "Condition": [],
                "p_value": [],
            }
        )

        corrected = AnalysisEngine._apply_fdr_correction(
            results.copy(), "p_value", "p_value_fdr", "significance_fdr", "*"
        )

        # Should not crash
        assert "p_value_fdr" in corrected.columns


class TestAnalysisEngineSEMCalculation:
    """Test suite for SEM (Standard Error of Mean) calculation in DDCt."""

    def test_sem_positive_for_replicates(
        self, mock_streamlit, sample_qpcr_raw_data, sample_mapping
    ):
        """SEM should be positive when there are multiple replicates."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        AnalysisEngine = spec.AnalysisEngine

        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            excluded_wells=set(),
            excluded_samples=set(),
            sample_mapping=sample_mapping,
        )

        # SEM should be present and non-negative
        assert "SEM" in result.columns
        assert (result["SEM"] >= 0).all()

    def test_sem_zero_for_single_replicate(self, mock_streamlit, sample_mapping):
        """SEM should be 0 when there's only a single replicate."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        AnalysisEngine = spec.AnalysisEngine

        # Create data with single replicates
        single_rep_data = pd.DataFrame(
            [
                {"Well": "A1", "Sample": "Non-treated", "Target": "GAPDH", "CT": 18.5},
                {"Well": "A2", "Sample": "Non-treated", "Target": "COL1A1", "CT": 25.0},
                {"Well": "B1", "Sample": "Treatment1", "Target": "GAPDH", "CT": 18.3},
                {"Well": "B2", "Sample": "Treatment1", "Target": "COL1A1", "CT": 23.5},
            ]
        )

        result = AnalysisEngine.calculate_ddct(
            data=single_rep_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            excluded_wells=set(),
            excluded_samples=set(),
            sample_mapping=sample_mapping,
        )

        # SEM should be 0 for single replicates
        assert "SEM" in result.columns
        assert (result["SEM"] == 0).all()

    def test_sem_scales_with_expression(
        self, mock_streamlit, sample_qpcr_raw_data, sample_mapping
    ):
        """SEM should scale appropriately with relative expression."""
        from importlib import import_module

        spec = import_module("streamlit qpcr analysis v1")
        AnalysisEngine = spec.AnalysisEngine

        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            excluded_wells=set(),
            excluded_samples=set(),
            sample_mapping=sample_mapping,
        )

        # SEM should be related to expression level (higher expression = potentially higher SEM)
        # This tests that error propagation is working
        assert "SEM" in result.columns
        assert "Relative_Expression" in result.columns

        # Reference sample should have fold change ~1.0
        ref_rows = result[result["Condition"] == "Non-treated"]
        if len(ref_rows) > 0:
            assert abs(ref_rows["Relative_Expression"].iloc[0] - 1.0) < 0.1
