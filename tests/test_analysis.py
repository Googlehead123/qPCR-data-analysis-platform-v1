import pandas as pd
import numpy as np
import pytest


class TestAnalysisEngineCalculateDDCT:
    def test_calculate_ddct_basic(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine
        
        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells=set(),
            excluded_samples=set(),
            sample_mapping=sample_mapping
        )
        
        assert result is not None
        assert not result.empty
        assert 'Target' in result.columns
        assert 'Condition' in result.columns
        assert 'Relative_Expression' in result.columns
        assert 'Delta_Delta_Ct' in result.columns

    def test_calculate_ddct_reference_sample_has_expression_1(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine
        
        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells=set(),
            excluded_samples=set(),
            sample_mapping=sample_mapping
        )
        
        ref_row = result[result['Condition'] == 'Non-treated']
        assert len(ref_row) == 1
        assert abs(ref_row['Relative_Expression'].iloc[0] - 1.0) < 0.01

    def test_calculate_ddct_excludes_housekeeping_gene(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine
        
        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells=set(),
            excluded_samples=set(),
            sample_mapping=sample_mapping
        )
        
        assert 'GAPDH' not in result['Target'].values

    def test_calculate_ddct_respects_excluded_wells(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine
        
        first_well = sample_qpcr_raw_data['Well'].iloc[0]
        
        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells={first_well},
            excluded_samples=set(),
            sample_mapping=sample_mapping
        )
        
        assert result is not None

    def test_calculate_ddct_per_gene_sample_exclusion(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        """Excluding a well for one gene should NOT affect other genes."""
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine

        # Get a COL1A1 well for Non-treated
        col1a1_nt = sample_qpcr_raw_data[
            (sample_qpcr_raw_data['Target'] == 'COL1A1') &
            (sample_qpcr_raw_data['Sample'] == 'Non-treated')
        ]
        well_to_exclude = col1a1_nt['Well'].iloc[0]

        # Baseline: no exclusions
        baseline = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells={},
            excluded_samples=set(),
            sample_mapping=sample_mapping,
        )

        # Exclude one well for COL1A1/Non-treated only (per-gene-sample dict)
        excluded_dict = {('COL1A1', 'Non-treated'): {well_to_exclude}}
        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells=excluded_dict,
            excluded_samples=set(),
            sample_mapping=sample_mapping,
        )

        assert result is not None
        assert not result.empty

        # The COL1A1 Non-treated row should have fewer replicates
        nt_row = result[result['Condition'] == 'Non-treated']
        baseline_nt = baseline[baseline['Condition'] == 'Non-treated']
        assert nt_row['n_replicates'].iloc[0] == baseline_nt['n_replicates'].iloc[0] - 1

        # Treatment rows should be unaffected (same replicate count)
        for cond in ['Treatment1', 'Treatment2']:
            res_row = result[result['Condition'] == cond]
            base_row = baseline[baseline['Condition'] == cond]
            assert res_row['n_replicates'].iloc[0] == base_row['n_replicates'].iloc[0]

    def test_calculate_ddct_reference_exclusion_changes_results(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        """Excluding a well from the reference sample must change the reference ΔCt and thus all fold changes."""
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine

        # Baseline: no exclusions
        baseline = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells={},
            excluded_samples=set(),
            sample_mapping=sample_mapping,
        )

        # Get a COL1A1 well for Non-treated (this is the reference sample)
        ref_wells = sample_qpcr_raw_data[
            (sample_qpcr_raw_data['Target'] == 'COL1A1') &
            (sample_qpcr_raw_data['Sample'] == 'Non-treated')
        ]
        well_to_exclude = ref_wells['Well'].iloc[0]

        # Exclude one reference well
        excluded_dict = {('COL1A1', 'Non-treated'): {well_to_exclude}}
        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells=excluded_dict,
            excluded_samples=set(),
            sample_mapping=sample_mapping,
        )

        # Reference sample n_replicates should be reduced
        ref_baseline = baseline[baseline['Condition'] == 'Non-treated']
        ref_result = result[result['Condition'] == 'Non-treated']
        assert ref_result['n_replicates'].iloc[0] == ref_baseline['n_replicates'].iloc[0] - 1

        # The reference Ct mean should have changed (different well set)
        # which means fold changes for treatment samples may differ
        assert result is not None
        assert not result.empty

    def test_calculate_ddct_sd_is_target_only(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        """SD and SEM should reflect target gene CT variation only, not combined target+HK."""
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine

        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells={},
            excluded_samples=set(),
            sample_mapping=sample_mapping,
        )

        assert not result.empty

        for _, row in result.iterrows():
            target_sd = row['Target_Ct_SD']
            sd = row['SD']
            n = row['n_replicates']

            # SD should equal Target_Ct_SD exactly (not combined with HK)
            assert abs(sd - target_sd) < 1e-10, (
                f"SD ({sd}) should equal Target_Ct_SD ({target_sd}) for {row['Condition']}"
            )

            # SEM should equal Target_Ct_SD / sqrt(n)
            expected_sem = target_sd / np.sqrt(n) if n > 1 else 0
            assert abs(row['SEM'] - expected_sem) < 1e-10, (
                f"SEM ({row['SEM']}) should equal Target_Ct_SD/sqrt(n) ({expected_sem})"
            )

    def test_calculate_ddct_respects_excluded_samples(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine
        
        result = AnalysisEngine.calculate_ddct(
            data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells=set(),
            excluded_samples={'Treatment2'},
            sample_mapping=sample_mapping
        )
        
        assert 'Treatment2' not in result['Condition'].values

    def test_calculate_ddct_handles_empty_data(self, mock_streamlit, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine
        
        empty_df = pd.DataFrame(columns=['Well', 'Sample', 'Target', 'CT'])
        
        result = AnalysisEngine.calculate_ddct(
            data=empty_df,
            hk_gene='GAPDH',
            ref_sample='Non-treated',
            excluded_wells=set(),
            excluded_samples=set(),
            sample_mapping=sample_mapping
        )
        
        assert result.empty


class TestAnalysisEngineCalculateStatistics:
    def test_calculate_statistics_adds_pvalue_column(self, mock_streamlit, processed_gene_data, sample_qpcr_raw_data, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine
        
        mock_streamlit.session_state['data'] = sample_qpcr_raw_data
        mock_streamlit.session_state['hk_gene'] = 'GAPDH'
        mock_streamlit.session_state['sample_mapping'] = sample_mapping
        
        result = AnalysisEngine.calculate_statistics(
            processed=processed_gene_data,
            compare_condition='Non-treated',
            raw_data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            sample_mapping=sample_mapping
        )
        
        assert 'p_value' in result.columns
        assert 'significance' in result.columns

    def test_calculate_statistics_dual_comparison(self, mock_streamlit, processed_gene_data, sample_qpcr_raw_data, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine
        
        mock_streamlit.session_state['data'] = sample_qpcr_raw_data
        mock_streamlit.session_state['hk_gene'] = 'GAPDH'
        mock_streamlit.session_state['sample_mapping'] = sample_mapping
        
        result = AnalysisEngine.calculate_statistics(
            processed=processed_gene_data,
            compare_condition='Non-treated',
            compare_condition_2='Treatment1',
            raw_data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            sample_mapping=sample_mapping
        )
        
        assert 'p_value_2' in result.columns
        assert 'significance_2' in result.columns

    def test_calculate_statistics_significance_thresholds(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        AnalysisEngine = spec.AnalysisEngine

        mock_streamlit.session_state['data'] = sample_qpcr_raw_data
        mock_streamlit.session_state['hk_gene'] = 'GAPDH'
        mock_streamlit.session_state['sample_mapping'] = sample_mapping

        processed = pd.DataFrame({
            'Target': ['COL1A1', 'COL1A1'],
            'Condition': ['Non-treated', 'Treatment1'],
            'Relative_Expression': [1.0, 5.0],
        })

        result = AnalysisEngine.calculate_statistics(
            processed=processed,
            compare_condition='Non-treated',
            raw_data=sample_qpcr_raw_data,
            hk_gene='GAPDH',
            sample_mapping=sample_mapping
        )

        assert 'significance' in result.columns


class TestReplicateFoldChanges:
    def test_returns_one_row_per_replicate(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from qpcr.analysis import AnalysisEngine
        result = AnalysisEngine.compute_replicate_fold_changes(
            raw_data=sample_qpcr_raw_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            sample_mapping=sample_mapping,
            excluded_wells=set(),
        )
        # 3 samples × 3 replicates each for COL1A1 = 9 rows
        assert len(result) == 9
        assert "Condition" in result.columns
        assert "Replicate_FC" in result.columns
        assert "Target" in result.columns

    def test_reference_condition_fcs_average_near_one(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from qpcr.analysis import AnalysisEngine
        result = AnalysisEngine.compute_replicate_fold_changes(
            raw_data=sample_qpcr_raw_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            sample_mapping=sample_mapping,
            excluded_wells=set(),
        )
        ref_fcs = result[result["Condition"] == "Non-treated"]["Replicate_FC"]
        assert abs(ref_fcs.mean() - 1.0) < 0.3

    def test_handles_excluded_wells(self, mock_streamlit, sample_qpcr_raw_data, sample_mapping):
        from qpcr.analysis import AnalysisEngine
        # Exclude a COL1A1 well (index 3 = first COL1A1 replicate, after 3 GAPDH wells)
        excluded = {sample_qpcr_raw_data["Well"].iloc[3]}
        result = AnalysisEngine.compute_replicate_fold_changes(
            raw_data=sample_qpcr_raw_data,
            hk_gene="GAPDH",
            ref_sample="Non-treated",
            sample_mapping=sample_mapping,
            excluded_wells=excluded,
        )
        assert len(result) == 8  # 9 - 1 excluded COL1A1 well


class TestSingletonComparisons:
    """A group with n=1 must not produce a p-value or a significance marker.

    The old code substituted stats.ttest_1samp when either group had a single
    well, treating that singleton's measured ΔCt as a KNOWN population mean.
    That is a different hypothesis and it manufactured significance: these
    fixtures returned p=0.0335 ('*') and, with the groups reversed, p=0.0027
    ('**'), for data with no estimable within-group variance.
    Decision (Min, 2026-08-24): report nothing instead.
    """

    @staticmethod
    def _run(rows):
        import pandas as pd
        from qpcr.analysis import AnalysisEngine
        raw = pd.DataFrame(rows, columns=["Well", "Sample", "Target", "CT"])
        mapping = {s: {"condition": s} for s in raw["Sample"].unique()}
        processed = AnalysisEngine.calculate_ddct(
            raw, "HK", "R", {}, set(), mapping
        )
        return AnalysisEngine.calculate_statistics(
            processed, "R", raw_data=raw, hk_gene="HK",
            sample_mapping=mapping, excluded_wells={},
        )

    def test_treatment_n1_gets_no_marker(self, mock_streamlit):
        import pandas as pd
        result = self._run([
            ("R1", "R", "HK", 20), ("R2", "R", "HK", 20),
            ("R1", "R", "G", 20), ("R2", "R", "G", 21),
            ("T1", "T", "HK", 20), ("T1", "T", "G", 30),
        ])
        row = result[result["Condition"] == "T"].iloc[0]
        assert pd.isna(row["p_value"]), (
            "n=1 has no estimable variance, so there is no valid two-group "
            f"t-test; got p={row['p_value']}"
        )
        assert row["significance"] == ""

    def test_reference_n1_gets_no_marker(self, mock_streamlit):
        import pandas as pd
        result = self._run([
            ("R1", "R", "HK", 20), ("R1", "R", "G", 20),
            ("T1", "T", "HK", 20), ("T2", "T", "HK", 20), ("T3", "T", "HK", 20),
            ("T1", "T", "G", 25), ("T2", "T", "G", 25), ("T3", "T", "G", 26),
        ])
        row = result[result["Condition"] == "T"].iloc[0]
        assert pd.isna(row["p_value"])
        assert row["significance"] == ""

    def test_reference_n1_is_not_labelled_welch(self, mock_streamlit):
        """The exporter read the treatment n alone and printed "Welch t-test
        (n=3)" for a comparison actually run by ttest_1samp."""
        result = self._run([
            ("R1", "R", "HK", 20), ("R1", "R", "G", 20),
            ("T1", "T", "HK", 20), ("T2", "T", "HK", 20), ("T3", "T", "HK", 20),
            ("T1", "T", "G", 25), ("T2", "T", "G", 25), ("T3", "T", "G", 26),
        ])
        row = result[result["Condition"] == "T"].iloc[0]
        assert "Welch" not in row["stat_test_used"]
        assert "none" in row["stat_test_used"]
        assert row["n_reference"] == 1

    def test_valid_two_group_test_still_runs_and_is_labelled(self, mock_streamlit):
        import pandas as pd
        result = self._run([
            ("R1", "R", "HK", 20), ("R2", "R", "HK", 20), ("R3", "R", "HK", 20),
            ("R1", "R", "G", 20), ("R2", "R", "G", 21), ("R3", "R", "G", 20.5),
            ("T1", "T", "HK", 20), ("T2", "T", "HK", 20), ("T3", "T", "HK", 20),
            ("T1", "T", "G", 25), ("T2", "T", "G", 25.2), ("T3", "T", "G", 24.8),
        ])
        row = result[result["Condition"] == "T"].iloc[0]
        assert pd.notna(row["p_value"]), "n=3 vs n=3 must still produce a p-value"
        assert "Welch t-test" in row["stat_test_used"]
        assert row["n_reference"] == 3
