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
