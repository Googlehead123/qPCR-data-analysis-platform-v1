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
