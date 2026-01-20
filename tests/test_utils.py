import pytest


class TestNaturalSortKey:
    def test_natural_sort_key_basic_numbers(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        samples = ['Sample10', 'Sample2', 'Sample1', 'Sample20']
        sorted_samples = sorted(samples, key=natural_sort_key)
        
        assert sorted_samples == ['Sample1', 'Sample2', 'Sample10', 'Sample20']

    def test_natural_sort_key_mixed_content(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        samples = ['A10B2', 'A2B1', 'A1B10']
        sorted_samples = sorted(samples, key=natural_sort_key)
        
        assert sorted_samples == ['A1B10', 'A2B1', 'A10B2']

    def test_natural_sort_key_no_numbers(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        samples = ['Charlie', 'Alpha', 'Bravo']
        sorted_samples = sorted(samples, key=natural_sort_key)
        
        assert sorted_samples == ['Alpha', 'Bravo', 'Charlie']

    def test_natural_sort_key_case_insensitive(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        samples = ['sample1', 'SAMPLE2', 'Sample3']
        sorted_samples = sorted(samples, key=natural_sort_key)
        
        assert sorted_samples == ['sample1', 'SAMPLE2', 'Sample3']

    def test_natural_sort_key_empty_string(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        result = natural_sort_key('')
        assert result == ['']

    def test_natural_sort_key_none_handling(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        natural_sort_key = spec.natural_sort_key
        
        result = natural_sort_key(None)
        assert isinstance(result, list)
