import io
import pandas as pd
import pytest


class TestQPCRParserDetectFormat:
    def test_detect_format1_with_well_position(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        df = pd.DataFrame([
            ['Block Type', '', '', ''],
            ['Well Position', 'Sample Name', 'Target Name', 'CT'],
            ['A1', 'Sample1', 'GAPDH', '18.5']
        ])
        
        fmt, idx = QPCRParser.detect_format(df)
        assert fmt == 'format1'
        assert idx == 1

    def test_detect_format1_with_well_and_sample_name(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        df = pd.DataFrame([
            ['Well', 'Well Position', 'Sample Name', 'Target Name', 'CT'],
            ['A1', 'A1', 'Sample1', 'GAPDH', '18.5']
        ])
        
        fmt, idx = QPCRParser.detect_format(df)
        assert fmt == 'format1'
        assert idx == 0

    def test_detect_unknown_format(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        df = pd.DataFrame([
            ['Random', 'Data', 'Here'],
            ['1', '2', '3']
        ])
        
        fmt, idx = QPCRParser.detect_format(df)
        assert fmt == 'unknown'
        assert idx == 0


class TestQPCRParserParseFormat1:
    def test_parse_format1_standard_columns(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        df = pd.DataFrame([
            ['Well Position', 'Well', 'Sample Name', 'Target Name', 'CT'],
            ['A1', 'A1', 'Non-treated', 'GAPDH', '18.52'],
            ['A2', 'A2', 'Non-treated', 'GAPDH', '18.48'],
            ['A3', 'A3', 'Non-treated', 'COL1A1', '25.12']
        ])
        
        result = QPCRParser.parse_format1(df, 0)
        
        assert result is not None
        assert len(result) == 3
        assert list(result.columns) == ['Well', 'Sample', 'Target', 'CT']
        assert result['CT'].dtype == float

    def test_parse_format1_missing_ct_column_returns_none(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        df = pd.DataFrame([
            ['Well Position', 'Sample Name', 'Target Name', 'Value'],
            ['A1', 'Sample1', 'GAPDH', '18.5']
        ])
        
        result = QPCRParser.parse_format1(df, 0)
        assert result is None

    def test_parse_format1_insufficient_columns_returns_none(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        df = pd.DataFrame([
            ['Col1', 'Col2', 'Col3'],
            ['A', 'B', 'C']
        ])
        
        result = QPCRParser.parse_format1(df, 0)
        assert result is None

    def test_parse_format1_filters_nan_ct_values(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        df = pd.DataFrame([
            ['Well Position', 'Sample Name', 'Target Name', 'CT'],
            ['A1', 'Sample1', 'GAPDH', '18.5'],
            ['A2', 'Sample1', 'GAPDH', 'Undetermined'],
            ['A3', 'Sample1', 'GAPDH', '18.7']
        ])
        
        result = QPCRParser.parse_format1(df, 0)
        
        assert len(result) == 2
        assert 'Undetermined' not in result['CT'].values


class TestQPCRParserParse:
    def test_parse_valid_format1_csv(self, mock_streamlit, format1_csv_content):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        file = io.StringIO(format1_csv_content)
        result = QPCRParser.parse(file)
        
        assert result is not None
        assert 'Well' in result.columns
        assert 'Sample' in result.columns
        assert 'Target' in result.columns
        assert 'CT' in result.columns
        assert len(result) == 12

    def test_parse_returns_none_on_unknown_format(self, mock_streamlit):
        from importlib import import_module
        spec = import_module('streamlit qpcr analysis v1')
        QPCRParser = spec.QPCRParser
        
        csv_content = """RandomCol1,RandomCol2,RandomCol3
A1,Sample1,18.5
"""
        file = io.StringIO(csv_content)
        result = QPCRParser.parse(file)
        
        assert result is None
