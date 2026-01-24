"""
Tests for ReportGenerator PPT generation functionality.

Tests cover:
- Title slide generation with COSMAX theme
- Gene slide generation with/without method text
- create_presentation() with mock data
- Edge cases (empty graphs, long method text)
"""

import io
from unittest.mock import patch, MagicMock
import pandas as pd
import pytest


# Create a minimal PNG bytes for mocking to_image
MOCK_PNG_BYTES = (
    b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
    b"\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\x0cIDATx\x9cc\xf8\x0f\x00"
    b"\x00\x01\x01\x00\x05\x18\xd8N\x00\x00\x00\x00IEND\xaeB`\x82"
)


class TestReportGeneratorTitleSlide:
    """Tests for title slide generation with COSMAX theme."""

    def test_title_slide_has_logo_placeholder(self, mock_streamlit):
        """Title slide should have a logo placeholder text box."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)

        analysis_params = {
            "Efficacy_Type": "Anti-Aging",
            "Housekeeping_Gene": "GAPDH",
            "Reference_Sample": "Non-treated",
        }

        ReportGenerator._add_title_slide(prs, analysis_params)

        slide = prs.slides[0]
        # Title slide should have shapes (logo, title, subtitle, info, accent bar)
        assert len(slide.shapes) >= 5

        # Check for logo placeholder text
        logo_found = False
        for shape in slide.shapes:
            if hasattr(shape, "text_frame"):
                if "[Add Logo Here]" in shape.text_frame.text:
                    logo_found = True
                    break
        assert logo_found, "Logo placeholder not found on title slide"

    def test_title_slide_has_main_title(self, mock_streamlit):
        """Title slide should have main title 'qPCR Gene Expression Analysis'."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)

        analysis_params = {
            "Efficacy_Type": "Anti-Aging",
            "Housekeeping_Gene": "GAPDH",
            "Reference_Sample": "Non-treated",
        }

        ReportGenerator._add_title_slide(prs, analysis_params)

        slide = prs.slides[0]
        title_found = False
        for shape in slide.shapes:
            if hasattr(shape, "text_frame"):
                if "qPCR Gene Expression Analysis" in shape.text_frame.text:
                    title_found = True
                    break
        assert title_found, "Main title not found on title slide"

    def test_title_slide_has_red_accent_bar(self, mock_streamlit):
        """Title slide should have COSMAX red accent bar at bottom."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches
        from pptx.enum.shapes import MSO_SHAPE_TYPE

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)

        analysis_params = {
            "Efficacy_Type": "Anti-Aging",
            "Housekeeping_Gene": "GAPDH",
            "Reference_Sample": "Non-treated",
        }

        ReportGenerator._add_title_slide(prs, analysis_params)

        slide = prs.slides[0]
        # Check for rectangle shape (accent bar)
        rect_found = False
        for shape in slide.shapes:
            if shape.shape_type == MSO_SHAPE_TYPE.AUTO_SHAPE:
                rect_found = True
                # Verify it has COSMAX red color (0xEA, 0x1D, 0x22)
                if hasattr(shape, "fill") and shape.fill.fore_color:
                    rgb = shape.fill.fore_color.rgb
                    assert rgb == (0xEA, 0x1D, 0x22) or str(rgb) == "EA1D22", (
                        f"Accent bar color should be COSMAX red, got {rgb}"
                    )
                break
        assert rect_found, "Red accent bar not found on title slide"


class TestReportGeneratorGeneSlide:
    """Tests for gene slide generation."""

    def test_gene_slide_with_method_text(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """Gene slide should display method text when provided."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        # Create a simple mock figure
        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)
        blank_layout = prs.slide_layouts[6]

        method_text = "HDF cells treated with test compound for 24h"

        # Mock to_image to avoid kaleido dependency
        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            ReportGenerator._add_gene_slide(
                prs, blank_layout, "COL1A1", fig, processed_gene_data, method_text
            )

        slide = prs.slides[0]
        method_found = False
        for shape in slide.shapes:
            if hasattr(shape, "text_frame"):
                if method_text in shape.text_frame.text:
                    method_found = True
                    break
        assert method_found, "Method text not found on gene slide"

    def test_gene_slide_without_method_text(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """Gene slide should work correctly without method text."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)
        blank_layout = prs.slide_layouts[6]

        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            ReportGenerator._add_gene_slide(
                prs, blank_layout, "COL1A1", fig, processed_gene_data, ""
            )

        slide = prs.slides[0]
        # Should have gene title
        title_found = False
        for shape in slide.shapes:
            if hasattr(shape, "text_frame"):
                if "COL1A1 Expression" in shape.text_frame.text:
                    title_found = True
                    break
        assert title_found, "Gene title not found on slide"

    def test_gene_slide_has_statistics_sidebar(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """Gene slide should have statistics sidebar when gene_data is provided."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)
        blank_layout = prs.slide_layouts[6]

        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            ReportGenerator._add_gene_slide(
                prs, blank_layout, "COL1A1", fig, processed_gene_data, ""
            )

        slide = prs.slides[0]
        stats_found = False
        for shape in slide.shapes:
            if hasattr(shape, "text_frame"):
                if "Statistics" in shape.text_frame.text:
                    stats_found = True
                    break
        assert stats_found, "Statistics sidebar not found on gene slide"

    def test_gene_slide_has_red_accent_bar(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """Gene slide should have COSMAX red accent bar at bottom."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches
        from pptx.enum.shapes import MSO_SHAPE_TYPE
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)
        blank_layout = prs.slide_layouts[6]

        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            ReportGenerator._add_gene_slide(
                prs, blank_layout, "COL1A1", fig, processed_gene_data, ""
            )

        slide = prs.slides[0]
        rect_found = False
        for shape in slide.shapes:
            if shape.shape_type == MSO_SHAPE_TYPE.AUTO_SHAPE:
                rect_found = True
                break
        assert rect_found, "Red accent bar not found on gene slide"


class TestReportGeneratorCreatePresentation:
    """Tests for create_presentation() method."""

    def test_create_presentation_returns_bytes(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """create_presentation should return bytes (BytesIO content)."""
        from importlib import import_module
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        graphs = {"COL1A1": fig}
        processed_data = {"COL1A1": processed_gene_data}
        analysis_params = {
            "Efficacy_Type": "Anti-Aging",
            "Housekeeping_Gene": "GAPDH",
            "Reference_Sample": "Non-treated",
        }

        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            result = ReportGenerator.create_presentation(
                graphs=graphs,
                processed_data=processed_data,
                analysis_params=analysis_params,
                layout="one_per_slide",
                include_title_slide=True,
                include_summary=False,
                method_text="",
            )

        assert isinstance(result, bytes)
        assert len(result) > 0

    def test_create_presentation_with_method_text(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """create_presentation should pass method_text to gene slides."""
        from importlib import import_module
        from pptx import Presentation
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        graphs = {"COL1A1": fig}
        processed_data = {"COL1A1": processed_gene_data}
        analysis_params = {
            "Efficacy_Type": "Anti-Aging",
            "Housekeeping_Gene": "GAPDH",
            "Reference_Sample": "Non-treated",
        }
        method_text = "Test method description"

        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            result = ReportGenerator.create_presentation(
                graphs=graphs,
                processed_data=processed_data,
                analysis_params=analysis_params,
                layout="one_per_slide",
                include_title_slide=True,
                include_summary=False,
                method_text=method_text,
            )

        # Load the presentation and verify method text is present
        prs = Presentation(io.BytesIO(result))
        method_found = False
        for slide in prs.slides:
            for shape in slide.shapes:
                if hasattr(shape, "text_frame"):
                    if method_text in shape.text_frame.text:
                        method_found = True
                        break
            if method_found:
                break
        assert method_found, "Method text not found in generated presentation"

    def test_create_presentation_slide_count(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """create_presentation should create correct number of slides."""
        from importlib import import_module
        from pptx import Presentation
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        graphs = {"COL1A1": fig}
        processed_data = {"COL1A1": processed_gene_data}
        analysis_params = {
            "Efficacy_Type": "Anti-Aging",
            "Housekeeping_Gene": "GAPDH",
            "Reference_Sample": "Non-treated",
        }

        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            result = ReportGenerator.create_presentation(
                graphs=graphs,
                processed_data=processed_data,
                analysis_params=analysis_params,
                layout="one_per_slide",
                include_title_slide=True,
                include_summary=False,
                method_text="",
            )

        prs = Presentation(io.BytesIO(result))
        # Should have: 1 title + 1 gene = 2 slides (summary disabled)
        assert len(prs.slides) == 2


class TestReportGeneratorEdgeCases:
    """Tests for edge cases in PPT generation."""

    def test_gene_slide_with_empty_gene_data(self, mock_streamlit, graph_settings):
        """Gene slide should handle empty gene_data gracefully."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        # Create a simple mock figure
        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)
        blank_layout = prs.slide_layouts[6]

        empty_df = pd.DataFrame()

        # Should not raise an exception
        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            ReportGenerator._add_gene_slide(
                prs, blank_layout, "TEST_GENE", fig, empty_df, ""
            )

        slide = prs.slides[0]
        assert slide is not None

    def test_gene_slide_with_none_gene_data(self, mock_streamlit, graph_settings):
        """Gene slide should handle None gene_data gracefully."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)
        blank_layout = prs.slide_layouts[6]

        # Should not raise an exception
        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            ReportGenerator._add_gene_slide(
                prs, blank_layout, "TEST_GENE", fig, None, ""
            )

        slide = prs.slides[0]
        assert slide is not None

    def test_gene_slide_with_long_method_text(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """Gene slide should handle long method text without breaking layout."""
        from importlib import import_module
        from pptx import Presentation
        from pptx.util import Inches
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)
        blank_layout = prs.slide_layouts[6]

        # Very long method text
        long_method_text = (
            "HDF cells were seeded at 1x10^5 cells/well in 6-well plates and "
            "treated with test compound at various concentrations (0.1, 1, 10, "
            "100 μM) for 24 hours. RNA was extracted using TRIzol reagent and "
            "reverse transcribed using SuperScript IV."
        )

        # Should not raise an exception
        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            ReportGenerator._add_gene_slide(
                prs, blank_layout, "COL1A1", fig, processed_gene_data, long_method_text
            )

        slide = prs.slides[0]
        assert slide is not None
        # Verify method text is present (may be truncated)
        method_found = False
        for shape in slide.shapes:
            if hasattr(shape, "text_frame"):
                if "HDF cells" in shape.text_frame.text:
                    method_found = True
                    break
        assert method_found, "Long method text not found on slide"

    def test_create_presentation_with_empty_graphs(self, mock_streamlit):
        """create_presentation should handle empty graphs dict."""
        from importlib import import_module
        from pptx import Presentation
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        graphs = {}
        processed_data = {}
        analysis_params = {
            "Efficacy_Type": "Anti-Aging",
            "Housekeeping_Gene": "GAPDH",
            "Reference_Sample": "Non-treated",
        }

        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            result = ReportGenerator.create_presentation(
                graphs=graphs,
                processed_data=processed_data,
                analysis_params=analysis_params,
                layout="one_per_slide",
                include_title_slide=True,
                include_summary=True,
                method_text="",
            )

        prs = Presentation(io.BytesIO(result))
        # Should have: 1 title + 0 genes + 1 summary = 2 slides
        assert len(prs.slides) == 2

    def test_create_presentation_without_title_slide(
        self, mock_streamlit, processed_gene_data, graph_settings
    ):
        """create_presentation should skip title slide when include_title_slide=False."""
        from importlib import import_module
        from pptx import Presentation
        import plotly.graph_objects as go

        spec = import_module("streamlit qpcr analysis v1")
        ReportGenerator = spec.ReportGenerator

        fig = go.Figure(data=[go.Bar(x=["A", "B"], y=[1, 2])])

        graphs = {"COL1A1": fig}
        processed_data = {"COL1A1": processed_gene_data}
        analysis_params = {
            "Efficacy_Type": "Anti-Aging",
            "Housekeeping_Gene": "GAPDH",
            "Reference_Sample": "Non-treated",
        }

        with patch.object(go.Figure, "to_image", return_value=MOCK_PNG_BYTES):
            result = ReportGenerator.create_presentation(
                graphs=graphs,
                processed_data=processed_data,
                analysis_params=analysis_params,
                layout="one_per_slide",
                include_title_slide=False,
                include_summary=False,
                method_text="",
            )

        prs = Presentation(io.BytesIO(result))
        # Should have: 0 title + 1 gene = 1 slide (summary disabled)
        assert len(prs.slides) == 1
