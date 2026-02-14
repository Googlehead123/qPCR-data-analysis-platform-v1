"""ReportGenerator & PPTGenerator — PowerPoint report generation.

ReportGenerator: Legacy PPTX with Cosmax branding
PPTGenerator: Korean-centric template-based PPTX with efficacy config
"""

import io
from datetime import datetime
from typing import Dict

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from qpcr.constants import PLOTLY_FONT_FAMILY, CM_TO_PX, CM_TO_EMU, EFFICACY_CONFIG


class ReportGenerator:
    SLIDE_WIDTH_INCHES = 13.333
    SLIDE_HEIGHT_INCHES = 7.5

    @staticmethod
    def _fig_to_image(fig: go.Figure, format: str = "png", scale: int = 2) -> bytes:
        """Convert Plotly figure to image bytes with proper error handling for Kaleido/Chrome."""
        import os

        chrome_paths = [
            "/usr/bin/chromium",
            "/usr/bin/chromium-browser",
            "/usr/bin/google-chrome",
            "/usr/bin/google-chrome-stable",
        ]
        for chrome_path in chrome_paths:
            if os.path.exists(chrome_path):
                os.environ["CHROME_PATH"] = chrome_path
                break

        try:
            return fig.to_image(format=format, scale=scale)
        except Exception as e:
            error_msg = str(e)
            if "Chrome" in error_msg or "chromium" in error_msg.lower():
                raise RuntimeError(
                    "Image export requires Chrome/Chromium. "
                    "On Streamlit Cloud, add 'chromium' to packages.txt. "
                    "Locally, install Chrome or run: plotly_get_chrome"
                ) from e
            raise

    @staticmethod
    def create_presentation(
        graphs: Dict[str, go.Figure],
        processed_data: Dict[str, pd.DataFrame],
        analysis_params: dict,
        layout: str = "one_per_slide",
        include_title_slide: bool = True,
        include_summary: bool = True,
        method_text: str = "",
    ) -> bytes:
        try:
            from pptx import Presentation
            from pptx.util import Inches, Emu
            from pptx.dml.color import RGBColor
        except ImportError:
            raise ImportError(
                "python-pptx is required. Install with: pip install python-pptx"
            )

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)

        blank_layout = prs.slide_layouts[6]

        if include_title_slide:
            ReportGenerator._add_title_slide(prs, analysis_params)

        for gene, fig in graphs.items():
            gene_data = processed_data.get(gene)
            ReportGenerator._add_gene_slide(
                prs, blank_layout, gene, fig, gene_data, method_text
            )

        if include_summary:
            ReportGenerator._add_summary_slide(
                prs, blank_layout, processed_data, analysis_params
            )

        output = io.BytesIO()
        prs.save(output)
        output.seek(0)
        return output.getvalue()

    @staticmethod
    def _add_title_slide(prs, analysis_params: dict):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN
        from pptx.enum.shapes import MSO_SHAPE
        from pptx.dml.color import RGBColor

        slide = prs.slides.add_slide(prs.slide_layouts[6])

        background = slide.background
        fill = background.fill
        fill.solid()
        fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

        logo_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(0.3), Inches(2), Inches(0.8)
        )
        logo_frame = logo_box.text_frame
        logo_frame.word_wrap = True
        logo_para = logo_frame.paragraphs[0]
        logo_para.text = "[Add Logo Here]"
        logo_para.font.size = Pt(12)
        logo_para.font.color.rgb = RGBColor(0xC1, 0xC6, 0xC7)

        title_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(2.5), Inches(12.33), Inches(1)
        )
        title_frame = title_box.text_frame
        title_para = title_frame.paragraphs[0]
        title_para.text = "qPCR Gene Expression Analysis"
        title_para.font.size = Pt(54)
        title_para.font.bold = True
        title_para.font.color.rgb = RGBColor(0x00, 0x00, 0x00)
        title_para.alignment = PP_ALIGN.CENTER

        subtitle_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(3.8), Inches(12.33), Inches(0.5)
        )
        subtitle_frame = subtitle_box.text_frame
        subtitle_para = subtitle_frame.paragraphs[0]
        efficacy = analysis_params.get("Efficacy_Type", "Analysis")
        subtitle_para.text = f"{efficacy} Efficacy Study"
        subtitle_para.font.size = Pt(32)
        subtitle_para.font.color.rgb = RGBColor(0x00, 0x00, 0x00)
        subtitle_para.alignment = PP_ALIGN.CENTER

        info_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(5.5), Inches(12.33), Inches(1)
        )
        info_frame = info_box.text_frame
        info_para = info_frame.paragraphs[0]
        date = analysis_params.get("Date", datetime.now().strftime("%Y-%m-%d"))
        hk = analysis_params.get("Housekeeping_Gene", "N/A")
        ref = analysis_params.get("Reference_Sample", "N/A")
        info_para.text = f"Date: {date}  |  HK Gene: {hk}  |  Reference: {ref}"
        info_para.font.size = Pt(14)
        info_para.font.color.rgb = RGBColor(0xC1, 0xC6, 0xC7)
        info_para.alignment = PP_ALIGN.CENTER

        accent_bar = slide.shapes.add_shape(
            MSO_SHAPE.RECTANGLE, Inches(0), Inches(7.0), Inches(13.333), Inches(0.5)
        )
        accent_bar.fill.solid()
        accent_bar.fill.fore_color.rgb = RGBColor(0xEA, 0x1D, 0x22)
        accent_bar.line.fill.background()

    @staticmethod
    def _add_gene_slide(
        prs,
        layout,
        gene: str,
        fig: go.Figure,
        gene_data: pd.DataFrame = None,
        method_text: str = "",
    ):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN
        from pptx.enum.shapes import MSO_SHAPE
        from pptx.dml.color import RGBColor

        slide = prs.slides.add_slide(layout)

        background = slide.background
        fill = background.fill
        fill.solid()
        fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

        title_box = slide.shapes.add_textbox(
            Inches(0.3), Inches(0.2), Inches(12.73), Inches(0.6)
        )
        title_frame = title_box.text_frame
        title_para = title_frame.paragraphs[0]
        title_para.text = f"{gene} Expression"
        title_para.font.size = Pt(28)
        title_para.font.bold = True
        title_para.font.color.rgb = RGBColor(0x00, 0x00, 0x00)
        title_para.alignment = PP_ALIGN.LEFT

        graph_top = Inches(0.9)
        if method_text:
            method_box = slide.shapes.add_textbox(
                Inches(0.3), Inches(0.7), Inches(12.73), Inches(0.4)
            )
            method_frame = method_box.text_frame
            method_frame.word_wrap = True
            method_para = method_frame.paragraphs[0]
            method_para.text = method_text
            method_para.font.size = Pt(11)
            method_para.font.color.rgb = RGBColor(0xC1, 0xC6, 0xC7)
            method_para.alignment = PP_ALIGN.LEFT
            graph_top = Inches(1.2)

        fig_copy = go.Figure(fig)
        fig_copy.update_layout(
            width=1000,
            height=550,
            margin=dict(l=60, r=60, t=60, b=80),
            font=dict(size=14, family=PLOTLY_FONT_FAMILY),
        )

        img_bytes = ReportGenerator._fig_to_image(fig_copy, format="png", scale=2)
        img_stream = io.BytesIO(img_bytes)

        left = Inches(0.5)
        top = graph_top
        width = Inches(9.5)
        slide.shapes.add_picture(img_stream, left, top, width=width)

        if gene_data is not None and not gene_data.empty:
            stats_box = slide.shapes.add_textbox(
                Inches(10.2), Inches(1.0), Inches(2.8), Inches(5.5)
            )
            stats_frame = stats_box.text_frame
            stats_frame.word_wrap = True

            header_para = stats_frame.paragraphs[0]
            header_para.text = "Statistics"
            header_para.font.size = Pt(14)
            header_para.font.bold = True

            for _, row in gene_data.iterrows():
                cond = row.get("Condition", "N/A")
                fc = row.get("Fold_Change", row.get("Relative_Expression", 0))
                pval = row.get("p_value", float("nan"))
                sig = row.get("significance", "")

                para = stats_frame.add_paragraph()
                para.text = f"{cond[:12]}"
                para.font.size = Pt(10)
                para.font.bold = True

                para2 = stats_frame.add_paragraph()
                fc_str = f"FC: {fc:.2f}" if pd.notna(fc) else "FC: N/A"
                p_str = f"p={pval:.3f}" if pd.notna(pval) else ""
                para2.text = f"  {fc_str} {sig}"
                para2.font.size = Pt(9)

        accent_bar = slide.shapes.add_shape(
            MSO_SHAPE.RECTANGLE, Inches(0), Inches(7.0), Inches(13.333), Inches(0.5)
        )
        accent_bar.fill.solid()
        accent_bar.fill.fore_color.rgb = RGBColor(0xEA, 0x1D, 0x22)
        accent_bar.line.fill.background()

    @staticmethod
    def _add_dual_gene_slide(prs, layout, gene_pairs: list, processed_data: dict):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN

        slide = prs.slides.add_slide(layout)

        positions = [
            (Inches(0.3), Inches(0.8), Inches(6.2)),
            (Inches(6.8), Inches(0.8), Inches(6.2)),
        ]

        for idx, (gene, fig) in enumerate(gene_pairs):
            if idx >= 2:
                break

            left, top, width = positions[idx]

            title_box = slide.shapes.add_textbox(left, Inches(0.2), width, Inches(0.5))
            title_frame = title_box.text_frame
            title_para = title_frame.paragraphs[0]
            title_para.text = f"{gene}"
            title_para.font.size = Pt(20)
            title_para.font.bold = True

            fig_copy = go.Figure(fig)
            fig_copy.update_layout(
                width=600,
                height=450,
                margin=dict(l=50, r=30, t=40, b=60),
                font=dict(size=11, family=PLOTLY_FONT_FAMILY),
            )

            img_bytes = ReportGenerator._fig_to_image(fig_copy, format="png", scale=2)
            img_stream = io.BytesIO(img_bytes)

            slide.shapes.add_picture(img_stream, left, top, width=width)

    @staticmethod
    def _add_grid_slide(prs, layout, graphs: dict):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN

        slide = prs.slides.add_slide(layout)

        title_box = slide.shapes.add_textbox(
            Inches(0.3), Inches(0.1), Inches(12.73), Inches(0.5)
        )
        title_frame = title_box.text_frame
        title_para = title_frame.paragraphs[0]
        title_para.text = "Gene Expression Overview"
        title_para.font.size = Pt(24)
        title_para.font.bold = True

        gene_list = list(graphs.items())
        n_genes = len(gene_list)

        if n_genes <= 2:
            cols, rows = 2, 1
        elif n_genes <= 4:
            cols, rows = 2, 2
        elif n_genes <= 6:
            cols, rows = 3, 2
        else:
            cols, rows = 4, 2

        cell_width = 12.5 / cols
        cell_height = 6.5 / rows

        for idx, (gene, fig) in enumerate(gene_list[: cols * rows]):
            col_idx = idx % cols
            row_idx = idx // cols

            left = Inches(0.4 + col_idx * cell_width)
            top = Inches(0.7 + row_idx * cell_height)

            fig_copy = go.Figure(fig)
            fig_copy.update_layout(
                width=350,
                height=280,
                margin=dict(l=40, r=20, t=35, b=40),
                font=dict(size=9, family=PLOTLY_FONT_FAMILY),
                title=dict(text=gene, font=dict(size=12, family=PLOTLY_FONT_FAMILY)),
            )

            img_bytes = ReportGenerator._fig_to_image(fig_copy, format="png", scale=2)
            img_stream = io.BytesIO(img_bytes)

            slide.shapes.add_picture(
                img_stream, left, top, width=Inches(cell_width - 0.2)
            )

    @staticmethod
    def _add_summary_slide(prs, layout, processed_data: dict, analysis_params: dict):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN

        slide = prs.slides.add_slide(layout)

        title_box = slide.shapes.add_textbox(
            Inches(0.3), Inches(0.2), Inches(12.73), Inches(0.6)
        )
        title_frame = title_box.text_frame
        title_para = title_frame.paragraphs[0]
        title_para.text = "Analysis Summary"
        title_para.font.size = Pt(28)
        title_para.font.bold = True

        all_results = (
            pd.concat(processed_data.values(), ignore_index=True)
            if processed_data
            else pd.DataFrame()
        )

        summary_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(1.0), Inches(6), Inches(5.5)
        )
        summary_frame = summary_box.text_frame
        summary_frame.word_wrap = True

        header = summary_frame.paragraphs[0]
        header.text = "Experimental Parameters"
        header.font.size = Pt(16)
        header.font.bold = True

        params_text = [
            f"Efficacy Type: {analysis_params.get('Efficacy_Type', 'N/A')}",
            f"Housekeeping Gene: {analysis_params.get('Housekeeping_Gene', 'N/A')}",
            f"Reference Condition: {analysis_params.get('Reference_Sample', 'N/A')}",
            f"Comparison Condition: {analysis_params.get('Compare_To', 'N/A')}",
            f"Genes Analyzed: {len(processed_data)}",
            f"Analysis Date: {analysis_params.get('Date', 'N/A')}",
        ]

        for text in params_text:
            para = summary_frame.add_paragraph()
            para.text = text
            para.font.size = Pt(12)

        if not all_results.empty:
            results_box = slide.shapes.add_textbox(
                Inches(7), Inches(1.0), Inches(5.8), Inches(5.5)
            )
            results_frame = results_box.text_frame
            results_frame.word_wrap = True

            header2 = results_frame.paragraphs[0]
            header2.text = "Key Findings"
            header2.font.size = Pt(16)
            header2.font.bold = True

            for gene, gene_df in processed_data.items():
                if gene_df.empty:
                    continue

                para = results_frame.add_paragraph()
                para.text = f"\n{gene}:"
                para.font.size = Pt(12)
                para.font.bold = True

                sig_results = gene_df[
                    gene_df.get("significance", pd.Series([""])).str.len() > 0
                ]
                if not sig_results.empty:
                    for _, row in sig_results.head(3).iterrows():
                        cond = row.get("Condition", "N/A")
                        fc = row.get("Fold_Change", row.get("Relative_Expression", 0))
                        sig = row.get("significance", "")
                        para2 = results_frame.add_paragraph()
                        para2.text = f"  {cond}: {fc:.2f}x {sig}"
                        para2.font.size = Pt(10)
                else:
                    para2 = results_frame.add_paragraph()
                    para2.text = "  No significant changes"
                    para2.font.size = Pt(10)


class PPTGenerator:
    NAVY_BLUE = (0, 0, 0)
    WHITE = (255, 255, 255)

    @staticmethod
    def _get_color_rgb(rgb_tuple):
        from pptx.dml.color import RGBColor

        return RGBColor(*rgb_tuple)

    @staticmethod
    def create_title_slide(prs, analysis_params):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN
        from pptx.dml.color import RGBColor

        slide = prs.slides.add_slide(prs.slide_layouts[6])

        background = slide.background
        fill = background.fill
        fill.solid()
        fill.fore_color.rgb = PPTGenerator._get_color_rgb(PPTGenerator.WHITE)

        title_box = slide.shapes.add_textbox(
            Inches(1), Inches(2.5), Inches(11.33), Inches(1)
        )
        tf = title_box.text_frame
        p = tf.paragraphs[0]
        p.text = "유전자 발현 분석 (Gene Expression Analysis)"
        p.font.size = Pt(40)
        p.font.bold = True
        p.font.name = "Malgun Gothic"
        p.alignment = PP_ALIGN.CENTER

        subtitle_box = slide.shapes.add_textbox(
            Inches(1), Inches(3.8), Inches(11.33), Inches(0.5)
        )
        tf = subtitle_box.text_frame
        p = tf.paragraphs[0]
        efficacy = analysis_params.get("Efficacy_Type", "Analysis")
        p.text = f"{efficacy} Efficacy Test"
        p.font.size = Pt(28)
        p.font.name = "Malgun Gothic"
        p.alignment = PP_ALIGN.CENTER

        info_box = slide.shapes.add_textbox(
            Inches(1), Inches(5.0), Inches(11.33), Inches(1)
        )
        tf = info_box.text_frame
        p = tf.paragraphs[0]
        date = analysis_params.get("Date", "")
        p.text = f"Date: {date}"
        p.font.size = Pt(18)
        p.font.name = "Arial"
        p.alignment = PP_ALIGN.CENTER

        logo_box = slide.shapes.add_textbox(
            Inches(11.5), Inches(6.5), Inches(1.5), Inches(0.5)
        )
        tf = logo_box.text_frame
        p = tf.paragraphs[0]
        p.text = "[Logo]"
        p.font.size = Pt(12)
        p.font.color.rgb = RGBColor(200, 200, 200)
        p.alignment = PP_ALIGN.RIGHT

    @staticmethod
    def create_gene_slide(prs, gene, fig, gene_data, analysis_params):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN
        from pptx.enum.shapes import MSO_SHAPE
        from pptx.dml.color import RGBColor

        slide = prs.slides.add_slide(prs.slide_layouts[6])

        background = slide.background
        fill = background.fill
        fill.solid()
        fill.fore_color.rgb = PPTGenerator._get_color_rgb(PPTGenerator.WHITE)

        title_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(0.3), Inches(9.0), Inches(0.8)
        )
        tf = title_box.text_frame
        p = tf.paragraphs[0]
        p.text = f"{gene} Expression"
        p.font.size = Pt(32)
        p.font.bold = True
        p.font.name = "Malgun Gothic"

        line = slide.shapes.add_shape(
            MSO_SHAPE.RECTANGLE, Inches(0.5), Inches(1.1), Inches(9.0), Inches(0.05)
        )
        line.fill.solid()
        line.fill.fore_color.rgb = PPTGenerator._get_color_rgb(PPTGenerator.NAVY_BLUE)
        line.line.fill.background()

        try:
            img_bytes = fig.to_image(format="png", scale=2, width=1200, height=900)
            image_stream = io.BytesIO(img_bytes)
            slide.shapes.add_picture(
                image_stream,
                Inches(0.5),
                Inches(1.5),
                width=Inches(6.0),
                height=Inches(4.0),
            )
        except Exception as e:
            err_box = slide.shapes.add_textbox(
                Inches(0.5), Inches(1.5), Inches(6.0), Inches(4.0)
            )
            err_box.text = f"Graph Error: {str(e)}"

        if gene_data is not None and not gene_data.empty:
            rows = len(gene_data) + 1
            cols = 4
            table_height = Inches(0.3 * rows)
            table_shape = slide.shapes.add_table(
                rows, cols, Inches(0.5), Inches(5.8), Inches(6.0), table_height
            )
            table = table_shape.table

            headers = ["시료명 (Sample)", "Condition", "Fold Change", "P-value"]
            for i, h in enumerate(headers):
                cell = table.cell(0, i)
                cell.text = h
                cell.text_frame.paragraphs[0].font.size = Pt(10)
                cell.text_frame.paragraphs[0].font.bold = True

            for i, row in enumerate(gene_data.itertuples()):
                r = i + 1
                sample_val = getattr(row, "Original_Sample", getattr(row, "Sample", ""))
                cond_val = getattr(row, "Condition", "")
                fc_val = getattr(
                    row, "Fold_Change", getattr(row, "Relative_Expression", 0)
                )
                pval = getattr(row, "p_value", None)
                sig = getattr(row, "significance", "")

                table.cell(r, 0).text = str(sample_val)
                table.cell(r, 1).text = str(cond_val)
                table.cell(r, 2).text = f"{fc_val:.2f}" if pd.notna(fc_val) else "-"
                table.cell(r, 3).text = f"{pval:.4f} {sig}" if pd.notna(pval) else "-"

                for j in range(4):
                    table.cell(r, j).text_frame.paragraphs[0].font.size = Pt(9)

        meta_box = slide.shapes.add_textbox(
            Inches(7.0), Inches(1.5), Inches(2.5), Inches(5.0)
        )
        tf = meta_box.text_frame
        tf.word_wrap = True

        def add_meta_line(text, bold=False):
            p = tf.add_paragraph()
            p.text = text
            p.font.size = Pt(11)
            p.font.bold = bold

        add_meta_line("Analysis Parameters", bold=True)
        add_meta_line(f"HK Gene: {analysis_params.get('Housekeeping_Gene', '-')}")
        add_meta_line(f"Ref Sample: {analysis_params.get('Reference_Sample', '-')}")
        add_meta_line(f"Compare: {analysis_params.get('Compare_To', '-')}")
        add_meta_line("")
        add_meta_line("통계적 유의성 (Significance)", bold=True)
        add_meta_line("* p < 0.05")
        add_meta_line("** p < 0.01")
        add_meta_line("*** p < 0.001")

        logo_box = slide.shapes.add_textbox(
            Inches(11.5), Inches(6.5), Inches(1.5), Inches(0.5)
        )
        tf = logo_box.text_frame
        p = tf.paragraphs[0]
        p.text = "[Logo]"
        p.font.size = Pt(12)
        p.font.color.rgb = RGBColor(200, 200, 200)
        p.alignment = PP_ALIGN.RIGHT

    @staticmethod
    def _delete_slide(prs, index):
        slides_list = prs.slides._sldIdLst
        rId = slides_list[index].attrib[
            '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
        ]
        prs.part.drop_rel(rId)
        slides_list.remove(slides_list[index])

    @staticmethod
    def _safe_copy_slide(prs, source_index):
        """Copy a slide by duplicating shapes via the spTree proxy.

        Uses python-pptx's spTree property to ensure shapes are properly
        recognized after copy. Copies relationships for images/media,
        then deep-copies each shape element from the source spTree.
        """
        import copy

        source = prs.slides[source_index]
        new_slide = prs.slides.add_slide(source.slide_layout)

        rid_map = {}
        for rid, rel in source.part.rels.items():
            rel_type = str(rel.reltype)
            if 'slideLayout' in rel_type:
                for new_rid, new_rel in new_slide.part.rels.items():
                    if 'slideLayout' in str(new_rel.reltype):
                        rid_map[rid] = new_rid
                        break
                continue
            if 'notesSlide' in rel_type:
                continue
            try:
                if rel.is_external:
                    new_rid = new_slide.part.rels.get_or_add_ext_rel(
                        rel.reltype, rel.target_ref
                    )
                else:
                    new_rid = new_slide.part.relate_to(
                        rel.target_part, rel.reltype
                    )
                rid_map[rid] = new_rid
            except Exception:
                pass

        source_spTree = source._element.spTree
        target_spTree = new_slide._element.spTree

        shape_tags = ['}sp', '}pic', '}cxnSp', '}grpSp']
        for child in list(target_spTree):
            if any(child.tag.endswith(t) for t in shape_tags):
                target_spTree.remove(child)

        for child in source_spTree:
            if any(child.tag.endswith(t) for t in shape_tags):
                new_child = copy.deepcopy(child)
                for elem in new_child.iter():
                    for attr_key in list(elem.attrib.keys()):
                        val = elem.attrib[attr_key]
                        if val in rid_map:
                            elem.attrib[attr_key] = rid_map[val]
                target_spTree.append(new_child)

        return new_slide

    @staticmethod
    def _move_slide_to_end(prs, slide_index):
        slides_list = prs.slides._sldIdLst
        el = slides_list[slide_index]
        slides_list.remove(el)
        slides_list.append(el)

    @staticmethod
    def _populate_gene_slide(prs, slide, gene, fig, gene_data, analysis_params, graph_settings):
        """Populate a template-duplicated gene slide with gene name and graph image."""
        from pptx.util import Inches, Pt, Emu

        for shape in slide.shapes:
            if not shape.has_text_frame:
                continue
            full_text = shape.text_frame.text

            # TextBox 17: "효능 평가" → "{gene} 효능 평가"
            if "효능 평가" in full_text and "Results" not in full_text and "有" not in full_text:
                for para in shape.text_frame.paragraphs:
                    for run in para.runs:
                        if "효능" in run.text:
                            run.text = run.text.replace("효능 평가", f"{gene} 효능 평가")
                            break
                    break

            # TextBox 14: "Results: 효능 有/無" → update with gene result
            elif "Results" in full_text or "효능 有" in full_text:
                has_efficacy = False
                if gene_data is not None and not gene_data.empty:
                    sig_vals = gene_data.get("significance", pd.Series())
                    has_efficacy = sig_vals.str.len().gt(0).any() if not sig_vals.empty else False
                result_text = f"Results: 효능 {'有' if has_efficacy else '無'}"
                for para in shape.text_frame.paragraphs:
                    runs = list(para.runs)
                    for i, run in enumerate(runs):
                        if i == 0:
                            run.text = result_text
                        else:
                            run.text = ""
                    break

            # TextBox 2: experiment details (top-right area) — update with params
            elif (shape.left is not None and shape.top is not None
                  and shape.left > 3500000 and shape.top < 500000):
                tf = shape.text_frame
                hk = analysis_params.get("Housekeeping_Gene", "-")
                ref = analysis_params.get("Reference_Sample", "-")
                cmp = analysis_params.get("Compare_To", "-")
                efficacy_t = analysis_params.get("Efficacy_Type", "-")
                for para in tf.paragraphs:
                    for run in para.runs:
                        run.text = ""
                if tf.paragraphs and tf.paragraphs[0].runs:
                    tf.paragraphs[0].runs[0].text = (
                        f"{efficacy_t} | HK: {hk} | Ref: {ref} | vs: {cmp}"
                    )

        # Add graph image
        gs = graph_settings or {}
        w_cm = gs.get("figure_width", 28)
        h_cm = gs.get("figure_height", 16)

        fig_copy = go.Figure(fig)
        fig_copy.update_layout(
            width=max(int(w_cm * CM_TO_PX), 800),
            height=max(int(h_cm * CM_TO_PX), 500),
        )

        try:
            img_bytes = ReportGenerator._fig_to_image(fig_copy, format="png", scale=2)
            img_stream = io.BytesIO(img_bytes)

            w_emu = int(w_cm * CM_TO_EMU)
            h_emu = int(h_cm * CM_TO_EMU)
            max_w_emu = int(11.5 * 914400)
            max_h_emu = int(5.2 * 914400)

            if w_emu > max_w_emu:
                scale_f = max_w_emu / w_emu
                w_emu = max_w_emu
                h_emu = int(h_emu * scale_f)
            if h_emu > max_h_emu:
                scale_f = max_h_emu / h_emu
                h_emu = max_h_emu
                w_emu = int(w_emu * scale_f)

            slide_w = int(prs.slide_width)
            left = max(0, int((slide_w - w_emu) / 2))
            top = int(1.0 * 914400)

            slide.shapes.add_picture(
                img_stream, left, top, width=Emu(w_emu), height=Emu(h_emu)
            )
        except Exception as e:
            err_box = slide.shapes.add_textbox(
                Inches(0.5), Inches(1.5), Inches(11.0), Inches(4.0)
            )
            err_box.text_frame.paragraphs[0].text = f"Graph Error: {str(e)}"
            err_box.text_frame.paragraphs[0].font.size = Pt(14)

    @staticmethod
    def generate_presentation(graphs, processed_data, analysis_params, graph_settings=None):
        try:
            from pptx import Presentation
            from pptx.util import Inches, Emu
        except ImportError:
            st.error("python-pptx not installed")
            return None

        import os

        script_dir = os.path.dirname(os.path.abspath(__file__))
        template_path = os.path.join(
            os.path.dirname(script_dir), "251215 효능평가 결과 TEMPLATE.pptx"
        )
        use_template = os.path.isfile(template_path)

        gs = graph_settings or {}
        gene_list = list(graphs.keys())
        n_genes = len(gene_list)

        if use_template and n_genes > 0:
            prs = Presentation(template_path)

            # Template slides: [0]=Title, [1]=Description, [2]=Gene template, [3]=Closing
            # Step 1: Update title slide — replace "소재" with efficacy type
            efficacy = analysis_params.get("Efficacy_Type", "소재")
            title_slide = prs.slides[0]
            for shape in title_slide.shapes:
                if shape.has_text_frame:
                    for para in shape.text_frame.paragraphs:
                        for run in para.runs:
                            if "소재" in run.text:
                                run.text = run.text.replace("소재", efficacy)

            # Step 2: Copy gene template (index 2) BEFORE deleting to avoid
            # ZIP part-name conflicts (delete frees a name that add_slide reuses)
            for _ in range(n_genes - 1):
                PPTGenerator._safe_copy_slide(prs, 2)
            # Now: [0]=Title, [1]=Desc, [2]=Gene, [3]=Closing, [4..N+2]=GeneCopies

            # Step 3: Delete slide 1 (description overview)
            PPTGenerator._delete_slide(prs, 1)
            # Now: [0]=Title, [1]=Gene, [2]=Closing, [3..N+1]=GeneCopies

            # Step 4: Move closing slide (index 2) to end
            PPTGenerator._move_slide_to_end(prs, 2)
            # Now: [0]=Title, [1]=Gene, [2..N]=GeneCopies, [N+1]=Closing

            # Step 5: Populate each gene slide with gene name and graph
            for i, gene in enumerate(gene_list):
                slide = prs.slides[1 + i]
                fig = graphs[gene]
                gene_data = processed_data.get(gene)
                PPTGenerator._populate_gene_slide(
                    prs, slide, gene, fig, gene_data, analysis_params, gs
                )
        else:
            # Fallback: no template available
            prs = Presentation()
            prs.slide_width = Emu(12192000)
            prs.slide_height = Emu(6858000)
            PPTGenerator.create_title_slide(prs, analysis_params)

            for gene in gene_list:
                fig = graphs[gene]
                gene_data = processed_data.get(gene)
                PPTGenerator.create_gene_slide(
                    prs, gene, fig, gene_data, analysis_params
                )

        output = io.BytesIO()
        prs.save(output)
        output.seek(0)
        return output.getvalue()
