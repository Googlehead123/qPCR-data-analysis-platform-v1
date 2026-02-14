"""GraphGenerator — Plotly bar chart visualizations for qPCR results.

Creates per-gene relative expression bar charts with error bars,
significance annotations, and customizable styling.
"""

import textwrap

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from qpcr.constants import DEFAULT_GROUP_COLORS, PLOTLY_FONT_FAMILY, CM_TO_PX


class GraphGenerator:
    @staticmethod
    def _wrap_text(text: str, width: int = 15) -> str:
        """Wrap text for x-axis labels using <br> for Plotly compatibility"""
        wrapped = textwrap.fill(text, width=width)
        return wrapped.replace("\n", "<br>")

    @staticmethod
    def create_gene_graph(
        data: pd.DataFrame,
        gene: str,
        settings: dict,
        efficacy_config: dict = None,
        sample_order: list = None,
        per_sample_overrides: dict = None,
        condition_colors: dict = None,
        display_gene_name: str = None,
        ref_line_value: float = None,
        ref_line_label: str = None,
    ) -> go.Figure:
        """Create individual graph for each gene with proper data handling"""

        # Guard against empty data
        if data is None or data.empty:
            fig = go.Figure()
            fig.add_annotation(text="No data available", showarrow=False)
            return fig

        gene_data = data.copy()

        if "Relative_Expression" not in gene_data.columns:
            if "Fold_Change" in gene_data.columns:
                gene_data["Relative_Expression"] = gene_data["Fold_Change"]
            else:
                st.error(
                    f"Missing Relative_Expression or Fold_Change column for {gene}"
                )
                fig = go.Figure()
                fig.add_annotation(
                    text=f"Missing data columns for {gene}", showarrow=False
                )
                return fig

        if "SEM" not in gene_data.columns:
            gene_data["SEM"] = 0

        # Use sample_order from mapping and deduplicate conditions while preserving order
        if sample_order:
            mapping = st.session_state.get("sample_mapping", {})
            condition_order = []
            seen_conditions = set()

            for sample in sample_order:
                if mapping.get(sample, {}).get("include", True):
                    cond = mapping.get(sample, {}).get("condition", sample)
                    if (
                        cond in gene_data["Condition"].unique()
                        and cond not in seen_conditions
                    ):
                        condition_order.append(cond)
                        seen_conditions.add(cond)

            for cond in gene_data["Condition"].unique():
                if cond not in seen_conditions:
                    condition_order.append(cond)
                    seen_conditions.add(cond)

            gene_data["Condition"] = pd.Categorical(
                gene_data["Condition"], categories=condition_order, ordered=True
            )
            gene_data = gene_data.sort_values("Condition")
        else:
            gene_data = gene_data.sort_values("Condition")

        gene_data_indexed = gene_data.reset_index(drop=True)

        condition_names = gene_data_indexed["Condition"].tolist()
        n_bars = len(gene_data_indexed)

        bar_colors = []

        for idx, row in gene_data_indexed.iterrows():
            condition = row["Condition"]
            group = row.get("Group", "Treatment")

            custom_key = f"{gene}_{condition}"
            if custom_key in settings.get("bar_colors_per_sample", {}):
                bar_colors.append(settings["bar_colors_per_sample"][custom_key])
            elif condition_colors and condition in condition_colors:
                bar_colors.append(condition_colors[condition])
            elif group in DEFAULT_GROUP_COLORS:
                bar_colors.append(DEFAULT_GROUP_COLORS[group])
            else:
                default_color = settings.get("bar_colors", {}).get(gene, "#D3D3D3")
                bar_colors.append(default_color)

        fig = go.Figure()

        # Error bars - use SEM or SD based on user preference
        error_bar_type = st.session_state.get("error_bar_type", "sem")
        error_col = "SD" if error_bar_type == "sd" else "SEM"
        if error_col not in gene_data_indexed.columns:
            error_col = "SEM"
        error_array = gene_data_indexed[error_col].values

        gene_bar_settings = st.session_state.get(f"{gene}_bar_settings", {})

        show_error_global = settings.get("show_error", True)
        show_sig_global = settings.get("show_significance", True)

        error_visible_array = []

        for idx in range(n_bars):
            row = gene_data_indexed.iloc[idx]
            condition = row["Condition"]
            bar_key = f"{gene}_{condition}"

            bar_config = gene_bar_settings.get(
                bar_key, {"show_sig": True, "show_err": True}
            )

            if show_error_global and bar_config.get("show_err", True):
                error_visible_array.append(error_array[idx])
            else:
                error_visible_array.append(0)

        fig.add_trace(
            go.Bar(
                x=list(range(n_bars)),
                y=gene_data_indexed["Relative_Expression"],
                error_y=dict(
                    type="data",
                    array=error_visible_array,
                    arrayminus=[0] * n_bars,
                    visible=True,
                    thickness=2,
                    width=4,
                    color="rgba(0,0,0,0.5)",
                    symmetric=False,
                ),
                marker=dict(
                    color=bar_colors,
                    line=dict(
                        width=settings.get("marker_line_width", 1), color="black"
                    ),
                    opacity=settings.get("bar_opacity", 0.95),
                ),
                showlegend=False,
            )
        )

        max_y_value = gene_data_indexed["Relative_Expression"].max()
        if pd.isna(max_y_value) or max_y_value <= 0:
            max_y_value = 1.0  # Fallback for all-NaN or zero expression
        max_error = error_array.max() if len(error_array) > 0 else 0
        if pd.isna(max_error):
            max_error = 0
        y_max_auto = (
            max_y_value + max_error + (max_y_value * 0.15)
        )
        if ref_line_value is not None and pd.notna(ref_line_value):
            y_max_auto = max(y_max_auto, ref_line_value * 1.20)

        fixed_symbol_spacing = y_max_auto * 0.05

        # Add significance symbols
        for idx in range(n_bars):
            row = gene_data_indexed.iloc[idx]
            condition = row["Condition"]
            bar_key = f"{gene}_{condition}"
            bar_config = gene_bar_settings.get(
                bar_key, {"show_sig": True, "show_err": True}
            )

            sig_1 = row.get("significance", "")
            sig_2 = row.get("significance_2", "")

            bar_height = row["Relative_Expression"]
            error_bar_height = error_visible_array[idx]
            base_y_position = bar_height + error_bar_height

            asterisk_font_size = 16
            hashtag_font_size = 10

            if show_sig_global and bar_config.get("show_sig", True):
                symbols_to_show = []
                font_sizes = []

                if sig_1 in ["*", "**", "***"]:
                    symbols_to_show.append(sig_1)
                    font_sizes.append(asterisk_font_size)

                if sig_2 in ["#", "##", "###"]:
                    symbols_to_show.append(sig_2)
                    font_sizes.append(hashtag_font_size)

                if len(symbols_to_show) == 2:
                    fig.add_annotation(
                        x=idx,
                        y=base_y_position + (fixed_symbol_spacing * 0.2),
                        text=symbols_to_show[0],
                        showarrow=False,
                        font=dict(size=font_sizes[0], color="black", family=PLOTLY_FONT_FAMILY),
                        xref="x",
                        yref="y",
                        xanchor="center",
                        yanchor="bottom",
                    )

                    fig.add_annotation(
                        x=idx,
                        y=base_y_position
                        + (fixed_symbol_spacing * 0.2)
                        + fixed_symbol_spacing,
                        text=symbols_to_show[1],
                        showarrow=False,
                        font=dict(size=font_sizes[1], color="black", family=PLOTLY_FONT_FAMILY),
                        xref="x",
                        yref="y",
                        xanchor="center",
                        yanchor="bottom",
                    )

                elif len(symbols_to_show) == 1:
                    fig.add_annotation(
                        x=idx,
                        y=base_y_position + (fixed_symbol_spacing * 0.2),
                        text=symbols_to_show[0],
                        showarrow=False,
                        font=dict(size=font_sizes[0], color="black", family=PLOTLY_FONT_FAMILY),
                        xref="x",
                        yref="y",
                        xanchor="center",
                        yanchor="bottom",
                    )

        gene_label = display_gene_name if display_gene_name else gene
        y_label_html = f"Relative <b style='color:red;'>{gene_label}</b> Expression Level"

        y_axis_config = dict(
            title=dict(
                text=y_label_html,
                font=dict(size=settings.get(f"{gene}_ylabel_size", 14), family=PLOTLY_FONT_FAMILY),
                standoff=15,
            ),
            showgrid=False,
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor="black",
            showline=True,
            linewidth=1,
            linecolor="black",
            mirror=False,
            range=[0, y_max_auto],
            fixedrange=False,
        )

        if settings.get("y_log_scale"):
            y_axis_config["type"] = "log"
            y_axis_config.pop("range", None)

        if settings.get("y_min") is not None or settings.get("y_max") is not None:
            y_range = []
            y_range.append(settings.get("y_min", 0))
            y_range.append(settings.get("y_max", y_max_auto))
            y_axis_config["range"] = y_range

        gene_bar_gap = settings.get(f"{gene}_bar_gap", settings.get("bar_gap", 0.15))
        gene_margins = settings.get(
            f"{gene}_margins", {"l": 80, "r": 40, "t": 60, "b": 200}
        )
        gene_bg_color = settings.get(
            f"{gene}_bg_color", settings.get("plot_bgcolor", "#FFFFFF")
        )
        gene_tick_size = settings.get(f"{gene}_tick_size", 12)

        wrapped_labels = [
            GraphGenerator._wrap_text(str(cond), 15) for cond in condition_names
        ]

        # FIX-19: Dynamic bottom margin based on max label line count
        max_label_lines = max(
            (label.count("<br>") + 1 for label in wrapped_labels), default=1
        )
        dynamic_b_margin = 180 + max(0, max_label_lines - 1) * 22
        default_margins = gene_margins.copy()
        if default_margins.get("b", 200) < dynamic_b_margin:
            default_margins["b"] = dynamic_b_margin
        gene_margins = default_margins

        # P-VALUE LEGEND
        cmp_ref_name = st.session_state.get("analysis_cmp_condition", "")
        legend_ref_label = f" (vs {cmp_ref_name})" if cmp_ref_name else ""
        legend_text = f"<b>Significance{legend_ref_label}:</b>  * p<0.05  ** p<0.01  *** p<0.001"

        if (
            "significance_2" in gene_data_indexed.columns
            and gene_data_indexed["significance_2"].notna().any()
        ):
            cmp_ref_name_2 = st.session_state.get("analysis_cmp_condition_2", "")
            legend_ref_label_2 = f" (vs {cmp_ref_name_2})" if cmp_ref_name_2 else ""
            legend_text += (
                f"<br><b>2nd Comparison{legend_ref_label_2}:</b>  # p<0.05  ## p<0.01  ### p<0.001"
            )

        fig_h_px = int(settings.get("figure_height", 16) * CM_TO_PX)
        b_margin_px = gene_margins.get("b", 200)
        t_margin_px = gene_margins.get("t", 60)
        plot_h_px = max(fig_h_px - b_margin_px - t_margin_px, 100)

        # Dynamic legend y: place below x-axis labels with gap
        label_px_est = max_label_lines * 18 + 10
        legend_y = -((label_px_est / plot_h_px) + 0.06)

        fig.update_layout(
            title="",
            xaxis=dict(
                title="",
                showgrid=False,
                zeroline=False,
                tickmode="array",
                tickvals=list(range(n_bars)),
                ticktext=wrapped_labels,
                tickfont=dict(size=gene_tick_size, family=PLOTLY_FONT_FAMILY),
                tickangle=0,
                ticks="outside",
                ticklen=8,
                tickcolor="rgba(0,0,0,0)",
                showline=False,
                mirror=False,
                side="bottom",
                range=[-0.5, n_bars - 0.5],
            ),
            yaxis=y_axis_config,
            template=settings.get("color_scheme", "plotly_white"),
            font=dict(size=settings.get("font_size", 14), family=PLOTLY_FONT_FAMILY),
            height=fig_h_px,
            width=int(settings.get("figure_width", 28) * CM_TO_PX),
            bargap=gene_bar_gap,
            showlegend=settings.get("show_legend", False),
            plot_bgcolor=gene_bg_color,
            paper_bgcolor="#FFFFFF",
            margin=dict(
                l=gene_margins.get("l", 80),
                r=gene_margins.get("r", 40),
                t=gene_margins.get("t", 60),
                b=b_margin_px,
            ),
        )

        fig.add_annotation(
            text=legend_text,
            xref="paper",
            yref="paper",
            x=1.0,
            y=legend_y,
            xanchor="right",
            yanchor="top",
            showarrow=False,
            font=dict(size=10, color="#888888", family=PLOTLY_FONT_FAMILY),
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="#DDDDDD",
            borderwidth=1,
            borderpad=3,
        )

        if ref_line_value is not None and pd.notna(ref_line_value):
            max_annotation_y = 0
            for _idx in range(n_bars):
                _row = gene_data_indexed.iloc[_idx]
                _bar_h = _row["Relative_Expression"]
                _err_h = error_visible_array[_idx] if _idx < len(error_visible_array) else 0
                _top_y = _bar_h + _err_h + (fixed_symbol_spacing * 1.2)
                if _top_y > max_annotation_y:
                    max_annotation_y = _top_y

            ann_position = "top right"
            if ref_line_value > max_annotation_y * 0.85:
                ann_position = "bottom right"

            fig.add_hline(
                y=ref_line_value,
                line_dash="dash",
                line_color="rgba(120, 120, 120, 0.6)",
                line_width=1.5,
                annotation_text=ref_line_label or f"{ref_line_value:.2f}",
                annotation_position=ann_position,
                annotation_font=dict(size=10, color="#666666", family=PLOTLY_FONT_FAMILY),
            )

        return fig
