"""GraphGenerator — Plotly bar chart visualizations for qPCR results.

Creates per-gene relative expression bar charts with error bars,
significance annotations, and customizable styling.
"""

import textwrap

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from qpcr.constants import PLOTLY_FONT_FAMILY, CM_TO_PX


def _darken_hex(hex_color: str, factor: float = 0.3) -> str:
    """Darken a hex color by the given factor (0-1)."""
    hex_color = hex_color.lstrip("#")
    if len(hex_color) != 6:
        return "#666666"
    r, g, b = int(hex_color[:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    r = int(r * (1 - factor))
    g = int(g * (1 - factor))
    b = int(b * (1 - factor))
    return f"#{r:02x}{g:02x}{b:02x}"


class GraphGenerator:
    @staticmethod
    def _wrap_text(text: str, width: int = 15) -> str:
        """Wrap text for x-axis labels using <br> for Plotly compatibility"""
        wrapped = textwrap.fill(text, width=width)
        return wrapped.replace("\n", "<br>")

    @staticmethod
    def _auto_wrap_width(n_bars: int, fig_width_cm: float = 28) -> int:
        """Calculate optimal wrap width based on number of bars and figure width.

        Adapts like Excel cell wrapping: more bars -> narrower wrap.
        Minimum 10 chars to keep labels readable even with 20+ bars.
        """
        px_per_bar = (fig_width_cm * 37.8) / max(n_bars, 1)
        chars_per_bar = max(int(px_per_bar / 7), 10)
        return min(chars_per_bar, 30)

    @staticmethod
    def create_gene_graph(
        data: pd.DataFrame,
        gene: str,
        settings: dict,
        efficacy_config: dict = None,
        sample_order: list = None,
        per_sample_overrides: dict = None,
        display_gene_name: str = None,
        ref_line_value: float = None,
        ref_line_label: str = None,
        show_data_points: bool = False,
        replicate_data: pd.DataFrame = None,
        color_preset: str = None,
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

        # Import GRAPH_PRESETS only when needed (avoid circular import at module level)
        from qpcr.constants import GRAPH_PRESETS

        # ---- BAR COLORS ----
        _is_ref = [False] * n_bars
        if "Fold_Change" in gene_data_indexed.columns:
            for i, (_, r) in enumerate(gene_data_indexed.iterrows()):
                if abs(r.get("Fold_Change", 0) - 1.0) < 0.001:
                    _is_ref[i] = True

        bar_colors = ["#FFFFFF"] * n_bars
        if color_preset and color_preset != "Custom" and color_preset in GRAPH_PRESETS:
            _tone = GRAPH_PRESETS[color_preset]["color"]
            _ref_c = GRAPH_PRESETS[color_preset]["ref"]
            bar_colors = [_ref_c if _is_ref[i] else _tone for i in range(n_bars)]
        else:
            for i, (_, row) in enumerate(gene_data_indexed.iterrows()):
                condition = row["Condition"]
                custom_key = f"{gene}_{condition}"
                if custom_key in settings.get("bar_colors_per_sample", {}):
                    bar_colors[i] = settings["bar_colors_per_sample"][custom_key]

        fig = go.Figure()

        # Error bars - use fold-change domain asymmetric bars (Livak method)
        # Falls back to Ct-domain SD/SEM if FC columns are missing
        has_fc_errors = "FC_Error_Upper" in gene_data_indexed.columns and "FC_Error_Lower" in gene_data_indexed.columns
        error_bar_type = st.session_state.get("error_bar_type", "sem")

        if has_fc_errors:
            error_upper_array = gene_data_indexed["FC_Error_Upper"].fillna(0).values
            error_lower_array = gene_data_indexed["FC_Error_Lower"].fillna(0).values
        else:
            # Legacy fallback: Ct-domain SD/SEM (symmetric, less accurate)
            error_col = "SD" if error_bar_type == "sd" else "SEM"
            if error_col not in gene_data_indexed.columns:
                error_col = "SEM"
            error_upper_array = gene_data_indexed[error_col].fillna(0).values
            error_lower_array = gene_data_indexed[error_col].fillna(0).values

        gene_bar_settings = st.session_state.get(f"{gene}_bar_settings", {})

        show_error_global = settings.get("show_error", True)
        show_sig_global = settings.get("show_significance", True)

        error_visible_upper = []
        error_visible_lower = []

        for idx in range(n_bars):
            row = gene_data_indexed.iloc[idx]
            condition = row["Condition"]
            bar_key = f"{gene}_{condition}"

            bar_config = gene_bar_settings.get(
                bar_key, {"show_sig": True, "show_err": True}
            )

            if show_error_global and bar_config.get("show_err", True):
                error_visible_upper.append(error_upper_array[idx])
            else:
                error_visible_upper.append(0)
            error_visible_lower.append(0)  # top-only error bars

        fig.add_trace(
            go.Bar(
                x=list(range(n_bars)),
                y=gene_data_indexed["Relative_Expression"],
                error_y=dict(
                    type="data",
                    array=error_visible_upper,
                    arrayminus=error_visible_lower,
                    visible=True,
                    thickness=2,
                    width=6,
                    color="rgba(0,0,0,0.6)",
                    symmetric=False,
                ),
                marker=dict(
                    color=bar_colors,
                    line=dict(
                        width=settings.get("marker_line_width", 1), color="black"
                    ),
                    opacity=settings.get("bar_opacity", 0.85),
                ),
                showlegend=False,
            )
        )

        # Data point overlay (jittered scatter on top of bars)
        if show_data_points and replicate_data is not None and not replicate_data.empty:
            import hashlib
            scatter_x = []
            scatter_y = []
            scatter_colors = []
            for idx, condition in enumerate(condition_names):
                cond_replicates = replicate_data[replicate_data["Condition"] == condition]
                if cond_replicates.empty:
                    continue
                seed = int(hashlib.md5(f"{gene}_{condition}".encode()).hexdigest()[:8], 16)
                rng = np.random.RandomState(seed)
                n_pts = len(cond_replicates)
                jitter = rng.uniform(-0.15, 0.15, size=n_pts)
                scatter_x.extend([idx + j for j in jitter])
                scatter_y.extend(cond_replicates["Replicate_FC"].tolist())
                base_color = bar_colors[idx] if idx < len(bar_colors) else "#666666"
                scatter_colors.extend([_darken_hex(base_color, 0.3)] * n_pts)

            if scatter_x:
                fig.add_trace(go.Scatter(
                    x=scatter_x,
                    y=scatter_y,
                    mode="markers",
                    marker=dict(size=5, color=scatter_colors, opacity=0.65, line=dict(width=0)),
                    showlegend=False,
                    hoverinfo="y",
                ))

        max_y_value = gene_data_indexed["Relative_Expression"].max()
        if pd.isna(max_y_value) or max_y_value <= 0:
            max_y_value = 1.0  # Fallback for all-NaN or zero expression
        max_error = error_upper_array.max() if len(error_upper_array) > 0 else 0
        if pd.isna(max_error):
            max_error = 0
        y_max_auto = (
            max_y_value + max_error + (max_y_value * 0.15)
        )
        if ref_line_value is not None and pd.notna(ref_line_value):
            y_max_auto = max(y_max_auto, ref_line_value * 1.20)

        fixed_symbol_spacing = y_max_auto * 0.05

        if show_sig_global:
            # Direct mode: add significance symbols above each bar
            for idx in range(n_bars):
                row = gene_data_indexed.iloc[idx]
                condition = row["Condition"]
                bar_key = f"{gene}_{condition}"
                bar_config = gene_bar_settings.get(
                    bar_key, {"show_sig": True, "show_err": True}
                )

                # Get all significance values (3 comparisons)
                sig_1 = row.get("significance", "")
                sig_2 = row.get("significance_2", "")
                sig_3 = row.get("significance_3", "")

                bar_height = row["Relative_Expression"]
                error_bar_height = error_visible_upper[idx]
                base_y_position = bar_height + error_bar_height

                asterisk_font_size = 16
                hashtag_font_size = 10
                dagger_font_size = 10

                if show_sig_global:
                    symbols_to_show = []
                    font_sizes = []

                    if sig_1 in ["*", "**", "***"] and bar_config.get("show_sig_1", bar_config.get("show_sig", True)):
                        symbols_to_show.append(sig_1)
                        font_sizes.append(asterisk_font_size)

                    if sig_2 in ["#", "##", "###"] and bar_config.get("show_sig_2", bar_config.get("show_sig", True)):
                        symbols_to_show.append(sig_2)
                        font_sizes.append(hashtag_font_size)

                    if sig_3 in ["\u2020", "\u2020\u2020", "\u2020\u2020\u2020"] and bar_config.get("show_sig_3", bar_config.get("show_sig", True)):
                        symbols_to_show.append(sig_3)
                        font_sizes.append(dagger_font_size)

                    # Stack symbols with fixed absolute spacing
                    for si, (sym, fs) in enumerate(zip(symbols_to_show, font_sizes)):
                        y_pos = base_y_position + (fixed_symbol_spacing * 0.2) + (si * fixed_symbol_spacing)
                        fig.add_annotation(
                            x=idx,
                            y=y_pos,
                            text=sym,
                            showarrow=False,
                            font=dict(size=fs, color="black", family=PLOTLY_FONT_FAMILY),
                            xref="x",
                            yref="y",
                            xanchor="center",
                            yanchor="bottom",
                        )

        gene_label = display_gene_name if display_gene_name else gene
        y_label_html = f"<b>Relative <span style='color:red;'>{gene_label}</span> Expression Level</b>"

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
            linewidth=1.5,
            linecolor="#2C3E50",
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

        # Auto-scale figure width for many bars (minimum 1.4cm per bar)
        configured_width = settings.get("figure_width", 28)
        min_width_for_bars = n_bars * 1.4
        effective_fig_width = max(configured_width, min_width_for_bars)

        # Auto-reduce tick font for dense graphs (>12 bars)
        if n_bars > 12 and gene_tick_size > 9:
            gene_tick_size = max(9, gene_tick_size - (n_bars - 12) // 3)

        # X-axis label mode: Auto-wrap / Angled 45 / Angled 90 / Horizontal
        label_mode = settings.get("label_mode", "Auto-wrap")
        x_tick_angle = 0

        if label_mode == "Auto-wrap":
            wrap_w = GraphGenerator._auto_wrap_width(n_bars, effective_fig_width)
            wrapped_labels = [
                GraphGenerator._wrap_text(str(cond), wrap_w) for cond in condition_names
            ]
        elif label_mode == "Angled 45\u00b0":
            wrapped_labels = [str(c) for c in condition_names]
            x_tick_angle = -45
        elif label_mode == "Angled 90\u00b0":
            wrapped_labels = [str(c) for c in condition_names]
            x_tick_angle = -90
        else:  # Horizontal
            wrapped_labels = [str(c) for c in condition_names]

        # Dynamic bottom margin based on label mode and content
        if label_mode == "Auto-wrap":
            max_label_lines = max(
                (label.count("<br>") + 1 for label in wrapped_labels), default=1
            )
            dynamic_b_margin = 180 + max(0, max_label_lines - 1) * 22
        elif label_mode == "Angled 45\u00b0":
            max_label_len = max((len(str(c)) for c in condition_names), default=5)
            dynamic_b_margin = 140 + max_label_len * 4
            max_label_lines = max(1, dynamic_b_margin // 22)
        elif label_mode == "Angled 90\u00b0":
            max_label_len = max((len(str(c)) for c in condition_names), default=5)
            dynamic_b_margin = 120 + max_label_len * 6
            max_label_lines = max(1, dynamic_b_margin // 22)
        else:
            dynamic_b_margin = 140
            max_label_lines = 1

        default_margins = gene_margins.copy()
        if default_margins.get("b", 200) < dynamic_b_margin:
            default_margins["b"] = dynamic_b_margin
        gene_margins = default_margins

        # P-VALUE LEGEND - Support dual/triple comparison with reference names
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

        if (
            "significance_3" in gene_data_indexed.columns
            and gene_data_indexed["significance_3"].notna().any()
        ):
            cmp_ref_name_3 = st.session_state.get("analysis_cmp_condition_3", "")
            legend_ref_label_3 = f" (vs {cmp_ref_name_3})" if cmp_ref_name_3 else ""
            legend_text += (
                f"<br><b>3rd Comparison{legend_ref_label_3}:</b>  \u2020 p<0.05  \u2020\u2020 p<0.01  \u2020\u2020\u2020 p<0.001"
            )

        # Reserve extra bottom margin for the significance legend below labels
        legend_line_count = legend_text.count("<br>") + 1
        legend_extra_px = legend_line_count * 18 + 20

        fig_h_px = int(settings.get("figure_height", 16) * CM_TO_PX)
        b_margin_px = gene_margins.get("b", 200)
        t_margin_px = gene_margins.get("t", 60)
        plot_h_px = max(fig_h_px - b_margin_px - legend_extra_px - t_margin_px, 100)

        fig.update_layout(
            title="",
            xaxis=dict(
                title="",
                showgrid=False,
                zeroline=False,
                tickmode="array",
                tickvals=list(range(n_bars)),
                ticktext=wrapped_labels,
                tickfont=dict(size=gene_tick_size, family=PLOTLY_FONT_FAMILY, color="black"),
                tickangle=x_tick_angle,
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
            font=dict(size=settings.get("font_size", 14), family=PLOTLY_FONT_FAMILY, color="black"),
            height=fig_h_px,
            width=int(effective_fig_width * CM_TO_PX),
            bargap=gene_bar_gap,
            showlegend=settings.get("show_legend", False),
            plot_bgcolor=gene_bg_color,
            paper_bgcolor="#FFFFFF",
            margin=dict(
                l=gene_margins.get("l", 80),
                r=gene_margins.get("r", 40),
                t=t_margin_px,
                b=b_margin_px + legend_extra_px,
            ),
        )

        # Place significance legend below the plot, beneath x-axis labels
        label_px_est = max_label_lines * 18 + 10
        legend_y_frac = -((label_px_est + 8) / plot_h_px)
        fig.add_annotation(
            text=legend_text,
            xref="paper",
            yref="paper",
            x=1.0,
            y=legend_y_frac,
            xanchor="right",
            yanchor="top",
            showarrow=False,
            font=dict(size=9, color="#666666", family=PLOTLY_FONT_FAMILY),
            bgcolor="rgba(255,255,255,0.90)",
            bordercolor="#CCCCCC",
            borderwidth=1,
            borderpad=4,
        )

        if ref_line_value is not None and pd.notna(ref_line_value):
            max_annotation_y = 0
            for _idx in range(n_bars):
                _row = gene_data_indexed.iloc[_idx]
                _bar_h = _row["Relative_Expression"]
                _err_h = error_visible_upper[_idx] if _idx < len(error_visible_upper) else 0
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
