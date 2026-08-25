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

# Where a chart image ends up. PPTGenerator places it 9.11in wide on a 13.33in
# slide, which is what turns figure pixels into points the reader actually sees.
PPT_PLACEMENT_WIDTH_IN = 9.11
# Floor for the smallest chart text once placed on that slide.
MIN_SLIDE_PT = 10.0
# The smallest hardcoded font in this module (the captions and n= labels).
SMALLEST_CHART_FONT_PX = 9


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
    def _auto_wrap_width(n_bars: int, fig_width_cm: float = 28,
                         tick_px: float = 12) -> int:
        """Calculate optimal wrap width based on number of bars and figure width.

        Adapts like Excel cell wrapping: more bars -> narrower wrap.
        Minimum 10 chars to keep labels readable even with 20+ bars.

        ``tick_px`` is the RENDERED tick font size. The character width used to
        be a bare /7, which silently assumed the old 12px tick font; once the
        fonts were scaled for slide legibility the labels grew ~1.8x wider than
        the wrap allowed for and adjacent labels ran into each other
        ("Test article 1 ppmTest article 10 ppm"). 0.58 em is a reasonable mean
        advance for this font stack at these sizes.
        """
        px_per_bar = (fig_width_cm * 37.8) / max(n_bars, 1)
        char_px = max(float(tick_px) * 0.58, 4.0)
        chars_per_bar = max(int(px_per_bar / char_px), 10)
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
        ref_condition: str = None,
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
        # Reference condition — match by condition name (exact), FC fallback (single match only)
        _is_ref = [False] * n_bars
        if ref_condition:
            for i, (_, r) in enumerate(gene_data_indexed.iterrows()):
                if r.get("Condition") == ref_condition:
                    _is_ref[i] = True
        elif "Fold_Change" in gene_data_indexed.columns:
            _fc_matches = [
                i for i, (_, r) in enumerate(gene_data_indexed.iterrows())
                if abs(r.get("Fold_Change", 0) - 1.0) < 0.001
            ]
            if len(_fc_matches) == 1:
                _is_ref[_fc_matches[0]] = True

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

        # Error bars — selectable mode. Every mode plots in the FOLD-CHANGE
        # domain, because that is the axis: the 'sem' and 'sd' modes used to put
        # the raw Ct-domain SEM/SD onto a fold-change bar, which is dimensionally
        # meaningless and inverts the comparison (a bar at FC 4.59 was given its
        # Ct SEM of 0.23 while its true fold-domain interval is 0.79, so the tall
        # bar drew a TIGHTER whisker than a short one).
        # error_caption documents which is shown (reviewers require this to be
        # stated) and now also states that the spread is of the target
        # replicates only — the housekeeping SD is deliberately not propagated,
        # so calling it plain "Livak" overstated it.
        cols = set(gene_data_indexed.columns)
        error_bar_mode = (settings.get("error_bar_mode")
                          or st.session_state.get("error_bar_mode")
                          or "livak_sd")
        if error_bar_mode == "ci95" and {"FC_CI_Upper", "FC_CI_Lower"} <= cols:
            _eu, _el = "FC_CI_Upper", "FC_CI_Lower"
            error_caption = "Error bars: 95% CI (fold-change domain, target replicates)"
        elif error_bar_mode == "sem" and {"FC_SEM_Upper", "FC_SEM_Lower"} <= cols:
            _eu, _el = "FC_SEM_Upper", "FC_SEM_Lower"
            error_caption = "Error bars: ±SEM of target replicates (fold-change domain)"
        elif {"FC_Error_Upper", "FC_Error_Lower"} <= cols:
            # 'sd' resolves here too: the Livak bars ARE the SD, transformed.
            _eu, _el = "FC_Error_Upper", "FC_Error_Lower"
            error_caption = "Error bars: ±SD of target replicates (fold-change domain)"
        else:
            # Last resort for a frame without the transformed columns (e.g. an
            # older stored result). State the domain honestly.
            _col = "SD" if error_bar_mode == "sd" else "SEM"
            _eu = _el = _col if _col in cols else "SEM"
            error_caption = f"Error bars: ±{_eu} (Ct domain — not fold-change scaled)"
        # NaN is kept, NOT filled with 0. calculate_ddct stores undefined spread
        # (n=1: no estimable variance) as NaN; filling it with 0 drew a
        # zero-length bar, which reads as "measured, perfectly precise" — the
        # opposite of the truth. Plotly serialises NaN in an error array to
        # null and omits that point's bar, leaving the others intact.
        error_upper_array = gene_data_indexed[_eu].values
        error_lower_array = gene_data_indexed[_el].values

        # Every trace used to be hardcoded showlegend=False, so the "Legend"
        # toggle switched on a legend with nothing in it. Name the traces and let
        # the setting drive their visibility instead.
        show_legend = bool(settings.get("show_legend", False))

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
                # None, not 0: a 0 still draws the cap, so a bar the user
                # switched off looked like a measured zero spread.
                error_visible_upper.append(None)
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
                name="Relative expression",
                showlegend=show_legend,
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
                    name="Replicates",
                    showlegend=show_legend,
                    hoverinfo="y",
                ))

        # ---- Font scale for the output medium -------------------------------
        # Chart text is sized in FIGURE PIXELS but read on a slide. The PPT
        # writer places the image 9.11in wide, so on a 28cm (1058px) figure one
        # figure pixel is 0.620pt: the 9px captions landed at 5.6pt and the 12px
        # ticks at 7.4pt, while the slide's own chrome is 11-16pt. The least
        # important text on the slide was three times the size of the axis
        # labels.
        #
        # One scale is applied to every chart font so the SMALLEST text clears a
        # 10pt floor on the slide (decision: Min, 2026-08-24, item 7a). Scaling
        # rather than clamping each size individually is deliberate — clamping
        # would raise the 9px captions and the 12px ticks to the same value and
        # flatten the hierarchy. The size sliders keep working; they now set
        # relative size, which is what they were always really doing.
        _cfg_w = settings.get("figure_width", 28)
        _eff_w_cm = max(_cfg_w, n_bars * 1.4)
        _fig_w_px = max(int(_eff_w_cm * CM_TO_PX), 1)
        font_scale = max(
            1.0,
            (MIN_SLIDE_PT * _fig_w_px)
            / (PPT_PLACEMENT_WIDTH_IN * 72.0 * SMALLEST_CHART_FONT_PX),
        )

        def _fs(px) -> int:
            """Scale a figure-pixel font size for the output medium.

            Rounds UP: rounding 9 * 1.792 = 16.13 down to 16 px lands at 9.92pt
            and misses the floor by 0.08pt, which would make the guarantee a
            near-miss rather than a guarantee.
            """
            return max(1, int(np.ceil(float(px) * font_scale - 1e-9)))

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

        # Stacked significance symbols are POSITIONED LATER, once plot_h_px is
        # known (see "emit the deferred significance symbols" below). The step
        # used to be y_max_auto * 0.05, a fraction of the DATA range, which at
        # the default 28x16cm geometry is ~13.5 px — smaller than the 16 px
        # asterisk it has to clear, so a dagger's stem was drawn straight
        # through the middle asterisk of a "***" and corrupted it. A glyph must
        # be cleared in PIXELS, so the step has to come from the font sizes.
        pending_sig = []

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
                # None (bar switched off) and NaN (undefined spread at n=1) both
                # mean "no error bar is drawn", so the marker sits at the bar
                # top. Only the LAYOUT coerces to 0 — the plotted array keeps
                # None/NaN so Plotly omits those bars.
                _err = error_visible_upper[idx]
                error_bar_height = (
                    0.0 if _err is None or pd.isna(_err) else float(_err)
                )
                base_y_position = bar_height + error_bar_height

                asterisk_font_size = _fs(16)
                hashtag_font_size = _fs(10)
                dagger_font_size = _fs(10)

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

                    if symbols_to_show:
                        pending_sig.append(
                            (idx, base_y_position,
                             list(zip(symbols_to_show, font_sizes)))
                        )

        gene_label = display_gene_name if display_gene_name else gene
        y_label_html = f"<b>Relative <span style='color:red;'>{gene_label}</span> Expression Level</b>"

        y_axis_config = dict(
            title=dict(
                text=y_label_html,
                font=dict(size=_fs(settings.get(f"{gene}_ylabel_size", 14)),
                          family=PLOTLY_FONT_FAMILY),
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

        is_log = bool(settings.get("y_log_scale"))
        if is_log:
            y_axis_config["type"] = "log"
            y_axis_config.pop("range", None)

        if settings.get("y_min") is not None or settings.get("y_max") is not None:
            y_lo = settings.get("y_min", 0)
            y_hi = settings.get("y_max", y_max_auto)
            if is_log:
                # Plotly reads a log axis's range as EXPONENTS, so passing the raw
                # bounds through turned a 1-100 request into 10^1-10^100 and drew a
                # blank chart. Convert, and fall back to auto when the bounds are
                # not expressible on a log axis rather than showing a wrong one.
                if y_lo is not None and y_lo > 0 and y_hi is not None and y_hi > y_lo:
                    y_axis_config["range"] = [np.log10(y_lo), np.log10(y_hi)]
            else:
                y_axis_config["range"] = [y_lo, y_hi]

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
            wrap_w = GraphGenerator._auto_wrap_width(
                n_bars, effective_fig_width, _fs(gene_tick_size)
            )
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

        # Bottom margin, measured from the LABEL BLOCK rather than a floor.
        # The old floors (180 / 140+len*4 / 120+len*6) were set independently of
        # the text they had to hold, and they over-reserved badly: on a default
        # 28x16cm figure the bottom margin resolved to 238px plus a 38px legend
        # reserve on a 604px figure, leaving a 268px plot area \u2014 44% of the
        # image, with 40% of every exported PNG blank (measured: ink rows
        # 117-849 of 1208, so 9.7% blank on top and 29.6% on the bottom). On the
        # slide that became a 1.54in empty band between the chart and the
        # "Results: \ud6a8\ub2a5" line.
        #
        # Now: one line of ticks costs one line of ticks. Measured after the
        # change on a 5-condition default figure: bottom blank 29.6% -> 0.0% for
        # Auto-wrap and Horizontal, and the plot area went from 44% to 75% of
        # the image.
        #
        # The angled modes need the label's PROJECTED height, so they still
        # scale with length: sin(45 deg) ~= 0.71 of an estimated glyph advance,
        # and the full advance at 90 deg. Those two keep ~9% bottom slack
        # (down from 12-17%) because the real projected height depends on glyph
        # metrics this code cannot measure without rendering first. That slack
        # is deliberate \u2014 under-reserving CLIPS the labels, which is worse than
        # white space. If it ever matters, xaxis automargin would let Plotly
        # measure it properly, at the cost of the caption anchoring below, which
        # is expressed as a fraction of the plot height.
        _line_px = _fs(gene_tick_size) * 1.45  # tick line height incl. leading
        _pad_px = _fs(10)                      # breathing room under the axis
        if label_mode == "Auto-wrap":
            max_label_lines = max(
                (label.count("<br>") + 1 for label in wrapped_labels), default=1
            )
            dynamic_b_margin = int(max_label_lines * _line_px + _pad_px)
        elif label_mode == "Angled 45\u00b0":
            max_label_len = max((len(str(c)) for c in condition_names), default=5)
            _glyph_px = _fs(gene_tick_size) * 0.50
            dynamic_b_margin = int(max_label_len * _glyph_px * 0.71 + _pad_px)
            max_label_lines = 1
        elif label_mode == "Angled 90\u00b0":
            max_label_len = max((len(str(c)) for c in condition_names), default=5)
            _glyph_px = _fs(gene_tick_size) * 0.50
            dynamic_b_margin = int(max_label_len * _glyph_px + _pad_px)
            max_label_lines = 1
        else:  # Horizontal
            dynamic_b_margin = int(_line_px + _pad_px)
            max_label_lines = 1

        # Use the measured value, do not merely floor it. The old code only ever
        # RAISED b towards dynamic_b_margin, so the {"b": 200} default won
        # whenever the labels needed less — which is every ordinary chart. An
        # explicit per-gene margin from the operator still wins over both.
        _has_explicit_margins = settings.get(f"{gene}_margins") is not None
        default_margins = gene_margins.copy()
        if not _has_explicit_margins:
            default_margins["b"] = dynamic_b_margin
            # No chart title is drawn (the gene name lives on the y-axis), so
            # the 60px top reserve was holding space for nothing. The
            # significance stack sits INSIDE the plot area — the y-range is
            # extended to contain it — so it does not need margin either.
            default_margins["t"] = int(_pad_px * 1.6)
        elif default_margins.get("b", 200) < dynamic_b_margin:
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
        # Scaled with the legend's own font, which the 18px constant predated.
        legend_extra_px = int(legend_line_count * _fs(9) * 1.5 + _pad_px)

        fig_h_px = int(settings.get("figure_height", 16) * CM_TO_PX)
        b_margin_px = gene_margins.get("b", 200)
        t_margin_px = gene_margins.get("t", 60)
        plot_h_px = max(fig_h_px - b_margin_px - legend_extra_px - t_margin_px, 100)

        # Emit the deferred significance symbols. The step between stacked
        # symbols is derived from the FONT SIZES in px and converted to data
        # units through the plot height, so each glyph clears the one below it
        # whatever the figure geometry or the data range. 1.25x the larger of the
        # two adjacent sizes gives a normal line's worth of leading.
        # Top of the tallest significance stack, in DATA units. Initialised
        # outside the block because the reference-line placement below reads it
        # whether or not any symbols were drawn.
        _stack_top = 0.0
        if pending_sig:
            _data_per_px = (y_max_auto / plot_h_px) if plot_h_px else 0.0
            for _idx, _base_y, _syms in pending_sig:
                _y = _base_y + max(_syms[0][1], 8) * 0.35 * _data_per_px
                _prev_fs = None
                _sym_fs = _syms[0][1]
                for _sym, _sym_fs in _syms:
                    if _prev_fs is not None:
                        _y += max(_prev_fs, _sym_fs) * 1.25 * _data_per_px
                    fig.add_annotation(
                        x=_idx,
                        y=_y,
                        text=_sym,
                        showarrow=False,
                        font=dict(size=_sym_fs, color="black",
                                  family=PLOTLY_FONT_FAMILY),
                        xref="x",
                        yref="y",
                        xanchor="center",
                        yanchor="bottom",
                    )
                    _prev_fs = _sym_fs
                # The glyph itself occupies height above its anchor.
                _stack_top = max(_stack_top, _y + _sym_fs * 1.2 * _data_per_px)

            # Give the stack room. Without this a taller stack is clipped by the
            # axis range, which was computed before the symbols were placed.
            # Only when the operator has not pinned the bounds themselves.
            if (settings.get("y_min") is None and settings.get("y_max") is None
                    and not is_log and _stack_top > y_max_auto):
                y_max_auto = _stack_top
                y_axis_config["range"] = [0, y_max_auto]

        fig.update_layout(
            title="",
            xaxis=dict(
                title="",
                showgrid=False,
                zeroline=False,
                tickmode="array",
                tickvals=list(range(n_bars)),
                ticktext=wrapped_labels,
                tickfont=dict(size=_fs(gene_tick_size), family=PLOTLY_FONT_FAMILY,
                              color="black"),
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
            font=dict(size=_fs(settings.get("font_size", 14)),
                      family=PLOTLY_FONT_FAMILY, color="black"),
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

        # Place significance legend below the plot, beneath x-axis labels.
        # label_px_est used a hardcoded 18px line height that predated the font
        # scaling, so once the ticks grew the captions were anchored INSIDE the
        # label block and the second line of a wrapped label ("...100 / ppm")
        # printed straight through them. _line_px is the same value the bottom
        # margin is computed from, so the two cannot drift apart again.
        label_px_est = max_label_lines * _line_px + _pad_px
        legend_y_frac = -((label_px_est + _pad_px * 0.5) / plot_h_px)
        fig.add_annotation(
            text=legend_text,
            xref="paper",
            yref="paper",
            x=1.0,
            y=legend_y_frac,
            xanchor="right",
            yanchor="top",
            showarrow=False,
            font=dict(size=_fs(9), color="#666666", family=PLOTLY_FONT_FAMILY),
            bgcolor="rgba(255,255,255,0.90)",
            bordercolor="#CCCCCC",
            borderwidth=1,
            borderpad=4,
        )

        if ref_line_value is not None and pd.notna(ref_line_value):
            # Tallest drawn element, in DATA units, so the label can sit clear of
            # it. This used to add `fixed_symbol_spacing * 1.2` — a name that has
            # not existed since the significance glyphs moved to `pending_sig`,
            # so EVERY chart with a reference line raised NameError. The memo
            # swallowed it (`_memoized_gene_figure`) and the gene then never
            # reached `st.session_state["graphs"]`, which is what the PPT and
            # image writers read: the gene vanished from the deck behind a
            # transient "chart unavailable" warning. `_stack_top` is the real
            # symbol-stack top computed above, already in these units.
            max_annotation_y = _stack_top
            for _idx in range(n_bars):
                _row = gene_data_indexed.iloc[_idx]
                _bar_h = _row["Relative_Expression"]
                _err_h = error_visible_upper[_idx] if _idx < len(error_visible_upper) else 0
                _top_y = _bar_h + _err_h
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
                annotation_font=dict(size=_fs(10), color="#666666",
                                     family=PLOTLY_FONT_FAMILY),
            )

        # Optional per-bar sample size (n=) — publication convention
        if settings.get("show_n") and "n_replicates" in gene_data_indexed.columns:
            for _i in range(n_bars):
                _n = gene_data_indexed.iloc[_i].get("n_replicates")
                if pd.notna(_n):
                    fig.add_annotation(
                        x=_i, y=0, yref="y", text=f"n={int(_n)}",
                        showarrow=False, yshift=9,
                        font=dict(size=_fs(9), color="#555555",
                                  family=PLOTLY_FONT_FAMILY),
                    )

        # Error-bar caption — always state what the error bars represent
        if settings.get("show_error", True):
            fig.add_annotation(
                text=error_caption, xref="paper", yref="paper",
                x=0.0, y=legend_y_frac, xanchor="left", yanchor="top",
                showarrow=False,
                # #666666 to match the significance legend beside it. #888888
                # measured 3.54:1 on white, below the 4.5:1 AA floor, while its
                # immediate neighbour sat at 5.74:1 — two adjacent captions in
                # two different greys, one of them failing.
                font=dict(size=_fs(9), color="#666666", family=PLOTLY_FONT_FAMILY),
            )

        return fig
