"""AnalysisEngine — DDCt calculations and statistical tests.

Contains calculate_ddct, calculate_statistics, and FDR correction.
run_full_analysis remains in the main Streamlit app (deep session_state coupling).
"""

import warnings

import numpy as np
import pandas as pd
import streamlit as st
from scipy import stats


class AnalysisEngine:
    @staticmethod
    def calculate_ddct(
        data: pd.DataFrame,
        hk_gene: str,
        ref_sample: str,
        excluded_wells,
        excluded_samples: set,
        sample_mapping: dict,
    ) -> pd.DataFrame:
        """Gene-by-gene DDCt calculation with housekeeping normalization.

        Args:
            excluded_wells: Either a dict {(gene, sample): set_of_wells} for
                per-gene-sample exclusions, or a flat set of well IDs for
                backward compatibility.
        """

        # Normalize excluded_wells: support both dict (per-gene-sample) and flat set
        if isinstance(excluded_wells, dict):
            excluded_wells_dict = excluded_wells
        else:
            excluded_wells_dict = None
            data = data[~data["Well"].isin(excluded_wells or set())].copy()

        # Filter excluded samples
        data = data[~data["Sample"].isin(excluded_samples)].copy()

        # Apply sample name mapping
        data["Condition"] = data["Sample"].map(
            lambda x: sample_mapping.get(x, {}).get("condition", x)
        )
        data["Group"] = data["Sample"].map(
            lambda x: sample_mapping.get(x, {}).get("group", "Treatment")
        )

        results = []
        # FIX-06: Accumulate warnings for skipped genes/conditions
        _skipped_warnings = []

        # Process each target gene separately (exclude housekeeping)
        for target in data["Target"].unique():
            if target.upper() == hk_gene.upper():
                continue

            target_data = data[data["Target"] == target]

            for condition in target_data["Condition"].unique():
                cond_data = target_data[target_data["Condition"] == condition]

                if len(cond_data) == 0:
                    continue

                # Per-gene-sample well exclusion: filter target wells
                if excluded_wells_dict is not None:
                    sample_name = cond_data["Sample"].iloc[0]
                    target_excluded = excluded_wells_dict.get((target, sample_name), set())
                    cond_data = cond_data[~cond_data["Well"].isin(target_excluded)]
                    if len(cond_data) == 0:
                        _skipped_warnings.append(
                            f"Gene '{target}', condition '{condition}': all wells excluded"
                        )
                        continue

                hk_data = data[
                    (data["Condition"] == condition) & (data["Target"] == hk_gene)
                ]

                # Per-gene-sample well exclusion: filter HK wells using HK gene key
                if excluded_wells_dict is not None:
                    hk_sample_name = hk_data["Sample"].iloc[0] if len(hk_data) > 0 else None
                    if hk_sample_name:
                        hk_excluded = excluded_wells_dict.get((hk_gene, hk_sample_name), set())
                        hk_data = hk_data[~hk_data["Well"].isin(hk_excluded)]

                if len(hk_data) == 0:
                    _skipped_warnings.append(
                        f"Gene '{target}', condition '{condition}': no HK gene ('{hk_gene}') data"
                    )
                    continue

                target_ct_values = cond_data["CT"].values
                hk_ct_values = hk_data["CT"].values

                if len(target_ct_values) == 0 or len(hk_ct_values) == 0:
                    continue

                target_ct_mean = target_ct_values.mean()
                hk_ct_mean = hk_ct_values.mean()

                if np.isnan(target_ct_mean) or np.isnan(hk_ct_mean):
                    _skipped_warnings.append(
                        f"Gene '{target}', condition '{condition}': NaN CT mean (target={target_ct_mean}, HK={hk_ct_mean})"
                    )
                    continue

                delta_ct = target_ct_mean - hk_ct_mean

                # Get reference DCt (ref_sample) — must also respect exclusions
                ref_target = target_data[target_data["Condition"] == ref_sample]
                ref_hk = data[
                    (data["Condition"] == ref_sample) & (data["Target"] == hk_gene)
                ]

                # Apply per-gene-sample exclusion to reference wells too
                if excluded_wells_dict is not None and len(ref_target) > 0:
                    ref_sample_name = ref_target["Sample"].iloc[0]
                    ref_target_excluded = excluded_wells_dict.get((target, ref_sample_name), set())
                    ref_target = ref_target[~ref_target["Well"].isin(ref_target_excluded)]

                if excluded_wells_dict is not None and len(ref_hk) > 0:
                    ref_hk_sample_name = ref_hk["Sample"].iloc[0]
                    ref_hk_excluded = excluded_wells_dict.get((hk_gene, ref_hk_sample_name), set())
                    ref_hk = ref_hk[~ref_hk["Well"].isin(ref_hk_excluded)]

                if len(ref_target) > 0 and len(ref_hk) > 0:
                    ref_delta_ct = ref_target["CT"].mean() - ref_hk["CT"].mean()
                else:
                    if condition == ref_sample:
                        ref_delta_ct = delta_ct
                    else:
                        continue

                ddct = delta_ct - ref_delta_ct
                # Guard against overflow: clamp ddct to prevent inf/0 in 2^(-ddct)
                ddct_clamped = np.clip(ddct, -50, 50)
                rel_expr = 2 ** (-ddct_clamped)

                target_sd = target_ct_values.std() if len(target_ct_values) > 1 else 0
                hk_sd = hk_ct_values.std() if len(hk_ct_values) > 1 else 0
                n_target = len(target_ct_values)
                n_hk = len(hk_ct_values)

                target_sem = target_sd / np.sqrt(n_target) if n_target > 1 else 0

                sd = target_sd
                sem = target_sem

                original_sample = cond_data["Sample"].iloc[0]
                group = sample_mapping.get(original_sample, {}).get(
                    "group", "Treatment"
                )

                # FIX-17: Return NaN instead of 0 when mean CT is non-positive
                target_cv = (
                    (target_sd / target_ct_mean * 100) if target_ct_mean > 0 else np.nan
                )
                hk_cv = (hk_sd / hk_ct_mean * 100) if hk_ct_mean > 0 else np.nan

                results.append(
                    {
                        "Target": target,
                        "Condition": condition,
                        "Original_Sample": original_sample,
                        "Group": group,
                        "n_replicates": n_target,
                        "n_hk_replicates": n_hk,
                        "Target_Ct_Mean": target_ct_mean,
                        "Target_Ct_SD": target_sd,
                        "Target_Ct_CV%": target_cv,
                        "HK_Ct_Mean": hk_ct_mean,
                        "HK_Ct_SD": hk_sd,
                        "HK_Ct_CV%": hk_cv,
                        "Delta_Ct": delta_ct,
                        "Delta_Delta_Ct": ddct,
                        "Relative_Expression": rel_expr,
                        "SEM": sem,
                        "SD": sd,
                        "Fold_Change": rel_expr,
                    }
                )

        result_df = pd.DataFrame(results)
        # FIX-06: Attach warnings as attribute for caller to display
        result_df.attrs["_skipped_warnings"] = _skipped_warnings
        return result_df

    @staticmethod
    def calculate_statistics(
        processed: pd.DataFrame,
        compare_condition: str,
        compare_condition_2: str = None,
        raw_data: pd.DataFrame = None,
        hk_gene: str = None,
        sample_mapping: dict = None,
        ttest_type: str = "welch",
        excluded_wells=None,
    ) -> pd.DataFrame:
        """T-test comparing each condition to compare_condition (and optionally compare_condition_2).

        Args:
            ttest_type: "welch" for Welch's t-test (unequal variance), "student" for Student's t-test (equal variance)
            excluded_wells: Either a dict {(gene, sample): set_of_wells} for
                per-gene-sample exclusions, or a flat set of well IDs, or None.
        """

        # Use session_state fallbacks
        raw_data = raw_data if raw_data is not None else st.session_state.get("data")
        hk_gene = hk_gene if hk_gene is not None else st.session_state.get("hk_gene")
        sample_mapping = (
            sample_mapping
            if sample_mapping is not None
            else st.session_state.get("sample_mapping", {})
        )

        # Normalize excluded_wells
        if isinstance(excluded_wells, dict):
            excluded_wells_dict = excluded_wells
        else:
            excluded_wells_dict = None
            if excluded_wells and raw_data is not None:
                raw_data = raw_data[~raw_data["Well"].isin(excluded_wells)].copy()

        if raw_data is None or hk_gene is None:
            return processed

        results = processed.copy()
        results["p_value"] = np.nan
        results["significance"] = ""
        # FIX-16: Track one-sample t-test usage
        _onesamp_warnings = []
        _stats_skipped = []

        if compare_condition_2:
            results["p_value_2"] = np.nan
            results["significance_2"] = ""

        for target in results["Target"].unique():
            if pd.isna(target):
                continue

            t_rows = raw_data[raw_data["Target"] == target].copy()
            if t_rows.empty:
                continue

            # Per-gene-sample well exclusion for target gene
            if excluded_wells_dict is not None:
                exclude_mask = t_rows.apply(
                    lambda r: r["Well"] in excluded_wells_dict.get((target, r["Sample"]), set()),
                    axis=1,
                )
                t_rows = t_rows[~exclude_mask]
                if t_rows.empty:
                    continue

            t_rows["Condition"] = t_rows["Sample"].map(
                lambda s: sample_mapping.get(s, {}).get("condition", s)
            )

            hk_rows = raw_data[raw_data["Target"] == hk_gene].copy()

            # Per-gene-sample well exclusion for HK gene
            if excluded_wells_dict is not None:
                hk_exclude_mask = hk_rows.apply(
                    lambda r: r["Well"] in excluded_wells_dict.get((hk_gene, r["Sample"]), set()),
                    axis=1,
                )
                hk_rows = hk_rows[~hk_exclude_mask]

            hk_rows["Condition"] = hk_rows["Sample"].map(
                lambda s: sample_mapping.get(s, {}).get("condition", s)
            )
            hk_means = hk_rows.groupby("Condition")["CT"].mean().to_dict()

            rel_expr = {}
            for cond, grp in t_rows.groupby("Condition"):
                hk_mean = hk_means.get(cond, np.nan)
                if np.isnan(hk_mean) or hk_mean <= 0:
                    _stats_skipped.append(f"{target}/{cond}: invalid HK mean ({hk_mean})")
                    continue
                rel_expr[cond] = 2 ** (-(grp["CT"].values - hk_mean))

            # FIRST COMPARISON: compare_condition
            ref_vals = rel_expr.get(compare_condition, np.array([]))
            if ref_vals.size >= 1:
                for cond, vals in rel_expr.items():
                    if cond == compare_condition or vals.size == 0:
                        continue

                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", RuntimeWarning)
                            equal_var = ttest_type == "student"
                            if ref_vals.size >= 2 and vals.size >= 2:
                                _, p_val = stats.ttest_ind(
                                    ref_vals, vals, equal_var=equal_var
                                )
                            elif vals.size == 1 and ref_vals.size >= 2:
                                _onesamp_warnings.append(f"{target}/{cond} (n=1 vs ref n={ref_vals.size})")
                                _, p_val = stats.ttest_1samp(ref_vals, vals[0])
                            elif ref_vals.size == 1 and vals.size >= 2:
                                _onesamp_warnings.append(f"{target}/{cond} (ref n=1 vs n={vals.size})")
                                _, p_val = stats.ttest_1samp(vals, ref_vals[0])
                            else:
                                p_val = np.nan
                    except (ValueError, TypeError):
                        p_val = np.nan

                    mask = (results["Target"] == target) & (
                        results["Condition"] == cond
                    )
                    results.loc[mask, "p_value"] = p_val

                    if not np.isnan(p_val):
                        if p_val < 0.001:
                            results.loc[mask, "significance"] = "***"
                        elif p_val < 0.01:
                            results.loc[mask, "significance"] = "**"
                        elif p_val < 0.05:
                            results.loc[mask, "significance"] = "*"

            # SECOND COMPARISON: compare_condition_2 (if provided)
            if compare_condition_2:
                ref_vals_2 = rel_expr.get(compare_condition_2, np.array([]))
                if ref_vals_2.size >= 1:
                    for cond, vals in rel_expr.items():
                        if cond == compare_condition_2 or vals.size == 0:
                            continue

                        try:
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", RuntimeWarning)
                                equal_var = ttest_type == "student"
                                if ref_vals_2.size >= 2 and vals.size >= 2:
                                    _, p_val_2 = stats.ttest_ind(
                                        ref_vals_2, vals, equal_var=equal_var
                                    )
                                elif vals.size == 1 and ref_vals_2.size >= 2:
                                    _, p_val_2 = stats.ttest_1samp(ref_vals_2, vals[0])
                                elif ref_vals_2.size == 1 and vals.size >= 2:
                                    _, p_val_2 = stats.ttest_1samp(vals, ref_vals_2[0])
                                else:
                                    p_val_2 = np.nan
                        except (ValueError, TypeError):
                            p_val_2 = np.nan

                        mask = (results["Target"] == target) & (
                            results["Condition"] == cond
                        )
                        results.loc[mask, "p_value_2"] = p_val_2

                        if not np.isnan(p_val_2):
                            if p_val_2 < 0.001:
                                results.loc[mask, "significance_2"] = "###"
                            elif p_val_2 < 0.01:
                                results.loc[mask, "significance_2"] = "##"
                            elif p_val_2 < 0.05:
                                results.loc[mask, "significance_2"] = "#"

        results = AnalysisEngine._apply_fdr_correction(
            results, "p_value", "p_value_fdr", "significance_fdr", "*"
        )
        if compare_condition_2:
            results = AnalysisEngine._apply_fdr_correction(
                results, "p_value_2", "p_value_fdr_2", "significance_fdr_2", "#"
            )

        # FIX-16: Attach one-sample t-test warnings for caller to display
        results.attrs["_onesamp_warnings"] = _onesamp_warnings
        results.attrs["_stats_skipped_warnings"] = _stats_skipped
        return results

    @staticmethod
    def _apply_fdr_correction(
        results: pd.DataFrame, p_col: str, fdr_col: str, sig_col: str, marker: str
    ) -> pd.DataFrame:
        valid_mask = results[p_col].notna()
        valid_pvals = results.loc[valid_mask, p_col].values

        if len(valid_pvals) == 0:
            results[fdr_col] = np.nan
            results[sig_col] = ""
            return results

        n = len(valid_pvals)
        sorted_idx = np.argsort(valid_pvals)
        sorted_pvals = valid_pvals[sorted_idx]

        fdr_vals = np.zeros(n)
        for i in range(n):
            rank = i + 1
            fdr_vals[i] = sorted_pvals[i] * n / rank

        for i in range(n - 2, -1, -1):
            fdr_vals[i] = min(fdr_vals[i], fdr_vals[i + 1])

        fdr_vals = np.clip(fdr_vals, 0, 1)

        unsorted_fdr = np.zeros(n)
        unsorted_fdr[sorted_idx] = fdr_vals

        results[fdr_col] = np.nan
        results.loc[valid_mask, fdr_col] = unsorted_fdr

        results[sig_col] = ""
        for i, (idx, fdr) in enumerate(zip(results.index[valid_mask], unsorted_fdr)):
            if fdr < 0.001:
                results.loc[idx, sig_col] = marker * 3
            elif fdr < 0.01:
                results.loc[idx, sig_col] = marker * 2
            elif fdr < 0.05:
                results.loc[idx, sig_col] = marker

        return results
