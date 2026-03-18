"""QualityControl — QC checks, outlier detection, and triplicate statistics.

Provides Grubbs test, plate heatmap, replicate stats, and well-level diagnostics.
No Streamlit dependency — all methods are pure computation.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy import stats
from typing import Tuple


class QualityControl:
    CT_HIGH_THRESHOLD = 35.0
    CT_LOW_THRESHOLD = 10.0
    CV_THRESHOLD = 0.05
    HK_VARIATION_THRESHOLD = 1.0
    GRUBBS_ALPHA = 0.05

    @staticmethod
    def get_triplicate_data(
        data: pd.DataFrame, excluded_wells: set = None
    ) -> pd.DataFrame:
        """
        Build a comprehensive triplicate-level view of all CT values.
        Returns DataFrame with one row per well, grouped by Sample+Target.
        """
        if data is None or data.empty:
            return pd.DataFrame()

        excluded_wells = excluded_wells or set()

        df = data[["Well", "Sample", "Target", "CT"]].copy()
        df["Excluded"] = df["Well"].isin(excluded_wells)

        group_stats = (
            df[~df["Excluded"]]
            .groupby(["Sample", "Target"])
            .agg({"CT": ["mean", "std", "count", "min", "max", list], "Well": list})
            .reset_index()
        )

        group_stats.columns = [
            "Sample",
            "Target",
            "Mean_CT",
            "SD",
            "n",
            "Min_CT",
            "Max_CT",
            "CT_Values",
            "Wells",
        ]

        # FIX-17: Return NaN instead of 0 when mean CT is non-positive (invalid)
        group_stats["CV_pct"] = np.where(
            group_stats["Mean_CT"] > 0,
            (group_stats["SD"] / group_stats["Mean_CT"]) * 100,
            np.nan,
        )
        group_stats["Range"] = group_stats["Max_CT"] - group_stats["Min_CT"]

        def get_health_status(row):
            issues = []
            severity = "ok"

            if row["n"] < 2:
                issues.append("Low n")
                severity = "warning"
            if pd.notna(row["CV_pct"]) and row["CV_pct"] > QualityControl.CV_THRESHOLD * 100:
                issues.append(f"High CV ({row['CV_pct']:.1f}%)")
                severity = "warning"
            if row["Mean_CT"] > QualityControl.CT_HIGH_THRESHOLD:
                issues.append("High CT")
                severity = "warning"
            if row["Mean_CT"] < QualityControl.CT_LOW_THRESHOLD:
                issues.append("Low CT")
                severity = "warning"
            if row["Range"] > 1.0 and row["n"] >= 2:
                issues.append(f"High range ({row['Range']:.2f})")
                severity = "error" if row["Range"] > 2.0 else "warning"

            if row["n"] >= 3:
                ct_vals = np.array(row["CT_Values"])
                is_outlier, _ = QualityControl.grubbs_test(
                    ct_vals, QualityControl.GRUBBS_ALPHA
                )
                if is_outlier:
                    issues.append("Has outlier")
                    severity = "error"

            if not issues:
                return "OK", "ok"
            return "; ".join(issues), severity

        group_stats[["Status", "Severity"]] = group_stats.apply(
            lambda row: pd.Series(get_health_status(row)), axis=1
        )

        group_stats["Mean_CT"] = group_stats["Mean_CT"].round(2)
        group_stats["SD"] = group_stats["SD"].round(3)
        group_stats["CV_pct"] = group_stats["CV_pct"].round(1)
        group_stats["Range"] = group_stats["Range"].round(2)

        return group_stats

    @staticmethod
    def get_wells_for_triplicate(
        data: pd.DataFrame, sample: str, target: str
    ) -> pd.DataFrame:
        """
        Get all wells for a specific Sample+Target combination with detailed info.
        """
        if data is None or data.empty:
            return pd.DataFrame()

        wells = data[(data["Sample"] == sample) & (data["Target"] == target)].copy()

        if wells.empty:
            return pd.DataFrame()

        mean_ct = wells["CT"].mean()
        std_ct = wells["CT"].std() if len(wells) > 1 else 0

        wells["Is_Outlier"] = False
        if len(wells) >= 3:
            ct_vals = wells["CT"].values
            is_outlier, outlier_idx = QualityControl.grubbs_test(
                ct_vals, QualityControl.GRUBBS_ALPHA
            )
            if is_outlier and outlier_idx >= 0:
                wells.iloc[outlier_idx, wells.columns.get_loc("Is_Outlier")] = True

        wells["Deviation"] = (wells["CT"] - mean_ct).round(3)
        wells["Z_Score"] = np.where(
            std_ct > 0, (wells["CT"] - mean_ct) / std_ct, np.nan
        ).round(2)

        def well_status(row):
            issues = []
            if row["CT"] > QualityControl.CT_HIGH_THRESHOLD:
                issues.append("High CT")
            if row["CT"] < QualityControl.CT_LOW_THRESHOLD:
                issues.append("Low CT")
            if row["Is_Outlier"]:
                issues.append("Grubbs outlier")
            if pd.notna(row["Z_Score"]) and abs(row["Z_Score"]) > 2:
                issues.append(f"Z={row['Z_Score']:.1f}")
            return "; ".join(issues) if issues else "OK"

        wells["Well_Status"] = wells.apply(well_status, axis=1)

        return wells[
            [
                "Well",
                "Sample",
                "Target",
                "CT",
                "Deviation",
                "Z_Score",
                "Is_Outlier",
                "Well_Status",
            ]
        ]

    @staticmethod
    def get_qc_summary_stats(data: pd.DataFrame, excluded_wells: set = None) -> dict:
        """
        Calculate overall QC summary statistics for the dataset.
        """
        if data is None or data.empty:
            return {}

        excluded_wells = excluded_wells or set()
        active_data = data[~data["Well"].isin(excluded_wells)]

        total_wells = len(data)
        excluded_count = len(data[data["Well"].isin(excluded_wells)])
        active_wells = total_wells - excluded_count

        ct_mean = active_data["CT"].mean()
        ct_std = active_data["CT"].std()
        ct_min = active_data["CT"].min()
        ct_max = active_data["CT"].max()

        high_ct_count = len(
            active_data[active_data["CT"] > QualityControl.CT_HIGH_THRESHOLD]
        )
        low_ct_count = len(
            active_data[active_data["CT"] < QualityControl.CT_LOW_THRESHOLD]
        )

        triplicate_stats = QualityControl.get_triplicate_data(data, excluded_wells)
        if not triplicate_stats.empty:
            total_triplicates = len(triplicate_stats)
            healthy_triplicates = len(
                triplicate_stats[triplicate_stats["Status"] == "OK"]
            )
            warning_triplicates = len(
                triplicate_stats[triplicate_stats["Severity"] == "warning"]
            )
            error_triplicates = len(
                triplicate_stats[triplicate_stats["Severity"] == "error"]
            )
            avg_cv = triplicate_stats["CV_pct"].mean()
            max_cv = triplicate_stats["CV_pct"].max()
        else:
            total_triplicates = healthy_triplicates = warning_triplicates = (
                error_triplicates
            ) = 0
            avg_cv = max_cv = np.nan

        n_samples = active_data["Sample"].nunique()
        n_genes = active_data["Target"].nunique()

        return {
            "total_wells": total_wells,
            "excluded_wells": excluded_count,
            "active_wells": active_wells,
            "ct_mean": round(ct_mean, 2) if pd.notna(ct_mean) else 0,
            "ct_std": round(ct_std, 2) if pd.notna(ct_std) else 0,
            "ct_min": round(ct_min, 2) if pd.notna(ct_min) else 0,
            "ct_max": round(ct_max, 2) if pd.notna(ct_max) else 0,
            "high_ct_count": high_ct_count,
            "low_ct_count": low_ct_count,
            "total_triplicates": total_triplicates,
            "healthy_triplicates": healthy_triplicates,
            "warning_triplicates": warning_triplicates,
            "error_triplicates": error_triplicates,
            "avg_cv_pct": round(avg_cv, 1) if pd.notna(avg_cv) else 0,
            "max_cv_pct": round(max_cv, 1) if pd.notna(max_cv) else 0,
            "n_samples": n_samples,
            "n_genes": n_genes,
            "health_score": round(healthy_triplicates / total_triplicates * 100, 1)
            if total_triplicates > 0
            else 0,
        }

    @staticmethod
    def suggest_exclusions(
        data: pd.DataFrame,
        sample: str,
        target: str,
        excluded_wells: set = None,
        strategy: str = "outlier",
    ) -> list:
        """
        Suggest wells to exclude based on different strategies.

        Strategies:
        - 'outlier': Exclude statistical outliers (Grubbs test)
        - 'worst': Exclude the well with highest deviation from mean
        - 'keep_best_2': Keep the 2 closest values, exclude others
        """
        excluded_wells = excluded_wells or set()
        wells_df = QualityControl.get_wells_for_triplicate(data, sample, target)

        if wells_df.empty:
            return []

        active_wells = wells_df[~wells_df["Well"].isin(excluded_wells)]

        if len(active_wells) < 2:
            return []

        suggestions = []

        if strategy == "outlier":
            outliers = active_wells[active_wells["Is_Outlier"]]
            suggestions = outliers["Well"].tolist()

        elif strategy == "worst":
            if len(active_wells) > 2:
                worst_idx = active_wells["Deviation"].abs().idxmax()
                suggestions = [active_wells.loc[worst_idx, "Well"]]

        elif strategy == "keep_best_2":
            if len(active_wells) > 2:
                median_ct = active_wells["CT"].median()
                active_wells_sorted = active_wells.copy()
                active_wells_sorted["Dist_to_Median"] = abs(
                    active_wells_sorted["CT"] - median_ct
                )
                active_wells_sorted = active_wells_sorted.sort_values("Dist_to_Median")
                to_exclude = active_wells_sorted.iloc[2:]["Well"].tolist()
                suggestions = to_exclude

        return suggestions

    @staticmethod
    def grubbs_test(values: np.ndarray, alpha: float = 0.05) -> Tuple[bool, int]:
        n = len(values)
        if n < 3:
            return False, -1

        mean_val = np.mean(values)
        std_val = np.std(values, ddof=1)

        if std_val == 0:
            return False, -1

        g_scores = np.abs(values - mean_val) / std_val
        max_idx = np.argmax(g_scores)
        g_stat = g_scores[max_idx]

        t_crit = stats.t.ppf(1 - alpha / (2 * n), n - 2)
        g_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t_crit**2 / (n - 2 + t_crit**2))

        return g_stat > g_crit, int(max_idx)

    @staticmethod
    def detect_outliers(data: pd.DataFrame, hk_gene: str = None) -> pd.DataFrame:
        if data is None or data.empty:
            return pd.DataFrame()

        qc_df = data[["Well", "Sample", "Target", "CT"]].copy()

        ct_high = qc_df["CT"] > QualityControl.CT_HIGH_THRESHOLD
        ct_low = qc_df["CT"] < QualityControl.CT_LOW_THRESHOLD

        high_ct_issue = f"CT > {QualityControl.CT_HIGH_THRESHOLD} (low expression)"
        low_ct_issue = f"CT < {QualityControl.CT_LOW_THRESHOLD} (unusually high)"

        qc_df["Issues"] = pd.DataFrame(
            {
                "high": ct_high.map(lambda x: high_ct_issue if x else ""),
                "low": ct_low.map(lambda x: low_ct_issue if x else ""),
            }
        ).apply(lambda row: "; ".join([x for x in row if x]) or "OK", axis=1)

        qc_df["Severity"] = np.where(ct_high | ct_low, "warning", "ok")
        qc_df["Flagged"] = ct_high | ct_low

        cv_stats = (
            data.groupby(["Sample", "Target"])["CT"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        cv_stats["cv"] = np.where(
            (cv_stats["mean"] > 0) & (cv_stats["count"] > 1),
            cv_stats["std"] / cv_stats["mean"],
            np.nan,
        )
        high_cv_groups = cv_stats[cv_stats["cv"] > QualityControl.CV_THRESHOLD][
            ["Sample", "Target", "cv"]
        ]

        if not high_cv_groups.empty:
            qc_df = qc_df.merge(high_cv_groups, on=["Sample", "Target"], how="left")
            has_high_cv = qc_df["cv"].notna()

            cv_issue = qc_df["cv"].apply(
                lambda x: f"CV={x:.1%} (high variability)" if pd.notna(x) else ""
            )
            qc_df.loc[has_high_cv, "Issues"] = qc_df.loc[has_high_cv].apply(
                lambda r: f"{r['Issues']}; {cv_issue[r.name]}"
                if r["Issues"] != "OK"
                else cv_issue[r.name],
                axis=1,
            )
            qc_df.loc[has_high_cv, "Severity"] = "warning"
            qc_df.loc[has_high_cv, "Flagged"] = True
            qc_df = qc_df.drop(columns=["cv"])

        grubbs_outliers = set()
        for (sample, target), group in data.groupby(["Sample", "Target"]):
            if len(group) >= 3:
                ct_vals = group["CT"].values
                is_outlier, outlier_idx = QualityControl.grubbs_test(
                    ct_vals, QualityControl.GRUBBS_ALPHA
                )
                if is_outlier:
                    outlier_well = group.iloc[outlier_idx]["Well"]
                    grubbs_outliers.add(outlier_well)

        if grubbs_outliers:
            grubbs_mask = qc_df["Well"].isin(grubbs_outliers)

            def add_grubbs_issue(current_issue):
                grubbs_issue = "Grubbs outlier"
                return (
                    f"{current_issue}; {grubbs_issue}"
                    if current_issue != "OK"
                    else grubbs_issue
                )

            qc_df.loc[grubbs_mask, "Issues"] = qc_df.loc[grubbs_mask, "Issues"].apply(
                add_grubbs_issue
            )
            qc_df.loc[grubbs_mask, "Severity"] = "error"
            qc_df.loc[grubbs_mask, "Flagged"] = True

        if hk_gene:
            hk_data = data[data["Target"] == hk_gene]
            if not hk_data.empty:
                hk_by_sample = hk_data.groupby("Sample")["CT"].mean()
                overall_hk_mean = hk_by_sample.mean()

                deviations = (hk_by_sample - overall_hk_mean).abs()
                flagged_samples = deviations[
                    deviations > QualityControl.HK_VARIATION_THRESHOLD
                ]

                if not flagged_samples.empty:
                    deviation_map = flagged_samples.to_dict()
                    hk_mask = (qc_df["Target"] == hk_gene) & (
                        qc_df["Sample"].isin(flagged_samples.index)
                    )

                    def add_hk_issue(row):
                        dev = deviation_map.get(row["Sample"], 0)
                        hk_issue = f"HK deviation={dev:.2f}"
                        return (
                            f"{row['Issues']}; {hk_issue}"
                            if row["Issues"] != "OK"
                            else hk_issue
                        )

                    qc_df.loc[hk_mask, "Issues"] = qc_df.loc[hk_mask].apply(
                        add_hk_issue, axis=1
                    )
                    qc_df.loc[hk_mask, "Severity"] = "error"
                    qc_df.loc[hk_mask, "Flagged"] = True

        return qc_df

    @staticmethod
    def get_replicate_stats(data: pd.DataFrame) -> pd.DataFrame:
        if data is None or data.empty:
            return pd.DataFrame()

        rep_stats = (
            data.groupby(["Sample", "Target"])["CT"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        rep_stats.columns = ["Sample", "Target", "Mean CT", "SD", "n"]

        rep_stats["SD"] = rep_stats["SD"].fillna(0)
        rep_stats["CV%"] = np.where(
            rep_stats["Mean CT"] > 0, (rep_stats["SD"] / rep_stats["Mean CT"]) * 100, np.nan
        )

        rep_stats["Status"] = np.select(
            [rep_stats["Mean CT"] < 10, rep_stats["Mean CT"] > 35, rep_stats["CV%"] > 5],
            ["Check Signal", "Low Expression", "High CV"],
            default="OK",
        )

        rep_stats["Mean CT"] = rep_stats["Mean CT"].round(2)
        rep_stats["SD"] = rep_stats["SD"].round(3)
        rep_stats["CV%"] = rep_stats["CV%"].round(1)
        rep_stats["n"] = rep_stats["n"].astype(int)

        return rep_stats[["Sample", "Target", "n", "Mean CT", "SD", "CV%", "Status"]]

    @staticmethod
    def create_plate_heatmap(
        data: pd.DataFrame, value_col: str = "CT", excluded_wells: set = None
    ) -> go.Figure:
        if data is None or data.empty:
            return go.Figure()

        excluded_wells = excluded_wells or set()

        rows = list("ABCDEFGH")
        cols = list(range(1, 13))

        plate_values = np.full((8, 12), np.nan)
        plate_text = [["" for _ in range(12)] for _ in range(8)]
        plate_colors = [[0 for _ in range(12)] for _ in range(8)]

        for _, row in data.iterrows():
            well = row["Well"]
            if len(well) >= 2:
                well_row = well[0].upper()
                try:
                    well_col = int(well[1:])
                except ValueError:
                    continue

                if well_row in rows and 1 <= well_col <= 12:
                    r_idx = rows.index(well_row)
                    c_idx = well_col - 1

                    ct_val = row[value_col]
                    plate_values[r_idx, c_idx] = ct_val

                    sample_short = str(row["Sample"])[:10]
                    target_short = str(row["Target"])[:8]
                    excluded_marker = " [X]" if well in excluded_wells else ""
                    plate_text[r_idx][c_idx] = (
                        f"{well}{excluded_marker}<br>{sample_short}<br>{target_short}<br>CT: {ct_val:.1f}"
                    )

                    if well in excluded_wells:
                        plate_colors[r_idx][c_idx] = -1

        fig = go.Figure(
            data=go.Heatmap(
                z=plate_values,
                x=[str(c) for c in cols],
                y=rows,
                text=plate_text,
                hoverinfo="text",
                colorscale=[[0, "#2ecc71"], [0.5, "#f1c40f"], [1, "#e74c3c"]],
                zmin=15,
                zmax=40,
                colorbar=dict(title="CT Value"),
            )
        )

        for r_idx, row_letter in enumerate(rows):
            for c_idx, col_num in enumerate(cols):
                well_name = f"{row_letter}{col_num}"
                if well_name in excluded_wells:
                    fig.add_annotation(
                        x=str(col_num),
                        y=row_letter,
                        text="X",
                        showarrow=False,
                        font=dict(size=20, color="red"),
                        opacity=0.8,
                    )

        fig.update_layout(
            title="96-Well Plate Overview",
            xaxis=dict(title="Column", side="top", dtick=1),
            yaxis=dict(title="Row", autorange="reversed"),
            height=400,
            width=800,
        )

        return fig
