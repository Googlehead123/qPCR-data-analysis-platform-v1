"""QualityControl — QC checks, outlier detection, and triplicate statistics.

Provides Grubbs test, replicate stats, and well-level diagnostics.
No Streamlit dependency — all methods are pure computation.
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import Tuple


class QualityControl:
    # IMMUTABLE DEFAULTS. Do NOT reassign these at runtime.
    #
    # The QC Settings widgets used to write straight onto these class
    # attributes. A class attribute is PROCESS state while a Streamlit session
    # is not: one container runs one Python process and every browser session
    # shares module and class state, so one user's thresholds silently became
    # every concurrent user's thresholds — and persisted across a refresh, so a
    # fresh session inherited them instead of these documented defaults. On a
    # validated instrument tool that is also a provenance problem:
    # build_provenance would record thresholds the operator never set.
    #
    # They feed QC DISPLAY and the provenance record, not exclusion —
    # apply_auto_qc passes its own sd_threshold explicitly to
    # auto_select_replicates — so no published fold change moved. That is what
    # makes this Major rather than Critical.
    #
    # Callers now pass a per-session `thresholds` dict; these remain the
    # fallback. See resolve_thresholds.
    CT_HIGH_THRESHOLD = 35.0
    CT_LOW_THRESHOLD = 10.0
    CV_THRESHOLD = 0.05
    HK_VARIATION_THRESHOLD = 1.0

    THRESHOLD_KEYS = ("ct_high", "ct_low", "cv", "hk_variation")

    GRUBBS_ALPHA = 0.05

    @staticmethod
    def resolve_thresholds(thresholds: dict = None) -> dict:
        """Fill a partial per-session threshold dict from the class defaults.

        Accepts None so every method keeps working for callers that do not care
        (the tests, and any pure-computation use), while the app threads the
        session's own values through.
        """
        base = {
            "ct_high": QualityControl.CT_HIGH_THRESHOLD,
            "ct_low": QualityControl.CT_LOW_THRESHOLD,
            "cv": QualityControl.CV_THRESHOLD,
            "hk_variation": QualityControl.HK_VARIATION_THRESHOLD,
        }
        if not thresholds:
            return base
        for key in QualityControl.THRESHOLD_KEYS:
            value = thresholds.get(key)
            if value is not None:
                try:
                    base[key] = float(value)
                except (TypeError, ValueError):
                    pass
        return base

    @staticmethod
    def excluded_mask(df: pd.DataFrame, excluded_wells) -> pd.Series:
        """Boolean mask of the rows excluded by ``excluded_wells``.

        Accepts a flat set of well IDs OR the per-(gene, sample) dict the app
        actually stores. The dict form is what makes this correct: ``Well`` is the
        raw plate coordinate with no per-file disambiguation, and the app's data
        frame is a concat of every uploaded file, so "A1" recurs once per plate
        (and once per target on a multiplex plate). Every QC *display* path used
        to flatten the dict to a bare set of well IDs, which then excluded that
        coordinate for EVERY gene and sample sharing it — so the QC screen showed
        a housekeeping triplicate at n=2 while ΔΔCt used all three wells. Two
        plates uploaded together is the normal workflow here, so this was live.
        """
        if df is None or len(df) == 0:
            return pd.Series(dtype=bool, index=getattr(df, "index", None))
        if not excluded_wells:
            return pd.Series(False, index=df.index)
        if isinstance(excluded_wells, dict):
            if not {"Target", "Sample"} <= set(df.columns):
                flat = set()
                for wells in excluded_wells.values():
                    flat |= set(wells or ())
                return df["Well"].isin(flat)
            return pd.Series(
                [
                    well in excluded_wells.get((target, sample), ())
                    for target, sample, well in zip(
                        df["Target"], df["Sample"], df["Well"]
                    )
                ],
                index=df.index,
            )
        return df["Well"].isin(excluded_wells)

    @staticmethod
    def get_triplicate_data(
        data: pd.DataFrame, excluded_wells: set = None, thresholds: dict = None
    ) -> pd.DataFrame:
        """
        Build a comprehensive triplicate-level view of all CT values.
        Returns DataFrame with one row per well, grouped by Sample+Target.
        """
        _qct = QualityControl.resolve_thresholds(thresholds)
        if data is None or data.empty:
            return pd.DataFrame()

        excluded_wells = excluded_wells or set()

        df = data[["Well", "Sample", "Target", "CT"]].copy()
        df["Excluded"] = QualityControl.excluded_mask(df, excluded_wells)

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

        # Restore the (Sample, Target) groups that lost every well to exclusions.
        # Grouping only the survivors dropped them entirely instead of leaving an
        # n=0 row that trips the "Low n" check below — so excluding a whole bad
        # triplicate RAISED the health score (2 groups / 50% healthy became
        # 1 group / 100%).
        all_groups = df[["Sample", "Target"]].drop_duplicates()
        if len(all_groups) > len(group_stats):
            group_stats = all_groups.merge(
                group_stats, on=["Sample", "Target"], how="left"
            ).reset_index(drop=True)
            group_stats["n"] = group_stats["n"].fillna(0).astype(int)
            group_stats["CT_Values"] = group_stats["CT_Values"].apply(
                lambda v: v if isinstance(v, list) else []
            )
            group_stats["Wells"] = group_stats["Wells"].apply(
                lambda v: v if isinstance(v, list) else []
            )

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
            if pd.notna(row["CV_pct"]) and row["CV_pct"] > _qct["cv"] * 100:
                issues.append(f"High CV ({row['CV_pct']:.1f}%)")
                severity = "warning"
            if row["Mean_CT"] > _qct["ct_high"]:
                issues.append("High CT")
                severity = "warning"
            if row["Mean_CT"] < _qct["ct_low"]:
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
        data: pd.DataFrame, sample: str, target: str, thresholds: dict = None
    ) -> pd.DataFrame:
        """
        Get all wells for a specific Sample+Target combination with detailed info.
        """
        _qct = QualityControl.resolve_thresholds(thresholds)
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
            if row["CT"] > _qct["ct_high"]:
                issues.append("High CT")
            if row["CT"] < _qct["ct_low"]:
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
    def get_qc_summary_stats(data: pd.DataFrame, excluded_wells: set = None,
                             thresholds: dict = None) -> dict:
        """
        Calculate overall QC summary statistics for the dataset.
        """
        _qct = QualityControl.resolve_thresholds(thresholds)
        if data is None or data.empty:
            return {}

        excluded_wells = excluded_wells or set()
        _excl_mask = QualityControl.excluded_mask(data, excluded_wells)
        active_data = data[~_excl_mask]

        total_wells = len(data)
        excluded_count = int(_excl_mask.sum())
        active_wells = total_wells - excluded_count

        ct_mean = active_data["CT"].mean()
        ct_std = active_data["CT"].std()
        ct_min = active_data["CT"].min()
        ct_max = active_data["CT"].max()

        high_ct_count = len(
            active_data[active_data["CT"] > _qct["ct_high"]]
        )
        low_ct_count = len(
            active_data[active_data["CT"] < _qct["ct_low"]]
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
    def find_high_sd_outliers(
        data: pd.DataFrame,
        excluded_wells,
        sd_threshold: float = 0.5,
        gene_filter: str = None,
    ) -> list:
        """Find the worst replicate in each gene-sample group exceeding SD threshold."""
        if data is None or data.empty:
            return []

        df = data.copy()

        # Scope the exclusion to its own (gene, sample) — see excluded_mask.
        # Flattening removed the coordinate from every group sharing it, which
        # hid genuine high-SD groups from this suggestion list.
        if excluded_wells:
            _mask = QualityControl.excluded_mask(df, excluded_wells)
            if _mask.any():
                df = df[~_mask]

        if gene_filter:
            df = df[df["Target"] == gene_filter]

        suggestions = []
        for (target, sample), group in df.groupby(["Target", "Sample"]):
            if len(group) < 3:
                continue

            ct_values = group["CT"].values
            sd = np.std(ct_values, ddof=1)

            if sd <= sd_threshold:
                continue

            mean_ct = np.mean(ct_values)
            deviations = np.abs(ct_values - mean_ct)
            worst_idx = np.argmax(deviations)
            worst_row = group.iloc[worst_idx]

            suggestions.append({
                "Target": target,
                "Sample": sample,
                "Well": worst_row["Well"],
                "CT": round(float(worst_row["CT"]), 2),
                "deviation": round(float(worst_row["CT"] - mean_ct), 3),
                "group_sd": round(float(sd), 3),
                "group_mean": round(float(mean_ct), 2),
                "n_replicates": len(group),
            })

        return suggestions

    @staticmethod
    def auto_select_replicates(
        data: pd.DataFrame,
        sd_threshold: float = 0.3,
        excluded_wells=None,
    ) -> Tuple[dict, list]:
        """Automatic best-2-of-3 triplicate QC by minimum pairwise SD.

        For each (Target, Sample) group:
          * SD (ddof=1) of the replicate CTs <= threshold  -> keep all ('ok',
            not surfaced).
          * otherwise drop the replicate farthest from the mean, iteratively,
            until SD <= threshold or only 2 remain. For a triplicate this leaves
            the two closest CT values — provably the minimum-pairwise-SD pair
            (dropping the farther extreme == keeping the closer adjacent pair).
          * if even the best pair still exceeds the threshold, keep those 2 and
            flag the group 'unresolved'. A whole condition is never dropped.

        Args:
            data: long-format frame with Well, Sample, Target, CT.
            sd_threshold: max acceptable replicate CT SD (default 0.3).
            excluded_wells: optional dict {(gene, sample): set} or flat set of
                wells already removed — ignored here, not reconsidered.

        Returns:
            (exclusions, audit)
              exclusions: {(target, sample): set(dropped_wells)} for groups where
                  wells were dropped.
              audit: list of per-group dicts for every trimmed/unresolved group,
                  each with Target, Sample, n, orig_sd, final_sd, kept_wells,
                  dropped_wells, status.
        """
        exclusions: dict = {}
        audit: list = []
        if data is None or getattr(data, "empty", True):
            return exclusions, audit

        df = data[["Well", "Sample", "Target", "CT"]].copy()
        # Scope the pre-exclusion to its own (gene, sample) via excluded_mask.
        # Flattening to bare well IDs dropped that plate coordinate from EVERY
        # group sharing it, so auto-QC evaluated an unrelated gene's triplicate
        # on the wrong wells: a group carrying a 10-Ct outlier was reported at
        # SD 0.00 and left untrimmed, while DDCt kept using the outlier.
        if excluded_wells:
            _pre_mask = QualityControl.excluded_mask(df, excluded_wells)
            if _pre_mask.any():
                df = df[~_pre_mask]
        # Only real numeric CTs participate in pair selection.
        df = df[pd.to_numeric(df["CT"], errors="coerce").notna()]

        def _sd(vals):
            return float(np.std(vals, ddof=1)) if len(vals) >= 2 else float("nan")

        for (target, sample), group in df.groupby(["Target", "Sample"]):
            wells = list(group["Well"])
            cts = [float(c) for c in group["CT"]]
            n = len(cts)
            if n < 2:
                continue

            orig_sd = _sd(cts)
            if not np.isfinite(orig_sd) or orig_sd <= sd_threshold:
                continue  # 'ok' — nothing to do

            keep_idx = list(range(n))
            dropped_idx = []
            while len(keep_idx) > 2:
                vals = np.array([cts[i] for i in keep_idx])
                if _sd(vals) <= sd_threshold:
                    break
                mean = vals.mean()
                far = int(np.argmax(np.abs(vals - mean)))
                dropped_idx.append(keep_idx.pop(far))

            final_sd = _sd([cts[i] for i in keep_idx])
            status = "trimmed" if final_sd <= sd_threshold else "unresolved"
            dropped_wells = {wells[i] for i in dropped_idx}
            if dropped_wells:
                exclusions[(target, sample)] = dropped_wells
            audit.append({
                "Target": target,
                "Sample": sample,
                "n": n,
                "orig_sd": round(orig_sd, 3),
                "final_sd": round(final_sd, 3) if np.isfinite(final_sd) else None,
                "kept_wells": [wells[i] for i in keep_idx],
                "dropped_wells": sorted(dropped_wells),
                "status": status,
            })

        return exclusions, audit

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
    def detect_outliers(data: pd.DataFrame, hk_gene: str = None,
                        thresholds: dict = None) -> pd.DataFrame:
        _qct = QualityControl.resolve_thresholds(thresholds)
        if data is None or data.empty:
            return pd.DataFrame()

        qc_df = data[["Well", "Sample", "Target", "CT"]].copy()

        ct_high = qc_df["CT"] > _qct["ct_high"]
        ct_low = qc_df["CT"] < _qct["ct_low"]

        high_ct_issue = f"CT > {_qct["ct_high"]} (low expression)"
        low_ct_issue = f"CT < {_qct["ct_low"]} (unusually high)"

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
        high_cv_groups = cv_stats[cv_stats["cv"] > _qct["cv"]][
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
                    deviations > _qct["hk_variation"]
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
    def get_replicate_stats(data: pd.DataFrame,
                            thresholds: dict = None) -> pd.DataFrame:
        _qct = QualityControl.resolve_thresholds(thresholds)
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

        # Use the configured thresholds, not literals. These were hardcoded to
        # 10 / 35 / 5, so the QC Settings sliders moved the summary counts but
        # never this table's Status column — raising CT_HIGH to 40 took the
        # "high CT" count to zero while these rows still read "Low Expression".
        # (_cached_replicate_stats also puts the threshold tuple in its cache key,
        # which advertised a dependency that did not exist until now.)
        rep_stats["Status"] = np.select(
            [
                rep_stats["Mean CT"] < _qct["ct_low"],
                rep_stats["Mean CT"] > _qct["ct_high"],
                rep_stats["CV%"] > _qct["cv"] * 100,
            ],
            ["Check Signal", "Low Expression", "High CV"],
            default="OK",
        )

        rep_stats["Mean CT"] = rep_stats["Mean CT"].round(2)
        rep_stats["SD"] = rep_stats["SD"].round(3)
        rep_stats["CV%"] = rep_stats["CV%"].round(1)
        rep_stats["n"] = rep_stats["n"].astype(int)

        return rep_stats[["Sample", "Target", "n", "Mean CT", "SD", "CV%", "Status"]]
