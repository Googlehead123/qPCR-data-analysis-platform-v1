"""Export functions for qPCR analysis results.

Provides multi-sheet Excel export with raw data, gene-by-gene analysis,
QC report, fold-change matrix, and summary.
"""

import io
from typing import Dict

import numpy as np
import pandas as pd


def export_to_excel(
    raw_data: pd.DataFrame,
    processed_data: Dict[str, pd.DataFrame],
    params: dict,
    mapping: dict,
    qc_stats: dict = None,
    replicate_stats: pd.DataFrame = None,
) -> bytes:
    """Export comprehensive Excel with gene-by-gene sheets, QC report, and FC matrix.

    Args:
        raw_data: Raw CT data.
        processed_data: Dict of gene -> DataFrame with analysis results.
        params: Analysis parameters dict.
        mapping: Sample mapping dict.
        qc_stats: Optional QC summary stats dict from QualityControl.get_qc_summary_stats().
        replicate_stats: Optional replicate stats DataFrame from QualityControl.get_replicate_stats().
    """
    output = io.BytesIO()

    with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
        # Parameters sheet
        pd.DataFrame([params]).to_excel(
            writer, sheet_name="Analysis_Parameters", index=False
        )

        # Sample mapping sheet
        pd.DataFrame([{"Original": k, **v} for k, v in mapping.items()]).to_excel(
            writer, sheet_name="Sample_Mapping", index=False
        )

        # Raw data (include mapped Condition column reflecting sample mapping)
        raw_export = raw_data.copy()
        if mapping:
            raw_export["Condition"] = raw_export["Sample"].map(
                lambda x: mapping.get(x, {}).get("condition", x)
            )
        else:
            raw_export["Condition"] = raw_export["Sample"]
        raw_export = (
            raw_export[["Well", "Sample", "Condition", "Target", "CT", "Source_File"]]
            if "Source_File" in raw_export.columns
            else raw_export
        )
        raw_export.to_excel(writer, sheet_name="Raw_Data", index=False)

        # Gene-by-gene calculations with statistical test method column
        for gene, gene_data in processed_data.items():
            sheet_name = f"{gene}_Analysis"[:31]
            gene_export = gene_data.copy()
            # Add statistical test method column
            if "p_value" in gene_export.columns:
                ttest_type = params.get("ttest_type", "welch")
                gene_export["Stat_Test"] = ""
                for idx, row in gene_export.iterrows():
                    n_rep = row.get("n_replicates", 0)
                    if pd.notna(row.get("p_value")) and row.get("p_value") is not np.nan:
                        if n_rep >= 2:
                            gene_export.loc[idx, "Stat_Test"] = (
                                f"{'Welch' if ttest_type == 'welch' else 'Student'} t-test (n={n_rep})"
                            )
                        elif n_rep == 1:
                            gene_export.loc[idx, "Stat_Test"] = f"One-sample t-test (n={n_rep})"
                        else:
                            gene_export.loc[idx, "Stat_Test"] = "N/A"
            gene_export.to_excel(writer, sheet_name=sheet_name, index=False)

        # Summary sheet
        if processed_data:
            non_empty = [df for df in processed_data.values() if not df.empty]
            if non_empty:
                all_data = pd.concat(non_empty, ignore_index=True)
                agg_dict = {"Relative_Expression": ["mean", "std", "count"]}
                if "p_value" in all_data.columns:
                    agg_dict["p_value"] = "min"
                group_cols = ["Target"]
                if "Group" in all_data.columns:
                    group_cols.append("Group")
                summary = (
                    all_data.groupby(group_cols)
                    .agg(agg_dict)
                    .round(4)
                )
                summary.to_excel(writer, sheet_name="Summary")

        # --- NEW: Fold Change Matrix (pivot table) ---
        if processed_data:
            non_empty = [df for df in processed_data.values() if not df.empty]
            if non_empty:
                all_data = pd.concat(non_empty, ignore_index=True)
                if "Fold_Change" in all_data.columns and "Condition" in all_data.columns:
                    fc_matrix = all_data.pivot_table(
                        values="Fold_Change",
                        index="Target",
                        columns="Condition",
                        aggfunc="first",
                    )
                    fc_matrix = fc_matrix.round(4)
                    fc_matrix.to_excel(writer, sheet_name="FC_Matrix")

        # --- NEW: QC Report sheet ---
        _write_qc_sheet(writer, qc_stats, replicate_stats)

    return output.getvalue()


def _write_qc_sheet(
    writer,
    qc_stats: dict = None,
    replicate_stats: pd.DataFrame = None,
):
    """Write QC Report sheet with summary stats and replicate-level CV/outlier info."""
    rows = []

    if qc_stats:
        rows.append({"Metric": "Total Wells", "Value": qc_stats.get("total_wells", "")})
        rows.append({"Metric": "Excluded Wells", "Value": qc_stats.get("excluded_wells", "")})
        rows.append({"Metric": "Active Wells", "Value": qc_stats.get("active_wells", "")})
        rows.append({"Metric": "CT Mean", "Value": qc_stats.get("ct_mean", "")})
        rows.append({"Metric": "CT SD", "Value": qc_stats.get("ct_std", "")})
        rows.append({"Metric": "CT Range", "Value": f"{qc_stats.get('ct_min', '')}-{qc_stats.get('ct_max', '')}"})
        rows.append({"Metric": "High CT Wells (>35)", "Value": qc_stats.get("high_ct_count", "")})
        rows.append({"Metric": "Low CT Wells (<10)", "Value": qc_stats.get("low_ct_count", "")})
        rows.append({"Metric": "Total Triplicates", "Value": qc_stats.get("total_triplicates", "")})
        rows.append({"Metric": "Healthy Triplicates", "Value": qc_stats.get("healthy_triplicates", "")})
        rows.append({"Metric": "Warning Triplicates", "Value": qc_stats.get("warning_triplicates", "")})
        rows.append({"Metric": "Error Triplicates", "Value": qc_stats.get("error_triplicates", "")})
        rows.append({"Metric": "Avg CV%", "Value": qc_stats.get("avg_cv_pct", "")})
        rows.append({"Metric": "Max CV%", "Value": qc_stats.get("max_cv_pct", "")})
        rows.append({"Metric": "Health Score (%)", "Value": qc_stats.get("health_score", "")})

    if rows:
        pd.DataFrame(rows).to_excel(writer, sheet_name="QC_Report", index=False, startrow=0)

    if replicate_stats is not None and not replicate_stats.empty:
        start_row = len(rows) + 3 if rows else 0
        replicate_stats.to_excel(
            writer, sheet_name="QC_Report", index=False, startrow=start_row
        )
