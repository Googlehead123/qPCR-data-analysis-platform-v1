"""QPCRParser — CSV parsing for QuantStudio format qPCR data.

Supports Format 1 (Well Position) and Format 2 (Well/Sample Name/Cт) layouts
with automatic detection, encoding fallback, and file size validation.
"""

import streamlit as st
import pandas as pd


class QPCRParser:
    MAX_FILE_SIZE_MB = 50

    @staticmethod
    def detect_format(df):
        for idx, row in df.iterrows():
            if len(row) == 0 or row.isna().all():
                continue

            row_str = " ".join(row.astype(str).values)
            if "Well Position" in row_str:
                return "format1", idx
            elif len(row) > 0 and row.iloc[0] == "Well" and "Sample Name" in row_str:
                return (
                    "format2" if "Cт" in row_str or "ΔCт" in row_str else "format1"
                ), idx
        return "unknown", 0

    @staticmethod
    def parse_format1(df, start):
        df = df.iloc[start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)

        # FIX-10: Strip whitespace from column names (matching Format 2 behavior)
        df.columns = [str(c).strip() if pd.notna(c) else c for c in df.columns]

        if len(df.columns) < 4:
            st.error("CSV must have at least 4 columns (Well, Sample, Target, CT)")
            return None

        well_col = next(
            (c for c in ["Well Position", "Well"] if c in df.columns), df.columns[0]
        )

        # Case-insensitive CT column detection
        ct_col = next(
            (c for c in df.columns if str(c).upper() in ["CT", "CТ"] or str(c) == "Cт"),
            None,
        )

        if not ct_col:
            st.error("CT column not found. Expected column named 'CT', 'Ct', or 'Cт'.")
            return None

        sample_col = (
            df.get("Sample Name", df.iloc[:, 2])
            if "Sample Name" in df.columns
            else df.iloc[:, 2]
        )
        target_col = (
            df.get("Target Name", df.iloc[:, 3])
            if "Target Name" in df.columns
            else df.iloc[:, 3]
        )

        parsed = pd.DataFrame(
            {
                "Well": df[well_col],
                "Sample": sample_col,
                "Target": target_col,
                "CT": pd.to_numeric(df[ct_col], errors="coerce"),
            }
        )

        # Count invalid CT values for user feedback
        invalid_ct_count = parsed["CT"].isna().sum()
        if invalid_ct_count > 0:
            st.info(
                f"Note: {invalid_ct_count} rows with invalid/undetermined CT values were filtered out."
            )

        pre_filter_count = len(parsed.dropna(subset=["CT"]))
        result = parsed.dropna(subset=["CT"]).query("Sample.notna() & Target.notna()")
        dropped_sample_target = pre_filter_count - len(result)
        if dropped_sample_target > 0:
            st.info(
                f"Note: {dropped_sample_target} rows with missing Sample or Target names were filtered out."
            )

        if result.empty:
            st.warning(
                "No valid data rows found after filtering. Check that your file contains valid CT values and Sample/Target names."
            )

        return result

    @staticmethod
    def parse_format2(df, start):
        try:
            df = df.iloc[start:].reset_index(drop=True)
            df.columns = df.iloc[0]
            df = df.iloc[1:].reset_index(drop=True)

            # Try to find columns with flexible matching
            well_col = next((c for c in df.columns if str(c).strip() == "Well"), None)
            sample_col = next((c for c in df.columns if "Sample" in str(c)), None)
            target_col = next((c for c in df.columns if "Target" in str(c)), None)
            ct_col = next(
                (c for c in df.columns if str(c).strip() in ["Cт", "CT", "Ct"]), None
            )

            if not all([well_col, sample_col, target_col, ct_col]):
                missing = []
                if not well_col:
                    missing.append("Well")
                if not sample_col:
                    missing.append("Sample Name")
                if not target_col:
                    missing.append("Target Name")
                if not ct_col:
                    missing.append("CT/Cт")
                st.error(
                    f"Format2 parsing failed: Missing columns: {', '.join(missing)}"
                )
                return None

            parsed = pd.DataFrame(
                {
                    "Well": df[well_col],
                    "Sample": df[sample_col],
                    "Target": df[target_col],
                    "CT": pd.to_numeric(df[ct_col], errors="coerce"),
                }
            )

            pre_filter_count = len(parsed.dropna(subset=["CT"]))
            result = parsed.dropna(subset=["CT"]).query(
                "Sample.notna() & Target.notna()"
            )
            dropped_sample_target = pre_filter_count - len(result)
            if dropped_sample_target > 0:
                st.info(
                    f"Note: {dropped_sample_target} rows with missing Sample or Target names were filtered out."
                )
            if result.empty:
                st.warning(
                    "No valid data rows found after filtering. Check CT values and Sample/Target names."
                )
            return result

        except KeyError as e:
            st.error(f"Format2 parsing failed: Column not found - {e}")
            return None
        except Exception as e:
            st.error(f"Format2 parsing error: {e}")
            return None

    @staticmethod
    def parse(file):
        try:
            file.seek(0, 2)
            file_size_mb = file.tell() / (1024 * 1024)
            file.seek(0)

            if file_size_mb > QPCRParser.MAX_FILE_SIZE_MB:
                st.error(
                    f"File too large ({file_size_mb:.1f} MB). Maximum size is {QPCRParser.MAX_FILE_SIZE_MB} MB."
                )
                return None

            df = None
            # FIX-20: Extended encoding fallback chain with utf-16 support
            for enc in ["utf-8", "utf-16", "utf-16-le", "latin-1", "cp1252"]:
                try:
                    df = pd.read_csv(
                        file, encoding=enc, low_memory=False, skip_blank_lines=False
                    )
                    break
                except (UnicodeDecodeError, UnicodeError):
                    file.seek(0)
                    continue

            if df is None:
                return None

            fmt, start = QPCRParser.detect_format(df)
            return (
                QPCRParser.parse_format1(df, start)
                if fmt == "format1"
                else QPCRParser.parse_format2(df, start)
                if fmt == "format2"
                else None
            )
        except Exception as e:
            st.error(f"Parse error: {e}")
            return None
