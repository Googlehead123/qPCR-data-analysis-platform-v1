"""QPCRParser — CSV parsing for QuantStudio format qPCR data.

Supports Format 1 (Well Position) and Format 2 (Well/Sample Name/Cт) layouts
with automatic detection, encoding fallback, and file size validation.
"""

import streamlit as st
import pandas as pd


class QPCRParser:
    MAX_FILE_SIZE_MB = 50

    # Column names that identify each field, lower-cased and stripped.
    _WELL_NAMES = ("well", "well position")
    _SAMPLE_NAMES = ("sample", "sample name")
    _TARGET_NAMES = ("target", "target name")
    _CT_NAMES = ("ct", "cт", "cq")

    @staticmethod
    def header_in_columns(names) -> bool:
        """True if ``names`` is itself the data header rather than a preamble row.

        ``detect_format`` only ever scanned ROWS, but pandas consumes line 1 of
        the file as ``df.columns`` — so a CSV whose real header IS line 1 (no
        instrument preamble, which is what a hand-made or re-exported file looks
        like, including this repo's own test_qpcr.csv) could never be detected
        and was rejected outright.
        """
        low = [str(n).strip().lower() for n in names]
        return (
            any(n in QPCRParser._WELL_NAMES for n in low)
            and any(n in QPCRParser._SAMPLE_NAMES for n in low)
            and any(n in QPCRParser._TARGET_NAMES for n in low)
            and any(n in QPCRParser._CT_NAMES for n in low)
        )

    @staticmethod
    def detect_format(df):
        if QPCRParser.header_in_columns(df.columns):
            # Signal to parse() that the header must be pushed back into the
            # frame as row 0 before the format parsers see it.
            return "format1_header_row", 0
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
            (c for c in df.columns
             if str(c).strip().upper() in ["CT", "CТ", "CQ"] or str(c).strip() == "Cт"),
            None,
        )

        if not ct_col:
            st.error("CT column not found. Expected a column named 'CT', 'Ct', 'Cq' or 'Cт'.")
            return None

        # Look the columns up by name. This used to fall back to df.iloc[:, 2]
        # and df.iloc[:, 3] with no validation and no message, so a header of
        # "Well Position,Sample,Target,CT" — which selects this format — silently
        # produced Sample=<the target name> and Target=<the CT value>, reported
        # "2 wells parsed", and listed CT values as gene names.
        sample_col_name = next(
            (c for c in df.columns if str(c).strip() in ("Sample Name", "Sample")), None
        )
        target_col_name = next(
            (c for c in df.columns if str(c).strip() in ("Target Name", "Target")), None
        )
        if sample_col_name is None or target_col_name is None:
            missing = [
                label
                for label, col in (("Sample Name", sample_col_name),
                                   ("Target Name", target_col_name))
                if col is None
            ]
            st.error(
                f"Column(s) not found: {', '.join(missing)}. Found columns: "
                f"{', '.join(str(c) for c in df.columns[:10])}. Nothing was loaded "
                f"from this file."
            )
            return None
        sample_col = df[sample_col_name]
        target_col = df[target_col_name]

        parsed = pd.DataFrame(
            {
                "Well": df[well_col].astype(str).str.strip(),
                "Sample": sample_col.astype(str).str.strip(),
                "Target": target_col.astype(str).str.strip(),
                "CT": pd.to_numeric(df[ct_col], errors="coerce"),
            }
        )
        parsed["Sample"] = parsed["Sample"].replace({"nan": pd.NA, "None": pd.NA, "": pd.NA})
        parsed["Target"] = parsed["Target"].replace({"nan": pd.NA, "None": pd.NA, "": pd.NA})
        parsed["Well"] = parsed["Well"].replace({"nan": pd.NA, "None": pd.NA, "": pd.NA})

        # Count invalid CT values for user feedback
        invalid_ct_count = parsed["CT"].isna().sum()
        if invalid_ct_count > 0:
            st.info(
                f"Note: {invalid_ct_count} rows with invalid/undetermined CT values were filtered out."
            )

        pre_filter_count = len(parsed.dropna(subset=["CT"]))
        result = parsed.dropna(subset=["CT", "Sample", "Target"])
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
                (c for c in df.columns if str(c).strip() in ["Cт", "CT", "Ct", "Cq", "CQ"]), None
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
                    missing.append("CT/Cq/Cт")
                st.error(
                    f"Format2 parsing failed: Missing columns: {', '.join(missing)}"
                )
                return None

            parsed = pd.DataFrame(
                {
                    "Well": df[well_col].astype(str).str.strip(),
                    "Sample": df[sample_col].astype(str).str.strip(),
                    "Target": df[target_col].astype(str).str.strip(),
                    "CT": pd.to_numeric(df[ct_col], errors="coerce"),
                }
            )
            parsed["Sample"] = parsed["Sample"].replace({"nan": pd.NA, "None": pd.NA, "": pd.NA})
            parsed["Target"] = parsed["Target"].replace({"nan": pd.NA, "None": pd.NA, "": pd.NA})
            parsed["Well"] = parsed["Well"].replace({"nan": pd.NA, "None": pd.NA, "": pd.NA})

            # Report dropped CT rows the way parse_format1 does. "Undetermined"
            # is in na_values, so a genuine no-amplification well is coerced to
            # NaN and the ROW IS DELETED — leaving a triplicate silently at n=2
            # with SD/SEM/p-values computed as if all three amplified, and no
            # trace in the QC tab or the provenance record.
            invalid_ct_count = int(parsed["CT"].isna().sum())
            if invalid_ct_count > 0:
                st.info(
                    f"Note: {invalid_ct_count} rows with invalid/undetermined CT "
                    f"values were filtered out."
                )

            pre_filter_count = len(parsed.dropna(subset=["CT"]))
            result = parsed.dropna(subset=["CT", "Sample", "Target"])
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

            # Encoding fallback chain. Decoding without an exception is NOT
            # proof of the right codec: utf-16-le happily decodes CP949 bytes
            # (what Excel on a Korean Windows writes for plain "CSV (Comma
            # delimited)") into mojibake column names, and latin-1 decodes any
            # byte string at all. So an encoding is accepted only once the layout
            # is recognisable, and the first merely-decodable result is kept as a
            # fallback so the error message can still show what was read.
            # utf-8-sig leads so a BOM is stripped instead of being glued to the
            # first column name.
            df = None
            first_readable = None
            for enc in ["utf-8-sig", "utf-8", "utf-16", "utf-16-le",
                        "cp949", "euc-kr", "latin-1"]:
                try:
                    candidate = pd.read_csv(
                        file, encoding=enc, low_memory=False, skip_blank_lines=False,
                        keep_default_na=False,
                        na_values=["", "NA", "N/A", "NaN", "#N/A", "#N/A N/A",
                                   "#NA", "-NaN", "-nan", "nan", "<NA>",
                                   "Undetermined", "undetermined"],
                    )
                except (UnicodeDecodeError, UnicodeError):
                    file.seek(0)
                    continue
                except Exception:
                    file.seek(0)
                    continue
                # read_csv leaves the handle at EOF; rewind for the next attempt.
                file.seek(0)
                if first_readable is None:
                    first_readable = candidate
                if QPCRParser.detect_format(candidate)[0] != "unknown":
                    df = candidate
                    break

            if df is None:
                df = first_readable
            if df is None:
                st.error(
                    "Could not decode this file with any supported encoding "
                    "(UTF-8, UTF-16, CP949/EUC-KR). Re-save it as CSV UTF-8."
                )
                return None

            fmt, start = QPCRParser.detect_format(df)
            if fmt == "unknown":
                # This returned None with no message at all, and the caller only
                # reports on success or on empty-but-parsed — so an unreadable
                # file was completely silent AND the previous file's data stayed
                # live behind the new filename.
                st.error(
                    "Could not recognise this file's layout. Expected a "
                    "QuantStudio/StepOne export with a 'Well Position' header, or "
                    "a block whose header row contains Well / Sample Name / "
                    "Target Name / Cт. Two things to check: the header row is "
                    "intact, and the file is saved as UTF-8 rather than "
                    "ANSI/CP949 (in Excel: 'CSV UTF-8', not plain 'CSV'). "
                    f"Columns read: {', '.join(str(c) for c in list(df.columns)[:8])}"
                )
                return None
            if fmt == "format1_header_row":
                # Put the header back as row 0 so parse_format1 — which promotes
                # row 0 to the column names — works unchanged.
                df = pd.concat(
                    [pd.DataFrame([list(df.columns)], columns=df.columns), df],
                    ignore_index=True,
                )
                return QPCRParser.parse_format1(df, 0)
            return (
                QPCRParser.parse_format1(df, start)
                if fmt == "format1"
                else QPCRParser.parse_format2(df, start)
            )
        except Exception as e:
            st.error(f"Parse error: {e}")
            return None
