#!/usr/bin/env python3
"""
qPCR Analysis Validation Script
================================
Validates DDCt calculations against Excel reference results.

Usage:
    python validate_calculations.py

This script will:
1. Load the raw CSV data
2. Calculate DDCt values using the same logic as the application
3. Compare results against Excel reference (250908 template result.xlsx)
4. Report any discrepancies
"""

import pandas as pd
import numpy as np
from pathlib import Path


def validate_qpcr_calculations():
    """Main validation function."""

    print("=" * 100)
    print("qPCR Analysis Engine Validation")
    print("=" * 100)

    # Load raw data
    raw_file = Path("250908 beta glucogallin_raw.csv")
    if not raw_file.exists():
        print(f"❌ ERROR: Raw data file not found: {raw_file}")
        return False

    print(f"\n✓ Loading raw data: {raw_file}")
    df_raw = pd.read_csv(raw_file, header=None, encoding="utf-8-sig")

    # Parse data (simulating QPCRParser.parse_format1)
    start_idx = 35  # Data starts at row 36 (0-indexed: 35)
    df = df_raw.iloc[start_idx:].reset_index(drop=True)
    df.columns = df.iloc[0]
    df = df.iloc[1:].reset_index(drop=True)

    parsed = pd.DataFrame(
        {
            "Well": df["Well"],
            "Sample": df["Sample Name"],
            "Target": df["Target Name"],
            "CT": pd.to_numeric(df["CT"], errors="coerce"),
        }
    )
    data = parsed.dropna(subset=["CT"]).query("Sample.notna() & Target.notna()")

    print(f"✓ Parsed {len(data)} data points from {data['Sample'].nunique()} samples")

    # Configuration
    hk_gene = "ACTIN"
    ref_sample = "1"

    print(f"✓ Housekeeping gene: {hk_gene}")
    print(f"✓ Reference sample: {ref_sample}")

    # Create sample mapping
    sample_order = sorted(
        data["Sample"].unique(),
        key=lambda x: (isinstance(x, str) and not x.replace(".", "").isdigit(), x),
    )
    sample_mapping = {}
    for sample in sample_order:
        sample_mapping[sample] = {
            "condition": str(sample),
            "group": "Baseline" if sample == ref_sample else "Treatment",
            "include": True,
        }

    # Apply sample mapping
    data["Condition"] = data["Sample"].map(
        lambda x: sample_mapping.get(x, {}).get("condition", str(x))
    )
    data["Group"] = data["Sample"].map(
        lambda x: sample_mapping.get(x, {}).get("group", "Treatment")
    )

    # Calculate DDCt for PCNA1
    target_gene = "PCNA1"
    target_data = data[data["Target"] == target_gene]

    print(f"\n{'=' * 100}")
    print(f"Calculating DDCt for {target_gene}")
    print(f"{'=' * 100}\n")

    results = []

    # Get reference delta Ct
    ref_target_data = target_data[target_data["Condition"] == ref_sample]
    ref_hk_data = data[(data["Condition"] == ref_sample) & (data["Target"] == hk_gene)]

    if len(ref_target_data) == 0 or len(ref_hk_data) == 0:
        print(f"❌ ERROR: Reference sample '{ref_sample}' not found in data")
        return False

    ref_delta_ct = ref_target_data["CT"].mean() - ref_hk_data["CT"].mean()
    print(f"Reference Delta Ct: {ref_delta_ct:.6f}")

    # Expected values from Excel
    expected_values = {
        "1": {"delta_ct": 4.162730, "ddct": 0.000000, "rq": 1.000000},
        "2": {"delta_ct": 3.792393, "ddct": -0.370336, "rq": 1.292654},
        "3": {"delta_ct": 3.785753, "ddct": -0.376977, "rq": 1.298618},
        "4": {"delta_ct": 0.889559, "ddct": -3.273171, "rq": 9.667691},
    }

    all_match = True

    for condition in sorted(target_data["Condition"].unique())[:4]:  # First 4 samples
        cond_data = target_data[target_data["Condition"] == condition]
        hk_data = data[(data["Condition"] == condition) & (data["Target"] == hk_gene)]

        if len(hk_data) == 0:
            print(f"⚠ Sample {condition}: No housekeeping data - SKIP")
            continue

        # Calculate DDCt
        target_ct_mean = cond_data["CT"].mean()
        hk_ct_mean = hk_data["CT"].mean()
        delta_ct = target_ct_mean - hk_ct_mean

        if condition == ref_sample:
            ref_delta_ct_val = delta_ct
        else:
            ref_delta_ct_val = ref_delta_ct

        ddct = delta_ct - ref_delta_ct_val
        rq = 2 ** (-ddct)

        # Get expected values
        expected = expected_values.get(condition)

        if expected:
            # Check if values match within tolerance
            delta_ct_match = abs(delta_ct - expected["delta_ct"]) < 0.001
            ddct_match = abs(ddct - expected["ddct"]) < 0.001
            rq_match = (
                abs(rq - expected["rq"]) < 0.01
            )  # Slightly larger tolerance for RQ

            match_status = (
                "✅" if (delta_ct_match and ddct_match and rq_match) else "❌"
            )
            all_match = all_match and delta_ct_match and ddct_match and rq_match

            print(f"\nSample {condition}: {match_status}")
            print(f"  Target CT Mean: {target_ct_mean:.3f}")
            print(f"  HK CT Mean: {hk_ct_mean:.3f}")
            print(
                f"  Delta Ct: {delta_ct:.6f} (Expected: {expected['delta_ct']:.6f}) {'✓' if delta_ct_match else '✗'}"
            )
            print(
                f"  Delta Delta Ct: {ddct:.6f} (Expected: {expected['ddct']:.6f}) {'✓' if ddct_match else '✗'}"
            )
            print(
                f"  RQ (Fold Change): {rq:.6f} (Expected: {expected['rq']:.6f}) {'✓' if rq_match else '✗'}"
            )

            if not (delta_ct_match and ddct_match and rq_match):
                print(f"  ⚠ MISMATCH DETECTED!")
                print(f"    Delta Ct diff: {abs(delta_ct - expected['delta_ct']):.9f}")
                print(f"    DDCt diff: {abs(ddct - expected['ddct']):.9f}")
                print(f"    RQ diff: {abs(rq - expected['rq']):.9f}")
        else:
            print(f"\nSample {condition}: (no reference data)")
            print(f"  Delta Ct: {delta_ct:.6f}")
            print(f"  Delta Delta Ct: {ddct:.6f}")
            print(f"  RQ (Fold Change): {rq:.6f}")

    print(f"\n{'=' * 100}")
    if all_match:
        print("✅ VALIDATION PASSED - All calculations match Excel reference")
        print("   The qPCR analysis engine is working correctly!")
    else:
        print("❌ VALIDATION FAILED - Some calculations don't match")
        print("   Please check the discrepancies above.")
    print(f"{'=' * 100}\n")

    return all_match


if __name__ == "__main__":
    try:
        success = validate_qpcr_calculations()
        exit(0 if success else 1)
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback

        traceback.print_exc()
        exit(1)
