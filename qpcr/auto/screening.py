"""Deterministic pre-analysis data-health screening.

`screen_data` inspects raw parsed qPCR data and returns a structured report of
issues a scientist should resolve before trusting the analysis — missing
columns, an undetectable reference gene, thin replicate groups, non-detects, and
duplicate wells. Pure function: no Streamlit, no global state.
"""

from __future__ import annotations

import pandas as pd

REQUIRED_COLUMNS = ("Well", "Sample", "Target", "CT")


def _issue(level: str, message: str) -> dict:
    return {"level": level, "message": message}


def screen_data(data: pd.DataFrame, hk_gene: str | None,
                min_replicates: int = 3) -> dict:
    """Screen raw qPCR data for common problems.

    Args:
        data: parsed long-format frame with Well/Sample/Target/CT columns.
        hk_gene: the chosen reference (housekeeping) gene, or None.
        min_replicates: recommended replicate count per (target, sample).

    Returns:
        {"ok": bool, "issues": [{level, message}], "summary": {...}}
        `ok` is False only when a blocking (error-level) issue is present.
    """
    issues: list[dict] = []

    if data is None or len(data) == 0:
        return {"ok": False, "issues": [_issue("error", "No data loaded.")],
                "summary": {}}

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in data.columns]
    if missing_cols:
        issues.append(_issue("error", f"Missing required column(s): {', '.join(missing_cols)}."))
        return {"ok": False, "issues": issues, "summary": {"n_rows": int(len(data))}}

    targets = sorted(str(t) for t in data["Target"].dropna().unique())
    samples = sorted(str(s) for s in data["Sample"].dropna().unique())

    # Reference gene presence (case-insensitive)
    hk_present = bool(hk_gene) and any(str(t).upper() == str(hk_gene).upper() for t in targets)
    if not hk_gene:
        issues.append(_issue("error", "No reference (housekeeping) gene selected."))
    elif not hk_present:
        issues.append(_issue("error",
            f"Reference gene '{hk_gene}' not found among targets ({', '.join(targets)})."))

    # Non-detects (NaN CT)
    n_nondetect = int(data["CT"].isna().sum())
    if n_nondetect:
        by_target = (
            data[data["CT"].isna()].groupby("Target").size().sort_values(ascending=False)
        )
        worst = "; ".join(f"{t}: {int(n)}" for t, n in by_target.head(3).items())
        issues.append(_issue("warning",
            f"{n_nondetect} non-detect / undetermined well(s) (excluded from means). {worst}"))

    # Replicate depth per (target, sample), on detected wells only
    detected = data[data["CT"].notna()]
    thin = []
    for (target, sample), grp in detected.groupby(["Target", "Sample"]):
        n = len(grp)
        if n < 2:
            thin.append((str(target), str(sample), n, "error"))
        elif n < min_replicates:
            thin.append((str(target), str(sample), n, "warning"))
    if thin:
        n_single = sum(1 for *_x, lvl in thin if lvl == "error")
        n_low = sum(1 for *_x, lvl in thin if lvl == "warning")
        if n_single:
            issues.append(_issue("warning",
                f"{n_single} target/sample group(s) have a single replicate "
                "(no SD / weak statistics)."))
        if n_low:
            issues.append(_issue("info",
                f"{n_low} group(s) have fewer than {min_replicates} replicates."))

    # Duplicate well IDs within the same (sample, target)
    dup = detected.duplicated(subset=["Sample", "Target", "Well"]).sum()
    if dup:
        issues.append(_issue("info", f"{int(dup)} duplicate well entries detected."))

    if not issues:
        issues.append(_issue("info", "No problems detected — data looks clean."))

    ok = not any(i["level"] == "error" for i in issues)
    summary = {
        "n_wells": int(data["Well"].nunique()) if "Well" in data else 0,
        "n_rows": int(len(data)),
        "n_targets": len(targets),
        "n_samples": len(samples),
        "hk_gene": hk_gene,
        "hk_present": hk_present,
        "n_nondetect": n_nondetect,
        "thin_groups": len(thin),
    }
    return {"ok": ok, "issues": issues, "summary": summary}
