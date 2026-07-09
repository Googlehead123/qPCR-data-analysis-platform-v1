"""Transparent statistical-test recommendation.

Given the per-replicate values of two or more groups, `recommend_test` runs the
standard decision tree — normality (Shapiro-Wilk) then equal-variance (Levene) —
and returns which test is appropriate *and the reasoning*, so the choice is
visible rather than hidden. It ADVISES; it does not run the test or change any
p-value the analysis pipeline already computed.

Decision tree
-------------
2 groups:
    normal + equal variance      -> Student's t-test
    normal + unequal variance    -> Welch's t-test
    non-normal                   -> Mann-Whitney U
>2 groups:
    normal                       -> one-way ANOVA (+ Tukey HSD post-hoc)
    non-normal                   -> Kruskal-Wallis (+ Dunn post-hoc)

Normality is only assessable at n>=3 per group; with smaller groups the test
defaults to the parametric branch but flags the low power.
"""

from __future__ import annotations

import numpy as np


def _clean(arr) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    return a[~np.isnan(a)]


def recommend_test(groups, alpha: float = 0.05) -> dict:
    """Recommend a significance test for a set of group value-arrays.

    Args:
        groups: list of 1D array-likes (per-replicate values, ideally ΔCt).
        alpha: significance threshold for the normality/variance screens.

    Returns:
        dict with keys: test, reason, n_groups, group_sizes, normal
        (True/False/None), shapiro_p (list), levene_p (float|None),
        equal_var (bool|None), followup (str|None).
    """
    from scipy import stats

    clean = [_clean(g) for g in (groups or [])]
    clean = [g for g in clean if g.size > 0]
    sizes = [int(g.size) for g in clean]
    n_groups = len(clean)

    base = {"n_groups": n_groups, "group_sizes": sizes, "normal": None,
            "shapiro_p": [], "levene_p": None, "equal_var": None, "followup": None}

    if n_groups < 2:
        return {**base, "test": None,
                "reason": "Need at least two non-empty groups to compare."}

    # --- Normality (Shapiro-Wilk), only where n>=3 and values vary ---
    testable = [g for g in clean if g.size >= 3 and np.ptp(g) > 0]
    shapiro_p = []
    for g in testable:
        try:
            shapiro_p.append(float(stats.shapiro(g).pvalue))
        except Exception:
            pass
    if len(testable) < n_groups or not shapiro_p:
        normal = None  # can't fully assess (small n or constant group)
    else:
        normal = all(p > alpha for p in shapiro_p)

    # --- Equal variance (Levene), needs each group to have >=2 values ---
    levene_p = None
    equal_var = None
    if all(g.size >= 2 for g in clean):
        try:
            levene_p = float(stats.levene(*clean).pvalue)
            equal_var = levene_p > alpha
        except Exception:
            levene_p, equal_var = None, None

    small_n_note = ""
    if any(s < 3 for s in sizes):
        small_n_note = " (n<3 in a group — normality not assessable; interpret with caution)"

    if n_groups == 2:
        if normal is False:
            test = "Mann-Whitney U test"
            reason = "Non-normal group(s) (Shapiro p<=%.2g) — use a rank-based test." % alpha
        else:
            if equal_var is False:
                test = "Welch's t-test"
                reason = ("Approximately normal but unequal variances "
                          f"(Levene p={levene_p:.3g})." if levene_p is not None
                          else "Approximately normal; unequal variances assumed (Welch is safe default).")
            else:
                test = "Student's t-test"
                reason = ("Approximately normal with equal variances"
                          + (f" (Levene p={levene_p:.3g})." if levene_p is not None else "."))
            reason += small_n_note
        followup = None
    else:  # >2 groups
        if normal is False:
            test = "Kruskal-Wallis test"
            reason = "Non-normal group(s) across >2 conditions — use a rank-based omnibus test."
            followup = "Dunn's test (post-hoc, pairwise)"
        else:
            test = "One-way ANOVA"
            reason = "Approximately normal across >2 conditions." + small_n_note
            followup = "Tukey HSD (post-hoc, pairwise)"

    return {**base, "test": test, "reason": reason, "normal": normal,
            "shapiro_p": shapiro_p, "levene_p": levene_p, "equal_var": equal_var,
            "followup": followup}
