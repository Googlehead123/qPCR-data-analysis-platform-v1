"""p-values must reproduce the hand calculation: two-tailed Student's t-test.

Min reported app p-values HIGHER than the manual calculation for some
comparisons. Root cause: the statistics defaulted to Welch's t-test. With equal
group sizes Welch has the same t statistic as Student's but fewer degrees of
freedom, so its p is >= Student's — and strictly greater whenever the two
groups' variances differ (and t != 0), which is exactly "some p-values
higher". On the six
plates in test_data/, Welch exceeded Student in 90/90 comparisons (median 2x,
max 156x) and moved the significance marker in 21 of them.
Decision (Min, 2026-10-23): two-tailed Student's t-test is the default; Welch
stays selectable.

Every expected value here is derived by hand from the Ct table below — simple
means, STDEV with n-1, SEM = SD/sqrt(n), pooled variance, t-distribution tail
via scipy.special.stdtr — never by calling the function under test or
scipy.stats.ttest_ind. The pinned literals are the same numbers Excel's
T.TEST(control, sample, 2, 2) returns for these ΔCt values.
"""
import numpy as np
import pandas as pd
import pytest
from scipy import special

# Ct values: 3 wells per gene per condition. "None" is the negative control,
# "B" a specified comparison condition. A's target replicates are deliberately
# wider than the control's, which is the case where Welch and Student diverge.
CT = {
    #        HK (ACTIN)             target G
    "None": ([18.00, 18.10, 17.90], [25.00, 25.05, 24.95]),
    "A":    ([18.20, 18.00, 18.10], [23.00, 23.60, 22.70]),
    "B":    ([18.00, 18.00, 18.00], [24.00, 24.20, 24.10]),
}

# Per-well ΔCt = target Ct - mean HK Ct of the same condition.
DCT = {c: np.array(g) - np.mean(hk) for c, (hk, g) in CT.items()}


def _raw():
    rows = []
    for cond, (hk, g) in CT.items():
        for i, ct in enumerate(hk, 1):
            rows.append((f"{cond}-HK{i}", cond, "ACTIN", ct))
        for i, ct in enumerate(g, 1):
            rows.append((f"{cond}-G{i}", cond, "G", ct))
    return pd.DataFrame(rows, columns=["Well", "Sample", "Target", "CT"])


def _student_p(ref, cond):
    """Pooled-variance two-tailed t-test, written out (Excel T.TEST type 2)."""
    ref, cond = np.asarray(ref, float), np.asarray(cond, float)
    n1, n2 = ref.size, cond.size
    sp2 = ((n1 - 1) * ref.var(ddof=1) + (n2 - 1) * cond.var(ddof=1)) / (n1 + n2 - 2)
    t = (ref.mean() - cond.mean()) / np.sqrt(sp2 * (1 / n1 + 1 / n2))
    return 2 * special.stdtr(n1 + n2 - 2, -abs(t))


def _welch_p(ref, cond):
    """Welch-Satterthwaite two-tailed t-test, written out (Excel T.TEST type 3)."""
    ref, cond = np.asarray(ref, float), np.asarray(cond, float)
    v1, v2 = ref.var(ddof=1) / ref.size, cond.var(ddof=1) / cond.size
    t = (ref.mean() - cond.mean()) / np.sqrt(v1 + v2)
    df = (v1 + v2) ** 2 / (v1 ** 2 / (ref.size - 1) + v2 ** 2 / (cond.size - 1))
    return 2 * special.stdtr(df, -abs(t))


def _run(compare, excluded_wells=None, **kw):
    from qpcr.analysis import AnalysisEngine
    raw = _raw()
    mapping = {c: {"condition": c} for c in CT}
    excluded_wells = excluded_wells or {}
    processed = AnalysisEngine.calculate_ddct(
        raw, "ACTIN", "None", excluded_wells, set(), mapping
    )
    stats_df = AnalysisEngine.calculate_statistics(
        processed, compare, raw_data=raw, hk_gene="ACTIN",
        sample_mapping=mapping, excluded_wells=excluded_wells, **kw,
    )
    return stats_df.set_index("Condition")


class TestSummaryStatisticsMatchHandCalculation:
    def test_mean_sd_sem_and_fold_change(self, mock_streamlit):
        row = _run("None").loc["A"]
        assert row["n_replicates"] == 3
        assert row["Target_Ct_Mean"] == pytest.approx(23.1)            # (23.0+23.6+22.7)/3
        assert row["HK_Ct_Mean"] == pytest.approx(18.1)
        assert row["Target_Ct_SD"] == pytest.approx(np.sqrt(0.21))     # STDEV, n-1
        assert row["SD"] == pytest.approx(np.sqrt(0.21))
        assert row["SEM"] == pytest.approx(np.sqrt(0.21) / np.sqrt(3))
        assert row["Delta_Ct"] == pytest.approx(5.0)                   # 23.1 - 18.1
        assert row["Delta_Delta_Ct"] == pytest.approx(-2.0)            # 5.0 - 7.0
        assert row["Fold_Change"] == pytest.approx(4.0)                # 2^2


class TestDefaultIsTwoTailedStudent:
    def test_default_vs_negative_control_matches_excel_ttest_type2(self, mock_streamlit):
        row = _run("None").loc["A"]
        assert row["p_value"] == pytest.approx(_student_p(DCT["None"], DCT["A"]), rel=1e-9)
        assert row["p_value"] == pytest.approx(0.0016784148735, rel=1e-9)
        assert row["significance"] == "**"
        assert row["stat_test_used"] == "Student t-test (n=3 vs ref n=3)"

    def test_default_is_not_welch(self, mock_streamlit):
        """The regression itself: Welch gives 0.0162 ('*') on these wells."""
        row = _run("None").loc["A"]
        assert row["p_value"] < _welch_p(DCT["None"], DCT["A"]) / 5

    @pytest.mark.parametrize("ref", ["None", "B"])
    def test_any_reference_condition_uses_student(self, mock_streamlit, ref):
        """'Condition vs None, or vs a specified condition' — same test either way."""
        row = _run(ref).loc["A"]
        assert row["p_value"] == pytest.approx(_student_p(DCT[ref], DCT["A"]), rel=1e-9)

    def test_second_comparison_uses_student(self, mock_streamlit):
        row = _run("None", compare_condition_2="B").loc["A"]
        assert row["p_value_2"] == pytest.approx(_student_p(DCT["B"], DCT["A"]), rel=1e-9)
        assert row["p_value_2"] == pytest.approx(0.0153232158254, rel=1e-9)
        assert row["stat_test_used_2"].startswith("Student t-test")


class TestQcReducedReplicates:
    def test_two_wells_after_qc_use_df3(self, mock_streamlit):
        """QC drops A's 23.60 well: 2 wells vs 3, pooled df = 2 + 3 - 2 = 3."""
        excluded = {("G", "A"): {"A-G2"}}
        row = _run("None", excluded_wells=excluded).loc["A"]
        assert row["n_replicates"] == 2
        assert row["Target_Ct_Mean"] == pytest.approx(22.85)
        assert row["Target_Ct_SD"] == pytest.approx(np.sqrt(0.045))    # (23.0, 22.7)
        assert row["SEM"] == pytest.approx(np.sqrt(0.045) / np.sqrt(2))
        assert row["p_value"] == pytest.approx(
            _student_p(DCT["None"], DCT["A"][[0, 2]]), rel=1e-9
        )
        assert row["p_value"] == pytest.approx(0.000313799085425, rel=1e-9)
        assert row["stat_test_used"] == "Student t-test (n=2 vs ref n=3)"

    def test_excluded_hk_well_moves_the_test_baseline(self, mock_streamlit):
        """An HK exclusion must change the ΔCt the test sees, as it does the bar."""
        excluded = {("ACTIN", "A"): {"A-HK1"}}       # drop HK 18.20 -> mean 18.05
        row = _run("None", excluded_wells=excluded).loc["A"]
        dct_a = np.array(CT["A"][1]) - 18.05
        assert row["HK_Ct_Mean"] == pytest.approx(18.05)
        assert row["p_value"] == pytest.approx(_student_p(DCT["None"], dct_a), rel=1e-9)


class TestTestSelection:
    def test_welch_remains_available(self, mock_streamlit):
        row = _run("None", ttest_type="welch").loc["A"]
        assert row["p_value"] == pytest.approx(_welch_p(DCT["None"], DCT["A"]), rel=1e-9)
        assert row["stat_test_used"].startswith("Welch t-test")

    def test_unknown_test_name_is_rejected(self, mock_streamlit):
        """Anything but 'student' used to fall through to Welch silently."""
        with pytest.raises(ValueError, match="ttest_type"):
            _run("None", ttest_type="studnet")

    def test_default_constant_is_student(self, mock_streamlit):
        from qpcr.analysis import DEFAULT_TTEST_TYPE, TTEST_TYPES
        assert DEFAULT_TTEST_TYPE == "student"
        assert set(TTEST_TYPES) == {"student", "welch"}
