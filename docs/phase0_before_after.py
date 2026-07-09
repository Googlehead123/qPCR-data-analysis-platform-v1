import numpy as np
from scipy import stats
np.set_printoptions(precision=4)

print("="*74)
print("PHASE 0 CORRECTNESS FIXES — BEFORE/AFTER (magnitude sanity-check)")
print("="*74)

# ---------------------------------------------------------------------------
# FIX 1: Pooled-Condition well exclusion (calculate_ddct)
# ---------------------------------------------------------------------------
print("\n[1] Pooled-Condition well exclusion")
print("-"*74)
print("Scenario: Condition 'Treatment' pools raw samples D1 + D2 (donor pool).")
print("COL1A1 wells: D1=[24.0,24.0,24.0]  D2=[24.0,24.0, 30.0(bad, excluded)]")
print("HK (GAPDH) mean = 20.0; reference condition COL1A1 ΔCt = 5.0 (fold=1.0)\n")
hk = 20.0
ref_dct = 5.0  # reference ΔCt
# BEFORE: exclusion keyed to first sample only -> B6 (D2) NOT dropped
before_ct = np.array([24.0,24.0,24.0, 24.0,24.0,30.0])
# AFTER: row-wise exclusion -> B6 dropped
after_ct  = np.array([24.0,24.0,24.0, 24.0,24.0])
for label, ct in [("BEFORE (bug)", before_ct), ("AFTER  (fix)", after_ct)]:
    tmean = ct.mean(); dct = tmean - hk; ddct = dct - ref_dct; fc = 2**(-ddct)
    print(f"  {label}: target Ct mean={tmean:.3f}  ΔCt={dct:.3f}  ΔΔCt={ddct:.3f}  fold-change={fc:.3f}")
fc_b = 2**(-((before_ct.mean()-hk)-ref_dct)); fc_a = 2**(-((after_ct.mean()-hk)-ref_dct))
print(f"  => fold-change shifts {fc_b:.3f} -> {fc_a:.3f}  ({100*(fc_a-fc_b)/fc_b:+.1f}%). "
      f"The excluded outlier is now actually removed.")

# ---------------------------------------------------------------------------
# FIX 2: t-test on ΔCt (≈normal) instead of 2^-ΔCt (right-skewed)
# ---------------------------------------------------------------------------
print("\n[2] Significance test domain: ΔCt vs 2^-ΔCt")
print("-"*74)
print("Per-replicate ΔCt (target Ct - HK mean) for reference vs treatment groups.")
print("p-value BEFORE = ttest on 2^-ΔCt (fold-change);  AFTER = ttest on ΔCt.\n")

rng = np.random.RandomState(0)  # fixed seed, reproducible
scenarios = [
    ("Modest ~2x, tight",      np.array([5.0, 5.1, 4.9]),  np.array([4.0, 4.1, 3.9])),
    ("Strong ~8x, some spread",np.array([5.0, 5.3, 4.7]),  np.array([2.0, 2.6, 1.4])),
    ("Down ~0.3x",             np.array([5.0, 5.1, 4.9]),  np.array([6.7, 6.8, 6.6])),
    ("Noisy, borderline",      np.array([5.0, 5.8, 4.2]),  np.array([3.8, 4.9, 2.7])),
]
print(f"  {'scenario':<26}{'p BEFORE (2^-ΔCt)':>18}{'p AFTER (ΔCt)':>16}{'Δp':>10}")
for name, ref_dct_arr, cond_dct_arr in scenarios:
    ref_fc = 2**(-ref_dct_arr); cond_fc = 2**(-cond_dct_arr)
    _, p_before = stats.ttest_ind(ref_fc, cond_fc, equal_var=False)
    _, p_after  = stats.ttest_ind(ref_dct_arr, cond_dct_arr, equal_var=False)
    print(f"  {name:<26}{p_before:>18.4f}{p_after:>16.4f}{p_after-p_before:>+10.4f}")
print("\n  Note: two-sided p only; both use Welch. ΔCt better meets the t-test's")
print("  normality/variance assumptions, so AFTER is the statistically valid p-value.")
