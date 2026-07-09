# Phase 0 correctness fixes — before/after

Sanity-check of the two numbers-changing correctness fixes in this PR. Reproduce
with `python docs/phase0_before_after.py` (deterministic, fixed seed).

## Fix 1 — Pooled-Condition well exclusion (`calculate_ddct`)

**Bug:** exclusions were looked up via the *first* raw sample in a Condition
(`cond_data["Sample"].iloc[0]`), so when a Condition pools 2+ raw samples, an
excluded well on any other sample was silently kept in the ΔΔCt point estimate —
while the p-value and scatter (computed row-wise) *did* drop it. The bar could
disagree with its own statistics.

**Scenario:** Condition `Treatment` pools donors D1 + D2. COL1A1 wells
`D1=[24,24,24]`, `D2=[24, 24, 30(bad, excluded)]`; HK mean 20; reference ΔCt 5.

| | target Ct mean | ΔCt | ΔΔCt | fold-change |
|---|---|---|---|---|
| BEFORE (bug) | 25.000 | 5.000 | 0.000 | **1.000** |
| AFTER (fix) | 24.000 | 4.000 | −1.000 | **2.000** |

The excluded outlier is now actually removed → the bar matches its p-value/scatter.
Magnitude depends on how bad the excluded well is; here a single outlier flips the
result from "no change" to "2-fold up."

## Fix 2 — Significance test on ΔCt, not 2^-ΔCt

**Bug:** the t-test ran on the exponentiated fold-change `2^-ΔCt` (right-skewed),
violating the test's normality/variance assumptions. It now runs on `ΔCt`
(≈normal). Two-sided p only; both use Welch by default.

| scenario | p BEFORE (2^-ΔCt) | p AFTER (ΔCt) | Δp |
|---|---|---|---|
| Modest ~2×, tight | 0.0017 | 0.0003 | −0.0014 |
| Strong ~8×, some spread | 0.0640 | 0.0048 | −0.0592 |
| Down ~0.3× | 0.0017 | 0.0000 | −0.0016 |
| Noisy, borderline | 0.2790 | 0.2078 | −0.0712 |

Note the "Strong ~8×" row: the skewed fold-change inflated the variance enough to
**hide** a real effect (p=0.064, "n.s.") that the correct ΔCt test resolves as
significant (p=0.005). Significance stars on graphs/exports will change accordingly.

## Impact summary
- Point estimates change **only** when a pooled Condition has an excluded well on a
  non-first sample (otherwise identical).
- p-values shift for essentially all multi-replicate comparisons; the direction is
  toward correctly detecting effects that the skew previously masked.
- No test asserted specific p-value magnitudes, so existing tests remain valid;
  new regression tests lock in the exclusion behaviour.
