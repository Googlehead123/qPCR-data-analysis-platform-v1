"""Rule/template-based interpretation of ΔΔCt results.

Turns the completed, deterministic analysis into a reviewed narrative: for each
gene and treatment condition, classify the fold-change direction, whether it is
statistically significant, and — when the efficacy category defines an expected
direction — whether the result matches expectation. Produces bilingual
(English + Korean) templated sentences. Pure function; narrates existing numbers
only, never recomputes them.
"""

from __future__ import annotations

import pandas as pd

_DIR_EN = {"up": "upregulation", "down": "downregulation", "flat": "no change"}
_DIR_KO = {"up": "증가", "down": "감소", "flat": "변화 없음"}


def _direction(fc: float, tol: float = 0.05) -> str:
    if fc is None or pd.isna(fc):
        return "flat"
    if fc > 1 + tol:
        return "up"
    if fc < 1 - tol:
        return "down"
    return "flat"


def _expected_for(gene: str, expected_direction: dict | None) -> str | None:
    if not expected_direction:
        return None
    if gene in expected_direction:
        return str(expected_direction[gene]).lower()
    low = {str(k).lower(): str(v).lower() for k, v in expected_direction.items()}
    return low.get(str(gene).lower())


def interpret_gene(gene: str, gene_df: pd.DataFrame, expected: str | None = None,
                   ref_condition: str | None = None) -> dict:
    """Interpret one gene's result frame. Returns per-condition rows + narrative."""
    if gene_df is None or gene_df.empty:
        return {"gene": gene, "expected": expected, "rows": [],
                "narrative_en": f"{gene}: no data.", "narrative_ko": f"{gene}: 데이터 없음."}

    fc_col = "Fold_Change" if "Fold_Change" in gene_df.columns else "Relative_Expression"
    rows, en_parts, ko_parts = [], [], []

    for _, r in gene_df.iterrows():
        cond = str(r.get("Condition", ""))
        fc = r.get(fc_col)
        fc = float(fc) if pd.notna(fc) else None
        # Identify the reference/baseline row and skip narrating it.
        is_ref = (ref_condition is not None and cond == ref_condition) or (
            ref_condition is None and fc is not None and abs(fc - 1.0) < 1e-6
        )
        if is_ref:
            continue

        p = r.get("p_value")
        p = float(p) if pd.notna(p) else None
        sig = str(r.get("significance", "") or "")
        significant = bool(sig.strip()) or (p is not None and p < 0.05)
        direction = _direction(fc)

        if not significant:
            verdict = "n.s."
        elif expected is None:
            verdict = f"significant {direction}"
        elif direction == expected:
            verdict = "as expected"
        elif direction == "flat":
            verdict = "n.s."
        else:
            verdict = "opposite to expected"

        rows.append({
            "condition": cond, "fold_change": fc, "p_value": p,
            "significance": sig, "direction": direction, "verdict": verdict,
        })

        fc_txt = f"{fc:.2f}-fold" if fc is not None else "n/a"
        p_txt = f"p={p:.3g}" if p is not None else "p=n/a"
        exp_en = ""
        exp_ko = ""
        if expected is not None and significant:
            if verdict == "as expected":
                exp_en = f" — matches the expected {_DIR_EN[expected]}"
                exp_ko = f" (기대 방향 {_DIR_KO[expected]}과 일치)"
            elif verdict == "opposite to expected":
                exp_en = f" — OPPOSITE to the expected {_DIR_EN[expected]}"
                exp_ko = f" (기대 방향 {_DIR_KO[expected]}과 반대)"
        en_parts.append(
            f"{cond}: {fc_txt} {_DIR_EN[direction]} "
            f"({'significant, ' if significant else 'n.s., '}{p_txt}){exp_en}."
        )
        ko_parts.append(
            f"{cond}: {fc_txt} {_DIR_KO[direction]} "
            f"({'유의, ' if significant else '비유의, '}{p_txt}){exp_ko}."
        )

    head_en = f"{gene}: " + (" ".join(en_parts) if en_parts else "no treatment conditions to interpret.")
    head_ko = f"{gene}: " + (" ".join(ko_parts) if ko_parts else "해석할 처리 조건 없음.")
    return {"gene": gene, "expected": expected, "rows": rows,
            "narrative_en": head_en, "narrative_ko": head_ko}


def interpret_results(processed_data: dict, expected_direction: dict | None = None,
                      ref_condition: str | None = None) -> list[dict]:
    """Interpret every gene in a processed_data mapping. Returns a list of
    per-gene interpretation dicts (see `interpret_gene`)."""
    out = []
    for gene, gene_df in (processed_data or {}).items():
        expected = _expected_for(gene, expected_direction)
        out.append(interpret_gene(gene, gene_df, expected=expected, ref_condition=ref_condition))
    return out
