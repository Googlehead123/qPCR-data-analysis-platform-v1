"""MIQE-style reporting checklist.

Builds a Markdown reporting appendix from an analysis provenance record. The
items this tool can determine (quantification method, reference gene, calibrator,
statistics, exclusions) are auto-filled; wet-lab / assay-validation items that
live outside the tool are listed as unchecked boxes for the user to complete.
Pure function — text only, no I/O.

Reference: Bustin et al., MIQE guidelines (Clin Chem 2009); MIQE 2.0 (2025).
"""

from __future__ import annotations


def build_miqe_checklist(prov: dict) -> str:
    """Return a Markdown MIQE-style checklist from a provenance dict
    (as produced by the app's build_provenance)."""
    prov = prov or {}
    cmps = ", ".join(prov.get("comparison_conditions") or []) or "—"
    auto = [
        f"- **Quantification method:** {prov.get('method', 'relative (Livak 2^-ΔΔCt)')}",
        f"- **Reference (housekeeping) gene:** {prov.get('reference_gene') or '—'}",
        f"- **Calibrator / reference condition:** {prov.get('reference_condition') or '—'}",
        f"- **Comparison condition(s):** {cmps}",
        f"- **Statistical test:** {prov.get('statistical_test', '—')}",
        f"- **Multiple-testing correction:** {prov.get('fdr_correction', 'Benjamini-Hochberg')}",
        f"- **Genes analysed:** {prov.get('n_genes', '—')}; "
        f"**biological/technical samples:** {prov.get('n_samples', '—')}",
        f"- **Excluded wells (QC):** {prov.get('excluded_wells_count', 0)} "
        "(itemised in the provenance record)",
    ]
    todo = [
        "Sample source, handling, and storage conditions",
        "RNA extraction method; yield and purity (A260:A280, A260:A230)",
        "RNA integrity (RIN / rRNA ratio) and how assessed",
        "DNase treatment; genomic-DNA contamination assessment (−RT control)",
        "Reverse transcription: enzyme, priming strategy, input amount, conditions",
        "Target information: gene symbol, amplicon location/length, primer sequences",
        "Primer concentrations, annealing Tm, and in-silico specificity",
        "qPCR chemistry (SYBR Green / hydrolysis probe), instrument, and thermal protocol",
        "Amplification efficiency per assay (standard-curve slope, R², efficiency %)",
        "Specificity confirmation (melt-curve / gel / sequencing)",
        "No-template control (NTC) and no-RT control results",
        "Number of biological vs. technical replicates (state explicitly)",
    ]
    lines = [
        "# MIQE-style reporting checklist (RT-qPCR)",
        f"_Generated: {prov.get('generated', '—')} · {prov.get('software', 'qPCR Analysis Suite')}_",
        "",
        "## Auto-documented by this analysis",
        *auto,
        "",
        "## To complete (assay / wet-lab details not captured by this tool)",
        *[f"- [ ] {item}" for item in todo],
        "",
        "_This checklist supports MIQE-compliant reporting; it does not replace "
        "assay validation. Complete the unchecked items before submission._",
    ]
    return "\n".join(lines)
