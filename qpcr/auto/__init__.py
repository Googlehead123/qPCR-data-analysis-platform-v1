"""Deterministic Auto-Analyze engine for the qPCR platform.

Pure, reproducible automation (NO external LLM/API): pre-analysis data screening,
transparent statistical-test recommendation, and rule/template-based
interpretation of results against each efficacy category's expected direction.

Every function here is a pure function of its inputs so it can be unit-tested and
gives identical output across machines. The deterministic ΔΔCt/stats pipeline
remains authoritative; this module advises, screens, and narrates — it never
mutates results.
"""

from qpcr.auto.screening import screen_data
from qpcr.auto.stats_advisor import recommend_test
from qpcr.auto.interpret import interpret_results, interpret_gene
from qpcr.auto.miqe import build_miqe_checklist

__all__ = ["screen_data", "recommend_test", "interpret_results", "interpret_gene",
           "build_miqe_checklist"]
