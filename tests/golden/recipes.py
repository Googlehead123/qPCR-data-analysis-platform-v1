"""Recipes: one scripted pass through the app, from CSV upload to every export.

Synthetic recipes and their CSVs are generated here (no client data) and run in
CI. Local recipes are JSON files under ``<test-data root>/golden/recipes`` and
are never committed.
"""

import hashlib
import json
from dataclasses import asdict, dataclass, field

from tests.golden import paths


def opaque_id(name: str) -> str:
    """What output shows instead of a local recipe's name (a client study code)."""
    return hashlib.sha256(name.encode("utf-8")).hexdigest()[:8]


def _app_default_ttest_type() -> str:
    import conftest  # noqa: F401 - installs the Streamlit mock before qpcr imports it
    from qpcr.analysis import DEFAULT_TTEST_TYPE

    return DEFAULT_TTEST_TYPE


@dataclass(kw_only=True)
class Recipe:
    name: str
    # Paths relative to the test-data root, or names in SYNTHETIC_CSVS.
    csv_files: list[str]
    # None = use the app's suggestion / detection.
    efficacy: str | None = None
    hk_gene: str | None = None
    conditions: dict[str, str] = field(default_factory=dict)
    include: dict[str, bool] = field(default_factory=dict)
    # None = keep the natural order the upload sets.
    sample_order: list[str] | None = None
    ref_sample: str
    cmp_samples: list[str | None]
    ttest_type: str = field(default_factory=_app_default_ttest_type)
    # None = every auto-QC decision overridden (all replicates restored).
    auto_qc_threshold: float | None = 0.3
    # Each entry is [gene, sample, well].
    manual_exclusions: list[list[str]] = field(default_factory=list)
    concentration: str = ""
    gene_display_names: dict[str, str] = field(default_factory=dict)
    # Keys exactly as resolve_gene_settings reads them.
    graph_settings: dict = field(default_factory=dict)
    image_export: dict = field(
        default_factory=lambda: {"fmt": "PNG", "width": 1200, "height": 800}
    )

    def __post_init__(self):
        if len(self.cmp_samples) != 3:
            raise ValueError(f"{self.label}: cmp_samples must have length 3")
        if any(len(e) != 3 for e in self.manual_exclusions):
            raise ValueError(f"{self.label}: manual_exclusions entries are [gene, sample, well]")

    @property
    def is_local(self) -> bool:
        """True when any CSV comes from the test-data root rather than this module."""
        return any(f not in SYNTHETIC_CSVS for f in self.csv_files)

    @property
    def label(self) -> str:
        """The name for a synthetic recipe; an opaque id for a local one."""
        return f"local {opaque_id(self.name)}" if self.is_local else self.name

    def to_json(self) -> str:
        return json.dumps(asdict(self), ensure_ascii=False, indent=2) + "\n"

    @classmethod
    def from_json(cls, text: str) -> "Recipe":
        return cls(**json.loads(text))


# ---- synthetic plates ---------------------------------------------------------
# Hand-picked Ct values (3 replicates each). Every triplicate sits within the
# 0.30 auto-QC SD except where a feature is planted on purpose.

_ACTIN_A = {
    "S1": (18.02, 18.10, 17.95), "S2": (18.20, 18.12, 18.25), "S3": (17.90, 18.00, 17.96),
    "S4": (18.05, 18.15, 18.08), "S5": (18.11, 18.03, 18.18), "S6": (18.30, 18.22, 18.27),
}
_COL1A1_A = {
    "S1": (25.10, 25.02, 25.18), "S2": (24.05, 23.98, 24.12),
    # Active B: every well between 35 and 38 (high-Ct flag territory).
    "S3": (36.30, 36.50, 36.70),
    "S4": (24.52, 24.60, 24.48), "S5": (24.40, 24.55, 24.47), "S6": (30.10, 30.25, 30.05),
}
_MMP1_A = {
    "S1": (22.00, 22.08, 21.95),
    # Active A: the third well sits +0.8 Ct off the other two -> trimmed at SD 0.30.
    "S2": (23.00, 23.10, 23.85),
    # Active B: B1 (first) is the manually excluded well; same ID as COL1A1's B1.
    "S3": (21.40, 21.00, 21.10),
    "S4": (22.60, 22.52, 22.68), "S5": (22.45, 22.58, 22.50), "S6": (27.00, 27.20, 26.90),
}

_ACTIN_B = {"S1": (17.52, 17.60, 17.48), "S2": (17.55, 17.47, 17.62), "S10": (17.40, 17.51, 17.45)}
_COL1A1_B = {"S1": (26.02, 25.95, 26.10), "S2": (25.31, 25.25, 25.40), "S10": (24.60, 24.68, 24.55)}
_MMP1_B = {
    "S1": (23.00, 23.07, 22.95),
    # Above the 0.30 SD: auto-QC trims it on upload, and the recipe overrides that.
    "S2": (23.55, 23.62, 24.40),
    "S10": (24.20, 24.12, 24.27),
}


def _plate_rows(genes: dict[str, dict[str, tuple]]) -> list[tuple[str, str, str, float]]:
    """(well, sample, target, ct) rows; each sample owns a 3-well block per gene.

    Genes go two abreast (columns 1-6 and 7-12), two samples per row, so the
    same sample and gene slot gives the same well ID on every plate.
    """
    rows = []
    for g, (target, per_sample) in enumerate(genes.items()):
        rows_per_gene_pair = (len(per_sample) + 1) // 2
        for i, (sample, cts) in enumerate(per_sample.items()):
            row = "ABCDEFGH"[(g // 2) * rows_per_gene_pair + i // 2]
            col0 = (i % 2) * 3 + (g % 2) * 6 + 1
            for k, ct in enumerate(cts):
                rows.append((f"{row}{col0 + k}", sample, target, ct))
    return rows


def _bare_csv(genes) -> bytes:
    """The ``test_qpcr.csv`` layout: a bare Well,Sample,Target,Ct header on line 1."""
    lines = ["Well,Sample,Target,Ct"]
    lines += [f"{w},{s},{t},{ct:.2f}" for w, s, t, ct in _plate_rows(genes)]
    return ("\n".join(lines) + "\n").encode("utf-8")


def _quantstudio_csv(genes, eds_name: str) -> bytes:
    """QuantStudio export: instrument preamble, then the 'Well Position' header."""
    lines = [
        "Block Type,96-Well Block (0.2mL),,,,,,",
        f"Experiment File Name,{eds_name},,,,,,",
        "Experiment Run End Time,2026-01-01 09:00:00,,,,,,",
        ",,,,,,,",
        "Well,Well Position,Sample Name,Target Name,Task,Reporter,Quencher,CT",
    ]
    for n, (w, s, t, ct) in enumerate(_plate_rows(genes), 1):
        lines.append(f"{n},{w},{s},{t},UNKNOWN,FAM,NFQ-MGB,{ct:.3f}")
    return ("\n".join(lines) + "\n").encode("utf-8")


SYNTHETIC_CSVS = {
    "synthetic_a_plate1.csv": lambda: _bare_csv({"COL1A1": _COL1A1_A, "ACTIN": _ACTIN_A}),
    # The companion plate carries no reference gene, and reuses the well IDs.
    "synthetic_a_plate2.csv": lambda: _bare_csv({"MMP1": _MMP1_A}),
    "synthetic_b.csv": lambda: _quantstudio_csv(
        {"ACTIN": _ACTIN_B, "COL1A1": _COL1A1_B, "MMP1": _MMP1_B}, "synthetic_b.eds"
    ),
}


def synthetic_csv(name: str) -> bytes:
    return SYNTHETIC_CSVS[name]()


def synthetic_recipes() -> list[Recipe]:
    return [
        Recipe(
            name="synthetic_a",
            csv_files=["synthetic_a_plate1.csv", "synthetic_a_plate2.csv"],
            efficacy="광노화",
            hk_gene="ACTIN",
            conditions={
                "S1": "Non-treated", "S2": "Active A", "S3": "Active B",
                "S4": "Active C", "S5": "Active C", "S6": "S6",
            },
            include={"S1": True, "S2": True, "S3": True, "S4": True, "S5": True, "S6": False},
            sample_order=["S3", "S1", "S2", "S4", "S5", "S6"],
            ref_sample="S1",
            cmp_samples=["S1", "S2", "S3"],
            ttest_type="student",
            auto_qc_threshold=0.3,
            manual_exclusions=[["MMP1", "S3", "B1"]],
            concentration="0.1%",
            gene_display_names={"MMP1": "MMP-1"},
            graph_settings={"error_bar_mode": "sem", "COL1A1_y_max": 3.0},
        ),
        Recipe(
            name="synthetic_b",
            csv_files=["synthetic_b.csv"],
            conditions={"S1": "Non-treated", "S2": "Low dose", "S10": "High dose"},
            include={"S1": True, "S2": True, "S10": True},
            sample_order=["S10", "S1", "S2"],
            ref_sample="S1",
            cmp_samples=["S1", None, None],
            ttest_type="welch",
            auto_qc_threshold=None,
        ),
    ]


def load_local_recipes() -> list[Recipe]:
    """Every ``*.json`` under ``<test-data root>/golden/recipes``, by name; [] if none."""
    rdir = paths.recipes_dir()
    if rdir is None or not rdir.is_dir():
        return []
    return [Recipe.from_json(p.read_text(encoding="utf-8")) for p in sorted(rdir.glob("*.json"))]
