"""Check that the running interpreter matches the locked ``requirements.txt``.

Golden outputs depend on exact numpy / pandas / plotly / kaleido builds, so any
capture or compare must run inside the locked ``.venv``.
"""

import importlib.metadata
from pathlib import Path

from packaging.requirements import InvalidRequirement, Requirement

_MISSING = "<missing>"


def lock_mismatches(requirements: Path = Path("requirements.txt")) -> list[str]:
    """Return one ``"pkg: installed X != locked Y"`` entry per drifted pin.

    Pins whose environment marker does not apply here (``colorama`` on Linux)
    are skipped. A missing package is reported as installed ``<missing>``.
    """
    mismatches: list[str] = []
    for raw in Path(requirements).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        try:
            req = Requirement(line)
        except InvalidRequirement:
            continue
        if req.marker is not None and not req.marker.evaluate():
            continue
        specs = list(req.specifier)
        if len(specs) != 1 or specs[0].operator != "==":
            continue
        locked = specs[0].version
        try:
            installed = importlib.metadata.version(req.name)
        except importlib.metadata.PackageNotFoundError:
            installed = _MISSING
        if installed != locked:
            mismatches.append(f"{req.name}: installed {installed} != locked {locked}")
    return mismatches
