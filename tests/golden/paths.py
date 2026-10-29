"""Where golden material lives.

``test_data/`` is gitignored client material in the main checkout. Worktrees
do not have it, so the root is found through the repository's common git dir.
Nothing derived from it is ever committed; synthetic goldens live in the repo.
"""

import os
import subprocess
from pathlib import Path

SYNTH_GOLDENS = Path("tests/golden/synthetic")


def test_data_root() -> Path | None:
    """``$QPCR_TEST_DATA``, else ``<parent of git common dir>/test_data``; None if neither exists."""
    override = os.environ.get("QPCR_TEST_DATA")
    if override:
        return Path(override)
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--git-common-dir"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None
    if not out:
        return None
    candidate = Path(out).resolve().parent / "test_data"
    return candidate if candidate.is_dir() else None


def recipes_dir() -> Path | None:
    root = test_data_root()
    return root / "golden" / "recipes" if root else None


def local_out(baseline_sha: str) -> Path | None:
    root = test_data_root()
    return root / "golden" / "out" / baseline_sha if root else None
