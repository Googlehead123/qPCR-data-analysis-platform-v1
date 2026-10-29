"""Capture and compare golden outputs.

    python -m tests.golden.capture [--update] [--local] [--recipe NAME]
                                   [--init-local] [--allow ALLOWLIST.json ...]

Without ``--update`` it compares and exits 1 on any change the allowlists do
not permit. Synthetic goldens live in ``tests/golden/synthetic``; local ones
(client plates, real rendering) only under the test-data root, whose diffs are
written to ``<root>/golden/diffs/`` and never printed.
"""

import argparse
import difflib
import fnmatch
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from tests.golden import env, paths, recipes
from tests.golden.recipes import Recipe

REPO = Path(__file__).resolve().parents[2]
SYNTH_DIR = REPO / paths.SYNTH_GOLDENS
MAX_LINES_PER_PART = 40
_PLATE_SUFFIX = re.compile(r"(?:-\d+|_data| ?\d+)$")


# ---- golden files ---------------------------------------------------------------


def part_filename(key: str, part: str) -> str:
    return f"{key}__{part.replace('/', '__')}.txt"


def _split_filename(name: str) -> tuple[str, str]:
    key, _, rest = name[: -len(".txt")].partition("__")
    return key, rest.replace("__", "/")


def _is_under(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def write_golden(recipe: Recipe, outputs: dict[str, dict[str, str]], dest: Path) -> None:
    """Replace ``dest``'s ``*.txt`` with one file per part.

    Raises ``ValueError`` before touching anything if a recipe built from
    test-data CSVs would be written into the repository's ``tests/`` tree.
    """
    if recipe.is_local and _is_under(dest, REPO / "tests"):
        raise ValueError(f"refusing to write local recipe {recipe.name!r} output under tests/")
    files = {}
    for key, parts in outputs.items():
        for part, text in parts.items():
            name = part_filename(key, part)
            if _split_filename(name) != (key, part) or name in files:
                raise ValueError(f"part name does not round-trip as a file name: {key}/{part}")
            files[name] = text
    dest.mkdir(parents=True, exist_ok=True)
    for old in dest.glob("*.txt"):
        old.unlink()
    for name, text in files.items():
        (dest / name).write_bytes(text.encode("utf-8"))


def load_golden(src: Path) -> dict[str, dict[str, str]]:
    golden: dict[str, dict[str, str]] = defaultdict(dict)
    if src.is_dir():
        for path in sorted(src.glob("*.txt")):
            key, part = _split_filename(path.name)
            golden[key][part] = path.read_bytes().decode("utf-8")
    return dict(golden)


# ---- comparison -------------------------------------------------------------------


def _allowed(name: str, line: str, allow: list[dict]) -> bool:
    return any(
        fnmatch.fnmatchcase(name, entry["part"]) and re.search(entry["line"], line)
        for entry in allow
    )


def compare(actual, golden, allow: list[dict]) -> list[str]:
    """Every changed line no allowlist entry permits, as ``"<key>/<part>: -|+<line>"``.

    An entry ``{"part": "<glob>", "line": "<regex>"}`` permits a removed or added
    line only when the glob matches ``<key>/<part>`` and the regex matches the
    line. A part that appears or disappears is always a violation.
    """
    violations = []
    for key in sorted(set(actual) | set(golden)):
        a_parts, g_parts = actual.get(key, {}), golden.get(key, {})
        for part in sorted(set(a_parts) | set(g_parts)):
            name = f"{key}/{part}"
            if part not in a_parts:
                violations.append(f"{name}: part missing")
                continue
            if part not in g_parts:
                violations.append(f"{name}: part added")
                continue
            if a_parts[part] == g_parts[part]:
                continue
            g_lines, a_lines = g_parts[part].splitlines(), a_parts[part].splitlines()
            matcher = difflib.SequenceMatcher(None, g_lines, a_lines, autojunk=False)
            for tag, i1, i2, j1, j2 in matcher.get_opcodes():
                if tag == "equal":
                    continue
                for line in g_lines[i1:i2]:
                    if not _allowed(name, line, allow):
                        violations.append(f"{name}: -{line}")
                for line in a_lines[j1:j2]:
                    if not _allowed(name, line, allow):
                        violations.append(f"{name}: +{line}")
    return violations


def _by_part(violations: list[str]) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = defaultdict(list)
    for v in violations:
        grouped[v.split(": ", 1)[0]].append(v)
    return grouped


# ---- local (test-data) material -------------------------------------------------


def _golden_root() -> Path | None:
    root = paths.test_data_root()
    return root / "golden" if root is not None and root.is_dir() else None


def local_baseline_dir() -> Path | None:
    """``local_out(<sha in <root>/golden/BASELINE>)``, or None if there is no baseline."""
    groot = _golden_root()
    marker = groot / "BASELINE" if groot else None
    if marker is None or not marker.is_file():
        return None
    return paths.local_out(marker.read_text().strip())


def local_skip_reason() -> str | None:
    """Why the local goldens cannot run here, or None if they can."""
    root = paths.test_data_root()
    if root is None or not root.is_dir():
        where = root if root is not None else "<git common dir>/test_data"
        return f"test_data root not found ({where}); set QPCR_TEST_DATA"
    if not recipes.load_local_recipes():
        return "no local recipes under test_data/golden/recipes; run capture --init-local"
    baseline = local_baseline_dir()
    if baseline is None or not baseline.is_dir():
        return "no local baseline (test_data/golden/BASELINE); run capture --update --local"
    from qpcr.export_utils import _find_system_browser

    if _find_system_browser() is None:
        return "Chrome not found; local goldens need real Kaleido rendering"
    return None


def write_local_diff(recipe: Recipe, violations: list[str]) -> Path:
    out = _golden_root() / "diffs" / f"{recipe.name}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(v + "\n" for v in violations), encoding="utf-8")
    return out


def plate_set_stem(stem: str) -> str:
    """Strip trailing plate numbers / ``_data`` until none is left: 'X 1', 'X2_data-1' -> 'X'."""
    while True:
        shorter = _PLATE_SUFFIX.sub("", stem)
        if shorter == stem or not shorter.strip():
            return stem.strip()
        stem = shorter


def _merge_override(base: dict, override: dict) -> dict:
    merged = dict(base)
    for field_name, value in override.items():
        if field_name not in merged:
            raise ValueError(f"recipe_overrides.json: unknown recipe field {field_name!r}")
        if isinstance(merged[field_name], dict) and isinstance(value, dict):
            merged[field_name] = {**merged[field_name], **value}
        else:
            merged[field_name] = value
    return merged


def init_local() -> int:
    from tests.golden import driver

    groot = _golden_root()
    if groot is None:
        print(local_skip_reason())
        return 2
    root = groot.parent
    groups: dict[str, list[str]] = defaultdict(list)
    for csv_path in sorted(root.glob("*.csv")):
        groups[plate_set_stem(csv_path.stem)].append(csv_path.name)
    overrides_path = groot / "recipe_overrides.json"
    overrides = json.loads(overrides_path.read_text(encoding="utf-8")) if overrides_path.is_file() else {}

    rdir = paths.recipes_dir()
    rdir.mkdir(parents=True, exist_ok=True)
    written = applied = no_hk = 0
    for stem, files in sorted(groups.items()):
        found = driver.detect_defaults(files)
        ref = found["samples"][0]
        base = Recipe(
            name=stem, csv_files=files, efficacy=found["efficacy"], hk_gene=found["hk_gene"],
            conditions={s: s for s in found["samples"]},
            include={s: True for s in found["samples"]},
            ref_sample=ref, cmp_samples=[ref, None, None], ttest_type=found["ttest_type"],
        )
        fields = json.loads(base.to_json())
        if stem in overrides:
            fields = _merge_override(fields, overrides[stem])
            applied += 1
        recipe = Recipe(**fields)
        no_hk += recipe.hk_gene is None
        (rdir / f"{recipe.name}.json").write_text(recipe.to_json(), encoding="utf-8")
        written += 1
    unused = len(set(overrides) - set(groups))
    # Counts only: recipe names are client study codes.
    print(f"init-local: {written} recipe(s) from {sum(map(len, groups.values()))} CSV(s) "
          f"written to {rdir}; {applied} override set(s) applied, {unused} unused; "
          f"{no_hk} need an hk_gene override")
    return 0


# ---- CLI ------------------------------------------------------------------------


def _head_sha() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True, cwd=str(REPO)
    ).stdout.strip()


def _load_allow(files: list[str]) -> list[dict]:
    allow = []
    for f in files:
        entries = json.loads(Path(f).read_text(encoding="utf-8"))
        for e in entries:
            if not {"part", "line"} <= set(e) <= {"part", "line", "note"}:
                raise ValueError(f"{f}: allowlist entries are {{part, line[, note]}}: {e}")
            re.compile(e["line"])
        allow.extend(entries)
    return allow


def _update(recs: list[Recipe], local: bool) -> int:
    from tests.golden import driver

    if local:
        sha = _head_sha()
        out_root = paths.local_out(sha)
        for r in recs:
            write_golden(r, driver.produce_outputs(r, real_render=True), out_root / r.name)
        (_golden_root() / "BASELINE").write_text(sha + "\n")
        print(f"update --local: {len(recs)} recipe(s) captured to {out_root}")
        return 0
    for r in recs:
        write_golden(r, driver.produce_outputs(r, real_render=False), SYNTH_DIR / r.name)
        print(f"update: {r.name} -> {(SYNTH_DIR / r.name).relative_to(REPO)}")
    return 0


def _compare(recs: list[Recipe], local: bool, allow: list[dict]) -> int:
    from tests.golden import driver

    base = local_baseline_dir() if local else SYNTH_DIR
    if base is None:
        print(local_skip_reason())
        return 2
    failed = 0
    for i, r in enumerate(recs, 1):
        golden = load_golden(base / r.name)
        violations = compare(driver.produce_outputs(r, real_render=local), golden, allow)
        if not golden:
            violations = ["<recipe>: no goldens captured"] + violations
        failed += bool(violations)
        if local:
            if violations:
                path = write_local_diff(r, violations)
                print(f"[local #{i}] {len(violations)} changed line(s) in "
                      f"{len(_by_part(violations))} part(s); diff in {path.parent}")
            else:
                print(f"[local #{i}] OK")
            continue
        if not violations:
            print(f"{r.name}: OK")
            continue
        print(f"{r.name}: {len(violations)} changed line(s)")
        for part, lines in sorted(_by_part(violations).items()):
            print(f"  {part}: {len(lines)}")
            for line in lines[:MAX_LINES_PER_PART]:
                print(f"    {line.split(': ', 1)[1]}")
            if len(lines) > MAX_LINES_PER_PART:
                print(f"    ... {len(lines) - MAX_LINES_PER_PART} more")
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m tests.golden.capture")
    parser.add_argument("--update", action="store_true", help="write goldens instead of comparing")
    parser.add_argument("--local", action="store_true", help="test-data recipes, real rendering")
    parser.add_argument("--recipe", help="only this recipe")
    parser.add_argument("--init-local", action="store_true", help="(re)generate local recipes")
    parser.add_argument("--allow", nargs="+", action="extend", default=[], metavar="ALLOWLIST.json")
    args = parser.parse_args(argv)

    mismatches = env.lock_mismatches(REPO / "requirements.txt")
    if mismatches:
        print("Installed packages differ from requirements.txt (use .venv):")
        print("\n".join(f"  {m}" for m in mismatches))
        return 2
    if args.init_local:
        return init_local()

    if args.local:
        if _golden_root() is None:
            print(local_skip_reason())
            return 2
        if args.update and args.recipe:
            print("--update --local captures every local recipe; drop --recipe")
            return 2
        recs = recipes.load_local_recipes()
    else:
        recs = recipes.synthetic_recipes()
    if args.recipe:
        recs = [r for r in recs if r.name == args.recipe]
    if not recs:
        print("no recipes selected")
        return 2
    if args.update:
        return _update(recs, args.local)
    return _compare(recs, args.local, _load_allow(args.allow))


if __name__ == "__main__":
    sys.exit(main())
