#!/usr/bin/env python3
"""Strip non-manuscript legacy hospital-rank diagnostics before public release.

This tool edits only released aggregate JSON files. It never reads controlled
source data and preserves all primary and controlled-thinning quantities used
by the current manuscript.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
THINNING_FILES = tuple(
    RESULTS / "thinning" / f"eicu_v4_controlled_thinning_{hours}h.json"
    for hours in (12, 24, 36, 48)
)
AUDIT_FILE = RESULTS / "v6_targeted_reanalysis_audit.json"
RANKING_REFERENCE_KEYS = {
    "hospital_count_for_rank_diagnostic",
    "minimum_reference_episodes_per_rank_hospital",
}
RANKING_EFFECT_KEYS = {
    "hospital_raw_vs_thinned_failure_rank_spearman",
    "hospital_changed_failure_quartile_proportion",
}


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: dict[str, object]) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def legacy_ranking_keys() -> list[str]:
    found: list[str] = []
    for path in THINNING_FILES:
        payload = read_json(path)
        for section, keys in (
            ("reference_cohort", RANKING_REFERENCE_KEYS),
            ("effect_decomposition", RANKING_EFFECT_KEYS),
            ("interpretation", {"hospital_rank"}),
        ):
            present = set(payload.get(section, {})).intersection(keys)
            found.extend(f"{path.name}:{section}.{key}" for key in sorted(present))
    audit = read_json(AUDIT_FILE)
    for hours, effects in audit.get("controlled_thinning", {}).items():
        present = set(effects).intersection(RANKING_EFFECT_KEYS)
        found.extend(f"{AUDIT_FILE.name}:controlled_thinning.{hours}.{key}" for key in sorted(present))
    return found


def strip_legacy_ranking() -> int:
    removed = 0
    for path in THINNING_FILES:
        payload = read_json(path)
        for section, keys in (
            ("reference_cohort", RANKING_REFERENCE_KEYS),
            ("effect_decomposition", RANKING_EFFECT_KEYS),
            ("interpretation", {"hospital_rank"}),
        ):
            container = payload.get(section, {})
            for key in keys:
                if key in container:
                    container.pop(key)
                    removed += 1
        write_json(path, payload)
    audit = read_json(AUDIT_FILE)
    for effects in audit.get("controlled_thinning", {}).values():
        for key in RANKING_EFFECT_KEYS:
            if key in effects:
                effects.pop(key)
                removed += 1
    write_json(AUDIT_FILE, audit)
    return removed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="remove legacy ranking diagnostics")
    parser.add_argument("--check", action="store_true", help="fail if a legacy ranking diagnostic remains")
    args = parser.parse_args()
    if args.write == args.check:
        parser.error("specify exactly one of --write or --check")
    if args.write:
        print(json.dumps({"removed": strip_legacy_ranking()}, indent=2))
        return
    found = legacy_ranking_keys()
    if found:
        raise SystemExit("legacy hospital-ranking diagnostics remain:\n" + "\n".join(found))
    print(json.dumps({"status": "PASS", "legacy_hospital_ranking_diagnostics": 0}, indent=2))


if __name__ == "__main__":
    main()
