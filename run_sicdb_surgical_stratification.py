#!/usr/bin/env python3
"""Prespecified SICdb perioperative-structure transport analysis.

SICdb exposes SurgicalAdmissionType as Unknown, Urgent Surgery, Elective
Surgery, or No Surgery. It does not expose a reliable cardiac-surgery subtype
in this field, so this analysis does not infer one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_interval_aki_primary import load_sicdb
from run_remaining_statistical_modules import category_summary, first_episode


LABELS = {
    "3124": "unknown surgery status",
    "3125": "urgent surgery",
    "3126": "elective surgery",
    "3127": "no surgery",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    spells, labs = load_sicdb()
    groups: dict[str, list] = {label: [] for label in LABELS.values()}
    groups["any surgery"] = []
    for identifier, spell in spells.items():
        episode = first_episode(spell, labs.get(identifier, []), "standard")
        if episode is None or not episode.potential48:
            continue
        label = LABELS.get(str(spell.admission_type), "unmapped")
        groups.setdefault(label, []).append(episode)
        if label in {"urgent surgery", "elective surgery"}:
            groups["any surgery"].append(episode)
    output = {
        "source": "SICdb v1.0.8",
        "scope": "Perioperative-structure transport description of the frozen primary phenotype; no outcomes or causal model.",
        "surgical_admission_type_mapping": LABELS,
        "strata": {label: category_summary(episodes) for label, episodes in groups.items()},
        "interpretation_constraint": "SICdb SurgicalAdmissionType supports elective, urgent, and no-surgery strata only. Cardiac versus other surgery is not inferred from this variable.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / "sicdb_surgical_admission_stratification.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
