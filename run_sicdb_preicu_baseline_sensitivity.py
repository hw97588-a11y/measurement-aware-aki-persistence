#!/usr/bin/env python3
"""SICdb ICU-pre-admission creatinine baseline sensitivity analysis."""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

from run_interval_aki_primary import first_interval_episode, load_sicdb
from run_mimic_baseline_sensitivity import summarize_strategy


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    spells, labs = load_sicdb()
    rolling = {
        case_id: episode
        for case_id, spell in spells.items()
        if (episode := first_interval_episode("SICdb v1.0.8", spell, labs.get(case_id, []))) is not None
    }
    preicu_medians = {
        case_id: statistics.median([value for time, value in values if time < 0])
        for case_id, values in labs.items()
        if any(time < 0 for time, _ in values)
    }
    output = {
        "source": "SICdb v1.0.8",
        "scope": "ICU-pre-admission measured-creatinine baseline sensitivity only; no outcomes or causal models",
        "rolling_reference": {"first_aki_transition_episodes": len(rolling), "episodes_with_48h_potential": sum(episode.potential48 for episode in rolling.values())},
        "preicu_median_strategy": summarize_strategy("preicu_measured_creatinine_median", preicu_medians, spells, labs, rolling),
        "interpretation_constraint": "The source provides ICU-pre-admission laboratory offsets, but their clinical setting can include perioperative measurements; this is a transport sensitivity rather than a preferred universal baseline.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / "sicdb_preicu_baseline_sensitivity_summary.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
