#!/usr/bin/env python3
"""Compare MIMIC hospital-observable versus ICU-only phenotype windows."""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

from run_interval_aki_primary import load_mimic, summarize


def compact(summary: dict[str, object]) -> dict[str, object]:
    return {
        "cohort": summary["cohort"],
        "primary_population_categories": summary["primary_population_categories"],
        "primary_not_definitely_classifiable": summary["primary_not_definitely_classifiable"],
        "primary_interval_indeterminate": summary["primary_interval_indeterminate"],
        "primary_right_censored_unresolved": summary["primary_right_censored_unresolved"],
        "interval_widths": summary["interval_widths"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    spells, labs = load_mimic()
    full, _ = summarize("MIMIC-IV v3.1 hospital-observable", spells, labs)
    icu_only_spells = {
        hadm: replace(spell, end=min(spell.end, float(spell.extra["continuous_icu_end_minutes"])))
        for hadm, spell in spells.items()
    }
    icu_only, _ = summarize("MIMIC-IV v3.1 ICU-only", icu_only_spells, labs)
    output = {
        "source": "MIMIC-IV v3.1",
        "scope": "source-observation-window sensitivity only; no outcomes or causal models",
        "hospital_observable_window": compact(full),
        "icu_only_window": compact(icu_only),
        "interpretation_constraint": "A shorter ICU-only window creates additional administrative censoring and is not assumed equivalent to an absence of AKI or recovery.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / "mimic_observation_window_sensitivity_summary.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
