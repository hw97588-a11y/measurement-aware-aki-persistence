#!/usr/bin/env python3
"""NDT final-gate sensitivity to long gaps in the observed AKI-state chain.

The frozen v4 primary results are not overwritten.  This additional post hoc
analysis asks how many definite-persistent classifications still have an
observed positive-state chain with adjacent creatinine gaps no longer than 24
or 36 hours.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from interval_aki_v4_engine import (
    M24,
    M48,
    classify_first_episode,
    load_source,
)
from run_interval_aki_primary import Spell, deduplicate, wilson


def _fast_source(database: str):
    if database == "eicu":
        # Numerically reconciled with the v4 loader while avoiding a Python
        # parse of every non-creatinine row in the 2.3-Gb laboratory table.
        from run_v4_controlled_thinning import fast_eicu_source

        return fast_eicu_source()
    return load_source(database)


def persistence_support_chain(
    spell: Spell,
    values: list[tuple[float, float]],
    onset_upper: float,
    baseline: float,
    recovery_lower: float | None,
) -> list[tuple[float, float]]:
    """Return observed AKI-state values through the persistence-support time."""
    recovery_limit = min(baseline + 0.3, 1.5 * baseline)
    series = [(time, value) for time, value in deduplicate(values) if onset_upper <= time <= spell.end]
    if recovery_lower is not None:
        support_end = recovery_lower
    else:
        support_end = max((time for time, value in series if value >= recovery_limit), default=onset_upper)
    return [
        (time, value)
        for time, value in series
        if time <= support_end and value >= recovery_limit
    ]


def maximum_adjacent_gap(chain: list[tuple[float, float]]) -> float:
    """Maximum adjacent observed-state gap in minutes; zero for one value."""
    return max((chain[index][0] - chain[index - 1][0] for index in range(1, len(chain))), default=0.0)


def analyse(database: str) -> dict[str, object]:
    source, spells, labs = _fast_source(database)
    categories: Counter[str] = Counter()
    persistent_max_gaps: list[float] = []
    persistent_chain_lengths: list[int] = []

    for identifier, spell in spells.items():
        episode = classify_first_episode(source, spell, labs.get(identifier, []))
        if episode is None or not episode.coverage_48h:
            continue
        categories[episode.category] += 1
        if episode.category != "definite_persistent":
            continue
        chain = persistence_support_chain(
            spell,
            labs.get(identifier, []),
            episode.onset_upper,
            episode.baseline,
            episode.recovery_lower,
        )
        persistent_chain_lengths.append(len(chain))
        persistent_max_gaps.append(maximum_adjacent_gap(chain))

    denominator = sum(categories.values())
    original_persistent = categories["definite_persistent"]
    original_transient = categories["definite_transient"]
    original_indeterminate = categories["interval_indeterminate"] + categories["right_censored_unresolved"]
    if len(persistent_max_gaps) != original_persistent:
        raise AssertionError("Every definite-persistent episode must have one chain diagnostic")

    sensitivities = {}
    for hours in (24, 36):
        cutoff = hours * 60
        reclassified = sum(gap > cutoff for gap in persistent_max_gaps)
        supported = original_persistent - reclassified
        indeterminate = original_indeterminate + reclassified
        sensitivities[f"max_adjacent_gap_{hours}h"] = {
            "continuity_supported_definite_persistent": supported,
            "continuity_gap_indeterminate": wilson(reclassified, denominator),
            "total_classification_indeterminate_after_reclassification": wilson(indeterminate, denominator),
            "persistent_identified_set": {
                "lower_continuity_supported_persistent": round(supported / denominator, 6),
                "upper_not_definite_transient_unchanged": round(1 - original_transient / denominator, 6),
                "width": round(indeterminate / denominator, 6),
            },
        }

    gaps = np.asarray(persistent_max_gaps, dtype=float) / 60
    lengths = np.asarray(persistent_chain_lengths, dtype=float)
    return {
        "source": source,
        "analysis_status": "Additional post hoc methodological robustness analysis; frozen v4 primary results unchanged.",
        "primary_convention": "The indexed episode is considered ongoing until the first observed recovery meeting the fixed episode-ending rule.",
        "sensitivity_interpretation": "A long gap does not prove recovery; it removes continuity support for the definite-persistent lower bound under the selected maximum-gap rule.",
        "primary_48h_icu_coverage_denominator": denominator,
        "original_categories": dict(categories),
        "original_persistent_identified_set": {
            "lower_definite_persistent": round(original_persistent / denominator, 6),
            "upper_not_definite_transient": round(1 - original_transient / denominator, 6),
            "width": round(original_indeterminate / denominator, 6),
        },
        "definite_persistent_chain_diagnostics": {
            "episodes": original_persistent,
            "observed_state_measurements_median": round(float(np.median(lengths)), 3),
            "observed_state_measurements_p25_p75": [
                round(float(np.quantile(lengths, 0.25)), 3),
                round(float(np.quantile(lengths, 0.75)), 3),
            ],
            "maximum_adjacent_gap_hours_median": round(float(np.median(gaps)), 3),
            "maximum_adjacent_gap_hours_p25_p75": [
                round(float(np.quantile(gaps, 0.25)), 3),
                round(float(np.quantile(gaps, 0.75)), 3),
            ],
            "maximum_adjacent_gap_hours_p90": round(float(np.quantile(gaps, 0.90)), 3),
        },
        "sensitivity_results": sensitivities,
        "methods_note": "All originally definite-transient, interval-indeterminate and right-censored-unresolved classifications remain unchanged. Only originally definite-persistent episodes can be reclassified.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", choices=("mimic", "sicdb", "eicu"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    output = analyse(args.database)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / f"{args.database}_ndt_continuity_gap_sensitivity.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
