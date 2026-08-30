#!/usr/bin/env python3
"""Prespecified MIMIC baseline-creatinine phenotype sensitivity analysis.

This script compares phenotype reclassification only.  It never reads or
models death, KRT, or another clinical outcome.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import timedelta
from pathlib import Path

from run_interval_aki_primary import (
    M48,
    M7D,
    MIMIC,
    as_datetime,
    as_float,
    deduplicate,
    first_interval_episode,
    load_mimic,
    mimic_csv,
    valid_creatinine,
)


@dataclass
class StaticEpisode:
    category: str
    potential48: bool


def static_baseline_episode(end: float, values: list[tuple[float, float]], baseline: float) -> StaticEpisode | None:
    """First ICU-onset AKI episode when a fixed, prespecified baseline is used."""
    series = deduplicate(values)
    previous_positive = False
    last_nonaki: float | None = None
    for index, (time, creatinine) in enumerate(series):
        positive = creatinine - baseline >= .3 or creatinine / baseline >= 1.5
        if positive and not previous_positive and 0 <= time <= min(M7D, end) and last_nonaki is not None:
            recovery_limit = min(baseline + .3, baseline * 1.5)
            after = [(t, c) for t, c in series[index:] if t >= time and t <= end]
            recovery_idx = next((j for j, (_, value) in enumerate(after) if value < recovery_limit), None)
            if recovery_idx is not None:
                recovery_upper, _ = after[recovery_idx]
                recovery_lower, _ = after[recovery_idx - 1] if recovery_idx else (time, creatinine)
                lower = max(0., recovery_lower - time)
                upper = recovery_upper - last_nonaki
                category = "definite_transient" if upper <= M48 else "definite_persistent" if lower > M48 else "interval_indeterminate"
            else:
                lower = max(0., after[-1][0] - time) if after else 0.
                category = "definite_persistent" if lower > M48 else "right_censored_unresolved"
            return StaticEpisode(category, end >= time + M48)
        if not positive:
            last_nonaki = time
        previous_positive = positive
    return None


def historical_baselines(spells) -> dict[str, list[float]]:
    """Extract actual same-patient creatinine 7--365 days before admission."""
    by_subject: dict[str, list[str]] = defaultdict(list)
    for hadm, spell in spells.items():
        by_subject[str(spell.extra["subject_id"])].append(hadm)
    result: dict[str, list[float]] = defaultdict(list)
    with zipfile.ZipFile(MIMIC) as archive:
        raw = gzip.GzipFile(fileobj=archive.open("mimic-iv-3.1/hosp/labevents.csv.gz"))
        next(raw, b"")
        for line in raw:
            if b",50912," not in line and b",52546," not in line:
                continue
            row = next(csv.reader([line.decode("utf-8", errors="replace")]))
            if len(row) < 11 or row[10].strip().casefold() != "mg/dl":
                continue
            observed, value = as_datetime(row[6]), as_float(row[9])
            if observed is None or not valid_creatinine(value):
                continue
            for hadm in by_subject.get(row[1], []):
                admit = spells[hadm].extra["admit_dt"]
                if admit - timedelta(days=365) <= observed < admit - timedelta(days=7):
                    result[hadm].append(value)
        raw.close()
    return result


def summarize_strategy(name: str, baselines: dict[str, float], spells, labs, rolling) -> dict[str, object]:
    alternate: dict[str, StaticEpisode] = {}
    for hadm, baseline in baselines.items():
        episode = static_baseline_episode(spells[hadm].end, labs.get(hadm, []), baseline)
        if episode is not None:
            alternate[hadm] = episode
    categories = Counter(episode.category for episode in alternate.values())
    potential = {hadm: episode for hadm, episode in alternate.items() if episode.potential48}
    both = set(rolling).intersection(alternate)
    both_potential = {hadm for hadm in both if rolling[hadm].potential48 and alternate[hadm].potential48}
    concordant = sum(rolling[hadm].category == alternate[hadm].category for hadm in both_potential)
    return {
        "strategy": name,
        "spells_with_strategy_baseline": len(baselines),
        "first_aki_transition_episodes": len(alternate),
        "episodes_with_48h_potential": len(potential),
        "categories_all": dict(categories),
        "rolling_only_transition": len(set(rolling).difference(alternate)),
        "strategy_only_transition": len(set(alternate).difference(rolling)),
        "detected_by_both": len(both),
        "both_with_48h_potential": len(both_potential),
        "category_concordant_in_both_potential": concordant,
        "category_concordance_proportion": round(concordant / len(both_potential), 5) if both_potential else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    spells, labs = load_mimic()
    rolling = {
        hadm: episode
        for hadm, spell in spells.items()
        if (episode := first_interval_episode("MIMIC-IV v3.1", spell, labs.get(hadm, []))) is not None
    }
    historical = historical_baselines(spells)
    strategies = {
        "actual_preadmission_7to365d_median": {
            hadm: statistics.median(values) for hadm, values in historical.items() if values
        },
        "early_hospital_24h_minimum": {
            hadm: min(value for time, value in values if 0 <= time <= 24 * 60)
            for hadm, values in labs.items()
            if any(0 <= time <= 24 * 60 for time, _ in values)
        },
        "early_hospital_48h_minimum": {
            hadm: min(value for time, value in values if 0 <= time <= 48 * 60)
            for hadm, values in labs.items()
            if any(0 <= time <= 48 * 60 for time, _ in values)
        },
    }
    output = {
        "source": "MIMIC-IV v3.1",
        "scope": "baseline-definition sensitivity only; no outcomes, causal models, or selection based on clinical results",
        "rolling_reference": {"first_aki_transition_episodes": len(rolling), "episodes_with_48h_potential": sum(episode.potential48 for episode in rolling.values())},
        "strategies": [summarize_strategy(name, baselines, spells, labs, rolling) for name, baselines in strategies.items()],
        "interpretation_constraint": "early-hospital minima can include measurements after unobserved biological onset; they are a prespecified sensitivity definition, not a preferred causal baseline",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / "mimic_baseline_sensitivity_summary.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
