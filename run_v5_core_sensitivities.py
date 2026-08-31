#!/usr/bin/env python3
"""Corrected ICU-only robustness analyses for the NDT revision.

This script deliberately excludes mortality, inverse-observation weighting and
hospital ranking.  It uses the v4 continuous-ICU-spell phenotype, true patient
clusters in all three sources and two-stage hospital/patient resampling in
eICU.  Outputs contain aggregate results only.
"""

from __future__ import annotations

import argparse
import json
import statistics
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from interval_aki_v4_engine import (
    M48,
    M72,
    CATEGORY_ORDER,
    V4Episode,
    classify_first_episode,
    load_source,
)
from run_interval_aki_primary import M7D, Spell, deduplicate, load_mimic
from mimic_history_helper import historical_baselines
from run_v4_primary_inference import (
    _cluster_bootstrap,
    _two_stage_hospital_patient_bootstrap,
)


@dataclass(frozen=True)
class ResampleEpisode:
    category: str
    cluster_id: str
    hospital_id: str | None


def _resampling_summary(
    episodes: list[V4Episode] | list[ResampleEpisode],
    database: str,
    replicates: int,
    seed: int,
) -> dict[str, object]:
    counts = Counter(episode.category for episode in episodes)
    denominator = len(episodes)
    indeterminate = counts["interval_indeterminate"] + counts["right_censored_unresolved"]
    if database == "eicu":
        inference = _two_stage_hospital_patient_bootstrap(episodes, replicates, seed)
    else:
        inference = _cluster_bootstrap(episodes, replicates, seed)
    return {
        "episodes": denominator,
        "unique_patients": len({episode.cluster_id for episode in episodes}),
        "categories": {category: int(counts[category]) for category in CATEGORY_ORDER},
        "classification_indeterminate": {
            "n": int(indeterminate),
            "proportion": round(indeterminate / denominator, 6),
            "cluster_bootstrap_ci95": inference["percentile_ci"]["identified_set_width"],
        },
        "persistent_identified_set": {
            "lower_definite_persistent": round(counts["definite_persistent"] / denominator, 6),
            "upper_definite_persistent_plus_indeterminate": round(1 - counts["definite_transient"] / denominator, 6),
            "width": round(indeterminate / denominator, 6),
        },
    }


def _episode_map(
    source: str,
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    recovery_policy: str = "first_recovery",
) -> dict[str, V4Episode]:
    output: dict[str, V4Episode] = {}
    for identifier, spell in spells.items():
        episode = classify_first_episode(source, spell, labs.get(identifier, []), recovery_policy)
        if episode is not None:
            output[identifier] = episode
    return output


def _flow_counts(
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    episodes: dict[str, V4Episode],
) -> dict[str, object]:
    measurements = {
        identifier: len([1 for time, _ in deduplicate(labs.get(identifier, [])) if -M7D <= time <= min(M7D, spell.end)])
        for identifier, spell in spells.items()
    }
    no_measurement = sum(count == 0 for count in measurements.values())
    one_measurement = sum(count == 1 for count in measurements.values())
    at_least_two_no_episode = sum(count >= 2 and identifier not in episodes for identifier, count in measurements.items())
    if no_measurement + one_measurement + at_least_two_no_episode + len(episodes) != len(spells):
        raise AssertionError("Flow exclusions do not reconcile to the source-spell denominator")
    primary = {identifier: episode for identifier, episode in episodes.items() if episode.coverage_48h}
    return {
        "adult_first_continuous_icu_spells": len(spells),
        "unique_patients_in_source_spells": len({str(spell.extra.get("cluster_id", spell.identifier)) for spell in spells.values()}),
        "excluded_no_valid_central_creatinine": no_measurement,
        "excluded_only_one_valid_central_creatinine": one_measurement,
        "excluded_at_least_two_measurements_but_no_observed_qualifying_nonaki_to_aki_transition": at_least_two_no_episode,
        "first_observed_transition_defined_aki_episodes": len(episodes),
        "unique_patients_with_episode": len({episode.cluster_id for episode in episodes.values()}),
        "insufficient_48h_icu_observation_opportunity": sum(not episode.coverage_48h for episode in episodes.values()),
        "primary_48h_potential_icu_coverage_episodes": len(primary),
        "unique_patients_in_primary_population": len({episode.cluster_id for episode in primary.values()}),
        "interpretation": "The source cohort is restricted to first observed transition-defined creatinine AKI, not all prevalent or biological AKI.",
    }


def _median_iqr(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0, "median": None, "p25": None, "p75": None}
    array = np.asarray(values, dtype=float)
    return {
        "n": len(values),
        "median": round(float(np.median(array)), 3),
        "p25": round(float(np.quantile(array, 0.25)), 3),
        "p75": round(float(np.quantile(array, 0.75)), 3),
    }


def _age_group(age: float | None) -> str:
    """Return a reporting category that preserves source top-coding."""
    if age is None:
        return "Unknown"
    if age < 65:
        return "<65"
    if age < 75:
        return "65-74"
    if age < 85:
        return "75-84"
    if age < 90:
        return "85-89"
    return ">=90"


def _primary_characteristics(
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    rolling: dict[str, V4Episode],
) -> dict[str, object]:
    primary = {identifier: episode for identifier, episode in rolling.items() if episode.coverage_48h}
    ages = [episode.age for episode in primary.values() if episode.age is not None]
    baselines = [episode.baseline for episode in primary.values()]
    onset_days = [episode.onset_upper / (24 * 60) for episode in primary.values()]
    spell_hours = [spells[identifier].end / 60 for identifier in primary]
    measurement_counts: list[float] = []
    adjacent_intervals: list[float] = []
    retested_24h = 0
    for identifier, episode in primary.items():
        observed = [
            time for time, _ in deduplicate(labs.get(identifier, []))
            if 0 <= time <= spells[identifier].end
        ]
        measurement_counts.append(float(len(observed)))
        adjacent_intervals.extend((later - earlier) / 60 for earlier, later in zip(observed, observed[1:]))
        retested_24h += any(episode.onset_upper < time <= episode.onset_upper + 24 * 60 for time in observed)
    return {
        "episodes": len(primary),
        "unique_patients": len({episode.cluster_id for episode in primary.values()}),
        "age_years_or_top_coded_lower_bound": _median_iqr([float(value) for value in ages]),
        "age_group_distribution": dict(Counter(_age_group(episode.age) for episode in primary.values())),
        "sex": dict(Counter((episode.sex or "Unknown") for episode in primary.values())),
        "admission_type_top10": dict(Counter((episode.admission_type or "Unknown") for episode in primary.values()).most_common(10)),
        "baseline_creatinine_mg_dl": _median_iqr(baselines),
        "initial_aki_stage": dict(Counter(str(episode.initial_aki_stage) for episode in primary.values())),
        "onset_icu_day": _median_iqr(onset_days),
        "database_covered_critical_care_spell_hours": _median_iqr(spell_hours),
        "creatinine_measurements_during_covered_spell": _median_iqr(measurement_counts),
        "adjacent_creatinine_interval_hours": _median_iqr(adjacent_intervals),
        "retested_within_24h_after_first_positive": {
            "n": retested_24h,
            "denominator": len(primary),
            "proportion": round(retested_24h / len(primary), 6) if primary else None,
        },
    }


def _feasibility_summary(rolling: dict[str, V4Episode]) -> dict[str, object]:
    all_episodes = list(rolling.values())
    classifiable = sum(
        episode.coverage_48h and episode.category in {"definite_transient", "definite_persistent"}
        for episode in all_episodes
    )
    indeterminate = sum(
        episode.coverage_48h and episode.category in {"interval_indeterminate", "right_censored_unresolved"}
        for episode in all_episodes
    )
    insufficient = sum(not episode.coverage_48h for episode in all_episodes)
    denominator = len(all_episodes)
    if classifiable + indeterminate + insufficient != denominator:
        raise AssertionError("Phenotype-feasibility categories do not reconcile")
    return {
        "all_first_observed_index_aki_episodes": denominator,
        "with_48h_opportunity_and_uniquely_classifiable": {
            "n": classifiable,
            "proportion": round(classifiable / denominator, 6) if denominator else None,
        },
        "with_48h_opportunity_but_classification_indeterminate": {
            "n": indeterminate,
            "proportion": round(indeterminate / denominator, 6) if denominator else None,
        },
        "insufficient_48h_icu_observation_opportunity": {
            "n": insufficient,
            "proportion": round(insufficient / denominator, 6) if denominator else None,
        },
    }


def _one_episode_per_patient(primary: dict[str, V4Episode]) -> list[V4Episode]:
    """Retain one deterministic eligible source record per unique patient."""
    selected: dict[str, tuple[str, V4Episode]] = {}
    for identifier, episode in primary.items():
        current = selected.get(episode.cluster_id)
        if current is None or str(identifier) < current[0]:
            selected[episode.cluster_id] = (str(identifier), episode)
    return [item[1] for item in selected.values()]


def _fixed_baseline_episode(
    spell: Spell,
    values: list[tuple[float, float]],
    baseline: float,
) -> ResampleEpisode | None:
    series = [(time, value) for time, value in deduplicate(values) if -M7D <= time <= spell.end]
    index_search_end = min(M7D, float(spell.extra.get("index_search_end_minutes", spell.end)))
    previous_positive = False
    last_non_aki: float | None = None
    for index, (time, creatinine) in enumerate(series):
        positive = creatinine - baseline >= 0.3 or creatinine / baseline >= 1.5
        if positive and not previous_positive and 0 <= time <= index_search_end and last_non_aki is not None:
            recovery_limit = min(baseline + 0.3, 1.5 * baseline)
            post = series[index:]
            recovery_index = next((j for j, (_, value) in enumerate(post) if value < recovery_limit), None)
            if recovery_index is not None:
                recovery_upper = post[recovery_index][0]
                recovery_lower = post[recovery_index - 1][0] if recovery_index else time
                duration_lower = max(0.0, recovery_lower - time)
                duration_upper = recovery_upper - last_non_aki
                category = (
                    "definite_transient"
                    if duration_upper <= M48
                    else "definite_persistent"
                    if duration_lower > M48
                    else "interval_indeterminate"
                )
            else:
                last_positive = max((observed for observed, value in post if value >= recovery_limit), default=time)
                duration_lower = max(0.0, last_positive - time)
                category = "definite_persistent" if duration_lower > M48 else "right_censored_unresolved"
            if spell.end < time + M48:
                return None
            return ResampleEpisode(
                category=category,
                cluster_id=str(spell.extra.get("cluster_id", spell.identifier)),
                hospital_id=str(spell.hospital) if spell.hospital is not None else None,
            )
        if not positive:
            last_non_aki = time
        previous_positive = positive
    return None


def _fixed_baseline_map(
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    baselines: dict[str, float],
) -> dict[str, ResampleEpisode]:
    output: dict[str, ResampleEpisode] = {}
    for identifier, baseline in baselines.items():
        if identifier not in spells:
            continue
        episode = _fixed_baseline_episode(spells[identifier], labs.get(identifier, []), baseline)
        if episode is not None:
            output[identifier] = episode
    return output


def _mimic_baseline_sensitivity(
    source: str,
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    rolling: dict[str, V4Episode],
    replicates: int,
    seed: int,
) -> dict[str, object]:
    historical = historical_baselines(spells)
    strategies = {
        "measured_preadmission_7to365d_median": {
            identifier: statistics.median(values) for identifier, values in historical.items() if values
        },
        "first_24h_minimum": {
            identifier: min(value for time, value in values if 0 <= time <= 24 * 60)
            for identifier, values in labs.items()
            if any(0 <= time <= 24 * 60 for time, _ in values)
        },
        "first_48h_minimum": {
            identifier: min(value for time, value in values if 0 <= time <= 48 * 60)
            for identifier, values in labs.items()
            if any(0 <= time <= 48 * 60 for time, _ in values)
        },
    }
    rolling_primary = {identifier: episode for identifier, episode in rolling.items() if episode.coverage_48h}
    output = []
    for index, (name, baselines) in enumerate(strategies.items()):
        alternate = _fixed_baseline_map(spells, labs, baselines)
        shared = set(rolling_primary).intersection(alternate)
        output.append({
            "strategy": name,
            "spells_with_strategy_baseline": len(baselines),
            "summary": _resampling_summary(list(alternate.values()), "mimic", replicates, seed + index + 20),
            "shared_primary_episodes": len(shared),
            "category_concordance_among_shared_primary_episodes": round(
                sum(rolling_primary[identifier].category == alternate[identifier].category for identifier in shared) / len(shared),
                6,
            ) if shared else None,
        })
    return {"source": source, "strategies": output}


def _sicdb_preicu_baseline_sensitivity(
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    rolling: dict[str, V4Episode],
    replicates: int,
    seed: int,
) -> dict[str, object]:
    baselines = {
        identifier: statistics.median([value for time, value in values if time < 0])
        for identifier, values in labs.items()
        if any(time < 0 for time, _ in values)
    }
    alternate = _fixed_baseline_map(spells, labs, baselines)
    rolling_primary = {identifier: episode for identifier, episode in rolling.items() if episode.coverage_48h}
    shared = set(rolling_primary).intersection(alternate)
    return {
        "spells_with_preicu_measured_baseline": len(baselines),
        "summary": _resampling_summary(list(alternate.values()), "sicdb", replicates, seed + 30),
        "shared_primary_episodes": len(shared),
        "category_concordance_among_shared_primary_episodes": round(
            sum(rolling_primary[identifier].category == alternate[identifier].category for identifier in shared) / len(shared),
            6,
        ) if shared else None,
    }


def _mimic_window_sensitivity(replicates: int, seed: int) -> dict[str, object]:
    raw_spells, labs = load_mimic()
    hospital_spells = {
        identifier: replace(
            spell,
            extra={
                **spell.extra,
                "cluster_id": str(spell.extra["subject_id"]),
                "index_search_end_minutes": float(spell.extra["continuous_icu_end_minutes"]),
            },
        )
        for identifier, spell in raw_spells.items()
    }
    icu_spells = {
        identifier: replace(spell, end=float(spell.extra["continuous_icu_end_minutes"]))
        for identifier, spell in hospital_spells.items()
        if float(spell.extra["continuous_icu_end_minutes"]) > 0
    }
    hospital_episodes = _episode_map("MIMIC-IV v3.1 hospital-wide", hospital_spells, labs)
    icu_episodes = _episode_map("MIMIC-IV v3.1 ICU-only", icu_spells, labs)
    return {
        "icu_only": _resampling_summary(
            [episode for episode in icu_episodes.values() if episode.coverage_48h], "mimic", replicates, seed + 40
        ),
        "hospital_wide": _resampling_summary(
            [episode for episode in hospital_episodes.values() if episode.coverage_48h], "mimic", replicates, seed + 41
        ),
        "note": "ICU-only is the cross-database primary scope. The hospital-wide sensitivity keeps index-AKI onset restricted to the same continuous ICU spell and extends only recovery follow-up to hospital discharge or earlier in-hospital death.",
    }


def _thinning_reference_comparison(
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    rolling: dict[str, V4Episode],
) -> dict[str, object]:
    primary = {identifier: episode for identifier, episode in rolling.items() if episode.coverage_48h}
    reference = {
        identifier: episode
        for identifier, episode in primary.items()
        if sum(episode.onset_upper <= time <= episode.onset_upper + M72 for time, _ in deduplicate(labs.get(identifier, []))) >= 4
    }

    def describe(items: dict[str, V4Episode]) -> dict[str, object]:
        identifiers = list(items)
        baselines = np.asarray([items[identifier].baseline for identifier in identifiers], dtype=float)
        onset_days = np.asarray([items[identifier].onset_upper / (24 * 60) for identifier in identifiers], dtype=float)
        measurement_counts = np.asarray([
            sum(items[identifier].onset_upper <= time <= items[identifier].onset_upper + M72 for time, _ in deduplicate(labs.get(identifier, [])))
            for identifier in identifiers
        ], dtype=float)
        coverage_hours = np.asarray([spells[identifier].end / 60 for identifier in identifiers], dtype=float)
        admission = Counter((spells[identifier].admission_type or "Unknown") for identifier in identifiers)
        stage = Counter(items[identifier].initial_aki_stage for identifier in identifiers)

        def median_iqr(values: np.ndarray) -> dict[str, float]:
            return {
                "median": round(float(np.median(values)), 3),
                "p25": round(float(np.quantile(values, 0.25)), 3),
                "p75": round(float(np.quantile(values, 0.75)), 3),
            }

        return {
            "episodes": len(items),
            "unique_patients": len({episode.cluster_id for episode in items.values()}),
            "initial_aki_stage": {str(key): int(value) for key, value in sorted(stage.items())},
            "baseline_creatinine_mg_dl": median_iqr(baselines),
            "onset_icu_day": median_iqr(onset_days),
            "creatinine_measurements_onset_to_72h": median_iqr(measurement_counts),
            "icu_observation_duration_hours": median_iqr(coverage_hours),
            "admission_source_top10": dict(admission.most_common(10)),
        }

    return {
        "full_eicu_primary_cohort": describe(primary),
        "measurement_rich_thinning_reference_cohort": describe(reference),
        "selection_fraction": round(len(reference) / len(primary), 6),
        "constraint": "Controlled-thinning results apply to this fixed measurement-rich reference cohort and are not population failure-rate estimates.",
    }


def analyse(database: str, replicates: int, seed: int) -> dict[str, object]:
    source, spells, labs = load_source(database)
    rolling = _episode_map(source, spells, labs)
    primary = [episode for episode in rolling.values() if episode.coverage_48h]
    policies = {
        "first_recovery_primary": "first_recovery",
        "two_observed_recoveries_at_least_6h_apart": "two_recoveries_6h",
        "recovery_confirmed_24to48h_without_intervening_relapse": "confirmed_24to48h",
    }
    recovery = {}
    for index, (label, policy) in enumerate(policies.items()):
        episodes = _episode_map(source, spells, labs, policy)
        recovery[label] = _resampling_summary(
            [episode for episode in episodes.values() if episode.coverage_48h],
            database,
            replicates,
            seed + index,
        )
    strict = [
        episode
        for episode in primary
        if episode.onset_upper >= M48 and episode.n_first48h >= 2 and episode.no_identifiable_aki_first48h
    ]
    primary_map = {identifier: episode for identifier, episode in rolling.items() if episode.coverage_48h}
    stage1 = [episode for episode in primary if episode.initial_aki_stage == 1]
    stage2plus = [episode for episode in primary if episode.initial_aki_stage >= 2]
    output: dict[str, object] = {
        "source": source,
        "scope": "Corrected v5 ICU-only core robustness; mortality, IPW and hospital ranking excluded.",
        "flow": _flow_counts(spells, labs, rolling),
        "all_index_episode_phenotype_feasibility": _feasibility_summary(rolling),
        "primary_clinical_and_observation_characteristics": _primary_characteristics(spells, labs, rolling),
        "recovery_definition_sensitivity": recovery,
        "strict_icu_acquired_aki": _resampling_summary(strict, database, replicates, seed + 10),
        "initial_aki_stage_sensitivity": {
            "stage_1": _resampling_summary(stage1, database, replicates, seed + 11),
            "stage_2_or_3": _resampling_summary(stage2plus, database, replicates, seed + 12),
            "interpretation": "Descriptive clinical-gradient sensitivity; no multiplicity-adjusted subgroup claim was prespecified.",
        },
        "one_episode_per_unique_patient_sensitivity": {
            **_resampling_summary(_one_episode_per_patient(primary_map), database, replicates, seed + 13),
            "selection_rule": "One deterministic eligible source record per unique patient, selected by stable source identifier; this addresses repeated contribution but is not interpreted as chronologic ordering where absolute dates are unavailable.",
        },
    }
    if database == "mimic":
        output["baseline_creatinine_sensitivity"] = _mimic_baseline_sensitivity(
            source, spells, labs, rolling, replicates, seed
        )
        output["observation_window_sensitivity"] = _mimic_window_sensitivity(replicates, seed)
    elif database == "sicdb":
        output["preicu_measured_baseline_sensitivity"] = _sicdb_preicu_baseline_sensitivity(
            spells, labs, rolling, replicates, seed
        )
    elif database == "eicu":
        output["controlled_thinning_reference_comparison"] = _thinning_reference_comparison(spells, labs, rolling)
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", choices=("mimic", "sicdb", "eicu"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260830)
    args = parser.parse_args()
    output = analyse(args.database, args.bootstrap_replicates, args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / f"{args.database}_v5_core_sensitivities.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
