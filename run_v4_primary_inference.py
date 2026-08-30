#!/usr/bin/env python3
"""ICU-coverage primary phenotype, partial identification, and bootstrap inference."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from interval_aki_v4_engine import M24, M48, V4Episode, derive_episodes, load_source, primary_population, primary_summary


def _parameters(episodes: list[V4Episode]) -> tuple[int, int, int]:
    n = len(episodes)
    lower = sum(episode.category == "definite_persistent" for episode in episodes)
    upper = sum(episode.category != "definite_transient" for episode in episodes)
    return n, lower, upper


def _summarize_bootstrap(lower: np.ndarray, upper: np.ndarray, point_lower: float, point_upper: float) -> dict[str, object]:
    width = upper - lower
    deviations = np.maximum(np.abs(lower - point_lower), np.abs(upper - point_upper))
    critical = float(np.quantile(deviations, .975))
    return {
        "replicates": int(len(lower)),
        "percentile_ci": {
            "persistent_lower": [round(float(np.quantile(lower, .025)), 6), round(float(np.quantile(lower, .975)), 6)],
            "persistent_upper": [round(float(np.quantile(upper, .025)), 6), round(float(np.quantile(upper, .975)), 6)],
            "identified_set_width": [round(float(np.quantile(width, .025)), 6), round(float(np.quantile(width, .975)), 6)],
        },
        "simultaneous_bootstrap_confidence_region": {
            "method": "Conservative unstudentized max-deviation percentile bootstrap region",
            "critical_absolute_probability_deviation": round(critical, 6),
            "persistent_lower": [round(max(0.0, point_lower - critical), 6), round(min(1.0, point_lower + critical), 6)],
            "persistent_upper": [round(max(0.0, point_upper - critical), 6), round(min(1.0, point_upper + critical), 6)],
        },
    }


def _cluster_bootstrap(episodes: list[V4Episode], replicates: int, seed: int) -> dict[str, object]:
    """Patient/episode-cluster bootstrap for MIMIC and SICdb."""
    clusters: dict[str, list[V4Episode]] = defaultdict(list)
    for episode in episodes:
        clusters[episode.cluster_id].append(episode)
    arrays = [
        np.asarray([(1, int(item.category == "definite_persistent"), int(item.category != "definite_transient")) for item in group], dtype=int)
        for group in clusters.values()
    ]
    total = np.asarray([array.sum(axis=0) for array in arrays], dtype=int)
    n_units = len(arrays)
    rng = np.random.default_rng(seed)
    lower = np.empty(replicates)
    upper = np.empty(replicates)
    for index in range(replicates):
        selection = rng.integers(0, n_units, n_units)
        counts = total[selection].sum(axis=0)
        lower[index] = counts[1] / counts[0]
        upper[index] = counts[2] / counts[0]
    n, lo_count, hi_count = _parameters(episodes)
    output = _summarize_bootstrap(lower, upper, lo_count / n, hi_count / n)
    output["resampling"] = {"level": "patient/episode cluster", "clusters": n_units, "episodes": n}
    return output


def _two_stage_hospital_patient_bootstrap(episodes: list[V4Episode], replicates: int, seed: int) -> dict[str, object]:
    """Resample hospitals, then unique patients within each sampled hospital."""
    grouped: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for episode in episodes:
        grouped[str(episode.hospital_id)][episode.cluster_id].append(
            (int(episode.category == "definite_persistent"), int(episode.category != "definite_transient"))
        )
    hospitals = [
        [np.asarray(records, dtype=int) for records in patients.values()]
        for patients in grouped.values()
    ]
    n_hospitals = len(hospitals)
    rng = np.random.default_rng(seed)
    lower = np.empty(replicates)
    upper = np.empty(replicates)
    for index in range(replicates):
        selected_hospitals = rng.integers(0, n_hospitals, n_hospitals)
        total_n = total_lower = total_upper = 0
        for selected in selected_hospitals:
            patient_records = hospitals[selected]
            selected_patients = rng.integers(0, len(patient_records), len(patient_records))
            for patient_index in selected_patients:
                records = patient_records[patient_index]
                total_n += len(records)
                total_lower += int(records[:, 0].sum())
                total_upper += int(records[:, 1].sum())
        lower[index] = total_lower / total_n
        upper[index] = total_upper / total_n
    n, lo_count, hi_count = _parameters(episodes)
    output = _summarize_bootstrap(lower, upper, lo_count / n, hi_count / n)
    output["resampling"] = {
        "level": "hospital then unique patient within hospital; all eligible episodes for a sampled patient retained",
        "hospitals": n_hospitals,
        "unique_hospital_patient_clusters": sum(len(patients) for patients in hospitals),
        "unique_patients_global": len({episode.cluster_id for episode in episodes}),
        "episodes": n,
    }
    return output


def _threshold_curve_with_coverage(source: str, spells, labs) -> list[dict[str, object]]:
    """Threshold curve, retaining coverage end in a local identifier lookup."""
    rows = []
    # derive_episodes preserves source iteration but deliberately removes IDs;
    # recreate with a local direct call for an ID-safe aggregate summary.
    from interval_aki_v4_engine import classify_first_episode
    for hours in (24, 36, 48, 72, 96):
        episodes = []
        for identifier, spell in spells.items():
            episode = classify_first_episode(source, spell, labs.get(identifier, []), threshold_minutes=hours * 60)
            if episode is not None and spell.end >= episode.onset_upper + hours * 60:
                episodes.append(episode)
        n, lower, upper = _parameters(episodes)
        rows.append({
            "threshold_hours": hours,
            "episodes_with_potential_coverage": n,
            "persistent_lower": round(lower / n, 6) if n else None,
            "persistent_upper": round(upper / n, 6) if n else None,
            "identified_set_width": round((upper - lower) / n, 6) if n else None,
        })
    return rows


def _indeterminacy_decomposition(episodes: list[V4Episode], threshold: float = M48) -> dict[str, object]:
    output = defaultdict(int)
    onset_widths = []
    recovery_widths = []
    for episode in episodes:
        if episode.category == "right_censored_unresolved":
            output["no_observed_recovery_right_censored"] += 1
            continue
        if episode.category != "interval_indeterminate":
            continue
        assert episode.recovery_lower is not None and episode.recovery_upper is not None
        onset_contributes = episode.recovery_lower - episode.onset_upper <= threshold < episode.recovery_lower - episode.onset_lower
        recovery_contributes = episode.recovery_lower - episode.onset_upper <= threshold < episode.recovery_upper - episode.onset_upper
        if onset_contributes and recovery_contributes:
            output["both_onset_and_recovery_interval_contribute"] += 1
        elif onset_contributes:
            output["onset_interval_contributes"] += 1
        elif recovery_contributes:
            output["recovery_interval_contributes"] += 1
        else:
            output["joint_interval_crossing_not_attributable_to_single_margin"] += 1
        onset_widths.append(episode.onset_interval_width)
        recovery_widths.append(episode.recovery_interval_width or 0.0)
    return {
        "counts": dict(output),
        "interval_indeterminate_episodes": sum(value for key, value in output.items() if key != "no_observed_recovery_right_censored"),
        "onset_interval_width_hours_median": round(float(np.median(onset_widths)) / 60, 3) if onset_widths else None,
        "recovery_interval_width_hours_median": round(float(np.median(recovery_widths)) / 60, 3) if recovery_widths else None,
        "note": "Attribution asks whether the 48-h boundary can be crossed by onset or recovery uncertainty while holding the other endpoint at its closest observed bound. Categories may be joint by construction.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", choices=("mimic", "sicdb", "eicu"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260830)
    args = parser.parse_args()
    source, spells, labs = load_source(args.database)
    episodes = derive_episodes(source, spells, labs)
    primary = primary_population(episodes)
    inference = _two_stage_hospital_patient_bootstrap(primary, args.bootstrap_replicates, args.seed) if args.database == "eicu" else _cluster_bootstrap(primary, args.bootstrap_replicates, args.seed)
    summary = primary_summary(episodes)
    summary["monitoring_indeterminate"]["cluster_bootstrap_ci95"] = inference["percentile_ci"]["identified_set_width"]
    spell_clusters = [str(spell.extra.get("cluster_id", spell.identifier)) for spell in spells.values()]
    spell_cluster_counts = Counter(spell_clusters)
    episode_cluster_counts = Counter(episode.cluster_id for episode in episodes)
    primary_cluster_counts = Counter(episode.cluster_id for episode in primary)
    output = {
        "source": source,
        "scope": "v4 ICU-coverage primary phenotype, structural-censoring separation, threshold robustness, and cluster-appropriate partial-identification inference.",
        "denominator_definition": "First AKI episodes with database ICU coverage continuing through first positive creatinine plus 48 h; no subsequent creatinine measurement is required for inclusion.",
        "population_structure": {
            "adult_first_continuous_icu_spells": len(spells),
            "unique_patients_in_source_spells": len(spell_cluster_counts),
            "patients_with_multiple_source_spells": sum(count > 1 for count in spell_cluster_counts.values()),
            "first_aki_episodes": len(episodes),
            "unique_patients_with_first_aki_episode": len(episode_cluster_counts),
            "patients_with_multiple_first_aki_episode_admissions": sum(count > 1 for count in episode_cluster_counts.values()),
            "episodes_from_patients_with_multiple_first_aki_episode_admissions": sum(count for count in episode_cluster_counts.values() if count > 1),
            "primary_48h_episodes": len(primary),
            "unique_patients_in_primary_48h_population": len(primary_cluster_counts),
            "patients_with_multiple_primary_48h_episode_admissions": sum(count > 1 for count in primary_cluster_counts.values()),
        },
        "primary": summary,
        "partial_identification_inference": inference,
        "threshold_curve": _threshold_curve_with_coverage(source, spells, labs),
        "monitoring_indeterminacy_decomposition": _indeterminacy_decomposition(primary),
        "interpretation_constraint": "The identified set and bootstrap region describe classification and sampling uncertainty under the stated observation process. They are not estimates of a latent exact biological duration.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / f"{args.database}_v4_primary_inference.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
