#!/usr/bin/env python3
"""Fixed-reference controlled thinning of maximally observed eICU trajectories.

The reference is an observed, measurement-rich trajectory rather than a
biological gold standard. Every imposed schedule retains the same reference
AKI episodes as its denominator, separating non-detection from uncertainty
among the AKI episodes that remain detected.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

from interval_aki_v4_engine import M72, UNRESOLVED, classify_first_episode
from run_interval_aki_primary import EICU, M7D, Spell, as_float, deduplicate, is_adult, valid_creatinine


CATEGORIES = ["definite_transient", "definite_persistent", "interval_indeterminate", "right_censored_unresolved"]


def thin(values: list[tuple[float, float]], phase: float, interval: float) -> list[tuple[float, float]]:
    """Select one observed value closest to every phase-shifted bin centre."""
    selected: dict[int, tuple[float, float, float]] = {}
    for time, value in deduplicate(values):
        index = math.floor((time - phase) / interval)
        centre = phase + (index + .5) * interval
        candidate = (abs(time - centre), time, value)
        if index not in selected or candidate < selected[index]:
            selected[index] = candidate
    return [(time, value) for _, time, value in sorted(selected.values(), key=lambda row: row[1])]


def percentile_summary(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "monte_carlo_median": round(float(np.median(array)), 6),
        "monte_carlo_p025": round(float(np.quantile(array, .025)), 6),
        "monte_carlo_p975": round(float(np.quantile(array, .975)), 6),
    }


def quartiles(values: np.ndarray) -> np.ndarray:
    ranks = rankdata(values, method="average")
    return np.minimum(4, np.maximum(1, np.ceil(4 * ranks / len(values)).astype(int)))


def fast_eicu_source() -> tuple[str, dict[str, Spell], dict[str, list[tuple[float, float]]]]:
    """Read only the creatinine rows, streaming the protected gzip in C tools.

    The standard CSV loader parses every one of the roughly 2.3-Gb
    uncompressed laboratory rows in Python.  Here gzip and ripgrep perform
    the broad text filter, after which Python validates the 1.28M candidate
    creatinine records using the same units and time rules as the v4 loader.
    """
    by_health: dict[str, list[tuple[float, float, float, str, str, str, float | None, str, str, float | None, str]]] = defaultdict(list)
    with gzip.open(EICU / "patient.csv.gz", "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle):
            if not is_adult(row["age"]):
                continue
            hospital_offset, unit_end = as_float(row["hospitaladmitoffset"]), as_float(row["unitdischargeoffset"])
            if hospital_offset is None or unit_end is None or unit_end <= 0:
                continue
            start = -hospital_offset
            by_health[row["patienthealthsystemstayid"]].append((
                start, start + unit_end, as_float(row["unitvisitnumber"]) or 0, row["patientunitstayid"], row["hospitalid"],
                row["uniquepid"], as_float(row["age"]), row["gender"], row["unitadmitsource"], as_float(row["hospitaldischargeoffset"]), row["hospitaldischargestatus"],
            ))
    spells: dict[str, Spell] = {}
    unit_mapping: dict[str, tuple[str, float]] = {}
    for health, units in by_health.items():
        units.sort(key=lambda item: (item[0], item[2], item[3]))
        first_start, current_end, _, first_stay, hospital, unique_patient, age, sex, admission_type, hospital_end_offset, hospital_status = units[0]
        included = [first_stay]
        for unit_start, unit_end, _, stay, unit_hospital, *_ in units[1:]:
            if unit_hospital != hospital or unit_start > current_end + 4 * 60:
                break
            included.append(stay)
            current_end = max(current_end, unit_end)
        spells[health] = Spell(
            health,
            current_end - first_start,
            hospital=hospital,
            age=age,
            sex=sex,
            admission_type=admission_type,
            extra={
                "hospital_discharge_offset": hospital_end_offset,
                "hospital_discharge_status": hospital_status,
                "uniquepid": unique_patient,
            },
        )
        for unit_start, _, _, stay, unit_hospital, *_ in units:
            if stay in included and unit_hospital == hospital:
                unit_mapping[stay] = (health, unit_start - first_start)
    raw_units: dict[str, Counter] = defaultdict(Counter)
    labs: dict[str, list[tuple[float, float]]] = defaultdict(list)
    gzip_process = subprocess.Popen(["gzip", "-dc", str(EICU / "lab.csv.gz")], stdout=subprocess.PIPE)
    if gzip_process.stdout is None:
        raise RuntimeError("Cannot start gzip lab stream")
    filter_process = subprocess.Popen(["rg", "-i", ",1,creatinine,"], stdin=gzip_process.stdout, stdout=subprocess.PIPE)
    gzip_process.stdout.close()
    if filter_process.stdout is None:
        raise RuntimeError("Cannot start creatinine row filter")
    columns = ["labid", "patientunitstayid", "labresultoffset", "labtypeid", "labname", "labresult", "labresulttext", "labmeasurenamesystem", "labmeasurenameinterface", "labresultrevisedoffset"]
    for raw in filter_process.stdout:
        values = next(csv.reader([raw.decode("utf-8", errors="replace").rstrip("\r\n")]))
        if len(values) != len(columns):
            continue
        row = dict(zip(columns, values))
        mapping = unit_mapping.get(row["patientunitstayid"])
        if mapping is None:
            continue
        spell_id, unit_start = mapping
        hospital = spells[spell_id].hospital
        unit = row["labmeasurenamesystem"].strip().casefold().replace("μ", "µ")
        raw_units[str(hospital)][unit] += 1
        value = as_float(row["labresult"])
        if unit == "mg/dl":
            pass
        elif unit in {"µmol/l", "umol/l"}:
            value = value / 88.4 if value is not None else None
        else:
            continue
        offset = as_float(row["labresultoffset"])
        spell = spells[spell_id]
        source_offset = unit_start + offset if offset is not None else None
        if source_offset is not None and valid_creatinine(value) and -M7D <= source_offset <= spell.end:
            labs[spell_id].append((source_offset, value))
    if filter_process.wait() != 0 or gzip_process.wait() != 0:
        raise RuntimeError("Creatinine row stream did not finish cleanly")
    incompatible = {hospital for hospital, units in raw_units.items() if any(unit not in {"mg/dl", "µmol/l", "umol/l"} for unit in units)}
    if incompatible:
        spells = {identifier: spell for identifier, spell in spells.items() if str(spell.hospital) not in incompatible}
        labs = {identifier: values for identifier, values in labs.items() if identifier in spells}
    return "eICU-CRD v2.0", spells, labs


def source_reference():
    source, spells, labs = fast_eicu_source()
    reference = {}
    for identifier, spell in spells.items():
        episode = classify_first_episode(source, spell, labs.get(identifier, []))
        n72 = sum(
            episode is not None and episode.onset_upper <= time <= episode.onset_upper + M72
            for time, _ in deduplicate(labs.get(identifier, []))
        )
        if episode is not None and episode.coverage_48h and n72 >= 4:
            reference[identifier] = episode
    return source, spells, labs, reference


def index_episode_match_state(original, scheduled) -> str:
    """Match a thinned episode to the fixed index episode.

    The index episode is retained when the thinned onset interval overlaps the
    original onset interval or its first positive measurement occurs before
    the original observed recovery.  A first thinned AKI beginning only after
    observed recovery is a later recurrence and cannot substitute for the
    missed index episode.  Without observed recovery, the prespecified
    single-episode continuity convention is retained.
    """
    if scheduled is None:
        return "index_not_detected"
    overlap = max(original.onset_lower, scheduled.onset_lower) <= min(original.onset_upper, scheduled.onset_upper)
    before_observed_recovery = original.recovery_upper is None or scheduled.onset_upper < original.recovery_upper
    if overlap or before_observed_recovery:
        return "index_retained"
    return "index_not_retained_later_recurrence_detected"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument("--interval-hours", type=int, nargs="+", default=[12, 24, 36, 48])
    parser.add_argument("--hospital-minimum-reference-episodes", type=int, default=20)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    progress_path = args.output_dir / "controlled_thinning_progress.json"

    print("Loading protected eICU source and fixing the reference cohort...", flush=True)
    source, spells, labs, reference = source_reference()
    if not reference:
        raise RuntimeError("No high-density reference episodes")
    by_hospital: dict[str, list[str]] = defaultdict(list)
    for identifier, episode in reference.items():
        by_hospital[str(episode.hospital_id)].append(identifier)
    eligible_hospitals = sorted(hospital for hospital, identifiers in by_hospital.items() if len(identifiers) >= args.hospital_minimum_reference_episodes)
    hospital_index = {hospital: index for index, hospital in enumerate(eligible_hospitals)}
    hospital_sizes = np.asarray([len(by_hospital[hospital]) for hospital in eligible_hospitals], dtype=float)
    raw_failure = np.asarray([
        sum(reference[identifier].category in UNRESOLVED for identifier in by_hospital[hospital]) / len(by_hospital[hospital])
        for hospital in eligible_hospitals
    ], dtype=float)
    reference_counts = Counter(episode.category for episode in reference.values())
    rng = np.random.default_rng(args.seed)
    progress_path.write_text(json.dumps({"status": "reference_fixed", "episodes": len(reference), "hospitals": len(by_hospital), "rank_hospitals": len(eligible_hospitals)}, indent=2) + "\n")

    for hours in args.interval_hours:
        interval = hours * 60
        target = args.output_dir / f"eicu_v4_controlled_thinning_{hours}h.json"
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite {target}")
        print(f"Running {args.replicates} random phases at {hours} h...", flush=True)
        per_phase: dict[str, list[float]] = defaultdict(list)
        transition = Counter()
        phases = rng.uniform(0, interval, args.replicates)
        for replicate, phase in enumerate(phases, start=1):
            counts = Counter()
            hospital_failure = np.zeros(len(eligible_hospitals), dtype=float)
            for identifier, original in reference.items():
                scheduled = classify_first_episode(source, spells[identifier], thin(labs.get(identifier, []), float(phase), interval))
                ref_category = original.category
                hospital = str(original.hospital_id)
                match_state = index_episode_match_state(original, scheduled)
                if match_state == "index_not_detected":
                    state = "index_not_detected"
                    counts[state] += 1
                    failed = 1
                elif match_state == "index_not_retained_later_recurrence_detected":
                    state = match_state
                    counts[state] += 1
                    counts["later_recurrence_detected"] += 1
                    failed = 1
                elif not scheduled.coverage_48h:
                    state = "index_retained_without_potential_48h_coverage"
                    counts[state] += 1
                    failed = 1
                else:
                    state = scheduled.category
                    counts["retained_primary_eligible_aki"] += 1
                    counts[state] += 1
                    failed = int(state in UNRESOLVED)
                counts["retained_index_aki"] += int(match_state == "index_retained")
                counts["total_failure"] += failed
                transition[(ref_category, state)] += 1
                if hospital in hospital_index:
                    hospital_failure[hospital_index[hospital]] += failed
            total = len(reference)
            retained_primary = counts["retained_primary_eligible_aki"]
            per_phase["index_episode_retention"].append(counts["retained_index_aki"] / total)
            per_phase["primary_eligible_retention"].append(retained_primary / total)
            per_phase["conditional_monitoring_indeterminate_among_retained_primary_eligible"].append(
                (counts["interval_indeterminate"] + counts["right_censored_unresolved"]) / retained_primary if retained_primary else float("nan")
            )
            per_phase["total_phenotype_failure_fixed_reference_denominator"].append(counts["total_failure"] / total)
            standardized_failure = hospital_failure / hospital_sizes
            rho = spearmanr(raw_failure, standardized_failure).statistic
            per_phase["hospital_raw_vs_thinned_failure_rank_spearman"].append(float(rho))
            per_phase["hospital_changed_failure_quartile_proportion"].append(float((quartiles(raw_failure) != quartiles(standardized_failure)).mean()))
            if replicate % 25 == 0 or replicate == args.replicates:
                progress_path.write_text(json.dumps({
                    "status": "running", "interval_hours": hours, "completed_replicates": replicate,
                    "replicates": args.replicates, "reference_episodes": len(reference),
                }, indent=2) + "\n")
        transition_rows = []
        state_order = [
            "index_not_detected",
            "index_not_retained_later_recurrence_detected",
            "index_retained_without_potential_48h_coverage",
            *CATEGORIES,
        ]
        for ref_category in CATEGORIES:
            for state in state_order:
                count = transition[(ref_category, state)]
                if count:
                    transition_rows.append({
                        "reference_category": ref_category,
                        "thinned_category": state,
                        "mean_count_per_phase": round(count / args.replicates, 3),
                        "mean_proportion_within_reference_category": round(count / (args.replicates * reference_counts[ref_category]), 6),
                    })
        output = {
            "source": source,
            "scope": "Post-result controlled thinning robustness analysis. The reference trajectory is maximally observed, not a biological gold standard. Hospital quantities describe phenotype sensitivity, never care quality.",
            "reference_cohort": {
                "definition": "Fixed first creatinine-defined AKI episodes with potential ICU coverage through 48 h and at least four observed creatinines in the first 72 h.",
                "episodes": len(reference), "hospitals": len(by_hospital),
                "hospital_count_for_rank_diagnostic": len(eligible_hospitals),
                "minimum_reference_episodes_per_rank_hospital": args.hospital_minimum_reference_episodes,
                "reference_categories": dict(reference_counts),
            },
            "schedule": {
                "imposed_maximum_sampling_frequency_hours": hours,
                "replicates": args.replicates, "seed": args.seed,
                "random_phase": "One global phase per replicate, sampled uniformly from [0, interval).",
                "selection_rule": "One existing observation closest to each phase-shifted bin centre; no values are imputed.",
            },
            "effect_decomposition": {name: percentile_summary(values) for name, values in per_phase.items()},
            "transition_matrix": transition_rows,
            "interpretation": {
                "index_episode_retention": "Fixed reference index-AKI episodes still identified as the same episode after thinning; later recurrent AKI cannot substitute for a missed index episode.",
                "conditional_monitoring_indeterminate": "Indeterminate proportion only among thinned episodes that remain AKI and retain potential 48-h ICU coverage.",
                "total_phenotype_failure": "Fixed-reference denominator: index episode not detected, only a later recurrence detected, index detected after losing potential 48-h coverage, or index detected but monitoring-indeterminate.",
                "episode_matching": "The thinned onset interval had to overlap the reference onset interval or the thinned first positive had to precede the reference observed recovery; later recurrence was reported separately.",
                "hospital_rank": "Raw-versus-thinned failure-rate ranks are a Monte Carlo sensitivity diagnostic. They are not hospital performance ranks and are not shrinkage-adjusted.",
            },
        }
        target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"Finished {hours} h: {target}", flush=True)
    progress_path.write_text(json.dumps({"status": "complete", "interval_hours": args.interval_hours, "replicates": args.replicates, "reference_episodes": len(reference)}, indent=2) + "\n")


if __name__ == "__main__":
    main()
