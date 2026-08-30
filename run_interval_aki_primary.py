#!/usr/bin/env python3
"""Primary, non-destructive interval-censored AKI phenotype analysis.

This program implements ANALYSIS_LOCK_v1_20260830.md.  It reports only
phenotype distributions and interval widths.  It deliberately does not read
or analyse death, KRT, or other clinical-outcome associations.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import os
import statistics
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path


def _configured_path(environment_variable: str) -> Path | None:
    value = os.environ.get(environment_variable)
    return Path(value).expanduser() if value else None


def require_source_path(path: Path | None, environment_variable: str) -> Path:
    if path is None:
        raise RuntimeError(f"Set {environment_variable} to the governed source-data path before running this analysis.")
    if not path.exists():
        raise FileNotFoundError(f"{environment_variable} does not exist: {path}")
    return path


MIMIC = _configured_path("MIMIC_IV_PATH")
EICU = _configured_path("EICU_CRD_PATH")
SICDB = _configured_path("SICDB_PATH")
SICDB_MEMBER_ROOT = os.environ.get(
    "SICDB_MEMBER_ROOT",
    "salzburg-intensive-care-database-sicdb-a-freely-accessible-intensive-care-database-1.0.8",
)

M48 = 48 * 60
M72 = 72 * 60
M7D = 7 * 24 * 60


def as_float(value: object) -> float | None:
    try:
        output = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return output if math.isfinite(output) else None


def as_datetime(value: str) -> datetime | None:
    try:
        return datetime.fromisoformat(value)
    except (TypeError, ValueError):
        return None


def is_adult(value: object) -> bool:
    text = str(value).strip()
    return text.startswith(">") or ((age := as_float(text)) is not None and age >= 18)


def valid_creatinine(value: float | None) -> bool:
    return value is not None and 0.1 <= value <= 25.0


def percentile(values: list[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    return ordered[max(0, math.ceil(probability * len(ordered)) - 1)]


def rounded(value: float | None, digits: int = 2) -> float | None:
    return round(value, digits) if value is not None else None


def quantile_summary(values: list[float]) -> dict[str, object]:
    return {
        "n": len(values),
        "median_hours": rounded(statistics.median(values) / 60) if values else None,
        "p25_hours": rounded(percentile(values, .25) / 60),
        "p75_hours": rounded(percentile(values, .75) / 60),
    }


def deduplicate(values: list[tuple[float, float]]) -> list[tuple[float, float]]:
    by_time: dict[float, list[float]] = defaultdict(list)
    for offset, value in values:
        by_time[offset].append(value)
    return [(offset, statistics.median(results)) for offset, results in sorted(by_time.items())]


def wilson(successes: int, denominator: int) -> dict[str, float | int | None]:
    if denominator == 0:
        return {"n": 0, "denominator": 0, "proportion": None, "ci95_low": None, "ci95_high": None}
    z = 1.959963984540054
    proportion = successes / denominator
    base = 1 + z * z / denominator
    centre = (proportion + z * z / (2 * denominator)) / base
    radius = z * math.sqrt(proportion * (1 - proportion) / denominator + z * z / (4 * denominator * denominator)) / base
    return {
        "n": successes,
        "denominator": denominator,
        "proportion": round(proportion, 5),
        "ci95_low": round(max(0.0, centre - radius), 5),
        "ci95_high": round(min(1.0, centre + radius), 5),
    }


@dataclass
class Spell:
    identifier: str
    end: float
    hospital: str | None = None
    age: float | None = None
    sex: str | None = None
    admission_type: str | None = None
    extra: dict[str, object] = field(default_factory=dict)


@dataclass
class Episode:
    source: str
    hospital: str | None
    onset_lower: float
    onset_upper: float
    recovery_lower: float | None
    recovery_upper: float | None
    duration_lower: float | None
    duration_upper: float | None
    category: str
    potential48: bool
    baseline: float
    stage2plus: bool
    onset_day: int
    age: float | None
    sex: str | None
    admission_type: str | None
    n_72h: int


def first_interval_episode(
    source: str,
    spell: Spell,
    values: list[tuple[float, float]],
    recovery_horizon_after_onset: float | None = None,
) -> Episode | None:
    """Return the first qualifying AKI episode, retaining interval bounds only."""
    series = deduplicate(values)
    previous_positive = False
    last_non_aki: float | None = None
    for index, (time, creatinine) in enumerate(series):
        preceding = [(t, c) for t, c in series[:index] if 0 < time - t <= M7D]
        baseline48 = min((c for t, c in preceding if time - t <= M48), default=None)
        baseline7 = min((c for _, c in preceding), default=None)
        positive = bool(
            (baseline48 is not None and creatinine - baseline48 >= .3)
            or (baseline7 is not None and creatinine / baseline7 >= 1.5)
        )
        eligible_time = 0 <= time <= min(M7D, spell.end)
        if positive and not previous_positive and eligible_time and last_non_aki is not None and baseline7 is not None:
            recovery_limit = min(baseline7 + .3, 1.5 * baseline7)
            ascertainment_end = min(
                spell.end,
                time + recovery_horizon_after_onset if recovery_horizon_after_onset is not None else spell.end,
            )
            post_onset = [(t, c) for t, c in series[index:] if t >= time and t <= ascertainment_end]
            recovery_index = next((j for j, (_, value) in enumerate(post_onset) if value < recovery_limit), None)
            stage2plus = bool(
                creatinine / baseline7 >= 2.0
                or (baseline48 is not None and creatinine >= 4.0 and creatinine - baseline48 >= .3)
            )
            following72 = [(t, c) for t, c in post_onset if t <= time + M72]
            if recovery_index is not None:
                recovery_upper, _ = post_onset[recovery_index]
                recovery_lower, _ = post_onset[recovery_index - 1] if recovery_index > 0 else (time, creatinine)
                duration_lower = max(0.0, recovery_lower - time)
                duration_upper = recovery_upper - last_non_aki
                if duration_upper <= M48:
                    category = "definite_transient"
                elif duration_lower > M48:
                    category = "definite_persistent"
                else:
                    category = "interval_indeterminate"
            else:
                recovery_lower = recovery_upper = duration_upper = None
                last_unrecovered_time = post_onset[-1][0] if post_onset else time
                duration_lower = max(0.0, last_unrecovered_time - time)
                category = "definite_persistent" if duration_lower > M48 else "right_censored_unresolved"
            return Episode(
                source=source,
                hospital=spell.hospital,
                onset_lower=last_non_aki,
                onset_upper=time,
                recovery_lower=recovery_lower,
                recovery_upper=recovery_upper,
                duration_lower=duration_lower,
                duration_upper=duration_upper,
                category=category,
                potential48=spell.end >= time + M48,
                baseline=baseline7,
                stage2plus=stage2plus,
                onset_day=max(0, int(time // (24 * 60))),
                age=spell.age,
                sex=spell.sex,
                admission_type=spell.admission_type,
                n_72h=len(following72),
            )
        if not positive:
            last_non_aki = time
        previous_positive = positive
    return None


def summarize(source: str, spells: dict[str, Spell], labs: dict[str, list[tuple[float, float]]]) -> tuple[dict[str, object], list[Episode]]:
    episodes: list[Episode] = []
    for identifier, spell in spells.items():
        episode = first_interval_episode(source, spell, labs.get(identifier, []))
        if episode is not None:
            episodes.append(episode)
    all_categories = Counter(episode.category for episode in episodes)
    eligible = [episode for episode in episodes if episode.potential48]
    eligible_categories = Counter(episode.category for episode in eligible)
    unresolved = eligible_categories["interval_indeterminate"] + eligible_categories["right_censored_unresolved"]
    onset_widths = [episode.onset_upper - episode.onset_lower for episode in episodes]
    recovery_widths = [
        episode.recovery_upper - episode.recovery_lower
        for episode in episodes
        if episode.recovery_upper is not None and episode.recovery_lower is not None
    ]
    lower_widths = [episode.duration_lower for episode in episodes if episode.duration_lower is not None]
    upper_widths = [episode.duration_upper for episode in episodes if episode.duration_upper is not None]
    summary = {
        "source": source,
        "frozen_analysis_scope": {
            "no_outcomes_or_causal_models": True,
            "category_order": ["definite_transient", "definite_persistent", "interval_indeterminate", "right_censored_unresolved"],
            "primary_estimand": "among episodes with at least 48 h source-observable opportunity after first AKI-positive creatinine, proportion not definitely classifiable at 48 h = interval_indeterminate + right_censored_unresolved",
        },
        "cohort": {
            "spells": len(spells),
            "first_aki_episodes": len(episodes),
            "episodes_with_48h_potential_observation": len(eligible),
        },
        "all_episodes_categories": dict(all_categories),
        "primary_population_categories": dict(eligible_categories),
        "primary_not_definitely_classifiable": wilson(unresolved, len(eligible)),
        "primary_interval_indeterminate": wilson(eligible_categories["interval_indeterminate"], len(eligible)),
        "primary_right_censored_unresolved": wilson(eligible_categories["right_censored_unresolved"], len(eligible)),
        "interval_widths": {
            "onset_interval_width": quantile_summary(onset_widths),
            "recovery_interval_width": quantile_summary(recovery_widths),
            "duration_lower_bound": quantile_summary(lower_widths),
            "duration_upper_bound_observed_recovery_only": quantile_summary(upper_widths),
        },
        "high_density_episode_counts": {
            "at_least_4_creatinines_0_to_72h": sum(episode.n_72h >= 4 for episode in episodes),
            "at_least_6_creatinines_0_to_72h": sum(episode.n_72h >= 6 for episode in episodes),
        },
    }
    return summary, episodes


def mimic_csv(archive: zipfile.ZipFile, member: str):
    return csv.DictReader(io.TextIOWrapper(gzip.GzipFile(fileobj=archive.open(member)), encoding="utf-8", errors="replace", newline=""))


def load_mimic() -> tuple[dict[str, Spell], dict[str, list[tuple[float, float]]]]:
    with zipfile.ZipFile(require_source_path(MIMIC, "MIMIC_IV_PATH")) as archive:
        adults: dict[str, tuple[float | None, str]] = {}
        for row in mimic_csv(archive, "mimic-iv-3.1/hosp/patients.csv.gz"):
            if is_adult(row["anchor_age"]):
                adults[row["subject_id"]] = (as_float(row["anchor_age"]), row["gender"])
        admissions: dict[str, tuple[datetime, datetime, datetime | None, str]] = {}
        for row in mimic_csv(archive, "mimic-iv-3.1/hosp/admissions.csv.gz"):
            admit, discharge = as_datetime(row["admittime"]), as_datetime(row["dischtime"])
            death = as_datetime(row["deathtime"]) if row["deathtime"] else None
            if admit is not None and discharge is not None and discharge > admit:
                admissions[row["hadm_id"]] = (admit, discharge, death, row["admission_type"])
        units: dict[str, list[tuple[datetime, datetime, str]]] = defaultdict(list)
        for row in mimic_csv(archive, "mimic-iv-3.1/icu/icustays.csv.gz"):
            start, end = as_datetime(row["intime"]), as_datetime(row["outtime"])
            if row["subject_id"] in adults and start is not None and end is not None and row["hadm_id"] in admissions:
                units[row["hadm_id"]].append((start, end, row["subject_id"]))
        spells: dict[str, Spell] = {}
        starts: dict[str, datetime] = {}
        for hadm, records in units.items():
            records.sort()
            start, end, subject = records[0]
            for next_start, next_end, _ in records[1:]:
                if next_start <= end + timedelta(hours=4):
                    end = max(end, next_end)
                else:
                    break
            _, discharge, death, admission_type = admissions[hadm]
            source_end = death if death is not None and start < death < discharge else discharge
            if source_end <= start:
                continue
            age, sex = adults[subject]
            spells[hadm] = Spell(
                hadm,
                (source_end - start).total_seconds() / 60,
                age=age,
                sex=sex,
                admission_type=admission_type,
                extra={
                    "subject_id": subject,
                    "start_dt": start,
                    "admit_dt": admissions[hadm][0],
                    "continuous_icu_end_minutes": (end - start).total_seconds() / 60,
                    "hospital_death_dt": death,
                    "hospital_end_dt": source_end,
                },
            )
            starts[hadm] = start
        labs: dict[str, list[tuple[float, float]]] = defaultdict(list)
        raw = gzip.GzipFile(fileobj=archive.open("mimic-iv-3.1/hosp/labevents.csv.gz"))
        next(raw, b"")
        for line in raw:
            if b",50912," not in line and b",52546," not in line:
                continue
            row = next(csv.reader([line.decode("utf-8", errors="replace")]))
            if len(row) < 11:
                continue
            hadm = row[2]
            spell = spells.get(hadm)
            observed, value = as_datetime(row[6]), as_float(row[9])
            if spell is None or observed is None or row[10].strip().casefold() != "mg/dl" or not valid_creatinine(value):
                continue
            offset = (observed - starts[hadm]).total_seconds() / 60
            if -M7D <= offset <= min(M7D, spell.end):
                labs[hadm].append((offset, value))
        raw.close()
    return spells, labs


def eicu_csv(name: str):
    source = require_source_path(EICU, "EICU_CRD_PATH")
    return csv.DictReader(gzip.open(source / name, "rt", encoding="utf-8", errors="replace", newline=""))


def load_eicu() -> tuple[dict[str, Spell], dict[str, list[tuple[float, float]]]]:
    by_health: dict[str, list[tuple[float, float, float, str, str, str, float | None, str, str, float | None, str]]] = defaultdict(list)
    for row in eicu_csv("patient.csv.gz"):
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
        first_start, current_end, _, first_stay, hospital, unique_patient, age, sex, admit_source, hospital_end_offset, hospital_status = units[0]
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
            admission_type=admit_source,
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
    for row in eicu_csv("lab.csv.gz"):
        mapping = unit_mapping.get(row["patientunitstayid"])
        if mapping is None or row["labname"].strip().casefold() != "creatinine" or row["labtypeid"].strip() != "1":
            continue
        spell_id, unit_start = mapping
        hospital = spells[spell_id].hospital
        unit = row["labmeasurenamesystem"].strip().casefold().replace("μ", "µ")
        raw_units[hospital][unit] += 1
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
        if source_offset is not None and valid_creatinine(value) and -M7D <= source_offset <= min(M7D, spell.end):
            labs[spell_id].append((source_offset, value))
    incompatible = {hospital for hospital, unit_counts in raw_units.items() if any(unit not in {"mg/dl", "µmol/l", "umol/l"} for unit in unit_counts)}
    if incompatible:
        spells = {identifier: spell for identifier, spell in spells.items() if spell.hospital not in incompatible}
        labs = {identifier: entries for identifier, entries in labs.items() if identifier in spells}
    return spells, labs


def sicdb_csv(member: str):
    import subprocess
    source = require_source_path(SICDB, "SICDB_PATH")
    process = subprocess.Popen(["bsdtar", "-xOf", str(source), f"{SICDB_MEMBER_ROOT}/{member}"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if process.stdout is None:
        raise RuntimeError(f"Cannot open SICdb member {member}")
    return csv.DictReader(io.TextIOWrapper(gzip.GzipFile(fileobj=process.stdout), encoding="utf-8", errors="replace", newline="")), process


def load_sicdb() -> tuple[dict[str, Spell], dict[str, list[tuple[float, float]]]]:
    spells: dict[str, Spell] = {}
    reader, process = sicdb_csv("cases.csv.gz")
    for row in reader:
        if not is_adult(row["AgeOnAdmission"]):
            continue
        duration, offset = as_float(row["TimeOfStay"]), as_float(row["ICUOffset"])
        if duration is None or offset is None or duration <= 0:
            continue
        spells[row["CaseID"]] = Spell(row["CaseID"], duration / 60, age=as_float(row["AgeOnAdmission"]), admission_type=row["SurgicalAdmissionType"], extra={"icu_offset": offset})
    if process.wait() != 0:
        raise RuntimeError("Cannot read SICdb cases")
    labs: dict[str, list[tuple[float, float]]] = defaultdict(list)
    reader, process = sicdb_csv("laboratory.csv.gz")
    for row in reader:
        spell = spells.get(row["CaseID"])
        if spell is None or row["LaboratoryID"] not in {"367", "368"}:
            continue
        offset, value = as_float(row["Offset"]), as_float(row["LaboratoryValue"])
        if offset is None or not valid_creatinine(value):
            continue
        relative = (offset - float(spell.extra["icu_offset"])) / 60
        if -M7D <= relative <= min(M7D, spell.end):
            labs[spell.identifier].append((relative, value))
    if process.wait() != 0:
        raise RuntimeError("Cannot read SICdb laboratory data")
    return spells, labs


def save_json(path: Path, result: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite {path}")
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def save_eicu_model_input(path: Path, episodes: list[Episode]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite {path}")
    rows = [episode for episode in episodes if episode.potential48]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["hospitalid", "age_group", "sex", "admission_type", "baseline_creatinine", "stage2plus", "onset_icu_day", "not_definitely_classifiable"])
        writer.writeheader()
        for episode in rows:
            age_group = "unknown" if episode.age is None else f"{int(episode.age // 10) * 10}-{int(episode.age // 10) * 10 + 9}"
            writer.writerow({
                "hospitalid": episode.hospital,
                "age_group": age_group,
                "sex": episode.sex or "unknown",
                "admission_type": episode.admission_type or "unknown",
                "baseline_creatinine": round(episode.baseline, 4),
                "stage2plus": int(episode.stage2plus),
                "onset_icu_day": episode.onset_day,
                "not_definitely_classifiable": int(episode.category in {"interval_indeterminate", "right_censored_unresolved"}),
            })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", required=True, choices=("mimic", "sicdb", "eicu"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    loaders = {"mimic": ("MIMIC-IV v3.1", load_mimic), "sicdb": ("SICdb v1.0.8", load_sicdb), "eicu": ("eICU-CRD v2.0", load_eicu)}
    source, loader = loaders[args.database]
    spells, labs = loader()
    summary, episodes = summarize(source, spells, labs)
    summary["source_observability"] = {
        "MIMIC-IV v3.1": "hospital discharge or earlier in-hospital death; hospital labevents retained after ICU transfer",
        "SICdb v1.0.8": "source CaseID TimeOfStay relative to ICUOffset",
        "eICU-CRD v2.0": "end of first continuous same-hospital unit spell; laboratory table is keyed to unit stay",
    }[source]
    save_json(args.output_dir / f"{args.database}_interval_aki_primary_summary.json", summary)
    if args.database == "eicu":
        save_eicu_model_input(args.output_dir / "eicu_interval_aki_hospital_model_input.csv", episodes)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
