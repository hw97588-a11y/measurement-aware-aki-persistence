"""ICU-coverage, interval-censored creatinine-AKI phenotype engine (v4).

The module is intentionally independent of the locked v1--v3 analysis
outputs. It has no write side effects and retains no raw identifiers in its
summaries. It implements the frozen definitions in
PRESUBMISSION_METHODS_AMENDMENT_v4_20260830.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Callable

from run_interval_aki_primary import M48, M7D, Spell, deduplicate, load_eicu, load_mimic, load_sicdb, sicdb_csv


M24 = 24 * 60
M72 = 72 * 60
UNRESOLVED = {"interval_indeterminate", "right_censored_unresolved"}
CATEGORY_ORDER = ["definite_transient", "definite_persistent", "interval_indeterminate", "right_censored_unresolved"]


@dataclass(frozen=True)
class V4Episode:
    source: str
    cluster_id: str
    hospital_id: str | None
    onset_lower: float
    onset_upper: float
    duration_lower: float
    duration_upper: float | None
    category: str
    coverage_48h: bool
    structural_coverage_censored: bool
    baseline: float
    stage2plus: bool
    initial_aki_stage: int
    recovery_lower: float | None
    recovery_upper: float | None
    recovery_value: float | None
    onset_interval_width: float
    recovery_interval_width: float | None
    n_first48h: int
    no_identifiable_aki_first48h: bool
    age: float | None
    sex: str | None
    admission_type: str | None


def _rolling_positive(series: list[tuple[float, float]], index: int) -> tuple[bool, float | None, float | None]:
    time, creatinine = series[index]
    previous = [(t, c) for t, c in series[:index] if 0 < time - t <= M7D]
    baseline48 = min((c for t, c in previous if time - t <= M48), default=None)
    baseline7 = min((c for _, c in previous), default=None)
    positive = bool(
        (baseline48 is not None and creatinine - baseline48 >= 0.3)
        or (baseline7 is not None and creatinine / baseline7 >= 1.5)
    )
    return positive, baseline48, baseline7


def _stage(creatinine: float, baseline48: float | None, baseline7: float) -> int:
    ratio = creatinine / baseline7
    if ratio >= 3 or (creatinine >= 4.0 and baseline48 is not None and creatinine - baseline48 >= 0.3):
        return 3
    if ratio >= 2:
        return 2
    return 1


def _confirmation_index(post: list[tuple[float, float]], recovery_limit: float, policy: str) -> int | None:
    recovered = [value < recovery_limit for _, value in post]
    if policy == "first_recovery":
        return next((j for j, state in enumerate(recovered) if state), None)
    if policy == "two_recoveries_6h":
        for j in range(len(post) - 1):
            if recovered[j] and recovered[j + 1] and post[j + 1][0] - post[j][0] >= 6 * 60:
                return j
        return None
    if policy == "confirmed_24to48h":
        for j, (time, _) in enumerate(post):
            if not recovered[j]:
                continue
            for k in range(j + 1, len(post)):
                later_time, _ = post[k]
                if later_time - time > M48:
                    break
                if not recovered[k]:
                    break
                if later_time - time >= M24:
                    return j
        return None
    raise ValueError(f"Unknown recovery confirmation policy: {policy}")


def classify_first_episode(
    source: str,
    spell: Spell,
    values: list[tuple[float, float]],
    recovery_policy: str = "first_recovery",
    threshold_minutes: float = M48,
) -> V4Episode | None:
    """Classify first ICU AKI episode under ICU coverage, without KRT assumptions."""
    series = [(t, c) for t, c in deduplicate(values) if -M7D <= t <= spell.end]
    previous_positive = False
    last_non_aki: float | None = None
    for index, (time, creatinine) in enumerate(series):
        positive, baseline48, baseline7 = _rolling_positive(series, index)
        if positive and not previous_positive and 0 <= time <= min(M7D, spell.end) and last_non_aki is not None and baseline7 is not None:
            recovery_limit = min(baseline7 + 0.3, 1.5 * baseline7)
            post = series[index:]
            recovery_index = _confirmation_index(post, recovery_limit, recovery_policy)
            if recovery_index is not None:
                recovery_upper, recovery_value = post[recovery_index]
                recovery_lower = post[recovery_index - 1][0] if recovery_index else time
                duration_lower = max(0.0, recovery_lower - time)
                duration_upper = recovery_upper - last_non_aki
                if duration_upper <= threshold_minutes:
                    category = "definite_transient"
                elif duration_lower > threshold_minutes:
                    category = "definite_persistent"
                else:
                    category = "interval_indeterminate"
            else:
                recovery_lower = recovery_upper = recovery_value = duration_upper = None
                last_positive = max((t for t, value in post if value >= recovery_limit), default=time)
                duration_lower = max(0.0, last_positive - time)
                category = "definite_persistent" if duration_lower > threshold_minutes else "right_censored_unresolved"
            coverage = spell.end >= time + threshold_minutes
            stage = _stage(creatinine, baseline48, baseline7)
            pre_aki = any(_rolling_positive(series, j)[0] for j, (observed, _) in enumerate(series) if 0 <= observed < M48)
            cluster = str(spell.extra.get("cluster_id", spell.identifier))
            return V4Episode(
                source=source,
                cluster_id=cluster,
                hospital_id=str(spell.hospital) if spell.hospital is not None else None,
                onset_lower=last_non_aki,
                onset_upper=time,
                duration_lower=duration_lower,
                duration_upper=duration_upper,
                category=category,
                coverage_48h=coverage,
                structural_coverage_censored=not coverage,
                baseline=baseline7,
                stage2plus=stage >= 2,
                initial_aki_stage=stage,
                recovery_lower=recovery_lower,
                recovery_upper=recovery_upper,
                recovery_value=recovery_value,
                onset_interval_width=time - last_non_aki,
                recovery_interval_width=(recovery_upper - recovery_lower) if recovery_upper is not None and recovery_lower is not None else None,
                n_first48h=sum(0 <= observed <= M48 for observed, _ in series),
                no_identifiable_aki_first48h=not pre_aki,
                age=spell.age,
                sex=spell.sex,
                admission_type=spell.admission_type,
            )
        if not positive:
            last_non_aki = time
        previous_positive = positive
    return None


def _mimic_v4() -> tuple[dict[str, Spell], dict[str, list[tuple[float, float]]]]:
    spells, labs = load_mimic()
    icu_spells = {}
    for identifier, spell in spells.items():
        icu_end = float(spell.extra["continuous_icu_end_minutes"])
        if icu_end > 0:
            extra = dict(spell.extra)
            extra["cluster_id"] = str(spell.extra["subject_id"])
            icu_spells[identifier] = replace(spell, end=icu_end, extra=extra)
    return icu_spells, labs


def _sicdb_v4() -> tuple[dict[str, Spell], dict[str, list[tuple[float, float]]]]:
    spells, labs = load_sicdb()
    # The public cases table exposes the de-identified PatientID, allowing
    # true patient-level (rather than CaseID) resampling for SICdb.
    patient_by_case: dict[str, str] = {}
    reader, process = sicdb_csv("cases.csv.gz")
    for row in reader:
        if row["CaseID"] in spells and row.get("PatientID"):
            patient_by_case[row["CaseID"]] = row["PatientID"]
    if process.wait() != 0:
        raise RuntimeError("Cannot read SICdb PatientID mapping")
    revised = {
        identifier: replace(spell, extra={**spell.extra, "cluster_id": patient_by_case.get(identifier, identifier)})
        for identifier, spell in spells.items()
    }
    return revised, labs


def _eicu_v4() -> tuple[dict[str, Spell], dict[str, list[tuple[float, float]]]]:
    # Use the validated creatinine-only streaming reader for the 2.3-GB eICU
    # laboratory file; it applies the same unit and time rules as load_eicu.
    from run_v4_controlled_thinning import fast_eicu_source

    _, spells, labs = fast_eicu_source()
    revised = {
        identifier: replace(
            spell,
            extra={**spell.extra, "cluster_id": str(spell.extra["uniquepid"])},
        )
        for identifier, spell in spells.items()
    }
    return revised, labs


LOADERS: dict[str, tuple[str, Callable[[], tuple[dict[str, Spell], dict[str, list[tuple[float, float]]]]]]] = {
    "mimic": ("MIMIC-IV v3.1", _mimic_v4),
    "sicdb": ("SICdb v1.0.8", _sicdb_v4),
    "eicu": ("eICU-CRD v2.0", _eicu_v4),
}


def load_source(database: str) -> tuple[str, dict[str, Spell], dict[str, list[tuple[float, float]]]]:
    label, loader = LOADERS[database]
    spells, labs = loader()
    return label, spells, labs


def derive_episodes(
    source: str,
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    recovery_policy: str = "first_recovery",
    threshold_minutes: float = M48,
) -> list[V4Episode]:
    return [
        episode
        for identifier, spell in spells.items()
        if (episode := classify_first_episode(source, spell, labs.get(identifier, []), recovery_policy, threshold_minutes)) is not None
    ]


def primary_population(episodes: list[V4Episode]) -> list[V4Episode]:
    return [episode for episode in episodes if episode.coverage_48h]


def category_counts(episodes: list[V4Episode]) -> dict[str, int]:
    return {category: sum(episode.category == category for episode in episodes) for category in CATEGORY_ORDER}


def primary_summary(episodes: list[V4Episode]) -> dict[str, object]:
    primary = primary_population(episodes)
    cats = category_counts(primary)
    uncertain = cats["interval_indeterminate"] + cats["right_censored_unresolved"]
    return {
        "all_first_aki_episodes": len(episodes),
        "structural_coverage_censored_before_48h": sum(episode.structural_coverage_censored for episode in episodes),
        "primary_48h_icu_coverage_denominator": len(primary),
        "primary_categories": cats,
        "monitoring_indeterminate": {
            "n": uncertain,
            "denominator": len(primary),
            "proportion": round(uncertain / len(primary), 6) if primary else None,
        },
        "persistent_identified_set": {
            "lower_definite_persistent": round(cats["definite_persistent"] / len(primary), 6) if primary else None,
            "upper_not_definite_transient": round(1 - cats["definite_transient"] / len(primary), 6) if primary else None,
            "width_monitoring_indeterminate": round(uncertain / len(primary), 6) if primary else None,
        },
    }
