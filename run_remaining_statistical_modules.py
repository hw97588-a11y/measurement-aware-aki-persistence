#!/usr/bin/env python3
"""Complete the locked remaining statistical modules for interval-censored AKI.

This program implements STATISTICAL_COMPLETION_AMENDMENT_v1_2_20260830.md.
It reads the external source data without modifying them and writes only
source-level summaries.  It neither pools patient records across databases nor
fits a causal model of monitoring, AKI duration, mortality, or hospital care.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from run_interval_aki_primary import (
    M48,
    M7D,
    Spell,
    as_float,
    deduplicate,
    load_eicu,
    load_mimic,
    load_sicdb,
    wilson,
)


M24 = 24 * 60
M72 = 72 * 60
CATEGORY_ORDER = [
    "definite_transient",
    "definite_persistent",
    "interval_indeterminate",
    "right_censored_unresolved",
]
UNRESOLVED = {"interval_indeterminate", "right_censored_unresolved"}
RECOVERY_POLICIES = ("standard", "two_consecutive", "sustained_48h", "complete_1p2x")


@dataclass
class ExtendedEpisode:
    identifier: str
    category: str
    potential48: bool
    potential24: bool
    onset_lower: float
    onset_upper: float
    baseline: float
    stage2plus: bool
    n_0_48h: int
    no_identifiable_aki_first48h: bool
    current_creatinine: float
    previous_creatinine: float | None
    recovery_upper: float | None
    recovery_value: float | None


def rolling_positive(series: list[tuple[float, float]], index: int) -> tuple[bool, float | None, float | None]:
    """KDIGO creatinine positivity using values before the index measurement."""
    time, creatinine = series[index]
    previous = [(t, c) for t, c in series[:index] if 0 < time - t <= M7D]
    baseline48 = min((c for t, c in previous if time - t <= M48), default=None)
    baseline7 = min((c for _, c in previous), default=None)
    positive = bool(
        (baseline48 is not None and creatinine - baseline48 >= 0.3)
        or (baseline7 is not None and creatinine / baseline7 >= 1.5)
    )
    return positive, baseline48, baseline7


def confirmation_index(
    post: list[tuple[float, float]],
    recovery_limit: float,
    policy: str,
) -> int | None:
    """Return the first accepted recovery value under a locked policy."""
    recovered = [value < recovery_limit for _, value in post]
    if policy == "standard":
        return next((index for index, value in enumerate(recovered) if value), None)
    if policy == "two_consecutive":
        return next(
            (
                index
                for index in range(len(post) - 1)
                if recovered[index] and recovered[index + 1]
            ),
            None,
        )
    if policy == "sustained_48h":
        for index, ((time, _), is_recovered) in enumerate(zip(post, recovered, strict=True)):
            if not is_recovered:
                continue
            has_48h_confirmation = False
            relapse = False
            for later_time, later_value in post[index + 1 :]:
                if later_value >= recovery_limit:
                    relapse = True
                    break
                if later_time >= time + M48:
                    has_48h_confirmation = True
                    break
            if has_48h_confirmation and not relapse:
                return index
        return None
    if policy == "complete_1p2x":
        # The baseline is required for this criterion and is available only in
        # ``first_episode``; that caller resolves this policy before entering
        # the generic confirmation helper.
        raise RuntimeError("complete_1p2x must be handled in first_episode")
    raise ValueError(f"Unknown recovery policy: {policy}")


def first_episode(
    spell: Spell,
    values: list[tuple[float, float]],
    recovery_policy: str = "standard",
) -> ExtendedEpisode | None:
    """Find the first episode and classify its duration with a fixed policy."""
    series = deduplicate(values)
    previous_positive = False
    last_non_aki: float | None = None
    for index, (time, creatinine) in enumerate(series):
        positive, baseline48, baseline7 = rolling_positive(series, index)
        eligible_time = 0 <= time <= min(M7D, spell.end)
        if positive and not previous_positive and eligible_time and last_non_aki is not None and baseline7 is not None:
            recovery_limit = min(baseline7 + 0.3, 1.5 * baseline7)
            post = [(t, c) for t, c in series[index:] if t <= spell.end]
            if recovery_policy == "complete_1p2x":
                complete_limit = min(recovery_limit, 1.2 * baseline7)
                recovery_index = next((j for j, (_, value) in enumerate(post) if value < complete_limit), None)
            else:
                recovery_index = confirmation_index(post, recovery_limit, recovery_policy)

            if recovery_index is not None:
                recovery_upper, recovery_value = post[recovery_index]
                recovery_lower = post[recovery_index - 1][0] if recovery_index > 0 else time
                duration_lower = max(0.0, recovery_lower - time)
                duration_upper = recovery_upper - last_non_aki
                if duration_upper <= M48:
                    category = "definite_transient"
                elif duration_lower > M48:
                    category = "definite_persistent"
                else:
                    category = "interval_indeterminate"
            else:
                recovery_upper = None
                recovery_value = None
                known_positive_times = [t for t, value in post if value >= recovery_limit]
                last_known_positive = max(known_positive_times, default=time)
                duration_lower = max(0.0, last_known_positive - time)
                category = "definite_persistent" if duration_lower > M48 else "right_censored_unresolved"

            stage2plus = bool(
                creatinine / baseline7 >= 2.0
                or (baseline48 is not None and creatinine >= 4.0 and creatinine - baseline48 >= 0.3)
            )
            n_0_48h = sum(0 <= observed <= M48 for observed, _ in series)
            pre_aki = any(rolling_positive(series, j)[0] for j, (observed, _) in enumerate(series) if 0 <= observed < M48)
            previous_creatinine = series[index - 1][1] if index > 0 else None
            return ExtendedEpisode(
                identifier=spell.identifier,
                category=category,
                potential48=spell.end >= time + M48,
                potential24=spell.end >= time + M24,
                onset_lower=last_non_aki,
                onset_upper=time,
                baseline=baseline7,
                stage2plus=stage2plus,
                n_0_48h=n_0_48h,
                no_identifiable_aki_first48h=not pre_aki,
                current_creatinine=creatinine,
                previous_creatinine=previous_creatinine,
                recovery_upper=recovery_upper,
                recovery_value=recovery_value,
            )
        if not positive:
            last_non_aki = time
        previous_positive = positive
    return None


def category_summary(episodes: list[ExtendedEpisode]) -> dict[str, object]:
    categories = Counter(episode.category for episode in episodes)
    denominator = len(episodes)
    not_definite = categories["interval_indeterminate"] + categories["right_censored_unresolved"]
    return {
        "n": denominator,
        "categories": {category: int(categories[category]) for category in CATEGORY_ORDER},
        "not_definitely_classifiable": wilson(int(not_definite), denominator),
    }


def bounds(episodes: list[ExtendedEpisode]) -> dict[str, object]:
    summary = category_summary(episodes)
    n = int(summary["n"])
    cats = summary["categories"]
    if not n:
        return {"n": 0, "persistent_probability_bounds": None, "transient_probability_bounds": None}
    return {
        "n": n,
        "persistent_probability_bounds": {
            "lower_definite_persistent": round(cats["definite_persistent"] / n, 5),
            "upper_if_all_unresolved_are_persistent": round(1 - cats["definite_transient"] / n, 5),
        },
        "transient_probability_bounds": {
            "lower_definite_transient": round(cats["definite_transient"] / n, 5),
            "upper_if_all_unresolved_are_transient": round(1 - cats["definite_persistent"] / n, 5),
        },
        "interpretation": "Sharp descriptive bounds induced by the frozen category definitions; not confidence intervals for a latent biological duration.",
    }


def summarize_recovery_state(episodes: list[ExtendedEpisode]) -> dict[str, object]:
    observed = [episode for episode in episodes if episode.recovery_value is not None]
    complete = sum(episode.recovery_value < 1.2 * episode.baseline for episode in observed)
    partial = sum(1.2 * episode.baseline <= episode.recovery_value < 1.5 * episode.baseline for episode in observed)
    other = len(observed) - complete - partial
    return {
        "episodes_with_observed_primary_recovery": len(observed),
        "complete_recovery_below_1p2x_baseline": complete,
        "partial_recovery_1p2x_to_under_1p5x_baseline": partial,
        "other_primary_recovery_value": other,
        "complete_recovery_proportion_among_observed_recoveries": round(complete / len(observed), 5) if observed else None,
        "partial_recovery_proportion_among_observed_recoveries": round(partial / len(observed), 5) if observed else None,
        "interpretation": "Recovery state describes the first observed primary-recovery measurement and is not an additional duration category.",
    }


def observation_rows(spells: dict[str, Spell], labs: dict[str, list[tuple[float, float]]]) -> tuple[pd.DataFrame, dict[str, tuple[float, float | None, str]]]:
    """Create recurrent 24-h retesting opportunities and onset covariates."""
    rows: list[dict[str, object]] = []
    onset_features: dict[str, tuple[float, float | None, str]] = {}
    for identifier, spell in spells.items():
        series = deduplicate(labs.get(identifier, []))
        for index in range(1, len(series) - 1):
            time, creatinine = series[index]
            if time < 0 or time > M7D or spell.end < time + M24:
                continue
            previous_time, previous = series[index - 1]
            next_time, _ = series[index + 1]
            positive, _, _ = rolling_positive(series, index)
            day = "0" if time < M24 else "1" if time < 2 * M24 else "2" if time < 3 * M24 else "3_7"
            rows.append({
                "spell_id": identifier,
                "last_creatinine": creatinine,
                "delta_creatinine": creatinine - previous,
                "current_aki": "yes" if positive else "no",
                "icu_day": day,
                "next_test_24h": int(next_time - time <= M24),
            })
        episode = first_episode(spell, labs.get(identifier, []), "standard")
        if episode is not None and episode.previous_creatinine is not None:
            day = "0" if episode.onset_upper < M24 else "1" if episode.onset_upper < 2 * M24 else "2" if episode.onset_upper < 3 * M24 else "3_7"
            onset_features[identifier] = (episode.current_creatinine, episode.current_creatinine - episode.previous_creatinine, day)
    return pd.DataFrame(rows), onset_features


def model_observation_process(
    source: str,
    spells: dict[str, Spell],
    labs: dict[str, list[tuple[float, float]]],
    primary: list[ExtendedEpisode],
) -> dict[str, object]:
    data, onset_features = observation_rows(spells, labs)
    if data.empty:
        return {"run": False, "reason": "No recurrent observation opportunities"}
    last_sd = float(data["last_creatinine"].std(ddof=0))
    delta_sd = float(data["delta_creatinine"].std(ddof=0))
    if not math.isfinite(last_sd) or last_sd == 0 or not math.isfinite(delta_sd) or delta_sd == 0:
        return {"run": False, "reason": "No variance in an observation-process predictor"}
    last_mean = float(data["last_creatinine"].mean())
    delta_mean = float(data["delta_creatinine"].mean())
    data["last_creatinine_z"] = (data["last_creatinine"] - last_mean) / last_sd
    data["delta_creatinine_z"] = (data["delta_creatinine"] - delta_mean) / delta_sd
    formula = "next_test_24h ~ last_creatinine_z + delta_creatinine_z + C(current_aki) + C(icu_day)"
    model = smf.gee(
        formula,
        groups="spell_id",
        data=data,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Independence(),
    )
    result = model.fit(maxiter=100)
    confint = result.conf_int()
    effect_rows = []
    label_map = {
        "last_creatinine_z": "Last creatinine (per source SD)",
        "delta_creatinine_z": "Change from preceding creatinine (per source SD)",
        "C(current_aki)[T.yes]": "Current AKI-positive state (vs no)",
        "C(icu_day)[T.1]": "ICU day 1 (vs day 0)",
        "C(icu_day)[T.2]": "ICU day 2 (vs day 0)",
        "C(icu_day)[T.3_7]": "ICU day 3-7 (vs day 0)",
    }
    for term, label in label_map.items():
        if term in result.params.index:
            effect_rows.append({
                "term": label,
                "odds_ratio": round(float(math.exp(result.params[term])), 5),
                "ci95_low": round(float(math.exp(confint.loc[term, 0])), 5),
                "ci95_high": round(float(math.exp(confint.loc[term, 1])), 5),
                "p_value": round(float(result.pvalues[term]), 8),
            })

    records = []
    for episode in primary:
        features = onset_features.get(episode.identifier)
        if features is None:
            continue
        last, delta, day = features
        records.append({
            "last_creatinine_z": (last - last_mean) / last_sd,
            "delta_creatinine_z": (delta - delta_mean) / delta_sd,
            "current_aki": "yes",
            "icu_day": day,
            "not_definitely_classifiable": int(episode.category in UNRESOLVED),
        })
    episode_data = pd.DataFrame(records)
    if episode_data.empty:
        weighting = {"run": False, "reason": "No primary episode had a preceding creatinine for the locked observation model"}
    else:
        probability = np.asarray(result.predict(episode_data), dtype=float)
        stabilizer = float(np.mean(data["next_test_24h"]))
        weights = stabilizer / probability
        naive = float(episode_data["not_definitely_classifiable"].mean())
        weighted = float(np.average(episode_data["not_definitely_classifiable"], weights=weights))
        weighting = {
            "run": True,
            "primary_episodes_with_complete_observation_predictors": int(len(episode_data)),
            "unweighted_not_definitely_classifiable": round(naive, 5),
            "stabilized_inverse_observation_probability_weighted_sensitivity": round(weighted, 5),
            "predicted_retest_probability": {
                "minimum": round(float(np.min(probability)), 5),
                "p01": round(float(np.quantile(probability, .01)), 5),
                "median": round(float(np.median(probability)), 5),
                "p99": round(float(np.quantile(probability, .99)), 5),
                "maximum": round(float(np.max(probability)), 5),
            },
            "stabilized_weight": {
                "minimum": round(float(np.min(weights)), 5),
                "p01": round(float(np.quantile(weights, .01)), 5),
                "median": round(float(np.median(weights)), 5),
                "p99": round(float(np.quantile(weights, .99)), 5),
                "maximum": round(float(np.max(weights)), 5),
            },
            "interpretation": "Sensitivity only. Weighting balances observed first-episode features by the modeled probability of a 24-h retest; it does not recover an unobserved biological recovery time or replace the primary interval-censored estimand.",
        }
    interval = data["next_test_24h"].mean()
    return {
        "run": True,
        "design": "Patient-clustered GEE logistic recurrent observation model. Each eligible observed creatinine is an opportunity; outcome is another measured creatinine within 24 h.",
        "opportunities": int(len(data)),
        "spells": int(data["spell_id"].nunique()),
        "next_test_within_24h_proportion": round(float(interval), 5),
        "predictor_standard_deviations": {"last_creatinine_mg_dL": round(last_sd, 5), "delta_creatinine_mg_dL": round(delta_sd, 5)},
        "effects": effect_rows,
        "inverse_observation_weighted_sensitivity": weighting,
        "interpretation": "Associations demonstrate an informative observation process. They are neither clinical-treatment effects nor causal effects of testing frequency.",
    }


def analyze_source(source: str, loader) -> dict[str, object]:
    spells, labs = loader()
    standard = [episode for identifier, spell in spells.items() if (episode := first_episode(spell, labs.get(identifier, []), "standard")) is not None]
    primary = [episode for episode in standard if episode.potential48]
    sensitivity = {}
    for policy in RECOVERY_POLICIES[1:]:
        policy_episodes = [episode for identifier, spell in spells.items() if (episode := first_episode(spell, labs.get(identifier, []), policy)) is not None]
        sensitivity[policy] = category_summary([episode for episode in policy_episodes if episode.potential48])
    stage2 = [episode for episode in primary if episode.stage2plus]
    strict = [
        episode
        for episode in standard
        if episode.potential48
        and episode.onset_upper >= M48
        and episode.n_0_48h >= 2
        and episode.no_identifiable_aki_first48h
    ]
    no_early_terminal = [episode for episode in standard if episode.potential24]
    result = {
        "source": source,
        "scope": "Remaining pre-specified robustness, partial-identification, and observation-process modules; no causal clinical-effect analysis.",
        "standard_reference": category_summary(primary),
        "partial_identification_bounds": bounds(primary),
        "recovery_definition_sensitivity": sensitivity,
        "recovery_state_description": summarize_recovery_state(standard),
        "cohort_robustness": {
            "stage2plus_first_positive": category_summary(stage2),
            "strict_icu_acquired": category_summary(strict),
            "no_early_terminal_source_end_24h_all_episode_population": category_summary(no_early_terminal),
            "note": "The primary 48-h source-observable population necessarily has at least 24 h of source observability; the final 24-h restriction is therefore reported only for all first episodes.",
        },
        "observation_process": model_observation_process(source, spells, labs, primary),
        "data_processing": {
            "adult_first_continuous_icu_spells": len(spells),
            "first_aki_episodes": len(standard),
            "primary_48h_source_observable_episodes": len(primary),
        },
        "interpretation_constraint": "Databases were analysed separately. No output is a causal effect of observation frequency, phenotype duration, mortality, or hospital care.",
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", choices=("mimic", "sicdb", "eicu"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sources = {
        "mimic": ("MIMIC-IV v3.1", load_mimic),
        "sicdb": ("SICdb v1.0.8", load_sicdb),
        "eicu": ("eICU-CRD v2.0", load_eicu),
    }
    source, loader = sources[args.database]
    output = analyze_source(source, loader)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / f"{args.database}_remaining_statistical_modules.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
