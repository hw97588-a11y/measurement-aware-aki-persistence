#!/usr/bin/env python3
"""72-hour landmark construct validation for interval-censored AKI phenotypes.

Implements CONSTRUCT_VALIDATION_AMENDMENT_v1_1_20260830.md.  Categories are
ascertained only using values available by each episode's 72-hour landmark,
preventing post-landmark phenotype leakage.  Results are associative only.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from run_interval_aki_primary import M72, first_interval_episode, load_eicu, load_mimic


REFERENCE = "definite_transient"
CATEGORY_ORDER = ["definite_transient", "definite_persistent", "interval_indeterminate", "right_censored_unresolved"]


def age_group(value: float | None) -> str:
    return "unknown" if value is None else f"{int(value // 10) * 10}-{int(value // 10) * 10 + 9}"


def outcome_offset_mimic(spell) -> tuple[float | None, bool]:
    start = spell.extra["start_dt"]
    end = (spell.extra["hospital_end_dt"] - start).total_seconds() / 60
    death = spell.extra["hospital_death_dt"]
    death_offset = (death - start).total_seconds() / 60 if death is not None and death <= spell.extra["hospital_end_dt"] else None
    return end, death_offset is not None


def make_mimic_records() -> tuple[pd.DataFrame, Counter]:
    spells, labs = load_mimic()
    flow = Counter()
    rows = []
    for hadm, spell in spells.items():
        episode = first_interval_episode("MIMIC-IV v3.1", spell, labs.get(hadm, []), recovery_horizon_after_onset=M72)
        if episode is None:
            continue
        flow["first_aki_episode"] += 1
        landmark = episode.onset_upper + M72
        hospital_end, hospital_death = outcome_offset_mimic(spell)
        if hospital_end is None or hospital_end <= landmark:
            flow["not_in_hospital_at_landmark"] += 1
            continue
        death_dt = spell.extra["hospital_death_dt"]
        death_offset = (death_dt - spell.extra["start_dt"]).total_seconds() / 60 if death_dt is not None else None
        if death_offset is not None and death_offset <= landmark:
            flow["died_by_landmark"] += 1
            continue
        flow["landmark_eligible"] += 1
        rows.append({
            "category": episode.category,
            "death_after_landmark": int(hospital_death and death_offset is not None and death_offset > landmark),
            "age_group": age_group(episode.age),
            "sex": episode.sex or "unknown",
            "admission_type": episode.admission_type or "unknown",
            "baseline_creatinine": episode.baseline,
            "stage2plus": str(int(episode.stage2plus)),
            "onset_icu_day": str(episode.onset_day),
            "hospital": "single_center",
        })
    return pd.DataFrame(rows), flow


def make_eicu_records() -> tuple[pd.DataFrame, Counter]:
    spells, labs = load_eicu()
    flow = Counter()
    rows = []
    for spell_id, spell in spells.items():
        episode = first_interval_episode("eICU-CRD v2.0", spell, labs.get(spell_id, []), recovery_horizon_after_onset=M72)
        if episode is None:
            continue
        flow["first_aki_episode"] += 1
        landmark = episode.onset_upper + M72
        hospital_end = spell.extra.get("hospital_discharge_offset")
        status = str(spell.extra.get("hospital_discharge_status") or "").casefold()
        if not isinstance(hospital_end, (float, int)) or not math.isfinite(hospital_end) or hospital_end <= landmark:
            flow["not_in_hospital_at_landmark"] += 1
            continue
        death = status == "expired"
        if death and hospital_end <= landmark:
            flow["died_by_landmark"] += 1
            continue
        flow["landmark_eligible"] += 1
        rows.append({
            "category": episode.category,
            "death_after_landmark": int(death and hospital_end > landmark),
            "age_group": age_group(episode.age),
            "sex": episode.sex or "unknown",
            "admission_type": episode.admission_type or "unknown",
            "baseline_creatinine": episode.baseline,
            "stage2plus": str(int(episode.stage2plus)),
            "onset_icu_day": str(episode.onset_day),
            "hospital": str(episode.hospital),
        })
    return pd.DataFrame(rows), flow


def risk_table(data: pd.DataFrame) -> list[dict[str, object]]:
    result = []
    for category in CATEGORY_ORDER:
        subset = data[data["category"] == category]
        n, deaths = len(subset), int(subset["death_after_landmark"].sum())
        result.append({
            "category": category,
            "n": n,
            "deaths_after_landmark": deaths,
            "risk": round(deaths / n, 5) if n else None,
        })
    return result


def adjusted_risk_ratios(data: pd.DataFrame, clustered: bool) -> dict[str, object]:
    table = risk_table(data)
    deaths = {row["category"]: row["deaths_after_landmark"] for row in table}
    if any(deaths[category] < 50 for category in CATEGORY_ORDER):
        return {"run": False, "reason": "at least one phenotype category has fewer than 50 post-landmark deaths", "risk_table": table}
    formula = "death_after_landmark ~ C(category, Treatment(reference='definite_transient')) + C(age_group) + C(sex) + C(admission_type) + baseline_creatinine + C(stage2plus) + C(onset_icu_day)"
    model = smf.glm(formula, data=data, family=sm.families.Poisson())
    if clustered and data["hospital"].nunique() > 1:
        fitted = model.fit(cov_type="cluster", cov_kwds={"groups": data["hospital"]})
        covariance = "hospital-cluster robust"
    else:
        fitted = model.fit(cov_type="HC0")
        covariance = "HC0 robust"
    ci = fitted.conf_int()
    estimates = []
    for category in CATEGORY_ORDER[1:]:
        term = f"C(category, Treatment(reference='definite_transient'))[T.{category}]"
        estimates.append({
            "category_vs_definite_transient": category,
            "adjusted_risk_ratio": round(float(np.exp(fitted.params[term])), 5),
            "ci95_low": round(float(np.exp(ci.loc[term, 0])), 5),
            "ci95_high": round(float(np.exp(ci.loc[term, 1])), 5),
        })
    return {
        "run": True,
        "risk_table": table,
        "adjusted_model": "modified Poisson regression; associative construct validation only",
        "covariance": covariance,
        "n": int(len(data)),
        "post_landmark_deaths": int(data["death_after_landmark"].sum()),
        "estimates": estimates,
    }


def analyse(source: str) -> dict[str, object]:
    if source == "mimic":
        data, flow = make_mimic_records()
        result = adjusted_risk_ratios(data, clustered=False)
        source_label = "MIMIC-IV v3.1"
    else:
        data, flow = make_eicu_records()
        result = adjusted_risk_ratios(data, clustered=True)
        source_label = "eICU-CRD v2.0"
    return {
        "source": source_label,
        "design": "72-h after first AKI-positive creatinine landmark; phenotype labels use no measurements later than landmark; eligible only if alive and still in hospital at landmark",
        "flow": dict(flow),
        "model_result": result,
        "interpretation_constraint": "Association-based construct validation; not a causal effect of AKI duration or an evaluation of hospital care.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", choices=("mimic", "eicu"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = analyse(args.database)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / f"{args.database}_landmark_construct_validation.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
