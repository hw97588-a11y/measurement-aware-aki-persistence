#!/usr/bin/env python3
"""Patient-clustered ICU-only recurrent creatinine observation model.

This descriptive module demonstrates that repeat testing depends on already
observed history.  It intentionally does not apply inverse-observation
weighting to the persistence estimand.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from interval_aki_v4_engine import M24, M7D, _rolling_positive, load_source
from run_interval_aki_primary import deduplicate


def build_rows(spells, labs) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for identifier, spell in spells.items():
        patient_id = str(spell.extra.get("cluster_id", spell.identifier))
        series = [
            (time, value)
            for time, value in deduplicate(labs.get(identifier, []))
            if -M7D <= time <= min(M7D, spell.end)
        ]
        for index in range(1, len(series) - 1):
            time, creatinine = series[index]
            if time < 0 or spell.end < time + M24:
                continue
            previous = series[index - 1][1]
            next_time = series[index + 1][0]
            positive, _, _ = _rolling_positive(series, index)
            day = "0" if time < M24 else "1" if time < 2 * M24 else "2" if time < 3 * M24 else "3_7"
            rows.append({
                "patient_id": patient_id,
                "spell_id": identifier,
                "last_creatinine": creatinine,
                "delta_creatinine": creatinine - previous,
                "current_aki": "yes" if positive else "no",
                "icu_day": day,
                "next_test_24h": int(next_time - time <= M24),
            })
    return pd.DataFrame(rows)


def analyse(database: str) -> dict[str, object]:
    source, spells, labs = load_source(database)
    data = build_rows(spells, labs)
    if data.empty:
        raise RuntimeError("No eligible recurrent observation opportunities")
    last_mean = float(data["last_creatinine"].mean())
    last_sd = float(data["last_creatinine"].std(ddof=0))
    delta_mean = float(data["delta_creatinine"].mean())
    delta_sd = float(data["delta_creatinine"].std(ddof=0))
    data["last_creatinine_z"] = (data["last_creatinine"] - last_mean) / last_sd
    data["delta_creatinine_z"] = (data["delta_creatinine"] - delta_mean) / delta_sd
    formula = "next_test_24h ~ last_creatinine_z + delta_creatinine_z + C(current_aki) + C(icu_day)"
    model = smf.gee(
        formula,
        groups="patient_id",
        data=data,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Independence(),
    )
    result = model.fit(maxiter=100)
    confint = result.conf_int()
    effects = []
    labels = {
        "last_creatinine_z": "Last creatinine, per source SD",
        "delta_creatinine_z": "Change from preceding creatinine, per source SD",
        "C(current_aki)[T.yes]": "Current creatinine is AKI-positive versus not AKI-positive",
        "C(icu_day)[T.1]": "ICU day 1 versus day 0",
        "C(icu_day)[T.2]": "ICU day 2 versus day 0",
        "C(icu_day)[T.3_7]": "ICU day 3-7 versus day 0",
    }
    for term, label in labels.items():
        if term not in result.params:
            continue
        effects.append({
            "term": label,
            "odds_ratio": round(float(math.exp(result.params[term])), 6),
            "ci95_low": round(float(math.exp(confint.loc[term, 0])), 6),
            "ci95_high": round(float(math.exp(confint.loc[term, 1])), 6),
            "p_value": round(float(result.pvalues[term]), 10),
        })

    yes = data.copy()
    no = data.copy()
    yes["current_aki"] = "yes"
    no["current_aki"] = "no"
    probability_yes = float(np.mean(result.predict(yes)))
    probability_no = float(np.mean(result.predict(no)))
    crude = data.groupby("current_aki", observed=True)["next_test_24h"].agg(["count", "mean"])
    return {
        "source": source,
        "scope": "ICU-only recurrent observation process; patient-clustered GEE; no IPW and no causal monitoring effect.",
        "opportunities": int(len(data)),
        "unique_patients": int(data["patient_id"].nunique()),
        "continuous_icu_spells": int(data["spell_id"].nunique()),
        "crude_next_test_within_24h": {
            state: {"opportunities": int(row["count"]), "proportion": round(float(row["mean"]), 6)}
            for state, row in crude.iterrows()
        },
        "patient_clustered_gee": {
            "formula": formula,
            "working_correlation": "independence",
            "effects": effects,
        },
        "standardized_retest_probability": {
            "if_current_creatinine_aki_positive": round(probability_yes, 6),
            "if_current_creatinine_not_aki_positive": round(probability_no, 6),
            "risk_difference": round(probability_yes - probability_no, 6),
            "risk_ratio": round(probability_yes / probability_no, 6),
            "method": "Marginal standardization over the observed covariate distribution; point estimates are descriptive.",
        },
        "interpretation": "Repeat testing was associated with already observed creatinine history. The model cannot reconstruct unmeasured onset or recovery and was not used to correct the primary persistence bounds.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", choices=("mimic", "sicdb", "eicu"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = analyse(args.database)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / f"{args.database}_v5_observation_process.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
