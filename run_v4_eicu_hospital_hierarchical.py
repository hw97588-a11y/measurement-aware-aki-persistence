#!/usr/bin/env python3
"""Convergent eICU hospital-comparability model for the v4 phenotype.

This replaces the exploratory variational-Bayes fit.  It holds the
case-mix coefficients fixed at their patient-level logistic estimates, then
fits a hospital random-intercept variance by a one-dimensional Laplace
marginal likelihood.  This is deliberately a phenotype-comparability
analysis, not a hospital-performance model.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from patsy import dmatrix
from scipy.optimize import minimize_scalar
from scipy.special import expit
from sklearn.linear_model import LogisticRegression

from interval_aki_v4_engine import UNRESOLVED, derive_episodes, load_source, primary_population


def _age_group(value: float | None) -> str:
    return "unknown" if value is None else f"{int(value // 10) * 10}-{int(value // 10) * 10 + 9}"


def _mor(tau: float) -> float:
    return math.exp(math.sqrt(2 * tau * tau) * 0.67448975)


def prepare() -> pd.DataFrame:
    source, spells, labs = load_source("eicu")
    rows = []
    for episode in primary_population(derive_episodes(source, spells, labs)):
        rows.append({
            "hospital": str(episode.hospital_id),
            "age_group": _age_group(episode.age),
            "sex": episode.sex or "unknown",
            "admission_type": episode.admission_type or "unknown",
            "baseline_creatinine": episode.baseline,
            "stage2plus": str(int(episode.stage2plus)),
            "onset_day": str(min(7, int(episode.onset_upper // (24 * 60)))),
            "not_definitely_classifiable": int(episode.category in UNRESOLVED),
        })
    return pd.DataFrame(rows)


FORMULA = (
    "not_definitely_classifiable ~ C(age_group) + C(sex) + C(admission_type) + "
    "baseline_creatinine + C(stage2plus) + C(onset_day)"
)


def _mode_and_log_marginal(y: np.ndarray, eta: np.ndarray, variance: float) -> tuple[float, float, float]:
    """Return hospital posterior mode, conditional SD, and Laplace log likelihood."""
    mode = 0.0
    for _ in range(100):
        probability = expit(np.clip(eta + mode, -35, 35))
        gradient = float((y - probability).sum() - mode / variance)
        information = float((probability * (1 - probability)).sum() + 1 / variance)
        step = gradient / information
        mode += step
        if abs(step) < 1e-10:
            break
    probability = expit(np.clip(eta + mode, -35, 35))
    information = float((probability * (1 - probability)).sum() + 1 / variance)
    log_likelihood = float((y * (eta + mode) - np.logaddexp(0, eta + mode)).sum())
    # Normal-prior and Laplace constants cancel, leaving this variance-aware
    # marginal likelihood up to a common constant.
    log_marginal = log_likelihood - mode * mode / (2 * variance) - 0.5 * math.log(variance) - 0.5 * math.log(information)
    return mode, math.sqrt(1 / information), log_marginal


def _fit_random_intercept(y: np.ndarray, eta: np.ndarray, hospital_codes: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, dict[str, object]]:
    groups = [np.flatnonzero(hospital_codes == code) for code in np.unique(hospital_codes)]

    def objective(log_tau: float) -> float:
        variance = math.exp(2 * float(log_tau))
        return -sum(_mode_and_log_marginal(y[index], eta[index], variance)[2] for index in groups)

    fit = minimize_scalar(objective, bounds=(math.log(.01), math.log(5.0)), method="bounded", options={"xatol": 1e-6})
    tau = math.exp(float(fit.x))
    variance = tau * tau
    modes, sds = [], []
    for index in groups:
        mode, sd, _ = _mode_and_log_marginal(y[index], eta[index], variance)
        modes.append(mode)
        sds.append(sd)
    diagnostics = {
        "optimizer_success": bool(fit.success),
        "optimizer_message": str(fit.message),
        "optimizer_log_tau": round(float(fit.x), 8),
        "optimizer_objective": round(float(fit.fun), 8),
        "tau_at_search_boundary": bool(abs(float(fit.x) - math.log(.01)) < 1e-4 or abs(float(fit.x) - math.log(5.0)) < 1e-4),
    }
    return tau, np.asarray(modes), np.asarray(sds), diagnostics


def _fit_one(data: pd.DataFrame, minimum: int, draws: int, seed: int) -> tuple[dict[str, object], pd.DataFrame]:
    hospital_counts = data.groupby("hospital").size()
    eligible = hospital_counts[hospital_counts >= minimum].index
    subset = data[data["hospital"].isin(eligible)].copy().reset_index(drop=True)
    # A very weak ridge penalty avoids non-identifiability from sparse source
    # categories while keeping this purely a prespecified case-mix nuisance
    # model. It is not used to make a prediction claim.
    design = dmatrix(FORMULA.split("~", 1)[1], subset, return_type="dataframe")
    fixed = LogisticRegression(C=100.0, fit_intercept=False, solver="lbfgs", max_iter=2000)
    fixed.fit(design, subset["not_definitely_classifiable"])
    eta = np.asarray(design) @ fixed.coef_.ravel()
    y = subset["not_definitely_classifiable"].to_numpy(dtype=float)
    category = pd.Categorical(subset["hospital"])
    codes = category.codes
    hospitals = np.asarray(category.categories, dtype=str)
    tau, modes, sds, diagnostic = _fit_random_intercept(y, eta, codes)

    reference_eta = eta
    rows = []
    for code, hospital in enumerate(hospitals):
        part = subset.loc[codes == code]
        mode, sd = float(modes[code]), float(sds[code])
        standardized = float(expit(np.clip(reference_eta + mode, -35, 35)).mean())
        lower = float(expit(np.clip(reference_eta + mode - 1.96 * sd, -35, 35)).mean())
        upper = float(expit(np.clip(reference_eta + mode + 1.96 * sd, -35, 35)).mean())
        info = max(0.0, 1 / (sd * sd) - 1 / (tau * tau))
        reliability = (tau * tau * info) / (1 + tau * tau * info)
        rows.append({
            "hospital": hospital,
            "episodes": int(len(part)),
            "raw_monitoring_indeterminate": float(part["not_definitely_classifiable"].mean()),
            "eb_random_intercept": mode,
            "conditional_posterior_sd": sd,
            "case_mix_standardized_monitoring_indeterminate": standardized,
            "conditional_95_interval_low": lower,
            "conditional_95_interval_high": upper,
            "conditional_reliability": reliability,
        })
    hospital = pd.DataFrame(rows)
    hospital["standardized_rank"] = hospital["case_mix_standardized_monitoring_indeterminate"].rank(method="average")

    rng = np.random.default_rng(seed + minimum)
    draw_effects = rng.normal(modes, sds, size=(draws, len(hospitals)))
    draw_rates = expit(np.clip(reference_eta.mean() + draw_effects, -35, 35))
    ranks = pd.DataFrame(draw_rates).rank(axis=1, method="average").to_numpy()
    hospital["rank_p025"] = np.quantile(ranks, .025, axis=0)
    hospital["rank_p975"] = np.quantile(ranks, .975, axis=0)
    hospital["posterior_probability_highest_uncertainty_quartile"] = (ranks >= .75 * len(hospital)).mean(axis=0)

    rank_width = hospital["rank_p975"] - hospital["rank_p025"]
    summary = {
        "minimum_hospital_episodes": minimum,
        "episodes": int(len(subset)),
        "hospitals": int(len(hospital)),
        "random_intercept_sd_tau": round(tau, 6),
        "random_intercept_variance": round(tau * tau, 6),
        "median_odds_ratio": round(_mor(tau), 6),
        "median_conditional_reliability": round(float(hospital["conditional_reliability"].median()), 6),
        "median_raw_monitoring_indeterminate": round(float(hospital["raw_monitoring_indeterminate"].median()), 6),
        "median_case_mix_standardized_monitoring_indeterminate": round(float(hospital["case_mix_standardized_monitoring_indeterminate"].median()), 6),
        "median_posterior_rank_interval_width": round(float(rank_width.median()), 3),
        "hospitals_with_rank_interval_spanning_at_least_half": int((rank_width >= .5 * len(hospital)).sum()),
        "rank_interval_method": "Conditional normal posterior draws for empirical-Bayes random intercepts; fixed-effect and variance-component uncertainty not propagated, so intervals are diagnostic rather than definitive performance intervals.",
        "fit_diagnostics": diagnostic,
    }
    return summary, hospital


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--posterior-draws", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260830)
    parser.add_argument("--minimum-hospital-episodes", type=int, nargs="+", default=[20, 30, 50])
    parser.add_argument("--input-csv", type=Path, help="Audited cached v4-equivalent episode-level input; avoids rereading protected raw data.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.input_csv is None:
        data = prepare()
        input_provenance = "Directly derived from protected eICU source tables in this run."
    else:
        data = pd.read_csv(args.input_csv, dtype={"hospitalid": str})
        data = data.rename(columns={"hospitalid": "hospital", "onset_icu_day": "onset_day"})
        required = {"hospital", "age_group", "sex", "admission_type", "baseline_creatinine", "stage2plus", "onset_day", "not_definitely_classifiable"}
        missing = required.difference(data.columns)
        if missing:
            raise ValueError(f"Cached input lacks required columns: {sorted(missing)}")
        data["hospital"] = data["hospital"].astype(str)
        data["stage2plus"] = data["stage2plus"].astype(str)
        data["onset_day"] = data["onset_day"].astype(str)
        input_provenance = f"Audited cached eICU input: {args.input_csv}"
    output: dict[str, object] = {
        "source": "eICU-CRD v2.0",
        "scope": "v4 ICU-coverage phenotype comparability; this analysis makes no hospital-quality, treatment-quality, or causal testing-frequency claim.",
        "outcome": "Monitoring-indeterminate 48-hour persistence status among potential-coverage first creatinine-defined AKI episodes.",
        "input_provenance": input_provenance,
        "case_mix": ["age group", "sex", "unit admission source", "fixed episode baseline creatinine", "initial stage 2/3", "AKI onset ICU day"],
        "models": [],
    }
    for minimum in args.minimum_hospital_episodes:
        summary, hospital = _fit_one(data, minimum, args.posterior_draws, args.seed)
        filename = f"eicu_v4_hierarchical_comparability_ridge_laplace_min{minimum}.csv"
        target = args.output_dir / filename
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite {target}")
        hospital.to_csv(target, index=False)
        summary["hospital_file"] = filename
        output["models"].append(summary)
    suffix = "_".join(str(item) for item in args.minimum_hospital_episodes)
    target = args.output_dir / f"eicu_v4_hierarchical_comparability_ridge_laplace_summary_min{suffix}.json"
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
