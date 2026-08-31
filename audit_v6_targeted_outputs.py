#!/usr/bin/env python3
"""Reconcile corrected v6 aggregate outputs before manuscript regeneration."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
DATABASES = ("mimic", "sicdb", "eicu")
CATEGORIES = (
    "definite_transient",
    "definite_persistent",
    "interval_indeterminate",
    "right_censored_unresolved",
)


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def close(a: float, b: float, tolerance: float = 1e-6) -> bool:
    return abs(float(a) - float(b)) <= tolerance


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"name": name, "status": "PASS" if condition else "FAIL", "detail": detail})
        if not condition:
            raise AssertionError(f"{name}: {detail}")

    headline: dict[str, object] = {}
    for database in DATABASES:
        primary = load(RESULTS / "primary" / f"{database}_v4_primary_inference.json")
        sensitivity = load(RESULTS / "robustness" / f"{database}_v5_core_sensitivities.json")
        p = primary["primary"]
        cats = p["primary_categories"]
        denominator = p["primary_48h_icu_coverage_denominator"]
        indeterminate = cats["interval_indeterminate"] + cats["right_censored_unresolved"]
        identified = p["persistent_identified_set"]

        check(f"{database}: categories sum", sum(cats.values()) == denominator, f"{sum(cats.values())} = {denominator}")
        check(
            f"{database}: all episodes reconcile",
            p["all_first_aki_episodes"] == denominator + p["structural_coverage_censored_before_48h"],
            "all index episodes equal primary opportunity plus insufficient opportunity",
        )
        check(
            f"{database}: indeterminate numerator",
            p["monitoring_indeterminate"]["n"] == indeterminate,
            f"{p['monitoring_indeterminate']['n']} = {indeterminate}",
        )
        check(
            f"{database}: identified width",
            close(identified["width_monitoring_indeterminate"], indeterminate / denominator),
            "identified-set width equals classification-indeterminate proportion",
        )
        check(
            f"{database}: lower bound",
            close(identified["lower_definite_persistent"], cats["definite_persistent"] / denominator),
            "lower bound equals definite-persistent proportion",
        )
        check(
            f"{database}: upper bound",
            close(identified["upper_not_definite_transient"], 1 - cats["definite_transient"] / denominator),
            "upper bound equals one minus definite-transient proportion",
        )
        check(
            f"{database}: primary/core agreement",
            sensitivity["recovery_definition_sensitivity"]["first_recovery_primary"]["categories"] == cats,
            "primary categories agree with core-sensitivity source",
        )
        stage1 = sensitivity["initial_aki_stage_sensitivity"]["stage_1"]
        stage23 = sensitivity["initial_aki_stage_sensitivity"]["stage_2_or_3"]
        check(
            f"{database}: stage strata sum",
            stage1["episodes"] + stage23["episodes"] == denominator,
            "Stage 1 plus Stage 2/3 equals the primary denominator",
        )
        one = sensitivity["one_episode_per_unique_patient_sensitivity"]
        check(
            f"{database}: one episode per patient",
            one["episodes"] == one["unique_patients"],
            "deduplicated episode count equals unique-patient count",
        )
        feasibility = sensitivity["all_index_episode_phenotype_feasibility"]
        fsum = sum(
            feasibility[key]["n"]
            for key in (
                "with_48h_opportunity_and_uniquely_classifiable",
                "with_48h_opportunity_but_classification_indeterminate",
                "insufficient_48h_icu_observation_opportunity",
            )
        )
        check(
            f"{database}: feasibility categories sum",
            fsum == feasibility["all_first_observed_index_aki_episodes"],
            "three mutually exclusive feasibility categories reconcile",
        )
        headline[database] = {
            "all_index_episodes": p["all_first_aki_episodes"],
            "primary_48h_opportunity_episodes": denominator,
            "unique_patients_primary": primary["population_structure"]["unique_patients_in_primary_48h_population"],
            "categories": cats,
            "classification_indeterminate": p["monitoring_indeterminate"],
            "persistent_episode_proportion_identified_set": [
                identified["lower_definite_persistent"],
                identified["upper_not_definite_transient"],
            ],
            "feasibility": feasibility,
            "stage_1_indeterminate": stage1["classification_indeterminate"],
            "stage_2_or_3_indeterminate": stage23["classification_indeterminate"],
            "one_episode_per_patient_indeterminate": one["classification_indeterminate"],
        }

    thinning: dict[str, object] = {}
    total_failure = []
    retention = []
    for hours in (12, 24, 36, 48):
        data = load(RESULTS / "thinning" / f"eicu_v4_controlled_thinning_{hours}h.json")
        check(f"thinning {hours}h: reference size", data["reference_cohort"]["episodes"] == 9790, "fixed reference n=9,790")
        reference_counts = data["reference_cohort"]["reference_categories"]
        check(f"thinning {hours}h: reference categories", sum(reference_counts.values()) == 9790, "reference categories sum to fixed denominator")
        effects = data["effect_decomposition"]
        total_failure.append(effects["total_phenotype_failure_fixed_reference_denominator"]["monte_carlo_median"])
        retention.append(effects["index_episode_retention"]["monte_carlo_median"])
        later = sum(
            row["mean_count_per_phase"]
            for row in data["transition_matrix"]
            if row["thinned_category"] == "index_not_retained_later_recurrence_detected"
        )
        check(f"thinning {hours}h: recurrence separation", later > 0, "later recurrence is detected and reported separately")
        thinning[str(hours)] = effects
    check("thinning: total failure gradient", total_failure == sorted(total_failure), "total failure increases from 12 to 48 h")
    check("thinning: index retention gradient", retention == sorted(retention, reverse=True), "index-episode retention decreases from 12 to 48 h")

    output = {
        "status": "PASS",
        "analysis_version": "v6 targeted reanalysis / manuscript package v1.2",
        "implementation_changes": [
            "Post-onset creatinine follow-up retained through the applicable spell end, while index onset remained restricted to ICU day 0-7.",
            "SICdb coverage end re-anchored as (TimeOfStay - ICUOffset)/60 with non-positive intervals excluded.",
            "Controlled thinning anchored to the fixed index AKI episode; later recurrence cannot substitute for a missed index episode.",
        ],
        "headline": headline,
        "controlled_thinning": thinning,
        "checks": checks,
        "checks_passed": sum(item["status"] == "PASS" for item in checks),
        "checks_total": len(checks),
    }
    target = RESULTS / "v6_targeted_reanalysis_audit.json"
    target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"], "checks": f"{output['checks_passed']}/{output['checks_total']}", "headline": headline}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
