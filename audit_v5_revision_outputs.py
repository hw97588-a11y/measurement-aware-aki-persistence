#!/usr/bin/env python3
"""Fail-fast reconciliation checks for the corrected NDT revision outputs."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PRIMARY_RESULTS = ROOT / "results" / "primary"
ROBUSTNESS_RESULTS = ROOT / "results" / "robustness"
OBSERVATION_RESULTS = ROOT / "results" / "observation"
AUDIT_RESULT = ROOT / "results" / "v5_revision_output_audit.json"
EXPECTED = {
    "mimic": {
        "n": 10504,
        "categories": [3859, 4022, 1473, 1150],
        "ci": [0.241726, 0.258215],
    },
    "sicdb": {
        "n": 2472,
        "categories": [652, 925, 605, 290],
        "ci": [0.342798, 0.380551],
    },
    "eicu": {
        "n": 14599,
        "categories": [4532, 5104, 3334, 1629],
        "ci": [0.326502, 0.353454],
    },
}


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: object = None) -> None:
        checks.append({"name": name, "passed": bool(condition), "detail": detail})

    for database, expected in EXPECTED.items():
        primary = read(PRIMARY_RESULTS / f"{database}_v4_primary_inference.json")
        sensitivity = read(ROBUSTNESS_RESULTS / f"{database}_v5_core_sensitivities.json")
        observation = read(OBSERVATION_RESULTS / f"{database}_v5_observation_process.json")
        p = primary["primary"]
        cats = p["primary_categories"]
        ordered = [cats[key] for key in ("definite_transient", "definite_persistent", "interval_indeterminate", "right_censored_unresolved")]
        check(f"{database}: locked primary denominator", p["primary_48h_icu_coverage_denominator"] == expected["n"], p["primary_48h_icu_coverage_denominator"])
        check(f"{database}: locked primary categories", ordered == expected["categories"], ordered)
        check(f"{database}: category arithmetic", sum(ordered) == expected["n"], sum(ordered))
        check(f"{database}: width equals indeterminate proportion", abs(p["persistent_identified_set"]["width_monitoring_indeterminate"] - p["monitoring_indeterminate"]["proportion"]) < 1e-12)
        check(f"{database}: one cluster-bootstrap interval for width", p["monitoring_indeterminate"]["cluster_bootstrap_ci95"] == expected["ci"], p["monitoring_indeterminate"]["cluster_bootstrap_ci95"])
        check(f"{database}: primary/sensitivity identity", sensitivity["recovery_definition_sensitivity"]["first_recovery_primary"]["categories"] == cats)
        flow = sensitivity["flow"]
        flow_total = (
            flow["excluded_no_valid_central_creatinine"]
            + flow["excluded_only_one_valid_central_creatinine"]
            + flow["excluded_at_least_two_measurements_but_no_observed_qualifying_nonaki_to_aki_transition"]
            + flow["first_observed_transition_defined_aki_episodes"]
        )
        check(f"{database}: source flow reconciles", flow_total == flow["adult_first_continuous_icu_spells"], flow_total)
        check(f"{database}: structural censoring reconciles", flow["first_observed_transition_defined_aki_episodes"] - flow["excluded_structural_icu_coverage_end_before_48h"] == expected["n"])
        for label, result in sensitivity["recovery_definition_sensitivity"].items():
            check(f"{database}: recovery sensitivity {label} categories sum", sum(result["categories"].values()) == result["episodes"])
        strict = sensitivity["strict_icu_acquired_aki"]
        check(f"{database}: strict ICU-acquired categories sum", sum(strict["categories"].values()) == strict["episodes"])
        check(f"{database}: observation process uses patients", observation["unique_patients"] <= observation["continuous_icu_spells"])
        check(f"{database}: observation process contains no IPW", "inverse" not in json.dumps(observation).lower())

    eicu = read(PRIMARY_RESULTS / "eicu_v4_primary_inference.json")
    check("eICU: true unique-patient two-stage bootstrap", eicu["partial_identification_inference"]["resampling"]["level"].startswith("hospital then unique patient"))
    check("eICU: primary unique patient count", eicu["population_structure"]["unique_patients_in_primary_48h_population"] == 13896)
    check("eICU: repeated primary-admission patient count", eicu["population_structure"]["patients_with_multiple_primary_48h_episode_admissions"] == 589)
    thinning_compare = read(ROBUSTNESS_RESULTS / "eicu_v5_core_sensitivities.json")["controlled_thinning_reference_comparison"]
    check("eICU: fixed thinning reference count", thinning_compare["measurement_rich_thinning_reference_cohort"]["episodes"] == 9323)
    check("eICU: thinning reference is selected subset", 0 < thinning_compare["selection_fraction"] < 1, thinning_compare["selection_fraction"])

    failures = [item for item in checks if not item["passed"]]
    output = {"status": "FAIL" if failures else "PASS", "checks_passed": len(checks) - len(failures), "checks_total": len(checks), "checks": checks}
    AUDIT_RESULT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: output[key] for key in ("status", "checks_passed", "checks_total")}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
