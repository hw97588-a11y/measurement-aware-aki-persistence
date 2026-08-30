#!/usr/bin/env python3
"""Convert private controlled-thinning simulator files into public aggregates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REFERENCE = {0: "definite_transient", 1: "definite_persistent", 2: "interval_indeterminate", 3: "right_censored_unresolved"}
THINNED = {0: "not_detected", 1: "detected_without_potential_48h_coverage", 2: "definite_transient", 3: "definite_persistent", 4: "interval_indeterminate", 5: "right_censored_unresolved"}
METRICS = {
    "retained": "phenotype_retention",
    "primary_retained": "primary_eligible_retention",
    "uncertain": "conditional_monitoring_indeterminate_among_retained_primary_eligible",
    "failure": "total_phenotype_failure_fixed_reference_denominator",
    "rho": "hospital_raw_vs_thinned_failure_rank_spearman",
    "quartile_change": "hospital_changed_failure_quartile_proportion",
}


def summarize(values: pd.Series) -> dict[str, float]:
    return {
        "monte_carlo_median": round(float(values.median()), 6),
        "monte_carlo_p025": round(float(values.quantile(.025)), 6),
        "monte_carlo_p975": round(float(values.quantile(.975)), 6),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulation-dir", type=Path, required=True)
    parser.add_argument("--cache-metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    args = parser.parse_args()
    metadata = json.loads(args.cache_metadata.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_index = []
    for hours in (12, 24, 36, 48):
        metric = pd.read_csv(args.simulation_dir / f"metrics_{hours}.tsv", sep="\t")
        if len(metric) != args.replicates:
            raise ValueError(f"Expected {args.replicates} phases at {hours} h, found {len(metric)}")
        transition = pd.read_csv(args.simulation_dir / f"transitions_{hours}.tsv", sep="\t")
        rows = []
        for row in transition.itertuples(index=False):
            reference = REFERENCE[int(row.reference_category)]
            denominator = metadata["reference_categories"][reference] * args.replicates
            rows.append({
                "reference_category": reference,
                "thinned_category": THINNED[int(row.thinned_state)],
                "mean_count_per_phase": round(int(row.count) / args.replicates, 3),
                "mean_proportion_within_reference_category": round(int(row.count) / denominator, 6),
            })
        output = {
            "source": metadata["source"],
            "scope": "Post-result controlled-thinning robustness analysis. The maximally observed trajectory is a measurement-rich reference, not a biological gold standard. No hospital-quality or causal testing-frequency inference is made.",
            "reference_cohort": {
                "definition": "Fixed first creatinine-defined AKI episodes with potential ICU coverage through 48 h and at least four observed creatinines from first AKI positivity through 72 h after onset.",
                "episodes": metadata["episodes"], "hospitals": metadata["hospitals"],
                "hospital_count_for_rank_diagnostic": int(pd.read_csv(args.simulation_dir / "simulation_meta.tsv", sep="\t", header=None, index_col=0).loc["rank_hospitals", 1]),
                "minimum_reference_episodes_per_rank_hospital": 20,
                "reference_categories": metadata["reference_categories"],
            },
            "schedule": {
                "imposed_maximum_sampling_frequency_hours": hours,
                "replicates": args.replicates,
                "seed": 20260830,
                "random_phase": "One global phase per replicate, sampled uniformly from [0, interval).",
                "selection_rule": "One existing observation closest to each phase-shifted bin centre; no values are imputed.",
            },
            "effect_decomposition": {label: summarize(metric[column]) for column, label in METRICS.items()},
            "transition_matrix": rows,
            "interpretation": {
                "phenotype_retention": "Reference AKI episodes still detected as AKI after thinning. It is a retention measure, not biological sensitivity.",
                "conditional_monitoring_indeterminate": "Indeterminate proportion only among thinned episodes still detected as AKI and retaining potential 48-h ICU coverage.",
                "total_phenotype_failure": "Fixed-reference denominator: not detected, detected after losing potential 48-h coverage, or detected but monitoring-indeterminate.",
                "hospital_rank": "Raw-versus-thinned failure-rate ranks are an unshrunk Monte Carlo sensitivity diagnostic only. The separate hierarchical model governs the hospital-comparability inference.",
            },
        }
        target = args.output_dir / f"eicu_v4_controlled_thinning_{hours}h.json"
        if target.exists():
            raise FileExistsError(f"Refusing to overwrite {target}")
        target.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        output_index.append({"hours": hours, "file": target.name, **{label: output["effect_decomposition"][label]["monte_carlo_median"] for label in METRICS.values()}})
    index = {"source": metadata["source"], "scope": "Index of aggregate-only 500-phase controlled-thinning outputs.", "schedules": output_index}
    target = args.output_dir / "eicu_v4_controlled_thinning_index.json"
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite {target}")
    target.write_text(json.dumps(index, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(index, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
