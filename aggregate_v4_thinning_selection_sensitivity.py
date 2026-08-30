#!/usr/bin/env python3
"""Aggregate first/nearest/last observation selection sensitivity outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


METRICS = ["retained", "primary_retained", "uncertain", "failure", "rho", "quartile_change"]


def summary(series: pd.Series) -> dict[str, float]:
    return {"median": round(float(series.median()), 6), "p025": round(float(series.quantile(.025)), 6), "p975": round(float(series.quantile(.975)), 6)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulation-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    table = []
    for rule in ("nearest", "first", "last"):
        for hours in (12, 24, 36, 48):
            data = pd.read_csv(args.simulation_dir / f"metrics_{rule}_{hours}.tsv", sep="\t")
            record: dict[str, object] = {"selection_rule": rule, "hours": hours, "replicates": len(data)}
            for metric in METRICS:
                record[metric] = summary(data[metric])
            table.append(record)
    output = {
        "scope": "Post-result selection-rule sensitivity within a controlled-thinning experiment; all schedules use the same fixed 9,323-episode observed reference cohort and 500 random global phases.",
        "rules": {"nearest": "primary: nearest observed value to bin centre", "first": "sensitivity: first observed value in each bin", "last": "sensitivity: last observed value in each bin"},
        "metrics": {"retained": "AKI phenotype retention", "primary_retained": "retention with potential 48-h coverage", "uncertain": "conditional monitoring-indeterminate proportion", "failure": "fixed-reference total phenotype failure", "rho": "unshrunk hospital failure-rate rank Spearman diagnostic", "quartile_change": "unshrunk hospital failure-rate quartile change diagnostic"},
        "results": table,
        "constraint": "This does not establish an optimal monitoring schedule or a hospital-quality rank.",
    }
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output}")
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"records": len(table)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
