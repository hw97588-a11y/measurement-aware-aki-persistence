#!/usr/bin/env python3
"""Create a temporary, non-public compact trajectory cache for thinning.

The cache belongs in a private temporary directory and must not be committed
or shared. It contains de-identified intermediate trajectories only to make
the 500-replicate controlled-thinning computation feasible.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from run_interval_aki_primary import deduplicate
from run_v4_controlled_thinning import source_reference


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()
    if args.cache.exists() or args.metadata.exists():
        raise FileExistsError("Refusing to overwrite a temporary trajectory cache")
    source, spells, labs, reference = source_reference()
    hospitals = {hospital: index for index, hospital in enumerate(sorted({str(episode.hospital_id) for episode in reference.values()}))}
    hospital_counts = Counter(str(episode.hospital_id) for episode in reference.values())
    reference_counts = Counter(episode.category for episode in reference.values())
    category_codes = {"definite_transient": 0, "definite_persistent": 1, "interval_indeterminate": 2, "right_censored_unresolved": 3}
    with args.cache.open("w", encoding="utf-8") as handle:
        for identifier, episode in reference.items():
            trajectory = ";".join(f"{time:.8f},{value:.8f}" for time, value in deduplicate(labs[identifier]))
            handle.write(f"{hospitals[str(episode.hospital_id)]}\t{category_codes[episode.category]}\t{spells[identifier].end:.8f}\t{trajectory}\n")
    args.metadata.write_text(json.dumps({
        "source": source,
        "cache_scope": "Temporary protected intermediate file; do not distribute.",
        "episodes": len(reference), "hospitals": len(hospitals),
        "hospital_reference_counts": {str(hospitals[hospital]): count for hospital, count in hospital_counts.items()},
        "reference_categories": dict(reference_counts),
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"episodes": len(reference), "hospitals": len(hospitals), "cache": str(args.cache)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
