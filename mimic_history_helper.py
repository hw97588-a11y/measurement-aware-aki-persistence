#!/usr/bin/env python3
"""Read-only helper for measured preadmission MIMIC creatinine baselines."""

from __future__ import annotations

import csv
import gzip
import zipfile
from collections import defaultdict
from datetime import timedelta

from run_interval_aki_primary import MIMIC, as_datetime, as_float, require_source_path, valid_creatinine


def historical_baselines(spells) -> dict[str, list[float]]:
    """Extract same-patient creatinine values 7–365 days before admission."""
    by_subject: dict[str, list[str]] = defaultdict(list)
    for admission_id, spell in spells.items():
        by_subject[str(spell.extra["subject_id"])].append(admission_id)
    result: dict[str, list[float]] = defaultdict(list)
    with zipfile.ZipFile(require_source_path(MIMIC, "MIMIC_IV_PATH")) as archive:
        raw = gzip.GzipFile(fileobj=archive.open("mimic-iv-3.1/hosp/labevents.csv.gz"))
        next(raw, b"")
        for line in raw:
            if b",50912," not in line and b",52546," not in line:
                continue
            row = next(csv.reader([line.decode("utf-8", errors="replace")]))
            if len(row) < 11 or row[10].strip().casefold() != "mg/dl":
                continue
            observed, value = as_datetime(row[6]), as_float(row[9])
            if observed is None or not valid_creatinine(value):
                continue
            for admission_id in by_subject.get(row[1], []):
                admitted = spells[admission_id].extra["admit_dt"]
                if admitted - timedelta(days=365) <= observed < admitted - timedelta(days=7):
                    result[admission_id].append(value)
        raw.close()
    return result
