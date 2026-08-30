# Measurement-aware identification of 48-hour AKI persistence

This repository contains the phenotype specification, analysis code, synthetic
tests, and non-disclosive aggregate outputs supporting the study:

> Routine creatinine surveillance incompletely identifies acute kidney injury
> persistence at 48 hours

Repository: https://github.com/hw97588-a11y/measurement-aware-aki-persistence

The analysis evaluates whether irregular routine creatinine measurements allow
each creatinine-defined acute kidney injury (AKI) episode to be uniquely
classified as lasting no more than 48 hours or more than 48 hours. Episodes
whose feasible duration interval crosses 48 hours remain classification-
indeterminate; consequently, the population prevalence of persistent AKI is reported
as a lower and upper bound.

## Data access

No patient-level or episode-level source data are included in this repository.
The study used the following controlled or contributor-reviewed databases:

- MIMIC-IV v3.1, PhysioNet, DOI: `10.13026/kpb9-mt58`
- eICU Collaborative Research Database v2.0, PhysioNet, DOI:
  `10.13026/C2WM1R`
- Salzburg Intensive Care database v1.0.8, PhysioNet, DOI:
  `10.13026/8m72-6j83`

Researchers must obtain each dataset through its own access procedure and
comply with its training, credentialing, contributor-approval, and data-use
requirements. Do not commit extracted data, temporary trajectory caches,
patient identifiers, hospital identifiers, or episode-level model inputs.

## Repository contents

- `run_interval_aki_primary.py`: database-specific source mapping and common
  cleaning utilities.
- `interval_aki_v4_engine.py`: ICU-coverage interval phenotype engine.
- `run_v4_primary_inference.py`: primary categories, partial-identification
  bounds, threshold curves, and unique-patient-cluster bootstrap inference;
  eICU uses hospital-to-patient two-stage resampling.
- `run_v5_core_sensitivities.py`: corrected ICU-only recovery,
  intensive-care-acquired AKI, measured-baseline, observation-window and
  thinning-reference-cohort analyses.
- `mimic_history_helper.py`: read-only helper for measured MIMIC preadmission
  creatinine baselines used by the corrected sensitivity module.
- `run_v5_observation_process.py`: corrected ICU-only, true-patient-clustered
  recurrent retesting models. These models are descriptive and do not impute
  persistence status.
- `audit_v5_revision_outputs.py`: 47 cross-file, denominator, clustering and
  arithmetic checks used for the revised submission freeze.
- `run_ndt_continuity_gap_sensitivity.py`: 24- and 36-hour observed-positive-
  chain continuity stress tests.
- `cache_v4_thinning_reference.py`, `controlled_thinning_sim.cpp`, and
  `aggregate_v4_controlled_thinning.py`: fixed-reference, random-phase
  controlled thinning. The temporary cache is protected and must never be
  shared.
- `tests` are supplied as root-level `test_*.py` files for compatibility with
  the frozen scripts.
- `docs/`: locked phenotype and statistical specifications.
- `results/`: aggregate outputs supporting the primary tables and robustness
  analyses.

## Environment

Python 3.13 was used for the frozen analysis. Install dependencies in an
isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

The source loaders also require `gzip`, `rg`, and `bsdtar`. A C++17 compiler is
needed only for the accelerated controlled-thinning simulation.

Set paths without editing source code:

```bash
export MIMIC_IV_PATH=/path/to/mimic-iv-3.1.zip
export EICU_CRD_PATH=/path/to/eicu-crd-2.0
export SICDB_PATH=/path/to/sicdb-1.0.8.rar
```

The expected directory and archive structures match the original PhysioNet
releases. `SICDB_MEMBER_ROOT` can be overridden if the archive uses a different
top-level member name.

## Synthetic tests

The tests require no restricted data:

```bash
python3 -m unittest -v test_interval_aki_v4_engine.py
python3 -m unittest -v test_ndt_continuity_gap_sensitivity.py
```

They cover the 48-hour boundary, duplicate measurements, unresolved recovery,
recurrent AKI, recovery confirmation, and continuity-gap logic.

## Primary analysis

Run each database separately; outputs are never pooled at patient level:

```bash
python3 run_v4_primary_inference.py --database mimic --output-dir outputs/mimic
python3 run_v4_primary_inference.py --database sicdb --output-dir outputs/sicdb
python3 run_v4_primary_inference.py --database eicu --output-dir outputs/eicu
```

The primary denominator depends on survival and ICU database coverage through
48 hours after first AKI positivity, not on whether another creatinine was
actually measured.

## Controlled thinning

The fixed eICU reference cohort is cached locally only to make 500 random-phase
replicates computationally feasible. The cache contains protected trajectories
and is ignored by Git; store it in a private temporary directory and destroy it
according to the governing data-use agreement after use. Only aggregate JSON
outputs belong in `results/thinning/`.

The measurement-rich reference trajectory is not a biological gold standard,
and controlled thinning is not a simulated clinical testing policy. The code
reports phenotype retention, conditional indeterminacy, and total phenotype
failure without causal interpretation.

## Frozen central results

Among first AKI episodes with 48-hour potential ICU observation, persistence
status was classification-indeterminate in 25.0% of MIMIC-IV, 36.2% of SICdb,
and 34.0% of eICU episodes. Patient-cluster bootstrap 95% confidence intervals
were 24.2%–25.8%, 34.3%–38.1%, and 32.7%–35.3%, respectively. The identified sets for the prevalence of persistent AKI were
38.3%–63.3%, 37.4%–73.6%, and 35.0%–69.0%, respectively. In the fixed eICU
reference cohort, total phenotype failure was 54.9% with a 24-hour observation
grid and 78.9% with a 48-hour observation grid.

## Reproducibility boundaries

- Primary bounds are logical under the prespecified single-episode continuity
  convention; they are not assumption-free biological bounds.
- No missing creatinine is interpreted as recovery or persistence.
- Death, discharge, and ICU departure are not encoded as non-recovery.
- Kidney replacement therapy is not automatically encoded as persistent AKI
  because timing and capture are not transportable across the three sources.
- The final manuscript excludes the earlier mortality landmark, inverse-
  observation weighting and hospital-ranking modules. Those exploratory
  analyses are not part of the v1.1.0 submission freeze.
- MIMIC-IV and SICdb inference resamples unique patients. eICU inference first
  resamples hospitals and then unique patients within sampled hospitals,
  retaining all eligible episodes for each sampled patient.

## Authors and contact

Cheng Shen, Bohao Xue, and Jin Li. Correspondence: Jin Li,
`leesunny2015@163.com`.

## Licence and citation

Code is released under the MIT License. Please cite this repository using
`CITATION.cff` and cite all three source datasets under their exact versions
and persistent identifiers. `.zenodo.json` contains deposit-ready metadata;
the DOI field will be added after the authors create the Zenodo archive.
