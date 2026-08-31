# Measurement-aware identification of 48-hour AKI persistence

This repository contains the phenotype specification, analysis code, synthetic
tests and non-disclosive aggregate outputs supporting the study:

> Routine creatinine surveillance incompletely identifies acute kidney injury
> persistence at 48 hours

Repository: https://github.com/hw97588-a11y/measurement-aware-aki-persistence

## Release status and citation

Version 1.2.1 is the portable reproducibility release for the current
manuscript. It is intended to be tagged as `v1.2.1` only after all checks in
[`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md) pass on the exact
release commit. A version-specific Zenodo DOI must be added only after Zenodo
has archived that exact tagged release; do not reuse the DOI of an earlier
release. Until then, cite the repository URL and version shown in
[`CITATION.cff`](CITATION.cff).

The analysis evaluates whether irregular routine serum-creatinine measurements
allow each creatinine-defined acute kidney injury (AKI) episode to be uniquely
classified as lasting no more than 48 hours or more than 48 hours. Episodes
whose feasible duration interval crosses 48 hours remain
classification-indeterminate; consequently, the episode-level proportion
persisting beyond 48 hours is reported as a lower and upper bound.

## Data access and public-archive boundary

No patient-level or episode-level source data are included in this repository.
The study used the following controlled or contributor-reviewed databases:

- MIMIC-IV v3.1, PhysioNet, DOI: `10.13026/kpb9-mt58`
- eICU Collaborative Research Database v2.0, PhysioNet, DOI:
  `10.13026/C2WM1R`
- Salzburg Intensive Care database v1.0.8, PhysioNet, DOI:
  `10.13026/8m72-6j83`

Researchers must obtain each dataset through its own access procedure and
comply with its training, credentialing, contributor-approval and data-use
requirements. Never commit or upload extracted data, temporary trajectory
caches, patient identifiers, hospital identifiers, episode-level model inputs
or data-derived per-phase simulation files. The public archive contains only
source code, documentation and source-level aggregate outputs. See
[`docs/PUBLIC_ARCHIVE_POLICY.md`](docs/PUBLIC_ARCHIVE_POLICY.md) for the
release-boundary audit.

## Repository contents

- `run_interval_aki_primary.py`: source mapping, cleaning utilities and
  database loaders.
- `interval_aki_v4_engine.py`: ICU-coverage interval phenotype engine. Index
  AKI is searched during ICU days 0–7, while recovery follow-up continues to
  the end of the continuous database-covered critical-care spell.
- `run_v4_primary_inference.py`: primary categories, partial-identification
  bounds, threshold curves and unique-patient-cluster bootstrap inference;
  eICU uses hospital-to-patient two-stage resampling.
- `run_v5_core_sensitivities.py`: corrected ICU-only recovery,
  intensive-care-acquired AKI, measured-baseline, observation-window and
  thinning-reference-cohort analyses.
- `mimic_history_helper.py`: read-only helper for measured MIMIC preadmission
  creatinine baselines used by the core-sensitivity module.
- `run_v5_observation_process.py`: ICU-only, true-patient-clustered recurrent
  retesting models. These models are descriptive and do not impute persistence
  status.
- `run_ndt_continuity_gap_sensitivity.py`: 24- and 36-hour observed-unrecovered
  chain continuity stress tests.
- `cache_v4_thinning_reference.py`, `controlled_thinning_sim.cpp` and
  `aggregate_v4_controlled_thinning.py`: fixed-reference, random-phase,
  index-episode-anchored controlled thinning. A later recurrent AKI cannot
  substitute for a missed index episode. The temporary cache is protected and
  must never be shared.
- `audit_v6_targeted_outputs.py`: 44 cross-file, time-window, denominator,
  clustering, index-episode-matching and arithmetic checks over the released
  aggregate outputs.
- `prepare_public_aggregate_release.py`: removes legacy aggregate
  hospital-ranking diagnostics that are not part of the current manuscript
  before the public source archive is built.
- `test_*.py`: data-free synthetic tests.
- `docs/`: phenotype, statistical and release specifications.
- `results/`: non-disclosive aggregate outputs supporting retained manuscript
  tables and robustness analyses.

## Clean-room setup

The frozen execution environment was Python 3.13.12 with the versions pinned
in [`requirements.txt`](requirements.txt). Create an isolated environment from
a fresh clone:

```bash
git clone https://github.com/hw97588-a11y/measurement-aware-aki-persistence.git
cd measurement-aware-aki-persistence
git checkout v1.2.1
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The source loaders also require `gzip`, `rg` and `bsdtar`. A C++17 compiler is
needed only for the accelerated controlled-thinning simulation. See
[`docs/ENVIRONMENT_FREEZE_v1.2.1.md`](docs/ENVIRONMENT_FREEZE_v1.2.1.md) for
the exact reference environment and no-data verification commands.

## Data-free verification

These commands require no restricted data and should all succeed from a fresh
clone before any governed source data are configured:

```bash
python -m compileall -q .
python run_interval_aki_primary.py --help
python run_v4_primary_inference.py --help
python run_v5_core_sensitivities.py --help
python run_v5_observation_process.py --help
python run_v4_controlled_thinning.py --help
python run_ndt_continuity_gap_sensitivity.py --help
python cache_v4_thinning_reference.py --help
python aggregate_v4_controlled_thinning.py --help
python prepare_public_aggregate_release.py --check
python -m unittest -v test_interval_aki_v4_engine.py test_ndt_continuity_gap_sensitivity.py
python audit_v6_targeted_outputs.py
```

Twenty-four synthetic tests cover the 48-hour boundary, late index onset, duplicate
measurements, unresolved recovery, recurrent AKI, recovery confirmation,
SICdb endpoint re-anchoring, index-episode-anchored thinning, patient-specific
phase selection and continuity-gap logic. The aggregate audit records 44
passed reconciliation checks.

## Configure controlled data locally

Set data paths outside the repository; the runtime validates each configured
path before reading source data:

```bash
export MIMIC_IV_PATH=/secure/path/to/mimic-iv-3.1.zip
export EICU_CRD_PATH=/secure/path/to/eicu-crd-2.0
export SICDB_PATH=/secure/path/to/sicdb-1.0.8.rar
# Optional when the SICdb archive has a different top-level member directory:
export SICDB_MEMBER_ROOT=salzburg-intensive-care-database-sicdb-a-freely-accessible-intensive-care-database-1.0.8
```

The expected directory and archive structures match the original PhysioNet
releases. Never place restricted source files inside the repository. Use a
protected output directory outside the repository for all regenerated files.

## Primary analysis

Run each database separately; outputs are never pooled at patient level:

```bash
python run_v4_primary_inference.py --database mimic --output-dir /secure/output/mimic
python run_v4_primary_inference.py --database sicdb --output-dir /secure/output/sicdb
python run_v4_primary_inference.py --database eicu --output-dir /secure/output/eicu
```

The primary denominator depends on survival and ICU database coverage through
48 hours after first AKI positivity, not on whether another creatinine was
actually measured. Creatinine history is retained from seven days before ICU
entry; index AKI is searched only during ICU days 0–7; and recovery follow-up
continues from index onset to spell end. In SICdb, laboratory time and spell
end are both re-anchored to first ICU bed assignment.

## Controlled thinning

The fixed eICU reference cohort is cached locally only to make 500 random-phase
replicates computationally feasible. The cache contains protected trajectories
and is ignored by Git; store it in a private temporary directory and destroy it
according to the governing data-use agreement after use. Only aggregate JSON
outputs belong in `results/thinning/`.

The measurement-rich reference trajectory is not a biological gold standard,
and controlled thinning is not a simulated clinical testing policy. The code
reports index-episode retention, conditional indeterminacy and total phenotype
failure without causal interpretation. Later recurrent AKI is recorded
separately and cannot rescue a missed index episode.

## Frozen central results

Among first observed transition-defined AKI episodes with 48-hour potential ICU
observation, persistence status was classification-indeterminate in 21.3% of
MIMIC-IV, 31.8% of SICdb and 30.8% of eICU episodes. Cluster-respecting
bootstrap 95% confidence intervals were 20.6%–22.2%, 29.9%–33.4% and
29.5%–32.3%, respectively. Identified sets for the episode-level proportion
persisting beyond 48 hours were 41.2%–62.6%, 41.5%–73.2% and 37.8%–68.6%.
In the fixed 9,790-episode eICU reference cohort, total phenotype failure was
54.8% with a 24-hour observation grid and 77.5% with a 48-hour observation
grid.

## Reproducibility boundaries

- Primary bounds are logical under the prespecified single-episode continuity
  convention; they are not assumption-free biological bounds.
- No missing creatinine is interpreted as recovery or persistence.
- Death, discharge and ICU departure are not encoded as non-recovery.
- Kidney replacement therapy is not automatically encoded as persistent AKI
  because timing and capture are not transportable across the three sources.
- The final manuscript excludes the earlier mortality landmark,
  inverse-observation weighting and hospital-ranking modules.
- MIMIC-IV and SICdb inference resamples unique patients. eICU inference first
  resamples hospitals and then unique patients within sampled hospitals,
  retaining all eligible episodes for each sampled patient.

## Authors and contact

Cheng Shen, Bohao Xue and Jin Li. Correspondence: Jin Li,
`leesunny2015@163.com`.

## Licence and citation

Code is released under the MIT License. Please cite this repository using
[`CITATION.cff`](CITATION.cff) and cite all three source datasets under their
exact versions and persistent identifiers. When Zenodo has archived the exact
`v1.2.1` tag, add the version-specific DOI to the manuscript, release note and
repository metadata without changing the tagged source archive.
