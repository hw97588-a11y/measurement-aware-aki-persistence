# Aggregate output manifest

All files below contain database-level summaries only. No patient-level,
episode-level, or hospital-identifying rows are distributed.

## Primary inference

- `results/primary/mimic_v4_primary_inference.json`
- `results/primary/sicdb_v4_primary_inference.json`
- `results/primary/eicu_v4_primary_inference.json`

These files contain cohort counts, four phenotype categories, 48-hour
identified sets, cluster-respecting bootstrap intervals, threshold curves, and
aggregate decomposition of indeterminacy.

MIMIC-IV and SICdb use true unique-patient clusters. eICU uses hospital-to-
unique-patient two-stage resampling and reports globally unique patient and
hospital–patient cluster counts.

## Corrected core robustness

- `results/robustness/mimic_v5_core_sensitivities.json`
- `results/robustness/sicdb_v5_core_sensitivities.json`
- `results/robustness/eicu_v5_core_sensitivities.json`

These files contain recovery-confirmation, strict ICU-acquired AKI, measured-
baseline and observation-window sensitivities, plus the aggregate comparison
of the full eICU primary cohort with the measurement-rich thinning reference
cohort.

## Corrected observation process

- `results/observation/mimic_v5_observation_process.json`
- `results/observation/sicdb_v5_observation_process.json`
- `results/observation/eicu_v5_observation_process.json`

These files contain only aggregate opportunities, unique-patient counts,
patient-clustered generalized estimating equation parameters, and standardized
descriptive retesting probabilities. No inverse-observation weighting result is
part of the submission freeze.

## Audit bundled with version 1.2.1

- `results/v6_targeted_reanalysis_audit.json`

This file records 44 passed consistency checks covering distinct time windows,
SICdb endpoint alignment, true patient clusters, index-episode-anchored
thinning, flow identities and cross-file agreement. The final submission freeze
does not include the earlier mortality landmark, inverse-observation-weighting,
hospital-ranking or unanchored thinning outputs. Version 1.2.1 retains these
same non-disclosive aggregates and adds no patient-, episode-, hospital- or
phase-level records.

## Episode-continuity sensitivity

- `results/continuity/*_ndt_continuity_gap_sensitivity.json`
These files report only source-level counts and proportions under 24- and
36-hour maximum unrecovered-state-support gaps.

## Controlled thinning

- `results/thinning/eicu_v4_controlled_thinning_{12,24,36,48}h.json`
- `results/thinning/controlled_thinning_progress.json`

These files contain index-episode-anchored Monte Carlo summaries and aggregate
transition counts over 500 phases. Later recurrence cannot substitute for a
missed index episode. The protected trajectory cache and phase-level files are
not distributed. The public release also removes legacy aggregate
hospital-ranking diagnostics because they are not part of the retained
manuscript analyses; `prepare_public_aggregate_release.py --check` verifies
their absence.

## Manuscript tables

- `results/tables/Table_1_cohort_and_observation_characteristics.csv`
- `results/tables/Table_2_primary_persistence_results.csv`

The output-to-code-to-table mapping is documented in
`docs/RESULT_PROVENANCE_V1.2.md`.
