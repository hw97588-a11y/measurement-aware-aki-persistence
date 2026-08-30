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

## Revision audit

- `results/v5_revision_output_audit.json`

This file records 47 passed consistency checks. The final submission freeze does
not include the earlier mortality landmark, inverse-observation-weighting, or
hospital-ranking outputs.

## Episode-continuity sensitivity

- `results/continuity/*_ndt_continuity_gap_sensitivity.json`
- `results/continuity/ndt_final_gate_audit.json`

These files report only source-level counts and proportions under 24- and
36-hour maximum observed-positive-chain gaps.

## Controlled thinning

- `results/thinning/eicu_v4_controlled_thinning_{12,24,36,48}h.json`
- `results/thinning/eicu_v4_controlled_thinning_index.json`
- `results/thinning/eicu_v4_controlled_thinning_selection_rule_sensitivity.json`

These files contain Monte Carlo summaries and aggregated transition counts
over 500 phases. The protected trajectory cache and phase-level files are not
distributed.

## Manuscript tables

- `results/tables/Table_1_cohort_and_observation_characteristics.csv`
- `results/tables/Table_2_primary_persistence_results.csv`
- `results/tables/Table_S9_episode_continuity_sensitivity.csv`

The output-to-code-to-table mapping is documented in
`docs/RESULT_PROVENANCE_V1.1.md`.
