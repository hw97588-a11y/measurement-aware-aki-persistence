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
