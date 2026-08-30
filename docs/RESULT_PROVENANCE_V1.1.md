# Result provenance for the NDT v1.1 freeze

This manifest links each retained manuscript result to an aggregate source file
and the code that generates it. Restricted patient-, episode-, hospital-, and
trajectory-level data are not distributed.

| Manuscript output | Aggregate source | Analysis code | Presentation code |
|---|---|---|---|
| Main Table 1, cohort flow and observation characteristics | `results/primary/*_v4_primary_inference.json`; `results/robustness/*_v5_core_sensitivities.json` | `run_interval_aki_primary.py`; `run_v4_primary_inference.py`; `run_v5_core_sensitivities.py` | `results/tables/Table_1_cohort_and_observation_characteristics.csv`; submission document builder |
| Main Table 2, classification-indeterminate proportions and persistence bounds | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` | `results/tables/Table_2_primary_persistence_results.csv`; submission document builder |
| Main Figure 1, source flow and interval-classification concept | `results/robustness/*_v5_core_sensitivities.json`; phenotype specification | `run_v5_core_sensitivities.py`; `interval_aki_v4_engine.py` | submission figure builder |
| Main Figure 2, database-specific identified sets | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` | submission figure builder |
| Main Figure 3, controlled thinning | `results/thinning/eicu_v4_controlled_thinning_{12,24,36,48}h.json`; `results/thinning/eicu_v4_controlled_thinning_index.json` | `run_v4_controlled_thinning.py` | submission figure builder |
| Supplementary Tables S1–S3, denominator, inference, and source flow | `results/primary/*_v4_primary_inference.json`; `results/robustness/*_v5_core_sensitivities.json` | v4 primary and v5 sensitivity scripts | submission document builder |
| Supplementary Table S4, corrected observation-process model | `results/observation/*_v5_observation_process.json` | `run_v5_observation_process.py` | submission document builder |
| Supplementary Tables S5–S8, recovery, ICU-acquired, baseline/window, and thinning-cohort sensitivities | `results/robustness/*_v5_core_sensitivities.json` | `run_v5_core_sensitivities.py` | submission document builder |
| Supplementary Table S9 and Figure S3, episode-continuity convention | `results/continuity/*_ndt_continuity_gap_sensitivity.json` | `run_ndt_continuity_gap_sensitivity.py` | `results/tables/Table_S9_episode_continuity_sensitivity.csv`; submission builders |
| Supplementary Tables S10–S11 and Figure S4, thinning transfer and selection-rule sensitivities | `results/thinning/eicu_v4_controlled_thinning_*.json` | `run_v4_controlled_thinning.py` | submission builders |
| Supplementary Figures S1–S2, threshold curve and indeterminacy decomposition | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` | submission figure builder |

The retained submission does not use the earlier mortality-landmark,
inverse-observation-weighting, or centre-ranking modules. The machine-readable
reconciliation record is `results/v5_revision_output_audit.json` (47/47 checks
passed).

