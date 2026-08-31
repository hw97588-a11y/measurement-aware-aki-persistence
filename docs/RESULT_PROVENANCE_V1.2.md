# Result provenance for version 1.2.0

This manifest maps every retained manuscript table and figure to its
non-disclosive aggregate source and the script that creates that source. The
Word and plotting layers only format these frozen aggregates; they do not
re-estimate any phenotype or model.

| Manuscript output | Frozen aggregate source | Analysis script |
|---|---|---|
| Table 1, clinical and observation characteristics | `results/primary/*_v4_primary_inference.json`; `results/robustness/*_v5_core_sensitivities.json`; `results/observation/*_v5_observation_process.json` | `run_v4_primary_inference.py`; `run_v5_core_sensitivities.py`; `run_v5_observation_process.py` |
| Table 2, primary persistence results | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` |
| Figure 1, flow and interval classification | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` |
| Figure 2, identified sets | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` |
| Figure 3, controlled thinning | `results/thinning/eicu_v4_controlled_thinning_{12,24,36,48}h.json` | `run_v4_controlled_thinning.py`; `controlled_thinning_sim.cpp`; `aggregate_v4_controlled_thinning.py` |
| Supplementary Tables S1–S3 | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` |
| Supplementary Table S4 | `results/robustness/*_v5_core_sensitivities.json` | `run_v5_core_sensitivities.py` |
| Supplementary Table S5 and Figure S1 | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` |
| Supplementary Table S6 | `results/observation/*_v5_observation_process.json` | `run_v5_observation_process.py` |
| Supplementary Tables S7–S8 | `results/thinning/eicu_v4_controlled_thinning_{12,24,36,48}h.json`; `results/robustness/eicu_v5_core_sensitivities.json` | `run_v4_controlled_thinning.py`; `run_v5_core_sensitivities.py` |
| Supplementary Table S9 and Figure S3 | `results/continuity/*_ndt_continuity_gap_sensitivity.json` | `run_ndt_continuity_gap_sensitivity.py` |
| Supplementary Table S10 | `results/v6_targeted_reanalysis_audit.json` and unit-test logs | `audit_v6_targeted_outputs.py`; `test_interval_aki_v4_engine.py`; `test_ndt_continuity_gap_sensitivity.py` |
| Supplementary Figure S2 | `results/primary/*_v4_primary_inference.json` | `run_v4_primary_inference.py` |
| Supplementary Figure S4 | `results/robustness/*_v5_core_sensitivities.json` | `run_v5_core_sensitivities.py` |

The protected patient-, episode-, hospital- and phase-level caches are excluded
from this repository. `results/v6_targeted_reanalysis_audit.json` records 44
passed reconciliation checks over the retained aggregate outputs.
