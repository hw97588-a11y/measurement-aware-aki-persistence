# Database crosswalk

The loaders map each source to a common sequence of `(minutes from first ICU
entry, creatinine in mg/dL)` observations and one continuous ICU spell.

| Construct | MIMIC-IV v3.1 | SICdb v1.0.8 | eICU-CRD v2.0 |
|---|---|---|---|
| Patient/admission unit | `subject_id`, `hadm_id`, joined consecutive `icustays` | `PatientID`, `CaseID`, `ICUOffset` | `uniquePid`, `patientHealthSystemStayId`, joined consecutive `patientUnitStayId` |
| ICU boundary | `icustays.intime/outtime` | `cases.ICUOffset` and `TimeOfStay` | `patient.hospitalAdmitOffset` and `unitDischargeOffset` |
| Creatinine source | `hosp.labevents` item IDs 50912 or 52546 | `laboratory` IDs 367 or 368 | `lab`, `labName=creatinine`, `labTypeID=1` |
| Clinical time | `labevents.charttime` | `laboratory.Offset` relative to `ICUOffset` | `lab.labResultOffset` aligned across joined units |
| Units | mg/dL | source values used after distribution/unit audit | mg/dL; µmol/L or umol/L divided by 88.4 |
| Hospital cluster | single centre | not used for primary resampling | `hospitalId` |
| Patient bootstrap cluster | `subject_id` | `PatientID` | hospital, then true `uniquepid` within hospital; all eligible episodes for each sampled patient are retained |

Only central laboratory serum/plasma creatinine is used. Duplicate values at
the same timestamp are reduced to their median. Values outside 0.1–25 mg/dL
are excluded. eICU hospitals with incompatible creatinine unit interfaces are
excluded as a whole rather than repaired at patient level.

The cross-database primary observation scope ends at departure from the first
continuous ICU spell. MIMIC hospital-wide laboratory follow-up is a labelled
sensitivity analysis and is not substituted for the ICU-only primary result.
