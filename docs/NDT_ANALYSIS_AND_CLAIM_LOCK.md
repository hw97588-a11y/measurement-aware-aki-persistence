# NDT analysis and claim lock

Date: 2026-08-30

Target: *Nephrology Dialysis Transplantation*, Original Article.

## Central claim

Routine creatinine surveillance leaves the 48-hour persistence status of a
substantial proportion of creatinine-defined AKI episodes indeterminate.
Consequently, the population prevalence of persistent AKI is bounded rather than
uniquely point-estimated, and phenotype comparability deteriorates when the
recorded observation schedule is thinned.

## Required terminology

- Episode level: `definite transient`, `definite persistent`, and
  `classification-indeterminate under the observed creatinine schedule`.
- Population level: `identified set for persistent AKI prevalence`, `lower bound`,
  and `upper bound`.
- Primary bounds: `logical bounds under the prespecified episode-continuity
  convention`.
- Simulation: `controlled thinning`, `reduced-resolution observation
  schedule`, `phenotype retention`, and `total phenotype failure`.
- Database relationship: `cross-database replication under a harmonized
  phenotype specification`.

Do not use `misclassified`, `true prevalence`, `diagnostic sensitivity`,
`standardized clinical monitoring`, `external validation`, `hospital quality
ranking`, or causal language.

## Frozen analysis hierarchy

### Primary

For each database separately, among first creatinine-defined AKI episodes with
ICU coverage through first positivity plus 48 hours:

1. proportion classification-indeterminate;
2. lower and upper bounds for the prevalence of persistent AKI;
3. 2,000-replicate cluster-respecting sampling intervals/confidence region.

No patient-level cross-database pooling.

### Key secondary

1. confirmed recovery definitions;
2. strict ICU-acquired AKI;
3. baseline-creatinine strategies;
4. recurrent observation-process model using true unique-patient clusters;
5. eICU fixed-reference controlled thinning.

### Additional post hoc robustness

1. 24- and 36-hour maximum observed-positive-chain gap sensitivity;
2. 24–96-hour threshold curves;
3. controlled-thinning selection-rule sensitivity.

## Final gates

- Primary denominator depends only on survival and ICU database coverage, not
  on actual repeat testing: PASS.
- Structural censoring, monitoring indeterminacy and KRT interruption are not
  conflated: PASS, subject to source-specific KRT availability described in
  the supplement.
- Cluster-respecting partial-identification inference: PASS.
- Controlled thinning uses one fixed reference cohort and separates AKI
  non-retention from conditional indeterminacy: PASS.
- Controlled thinning uses 500 random phases at every schedule and no
  imputation: PASS.
- Primary engine synthetic boundary tests and reconciliation audit: PASS
  (15/15 unit tests; 47/47 revision checks).
- eICU inference resamples hospitals and then true unique patients within
  sampled hospitals: PASS.
- Single-episode continuity convention disclosed and 24/36-hour gap
  sensitivity completed: PASS.

## Results that must appear in the abstract

- Classification-indeterminate: 25.0% MIMIC-IV, 36.2% SICdb, 34.0% eICU.
- Identified sets for the prevalence of persistent AKI: 38.3%–63.3%,
  37.4%–73.6%, and 35.0%–69.0%, respectively.
- In eICU controlled thinning, fixed-reference total phenotype failure:
  54.9% under a 24-hour schedule and 78.9% under a 48-hour schedule.

Continuity-gap results belong in the Results robustness paragraph and
supplement, not the abstract unless requested by an editor or reviewer.

## Final interpretive boundary

This study concerns the computability and cross-centre comparability of
creatinine-defined 48-hour persistence under routine observation. It does not
test whether more frequent creatinine testing improves outcomes, whether the
48-hour clinical construct is biologically correct, or whether any hospital
provides higher-quality care.
