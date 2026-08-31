# NDT analysis and claim lock

Date: 2026-08-31

Target: *Nephrology Dialysis Transplantation*, Original Article.

## Central claim

Routine creatinine surveillance leaves the 48-hour persistence status of a
substantial proportion of creatinine-defined AKI episodes indeterminate.
Consequently, the episode-level proportion persisting beyond 48 hours is
bounded rather than uniquely point-estimated, and phenotype comparability
deteriorates when the recorded observation schedule is thinned.

## Required terminology

- Episode level: `definite transient`, `definite persistent`, and
  `classification-indeterminate under the observed creatinine schedule`.
- Population level: `identified set for the episode-level proportion persisting
  beyond 48 hours`, `lower bound`, and `upper bound`.
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
2. lower and upper bounds for the episode-level proportion persisting beyond
   48 hours;
3. 2,000-replicate cluster-respecting sampling intervals/confidence region.

No patient-level cross-database pooling.

### Key secondary

1. confirmed recovery definitions;
2. strict ICU-acquired AKI;
3. baseline-creatinine strategies;
4. recurrent observation-process model using true unique-patient clusters;
5. eICU fixed-reference controlled thinning.

### Additional post hoc robustness

1. 24- and 36-hour maximum unrecovered-state-support gap sensitivity;
2. 24–96-hour threshold curves;
3. initial Stage 1 versus Stage 2/3 and one-episode-per-patient sensitivity.

## Final gates

- Primary denominator depends only on survival and ICU database coverage, not
  on actual repeat testing: PASS.
- Insufficient 48-hour observation opportunity, monitoring indeterminacy and
  KRT interruption are not conflated: PASS, subject to source-specific KRT
  availability described in the supplement.
- Cluster-respecting partial-identification inference: PASS.
- Controlled thinning uses one fixed reference cohort, is anchored to the
  original index episode and separates index-episode non-retention from
  conditional indeterminacy: PASS.
- Controlled thinning uses 500 random phases at every schedule and no
  imputation: PASS.
- Baseline, index-search and recovery-follow-up windows are distinct; SICdb
  laboratory and coverage endpoints share the ICUOffset origin: PASS.
- Primary engine synthetic boundary tests and reconciliation audit: PASS
  (24/24 unit tests: 20 phenotype-engine and 4 continuity-support tests;
  44/44 final-revision checks).
- eICU inference resamples hospitals and then true unique patients within
  sampled hospitals: PASS.
- Single-episode continuity convention disclosed and 24/36-hour gap
  sensitivity completed: PASS.

## Results that must appear in the abstract

- Classification-indeterminate: 21.3% MIMIC-IV, 31.8% SICdb, 30.8% eICU.
- Identified sets for the episode-level proportion persisting beyond 48 hours:
  41.2%–62.6%, 41.5%–73.2%, and 37.8%–68.6%, respectively.
- In eICU controlled thinning, fixed-reference total phenotype failure:
  54.8% under a 24-hour grid and 77.5% under a 48-hour grid.

Continuity-gap results belong in the Results robustness paragraph and
supplement, not the abstract unless requested by an editor or reviewer.

## Final interpretive boundary

This study concerns the computability and cross-centre comparability of
creatinine-defined 48-hour persistence under routine observation. It does not
test whether more frequent creatinine testing improves outcomes, whether the
48-hour clinical construct is biologically correct, or whether any hospital
provides higher-quality care.
