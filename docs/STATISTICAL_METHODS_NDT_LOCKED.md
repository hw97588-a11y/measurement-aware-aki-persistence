# Locked statistical methods for the NDT manuscript

## Analysis unit, observation scope and primary estimand

The unit of analysis was the first creatinine-defined acute kidney injury
episode in an adult patient's first continuous intensive care spell. MIMIC-IV
v3.1, SICdb v1.0.8 and eICU-CRD v2.0 were analysed separately under a
harmonized phenotype specification; patient-level data were not pooled. The
cross-database primary observation scope was limited to creatinine measurements
recorded during the continuous intensive care spell.

The primary denominator comprised episodes for which database intensive care
coverage continued through 48 hours after the first AKI-positive creatinine.
Inclusion depended on survival and continued database coverage, not on whether
a subsequent creatinine was measured. Episodes without this opportunity were
reported as structurally censored and were not combined with monitoring-related
classification indeterminacy.

The primary estimand was the database-specific proportion of observable
episodes whose feasible duration interval crossed 48 hours and therefore
could not be uniquely classified as transient or persistent under the recorded
creatinine schedule.

## Episode construction and interval-defined categories

At each creatinine measurement, AKI positivity was assessed using a rolling
minimum creatinine over the preceding 48 hours for the absolute KDIGO criterion
and over the preceding seven days for the relative criterion. The first
transition from a non-AKI to an AKI-positive observed state during the first
seven intensive care days defined the index episode. The episode recovery
baseline was then fixed as the minimum available creatinine in the seven days
before first positivity. Recovery required an observed creatinine below both
the fixed baseline plus 0.3 mg/dL and 1.5 times the fixed baseline.

The true onset lay between the last observed non-AKI value and first observed
AKI-positive value. The true recovery lay between the last observed AKI-positive
value and first observed recovery value. The minimum and maximum feasible
durations were derived from these intervals. Episodes with a maximum duration
of 48 hours or less were definite transient; those with a minimum duration
strictly greater than 48 hours were definite persistent; and those whose
duration interval crossed 48 hours were interval-indeterminate. Episodes with
no observed recovery and insufficient observed AKI-state duration to establish
persistence were right-censored unresolved. The two latter groups constituted
classification indeterminacy in the primary estimand.

The primary classification used a prespecified single-episode continuity
convention: after first positivity, the index episode was considered ongoing
until the first observed value meeting the episode-ending recovery rule. Thus,
the population bounds are logical bounds under this convention; unobserved
recovery followed by recurrent AKI between measurements cannot be excluded.

## Partial identification and sampling uncertainty

For threshold c=48 hours, the lower bound for the prevalence of persistent AKI was the
proportion with minimum feasible duration greater than c. The upper bound was
the proportion whose maximum feasible duration could exceed c, equivalent to
one minus the definite-transient proportion. Their difference equalled the
classification-indeterminate proportion. The identified set and its width
quantified identification uncertainty.

Sampling uncertainty was assessed with 2,000 bootstrap replicates using a
fixed seed. MIMIC-IV and SICdb used unique-patient cluster resampling; SICdb
clusters were defined by its de-identified PatientID. eICU first resampled
hospitals and then true unique patients (`uniquepid`) within sampled hospitals,
retaining all eligible episodes belonging to each sampled patient. Percentile
intervals were reported for each bound and the identified-set width. A
conservative unstudentized max-deviation bootstrap region was additionally
constructed to cover both endpoints simultaneously. The same cluster-
respecting bootstrap interval was used whenever the classification-
indeterminate proportion or identical identified-set width was reported.

## Robustness analyses

Phenotyping was repeated using two consecutive recovery measurements separated
by at least six hours and using recovery confirmed by a subsequent non-AKI
measurement 24–48 hours later. Further analyses used stricter ICU-acquired AKI,
initial Stage 2/3 AKI, alternative observed baseline strategies where available,
and persistence thresholds from 24 to 96 hours.

As an additional post hoc analysis of the continuity convention, originally
definite-persistent episodes were required to have an observed AKI-state chain
with adjacent creatinine measurements no more than 24 or 36 hours apart through
the measurement supporting persistence. Episodes violating the selected rule
were reclassified as continuity-gap indeterminate. The rule does not assert
that recovery occurred; it removes support for bridging a long unobserved gap.

## Observation process

Repeat testing was analysed using patient-clustered generalized estimating
equation logistic models in which each eligible measured creatinine was a
recurrent opportunity and the outcome was another measurement within 24
hours. Predictors available at the opportunity included the most recent
creatinine, preceding change, current observed AKI state and intensive care
day. Model coefficients were reported as odds ratios, accompanied by
marginally standardized retesting probabilities and risk differences. These
descriptive models used the same ICU-only scope and true unique-patient
clustering as the primary analysis. They were not used to correct, impute, or
replace the interval-defined primary persistence estimand.

## Controlled thinning

Controlled thinning used a fixed eICU measurement-rich reference cohort of
9,323 first AKI episodes from 184 hospitals, each with potential 48-hour ICU
coverage and at least four observed creatinines from first positivity through
72 hours. This was a maximally observed reference trajectory, not a biological
gold standard. Existing measurements were reduced under observation-grid
intervals of 12, 24, 36 and 48 hours. For each grid, 500 global phases
were independently sampled from a uniform distribution over one grid
interval. The observation nearest each phase-shifted bin centre was retained;
no values were imputed. First- and last-observation selection rules were
examined as sensitivity analyses.

For the fixed reference denominator, analyses separately quantified phenotype
retention, conditional classification indeterminacy among retained and
primary-eligible AKI episodes, total phenotype failure (non-retention or
indeterminacy), and the full reference-to-thinned category transition matrix.

## Multiplicity and software

The study had one primary descriptive estimand. Robustness, observation-process
and controlled-thinning analyses were supportive or exploratory; no
multiplicity-adjusted confirmatory claims were made. Analyses
used Python 3.13 with NumPy, SciPy, pandas and statsmodels. Raw source data
remained read-only on the governed data volume. Fifteen phenotype and
continuity unit tests and 47 final-revision reconciliation checks passed.
