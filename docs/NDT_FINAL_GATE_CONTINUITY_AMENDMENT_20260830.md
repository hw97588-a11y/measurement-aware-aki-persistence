# NDT final-gate amendment: observed-positive-chain continuity sensitivity

Date: 2026-08-30

Status: additional post hoc methodological robustness analysis. This amendment
does not replace or alter the frozen v4 primary phenotype.

## Rationale

The v4 primary analysis uses a prespecified single-episode continuity
convention: after the first observed creatinine-defined AKI value, the indexed
episode is considered ongoing until the first observed creatinine value that
meets the fixed episode-ending recovery rule. Consequently, two AKI-state
measurements separated by a long unobserved interval could support a definite
persistent classification even though an unobserved recovery followed by a
recurrent AKI cannot be excluded.

The primary bounds must therefore be described as logical bounds under the
prespecified episode-continuity convention, not as assumption-free bounds.

## Fixed primary convention

- The AKI-detection baseline remains the rolling 48-hour/7-day baseline.
- The episode recovery baseline is fixed at the first AKI-positive
  measurement.
- The episode ends at the first observed creatinine below both recovery
  thresholds under the primary first-recovery rule.
- A duration of exactly 48 hours is transient; persistence requires a duration
  strictly greater than 48 hours.
- The primary denominator remains first AKI episodes with ICU coverage through
  the first AKI-positive measurement plus 48 hours. Inclusion never requires a
  subsequent creatinine measurement.

## Additional continuity-support sensitivity

Only episodes classified as definite persistent in the v4 primary analysis
are reconsidered. For each such episode, the observed creatinine measurements
at or above the fixed episode recovery threshold are followed from the first
AKI-positive measurement through the measurement that establishes a minimum
possible duration greater than 48 hours.

Two maximum-adjacent-gap rules are evaluated separately:

1. no adjacent observed AKI-state creatinine measurements more than 24 hours
   apart;
2. no adjacent observed AKI-state creatinine measurements more than 36 hours
   apart.

If the chain contains a gap longer than the selected limit, the episode is
reclassified as `continuity_gap_indeterminate` for that sensitivity analysis.
All definite transient, interval-indeterminate and right-censored unresolved
episodes retain their v4 classifications. Thus, this sensitivity analysis can
lower the persistent-prevalence lower bound and widen the identified set, but
cannot change its upper bound.

This rule is deliberately conservative. It does not assert that recovery
occurred during a long gap; it asks how much of the definite-persistent lower
bound remains supported when long unobserved intervals are not bridged by the
single-episode continuity convention.

## Outputs

For each database and each gap limit, report:

- the unchanged primary denominator;
- original and continuity-supported definite-persistent counts;
- the number and proportion reclassified for a continuity gap;
- total classification-indeterminate proportion after reclassification;
- the revised persistent-prevalence lower bound, unchanged upper bound and
  revised identified-set width;
- the distribution of the maximum adjacent gap among originally definite
  persistent episodes.

No mortality, hospital ranking or other outcome association is examined.
