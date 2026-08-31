# Version 1.2.1 — portable reproducibility release

## Summary

Version 1.2.1 is a release-quality reproducibility update for the
measurement-aware 48-hour AKI-persistence phenotype. It preserves the frozen
v1.2 numerical results and makes the public implementation and release
metadata match the documented workflow.

## Changes

- Uses environment-configured controlled-data locations rather than author
  machine paths, with clear failure messages when a required path is absent.
- Restores importability of every documented analysis entry point.
- Adds continuous-integration gates for compilation, every command-line
  `--help` path, synthetic tests and the aggregate-output reconciliation audit.
- Pins the reference Python execution environment and documents data-free
  clean-room verification.
- Updates the public archive boundary, source-output provenance, two-stage
  eICU hospital-to-unique-patient resampling description, citation metadata and
  Zenodo metadata.
- Removes legacy aggregate hospital-ranking diagnostics from the distributed
  controlled-thinning results; the retained manuscript does not use or
  interpret hospital ranking.
- Keeps only source code, documentation and non-disclosive aggregate outputs;
  restricted source data and trajectory caches remain excluded.

## Scientific scope preserved

The v1.2 corrections remain in force: index-AKI search is separate from
post-index follow-up; SICdb coverage is re-anchored to first ICU bed assignment;
and controlled thinning is anchored to the original index episode, so a later
recurrent AKI cannot replace a missed index episode. No primary numerical
result is changed by this release.

## Citation

The version-specific Zenodo DOI must be added only after Zenodo publishes the
archive created from the exact `v1.2.1` tag. Do not use the DOI of an earlier
release for this version.
