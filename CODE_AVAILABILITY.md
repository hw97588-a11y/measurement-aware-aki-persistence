# Code availability

This repository provides the common phenotype engine, source mappings,
statistical scripts, synthetic boundary tests, controlled-thinning code,
analysis specifications and non-disclosive aggregate outputs supporting the
manuscript.

The repository intentionally excludes all source patient-level data,
episode-level analytic files, hospital-specific model inputs, phase-level
simulation files and protected trajectory caches. Qualified researchers must
obtain MIMIC-IV v3.1, eICU-CRD v2.0 and SICdb v1.0.8 through their respective
PhysioNet access routes and reproduce the analysis locally under the applicable
data-use agreements.

The active version-control record is
https://github.com/hw97588-a11y/measurement-aware-aki-persistence. The
reproducibility release associated with the current submission is version
1.2.1. Its version-specific Zenodo DOI must be inserted only after Zenodo has
published an archive created from the exact `v1.2.1` tag; the DOI of an earlier
release must not be reused.

## Manuscript-ready wording after DOI publication

> **Code availability.** Analysis code, the phenotype specification, synthetic
> tests and non-disclosive aggregate outputs are available in the versioned
> Zenodo archive (version 1.2.1, DOI: **[insert version-specific DOI]**) and at
> https://github.com/hw97588-a11y/measurement-aware-aki-persistence. The
> controlled MIMIC-IV, eICU-CRD and SICdb source data are not redistributed and
> must be accessed through their respective authorization procedures.

Replace the bracketed DOI after confirming that it resolves to the exact public
archive. Do not claim that the controlled source data themselves are publicly
available through Zenodo or GitHub.
