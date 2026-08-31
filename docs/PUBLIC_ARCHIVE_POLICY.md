# Public archive policy

## What may be released

The public GitHub and Zenodo archives may contain only:

- source code and tests;
- source-mapping and phenotype documentation;
- the MIT licence, citation metadata and release notes; and
- non-disclosive source-level aggregate outputs in `results/`.

The aggregate outputs contain no direct patient, admission, ICU-stay,
episode, hospital or phase identifiers. They support the released tables and
figures as documented in `RESULT_PROVENANCE_V1.2.md`.

## What must never be released

Do not add, commit, attach or upload:

- raw MIMIC-IV, eICU-CRD or SICdb files, extracts or source archives;
- any patient-, admission-, encounter-, spell-, unit- or hospital-level table;
- trajectory caches, phase-level thinning files, temporary source-derived
  intermediate files or local result workspaces;
- database credentials, API keys, data-access records, login cookies or local
  environment files; or
- manually copied line-level examples from controlled data.

`.gitignore` is a safeguard, not a release decision. A maintainer must inspect
the staged file list and the generated release archive every time.

## Required pre-publication audit

From the exact release commit, run:

```bash
git status --short
git ls-files
git check-ignore -v .env data/example.csv protected/example.parquet || true
git archive --format=tar --prefix=measurement-aware-aki-persistence-v1.2.1/ HEAD \
  | tar -tf -
```

The Git working tree must be clean. Inspect the archive list for excluded
filenames and confirm that it contains no archives or raw-data extensions such
as `.zip`, `.rar`, `.gz`, `.parquet`, `.feather`, `.pkl`, `.sqlite`, `.dta`,
`.sas7bdat`, `.rds`, `.RData`, `.xlsx` or `.xls`. The only intended `.csv` and
`.json` files are the documented aggregate outputs under `results/` and
`.zenodo.json`.

For the final archive, record its SHA-256 digest and compare the Zenodo asset
after upload with the local file byte for byte or by matching checksum. See
`RELEASE_CHECKLIST.md` for the required release sequence.
