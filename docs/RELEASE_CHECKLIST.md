# Release checklist — version 1.2.1

This checklist is intentionally release-specific. Complete it from a clean
clone of the exact commit intended for tag `v1.2.1`; do not publish a GitHub
release or Zenodo record before every applicable item is complete.

## 1. Code and aggregate-output gates

- [ ] Confirm all runtime source paths are read from `MIMIC_IV_PATH`,
  `EICU_CRD_PATH`, `SICDB_PATH` and, when needed, `SICDB_MEMBER_ROOT`; no
  author-machine path is present in a published Python file.
- [ ] Run `python -m compileall -q .`.
- [ ] Run every `--help` command listed in `README.md`.
- [ ] Run `python -m unittest -v test_interval_aki_v4_engine.py test_ndt_continuity_gap_sensitivity.py`.
- [ ] Run `python audit_v6_targeted_outputs.py` and confirm `44/44` checks.
- [ ] Run `python prepare_public_aggregate_release.py --write`, then
  `python prepare_public_aggregate_release.py --check`, to remove legacy
  aggregate hospital-ranking diagnostics not used by the current manuscript.
- [ ] Confirm that CI for the exact release commit is green.
- [ ] Re-run the primary scripts in the governed environment or verify their
  stored aggregate outputs against the release commit under the documented
  data-use agreements.

## 2. Environment and documentation gates

- [ ] Confirm `requirements.txt` matches the actual clean-room execution
  environment described in `ENVIRONMENT_FREEZE_v1.2.1.md`.
- [ ] Confirm `README.md`, `CODE_AVAILABILITY.md`, `CITATION.cff`,
  `.zenodo.json`, `CHANGELOG.md` and `RELEASE_NOTES_v1.2.1.md` all state
  version `1.2.1`.
- [ ] Confirm all links target the canonical GitHub repository and do not cite
  an earlier release DOI as if it identified version 1.2.1.
- [ ] Confirm the output-to-code mapping in `RESULT_PROVENANCE_V1.2.md` and
  `AGGREGATE_OUTPUT_MANIFEST.md` matches every retained public result file.

## 3. Privacy and archive gates

- [ ] Complete `PUBLIC_ARCHIVE_POLICY.md`'s file-list audit.
- [ ] Inspect `git diff --cached --name-only` and `git ls-files` for governed
  data, caches, credentials and personal local paths.
- [ ] Build the source archive from the exact commit and inspect its contents.
- [ ] Record the archive SHA-256 digest in the GitHub release description and
  Zenodo upload notes.

## 4. Publish in the correct order

1. Commit the verified source and metadata; capture the commit SHA.
2. Create annotated Git tag `v1.2.1` at that SHA.
3. Push the commit and tag; wait for the GitHub Actions workflow to complete.
4. Create the GitHub release from that tag using `RELEASE_NOTES_v1.2.1.md`.
5. Create a **new version** from the previous Zenodo record or upload the exact
   tagged source archive. Use `.zenodo.json` as the metadata source and verify
   that the GitHub release URL, version and publication date are correct.
6. Publish Zenodo and record the new version-specific DOI, concept DOI, file
   name, file checksum and public URL.
7. Confirm that Zenodo's uploaded archive has the same checksum as the local
   release asset and that the DOI resolves publicly.
8. Insert the new version-specific DOI into the manuscript Data Availability
   Statement, cover letter where applicable, `CITATION.cff` and
   `CODE_AVAILABILITY.md`. If a metadata change is made after tag creation,
   create a new patch release rather than silently changing the source behind
   the archived tag.

## 5. Final record

Save a short `PUBLIC_ARCHIVE_VERIFICATION_v1.2.1.md` outside the public source
tree with the release commit SHA, tag, archive SHA-256, Zenodo DOI, concept DOI,
Zenodo file checksum, GitHub Actions run URL and date of verification.
