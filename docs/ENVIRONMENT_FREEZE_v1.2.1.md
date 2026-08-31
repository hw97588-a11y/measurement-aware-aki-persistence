# Frozen execution environment — version 1.2.1

## Reference interpreter and package versions

The final clean-room verification must use Python 3.13.12 and the exact pinned
packages in `requirements.txt`:

| Component | Reference version |
|---|---:|
| Python | 3.13.12 |
| NumPy | 2.4.3 |
| pandas | 2.3.3 |
| patsy | 1.0.2 |
| scikit-learn | 1.8.0 |
| SciPy | 1.17.1 |
| statsmodels | 0.14.6 |

The frozen analysis does not depend on scikit-learn for the retained primary
estimand, but it remains pinned because it was present in the execution
environment and may be used by supporting scripts.

## System utilities

The data-free tests require only Python and the packages above. Reproduction
against controlled source data additionally requires:

- `gzip` for streamed compressed eICU files;
- `rg` (ripgrep) for the eICU creatinine stream filter;
- `bsdtar` for streaming members of the SICdb archive; and
- a C++17 compiler only when running the accelerated controlled-thinning
  simulation.

Verify the data-free environment from a clean clone:

```bash
python --version
python -m pip check
python -m compileall -q .
python -m unittest -v test_interval_aki_v4_engine.py test_ndt_continuity_gap_sensitivity.py
python audit_v6_targeted_outputs.py
```

Then run every command listed under “Data-free verification” in `README.md`.
The command-line `--help` checks must pass before any restricted database path
is configured. The source data must remain outside the repository and outside
any public build or archive directory.
