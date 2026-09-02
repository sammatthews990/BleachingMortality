# Source code map

- `data/`: environmental acquisition, feature engineering, monitoring-data extraction and mortality harmonisation.
- `models/`: maintained fitting frameworks, including selected INLA components and retained BRMS/BRT candidates.
- `evaluation/`: validation summaries, diagnostics, outlier audits and controlled sensitivity tests.
- `lib/`: shared helpers, model registry access and diagnostic utilities.
- `pipeline/`: the readable manifest runner and repository-wide checks.

The production order is not alphabetical. Follow `config/pipeline.toml` or run `python src/pipeline/run_pipeline.py --list`.

Files beginning with `test_` in `evaluation/` are scientific candidate screens that fit models; they are not unit tests. Fast unit tests live in `tests/`.
