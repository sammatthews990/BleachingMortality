# Project conventions

## Paths and entry points

Run every command from the repository root. Production scripts use root-relative paths and should locate the root through a stable marker such as `config/model_registry.yml` when they may be sourced by Quarto.

`config/pipeline.toml` owns command order. If a production filename or prerequisite changes, update the manifest and `docs/pipeline.md` in the same change.

## Code placement

- `src/data/`: downloading, extraction, feature engineering and outcome harmonisation.
- `src/models/`: model fitting and ensemble construction.
- `src/evaluation/`: diagnostics, model comparisons, audits and sensitivity analyses.
- `src/lib/`: shared helpers without a top-level analysis side effect.
- `src/pipeline/`: orchestration and repository-wide checks.

Readable sequential code is preferred when an analysis is naturally linear. Create a function when logic is reused, independently testable, or clearer under a name; do not wrap every report step in a function.

## Names and outputs

- Use `snake_case` for files, variables and tabular fields, except retained programme names and established registry IDs.
- Use proportions in `[0, 1]` internally for cover, occurrence and mortality. State units at ingestion boundaries.
- Use `event_year` for bleaching events and retain the programme-specific source observation ID.
- Write generated tables, figures and fitted objects only beneath `output/`.
- Write rendered Quarto output only beneath `reports/`.
- Never make a generated file the only record of a scientific decision; decisions belong in `config/` or `docs/`.

## Reproducibility

- Set and record random seeds for stochastic fitting or resampling.
- Preprocessing must be learned inside each training fold and applied unchanged to its assessment fold.
- Save fold assignments, tuning results, metrics, predictions and convergence diagnostics with fitted models.
- Keep acquisition provenance, spatial/temporal coverage, transformations and missingness explicit.
- Do not replace missing measured data with a proxy without marking the value and its source.

## Model changes

New candidates do not become production models by overwriting outputs. Fit them to a new output directory, compare them on the fixed validation folds, document the decision, then update the registry deliberately.

Model comparison should report at least RMSE, MAE, predictive R-squared, occurrence Brier score where applicable, severe-event RMSE, false-extreme rate and event-level results. WAIC, DIC or AIC are supplementary because prediction of unseen events is the primary objective.

## Tests

Python extraction or transformation logic needs a focused unit test. R changes must pass the syntax check and should add a small deterministic check when logic can be isolated without fitting a full model. Data-contract changes require updates to both validation code and documentation.
