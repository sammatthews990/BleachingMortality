# Development history and line of sight

## Original analyses

The project began in two large notebooks:

- `archive/notebooks/01_DHW_Mortality_Exploratory.qmd` assembled and explored mortality and environmental predictors.
- `archive/notebooks/02_DHW_Mortality_Modelling.qmd` developed the earlier regression, BRMS and BRT candidates.

They are retained as historical evidence, not production entry points. Their extracted, maintained replacements now live in `src/data/`, `src/models/` and `src/evaluation/`.

## Pipeline maturation

The detailed chronological working record is preserved at `archive/plans/research_pipeline_implementation_plan.md`. Major steps were:

1. recreate environmental inputs, including PATMOS-x cloud rather than reusing historical values;
2. harmonise mortality outcomes with programme-specific survey timing and the restricted bleaching response;
3. establish fixed blocked validation and retain 2024 in cross-validation and production fitting;
4. restore preceding-year Acropora and leakage-safe prior thermal exposure;
5. recover MMP Secchi using IMOS and assess salinity, rainfall, wind and freshwater proxies;
6. compare beta, binomial, BRMS and BRT frameworks and add formal outlier assessment;
7. develop joint and cause-aware models to prevent known COTS/cyclone losses flattening the thermal relationship;
8. add spatial INLA effects, local-first DHW correction, prospective COTS state and independent BRT comparisons;
9. consolidate the selected model in the registry and numbered current reports.

## Current analysis chain

The five files in `analysis/current/` form the current evidence chain:

- `03_Best_Fit_Model.qmd`: persistent registry-driven model ledger and diagnostics;
- `04_COTS_ENSO_reef_control_assessment.qmd`: COTS intensity, ENSO screen and named controls;
- `05_spatial_early_bleaching_update_assessment.qmd`: leakage-controlled early-event update;
- `06_prospective_cots_and_inla_brt_assessment.qmd`: prospective COTS and residual-BRT decision;
- `07_standalone_brt_envelope_ensemble_assessment.qmd`: independent BRT and envelope comparison.

Other executable reports are retained under `analysis/investigations/`. They document useful negative results and sensitivities but do not define the current selected model.

Git history remains the authoritative record of files removed from the active tree. The cleanup did not alter fitted model specifications or registry results.
