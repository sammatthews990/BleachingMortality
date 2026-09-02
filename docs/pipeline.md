# Reproducible pipeline

`config/pipeline.toml` is the executable manifest. `src/pipeline/run_pipeline.py` prints and runs it without hiding the underlying commands.

## Before running

1. Install R, Quarto and Python 3.11 or newer.
2. Restore the project data bundle beneath `data/` using the paths in the pipeline manifest.
3. Install focused Python dependencies from `requirements/` when rerunning environmental acquisition.
4. Install the R packages used by the selected scripts, including INLA, `dplyr`, `readr`, `sf`, `Matrix`, `yaml`, `jsonlite`, `ggplot2`, `gbm`, `readxl`, `stringr` and `tidyr`.
5. Run the input check before an expensive fit.

```powershell
python src/pipeline/run_pipeline.py --profile current --check-inputs
python src/pipeline/run_pipeline.py --profile current --dry-run
```

## Current production profile

The `current` profile runs these stages in order:

1. `environment-check`: validate the assembled environmental table and its contract.
2. `outcomes`: rebuild the programme-specific restricted mortality outcomes.
3. `validation`: recreate fixed event-, reef- and region-blocked folds.
4. `pressure-context`: rebuild RRN pressure and the cause-labelled annual transition table.
5. `screening-models`: rebuild cached spatial INLA candidates used by the selected-model evidence and residual map.
6. `selected-model`: fit the cause-aware thermal/cyclone/COTS chain and select the prospective COTS state.
7. `machine-learning`: refit the residual and standalone BRT comparisons without promoting them automatically.
8. `diagnostics`: rebuild the current model ledger, major-miss register and interpretation figures.
9. `current-reports`: render the registry-driven current reports.

Run it with `python src/pipeline/run_pipeline.py --profile current`.

The selected-model stage is computationally expensive and reuses explicitly named cached fit files when the underlying scripts allow it. Never infer success from an existing output timestamp; inspect the stage log, validation tables and registry snapshot.

## Formal BRMS/BRT profile

The retained methodological comparison is separate:

```powershell
python src/pipeline/run_pipeline.py --profile formal-comparators
```

This rebuilds the formal BRT, inflated-beta BRMS, binomial sensitivity, summaries, ensemble, prediction uncertainty and outlier assessment. It must not update the production registry without a new matched comparison against the selected INLA composite.

## Environmental acquisition

Large remote-data steps are deliberately not part of the default current profile. They can be run directly:

```powershell
python -m pip install -r requirements/patmosx.txt
python src/data/fetch_patmosx_cloud.py
python src/data/fetch_dms_environmental_data.py
python src/data/validate_environmental_data.py
```

ERA5 weather uses `requirements/era5.txt`. Environmental scripts cache downloads beneath `data/cache/` and write processed tables and provenance beneath `data/processed/`.

## Outputs and provenance

`data/`, `output/` and `reports/` are local generated stores and are ignored by Git. A clean clone therefore needs the source data bundle or must rerun acquisition. The registry records the selected model identity and expected fit path, while model scripts write fold predictions, metrics, tuning, convergence and fit objects beneath `output/`.

When a production result changes, update together:

- `config/model_registry.yml` and, if needed, `config/model_terms.yml`;
- the relevant current Quarto analysis;
- `docs/model-status.md` and the development history;
- tests or data-contract checks;
- the pipeline manifest if execution order or prerequisites changed.
