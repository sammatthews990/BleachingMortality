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
7. `within-event-update`: build the locked 2025 initial assessment, harmonise audited six-event aerial surveys and validate the optional update separately from the initial forecast.
8. `machine-learning`: refit the residual and standalone BRT comparisons without promoting them automatically.
9. `diagnostics`: rebuild the current model ledger, major-miss register and interpretation figures.
10. `residual-gap-audit`: locate selected-model held-out error by event, region, reef, attributed cause and observation scale without fitting another model.
11. `experiment-ledger`: verify that every analysis has a concise finding and evidence trail.
12. `current-reports`: render the exploration index, registry-driven current reports and standalone diagnostics.

Run it with `python src/pipeline/run_pipeline.py --profile current`.

The selected-model stage is computationally expensive and reuses explicitly named cached fit files when the underlying scripts allow it. Never infer success from an existing output timestamp; inspect the stage log, validation tables and registry snapshot.

The residual audit can be refreshed without refitting the model after canonical
selected-model predictions are present:

```powershell
python src/pipeline/run_pipeline.py --profile model-gap-audit
quarto render analysis/investigations/remaining_residual_gap_audit.qmd
```

It uses only leave-one-event-out predictions from the registry-selected model.
Its squared-error shares describe where predictive error is concentrated; they
are diagnostics and cannot promote a new predictor or change the selected
model.

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

The IMOS thermal investigation is also outside the production profile. It
reconstructs an NOAA-equivalent DHW from bias-corrected 0.02-degree nighttime
SST and a SSTAARS-derived MMM before screening coral-threshold duration/peak,
standard Hobday MHW, and paired day/night terms. It also validates the cooling
contrast against the available hourly AIMS logger subset and constructs a
fixed-baseline NOAA CoralTemp empirical IDF sensitivity from 1986-2012 annual
block maxima:

```powershell
python -m pip install -r requirements/imos-thermal.txt
python src/pipeline/run_pipeline.py --profile imos-thermal-screen
```

The extractor caches SSTAARS harmonic and monthly percentile subsets plus
annual reef-event outputs. Use `--force` only when intentionally rebuilding
those generated parts. The screen is an investigation and cannot alter the
selected model without the formal promotion workflow.

## Within-event aerial update

The optional stage requires explicit 2025 ERA5 and PATMOS-x snapshots. Remote
acquisition remains outside the default profile; create them when absent:

```powershell
python src/data/extract_era5_weather.py --years 2025 --output data/processed/era5_weather_reef_year_2025.csv --part-dir data/processed/era5_weather_parts_2025
python src/data/fetch_patmosx_cloud.py --years 2025 --output data/processed/patmosx_cloud_reef_year_2025.csv --metadata data/processed/patmosx_cloud_metadata_2025.json
```

Then build the locked initial assessment and the distinct optional update:

```powershell
python src/pipeline/run_pipeline.py --stage within-event-update
```

The stage normally takes several minutes because the candidate coefficients are
refit across event-held-out and two spatial exclusion designs. The aerial
crosswalk, status-code exclusions and promotion boundary are documented in
`docs/aerial-early-bleaching-input-contract.md`; the retrospective locked 2025
boundary is in `docs/2025-initial-forecast-contract.md`. This stage does not
alter the canonical environmental table, model terms or selected-model registry.

## Standalone DHW occurrence diagnostic

The `dhw-diagnostic` profile compares raw NOAA and local-first
logger-adjusted DHW using the same simple Bernoulli mortality-occurrence model.
It excludes rows explicitly attributed to COTS, cyclone or flood, and reports
apparent, event-held-out, sector-held-out, pooled reef-blocked and
programme-held-out performance. It also fits separate descriptive curves for
all six events.

```powershell
python src/pipeline/run_pipeline.py --profile dhw-diagnostic --check-inputs
python src/pipeline/run_pipeline.py --profile dhw-diagnostic
quarto render analysis/investigations/bleaching_only_dhw_binomial_assessment.qmd
```

This profile is diagnostic only and cannot update the selected model registry.
The 2025 rows have no validated local-first correction, so zero adjustment is
recorded as unavailable rather than interpreted as evidence of zero NOAA bias.

## Outputs and provenance

`data/`, `output/` and `reports/` are local generated stores and are ignored by Git. A clean clone therefore needs the source data bundle or must rerun acquisition. The registry records the selected model identity and expected fit path, while model scripts write fold predictions, metrics, tuning, convergence and fit objects beneath `output/`.

When a production result changes, update together:

- `config/model_registry.yml` and, if needed, `config/model_terms.yml`;
- the relevant current Quarto analysis;
- `docs/model-status.md` and the development history;
- tests or data-contract checks;
- the pipeline manifest if execution order or prerequisites changed.
- `config/experiments.yml` with the experiment question, decision, concise finding and primary evidence.
