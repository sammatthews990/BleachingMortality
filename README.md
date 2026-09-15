# GBR coral bleaching mortality

This repository builds reef-level predictions of coral mortality associated with bleaching events on the Great Barrier Reef. It combines monitoring outcomes from LTMP, manta tow and MMP with thermal exposure, prior heat history, coral composition, water quality, cloud, currents, freshwater, cyclone and COTS information.

The repository now has one explicit production path. Historical notebooks and model screens are retained, but they are not mixed into the current execution path.

## Current result

The selected initial-forecast model is declared in `config/model_registry.yml`: `operational_rrn_raw_plus_manta_state`, a cause-aware INLA composite. Its primary leave-one-event-out metrics are RMSE 0.1485, predictive R-squared 0.2998, severe-event RMSE 0.3267 and false-extreme rate 0.0458.

The direct BRT remains useful for interpolation after an event state is represented, and the formal BRMS/BRT models remain available as methodological comparators. Neither replaces the selected initial forecast. See `docs/model-status.md` for the full decision.

## Repository map

```text
.
|-- AGENTS.md                 standing instructions for people and AI tools
|-- README.md                 project overview and quick start
|-- analysis/                 executable Quarto analyses
|   |-- exploration.qmd       fast registry-driven findings overview
|   |-- current/              reports supporting the selected model
|   `-- investigations/       retained scientific screens and sensitivities
|-- archive/                  superseded notebooks, plans and legacy assets
|-- config/                   model registry, terms and pipeline order
|-- data/                     local inputs and derived data; not committed
|-- docs/                     canonical methods and development documentation
|-- manuscript/               manuscript source
|-- output/                   generated models, predictions and figures
|-- reports/                  rendered Quarto output
|-- requirements/             focused Python dependency lists
|-- src/                      production code and retained model candidates
|   |-- data/                 acquisition, extraction and harmonisation
|   |-- evaluation/           diagnostics, comparisons and sensitivity tests
|   |-- lib/                  shared R helpers and registry utilities
|   |-- models/               fitting code for INLA, BRMS, BRT and ensembles
|   `-- pipeline/             orchestration and structural checks
`-- tests/                    Python unit and data-contract tests
```

## Quick start

Run all commands from the repository root.

```powershell
python src/pipeline/run_pipeline.py --list
python src/pipeline/run_pipeline.py --profile current --check-inputs
python src/pipeline/run_pipeline.py --profile current --dry-run
```

The dry run prints the exact commands used to rebuild the current selected model, its machine-learning comparators, diagnostics and current reports. Remove `--dry-run` only after the input check passes; the INLA fits are computationally expensive.

To run one named stage:

```powershell
python src/pipeline/run_pipeline.py --stage selected-model
```

To rebuild the locked 2025 initial-forecast assessment and validate the separate six-event aerial/RHIS within-event update:

```powershell
python src/pipeline/run_pipeline.py --stage within-event-update
```

To rerun the formal BRMS/BRT comparison separately:

```powershell
python src/pipeline/run_pipeline.py --profile formal-comparators
```

To compare raw and local-first logger-adjusted DHW in the standalone
bleaching-compatible binomial occurrence diagnostic:

```powershell
python src/pipeline/run_pipeline.py --profile dhw-diagnostic
```

## Tests

```powershell
python -m unittest discover -s tests -v
Rscript src/pipeline/check_r_syntax.R
```

Start with `quarto render analysis/exploration.qmd` for the fast experiment synthesis. Render the exploration page, five current model reports and the standalone DHW diagnostic with `quarto render`. Rendered files are written under `reports/` and are not versioned.

## Documentation

- `docs/pipeline.md` gives the ordered production workflow and data boundaries.
- `docs/methods.md` records the response, event timing, predictors, model structure and validation logic.
- `docs/aerial-early-bleaching-input-contract.md` defines the audited within-event aerial input and its promotion boundary.
- `analysis/investigations/bleaching_only_dhw_binomial_assessment.qmd` compares raw and logger-adjusted DHW without changing production selection.
- `docs/model-status.md` distinguishes the selected model from retained comparators and rejected candidates.
- `docs/2025-initial-forecast-contract.md` defines the retrospective locked-input 2025 assessment and its separation from the nowcast.
- `docs/validation.md` contains the blocked validation protocol.
- `docs/history.md` maps the original exploratory notebooks to the current pipeline.
- `docs/conventions.md` contains coding, naming and provenance rules.
- `docs/token-efficient-workflow.md` explains the summary-first experiment and conversation workflow.

Transient plans and handover notes belong in `.ai/` and are intentionally ignored. Reviewed decisions belong in `docs/` or `config/`.
