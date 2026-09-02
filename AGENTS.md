# Project: GBR coral bleaching mortality

## What this is

This project estimates reef-level coral mortality associated with bleaching events on the Great Barrier Reef. R is used for statistical modelling and reporting; Python is used mainly for environmental-data extraction and validation.

## Start here

- Read `README.md` for the human overview.
- Read `docs/pipeline.md` before running or changing the production pipeline.
- Read `docs/methods.md` before changing outcomes, predictors, validation, or likelihoods.
- Treat `config/model_registry.yml` as the source of truth for the selected model.
- Treat `config/pipeline.toml` as the source of truth for execution order.

## Conventions

- Follow `docs/conventions.md`.
- Run commands from the repository root; scripts use root-relative paths.
- Put reusable code in `src/`, executable Quarto analyses in `analysis/`, and tests in `tests/`.
- Keep generated data in `data/`, fitted objects and tables in `output/`, and rendered reports in `reports/`.
- Do not silently change the selected model, response definition, event window, validation folds, or registry metrics.
- Preserve the May cutoff: mortality is assessed at least two to three months after bleaching, and January to April surveys map to the preceding event year.
- Preserve the bleaching-mortality restriction. Known COTS and cyclone losses train separate cause layers; do not let them flatten the thermal response.
- Do not use event identity as an operational predictor.

## Current model contract

- The selected initial-forecast model is the registry entry `operational_rrn_raw_plus_manta_state`.
- It is a cause-aware INLA composite with Bernoulli occurrence and beta positive-magnitude components.
- The independent direct BRT is retained for within-event interpolation and nonlinear diagnostics, not as the initial forecast.
- BRMS/BRT formal models remain reproducible comparators; they are not the selected registry model.
- Leave-one-event-out performance is the primary selection evidence. Reef-blocked validation, severe-event RMSE, false-extreme rate, predictive R-squared, WAIC and DIC are supporting evidence.

## How to check work

- List the canonical pipeline: `python src/pipeline/run_pipeline.py --list`
- Check required inputs: `python src/pipeline/run_pipeline.py --profile current --check-inputs`
- Preview without fitting: `python src/pipeline/run_pipeline.py --profile current --dry-run`
- Run Python tests: `python -m unittest discover -s tests -v`
- Parse all R sources: `Rscript src/pipeline/check_r_syntax.R`
- Render the current report set: `quarto render`

A change is complete only when relevant tests pass and the registry, pipeline manifest, and documentation agree.

## Please do not

- Do not commit `data/`, `output/`, `reports/`, `.ai/`, caches, or local environments.
- Do not commit scratch plans from `.ai/`; promote reviewed conclusions to `docs/`.
- Do not edit archived notebooks as if they were production code.
- Do not overwrite fitted production artefacts without retaining provenance and validation outputs.
- Do not add a predictor solely because it improves an in-sample criterion.
- Do not rename registry model IDs or public helper functions without discussing the migration.
