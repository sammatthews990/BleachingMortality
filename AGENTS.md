# Project: GBR coral bleaching mortality

## What this is

This project estimates reef-level coral mortality associated with bleaching events on the Great Barrier Reef. R is used for statistical modelling and reporting; Python is used mainly for environmental-data extraction and validation.

## Start here

- Read `README.md` for the human overview.
- Read `config/experiments.yml` before searching analyses or outputs; it is the compact record of questions, findings and evidence.
- Read `docs/pipeline.md` before running or changing the production pipeline.
- Read `docs/methods.md` before changing outcomes, predictors, validation, or likelihoods.
- Treat `config/model_registry.yml` as the source of truth for the selected model.
- Treat `config/pipeline.toml` as the source of truth for execution order.

## Token-efficient working rules

- Load context in layers: this file, the experiment/model registries, then only the linked report or source needed for the task.
- Do not recursively read `data/`, `output/`, `reports/`, `archive/` or all investigation notebooks. Search filenames or symbols first with `rg`.
- Treat `analysis/exploration.qmd` as the fast human overview. It renders from `config/experiments.yml` without fitting models or reading large result tables.
- Every completed experiment must update one registry record in the same change: one short question, one decision status, one finding of at most 45 words, and links to its report and primary evidence.
- Prefer targeted checks while developing. Run full suites or model fits only when the change can affect them.
- In handovers and conversation, report decisions, changed paths and failed checks; link to logs or reports instead of copying them.
- Keep `.ai/` notes limited to the active task. Promote durable findings to the registry or `docs/`, then remove stale scratch context.

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
- Validate the experiment ledger: `Rscript src/evaluation/validate_experiment_registry.R`
- Render the lightweight findings report: `quarto render analysis/exploration.qmd`
- Render the current report set: `quarto render`

A change is complete only when relevant tests pass and the registry, pipeline manifest, and documentation agree.

## Please do not

- Do not commit `data/`, `output/`, `reports/`, `.ai/`, caches, or local environments.
- Do not commit scratch plans from `.ai/`; promote reviewed conclusions to `docs/`.
- Do not edit archived notebooks as if they were production code.
- Do not overwrite fitted production artefacts without retaining provenance and validation outputs.
- Do not add a predictor solely because it improves an in-sample criterion.
- Do not rename registry model IDs or public helper functions without discussing the migration.
