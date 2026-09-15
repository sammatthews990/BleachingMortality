# Token-efficient experiment workflow

The aim is to preserve scientific traceability without requiring each person or AI session to reread dozens of notebooks and large generated tables.

## Context loading order

1. Read `AGENTS.md` and the relevant section of `README.md`.
2. Read `config/experiments.yml` and `config/model_registry.yml` for prior decisions.
3. Open only the report and primary evidence linked by the relevant registry record.
4. Inspect source code, detailed outputs or raw data only when the question cannot be answered at the preceding layer.

Use targeted filename and symbol searches. Avoid broad dumps of `data/`, `output/`, rendered HTML, archived notebooks or model objects.

## Closing an experiment

An experiment is complete only when its code and detailed report are reproducible and its registry record is updated. Keep that record deliberately small:

- `question`: the decision the experiment addresses;
- `status`: whether it changed the production path;
- `finding`: one decision-ready sentence, no more than 45 words;
- `report`: the executable Quarto source;
- `evidence`: one or a few primary generated artefacts, not every output;
- `updated`: the date on which the summary was checked.

Run `Rscript src/evaluation/validate_experiment_registry.R` after adding or renaming an analysis. The validator checks coverage, uniqueness, required fields, length limits and report paths. Evidence files are reported as locally available or absent because generated outputs are intentionally not committed.

`analysis/exploration.qmd` renders the registry into a current synthesis and reproducible per-experiment blocks. It does not fit models or ingest large result tables, so it is safe to render frequently.

## Efficient conversations and handovers

- Refer to experiment IDs and file paths instead of restating whole reports.
- Share the smallest relevant error excerpt, not full console logs.
- State what changed since the previous result; do not repeat unchanged background.
- Start a new conversation at a major scientific boundary, with a short handover containing the active question, experiment IDs, changed paths and last checks.
- Put temporary reasoning in `.ai/`; put lasting scientific conclusions in the experiment registry or canonical documentation.
- Run the narrowest relevant validation first. Expensive refits should follow a dry run and occur only after the specification is frozen.

The registry is a navigation and decision layer, not a substitute for detailed methods, diagnostics or model artefacts.
