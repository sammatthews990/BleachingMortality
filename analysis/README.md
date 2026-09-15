# Executable analyses

`exploration.qmd` is the lightweight, reproducible findings index. It reads the compact experiment registry without rerunning models or loading large outputs. Read it before opening individual analyses.

`current/` contains the small report set that supports the registry-selected model. `investigations/` retains earlier model comparisons, environmental audits and sensitivity analyses. Every report in both folders must have a record in `config/experiments.yml`.

Quarto renders the current set and the standalone DHW diagnostic through the root `_quarto.yml`. Generated HTML belongs in `reports/` and is not versioned.
