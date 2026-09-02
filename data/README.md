# Local data store

Data are not committed. Restore the project data bundle here or rerun the acquisition scripts in `src/data/`.

- `raw/`: immutable downloaded/source files.
- `cache/`: restartable remote-data cache.
- `curated/`: manually reviewed inputs with provenance.
- `processed/`: model-ready tables, manifests and QA summaries.

See `docs/pipeline.md` and `config/pipeline.toml` for required production inputs.
