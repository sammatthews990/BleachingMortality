# Aerial early-bleaching input contract

## Purpose and model boundary

The aerial survey product is an explicitly timed **within-event update** input. It is not available to the selected initial forecast, is not part of the canonical environmental table, and does not change `config/model_terms.yml` or the selected model registry. Any future promotion requires matched event-held-out and spatially independent validation.

## Accepted sources

`src/data/build_aerial_early_bleaching.R` reads the historical Hughes reef-level files for 2016, 2017 and 2020 and the reviewed AIMS-GBRMPA files beneath `data/AerialBleach2022_25/` for 2022, 2024 and 2025. Raw files are immutable inputs. The builder records source paths and MD5 hashes in the generated audit.

The contemporary observations follow AIMS SOP 11. Valid bleaching categories describe the percentage of visible shallow living coral cover bleached:

| Category | Bleached cover |
|---:|---:|
| 0 | less than 1% |
| 1 | 1-10% |
| 2 | 11-30% |
| 3 | 31-60% |
| 4 | 61-90% |
| 5 | more than 90% |

Codes above 5 and their accompanying no-view, sediment, depth, or low/no-live-coral descriptions are observation-status codes. They must never be converted to severity. Historical `score` values are mapped through the reviewed `config/aerial_bleaching_score_crosswalk.csv`; the unexpected 2020 value 5 is flagged and top-coded to the documented historical greater-than-60% class rather than silently treated as a contemporary category.

## Harmonised fields and timing

The common continuous field is `aerial_bleached_cover_midpoint` (also exposed as `aerial_severity_01`). It is the midpoint of the reviewed historical or SOP11 percentage-cover bin. Contemporary non-integer mean categories are converted by piecewise-linear interpolation between SOP11 bin midpoints. The binary sensitivity `aerial_high_bleaching` retains the supplied historical indicator or a contemporary rounded category of at least 3. Native scores, scale, mapping description and review flags remain in the outputs.

For sources with dates, the builder also writes March 31 and April 30 availability fields and cutoff-specific values. Survey dates are absent from the supplied 2016, 2017 and 2020 files, so their timing remains explicitly unavailable. The strictly dated validation is therefore restricted to 2022, 2024 and 2025. The common event-time midpoint is tested across all six mortality-labelled events, while its historical timing limitation remains explicit.

Every output carries or inherits these operational restrictions:

- `product_mode = within_event_update`
- `use_in_initial_forecast = FALSE`
- `use_in_selected_model = FALSE`

## Generated outputs

Running `Rscript src/data/build_aerial_early_bleaching.R` writes:

- `data/processed/aerial_early_bleaching_observation.csv`: one row per native survey unit after removal of duplicate GIS-join rows;
- `data/processed/aerial_early_bleaching_reef_event.csv`: one audited row per `(ReefID, event_year)` for validation joins;
- `data/processed/aerial_early_bleaching_source_audit.csv`: row counts, usable/status counts, dates, cutoff availability, hashes and provenance;
- `data/processed/aerial_early_bleaching_conflict_audit.csv`: repeated reef identifiers, coordinate spans and score conflicts requiring review;
- `data/processed/aerial_early_bleaching_provenance.csv`: contract version, category definitions and method reference.

The builder fails if required source columns are missing, native-unit or reef-event keys remain duplicated, or a status code is retained as severity.

## Validation and promotion rule

`src/evaluation/test_spatial_early_bleaching_update.R` joins this table only to the separate early-event validation grid. It excludes every target sector or latitude block before constructing the spatial signal and holds out the target event when estimating update coefficients. The selected conditional mortality magnitude stays fixed; the tested update changes occurrence only.

The most constrained candidate standardises common aerial severity and April 30 RHIS burden within each event. Each signal is weighted by `effective_n / (effective_n + 3) * exp(-nearest_km / 450)`, single-source support is attenuated, and the reliability-weighted event mean is removed. One non-negative occurrence coefficient is estimated with the initial logit as an offset, no intercept, ridge precision 4, and no change to conditional magnitude. This prevents an event-wide bleaching level from becoming an event-identity proxy.

Promotion review requires a positive residual direction in at least four of all six events, lower aggregate RMSE, MAE and occurrence Brier score under both spatial designs, no material severe-RMSE or false-extreme penalty, gains across events, and lower RMSE and occurrence Brier score in the locked 2025 assessment under both designs. The constrained consensus clears the direction guard but has negligible skill and fails the aggregate and locked-2025 both-design guards. Until all guards pass, it remains off by default and may be reported only as a separate observed-condition layer. The locked initial-forecast boundary is documented in `docs/2025-initial-forecast-contract.md`.

## Method references

- [AIMS, *Standard Operating Procedure 11: Aerial Surveys of Coral Bleaching* (version 3, 2022)](https://www.aims.gov.au/sites/default/files/2022-06/AIMS_SOP11v3_Aerial-Surveys-Coral-Bleaching_202206.pdf), DOI `10.25845/n00q-z603`.
- [AIMS 2024 GBR aerial bleaching survey metadata and final report](https://apps.aims.gov.au/metadata/view/3fd37e52-12cf-4b59-8873-96b82bc955c9).
