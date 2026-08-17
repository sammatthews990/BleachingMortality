# GBR coral-bleaching mortality pipeline: ordered implementation plan

This plan records the improvements identified in the review of
`01_DHW_Mortality_Exploratory.qmd` and `02_DHW_Mortality_Modelling.qmd`.
Work is ordered so that predictor provenance and leakage are resolved before
models are refitted or interpreted.

## 1. Reconstruct event-year cloud cover — completed

- Source NOAA-18 AVHRR PATMOS-x v6.0 `cloud_fraction` from NOAA's public CDR
  archive, rather than importing the archived paper values.
- Compute the mean of all available ascending-pass observations from 1 January
  through 31 March for each event year. `cloudp_90` means the roughly 90-day
  JFM window; it is not a 90th percentile.
- Use the paper-matching raster convention (nearest latitude and first
  longitude cell centre at or east of the reef) for the primary predictor,
  while retaining ordinary nearest-neighbour and alternative pass aggregates
  as audit columns.
- Reconstruct 2016, 2017, 2020, 2022, and 2024 independently. Never fill a
  missing event year from another year's reef mean or a global median.
- Preserve file status and sampling coverage. The available paired days are
  70, 89, 89, 82, and 82 respectively; 2020 is predominantly preliminary and
  2024 is entirely preliminary in the current archive selection.
- Validate, but do not calibrate, against the archived Cheung et al. values.
  Across the authors' exact reef coordinates, the chosen reconstruction gives
  correlations of 0.9775, 0.9980, and 0.9979 and MAEs of 0.0074, 0.0037, and
  0.0061 for 2016, 2017, and 2020.

Implementation: `scripts/fetch_patmosx_cloud.py`,
`scripts/requirements-patmosx.txt`, and `tests/test_patmosx_cloud.py`.

## 2. Formalise data contracts and provenance

- Separate raw acquisition, reef registry, predictor construction, survey
  harmonisation, modelling tables, fitted models, and figures into explicit
  targets with immutable inputs and metadata manifests.
- Add schema, units, valid ranges, key uniqueness, missingness, spatial-match
  distance, and event-year coverage assertions at every hand-off.
- Pin R/Python dependencies and random seeds; render reports only from saved,
  versioned analysis products.

## 3. Harmonise outcomes and observation processes

- Define mortality estimands separately for manta tow, LTMP benthic surveys,
  and any other programmes; document denominators, timing, depth, taxonomic
  scope, and aggregation.
- Avoid treating heterogeneous proportions as interchangeable. Model survey
  programme and observation precision explicitly, with reef/event identifiers
  retained through every join.
- Audit zero/one handling and replace arbitrary nudging where an observation
  model can represent boundary values directly.

## 4. Build leakage-resistant validation

- Predefine spatially and temporally blocked outer resampling (held-out reefs,
  regions, and bleaching years) before model selection.
- Perform preprocessing, feature selection, tuning, and imputation inside each
  training fold. Keep the final event/year holdout untouched until the model
  specification is frozen.
- Report discrimination, calibration, proper scoring rules, absolute error,
  and uncertainty coverage by survey programme, region, year, and DHW range.

## 5. Reduce the model set to defensible hypotheses

- Start with a parsimonious hierarchical response curve centred on DHW and a
  small, preregistered set of physically motivated modifiers.
- Compare nonlinear or tree-based models only through the same nested blocked
  resampling. Treat feature importance and partial dependence as descriptive,
  not causal evidence.
- Quantify collinearity, interaction support, stability across folds/years,
  and sensitivity to outcome aggregation and predictor windows.

## 6. Quantify predictive and extrapolation uncertainty

- Propagate parameter, process, observation, and environmental-input
  uncertainty into reef-level predictions.
- Map environmental novelty and distance from the training domain; suppress or
  flag predictions outside supported DHW, cloud, water-clarity, current, depth,
  community, or geographic ranges.
- Calibrate prediction intervals and summarise uncertainty at reef, region,
  and GBR scales without presenting smoothed means as deterministic mortality.

## 7. Refactor the reports into a research-grade narrative

- Keep exploratory diagnostics separate from confirmatory modelling and final
  evaluation.
- Generate one canonical table/figure per claim from reusable functions;
  remove duplicated model fits, stale cached results, and contradictory labels.
- Lead with the most stable, out-of-sample findings and clearly distinguish
  association, prediction, and mechanistic interpretation.
