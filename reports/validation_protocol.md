# Prespecified validation protocol

This protocol defines how models in `02_DHW_Mortality_Modelling.qmd` will be
selected and evaluated after the outcome and environmental-data rebuilds. Fold
assignments are created before fitting by `scripts/build_validation_splits.R`.

## Data included

The primary predictor models use the independently reconstructed environmental
table, which covers the 2016, 2017, 2020, 2022, and 2024 bleaching events.
Programme-specific outcomes remain separate:

| Programme | Modelling rows | Reefs | Events |
|---|---:|---:|---|
| LTMP manta tow | 248 | 110 | 2016, 2017, 2020, 2022, 2024 |
| LTMP benthic | 199 | 80 | 2016, 2017, 2020, 2022, 2024 |
| MMP inshore benthic | 90 | 23 | 2017, 2022, 2024 |

The 1998 and 2002 outcomes remain available for DHW-only historical sensitivity
analyses, but they are not mixed into models requiring the modern environmental
predictor suite.

Environmental data are joined first by exact reef name and event year. Four MMP
rows use a documented reef-ID fallback because the North Keppel reef name differs
between registries. The join is checked not to duplicate or discard observations.

Pre-event community composition is reconstructed at the year of the observed
baseline cover. Exact reef-depth-year Acropora observations are preferred,
followed by earlier observations at the same reef and then a past-only
spatial-temporal estimate. Future composition observations are prohibited.
Both Acropora proportion and absolute Acropora cover are carried into every
fold. All validation rows now have an estimate, with source and latest source
year retained for sensitivity analysis.

## Use of the 2024 event

The 2024 event is too widespread, intense, and informative to exclude from the
final ensemble fit. It is therefore treated as a prespecified
leave-one-event-out fold rather than a permanent holdout. Its fold tests how a
model trained on earlier events transfers to 2024; the other event folds test
the reverse transfer when 2024 is part of training.

Candidate specifications are compared using predictions made while each event,
reef, or region is held out. After the full architecture is selected, the BRMS,
BRT, and ensemble components are refitted using every available event, including 2024. This improves
coverage of the observed DHW range but means there is no untouched final event;
performance claims must therefore be described as blocked cross-validation,
not external validation.

## Development-set validation schemes

Three complementary schemes are saved. They answer different questions and
should not be collapsed into a single randomly shuffled score.

1. `leave_one_event_out`: holds out an entire bleaching year. This is the main
   estimate of performance for a future event. Manta and LTMP have five folds;
   MMP has only 2017, 2022, and 2024 and is therefore much less stable.
2. `reef_blocked_5fold`: holds all observations from a reef together in one of
   five balanced folds. This tests transfer to reefs not used in fitting.
3. `region_blocked`: holds out Northern, Central, or Southern GBR as a whole
   where represented. This is a deliberately difficult geographic stress test,
   not the sole model-selection score.

The manifest reports the analysis and assessment sizes for every fold and
contains a passed leakage check for the relevant event, reef, or region grouping.

## Fold-local preprocessing

No predictor is globally imputed or scaled in the split-building script. For
each outer fold, the modelling code must:

1. estimate imputation values from the analysis rows only;
2. estimate centring, scaling, transformations, and factor levels from analysis
   rows only;
3. conduct feature selection and hyperparameter tuning using only inner splits
   of the analysis rows; and
4. apply the frozen preprocessing and fitted model to the assessment rows.

This order applies to all environmental and community predictors. Secchi depth
uses IMOS MODIS-Aqua Kd490 throughout. MMP centroid gaps caused by masked
nearshore pixels have been recovered from nearby valid water pixels with the
extraction radius and support count retained; no eReefs Kd490 values or global
median substitutions are used.

## Performance reporting

For each programme, report out-of-fold predictions and at least:

- MAE and RMSE for point accuracy, with predictive R2 on the same held-out rows;
- balanced MAE giving zero and positive observations equal weight, plus MAE for
  observations with at least 20% mortality;
- Brier score for occurrence of any mortality from the boundary component;
- mean predicted versus mean observed mortality;
- calibration intercept and slope where the assessment sample supports them;
- interval coverage and mean interval width for probabilistic models; and
- sample size alongside every result.

Also summarise errors by event, GBR region, and DHW range when a group contains
at least 10 observations. In-sample fit, ordinary random-fold CV, and internal
GBM CV must not be presented as evidence of generalisation.

Within-family explanatory summaries are also reported: occurrence-component
AICc and McFadden R2 plus positive-component deviance explained for the two-part
GLMs, and marginal/conditional Bayesian R2 for BRMS. AIC is not compared across
quasibinomial, BRT, and Bayesian model families.

## Saved products

- `data/processed/validation_rows_manta.{csv,rds}`
- `data/processed/validation_rows_ltmp.{csv,rds}`
- `data/processed/validation_rows_mmp.{csv,rds}`
- `data/processed/validation_split_manifest.csv`

The RDS files preserve dates and types for modelling. CSV files and the manifest
provide a human-readable audit trail.
