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

## 2. Add a clear modelling-table contract — completed

- Keep this boundary lightweight: no new workflow framework or wholesale
  notebook refactor.
- Validate the assembled environmental table before it is consumed: required
  columns, spatially unique reef-feature-year keys, the five event years,
  complete feature-year coverage, broad physical ranges, cloud provenance,
  and explicit missingness limits.
- Write a small JSON manifest containing the input checksum, dimensions,
  per-column missingness/ranges, cloud coverage by year, warnings, and errors.
- Retain known remote-sensing gaps as visible warnings: 0.42% for the annual
  SST/DHW-derived fields and 0.07% for the neighbourhood-recovered IMOS Secchi
  fields in the current table.
- Run the same validator automatically at the end of environmental assembly
  and independently with `python scripts/validate_environmental_data.py`.

Implementation: `scripts/validate_environmental_data.py`,
`tests/test_environmental_contract.py`, and the validation call in
`scripts/fetch_dms_environmental_data.py`.

## 2a. Recover masked MMP Secchi values — completed

- Use the IMOS MODIS-Aqua Kd490 product consistently across all event years. It
  provides the satellite-derived water-clarity signal needed for 2024 and is
  preferred here to the modelled eReefs chlorophyll/Kd490 fields, which are
  known to be less accurate for this application.
- Diagnose missing MMP values before imputation. All affected survey rows
  matched the intended environmental reef by exact name; the missing values
  occurred because reef centroids fell on masked land or shallow-reef pixels.
- Preserve every valid centroid estimate. Only for a missing centroid, compute
  `1.7 / mean(K_490)` across valid Q1 observations within 2 km, expanding to 3
  or 5 km only if needed. Record the radius, grid-cell count, observation count,
  and nearest valid-water distance.
- Recover all 17 missing MMP reef-year combinations (44 observation rows). All
  were supported within 2 km; the nearest valid water was 0.60--1.75 km away.
  MMP validation now has complete Secchi coverage (0 of 90 rows missing).

Implementation: `scripts/secchi_extraction.py`,
`scripts/recover_mmp_secchi.py`, `tests/test_secchi_extraction.py`, and the
Secchi extraction call in `scripts/fetch_dms_environmental_data.py`.

## 3. Harmonise outcomes and observation processes — completed

- Use proportional loss of total hard-coral cover relative to observed
  pre-event cover as the primary outcome, bounded to 0--1. Retain absolute loss
  in percentage points as the secondary outcome. Keep the former Gompertz
  growth-adjusted outcome and its component cover estimates as sensitivity
  variables rather than building them into the primary estimand.
- Preserve the established May-to-April event window. May--December surveys
  map to the current bleaching year and January--April surveys to the preceding
  bleaching year, allowing the required post-bleaching mortality assessment
  delay while avoiding attribution to the next event.
- Require the most recent baseline to be earlier than the outcome survey and no
  more than 24 months old. Retain exact baseline and outcome dates plus the
  interval in days.
- Keep LTMP manta tow, LTMP benthic, and MMP inshore benthic observations in
  separate primary-analysis tables at the finest supplied observation level.
  Any pooled hierarchical analysis remains a later sensitivity analysis.
- Retain genuine boundary values without beta-regression nudging. Cover gains
  map to zero mortality but are explicitly flagged; complete loss maps to one.
- Write programme-specific CSV and RDS tables plus a compact build manifest.

Implementation: `scripts/harmonise_mortality_outcomes.R` and the harmonisation
chunk in `01_DHW_Mortality_Exploratory.qmd`.

## 4. Build leakage-resistant validation — completed

- Treat 2024 as a prespecified leave-one-event-out fold rather than a permanent
  holdout. This tests transfer from earlier events to 2024 while allowing the
  final production model to learn from the most widespread, intense event.
- Use three explicit development-set checks: leave-one-event-out validation,
  five balanced reef-blocked folds, and leave-one-region-out stress tests.
  Keep all observations from a reef, event, or region together as appropriate.
- Save row-level assignments separately for manta, LTMP benthic, and MMP
  inshore data, plus a manifest giving every fold's analysis/assessment size
  and a passed leakage assertion.
- Join environmental predictors without multiplying observations: exact reef
  name is primary and a visible reef-ID fallback handles four North Keppel MMP
  rows. Missing predictors are retained rather than globally imputed.
- Reconstruct preceding-year Acropora composition against the actual baseline
  survey. Prefer exact reef-depth-year observations, then past observations at
  the same reef, then past-only spatial-temporal interpolation. Retain source,
  latest source year, and interpolation flags; prohibit future observations.
- Require all imputation, scaling, feature selection, and tuning to occur using
  analysis rows inside each fold. Report out-of-fold MAE, RMSE, calibration,
  interval coverage, and subgroup sample size; random-fold and in-sample scores
  are not generalisation evidence.
- Treat MMP validation as low-information: it contains 90 observations across
  23 reefs and only three event years. Secchi depth is now complete, although
  the recovered 2024 Stingaree estimate has only six valid Q1 observations and
  should be retained as a low-support sensitivity case.

Implementation: `scripts/build_validation_splits.R`,
`reports/validation_protocol.md`, and the validation-design loader in
`02_DHW_Mortality_Modelling.qmd`.

## 5. Establish interpretable benchmark models — completed

- Compare five nested benchmark hypotheses: linear DHW, nonlinear DHW,
  nonlinear DHW with optical modifiers, a physical-context extension, and a
  composition extension containing preceding-year Acropora proportion,
  absolute cover, and DHW--Acropora interaction. Use these to quantify the
  value of added context, not as the final limit on predictive complexity.
- Represent the boundary-heavy response with a two-part model: binomial
  occurrence of any cover loss and fractional-logit magnitude conditional on
  positive loss. Retain exact zero and one outcomes without nudging.
- Learn median imputation and standardisation from each analysis fold only.
  Compare candidates through the fixed event-, reef-, and region-blocked
  assignments. Record the simplest candidate within 2% of minimum RMSE for the
  explanatory comparison, but use minimum RMSE for the predictive benchmark.
- Select the physical-context predictive benchmark for LTMP benthic and manta,
  and linear DHW for MMP. The composition-rich GLM improves occurrence AICc and
  apparent fit but overfits blocked folds, especially in MMP; Acropora is
  retained through regularisation in BRMS/BRT rather than removed.
- Preserve the 2024 event-fold result: models trained only on prior events
  materially underpredict mean 2024 manta and MMP mortality. Refit the selected
  benchmarks on all modern observations, including 2024, without
  presenting that refit as independent validation.

Implementation: `scripts/fit_candidate_models.R` and the prespecified candidate
model section in `02_DHW_Mortality_Modelling.qmd` and
`reports/model_validation_report.qmd`.

## 6. Build and validate the formal BRMS--BRT learners — completed

- Refit a nonlinear, boundary-aware hierarchical BRMS learner using the
  canonical programme-specific outcomes, all attainable environmental
  predictors, event effects, and reef/region structure. Use regularising priors
  rather than removing accessible predictors solely for parsimony.
- Include preceding-year Acropora proportion and absolute Acropora cover in
  both formal frameworks, with a prespecified DHW--Acropora interaction in BRMS.
- Refit a regularised BRT learner with the same predictor contract. Tune tree
  number, interaction depth, learning rate, and minimum node size within the
  analysis portion of each blocked fold.
- Generate matched out-of-fold BRMS and BRT predictions for every event-, reef-,
  and region-blocked assessment row. Compare both learners directly with the
  benchmark under identical rows and metrics.
- Learn non-negative ensemble weights from cross-fitted base predictions so an
  assessment outcome never determines its own weight. Retain a base learner
  alone when ensembling does not improve blocked predictive accuracy.
- The final BRMS specification is zero-inflated beta for manta/LTMP and
  zero-one-inflated beta for MMP. Annual maximum DHW is nonlinear, MMP includes
  depth, and prespecified DHW interactions cover cloud, Secchi depth, current,
  and prior severe heat exposure. Held-out reefs never borrow fitted reef
  effects; event and region effects are included only when that event or region
  is known for the stated prediction target.
- The final BRT is a tuned two-part occurrence/magnitude learner. Its
  hyperparameters are selected inside each outer analysis fold and its
  production refit uses all retained modern observations.
- Complete 36 outer BRMS fits and 36 matched BRT fits. Targeted four-chain BRMS
  refits leave all folds with R-hat <= 1.01, bulk ESS >= 167, tail ESS >= 301,
  zero divergences, and zero tree-depth hits. The three production BRMS fits are
  also diagnostically clean.
- Retain both a BRMS--BRT-only blend and a three-member
  benchmark--BRMS--BRT stack. Meta-learner weights are cross-fitted by outer
  fold; production weights use only event- and reef-blocked predictions, never
  region stress tests.
- Route on multiple held-out diagnostics. RMSE remains primary for expected
  mortality, alongside predictive R2, balanced MAE, severe-event MAE, and the
  occurrence Brier score. Current routes are the benchmark--BRMS--BRT stack for
  LTMP and the BRMS--BRT ensemble for manta. MMP uses the BRMS--BRT ensemble
  provisionally: it is within 1% of the lowest RMSE and improves balanced MAE
  and occurrence Brier score while retaining both formal frameworks.
- Report occurrence AICc/McFadden R2 and positive deviance explained for the
  two-part GLMs, and marginal/conditional Bayesian R2 for BRMS. Do not compare
  AIC across quasibinomial, BRT, and Bayesian likelihood constructions.
- Refit every production component on all modern events, including 2024. GBR
  prediction surfaces are deliberately deferred until the uncertainty and
  environmental-novelty contract in Stage 7 is implemented.

Implementation: `scripts/formal_model_helpers.R`,
`scripts/fit_formal_brt.R`, `scripts/fit_formal_brms.R`,
`scripts/summarise_formal_models.R`, `scripts/build_formal_ensemble.R`, and the
rendered `reports/formal_model_report.html`.

## 6a. Add event context and separate forecasting from reconstruction — completed

- Reconstruct ten-year thermal history from annual NOAA CRW maximum DHW. Retain
  cumulative load above 4 DHW, exposure counts, the preceding-window maximum
  and mean, and signed novelty (current maximum minus the preceding ten-year
  maximum). Fit cumulative load and signed novelty as complementary predictors.
- Add an acute IMOS optical-extreme predictor from `1.7 / Q90(Kd490)` during
  January--March, using the same valid-water neighbourhood rules as mean
  Secchi. Treat it as a plume/turbidity proxy, not direct salinity.
- Keep eReefs salinity as a later sensitivity because the public v2 product
  ends during January 2024 and cannot consistently represent the full event.
- Expand the BRMS zero-boundary formula with heat load, signed novelty, acute
  optical conditions, event effects, and region effects. Retain zero-inflated
  beta for manta/LTMP and zero-one-inflated beta for MMP.
- Distinguish prediction targets. Future-event forecasts exclude event and
  region effects. Known-event reef mapping includes event and region effects;
  region-blocked stress tests include the event effect but exclude the held-out
  region effect.
- Compare Gaussian and logit-Gaussian bounded positive-magnitude BRT responses
  within inner folds. The bounded response wins 6 of 36 outer folds and is
  retained only where it improves original-scale expected-mortality RMSE.
- Complete targeted four-chain refits for marginal BRMS folds. All 36 outer
  fits now pass R-hat <= 1.01, bulk ESS >= 167, tail ESS >= 301, zero
  divergences, and zero tree-depth hits; all production fits are also clean.
- For known-event 2024 reconstruction, BRMS predicts programme means of 0.112
  (LTMP benthic), 0.171 (manta), and 0.246 (MMP), versus observations of 0.113,
  0.191, and 0.238. Leave-one-event-out forecasting remains harder, and 2022
  positive magnitude remains an explicit calibration gap.

Implementation: `scripts/secchi_extraction.py`,
`scripts/fetch_dms_environmental_data.py`, `scripts/formal_model_helpers.R`,
`scripts/fit_formal_brms.R`, `scripts/fit_formal_brt.R`,
`scripts/build_formal_ensemble.R`, and `reports/formal_model_report.qmd`.

## 6b. Test direct salinity and the 2024 optical reconstruction — completed

- Extract Q1 surface salinity for all 7,063 GBR registry features from the
  AIMS eReefs GBR4 v2 daily aggregation. Retain mean, Q10, minimum, days below
  26/28 PSU, and accumulated daily deficits below both thresholds. The output
  contains 35,315 reef-year rows for 2016, 2017, 2020, 2022, and 2024.
- Mark 2024 as partial (1--17 January; 17 of 91 requested days) and retain the
  documented 2022 river-forcing warning. Neither period is used as an ordinary
  complete-Q1 calibration event.
- Compare eReefs and IMOS at unique eReefs grid-cell/year combinations to avoid
  pseudo-replicating features that share a 4-km cell. Use complete, unflagged
  2016, 2017, and 2020 events for leave-one-event-out validation.
- Reject quantitative salinity imputation. Kd490-plus-spatial minimum-salinity
  predictions have only 0.34--0.41 PSU all-cell MAE, but miss the ecologically
  decisive flood tail: flood-cell MAE is 6.5--10.3 PSU, bias is always positive,
  and no held-out cell exposed below 28 PSU is predicted below 28 PSU.
- Retain an exploratory risk layer with Kd490-only, spatial-only, and combined
  probabilities kept separate. Kd490-only transfer is inconsistent and reverses
  in 2020. All three low-salinity cells in the partial 2024 eReefs record rank
  above the 99.5th percentile in the combined layer, but two rank below the
  15th percentile using Kd490 alone, showing that geography drives much of the
  apparent success.
- Do not add model-derived salinity to the formal mortality predictors. Keep
  the directly observed IMOS acute optical predictor and expose the 2024 risk
  layers only for sensitivity analysis, with
  `use_in_primary_mortality_model = FALSE`.

Implementation: `scripts/extract_ereefs_salinity.py`,
`reports/salinity_kd490_comparison.qmd`, the rendered
`reports/salinity_kd490_comparison.html`, and
`data/processed/ereefs_salinity_2024_optical_sensitivity.csv`.

## 6c. Represent wind and doldrum-like conditions — source audit completed

- Audit the linked IMOS Sentinel-1 SAR wind product across every Q1 2024 daily
  file. It provides calibrated 10-m neutral wind at approximately 1-km
  resolution but starts in October 2017 and retains satellite-swath gaps. Mean
  daily GBR-domain coverage is 3.9%; the median registry feature has seven
  observed days in Q1, so it cannot estimate calm duration.
- Extract calibrated, good-QC Metop-B scatterometer wind for all event years,
  aggregating observations within 30 km to daily medians. This produces useful
  offshore comparisons, but adequate coverage is limited to 60/199 LTMP and
  86/248 manta model rows and 0/90 MMP rows. Do not substitute distant offshore
  pixels for inshore wind.
- Extract the eReefs GBR4 daily `mean_wspeed` forcing for all 7,063 registry
  features and derive Q1 mean, Q10, minimum, fractions/days below 3 and 5 m/s,
  and longest calm spells. Coverage is complete for 2016, 2017, 2020, and 2022;
  2024 is explicitly partial through 17 January.
- Product agreement is strongest for broad mean wind, while event-specific
  lower-tail and calm-fraction agreement is weak. Do not use Metop or
  Sentinel-1 to impute missing February--March 2024 eReefs doldrum metrics.
- Before adding wind to BRMS/BRT, extract one continuous hourly 10-m wind
  reanalysis (ERA5 or matched ACCESS) for every event, validate its mean and
  low-wind tail against eReefs and good-QC Metop, then test a compact mean-plus-
  calm predictor set under the existing outer folds.

Implementation: `scripts/audit_sar_wind_coverage.py`,
`scripts/extract_metop_wind.py`, `scripts/extract_ereefs_wind.py`,
`reports/wind_doldrum_assessment.qmd`, and the rendered
`reports/wind_doldrum_assessment.html`.

## 7. Quantify predictive and extrapolation uncertainty — contract completed; mapping gated

- Relabel the saved BRMS Q05--Q95 limits correctly as posterior intervals for
  expected mortality (`posterior_epred`), not observation-level prediction
  intervals. Their held-out empirical coverage is only 20--54% across
  programmes and schemes.
- Add 90% split-conformal intervals calibrated with residuals from other outer
  folds of the same programme and scheme. Event-held-out marginal coverage is
  90% for LTMP, 84% for manta, and 81% for MMP. Report positive and severe rows
  separately so good zero-boundary coverage cannot hide high-loss failures.
- Add a severe-event upper guardrail calibrated only from mortality of at least
  20%. It reaches 91--93% event-held-out severe coverage for LTMP and manta but
  only 51% for MMP; MMP temporal uncertainty remains unresolved rather than
  being widened until it looks acceptable.
- Carry absolute BRMS--BRT disagreement as a local epistemic diagnostic. Do not
  treat agreement as proof of low uncertainty because both learners can omit
  the same event mechanism.
- Screen all 7,063 registry features against each programme's 15-variable
  environmental training domain using univariate range checks and nearest-row
  standardised distances. Only 18 registry features are within MMP support, so
  MMP is an inshore target and must not be mapped as a general GBR surface.
- Gate all production mortality surfaces because full-registry, past-only
  `prop_acropora_pre` and absolute `acropora_cover_pre` are not yet audited;
  MMP also needs depth. Every future prediction row must carry component
  predictions, route, interval, severe guardrail, disagreement, support flags,
  community-input provenance, and suppression status.

Implementation: `scripts/quantify_prediction_uncertainty.R`,
`output/prediction_uncertainty/`, `reports/prediction_uncertainty_report.qmd`,
and the rendered `reports/prediction_uncertainty_report.html`.

## 8. Refactor the reports into a research-grade narrative

- Keep exploratory diagnostics separate from confirmatory modelling and final
  evaluation.
- Generate one canonical table/figure per claim from reusable functions;
  remove duplicated model fits, stale cached results, and contradictory labels.
- Lead with the most stable, out-of-sample findings and clearly distinguish
  association, prediction, and mechanistic interpretation.
