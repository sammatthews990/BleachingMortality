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
  26/28/30 PSU, and accumulated daily deficits below all three thresholds. Use
  30 PSU as the primary freshwater-exposure screen, while retaining 28 and 26
  PSU as nested severe thresholds. The output
  contains 35,315 reef-year rows for 2016, 2017, 2020, 2022, and 2024.
- Mark 2024 as partial (1--17 January; 17 of 91 requested days) and retain the
  documented 2022 river-forcing warning. Neither period is used as an ordinary
  complete-Q1 calibration event.
- Compare eReefs and IMOS at unique eReefs grid-cell/year combinations to avoid
  pseudo-replicating features that share a 4-km cell. Use complete, unflagged
  2016, 2017, and 2020 events for leave-one-event-out validation.
- The 30-PSU screen yields 54 unique historical eReefs cell-events, compared
  with 27 below 28 PSU, substantially increasing the rare-event training set.
- Reject quantitative salinity imputation. Kd490-plus-spatial minimum-salinity
  predictions have only 0.34--0.41 PSU all-cell MAE, but miss the ecologically
  decisive flood tail: flood-cell MAE is 6.5--10.3 PSU, bias is always positive,
  and no held-out cell exposed below 28 PSU is predicted below 28 PSU.
- Retain an exploratory risk layer with Kd490-only, spatial-only, and combined
  probabilities kept separate. Kd490-only transfer is inconsistent and reverses
  in 2020. All three cells below 28 PSU in the partial 2024 eReefs record rank
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

## 6c. Represent wind and doldrum-like conditions — completed

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
- Extract hourly ERA5 10-m wind for every registry feature and event from the
  public AWS Earthmover ERA5 archive. This provides complete Q1 coverage for
  2016, 2017, 2020, 2022, and 2024, unlike the partial eReefs 2024 record.
- Validate ERA5 on unique eReefs grid cells. Event-level Spearman correlations
  are 0.91--0.95 for Q1 mean wind, 0.80--0.91 for Q10 wind, and 0.84--0.90 for
  the fraction below 3 m/s; mean-wind RMSE is 0.41--0.48 m/s. ERA5 therefore
  passes the product-agreement gate.
- Carry a compact ERA5 mean-wind plus calm-fraction pair into the next formal
  BRMS/BRT candidate refit. Do not include multiple strongly redundant calm
  summaries merely because they are available.

Implementation: `scripts/audit_sar_wind_coverage.py`,
`scripts/extract_metop_wind.py`, `scripts/extract_ereefs_wind.py`,
`reports/wind_doldrum_assessment.qmd`, and the rendered
`reports/wind_doldrum_assessment.html`.

## 6d. Test rainfall as freshwater-event context — completed

- Extract daily ERA5 precipitation for reef cells and their nearest ERA5 land
  cells, retaining December--March total rainfall and maximum 1-, 7-, and
  30-day accumulations. The same file contains the validated wind summaries
  and covers all 7,063 registry features in all five modern bleaching events.
- Use leave-one-event-out validation on unique eReefs cells from complete,
  unflagged 2016, 2017, and 2020 events. The combined rainfall, optical, wind,
  and spatial classifier ranks all four partial-record 2024 cells below 30 PSU
  above the 98th percentile of predicted risk.
- Do not convert those ranks to quantitative imputed salinity: no held-out cell
  below 30 or 28 PSU is reconstructed below its threshold, and flood-cell
  minimum-salinity predictions are biased high by 3.7--8.4 PSU. Retain direct
  ERA5 coastal rainfall as a candidate mortality predictor and test it only
  inside the existing event-held-out folds.
- Use the public ERA5 archive consistently. The tested Open-Meteo ERA5 mirror
  differed enough in rainfall to be rejected, while ACCESS operational forecast
  products do not provide a single public, consistent reanalysis for this
  historical extraction.
- In the next formal candidate refit, remove `dist_to_er_km`: it measures the
  reef-to-eReefs grid-match distance and is a quality-control field, not an
  ecological predictor. Test the true ERA5 land-cell distance together with a
  compact four-variable addition: coastal maximum 30-day rainfall, reef mean wind,
  reef calm fraction, and coastal distance. Judge retention entirely within
  the existing outer validation folds.

Implementation: `scripts/extract_era5_weather.py`,
`scripts/model_rainfall_salinity.R`,
`data/processed/era5_weather_reef_year.csv`,
`data/processed/ereefs_salinity_2024_rainfall_sensitivity.csv`,
`reports/rainfall_salinity_reconstruction.qmd`, and the rendered
`reports/rainfall_salinity_reconstruction.html`.

## 6e. Add ENSO event context — extraction completed; model sensitivity pending

- Download NOAA's seasonal Relative Oceanic Niño Index (RONI) and the Bureau
  of Meteorology's monthly Troup Southern Oscillation Index (SOI), caching the
  raw source files. Summarise RONI over DJF, JFM, and FMA and SOI over December
  to March for each bleaching-event year.
- Use continuous mean summer RONI as the only primary model candidate. Retain
  Strong/Moderate/Weak El Niño/La Niña labels for interpretation and SOI as an
  atmospheric cross-check; do not spend degrees of freedom on both indices or
  fit sparse categories as separate coefficients.
- The modern events cover Strong El Niño (2016), Neutral (2017 and 2020),
  Moderate La Niña (2022), and Weak El Niño (2024). The 2024 oceanic signal is
  not coupled in the December--March mean SOI. ENSO therefore distinguishes
  2022 from 2020 rather than supplying one common explanation for both low-
  mortality events.
- Treat the ENSO coefficient as a tightly regularised event-level sensitivity.
  There are only five modern events, so it is highly confounded with event
  identity and must not be promoted from apparent in-sample fit. Compare it
  only through leave-one-event-out predictions and do not include an event
  random effect in the same future-event prediction specification.

Implementation: `scripts/extract_enso_indices.py` and
`data/processed/enso_event_context.csv`.

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

## 7a. Diagnose ecological prediction failures and response family — completed

- Assess signed residuals only from leave-one-event-out predictions. Calibrate
  each event's 5th/95th residual tails and 90th absolute-error threshold using
  other events from the same programme and learner. Flags are ecological review
  targets and do not remove observations from model fitting.
- Retain all 96 flagged survey rows but provide a deduplicated review table of
  82 programme--reef--years. Forty unique reef-years are severe 2024 misses
  (5 LTMP, 25 manta, and 10 MMP); three are large overpredictions below 4 DHW.
  All three core learners underestimate 39 of the priority reef-years.
- Separate candidate 2024 mechanisms spatially. Northern/central severe misses
  coincide with high coastal rainfall, while southern misses combine higher
  DHW with high pre-event Acropora. Test both mechanisms in blocked prediction
  rather than forcing one steeper DHW curve to explain all regions.
- Fit a like-for-like Bayesian binomial sensitivity with 100 pseudo-trials and
  an observation-level random effect. This is not a literal colony-count model:
  the response is loss relative to estimated pre-event cover. The primary
  programme-specific zero-inflated/zero-one-inflated beta family has lower
  event-held-out RMSE in LTMP (0.133 versus 0.167), manta (0.141 versus 0.176),
  and MMP (0.243 versus 0.281).
- The binomial family raises manta's held-out 2024 mean prediction from 0.107 to
  0.168 (observed 0.191), but worsens its 2024 RMSE and overpredicts the low-
  mortality events. It barely changes the MMP 2024 miss (0.121 to 0.128 versus
  observed 0.238). It therefore does not solve the omitted-2024-mechanism
  problem.
- Keep the binomial fit as a sensitivity only. Its manta-2022 reef-variance
  parameter improved from R-hat 1.020 to 1.013 after an extended four-chain
  refit, but remains above the prespecified 1.01 gate. Effective sample sizes
  are adequate and there are no divergences or tree-depth hits; all other
  folds pass every diagnostic gate. The remaining failure is reported rather
  than hidden or sampled repeatedly until it happens to pass.

Implementation: `scripts/assess_prediction_outliers.R`,
`scripts/fit_binomial_family_sensitivity.R`,
`output/prediction_outliers/`, `reports/prediction_outlier_assessment.qmd`, and
the rendered `reports/prediction_outlier_assessment.html`.

## 7b. Restore an ecological high-DHW response - first candidate completed

- Preserve the original formal BRMS/BRT fits as benchmarks. Fit a separate
  shape-aware screen that retains preceding proportional Acropora as the
  susceptibility variable, replaces the collinear absolute Acropora term with
  total pre-event coral cover, adds DHW--Acropora effects to both mortality
  occurrence and positive magnitude, and prevents event identity from
  absorbing a response that must transfer to a future event.
- Use directional rather than outcome-forcing Bayesian priors. The BRT
  candidate constrains DHW and proportional Acropora to be non-decreasing,
  gives each reef-year equal fitting weight, and tunes against overall,
  event-balanced, and high-DHW error.
- The candidate is not a universal replacement. The shape-aware ensemble
  improves LTMP event-held-out RMSE from 0.131 to 0.128 and predictive R2 from
  0.388 to 0.408, but overpredicts the high-DHW/high-Acropora subset. It worsens
  manta RMSE from 0.141 to 0.149, so retain the existing manta benchmark until
  a shared ecological response with a programme-specific observation layer is
  fitted.
- The MMP shape-aware BRMS is a material improvement: event-held-out RMSE falls
  from 0.243 to 0.216, predictive R2 rises from -0.098 to 0.129, severe-loss
  MAE falls from 0.288 to 0.207, and high-DHW/high-Acropora bias is nearly zero
  (-0.010 observed minus predicted). With 2024 held out, its event mean rises
  from 0.121 in the benchmark BRMS to 0.203 versus 0.238 observed.
- In the all-event production scenario at 12 DHW, 90% proportional Acropora,
  and 55% total pre-event coral cover, the shape-aware BRMS predicts 0.562 LTMP
  and 0.609 MMP mortality. The shape-aware LTMP ensemble falls to 0.383 because
  its BRT component cannot extrapolate beyond LTMP's observed DHW range. Treat
  BRT weight as zero outside programme-specific training support rather than
  allowing a clamped tree response to flatten extrapolated mortality.
- All 13 shape-aware Bayesian event folds pass R-hat <= 1.01, effective sample
  size, divergence, and tree-depth gates after targeted four-chain refits.
  The implemented production route is programme-specific: shape-aware
  ensemble for LTMP within its observed DHW support and shape-aware BRMS beyond
  it; benchmark ensemble for manta within support and benchmark BRMS beyond it;
  and shape-aware BRMS for MMP throughout. At 12 DHW this routed predictor gives
  0.562 LTMP, 0.254 manta, and 0.609 MMP mortality in the dominant-Acropora
  scenario. Promotion remains conditional on spatial validation of the 2024
  compound-event mechanisms.

Implementation: `scripts/shape_aware_model_helpers.R`,
`scripts/fit_shape_aware_brms.R`, `scripts/fit_shape_aware_brt.R`,
`scripts/compare_shape_aware_models.R`, `output/shape_aware_models/`,
`reports/shape_aware_model_comparison.qmd`, and the rendered
`reports/shape_aware_model_comparison.html`.

## 7c. Joint observation and freshwater-event sensitivity - completed

- Stack all 537 LTMP, manta, and MMP observations while grouping 152 paired
  programme observations into shared reef-events. Fit a common ecological
  response with programme-specific boundary behaviour, complete-loss
  probability, mean offset, and precision. Keep disturbance codes as
  interpretation fields only and do not change selection logic.
- Add a compact ERA5 exposure block: nearest-coastal maximum 30-day rainfall,
  mean reef wind, calm fraction below 3 m/s, coastal distance, and DHW
  interactions with rainfall and calm fraction. Do not label ERA5 rainfall as
  measured salinity and do not use the partial January 2024 eReefs salinity flag
  as a mortality predictor.
- Reject both joint candidates for future unseen events. Relative to the
  current routed benchmark, event-held-out RMSE and severe MAE are worse for
  every programme. The all-event BRMS supports a positive rainfall main effect
  on positive mortality (estimate 0.381, 95% interval 0.173--0.596), but not a
  positive DHW--rainfall interaction (estimate -0.102, interval
  -0.270--0.065).
- Retain the joint compound BRT as a 2024 spatial sensitivity only. Against an
  otherwise identical joint BRT without weather, manta spatial RMSE improves
  from 0.161 to 0.151 and severe MAE from 0.179 to 0.165. LTMP severe MAE also
  improves, but overall error is unchanged; MMP worsens. The joint BRMS does
  not reproduce the manta gain.
- Weather raises predictions most at several wet/calm northern and central
  manta reefs, including MacGillivray, Swinger, South Direction, Hastings,
  Lizard Island NW, St Crispin, Mackay, and Eyrie. This is consistent with an
  additive freshwater-event mechanism but is framework-specific and does not
  justify promotion or a causal salinity claim.
- Preserve the existing programme-specific production route. Several
  two-chain joint BRMS screening folds exceed the R-hat gate, although all have
  adequate ESS and no divergences/tree-depth hits and the four-chain production
  fit passes. Because predictive promotion already fails, do not spend compute
  on targeted refits.

Implementation: `scripts/joint_compound_model_helpers.R`,
`scripts/fit_joint_compound_brms.R`, `scripts/fit_joint_compound_brt.R`,
`scripts/compare_joint_compound_models.R`,
`output/joint_compound_models/`, and
`reports/joint_compound_stress_assessment.qmd`.

## 7c. Joint observation and freshwater-event sensitivity - completed

- Stack all 537 LTMP, manta, and MMP observations while grouping 152 paired
  programme observations into shared reef-events. Fit a common ecological
  response with programme-specific boundary behaviour, complete-loss
  probability, mean offset, and precision. Keep disturbance codes as
  interpretation fields only and do not change selection logic.
- Add a compact ERA5 exposure block: nearest-coastal maximum 30-day rainfall,
  mean reef wind, calm fraction below 3 m/s, coastal distance, and DHW
  interactions with rainfall and calm fraction. Do not label ERA5 rainfall as
  measured salinity and do not use the partial January 2024 eReefs salinity flag
  as a mortality predictor.
- Reject both joint candidates for future unseen events. Relative to the
  current routed benchmark, event-held-out RMSE and severe MAE are worse for
  every programme. The all-event BRMS supports a positive rainfall main effect
  on positive mortality (estimate 0.381, 95% interval 0.173--0.596), but not a
  positive DHW--rainfall interaction (estimate -0.102, interval
  -0.270--0.065).
- Retain the joint compound BRT as a 2024 spatial sensitivity only. Against an
  otherwise identical joint BRT without weather, manta spatial RMSE improves
  from 0.161 to 0.151 and severe MAE from 0.179 to 0.165. LTMP severe MAE also
  improves, but overall error is unchanged; MMP worsens. The joint BRMS does
  not reproduce the manta gain.
- Weather raises predictions most at several wet/calm northern and central
  manta reefs, including MacGillivray, Swinger, South Direction, Hastings,
  Lizard Island NW, St Crispin, Mackay, and Eyrie. This is consistent with an
  additive freshwater-event mechanism but is framework-specific and does not
  justify promotion or a causal salinity claim.
- Preserve the existing programme-specific production route. Several
  two-chain joint BRMS screening folds exceed the R-hat gate, although all have
  adequate ESS and no divergences/tree-depth hits and the four-chain production
  fit passes. Because predictive promotion already fails, do not spend compute
  on targeted refits.

Implementation: `scripts/joint_compound_model_helpers.R`,
`scripts/fit_joint_compound_brms.R`, `scripts/fit_joint_compound_brt.R`,
`scripts/compare_joint_compound_models.R`,
`output/joint_compound_models/`, and
`reports/joint_compound_stress_assessment.qmd`.

## 7d. RRN water-colour and cyclone-pressure extension - completed

- Extract `wqc_freqcc12` and `cyc_maxHrs4mw` reproducibly from
  `data/GBRMPA_RRN_2025_AIMS.xlsx`. Map austral summer `202324` to event
  year 2024 and preserve both codes in the processed table. All 537 mortality
  observations match directly by GBRMPA `LABEL_ID`; no spatial imputation is
  required.
- Treat water colour as a delayed freshwater/plume proxy and cyclone wave
  exposure as a separate mechanical-disturbance predictor. Do not merge the
  cyclone variable into a salinity interaction and do not change disturbance
  inclusion or attribution logic.
- Validate water colour against complete 2016, 2017, 2020 and 2022 eReefs
  salinity fields. It improves rare below-30 PSU classification over ERA5
  rainfall alone (leave-one-event-out average precision 0.133 versus 0.101);
  water colour plus rainfall reaches 0.146, while IMOS Kd490 remains the best
  single classifier at 0.159. Combining all three gives the best continuous
  minimum-salinity RMSE (0.610 PSU).
- In mortality BRT ablations, water colour gives the lowest unseen-event RMSE
  among the tested joint variants (0.178 versus 0.181 for the core model), but
  the gain is small. Cyclone exposure gives the lowest overall reef-blocked
  2024 RMSE (0.171 versus 0.172), with the clearest programme gain in MMP.
  Rainfall/wind remains best for 2024 manta tow, and the complete
  weather-plus-RRN stack worsens overall 2024 RMSE to 0.176.
- The regularised all-event BRMS supports a positive additive water-colour
  association with positive mortality (estimate 0.310, 95% interval
  0.041--0.574), but not a DHW-by-water-colour interaction (0.022,
  -0.116--0.166) or positive cyclone-wave mean effect (-0.035,
  -0.164--0.092). The four-chain fit passes the diagnostic gate.
- This initial screen is superseded by the corrected matched comparison in
  section 7f. Retain rainfall, raw water colour and 4 m-wave hours together in
  the ensemble candidate, with IMOS Kd490/rainfall as the operational fallback
  during the one-wet-season WQC release lag. Promotion still depends on final
  ensemble validation rather than any single programme result.

Implementation: `scripts/extract_rrn_pressure_data.R`,
`scripts/evaluate_rrn_pressure_proxies.R`,
`scripts/compare_rrn_pressure_models.R`,
`scripts/fit_rrn_pressure_brms.R`,
`data/processed/rrn_pressure_reef_year.csv`,
`output/rrn_pressure_assessment/`, and
`reports/rrn_pressure_assessment.qmd`.

## 7e. Severe 2024 manta reef and community-composition audit - completed

- Extract all six RRN individual-pressure time series and all four metric
  families for Mackay, Linnet, Lizard Island NW, Eyrie and Swinger. Benchmark
  2023--24 pressures against all RRN reefs in the same one-degree latitude
  band and against each reef's preceding decade. Retain but flag the latest
  `cmex`, `cmfi` and `quad` metrics because they include summer 2024--25
  and would leak post-event information into a 2024 prediction.
- Flag Linnet and Eyrie for unusually high low-level COTS pressure: their IDW
  COTS values are at the 98.9th and 99.7th latitude-band percentiles,
  respectively. Do not automatically label these formal outbreaks and do not
  change the existing disturbance-selection logic.
- Flag Linnet, Eyrie and Swinger for coloured-water values at the 95th
  percentile of their own preceding decade. Across all 89 reef-blocked 2024
  manta observations, coloured water correlates 0.170 with the core-model
  residual; its upper quartile averages 5.7 percentage points more
  underprediction than the remainder.
- Flag Lizard Island NW's coral-sink index of 0.105 (15.6th local percentile)
  as a connectivity/recovery vulnerability, not an acute mortality exposure.
  Mackay has no convincing gridded RRN explanation, but its retained survey
  tooltip explicitly records bleaching, Cyclone Jasper and flood waters.
- Confirm that the manta model contains `prop_acropora_pre` for all five
  reefs. Linnet, Lizard NW and corrected Mackay use observed depth/year
  estimates; Eyrie and Swinger use past-only spatial-temporal interpolation.
  In the reef-blocked BRT folds, setting Acropora to 90% still leaves 23--44
  percentage points of underprediction. Acropora uncertainty therefore
  contributes but cannot explain the misses alone.
- Add a broader community-susceptibility/recovery-state predictor to the next
  candidate stage. It should represent other fast-growing or
  disturbance-sensitive taxa/growth forms and recent recovery state, rather
  than treating Acropora as the entire susceptible community. Carry provenance
  and uncertainty into prediction.

Implementation: `scripts/audit_rrn_target_reefs.R`,
`scripts/assess_manta_acropora_counterfactual.R`,
`output/rrn_target_reefs/`, and
`reports/rrn_target_reef_diagnostics.qmd`.

## 7f. Reef-lineage repair and expanded RRN pressure screen - completed

- Replace ambiguous cleaned-name crosswalks with exact-name-first matching.
  The former rule collapsed `MACKAY REEF` (`16-015`) and `MACKAY REEFS`
  (`15-024`) to the same key. Both 2024 Mackay survey programmes now resolve
  to `16-015`, and the pre-event Acropora estimate is the observed 0.111 rather
  than an interpolated value from the wrong reef.
- Carry the raw disturbance tooltip and explicit bleaching, cyclone, flood and
  COTS flags into the canonical outcome and validation rows. The 2024 Mackay
  tooltip now retains `coral bleaching, Cyclone Jasper, flood waters`.
- Extract observed and IDW-modelled COTS pressure alongside cyclone wave-hours
  and WQC. Add a leakage-safe within-reef WQC percentile based only on the
  preceding ten event years; retain raw WQC as a separate absolute-exposure
  measure.
- Compare core, weather, raw-WQC, reef-relative-WQC and COTS BRT candidates in
  matched leave-one-event-out and reef-blocked-2024 folds. Weather + raw WQC +
  cyclone waves gives the best overall reef-blocked-2024 RMSE (0.174) and
  improves unseen-event RMSE relative to weather alone (0.180 versus 0.185).
  The relative WQC candidate has some severe-case benefits but is weaker
  overall.
- Retain COTS in the research dataset but do not require it in the production
  route. Adding modelled COTS worsens blocked BRT prediction. In matched BRMS
  fits its mean effect is uncertain and its LOO gain is 0.1 (SE 1.6), despite a
  negligible Bayesian-R2 increase from 0.334 to 0.338.

Implementation: `scripts/harmonise_mortality_outcomes.R`,
`scripts/extract_rrn_pressure_data.R`, `scripts/build_validation_splits.R`,
`scripts/fit_joint_compound_brt.R`, `scripts/fit_rrn_pressure_brms.R`,
`scripts/compare_rrn_brms_models.R`, and
`reports/rrn_pressure_assessment.qmd`.

## 7g. Shared latent ecological mortality with programme observation layers - completed candidate

- Implement one ecological logit-mortality state per reef-event. Environmental
  and community predictors enter this state once. Manta anchors the scale;
  LTMP and MMP have separate observation bias, beta precision and zero/one
  probabilities.
- Reject programme-specific ecological response loadings because MMP has only
  one paired benthic reef-event and the expanded model fails the sampling gate.
  Retain a common ecological scale with regularised observation biases.
- Validate on matched leave-one-event-out and reef-blocked-2024 observations.
  The latent model improves unseen-event MMP RMSE from 0.324 to 0.271, but
  worsens LTMP and manta and trails the weather--RRN BRT overall. Retain it as
  an MMP calibration/ensemble candidate; do not replace the production route.
- Add a programme-free operational prediction method that excludes fitted
  reef effects and samples a new reef-event effect for GBR-wide uncertainty.

Implementation: `scripts/latent_mortality_model_helpers.R`,
`scripts/fit_latent_mortality_brms.R`,
`scripts/compare_latent_mortality_models.R`,
`output/latent_mortality/`, and
`reports/latent_and_annual_cover_assessment.qmd`.

## 7h. Annual coral-cover dynamics - decomposed BRT candidate completed

- Build 3,089 approximately annual transitions (240--550 days) from the raw
  manta, LTMP and MMP series. Retain positive, negative and zero change. Use
  exact-name-first reef lineage with explicit ambiguous-name rules and Mackay
  collision assertions.
- Predict both next cover and absolute percentage-point change. Include current
  cover, previous trajectory, interval, preceding Acropora, RRN DHW, WQC,
  cyclone waves and COTS. All disturbance types remain in scope because this
  target is net annual change rather than bleaching-attributed mortality.
- Benchmark persistence, a growth/environment linear model, a post-cover BRT
  and a change BRT. The growth model modestly improves temporal/forward RMSE;
  the change BRT improves unseen-reef RMSE and captures the mean 2024 decline,
  but not reef-level extremes. No candidate is ready for operational promotion.
- Fit a predictive decomposition in which gain is scaled by available space
  and loss by current coral cover. Use separate occurrence and bounded
  magnitude BRTs, retaining nonlinear DHW, Acropora, thermal history/novelty,
  Secchi/cloud, rainfall, calm, WQC, cyclone waves and COTS in the loss process.
- Transfer the restricted bleaching-mortality response into a source-balanced
  auxiliary loss learner. Apply the same event and reef exclusions as the
  annual validation fold so 2024 is not leaked into the novel-event test.
- The core decomposition is best when 2024 is absent (forward RMSE 10.05 pp).
  The modifier-rich model is best across unseen reefs (6.76 pp), while the
  mortality-augmented candidate is best for reef-blocked 2024 (8.59 pp) and
  predicts mean 2024 change as -4.52 pp versus -5.89 observed.
- Do not promote the candidate yet: severe 2024 losses average -21.58 pp but
  remain predicted as -11.14 pp. Residual underprediction remains associated
  with calm, rainfall and water colour, and Jasper-annotated reefs can have
  zero RRN cyclone-wave exposure.
- Next fit bounded beta magnitude and binomial occurrence models with
  programme-specific observation layers. Retain the decomposed BRTs as
  ensemble members and cross-fit any routing or stacking weights.

Implementation: `scripts/build_annual_coral_transitions.R`,
`scripts/fit_annual_coral_change_benchmarks.R`,
`scripts/fit_decomposed_annual_change.R`,
`scripts/summarise_decomposed_annual_change.R`,
`data/processed/annual_coral_transitions.*`, and
`output/annual_coral_change/`, `output/decomposed_annual_change/`, and
`reports/decomposed_annual_change_assessment.qmd`.

## 7i. Extreme-loss expert and distributional annual-change benchmark - completed candidate

- Add the RRN coral-sink connectivity index to the annual transition table as
  a recovery predictor. The remaining RRN individual-pressure fields were
  already present. Do not use the temporally composite RRN metrics as acute
  mortality predictors when their component variables are available.
- Fit a Bernoulli severe-loss gate (annual loss of at least 10 percentage
  points) and a severe-magnitude BRT. Train the tail with both severe annual
  transitions and source-balanced restricted bleaching-mortality rows, while
  applying the same reef/year exclusions as every validation fold.
- Blend the tail expert with the mortality-augmented decomposition. In
  reef-blocked 2024 validation, the soft gate improves RMSE from 8.59 to 8.14
  points and severe-loss MAE from 11.43 to 9.67. The aggressive gate improves
  severe MAE to 8.47 but raises false-extreme predictions among non-severe rows
  to 18.2%; retain it as a risk-sensitive scenario rather than the central
  prediction.
- Fit a distributionally appropriate GAM benchmark with binomial gain/loss
  occurrence and conditional beta magnitudes. Retain current cover, available
  space, Acropora, thermal history/novelty, Secchi/cloud, rainfall/calm, WQC,
  cyclone waves, COTS and the restricted-mortality transfer. Treat this as a
  transparent distributional benchmark, not a substitute for the formal
  Bayesian latent observation model. It is competitive but not promoted:
  unseen-reef RMSE is 6.89 pp (predictive R2 0.178), while reef-blocked 2024
  RMSE is 9.48 and its mean prediction for severe cases is -10.29 pp versus
  -21.58 observed. The severe underprediction is therefore not explained only
  by use of an unbounded response distribution.
- Prioritise measured MMP event salinity/TSS/CDOM/Secchi, Queensland river
  discharge, raw AIMS community and agents-of-mortality fields, cyclone
  best-track/rainfall exposure, 2024 aerial bleaching and GBR10 reef habitat.
  Derive event timing/duration rather than relying only on wet-season means.
- Carry the base, soft-gated central candidate and aggressive severe-risk
  scenario forward. Select any final blending weight and threshold with nested
  reef/year validation and an explicit asymmetric cost for missed extremes,
  while reporting the false-extreme rate as a separate constraint.

Implementation: `scripts/fit_extreme_loss_ensemble.R`,
`scripts/fit_beta_binomial_annual_gam.R`,
`output/extreme_loss_ensemble/`, `output/beta_binomial_annual/`, and
`reports/extreme_loss_ensemble_assessment.qmd`.

## 7j. Measured freshwater, cyclone tracks and upwelling relief - data layer completed

- Reconstruct the public AIMS water-quality record rather than reusing the
  existing eReefs salinity field. The API exposes MMP-only/all-programme,
  depth-average/depth-weighted and daily/hourly choices. The reproducible
  all-programme extraction contains 10 logger sites, 6,366 depth-averaged and
  3,868 depth-weighted sample rows across the 2015--2025 wet seasons. It
  supplies salinity, TSS, CDOM, Secchi, chlorophyll and nutrients, including
  943 depth-averaged observations below 30 PSU.
- Record the genuine observation gap: public daily AIMS logger downloads end
  by 25 October 2023. Do not treat the absence of a 2024 logger summary as a
  failed join and do not invent direct Jasper-period logger salinity.
- Link shallow (0--10 m) sample summaries to unique GBR registry points while
  retaining nearest-reef distance and 5/20 km proximity flags. This creates
  3,267 AIMS point-event calibration rows; 2,193 are within 5 km and 3,164 are
  within 20 km of a registry reef point. Exact unique survey-name lineage is
  used before spatial matching. River-mouth observations remain plume
  calibration data and are not relabelled as direct reef measurements.
- Build the corresponding 46,760-row all-GBR reef-event prediction grid for
  2016--2025. Join RRN coloured water, cyclone waves and COTS; ERA5 rainfall
  and wind; IMOS Kd490/Secchi; eReefs salinity; and the new BoM cyclone fields.
  Complete ERA5/eReefs predictor layers currently cover the five principal
  bleaching years 2016, 2017, 2020, 2022 and 2024.
- Parse the complete BoM Australian tropical-cyclone best-track database into
  reef-event proximity, nearest-fix wind/pressure, storms within 300 km and
  quadrant-aware gale-radius exposure. Jasper passes approximately 31 km from
  Mackay 16-015, 68 km from Swinger, 78 km from Mackay 15-024, and 124--132 km
  from Linnet/Lizard/Eyrie. Keep this forcing separate from rainfall/flooding
  and wave damage.
- Generalise freshwater with a two-part model: below-30-PSU occurrence and
  conditional magnitude/freshwater deficit. Calibrate the broad eReefs prior
  against AIMS, then predict from event rainfall, IMOS optical plume anomalies,
  cyclone forcing and river-to-reef connectivity. Keep 2024 in the training
  pool and validate with held-out spatial/catchment groups; report a no-2024
  sensitivity analysis rather than using 2024 only as a holdout.
- The first whole-reef-blocked pilot is complete. The hybrid occurrence model
  has average precision 0.339 versus prevalence 0.184, but conditional
  low-salinity magnitude RMSE remains 11.5 PSU with predictions biased 5.7 PSU
  too high. Do not pass these magnitude estimates into mortality yet. Add
  event timing and river-to-reef connectivity before expanding model
  complexity.
- Derive upwelling relief dynamically. First extract daily IMOS 0.02-degree L3S
  reef SST residuals relative to contemporaneous surrounding waters. Retain
  minimum residual, cool-day counts and cool-spell duration during peak heat.
  Then add a coarser subsurface forcing layer from vertical velocity, currents
  and temperature gradients, calibrated against Palm Passage/Myrmidon and
  other suitable AIMS/IMOS temperature records. Surface SST alone cannot detect
  all bottom intrusions.
- Enter predicted freshwater and cooling relief into the mortality ensemble as
  opposite-signed DHW interactions. Pass their uncertainty forward. Retain
  aerial bleaching as an optional within-event nowcast gate. Exclude the
  habitat map from this stage because it does not resolve either dynamic
  process.

Implementation: `scripts/extract_aims_water_quality.R`,
`scripts/extract_aims_logger_hourly.R`,
`scripts/extract_bom_cyclone_exposure.R`,
`scripts/build_aims_freshwater_calibration.R`,
`scripts/fit_aims_freshwater_proxy.R`,
`data/processed/aims_freshwater_training_points.csv`,
`data/processed/freshwater_predictor_grid_2016_2025.csv`,
`data/processed/bom_cyclone_reef_year.csv`, and
`reports/aims_cyclone_upwelling_assessment.qmd`.

## 7k. Repeated thermal exposure and susceptible composition - predictor contract updated

- Following McWilliam et al. (2026), calculate the number of completed years
  exceeding 6 DHW since 2016 and the years since the most recent >6-DHW event.
  Every history stops at the year before the modelled event; the response-year
  DHW can never contribute. Retain an explicit no-prior-event flag instead of
  treating the legacy value of 35 years as an observed recovery interval. Use
  an explicit 8+ year cap in model fitting while retaining the uncapped value
  for audit.
- Also retain an eight-year rolling event count for operational forecasts after
  the fixed 2016 origin becomes increasingly remote. Compare the fixed-origin
  and rolling definitions in validation; do not fit both together by default.
- The completed all-GBR reconstruction contains 35,315 feature-years across
  7,063 feature locations. For the 2024 prediction surface, 34.5% have no prior
  >6-DHW event since 2016 and 24.7% have at least three. The target Mackay,
  Linnet, Lizard, Eyrie and Swinger features have a seven-year interval since
  their last >6-DHW event, while prior-event counts vary from one to two.
- Allow current DHW to interact with event frequency and recovery interval in
  the formal BRMS/BRT, shape-aware, shared-latent and annual-loss candidates.
  This represents selective filtering after frequent events and renewed
  susceptibility following long recovery windows.
- Keep total pre-event Acropora in all models. Add tabular and staghorn
  Acropora as separate composition fields only when source observations support
  them. The local LTMP export contains only an aggregate `Acropora` category,
  so tabular cover must not be inferred from total Acropora. Request the raw
  point-classification/taxon-growth-form data and build a past-only tabular
  composition layer with the same provenance rules as total Acropora.

Implementation: `scripts/fetch_dms_environmental_data.py`,
`scripts/add_repeated_dhw_exposure.py`,
`scripts/audit_acropora_resolution.R`, `scripts/build_validation_splits.R`,
the formal/joint/latent/shape-aware model helpers,
`scripts/fit_beta_binomial_annual_gam.R`, and
`reports/repeated_exposure_and_composition_assessment.qmd`.

## 7l. R-INLA rapid Bayesian screening - first test completed

- Use a two-component INLA approximation for rapid model development:
  Bernoulli occurrence plus conditional beta mortality magnitude, with
  independent reef, event and region effects and fold-local preprocessing.
  Keep the exact zero-inflated/zero-one-inflated beta BRMS model as the formal
  endpoint.
- The first screen completed 56 blocked fits in 4.5 minutes (median 5.2 seconds
  per two-component fit). This is fast enough to compare a small prespecified
  set of ecological hypotheses under the full blocked-validation protocol.
- The McWilliam-style event-count/recovery history improves manta tow under
  reef and event blocking and improves 2024 manta prediction. It worsens LTMP
  and is mixed for MMP. Pooled reef-blocked RMSE is essentially unchanged
  (0.142 legacy versus 0.143 new history); retain both as candidate
  specifications rather than replacing the legacy history globally.
- Both formulations continue to underpredict observed mortality >=50% by about
  26 percentage points under reef blocking. The next INLA screen must therefore
  add operational exacerbating and ameliorating modifiers before any BRMS
  refit: freshwater/coloured water or rainfall-runoff, cyclone wave exposure,
  COTS, doldrum wind exposure, and dynamic cooling/upwelling.
- Rank candidates with programme-specific blocked RMSE/MAE, predictive R2,
  occurrence Brier score, severe-case bias and calibration. Use WAIC/DIC only
  as within-likelihood full-data support, not as a replacement for blocked
  prediction. Promote only stable predictive gains to the formal BRMS model.

Implementation: `scripts/fit_inla_history_screen.R`,
`output/inla_history_screen/`, and
`reports/inla_history_screen.qmd`.

## 7m. Low-DHW high-loss sensitivity and Lizard freshwater audit - completed

- Treat low-DHW/high-loss removal only as an influence test. In every
  reef-blocked fold, remove candidate rows from training while retaining the
  complete held-out assessment set. Compare 5%, 10% and 20% mortality
  thresholds below 4 DHW with the unchanged two-part Bernoulli/beta INLA
  model.
- The mechanical hypothesis is only partly supported. Excluding rows with
  mortality >=10% steepens the mean predicted difference from <4 to 8--<12
  DHW from 13.8 to 15.8 percentage points in Manta and 15.2 to 16.4 points in
  LTMP. It does not improve extreme prediction: Manta severe-case mean
  prediction falls from 27.4% to 26.3% and RMSE worsens from 34.0 to 35.7
  points; LTMP changes negligibly.
- Do not remove these rows from the production pipeline. Mackay, Swinger and
  other 2024 cases carry information about the compound-event response that
  the model needs to learn. The reconstructed MMP field has no rows below 4,
  although four retained Stingaree rows have legacy survey DHW of 3.98 and
  reconstructed NOAA DHW of 4.15. This threshold difference does not explain
  MMP severe high-DHW underprediction.
- Record the Lizard observation gap precisely. There is no continuous AIMS
  salinity logger near Lizard in the extracted public record; RM8 is the
  nearest at about 285 km and continuous coverage ends in October 2023.
  Discrete AIMS samples do exist around Lizard, Eyrie and Linnet during the
  2024 wet season.
- Use the discrete samples as compound-event evidence, not continuous
  exposure. Lizard Island north changed from 35.3 PSU in November to 33.8 PSU
  in February and MacGillivray from 35.2 to 34.1 PSU, with concurrent
  chlorophyll increases; Eyrie measured 34.26 PSU and Linnet 33.28 PSU in
  January. These snapshots support moderate freshening but cannot rule out a
  shorter pulse below 30 PSU.
- Retain rather than delete the northern misses, and add a measured-or-proxy
  freshwater interaction with DHW. Consider a robust contamination or
  heavy-tailed observation layer for remaining unexplained non-thermal loss,
  but do not allow it to suppress compound events that are forecast targets.

Implementation: `scripts/test_low_dhw_sensitivity.R`,
`scripts/audit_lizard_freshwater_context.R`,
`output/low_dhw_sensitivity/`, and
`reports/low_dhw_sensitivity_and_lizard_audit.qmd`.

## 7n. Zero-mortality influence and modifier-dependent DHW curves - completed

- Preserve zeros in the production estimand. They comprise 56% of Manta and
  LTMP and 34% of MMP, but are event structured: 2020 contains 92% Manta and
  87% LTMP zeros, compared with 27% and 37% in 2024. They encode the genuinely
  weak 2020/2022 responses and cannot be treated as generic contamination.
- Test zero leverage by retaining 100%, 75%, 50% and 25% of zeros in each
  training event while leaving every held-out zero untouched. Report raw
  downsampling, case-control prevalence-corrected predictions, and the
  conditional-positive beta prediction.
- Raw downsampling modestly raises severe predictions but worsens overall
  reef-blocked RMSE in every programme. With 25% of zeros retained, severe
  means rise only from 27.4% to 29.2% in Manta, 38.8% to 40.9% in LTMP and
  37.2% to 39.5% in MMP. The gain disappears after prevalence correction.
- Even assuming mortality occurrence is certain, the conditional beta
  magnitude predicts only 32.4%, 43.6% and 41.8% for severe Manta, LTMP and
  MMP cases, versus observed means of 57.0%, 61.4% and 64.4%. The severe-tail
  problem is therefore mainly magnitude/heterogeneity, not the presence of
  zeros.
- Compare production BRT and Bayesian beta-plus-boundary curves using the
  three strongest cross-framework, non-redundant modifiers: pre-event
  Acropora composition, DHW novelty and cloud fraction. Show programme-specific
  10th, 50th and 90th percentile profiles and clearly fade curves outside
  observed DHW support.
- Reparameterise Acropora before formal refitting. The existing Manta Bayesian
  model contains correlated Acropora proportion and absolute-cover terms with
  opposing effects, producing a counterintuitive high-composition curve. Test
  a single interpretable composition basis plus DHW interaction rather than
  retaining both unconstrained representations.
- Correct the MMP threshold audit. Four retained 2024 Stingaree observations
  have legacy survey DHW of 3.98, but reconstructed NOAA DHW of 4.15. No rows
  were removed; the apparent absence below 4 DHW is a product/threshold
  difference. Use <4.5 only for descriptive sensitivity bins and retain
  continuous reconstructed DHW in models.

Implementation: `scripts/test_zero_mortality_sensitivity.R`,
`scripts/build_brt_beta_modifier_curves.R`,
`scripts/audit_mmp_dhw_support.R`,
`output/zero_mortality_sensitivity/`, and
`reports/zero_mortality_and_modifier_curves.qmd`.

## 7o. Coherent Acropora and extreme-aware relative mortality - completed

- Retain all zeros and target the remaining conditional-magnitude deficit.
  Fit a standard two-part BRT plus a specialist gate for mortality of at least
  30% and a robust tail-magnitude learner. All comparisons use untouched
  reef-blocked folds and reef/event-balanced training weights.
- Use Acropora proportion plus total pre-event cover as the common ecological
  basis. Absolute Acropora cover is competitive for LTMP and MMP but weaker
  for Manta; fitting proportion and derived absolute cover together is not
  consistently better and recreates correlated, opposing terms.
- Retain the conservative tail learner as an operational ensemble candidate.
  Pooled RMSE improves from 0.157 to 0.153 and mean prediction for observed
  mortality of at least 50% rises from 0.281 to 0.325. The soft gate raises the
  severe mean to 0.369 with RMSE 0.154, but increases the false-extreme rate
  among observations below 30% from 4.1% to 8.3%.
- Treat the soft gate as an event-risk prediction and the aggressive gate only
  as an upper-risk sensitivity. In 2024, the soft gate raises severe means
  from 39.8% to 52.9% for LTMP, 31.4% to 35.2% for Manta and 23.9% to 35.0%
  for MMP. MMP severe loss remains substantially underestimated.
- Preserve the beta-plus-Bernoulli INLA/BRMS framework. Matched held-out
  predictions show that INLA remains the strongest single general model for
  LTMP and MMP, while the tail BRT contributes more to Manta severe cases.
  Equal INLA/BRT blending improves LTMP and Manta but weakens MMP relative to
  INLA alone; estimate programme-aware stacking weights only inside nested CV.
- Move modifiers jointly when diagnosing the ecological response. At 11.5 DHW
  within observed support, the high Acropora, high novelty, high freshwater/
  coloured-water and low-cloud profile predicts 35.9%, 46.6% and 51.2% for
  LTMP, Manta and MMP; the conservative tail version predicts 45.4%, 48.9%
  and 54.1%. Low-risk profiles remain below 7%.
- Carry the resulting structure into the formal BRMS shared ecological model:
  Bernoulli occurrence, conditional beta magnitude, a continuous tail mixture
  or robust event scale, DHW interactions for composition, novelty,
  freshwater and cloud, and programme-specific observation layers. Require
  reef-blocked, leave-event-out and forward-event validation.

Implementation: `scripts/fit_joint_extreme_relative_screen.R`,
`scripts/combine_relative_framework_predictions.R`,
`output/joint_extreme_relative/`, and
`reports/joint_extreme_relative_screen.qmd`.

## 7p. INLA spatial and spatio-temporal bounded model - completed screen

- Use INLA as the main bounded-response inference engine. The implemented
  screen jointly fits Bernoulli occurrence and conditional beta magnitude for
  LTMP, Manta tow and MMP, with shared ecological effects and programme/layer
  intercepts.
- Carry a persistent SPDE field into the operational GBR forecast. It improves
  reef-blocked RMSE from 0.154 to 0.147 and predictive R2 from 0.248 to 0.311.
  Its posterior median spatial range is about 223 km, although uncertainty is
  wide because only 139 reefs contribute to the field.
- Keep independent event-specific SPDE fields for retrospective diagnosis and
  event updating. They reach RMSE 0.145 and predictive R2 0.332 at new reefs
  within observed events, but their latent surface is not available for a pure
  pre-survey forecast.
- Do not promote the AR1 event field. Its adjacent-event correlation is
  negative (median -0.49), and the five events are irregularly spaced. The
  simpler independent event fields predict just as well.
- Retain the persistent-field plus RW1 DHW-slope model as a sensitivity. It is
  the best leave-event-out candidate (RMSE 0.212 versus 0.223 non-spatial), but
  all event-specific slope intervals cross zero and the direct repeat-exposure
  interaction is uncertain. This is not evidence of adaptation.
- Preserve the severe-tail and BRT work. Spatial modelling raises the mean
  prediction for observed mortality of at least 50% from 0.263 to about 0.31,
  still far below the observed 0.617. In 2024, reef-blocked spatial predictions
  explain 26.5--28.2% of variation but predict severe cases near 0.33 versus
  0.623 observed.
- In the formal model, estimate programme/component calibration loadings around
  the shared ecological response, use posterior sampling for the occurrence
  times magnitude expectation, and compare a pre-survey persistent-field map
  with an optional event-updated map.
- A prediction-stack audit found and corrected constrained random effects for
  prediction-only reef-event and held-event levels. Corrected caches use the
  `v2` prefix; the original caches are retained only for provenance.

Implementation: `scripts/fit_inla_spatiotemporal_screen.R`,
`output/inla_spatiotemporal/`, and
`reports/inla_spatiotemporal_screen.qmd`.

## 7q. Restore current speed and SST-shape candidates - completed screen

- Restore eReefs current speed to the operational ecological core. Its omission
  from the first spatio-temporal screen was a predictor-list refactor gap, not a
  model-selection decision. Fit both its main effect and DHW interaction.
- Carry summer SST skewness, summer SST excess kurtosis and wet-season median
  IMOS Chl-a as one regularised feature group. Do not interpret skewness and
  kurtosis independently because their correlation is approximately -0.80.
- Hold the selected latent structure fixed while testing these additions:
  persistent SPDE field, event-varying DHW slope, joint Bernoulli occurrence
  and conditional beta magnitude, and all previously retained ecological
  modifiers.
- Current speed improves pooled reef-blocked RMSE from 0.146 to 0.142 and
  leave-event-out RMSE from 0.212 to 0.202. Adding SST shape and Chl-a improves
  these further to 0.139 and 0.198; pooled reef-blocked predictive R2 rises
  from 0.327 to 0.384.
- The full feature group improves LTMP under both validation schemes and MMP
  for new reefs within observed events. Manta is nearly unchanged under
  reef-blocking. MMP leave-event-out RMSE worsens from 0.265 to 0.282, so the
  formal model must allow programme-specific, partially pooled feature slopes.
- Retain the group for prediction despite uncertain individual coefficients.
  In the full shared screen, current speed and DHW-by-current are negative and
  excess kurtosis has the largest SST-shape coefficient, but all 95% intervals
  include zero. The predictive gain is multivariate rather than a supported
  single-driver claim.
- Current speed is complete for all 537 rows, but post-January 2024 eReefs
  values are partly reef-climatological and should be interpreted as persistent
  hydrodynamic exposure. Mackay Reef SST/Chl-a features for 2016 and 2024 are
  fold-locally imputed; the other 304 joint reef-event keys are complete.
- Promote the full group to the formal partially pooled INLA model and retain
  the nonlinear BRT as an interaction learner. Do not automatically add DHW
  interactions for SST skewness, kurtosis or Chl-a unless nested validation
  demonstrates an additional gain.
- Implement partial pooling for the five prespecified DHW interactions:
  Acropora, novelty, freshwater/coloured water, cloud and current speed. Each
  keeps a common ecological coefficient plus shrunk deviations for occurrence
  and positive magnitude in LTMP, Manta and MMP.
- Partial pooling improves pooled reef-blocked RMSE from 0.139 to 0.137,
  predictive R2 from 0.384 to 0.405, and the severe-case mean prediction from
  0.331 to 0.354. It improves reef-blocked performance in all three programmes.
  Leave-event-out RMSE is 0.201 versus 0.198 for completely shared slopes: MMP
  improves, while LTMP and Manta weaken slightly.
- Prefer the partially pooled model for within-event GBR mapping and retain the
  completely shared full-feature model as a future-event forecast sensitivity.
  Most layer deviations shrink close to zero; the clearest exception is the
  novelty response in MMP positive mortality.

Implementation: expanded candidates in
`scripts/fit_inla_spatiotemporal_screen.R`, corrected `v2` caches in
`output/inla_spatiotemporal/`, and updated figures in
`reports/inla_spatiotemporal_screen.qmd`.

## 7r. INLA severe-tail, WQC threshold and local heat audit - completed screen

- Add a severe-tail BRT on top of the selected partially pooled INLA
  expectation. The BRT separately estimates the probability of at least 30%
  mortality and robust tail magnitude; corrections are upward-only so the
  tail expert cannot flatten or reduce the calibrated INLA prediction.
- Encode the RRN WQC signal on its native 0-1 scale using raw frequency, a
  hinge above 0.50 and a rolling ten-year sum of annual excess above 0.50.
  Include DHW interactions for current and cumulative excess. Keep the
  reef-relative prior-ten-year percentile as a distinct anomaly feature.
- Under reef-blocked validation, the conservative WQC-plus-track tail raises
  the severe-case mean prediction from 35.4% to 38.6% and reduces severe RMSE
  from 0.331 to 0.305, while pooled RMSE remains approximately 0.137. The soft
  gate raises the severe mean to 41.9% and reduces severe RMSE to 0.283, but
  increases pooled RMSE to 0.141 and false extremes.
- Do not use the tail learner as an unconditional replacement for INLA. When a
  complete event is held out, all tail variants worsen pooled RMSE even though
  the WQC-plus-track variant improves severe-case error. Use the conservative
  tail as an operational ensemble candidate and the soft gate as a high-risk
  scenario.
- Retain WQC in the severe-tail candidate but do not overstate the 0.50 hinge.
  High WQC is enriched for severe mortality and raw WQC is influential in the
  gate, yet adding threshold/cumulative terms does not improve blocked RMSE
  over the otherwise identical legacy tail model.
- Treat the current BoM best-track distance, wind and wind-by-distance index as
  provisional. They resolve a Jasper proximity signal at Mackay and Snapper,
  but a better cyclone data set is expected shortly and should replace this
  isolated feature block for Jasper and Kirrily before the production fit.
- Use the public AIMS daily endpoint to compare 2024 Lizard Island and North
  Direction reef-flat/slope loggers with the same NOAA CoralTemp grid used in
  the model. Seasonal mean differences are -0.04 to +0.07 degrees C but the
  January means are +0.08 to +0.44 degrees C. The signed January offset is
  0.37-1.94 degree-weeks and the positive-only offset is 0.87-2.05
  degree-weeks. These are offset integrals rather than formal extra DHW because
  they are not yet thresholded against a local MMM. Do not apply a blanket
  uplift; calculate logger-consistent HotSpots and DHW as the next diagnostic.
  Linnet has no AIMS logger coverage after 2018.

Implementation: `scripts/fit_inla_severe_tail_brt.R`,
`scripts/build_lizard_logger_diagnostic.py`,
`output/inla_severe_tail/`, and
`reports/inla_severe_tail_assessment.qmd`.

## 7s. Explained extreme misses and independent reef-event audit - in progress

- Keep a versioned reef-event explanation registry at
  `data/curated/extreme_miss_explanations.csv`. Evidence labels are diagnostic
  annotations only and do not exclude observations or change the existing
  disturbance-attribution logic.
- Report both observation-level and reef-event-level residual rankings.
  Snapper contributes three observations among the top five raw misses; the
  reef-event table collapses these only for interpretation and leaves the fit
  unchanged.
- Preserve the current explanations and actions: Snapper 2024 is Jasper
  flooding; Penrith 2017 has an LTMP storm note and 34.85 RRN damaging-wave
  hours from Debbie; Gannett Cay 2020 is consistent with COTS; Mackay 2024 is
  missing Jasper/freshwater forcing; Linnet 2024 combines a plausible local
  heat mismatch with rainfall and possible COTS effects; Opal 2016 requires
  preceding-cover uncertainty propagation; and Tobias Spit 2024 has 2 m
  floodwater/Jasper notes.
- Add `data/gbrPredsAdj_20262408.csv` as a COTS outbreak-probability predictor
  in the severe-tail screen. Gannett Cay 2020 has outbreak probability 0.934,
  independently supported by RRN IDW COTS pressure of 0.8393. Retain both
  until blocked validation establishes whether either is redundant.
- Replace the provisional cyclone fields when the forthcoming Jasper/Kirrily
  data arrive; preserve RRN wave hours in the interim because they clearly
  identify Penrith's Debbie exposure.

Implementation: `data/curated/extreme_miss_explanations.csv`, the monthly
logger-NOAA diagnostic, the COTS hindcast extension and reef-event prediction
output in `scripts/fit_inla_severe_tail_brt.R`, and the corresponding tables in
`reports/inla_severe_tail_assessment.qmd`.

## 7t. Cyclone/COTS INLA baseline, local DHW uplift and disease risk - completed screen

- Reconstruct local DHW from the AIMS logger temperatures using NOAA CRW v3.1
  MMM and the standard 84-day accumulation of HotSpots at least 1 degree C.
  The reconstructed NOAA peak agrees with the official daily product within
  0.18 DHW, supporting the implementation.
- Across November--April, logger peak DHW exceeds official NOAA by 4.26 and
  3.44 DHW at Lizard Island and by 5.38 and 1.05 DHW at North Direction.
  Carry a regional median sensitivity of +3.85 DHW, with depth-stratified
  means of +4.82 for reef flats and +2.25 for slopes. For the explicitly local
  2024 sensitivity, apply the flat mean to all six nearby reef exposures
  regardless of survey depth: NOAA is a surface-temperature product, so this
  test corrects the satellite exposure rather than the benthic survey depth.
  Do not extrapolate the fixed uplift beyond the named cluster.
- Add provisional cyclone-track proximity, a wind-distance index and the GBR
  COTS outbreak-probability hindcast to the central INLA candidate. This lowers
  reef-blocked RMSE from 0.1369 to 0.1357 and leave-event-out RMSE from 0.2008
  to 0.1973, while raising the severe-case mean from 35.4% to 35.9%. Promote
  this candidate as the INLA baseline for the severe-tail ensemble.
- Extract NOAA CRW disease outbreak risk over each November--April event
  window. Preserve negative product mask codes as structural non-applicability,
  not negative or missing risk. The product is applicable at 23 validation
  reef-years in 2020 and three in 2022, but none in 2016, 2017 or 2024.
- Disease intensity and duration do not improve the current models. Adding
  them to INLA worsens leave-event-out RMSE from 0.1973 to 0.2114 and raises
  false extremes; adding them to the tail BRT also weakens blocked prediction.
  Retain the metric for ecological attribution and sensitivity analysis, not
  the operational ensemble.
- Keep the cyclone source block replaceable. The forthcoming Jasper/Kirrily
  data should be evaluated with the same frozen validation splits and compared
  against both the provisional track fields and RRN damaging-wave hours.

Implementation: `scripts/build_lizard_logger_diagnostic.py`,
`scripts/extract_noaa_disease_risk.py`,
`tests/test_noaa_disease_risk.py`,
`scripts/fit_inla_spatiotemporal_screen.R`,
`scripts/fit_inla_severe_tail_brt.R`,
`data/processed/lizard_cluster_logger_dhw_uplift_2024.csv`,
`data/processed/noaa_disease_risk_validation.csv`, and
`reports/inla_severe_tail_assessment.qmd`.

## 7u. GBR-wide AIMS logger validation and local 2024 DHW sensitivity - completed screen

- Extract all public AIMS daily temperature-logger records for the 2016, 2017,
  2020, 2022 and 2024 bleaching events and collocate each series with NOAA
  CoralTemp SST, official DHW and the NOAA MMM. Reconstruct logger DHW with the
  same 84-day, HotSpot-at-least-1-degree definition used by NOAA.
- Accept a logger-event DHW comparison only when at least 70 daily
  observations occur in the 84-day window at its eligible peak. This yields
  143, 186, 202, 163 and 148 coverage-qualified GBR series-events respectively.
- Confirm a strong 2016 Lizard-region discrepancy. Linnet reaches 15.39 DHW
  versus NOAA 7.76 (+7.63), while Martin reaches 11.59 versus NOAA 7.87
  (+3.72). Both have complete summer records. The logger API has no 2016 North
  Direction series; its available record begins during the 2017 event.
- Preserve the discrepancy as event-specific and spatially heterogeneous.
  The GBR medians are -0.16, -0.22, -0.25, -0.47 and -1.50 DHW across the five
  events, while individual reefs can differ by several DHW in either
  direction. A universal correction would therefore be poorly calibrated.
- Test a transparent 2024 local sensitivity using +4.82 DHW at MacGillivray,
  Lizard NW, Eyrie, Martin, Linnet and North Direction. With the same INLA
  structure, reef-blocked GBR-wide RMSE improves from 0.1357 to 0.1349 and
  leave-event-out RMSE from 0.1973 to 0.1944. At the six target reefs,
  reef-blocked RMSE improves from 0.178 to 0.101 and mean bias from -0.147 to
  -0.010; when all 2024 observations are excluded from training, local RMSE
  improves from 0.349 to 0.104. The correction improves Eyrie and Linnet substantially but
  overpredicts Martin, so retain it as a sensitivity rather than production
  truth.
- Refit the severe-tail BRT on the adjusted INLA baseline. The conservative
  WQC/cyclone-track ensemble gives reef-blocked GBR-wide RMSE 0.134 versus
  0.136 with the unadjusted baseline, while the fully soft gate still adds too
  much positive bias. Keep the conservative gate as the operational candidate.
- Next replace the fixed local uplift with a spatial, event-specific NOAA-bias
  surface trained on coverage-qualified logger discrepancies. Validate it by
  leaving logger sites and entire events out, and propagate correction
  uncertainty into the mortality model. The Lizard automated weather station
  lists additional 2016 temperature series, but its measurement endpoint
  requires an AIMS API key and remains a documented follow-up.

Implementation: `scripts/build_aims_noaa_dhw_validation.py`, cached event
downloads in `data/raw/aims_temperature/events/`,
`data/processed/aims_logger_noaa_dhw_validation.csv`, the adjusted candidates
in `scripts/fit_inla_spatiotemporal_screen.R` and
`scripts/fit_inla_severe_tail_brt.R`, and
`reports/aims_noaa_dhw_validation.qmd`.

## 7v. Authenticated AIMS extraction and event-specific DHW correction - completed screen

- Store the AIMS Data Platform credential only in the Windows user environment
  as `AIMS_DATAPLATFORM_API_KEY`. The extraction script reads that variable at
  runtime; the credential is not written to code, data products or reports.
- Use `scripts/pull_aims_temperature_data.R` to retrieve the complete
  temperature deployment catalogue and all daily logger observations in the
  prespecified 2016, 2017, 2020, 2022 and 2024 event windows. Retain deployment
  depth and optionally add automated-station water temperatures. Raw responses
  are cached so routine rebuilds do not repeatedly download high-frequency
  data. Environment settings control years, refreshes and automated-site scope.
- Add the authenticated 2016 Lizard Island automated records. The shallow
  0.6 m and deep 10.1 m series have logger-minus-NOAA discrepancies of -0.68
  and -0.35 DHW. This confirms that the strong positive 2016 discrepancies at
  Linnet (+7.63) and Martin (+3.72) were spatially local rather than a uniform
  Lizard Island offset.
- Retain depth explicitly, but do not yet estimate a depth slope. Only five
  coverage-qualified site-events contain paired shallow and deep observations,
  with a median shallow-minus-deep discrepancy of 0.50 DHW and inconsistent
  event contrasts. Prefer
  shallow/flat records (0--5 m) for NOAA surface calibration and use deeper
  records where shallow records are unavailable; preserve the paired depth
  audit for future expansion.
- Replace the fixed local uplift with two partially pooled spatial products in
  `scripts/build_event_dhw_correction_layer.py`:
  (1) a historical-only rolling-prior surface for prospective prediction, and
  (2) a within-event update for retrospective or live-event recalibration.
  Both expose distance, effective logger count and information-availability
  fields. The selected ranges are 100 km for the historical surface and 200 km
  for the broad within-event offset.
- In blocked logger validation, the historical surface improves event-held-out
  discrepancy RMSE from 4.93 DHW with no correction to 4.73 DHW at the selected
  100 km range (best screened event-held-out RMSE 4.67 DHW at 25 km). Large
  individual logger anomalies dominate the remaining error, so correction
  uncertainty must be propagated rather than treating the surface as truth.
- Add historical and within-event correction candidates to the same INLA
  mortality structure. The historical surface slightly improves
  leave-one-event-out mortality RMSE from 0.1973 to 0.1942 and severe RMSE from
  0.4069 to 0.3985. Reef-blocked RMSE is effectively unchanged (0.1357 versus
  0.1368). The within-event surface raises exposure at the 2024 Lizard-region
  misses but gives leave-event-out RMSE 0.2000; it is not promoted as the
  prospective default.
- Preserve the ecological interpretation: corrected 2024 DHW rises by about
  0.96 at Linnet, 0.94 at Martin and 0.84 at Eyrie, yet severe mortality remains
  strongly underpredicted. Local NOAA mismatch explains part of the misses but
  cannot replace freshwater/WQC, Acropora composition, cyclone or COTS
  modifiers. Keep the former fixed +4.82 DHW uplift as a stress-test only.

Implementation: `scripts/pull_aims_temperature_data.R`,
`scripts/build_aims_noaa_dhw_validation.py`,
`scripts/build_event_dhw_correction_layer.py`, the rolling and within-event
candidates in `scripts/fit_inla_spatiotemporal_screen.R`,
`data/processed/aims_temperature_daily_with_depth.csv.gz`,
`data/processed/noaa_dhw_correction_layer_validation.csv`, outputs under
`output/dhw_correction/`, and
`reports/event_dhw_correction_assessment.qmd`.

## 7w. Local-first, peak-qualified DHW calibration - completed screen

- Replace the broad within-event surface with a measurement hierarchy. Use a
  coverage-qualified reef logger before an automated-station record at the same
  site-event, prefer shallow/flat reef series, and require matching normalised
  names plus separation below 5 km for a direct reef match. Otherwise
  interpolate from no more than four current-event sites inside a hard 25 km
  radius. Where current local support is absent, retain NOAA unchanged.
- Qualify records at the actual unconstrained logger-DHW peak. The previous
  rule could select an earlier peak merely because it had 70 observations in
  its 84-day window. In 2024 the Lizard relay-pole 5 m record consequently
  appeared to differ from NOAA by -2.37 DHW, although its later actual peak
  implied +2.18 DHW with only 66 days. The new rule excludes this incomplete
  record and retains the complete reef-flat comparison of +4.26 DHW.
- Do not combine a local measurement with a GBR event offset or historical
  surface. The previous approximately +0.55 to +0.96 DHW Lizard-cluster
  corrections resulted from pooling incompatible site records, smoothing over
  200 km, then adding negative event and historical terms. The revised local
  estimates are +0.33 at MacGillivray, +4.26 at Lizard NW, +3.52 at Eyrie,
  +3.84 at Martin, +3.62 at Linnet and +5.38 at North Direction. The first,
  second, fourth and sixth values are direct matched measurements; Eyrie and
  Linnet are locally interpolated.
- Retain ENSO as event context, not a GBR-wide correction. Event-median
  logger-minus-NOAA discrepancies are +0.34 in 2016, +0.24 in 2017, +0.29 in
  2020, -0.03 in 2022 and -0.59 DHW in 2024. These do not support a simple
  phase offset. Allow an ENSO-matched historical fallback only after at least
  two previous same-phase events at a site and at least 75% agreement in sign;
  none of the current events meets that prospective evidence rule.
- Select 25 km despite the marginally lower overall held-out discrepancy RMSE
  at 50 km: 25 km gives lower MAE (1.59 versus 1.60 DHW), lower extreme RMSE
  (4.06 versus 4.19 DHW), and better preserves the intended local estimand.
- Use the same local calibration hierarchy for historical fitting rows and the
  held-out/current event when those logger records are available. This primary
  leave-one-event-out candidate reduces INLA mortality RMSE from 0.1973 to
  0.1720, changes predictive R2 from -0.236 to 0.061, reduces severe-case RMSE
  from 0.4069 to 0.3562, and lowers the false-extreme rate from 0.113 to 0.087.
  Retain a current-event-only update as a severe-tail sensitivity: it gives
  RMSE 0.1827 and severe RMSE 0.3195, but raises the false-extreme rate to
  0.163. Local calibration remains incomplete without freshwater/WQC,
  composition, cyclone and COTS modifiers.
- Preserve the original NOAA exposure and store correction, source, distance,
  effective logger count and local variability in separate versioned fields.
  Propagate correction uncertainty in the final INLA model rather than treating
  the mean local update as error-free.

Implementation: `scripts/build_local_first_dhw_correction.py`, local-first and
operational-update candidates in `scripts/fit_inla_spatiotemporal_screen.R`,
`data/processed/noaa_dhw_correction_layer_local_first_validation.csv`, outputs
under `output/dhw_correction_local_first/`, tests in
`tests/test_local_first_dhw_correction.py`, and
`reports/local_first_dhw_correction_assessment.qmd`.

## 7x. Formal fully local-calibrated INLA model - active production stage

### Locked central-estimate contract

- Promote `persistent_rw1_cyclone_cots_local_first_dhw_partial_pool` from a
  sensitivity candidate to the central local-temperature model. Keep NOAA DHW
  immutable, retain the correction as a separate field, recompute all derived
  DHW hinges and novelty terms after correction, and apply the same hierarchy
  to fitting and prediction rows.
- This is a multi-event calibration, not a 2024 override. Direct/local support
  exists for 22 of 29 reef-events in 2016, 37 of 38 in 2017, 56 of 73 in 2020,
  43 of 61 in 2022 and 83 of 103 in 2024. Unsupported reefs retain NOAA as the
  point estimate; they do not receive a GBR event offset or inconsistent
  historical surface.
- Preserve the direct-measurement precedence rule: coverage-qualified reef
  logger, then automated station only if no reef logger exists, direct survey
  match only with normalised name agreement and distance below 5 km, otherwise
  no more than four sources inside 25 km. Record source, distance, effective
  source count, correction SD and uncertainty method in every prediction row.

### Validation interpretation

- Pooled leave-one-event-out RMSE improves from 0.1973 to 0.1720 and predictive
  R2 from -0.236 to 0.061. The gain is not temporally uniform: event-level RMSE
  improves by 0.011 in 2020, 0.014 in 2022 and 0.059 in 2024, but worsens by
  0.020 in 2016 and 0.036 in 2017.
- Do not tune the correction toward mortality to remove those early-event
  errors. Penrith is a labelled Cyclone Debbie loss, and Carter/Yonge show high
  mortality despite local evidence that NOAA overestimated heat. These rows
  indicate unresolved non-thermal processes rather than invalid logger
  measurements.
- Retain the current-event-only candidate as a severe-tail ensemble input. Its
  severe RMSE is 0.319 versus 0.356 for the central fully calibrated model, but
  its higher false-extreme rate precludes using it as the central estimate.

### Operational uncertainty - implemented first stage

- Estimate DHW-correction uncertainty without mortality outcomes. Direct
  corrections use replicate logger spread with a pooled 0.61 DHW floor. Local
  interpolation combines the weighted source spread with a 1.59 DHW held-out
  floor. Unsupported NOAA values retain zero mean correction with 3.47 DHW
  uncertainty from omitted-site validation.
- Translate DHW error to mortality with a non-negative, programme-specific
  first-order response slope estimated from the unchanged-training/current-
  event-update contrast. Combine it with the maximum programme-specific 90%
  conformal residual radius among other events. The resulting intervals attain
  93.9% pooled coverage for a nominal 90%, with mean total width 0.526 and mean
  DHW-only width 0.068. Manta 2024 coverage remains lower, so expose event and
  programme calibration alongside every operational map.
- Treat this as the operational interval implementation, not the final Bayesian
  measurement-error fit. The formal endpoint is 20--50 repeated, spatially
  correlated correction-layer imputations with INLA posterior predictions
  pooled using within- and between-imputation variance. This is required
  because corrected DHW enters nonlinear hinges, novelty and multiple
  interactions; a single linear measurement-error coefficient is insufficient.

### Updated miss register and ordered next improvements

1. **Implemented: separate thermal, cyclone and COTS pressure blocks.** The
   central INLA candidate is now
   `persistent_rw1_local_dhw_decomposed_hazards_freshwater_partial_pool`.
   It retains every observation, uses shared thermal effects, and adds distinct
   interval-aligned cyclone-wave/track, COTS-surface/hindcast and freshwater
   terms. This is an additive predictive decomposition, not a hard causal
   assignment from disturbance notes. Carter and Yonge exposed a material join
   error: their 2016 cover intervals span Cyclone Nathan, but the event-summer
   join gave zero wave hours. The new interval calculation assigns 68.72 h at
   Carter and 64.44 h at Yonge; Agincourt's baseline is after Nathan and stays
   unexposed.
2. **Implemented: continuous DHW-by-freshwater/coloured-water screen.** The
   new central candidate replaces the single freshwater composite with DHW
   interactions for coastal rainfall, current WQC frequency, reef-relative WQC
   percentile and the unthresholded rolling 10-year WQC sum. The previous
   `wqc_excess50_10yr_sum` is retained for diagnostics, but is not imposed as a
   universal threshold in the selected model. This candidate improves pooled
   leave-one-event-out RMSE from 0.172 to 0.158, predictive R2 from 0.060 to
   0.209, severe RMSE from 0.356 to 0.346 and false-extreme rate from 0.143 to
   0.060. It retains 2024 severe-tail performance (0.316 versus 0.310 severe
   RMSE for the prior local model).
3. **COTS probability--intensity combination tested; competing hazards remain
   required.** The hindcast represents the probability of exceeding 0.22
   COTS/tow, whereas the RRN interval maximum supplies continuous intensity.
   A candidate retains outbreak probability and replaces raw log intensity
   with $p\log[1+\max(I-0.22,0)]$. It improves event-held-out RMSE from 0.1576
   to 0.1556, predictive R2 from 0.211 to 0.231, severe RMSE from 0.3454 to
   0.3418 and false-extreme rate from 0.0566 to 0.0523, while reef-blocked RMSE
   is effectively unchanged. Gannett Cay 2020 has 9.91 COTS/tow and hindcast
   probability 0.934, but its held-out prediction remains only 0.055 versus
   0.566 observed. Keep the combined feature as a promising candidate, not a
   completed attribution. The next model must use labelled disturbance
   evidence to estimate competing thermal and COTS hazards and audit
   high-pressure/low-loss rows for interval timing and cover support.
4. **Operational uncertainty updated.** The matching current-event-local
   candidate supplies DHW sensitivity for the selected structure. The updated
   90% interval has 94.6% pooled coverage and mean width 0.495; the DHW-only
   component averages 0.115. Keep the repeated correction-layer imputation as
   the formal measurement-error endpoint.
5. **Replace the current cyclone proxy when the Jasper/Kirrily-resolving layer
   arrives.** Keep track distance, intensity, damaging-wave exposure and
   rainfall/runoff as separate mechanisms and repeat the unchanged blocked
   comparison.
6. **Use the current-event-only INLA and severe-tail BRT as constrained ensemble
   members.** Optimise weights against event-blocked RMSE, severe RMSE and
   false-extreme rate rather than allowing the severe-tail component to
   dominate the mean prediction.
7. **Resolve the 2016--2017 northern residuals.** Review survey notes, cover
   uncertainty, cyclone/COTS evidence and temperature records for Carter,
   Yonge, Agincourt and St Crispin. Their local corrections mostly reduce DHW,
   so another heat uplift is not supported.
8. **Refine composition only with complementary structure.** Retain preceding
   Acropora; add tabular Acropora and weedy-recovery assemblage information
   where attainable instead of another correlated total-cover measure.

Implementation: uncertainty fields in
`scripts/build_local_first_dhw_correction.py`, uncertainty provenance carried
through `scripts/fit_inla_spatiotemporal_screen.R`, operational propagation in
`scripts/propagate_local_dhw_uncertainty.py`, outputs under
`output/local_calibrated_operational/`, regression tests in
`tests/test_local_dhw_uncertainty.py`, and updated tables and figures in
`reports/local_first_dhw_correction_assessment.qmd`.

## 7v. Cause-aware COTS and cyclone competing hazards - selected

- Refit the thermal/freshwater INLA response after excluding only rows with
  explicit COTS, cyclone, storm or flood attribution. Unlabelled mortality is
  retained, preserving the restricted bleaching-mortality logic while stopping
  strongly supported non-thermal losses from flattening the thermal curve.
- Fit separate Bernoulli-beta COTS and cyclone hazards to the longer annual
  coral-cover transition record. COTS uses hindcast outbreak probability plus
  probability-weighted RRN excess density above 0.22 COTS/tow. Cyclone uses
  damaging-wave hours, wind-distance exposure and rainfall. Event-overlapping
  annual intervals and held-out reefs are excluded within their respective
  validation folds.
- Combine hazards on the response scale as
  $1-(1-T)(1-C)(1-S)$. The operational version activates COTS only above 0.22
  COTS/tow and cyclone only above 20 hours of waves exceeding 4 m. The ungated
  version improves the severe tail but creates too much unsupported background
  mortality and is rejected.
- Promote `cause_aware_exposure_gated_competing_hazards`. Against the prior
  selected model, leave-one-event-out RMSE improves from 0.1576 to 0.1507,
  predictive R2 from 0.211 to 0.279, severe RMSE from 0.345 to 0.326 and the
  false-extreme rate from 0.0566 to 0.0436. Reef-blocked RMSE is nearly
  retained (0.1398 versus 0.1390) and severe RMSE improves (0.284 versus
  0.308).
- Treat Gannett Cay and Penrith as cause-supported but magnitude-unresolved
  positive controls. Gannett's held-out prediction rises from 0.060 to 0.228;
  Penrith's rises from 0.077 to 0.206. Penrith joins Double Cone, Daydream and
  Shute as a Debbie case with more than 20 damaging-wave hours. Chinaman,
  Taylor and Rib Reef place their main COTS losses in the 2017-18 transitions;
  high-pressure/low-loss 2020 rows remain as timing and depleted-cover negative
  controls.
- Retest the 20-hour cyclone gate and conditional severity as soon as the
  improved Jasper/Kirrily dataset arrives. Next refine COTS with time since
  peak pressure and available starting cover, then rerun the major-miss table.

Implementation: `scripts/fit_cause_aware_competing_hazards.R`, outputs under
`output/cause_aware_competing_hazards/`, model selection in
`config/model_registry.yml`, and the complete assessment in
`reports/cause_aware_competing_hazard_assessment.qmd`.

## 7w. COTS timing, soft cyclone activation and La Nina residuals - completed

- Test event-year-available COTS timing and cover support using preceding
  three-year pressure, years since the prior outbreak, starting cover and
  pressure-by-cover interactions. Do not promote this block: event-held-out
  RMSE is 0.1508 versus 0.1507 for the simpler COTS model, and the two
  pressure-by-cover features correlate above 0.999 with their parent pressure
  terms. Starting cover remains in the COTS occurrence and magnitude models;
  develop the formal biomass constraint in the absolute-cover-change model.
- Replace the hard 20-hour cyclone switch with a logistic activation
  $\operatorname{logit}^{-1}[(W-20)/5]$. Promote
  `baseline_cots_logistic20_5_cyclone`: event-held-out RMSE improves from
  0.1507 to 0.1492, predictive R2 from 0.279 to 0.293, severe RMSE from 0.3259
  to 0.3254 and false-extreme rate stays 0.0436. Reef-blocked RMSE improves
  marginally from 0.13983 to 0.13979. Retest the activation and propagate its
  uncertainty when the replacement cyclone layer arrives.
- Audit the negative 2020/2022 predictive R2 values. These events have very low
  observed means (0.016 and 0.027) and many zeros, while model means are 0.064
  and 0.059. The problem is primarily a positive occurrence floor plus a small
  number of opposing outliers, not universally high magnitude error. Test
  ENSO/SOI and rainfall/runoff anomalies in the Bernoulli occurrence component
  rather than flattening the conditional DHW response.
- Retain Gannett, Penrith and Mackay as supported but unresolved focal cases.
  Add high-DHW zero-loss reefs in 2020 and low-DHW zero-loss reefs in 2022 as
  counter-controls for composition, cooling and occurrence calibration.
- Use `reports/_model_next_steps.md` as the canonical next-step list included
  by every current model report, preventing priorities from drifting between
  iterations.

Implementation: `scripts/test_cots_timing_cover_and_cyclone_soft_gate.R`,
outputs under `output/cots_timing_cover_cyclone_soft_gate/`, figures
`Fig-INLA-13` to `Fig-INLA-15` and `Fig-DATA-03`, and assessment report
`reports/cots_timing_lanina_residual_assessment.qmd`.

## 7x. Raw COTS intensity, La Nina occurrence state and sequential update - completed

- Compare the selected log COTS interval excess with untransformed interval
  density, a prospective current-year density plus time-since-peak state, and
  the same compact state fitted to absolute percentage-point cover loss. The
  compact predictors are complementary (largest absolute off-diagonal
  Spearman correlation below 0.42), so rejection is based on transfer rather
  than redundancy.
- Promote `cots_raw_interval_logistic20_5_cyclone`. Raw density improves
  event-held-out RMSE from 0.14922 to 0.14891 and predictive R2 from 0.2928 to
  0.2957; reef-blocked RMSE improves from 0.13979 to 0.13933. Severe RMSE
  weakens slightly from 0.3254 to 0.3275 and Gannett 2020 is predicted less
  well, so retain the log form as a severe-tail sensitivity and require a
  Gannett-positive-control check for any future COTS tail model.
- Do not promote current-year density plus time since peak or the
  absolute-cover-loss component: both worsen event transfer. The selected
  interval maximum remains retrospectively timed, so reconstruct an
  event-start COTS pressure nowcast from the latest survey, time since
  survey/peak, culling and hindcast probability before prospective mapping.
- Do not promote the RONI or SOI occurrence recalibrations. Coupled indices
  with reef-centred ERA5 rainfall and WQC anomaly improve 2022 but worsen 2020
  and suppress 2024. Five events and a coloured-water proxy cannot identify a
  transferable low-mortality state; true runoff/discharge remains a data gap.
- Retain the first-20%-of-reefs prevalence model only as a sequential update.
  On later reefs RMSE improves from 0.1554 to 0.1533 and severe RMSE from
  0.3655 to 0.3234. Next impose spatial separation and substitute aerial or
  rapid in-water prevalence where possible.
- Audit named controls with local-first logger support and survey context.
  Gannett remains a severe COTS magnitude miss; Havannah has positive direct
  logger correction; Linnet has programme disagreement; U/N 20-104 lacks a
  supported local correction; Pandora remains a high-DHW zero-loss control
  despite positive direct correction; and the named 2022 zero-loss reefs
  remain occurrence-floor controls.

Implementation: `scripts/test_cots_raw_enso_occurrence_controls.R`, outputs
under `output/cots_raw_enso_occurrence/`, figures `Fig-INLA-16` to
`Fig-INLA-18`, model registry version 4, and assessment report
`reports/04_COTS_ENSO_reef_control_assessment.qmd`.

## 7y. Spatially independent aerial/RHIS early-event update - completed, gated

- Convert RHIS morphology-specific bleaching to the requested 0--4 scale and
  calculate a cover- and percent-bleached-weighted community burden. Retain
  recently dead benthos only where bleaching is explicitly present. Use April
  30 as the primary common cutoff and March 31 as the strict timing
  sensitivity; never use later post-event mortality observations or cause
  labels.
- Treat Hughes aerial `bin.score` as a binary severe-bleaching event-time
  indicator for 2016, 2017 and 2020. Record the absence of survey dates in the
  supplied CSV as a blocking provenance limitation rather than assuming exact
  timing. Add 2022/2024 only when their scores and survey dates are supplied.
- Validate twice: fit update coefficients with the target event held out, and
  construct its spatial signal only after excluding the target sector or its
  full latitude block. Keep the selected initial forecast and its conditional
  magnitude fixed in the primary occurrence-only comparison.
- Do not promote the current candidates. RHIS-only updates worsen five-event
  RMSE and Brier score; aerial-only updates also fail. Aerial plus RHIS lowers
  RMSE from 0.1380 to 0.1338 and raises predictive R2 from 0.228 to 0.275 on
  the three available aerial events under sector exclusion, but severe RMSE
  worsens from 0.381 to 0.409. A two-part occurrence/magnitude update is worse.
- Avoid a larger early-data predictor block: sector-independent RHIS severity
  versus burden has Spearman correlation 0.912 and burden versus recently dead
  is 0.815. The next candidate is aerial severity plus one RHIS burden term.
- Keep the update disabled by default and separate from the initial forecast.
  Promotion requires improvement in overall RMSE, occurrence Brier score and
  severe RMSE under both event-held-out and spatially independent validation.

Implementation: `scripts/test_spatial_early_bleaching_update.R`, processed
sources `data/processed/rhis_early_bleaching_reef_event.csv` and
`data/processed/aerial_early_bleaching_reef_event.csv`, outputs under
`output/spatial_early_bleaching_update/`, figures `Fig-NOWCAST-01` to
`Fig-NOWCAST-03`, model registry version 5, and assessment report
`reports/05_spatial_early_bleaching_update_assessment.qmd`.

## 7z. Prospective COTS nowcast and cross-fitted residual BRT - completed

- Treat undated RRN event-season density as pre-bleaching pressure, following
  the data-owner interpretation. Freeze exact-date Manta and culling evidence
  at the end of February. Construct reef-event features from raw/log RRN
  excess, latest Manta density, the prior three-year Manta peak, time since
  that peak and the event-year outbreak hindcast probability.
- Audit the targeted Cull table over the preceding 365 days using total
  removals, dive effort, removals per dive and positive-dive fraction. Do not
  interpret removals as spatially standardised density. Culling terms slightly
  worsen both event-held-out and reef-blocked transfer and are retained only as
  ecological context.
- Promote `operational_rrn_raw_plus_manta_state`. Relative to the prior raw
  interval reference, leave-one-event-out RMSE improves from 0.14891 to
  0.14848 and predictive R-squared from 0.296 to 0.300; severe RMSE improves
  from 0.3275 to 0.3267. Reef-blocked RMSE is 0.13958. The raw RRN signal lifts
  Gannett's held-out prediction from the Manta-only 0.121 to 0.183. The log
  hybrid lifts it further but slightly worsens aggregate event transfer.
- Fit balanced and four-times severe-weighted residual BRTs inside nested
  leave-one-event-out and reef-blocked folds. Learn the INLA--BRT blend weight
  inside each training fold. Neither BRT passes promotion: balanced BRT
  event-held-out RMSE is 0.14923 and severe RMSE 0.33218; the severe-weighted
  version is 0.15019 and 0.33585. Both improve reef-blocked severe RMSE, which
  identifies within-event spatial structure rather than transferable
  new-event skill.
- Retain the residual BRT importance and partial-dependence plots as diagnostic
  evidence. Rainfall, the base INLA prediction, starting cover and SST
  skewness are its strongest signals, but they do not yet yield a safe
  prospective correction.

Implementation: `scripts/test_prospective_cots_nowcast.R`,
`scripts/test_inla_brt_residual_ensemble.R`, outputs under
`output/prospective_cots_nowcast/` and
`output/inla_brt_residual_ensemble/`, model registry version 6, and assessment
report `reports/06_prospective_cots_and_inla_brt_assessment.qmd`.

## 7aa. Independent full BRTs and environmental-envelope ensemble - completed

- Correct the scope of the machine-learning comparison. The preceding residual
  BRT estimated only errors left by INLA; it was not the independently fitted
  all-predictor BRT requested for model comparison. Fit two independent BRT
  candidates from the same operational predictor matrix: a direct bounded
  mortality BRT and a two-part occurrence times positive-magnitude BRT.
- Use the selected local-first DHW, composition and starting cover, thermal
  history and novelty, cloud/current/Secchi/SST shape, rainfall and coloured
  water, wind, cyclone, COTS, disease, depth, programme and region predictors.
  Apply the fixed outcome-free Spearman screen before fitting. Disease-risk
  applicability is removed because its absolute correlation with disease risk
  is 0.999; SST kurtosis and median chlorophyll remain excluded by the earlier
  complementarity decision in favour of SST skewness and Secchi.
- Tune interaction depth inside every outer fold. Evaluate both
  leave-one-event-out and five-fold reef-blocked transfer. Do not use event
  identity as an operational predictor.
- Define a fold-local environmental envelope from training predictors only:
  training-median imputation, standardisation, PCA retaining 90% variance
  (maximum eight components), whitening and ten-neighbour distance. Test hard,
  smooth 50% and full-applicability INLA--BRT weights. No observed mortality or
  residual enters the gate.
- Do not promote any BRT or envelope blend to the initial forecast. The direct
  BRT improves reef-blocked RMSE from 0.13958 to 0.13540, and the direct
  applicability blend reaches 0.12817, but leave-one-event-out RMSE worsens
  from 0.14848 for INLA to 0.17994 for the direct BRT. Severe RMSE worsens from
  0.32666 to 0.46502. The smooth direct blend is closer but still worse
  (RMSE 0.15079; severe RMSE 0.35539).
- Reject the proposed rule that BRT should dominate merely because a reef is
  inside the training environmental envelope. In event-held-out rows inside
  the envelope, INLA RMSE is 0.13541 versus 0.16916 for the direct BRT. This
  shows that represented covariate space does not guarantee representation of
  the event-level mortality process.
- Retain the direct BRT for nonlinear effect inspection, within-event
  interpolation after the event state is represented, and future spatially
  independent early-event updating. Preserve the two-part BRT as a likelihood
  sensitivity; it better controls low-event mortality in some years but
  severely underpredicts 2024 magnitude.

Implementation: `scripts/fit_standalone_brt_envelope_ensemble.R`, outputs
under `output/standalone_brt_envelope_ensemble/`, figures `Fig-FULLBRT-01` to
`Fig-FULLBRT-08`, model registry version 7, and assessment report
`reports/07_standalone_brt_envelope_ensemble_assessment.qmd`.

## 8. Refactor the reports into a research-grade narrative

### Addendum: SST distribution shape and median chlorophyll screen

- NOAA CoralTemp v3.1 now supplies acute-summer SST skewness and excess
  kurtosis, while IMOS MODIS-Aqua OC3 supplies wet-season median chlorophyll.
  The validation extraction has complete summer SST and wet-season chlorophyll
  coverage for all 304 reef-years.
- Q1 and wet-season chlorophyll are near duplicates (Spearman rho 0.97), so
  only wet-season median chlorophyll is retained as a direct candidate.
- The first fixed-settings blocked BRT/INLA screen finds no universal upgrade.
  Raw SST-shape/chlorophyll improves LTMP BRT reef-blocked RMSE from 0.133 to
  0.127 and 2024 RMSE from 0.121 to 0.107. Modifier PCA improves MMP BRT
  reef-blocked RMSE from 0.266 to 0.250. These gains do not remain stable under
  event-held-out validation, and manta's best INLA result remains the core
  model.
- Keep the SST shape and wet-season chlorophyll variables as grouped candidates.
  Test them jointly with freshwater, coloured water, cyclone, wind, COTS and
  cooling-relief variables under nested tuning before promoting them to the
  production BRMS/BRT ensemble.

- Keep exploratory diagnostics separate from confirmatory modelling and final
  evaluation.
- Generate one canonical table/figure per claim from reusable functions;
  remove duplicated model fits, stale cached results, and contradictory labels.
- Lead with the most stable, out-of-sample findings and clearly distinguish
  association, prediction, and mechanistic interpretation.
