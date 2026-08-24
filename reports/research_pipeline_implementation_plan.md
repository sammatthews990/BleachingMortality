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

## 8. Refactor the reports into a research-grade narrative

- Keep exploratory diagnostics separate from confirmatory modelling and final
  evaluation.
- Generate one canonical table/figure per claim from reusable functions;
  remove duplicated model fits, stale cached results, and contradictory labels.
- Lead with the most stable, out-of-sample findings and clearly distinguish
  association, prediction, and mechanistic interpretation.
