# Figure registry

## Fig-INLA-05_hinge_validation — model_comparison

Caption: Held-out model performance comparing thermal spline hinges at 4+8 and 6+8 DHW.
Interpretation: Tests whether moving the first change in thermal slope from 4 to 6 DHW improves transfer to unseen events and reefs.
Caveat: The alternative also moves the onset of DHW interactions with Acropora, novelty, freshwater and cloud; it is therefore a whole-response sensitivity test.
Model: dhw_hinge_sensitivity (INLA).
File: output/fig/Fig-INLA-05_hinge_validation_v001.png

## Fig-INLA-06_hinge_response — hinge_response

Caption: Binned held-out observed and predicted mortality for the 4+8 and 6+8 DHW spline candidates.
Interpretation: Shows where moving the first hinge changes practical predictions across the observed DHW range.
Caveat: This is a binned validation diagnostic rather than a covariate-adjusted posterior dose-response curve.
Model: dhw_hinge_sensitivity (INLA).
File: output/fig/Fig-INLA-06_hinge_response_v001.png

## Fig-INLA-07_spatial_residual_fields — spatial_residual_field

Caption: Posterior mean event-specific spatial residual fields for 2016, 2017, 2020, 2022 and 2024, with the selected model persistent field shown as the overall panel.
Interpretation: Maps residual geographic structure after measured covariates. Repeated red or blue areas suggest unresolved persistent mechanisms; event-only patterns suggest event-specific processes or observation gaps.
Caveat: The five event panels come from the event-spatial AR1 screening model, whereas Overall is the selected operational model persistent field. These are latent link-scale effects, not observed mortality or causal attribution.
Model: persistent_rw1_local_dhw_decomposed_hazards_freshwater_partial_pool (INLA).
File: output/fig/Fig-INLA-07_spatial_residual_fields_v001.png

## Fig-INLA-08_collinearity_reduction — model_comparison

Caption: Held-out performance of the selected full model and a candidate with redundant heat-history and inclusive cumulative-WQC terms removed.
Interpretation: Tests whether reducing genuine ecological redundancy preserves or improves transfer to unseen events and reefs.
Caveat: Structural DHW spline bases and their modifier interactions remain correlated by construction and are assessed jointly rather than pruned.
Model: collinearity_reduced_candidate (INLA).
File: output/fig/Fig-INLA-08_collinearity_reduction_v001.png

## Fig-BRT-03_operational_variable_importance — variable_importance

Caption: Relative influence of predictors in the full-data operational BRT fitted without event identity.
Interpretation: Shows which variables most often improve tree splits; correlated predictors can share influence.
Caveat: Descriptive full-data importance is not held-out skill and does not establish causal importance.
Model: operational_brt_collinearity_screened (BRT).
File: output/fig/Fig-BRT-03_operational_variable_importance_v001.png

## Fig-BRT-04_operational_bootstrap_pdp — partial_dependence

Caption: One-dimensional partial-dependence curves with 95% reef-event cluster-bootstrap intervals for the twelve most influential numeric operational BRT predictors.
Interpretation: Shows the fitted marginal shape and sampling stability of influential BRT relationships.
Caveat: Intervals quantify cluster-resampling uncertainty in the full-data fitted relationship; they are not prediction intervals or proof of causality.
Model: operational_brt_collinearity_screened (BRT).
File: output/fig/Fig-BRT-04_operational_bootstrap_pdp_v001.png

## Fig-INLA-04_event_specific_dhw_response — event_specific_dose_response

Caption: Event-specific DHW response curves from the full-data explanatory INLA model.
Interpretation: Shows how the fitted thermal response differs among observed bleaching events while holding the other fixed effects at their reference values.
Caveat: These curves are explanatory only: event effects are estimated using all data and do not represent prospective forecast skill.
Model: explanatory_event_dhw_interaction_full_data (INLA).
File: output/fig/Fig-INLA-04_event_specific_dhw_response_v001.png

## Fig-BRT-01_event_dhw_variable_importance — variable_importance

Caption: Relative influence from the full-data explanatory BRT with event-by-thermal terms.
Interpretation: Ranks predictor contribution within this fitted BRT; it is descriptive and not used for operational model selection.
Caveat: Relative influence is not a causal effect size and can be shared among correlated predictors.
Model: explanatory_event_dhw_interaction_full_data (BRT).
File: output/fig/Fig-BRT-01_event_dhw_variable_importance_v001.png

## Fig-BRT-02_event_dhw_partial_dependence — partial_dependence

Caption: One-dimensional BRT partial-dependence curves for the most influential numeric predictors.
Interpretation: Describes marginal patterns after averaging the other predictors across observed rows.
Caveat: Partial-dependence curves can be unreliable in sparsely supported predictor combinations; read with the event-specific INLA dose-response curves.
Model: explanatory_event_dhw_interaction_full_data (BRT).
File: output/fig/Fig-BRT-02_event_dhw_partial_dependence_v001.png

## Fig-INLA-09_explanatory_vs_dhw_only — model_comparison

Caption: Observed versus full-data fitted mortality for the explanatory event-DHW ecological model and a programme-layer compound model using only the local-first DHW spline.
Interpretation: Quantifies the additional apparent variation explained by event-specific ecological, disturbance and spatial structure relative to thermal exposure alone.
Caveat: These predictive R-squared values reuse the fitting data and known event/reef effects. They are explanatory apparent fit, not future-event forecast skill.
Model: explanatory_event_dhw_comparison (INLA).
File: output/fig/Fig-INLA-09_explanatory_vs_dhw_only_v001.png

## Fig-INLA-10_cots_probability_severity — model_comparison

Caption: Held-out comparison of the current additive COTS predictors and a complementary outbreak-probability plus probability-weighted-intensity representation.
Interpretation: Tests whether realised COTS intensity above the 0.22 COTS/tow outbreak threshold adds complementary information to the hindcast probability signal.
Caveat: Both terms use the hindcast outbreak definition: probability of density above 0.22 COTS/tow, and probability multiplied by log intensity excess above 0.22 COTS/tow.
Model: cots_probability_severity_candidate (INLA).
File: output/fig/Fig-INLA-10_cots_probability_severity_v001.png

## Fig-INLA-11_cause_aware_validation — model_comparison

Caption: Held-out performance of the current selected model, a cause-filtered thermal/freshwater INLA component, and its combination with independently trained COTS and cyclone hazards.
Interpretation: Tests whether cause-labelled annual losses improve transfer without flattening the thermal response.
Caveat: Cause labels are incomplete and annual cover changes may combine multiple processes; event-overlapping intervals are excluded from each leave-event-out cause fit.
Model: cause_aware_competing_hazards (INLA composite).
File: output/fig/Fig-INLA-11_cause_aware_validation_v001.png

## Fig-DATA-02_cots_cause_evidence — cause_evidence

Caption: Relative coral-cover loss against RRN COTS density for cause-labelled annual transitions from 2016 to 2021, highlighting Gannett Cay, Chinaman, Taylor and Rib reefs.
Interpretation: Shows the outbreak observations used to train the independent COTS hazard and the high-pressure/low-loss cases that constrain it.
Caveat: Cover-transition labels are observational and can contain multiple causes; the 0.22 line is the outbreak threshold, not a deterministic mortality threshold.
Model: cause_aware_cots_evidence (data).
File: output/fig/Fig-DATA-02_cots_cause_evidence_v001.png

## Fig-INLA-12_debbie_residual_audit — spatial_residual_audit

Caption: Cyclone Debbie track, 2017 annual coral losses at reefs with more than 20 hours of waves above 4 m, and held-out residuals before and after the competing cyclone hazard.
Interpretation: Tests whether the spatial underprediction around Debbie is supported by independent wave-exposure and cover-loss evidence, with Penrith as the bleaching-mortality positive control.
Caveat: The annual transition and bleaching-mortality panels use different response records; cyclone labels are observational and the track layer remains provisional.
Model: cause_aware_debbie_audit (INLA composite).
File: output/fig/Fig-INLA-12_debbie_residual_audit_v001.png

## Fig-INLA-13_cots_timing_cover_validation — model_comparison

Caption: Leave-one-event-out performance for COTS timing/starting-cover support and hard, ramped or logistic cyclone activation candidates.
Interpretation: Tests whether operational COTS history and a continuous cyclone activation improve event transfer without using future disturbance labels.
Caveat: The current cyclone exposure layer is provisional; the activation shape must be retested with the forthcoming track/intensity product.
Model: baseline_cots_logistic20_5_cyclone (INLA composite).
File: output/fig/Fig-INLA-13_cots_timing_cover_validation_v001.png

## Fig-DATA-03_cyclone_activation_functions — activation_function

Caption: Hard and continuous cyclone-hazard activation functions compared around the provisional 20-hour damaging-wave reference.
Interpretation: Makes the operational consequence of the cyclone gate explicit before the improved exposure layer is available.
Caveat: These are deterministic sensitivity functions and do not yet propagate cyclone-exposure uncertainty.
Model: baseline_cots_logistic20_5_cyclone (data).
File: output/fig/Fig-DATA-03_cyclone_activation_functions_v001.png

## Fig-INLA-14_lanina_residual_contributors — residual_contributors

Caption: Largest positive and negative reef-event contributors to squared prediction error in the 2020 and 2022 La Nina bleaching events.
Interpretation: Identifies the small number of reefs driving negative event-specific predictive R-squared and separates underprediction from overprediction.
Caveat: Residual attribution is hypothesis-generating; correlated environmental covariates and observation error prevent causal assignment.
Model: baseline_cots_logistic20_5_cyclone (INLA composite).
File: output/fig/Fig-INLA-14_lanina_residual_contributors_v001.png

## Fig-INLA-15_severe_residual_audit — severe_residual_audit

Caption: Held-out observed and predicted mortality for the largest residuals among observations with observed or predicted mortality of at least 50%.
Interpretation: Shows which severe misses remain after COTS/cyclone cause separation and timing/activation refinement.
Caveat: Repeated programme/depth rows can represent the same reef-event and should be interpreted with the reef-event audit table.
Model: baseline_cots_logistic20_5_cyclone (INLA composite).
File: output/fig/Fig-INLA-15_severe_residual_audit_v001.png

## Fig-INLA-16_cots_raw_timing_screen — model_comparison

Caption: Held-out comparison of log-compressed COTS intensity, raw interval density, current-year temporal state and an absolute-cover-loss biomass candidate.
Interpretation: Tests whether retaining the original density contrast and compact temporal support improves transfer without redundant interactions.
Caveat: The selected interval maximum is retrospectively timed and raw density slightly worsens the Gannett severe positive control.
Model: cots_raw_interval_logistic20_5_cyclone (INLA composite).
File: output/fig/Fig-INLA-16_cots_raw_timing_screen_v001.png

## Fig-INLA-17_lanina_occurrence_state — event_occurrence_state

Caption: Observed and held-out predicted event means for the current occurrence layer and continuous RONI/SOI plus rainfall/WQC candidates.
Interpretation: Shows that the tested climate-state corrections improve 2022 but worsen 2020 and suppress 2024.
Caveat: Only five events are available and WQC is a coloured-water/runoff proxy rather than measured discharge.
Model: cots_raw_interval_logistic20_5_cyclone (INLA occurrence calibration).
File: output/fig/Fig-INLA-17_lanina_occurrence_state_v001.png

## Fig-INLA-18_lanina_named_reef_audit — named_reef_residual_audit

Caption: Held-out residuals for named 2020 and 2022 positive and negative reef controls.
Interpretation: Connects the leading La Nina misses to local logger correction, composition, cooling and survey support.
Caveat: The audit is hypothesis-generating and repeated programme rows can describe the same reef-event.
Model: cots_raw_interval_logistic20_5_cyclone (INLA composite).
File: output/fig/Fig-INLA-18_lanina_named_reef_audit_v001.png

## Fig-NOWCAST-01_early_source_coverage — source coverage

Caption: Availability and prevalence of RHIS bleaching observations by April 30 for each event; aerial data are summarised separately because dates are not present in the supplied file.
Interpretation: Rapid in-water observations exist for all five events and can support a within-event occurrence update.
Caveat: April 30 is a common operational cutoff rather than a reef-specific heat stress peak; March 31 is tested as a stricter sensitivity.
Model: spatial_early_bleaching_update (INLA occurrence update).
File: output/fig/Fig-NOWCAST-01_early_source_coverage_v001.png

## Fig-NOWCAST-02_spatial_update_validation — held-out validation

Caption: Change in overall and severe-event RMSE after applying early aerial/RHIS bleaching indicators, with every target sector or latitude block excluded from construction of its update signal.
Interpretation: Negative values indicate that independent early-event observations improve the selected initial operational forecast.
Caveat: Aerial comparisons cover 2016, 2017 and 2020 only; 2022 and 2024 aerial scores have not yet been supplied.
Model: spatial_early_bleaching_update (INLA occurrence update).
File: output/fig/Fig-NOWCAST-02_spatial_update_validation_v001.png

## Fig-NOWCAST-03_event_transfer — event transfer

Caption: Event-specific RMSE change for spatially independent early-event occurrence updates. Each event is predicted from update coefficients fitted to other events.
Interpretation: The figure identifies whether gains are general or dominated by a single event such as 2020.
Caveat: With five events, event-level heterogeneity remains a major source of uncertainty; aerial candidates have only three event folds.
Model: spatial_early_bleaching_update (INLA occurrence update).
File: output/fig/Fig-NOWCAST-03_event_transfer_v001.png

## Fig-COTS-01_prospective_nowcast_validation — validation

Caption: Validation of raw/log prospective COTS pressure states with and without targeted culling context, compared with the retrospective interval reference.
Interpretation: Shows the predictive cost or gain from replacing future interval information with measurements available at event start.
Caveat: Cull removals are targeted evidence and intervention effort, not spatially representative density.
Model: prospective_cots_nowcast (INLA competing hazard).
File: output/fig/Fig-COTS-01_prospective_nowcast_validation_v001.png

## Fig-COTS-02_prospective_focal_controls — focal controls

Caption: Held-out mortality predictions at focal COTS controls using the retrospective interval signal and the best-ranked prospective event-start pressure state.
Interpretation: Makes Gannett and other outbreak-associated losses explicit positive controls rather than relying only on aggregate scores.
Caveat: Observed mortality can include simultaneous thermal or other hazards; culling coverage is targeted.
Model: prospective_cots_nowcast (INLA competing hazard).
File: output/fig/Fig-COTS-02_prospective_focal_controls_v001.png

## Fig-ENS-01_inla_brt_validation — validation

Caption: Nested event-held-out and reef-blocked validation of balanced and severe-weighted residual BRT corrections to the prospective explicit-hazard INLA composite.
Interpretation: Tests whether the machine-learning residual contribution improves severe losses without degrading overall calibration.
Caveat: Blend weights are learned within each training fold; the diagnostic full-data BRT is not used for skill estimates.
Model: inla_brt_residual_ensemble (INLA-BRT).
File: output/fig/Fig-ENS-01_inla_brt_validation_v001.png

## Fig-ENS-02_event_performance — event validation

Caption: Event-specific RMSE for prospective explicit-hazard INLA and its nested residual-BRT candidates.
Interpretation: Shows whether any aggregate or severe-tail gain is purchased by worsening the low-mortality 2020/2022 events.
Caveat: Only one severe observation occurs in 2020 and none in 2022, so event RMSE complements rather than replaces severe-tail assessment.
Model: inla_brt_residual_ensemble (INLA-BRT).
File: output/fig/Fig-ENS-02_event_performance_v001.png

## Fig-ENS-03_heldout_calibration — calibration

Caption: Observed versus predicted mortality for prospective explicit-hazard INLA and residual-BRT ensembles under leave-one-event-out validation.
Interpretation: Provides a direct visual check for underprediction of severe mortality and inflation of low-mortality events.
Caveat: Smooths are descriptive diagnostics and do not represent an additional fitted calibration layer.
Model: inla_brt_residual_ensemble (INLA-BRT).
File: output/fig/Fig-ENS-03_heldout_calibration_v001.png

## Fig-ENS-04_severe_controls — severe controls

Caption: Held-out predictions for Gannett and remaining absolute residuals of at least 0.30 under prospective INLA and residual-BRT ensembles.
Interpretation: Shows whether the severe-tail learner actually addresses focal misses rather than only improving aggregate scores.
Caveat: Focal observations may combine thermal, COTS, cyclone or freshwater processes; the figure is diagnostic rather than causal attribution.
Model: inla_brt_residual_ensemble (INLA-BRT).
File: output/fig/Fig-ENS-04_severe_controls_v001.png

## Fig-ENS-05_residual_brt_importance — variable importance

Caption: Relative influence in the diagnostic residual BRT fitted to cross-fitted prospective-INLA residuals.
Interpretation: Identifies which operational predictors contribute remaining nonlinear structure after explicit hazards.
Caveat: Full-data importance is interpretive; nested outer-fold predictions determine ensemble skill.
Model: inla_brt_residual_ensemble (BRT).
File: output/fig/Fig-ENS-05_residual_brt_importance_v001.png

## Fig-ENS-06_residual_brt_partial_dependence — partial dependence

Caption: Partial-dependence curves for the eight most influential residual-BRT predictors.
Interpretation: Shows the direction and shape of candidate nonlinear corrections remaining after explicit INLA hazards.
Caveat: These full-data marginal curves are descriptive and can be affected by interactions; promotion depends on nested held-out validation.
Model: inla_brt_residual_ensemble (BRT).
File: output/fig/Fig-ENS-06_residual_brt_partial_dependence_v001.png

## Fig-FULLBRT-01_validation — model validation

Caption: Matched validation of the operational INLA composite, direct and two-part standalone BRTs, and pre-specified environmental-envelope ensembles.
Interpretation: Separates new-event transfer from interpolation among reefs in represented events.
Caveat: Envelope weights use predictors only; none of the BRT candidates passes the event-held-out promotion safeguards unless explicitly stated.
Model: operational_rrn_raw_plus_manta_state (INLA-BRT comparison).
File: output/fig/Fig-FULLBRT-01_validation_v001.png

## Fig-FULLBRT-02_event_transfer — event transfer

Caption: Event-specific RMSE for operational INLA, standalone BRTs and smooth envelope ensembles.
Interpretation: Shows whether interpolation gains survive the distinct environmental states of 2016, 2017, 2020, 2022 and 2024.
Caveat: Only five events are available; event identity is never used as an operational predictor.
Model: operational_rrn_raw_plus_manta_state (INLA-BRT comparison).
File: output/fig/Fig-FULLBRT-02_event_transfer_v001.png

## Fig-FULLBRT-03_envelope_performance — applicability domain

Caption: INLA and standalone-BRT RMSE inside and outside each fold-specific multivariate environmental envelope.
Interpretation: Directly tests the hypothesis that BRT is the stronger learner within supported environmental space.
Caveat: Applicability is relative to each training fold and does not imply ecological causality.
Model: operational_rrn_raw_plus_manta_state (INLA-BRT comparison).
File: output/fig/Fig-FULLBRT-03_envelope_performance_v001.png

## Fig-FULLBRT-04_calibration — calibration

Caption: Observed versus predicted mortality under leave-one-event-out validation for INLA, standalone BRTs and the best fixed smooth envelope blend.
Interpretation: Displays mortality-floor and severe-tail calibration that aggregate metrics can obscure.
Caveat: Loess lines are descriptive and are not used to recalibrate predictions.
Model: operational_rrn_raw_plus_manta_state (INLA-BRT comparison).
File: output/fig/Fig-FULLBRT-04_calibration_v001.png

## Fig-FULLBRT-05_local_competence — local competence

Caption: Observation-level BRT absolute-error gain over INLA against training-fold environmental applicability.
Interpretation: Tests whether a denser supported environment is sufficient to justify greater BRT ensemble weight.
Caveat: The smoother is diagnostic only; the evaluated envelope weights are outcome-free and pre-specified.
Model: operational_rrn_raw_plus_manta_state (INLA-BRT comparison).
File: output/fig/Fig-FULLBRT-05_local_competence_v001.png

## Fig-FULLBRT-06_variable_importance — variable importance

Caption: Relative predictor influence for the full-data direct BRT and two-part occurrence and positive-magnitude components.
Interpretation: Shows how the independently fitted BRT frameworks allocate split improvement across operational predictors.
Caveat: Full-data importance is descriptive; correlated mechanisms can share influence and held-out validation governs model promotion.
Model: operational_rrn_raw_plus_manta_state (Standalone BRT).
File: output/fig/Fig-FULLBRT-06_variable_importance_v001.png

## Fig-FULLBRT-07_bootstrap_pdp — partial dependence

Caption: Mortality-scale partial-dependence curves with 95% reef-event cluster-bootstrap intervals for the nine most influential direct-BRT predictors.
Interpretation: Provides visual checks of nonlinear thresholds, slopes and unstable tails in the independently fitted direct BRT.
Caveat: These are marginal full-data diagnostics, not prediction intervals or causal effects.
Model: operational_rrn_raw_plus_manta_state (Standalone direct BRT).
File: output/fig/Fig-FULLBRT-07_bootstrap_pdp_v001.png

## Fig-FULLBRT-08_severe_controls — severe controls

Caption: Held-out predictions at severe and high-residual reef-events for INLA, standalone BRTs and the best fixed smooth envelope blend.
Interpretation: Shows whether an aggregate validation change addresses Gannett, Mackay, Penrith and other ecological controls.
Caveat: Repeated programme/depth observations are summarised to reef-event means.
Model: operational_rrn_raw_plus_manta_state (INLA-BRT comparison).
File: output/fig/Fig-FULLBRT-08_severe_controls_v001.png

## Fig-INLA-01_fixed_effect_posteriors — fixed_effect_posterior

Caption: Posterior means with 80% and 95% credible intervals for selected INLA fixed effects.
Interpretation: Shows direction and uncertainty of population-level link-scale effects.
Caveat: Coefficient magnitudes are not directly comparable across differently transformed predictors.
Model: operational_rrn_raw_plus_manta_state (INLA).
File: output/fig/Fig-INLA-01_fixed_effect_posteriors_v001.png

## Fig-INLA-02_heldout_calibration — heldout_calibration

Caption: Held-out observed versus predicted relative mortality by programme and bleaching event.
Interpretation: Assesses calibration of the selected operational model under event transfer.
Caveat: Predictions are cross-validated; residual structure may reflect unresolved observation and disturbance processes.
Model: operational_rrn_raw_plus_manta_state (INLA).
File: output/fig/Fig-INLA-02_heldout_calibration_v001.png

## Fig-INLA-03_empirical_event_dhw_response — empirical_dhw_response

Caption: Binned held-out observed and predicted mortality across DHW by programme and event.
Interpretation: Visualises event-specific response patterns that aggregate performance metrics can hide.
Caveat: This is an empirical diagnostic, not a covariate-adjusted partial dependence curve.
Model: operational_rrn_raw_plus_manta_state (INLA).
File: output/fig/Fig-INLA-03_empirical_event_dhw_response_v001.png

## Fig-DATA-01_predictor_correlation — predictor_correlation

Caption: Pairwise Pearson correlation matrix for every numeric fixed-effect predictor in the selected model.
Interpretation: Highlights structural redundancy among DHW spline bases, derived interactions and environmental predictors before interpreting individual coefficients.
Caveat: High correlation does not by itself justify removing an ecological mechanism; derived hinge and interaction terms are expected to be correlated and candidate removal must be judged by held-out prediction.
Model: operational_rrn_raw_plus_manta_state (data).
File: output/fig/Fig-DATA-01_predictor_correlation_v001.png

