# Nested cross-fitted BRT residual/severe-tail contribution after explicit
# thermal, cyclone and prospective COTS hazards.

suppressPackageStartupMessages({
  library(dplyr)
  library(gbm)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source('src/lib/model_registry.R')
source('src/lib/model_diagnostics.R')

root <- normalizePath('.', winslash = '/', mustWork = TRUE)
out_dir <- file.path(root, 'output', 'inla_brt_residual_ensemble')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260902L)

metric_summary <- function(rows) rows |> summarise(
  n = n(),
  rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
  mae = mean(abs(observed_mortality - predicted_mortality)),
  predictive_r2 = 1 - sum((observed_mortality-predicted_mortality)^2) /
    sum((observed_mortality-mean(observed_mortality))^2),
  severe_rmse = if (any(observed_mortality >= .5)) sqrt(mean(
    (observed_mortality[observed_mortality >= .5] -
       predicted_mortality[observed_mortality >= .5])^2
  )) else NA_real_,
  false_extreme_rate = mean(
    predicted_mortality[observed_mortality < .3] >= .3
  ),
  mean_bias = mean(observed_mortality - predicted_mortality),
  .groups = 'drop')

nowcast_dir <- file.path(root, 'output', 'prospective_cots_nowcast')
prospective_id <- read_lines(
  file.path(nowcast_dir, 'selected_prospective_candidate.txt')
)[[1]]
nowcast_predictions <- read_csv(
  file.path(nowcast_dir, 'cv_predictions.csv'), show_col_types = FALSE
)

context <- read_csv(
  file.path(root, 'output', 'explanatory_event_dhw',
            'event_dhw_brt_data.csv'), show_col_types = FALSE
) |>
  mutate(source_observation_id = as.character(source_observation_id)) |>
  distinct(source_observation_id, .keep_all = TRUE)

context_features <- c(
  'dhw_novelty10', 'dhw_events_since2016_n6',
  'dhw_years_since_last_n6_capped8', 'secc3m_p10', 'cloudp_90',
  'mcur_90', 'sst_summer_skewness', 'sst_summer_excess_kurtosis',
  'chla_wetseason_median', 'log_coastal_rain30',
  'wqc_prior10_percentile'
)

base <- nowcast_predictions |>
  filter(candidate == prospective_id) |>
  mutate(
    source_observation_id = as.character(source_observation_id),
    fold = as.character(fold),
    base_prediction = predicted_mortality,
    base_residual = observed_mortality - base_prediction,
    log1p_cyclone_wave_hours = log1p(pmax(cyclone_wave_hours, 0)),
    programme_factor = factor(programme_key)
  ) |>
  left_join(
    context |> select(source_observation_id, all_of(context_features)),
    by = 'source_observation_id', relationship = 'many-to-one'
  )

numeric_candidates <- c(
  'base_prediction', 'ann_maxdhw', 'pre_cover', 'acropora',
  'dhw_novelty10', 'dhw_events_since2016_n6',
  'dhw_years_since_last_n6_capped8', 'secc3m_p10', 'cloudp_90',
  'mcur_90', 'sst_summer_skewness', 'sst_summer_excess_kurtosis',
  'chla_wetseason_median', 'log_coastal_rain30',
  'wqc_prior10_percentile', 'log1p_cyclone_wave_hours',
  'cots_outbreak_probability', 'log1p_manta_latest_excess',
  'log1p_manta_peak_excess', 'manta_years_since_peak', 'lon', 'lat'
)
stopifnot(all(numeric_candidates %in% names(base)))

# Outcome-free, fixed collinearity screen. Priority order preserves the base
# prediction and direct mechanism measures before proxy alternatives.
cor_data <- base |>
  filter(scheme == 'leave_one_event_out') |>
  distinct(source_observation_id, .keep_all = TRUE) |>
  select(all_of(numeric_candidates)) |>
  mutate(across(everything(), as.numeric))
cor_matrix <- cor(cor_data, use = 'pairwise.complete.obs', method = 'spearman')
retained_numeric <- character()
screen_rows <- tibble()
for (candidate in numeric_candidates) {
  if (!length(retained_numeric)) {
    retained_numeric <- candidate
    screen_rows <- bind_rows(screen_rows, tibble(
      predictor = candidate, retained = TRUE,
      competing_predictor = NA_character_, absolute_correlation = NA_real_
    ))
    next
  }
  correlations <- abs(cor_matrix[candidate, retained_numeric])
  competing <- retained_numeric[[which.max(correlations)]]
  maximum <- max(correlations, na.rm = TRUE)
  keep <- !is.finite(maximum) || maximum <= .8
  if (keep) retained_numeric <- c(retained_numeric, candidate)
  screen_rows <- bind_rows(screen_rows, tibble(
    predictor = candidate, retained = keep,
    competing_predictor = competing, absolute_correlation = maximum
  ))
}
predictors <- c(retained_numeric, 'programme_factor')
write_csv(screen_rows, file.path(out_dir, 'collinearity_screen.csv'))
write_csv(
  as.data.frame(cor_matrix) |> as_tibble(rownames = 'predictor'),
  file.path(out_dir, 'predictor_correlations.csv')
)

prepare_fold <- function(training, assessment) {
  for (feature in retained_numeric) {
    replacement <- median(training[[feature]][is.finite(training[[feature]])],
                          na.rm = TRUE)
    if (!is.finite(replacement)) replacement <- 0
    training[[feature]][!is.finite(training[[feature]])] <- replacement
    assessment[[feature]][!is.finite(assessment[[feature]])] <- replacement
  }
  levels <- sort(unique(as.character(training$programme_factor)))
  training$programme_factor <- factor(training$programme_factor, levels = levels)
  assessment$programme_factor <- factor(
    as.character(assessment$programme_factor), levels = levels
  )
  list(training = training, assessment = assessment)
}

model_weights <- function(rows, severe_multiplier) {
  reef_n <- table(rows$ReefID)
  programme_n <- table(rows$programme_key)
  weight <- 1 / as.numeric(reef_n[rows$ReefID]) /
    as.numeric(programme_n[rows$programme_key])
  if (severe_multiplier > 1) {
    weight <- weight * ifelse(rows$observed_mortality >= .5,
                              severe_multiplier, 1)
  }
  weight / mean(weight)
}

fit_residual_brt <- function(training, assessment, severe_multiplier, seed) {
  prepared <- prepare_fold(training, assessment)
  training <- prepared$training
  assessment <- prepared$assessment
  set.seed(seed)
  fit <- gbm(
    as.formula(paste(
      'base_residual ~', paste(predictors, collapse = ' + ')
    )),
    data = training, distribution = 'gaussian',
    weights = model_weights(training, severe_multiplier),
    n.trees = 700, interaction.depth = 2, shrinkage = .02,
    n.minobsinnode = 12, bag.fraction = .7,
    train.fraction = 1, keep.data = TRUE, verbose = FALSE
  )
  prediction <- predict(
    fit, assessment, n.trees = 700, type = 'response'
  )
  list(fit = fit, prediction = prediction, assessment = assessment)
}

estimate_nested_alpha <- function(training, scheme, severe_multiplier, seed) {
  inner_groups <- if (scheme == 'leave_one_event_out') {
    sort(unique(training$event_year))
  } else sort(unique(training$fold))
  inner <- tibble()
  for (j in seq_along(inner_groups)) {
    group <- inner_groups[[j]]
    validation <- if (scheme == 'leave_one_event_out') {
      training |> filter(event_year == group)
    } else training |> filter(fold == group)
    analysis <- if (scheme == 'leave_one_event_out') {
      training |> filter(event_year != group)
    } else training |> filter(fold != group)
    fitted <- fit_residual_brt(
      analysis, validation, severe_multiplier, seed + j
    )
    inner <- bind_rows(inner, validation |>
      transmute(
        observed_mortality, base_prediction,
        residual_prediction = fitted$prediction
      ))
  }
  weights <- ifelse(inner$observed_mortality >= .5, severe_multiplier, 1)
  target <- inner$observed_mortality - inner$base_prediction
  denominator <- sum(weights * inner$residual_prediction^2)
  alpha <- if (denominator > 1e-10) sum(
    weights * inner$residual_prediction * target
  ) / denominator else 0
  pmin(pmax(alpha, 0), 1)
}

candidate_specs <- tribble(
  ~candidate, ~severe_multiplier,
  'INLA + balanced residual BRT', 1,
  'INLA + severe-weighted residual BRT', 4
)

ensemble_predictions <- tibble()
blend_weights <- tibble()
for (scheme_name in unique(base$scheme)) {
  scheme_rows <- base |> filter(scheme == scheme_name)
  for (fold_name in unique(scheme_rows$fold)) {
    assessment <- scheme_rows |> filter(fold == fold_name)
    training <- scheme_rows |> filter(fold != fold_name)
    for (i in seq_len(nrow(candidate_specs))) {
      spec <- candidate_specs[i, ]
      seed <- 20260902L + match(scheme_name, unique(base$scheme)) * 1000L +
        match(fold_name, unique(scheme_rows$fold)) * 20L + i
      alpha <- estimate_nested_alpha(
        training, scheme_name, spec$severe_multiplier, seed
      )
      fitted <- fit_residual_brt(
        training, assessment, spec$severe_multiplier, seed + 10L
      )
      correction <- alpha * fitted$prediction
      ensemble_predictions <- bind_rows(
        ensemble_predictions,
        assessment |>
          mutate(
            brt_residual_prediction = fitted$prediction,
            blend_weight = alpha,
            ensemble_correction = correction,
            predicted_mortality = pmin(pmax(
              base_prediction + correction, 0
            ), 1),
            residual = observed_mortality - predicted_mortality,
            candidate = spec$candidate,
            severe_multiplier = spec$severe_multiplier
          )
      )
      blend_weights <- bind_rows(blend_weights, tibble(
        candidate = spec$candidate, scheme = scheme_name,
        fold = fold_name, blend_weight = alpha,
        training_n = nrow(training), assessment_n = nrow(assessment)
      ))
    }
  }
}

baseline <- base |>
  mutate(candidate = 'Prospective INLA explicit hazards',
         predicted_mortality = base_prediction,
         residual = observed_mortality - predicted_mortality)
retrospective <- nowcast_predictions |>
  filter(candidate == 'selected_retrospective_raw_interval') |>
  mutate(
    source_observation_id = as.character(source_observation_id),
    fold = as.character(fold),
    candidate = 'Retrospective COTS INLA reference'
  )
all_predictions <- bind_rows(
  baseline, ensemble_predictions,
  retrospective |> select(any_of(names(baseline)))
)

comparison <- all_predictions |>
  group_by(candidate, scheme) |>
  metric_summary() |>
  arrange(scheme, rmse)
event_metrics <- all_predictions |>
  filter(scheme == 'leave_one_event_out') |>
  group_by(candidate, event_year) |>
  metric_summary()

base_metrics <- comparison |>
  filter(candidate == 'Prospective INLA explicit hazards')
base_events <- event_metrics |>
  filter(candidate == 'Prospective INLA explicit hazards')
eligibility <- comparison |>
  filter(str_detect(candidate, 'BRT')) |>
  select(candidate, scheme, rmse, severe_rmse, false_extreme_rate) |>
  pivot_wider(names_from = scheme,
              values_from = c(rmse, severe_rmse, false_extreme_rate)) |>
  left_join(
    event_metrics |>
      filter(event_year %in% c(2020L, 2022L), str_detect(candidate, 'BRT')) |>
      select(candidate, event_year, rmse) |>
      pivot_wider(names_from = event_year, values_from = rmse,
                  names_prefix = 'rmse_'), by = 'candidate'
  ) |>
  mutate(
    passes = rmse_leave_one_event_out <=
        base_metrics$rmse[base_metrics$scheme == 'leave_one_event_out'] &
      severe_rmse_leave_one_event_out <
        base_metrics$severe_rmse[base_metrics$scheme == 'leave_one_event_out'] &
      severe_rmse_reef_blocked_5fold <
        base_metrics$severe_rmse[base_metrics$scheme == 'reef_blocked_5fold'] &
      false_extreme_rate_leave_one_event_out <=
        base_metrics$false_extreme_rate[
          base_metrics$scheme == 'leave_one_event_out'
        ] &
      rmse_2020 <= base_events$rmse[base_events$event_year == 2020] + .001 &
      rmse_2022 <= base_events$rmse[base_events$event_year == 2022] + .001
  ) |>
  arrange(desc(passes), severe_rmse_leave_one_event_out,
          rmse_leave_one_event_out)
selected_ensemble <- if (any(eligibility$passes)) {
  eligibility$candidate[which(eligibility$passes)[[1]]]
} else 'Prospective INLA explicit hazards'

write_csv(all_predictions, file.path(out_dir, 'cv_predictions.csv'))
write_csv(comparison, file.path(out_dir, 'model_comparison.csv'))
write_csv(event_metrics, file.path(out_dir, 'event_metrics.csv'))
write_csv(blend_weights, file.path(out_dir, 'blend_weights.csv'))
write_csv(eligibility, file.path(out_dir, 'promotion_assessment.csv'))
write_lines(selected_ensemble, file.path(out_dir, 'selected_ensemble.txt'))

severe_audit <- all_predictions |>
  filter(
    scheme == 'leave_one_event_out',
    candidate %in% c('Prospective INLA explicit hazards',
                     candidate_specs$candidate),
    observed_mortality >= .5 |
      str_detect(str_to_lower(ReefName), 'gannett')
  ) |>
  group_by(candidate, ReefID, ReefName, event_year) |>
  summarise(
    observed_mortality = mean(observed_mortality),
    predicted_mortality = mean(predicted_mortality),
    residual = mean(residual),
    DHW = mean(ann_maxdhw),
    cots_probability = mean(cots_outbreak_probability),
    .groups = 'drop'
  ) |>
  arrange(desc(abs(residual)))
write_csv(severe_audit, file.path(out_dir, 'severe_and_gannett_audit.csv'))

# Full-data residual BRT is interpretation only; its response remains the
# cross-fitted prospective-INLA residual, not the in-sample residual.
diagnostic_candidate <- if (selected_ensemble ==
                            'Prospective INLA explicit hazards') {
  'INLA + severe-weighted residual BRT'
} else selected_ensemble
diagnostic_multiplier <- candidate_specs$severe_multiplier[
  match(diagnostic_candidate, candidate_specs$candidate)
]
diagnostic_rows <- base |>
  filter(scheme == 'leave_one_event_out') |>
  distinct(source_observation_id, .keep_all = TRUE)
diagnostic_fit <- fit_residual_brt(
  diagnostic_rows, diagnostic_rows, diagnostic_multiplier, 20261999L
)$fit
saveRDS(diagnostic_fit, file.path(out_dir, 'residual_brt_diagnostic.rds'))
influence <- summary(diagnostic_fit, n.trees = 700, plotit = FALSE) |>
  as_tibble() |>
  transmute(variable = var, relative_influence = rel.inf) |>
  arrange(desc(relative_influence))
write_csv(influence, file.path(out_dir, 'residual_brt_importance.csv'))

top_variables <- influence |>
  filter(variable != 'programme_factor') |>
  slice_head(n = 8) |>
  pull(variable)
pdp <- map_dfr(top_variables, function(variable) {
  variable_index <- match(variable, diagnostic_fit$var.names)
  partial <- gbm::plot.gbm(
    diagnostic_fit, i.var = variable_index, n.trees = 700,
    return.grid = TRUE
  )
  tibble(variable = variable, value = partial[[1]],
         residual_correction = partial[[2]])
})
write_csv(pdp, file.path(out_dir, 'residual_brt_partial_dependence.csv'))

# Registered figures.
validation_data <- comparison |>
  filter(candidate != 'Retrospective COTS INLA reference') |>
  select(candidate, scheme, rmse, severe_rmse, false_extreme_rate) |>
  pivot_longer(c(rmse, severe_rmse, false_extreme_rate),
               names_to = 'metric', values_to = 'value') |>
  mutate(
    metric = recode(metric, rmse = 'RMSE', severe_rmse = 'Severe RMSE',
                    false_extreme_rate = 'False-extreme rate'),
    scheme = recode(scheme, leave_one_event_out = 'Event held out',
                    reef_blocked_5fold = 'Reef blocked')
  )
p_validation <- ggplot(
  validation_data,
  aes(value, reorder(candidate, value), colour = scheme)
) + geom_point(size = 3, position = position_dodge(width = .5)) +
  facet_wrap(~metric, scales = 'free_x') +
  labs(
    title = 'Prospective INLA–BRT ensemble validation',
    subtitle = 'Nested blend weights; lower RMSE and false-extreme rate are better',
    x = NULL, y = NULL, colour = 'Validation'
  ) + theme_bw(base_size = 11) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_validation, 'Fig-ENS-01_inla_brt_validation', validation_data,
  'Nested event-held-out and reef-blocked validation of balanced and severe-weighted residual BRT corrections to the prospective explicit-hazard INLA composite.',
  'Tests whether the machine-learning residual contribution improves severe losses without degrading overall calibration.',
  'Blend weights are learned within each training fold; the diagnostic full-data BRT is not used for skill estimates.',
  'inla_brt_residual_ensemble', 'INLA-BRT', 'validation',
  'operational_candidate', root, TRUE,
  code_source = 'src/evaluation/test_inla_brt_residual_ensemble.R',
  width = 11.5, height = 6.5
)

event_plot_data <- event_metrics |>
  filter(candidate != 'Retrospective COTS INLA reference')
p_events <- ggplot(
  event_plot_data,
  aes(factor(event_year), rmse, colour = candidate, group = candidate)
) + geom_line(linewidth = .7) + geom_point(size = 2.4) +
  labs(
    title = 'Event transfer of the prospective INLA–BRT ensemble',
    subtitle = '2020 and 2022 are explicit low-mortality safeguards',
    x = 'Held-out event', y = 'RMSE', colour = 'Model'
  ) + theme_bw(base_size = 11) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_events, 'Fig-ENS-02_event_performance', event_plot_data,
  'Event-specific RMSE for prospective explicit-hazard INLA and its nested residual-BRT candidates.',
  'Shows whether any aggregate or severe-tail gain is purchased by worsening the low-mortality 2020/2022 events.',
  'Only one severe observation occurs in 2020 and none in 2022, so event RMSE complements rather than replaces severe-tail assessment.',
  'inla_brt_residual_ensemble', 'INLA-BRT', 'event validation',
  'operational_candidate', root, TRUE,
  code_source = 'src/evaluation/test_inla_brt_residual_ensemble.R',
  width = 10.5, height = 6.5
)

calibration_data <- all_predictions |>
  filter(
    scheme == 'leave_one_event_out',
    candidate %in% c('Prospective INLA explicit hazards',
                     candidate_specs$candidate)
  )
p_calibration <- ggplot(
  calibration_data,
  aes(predicted_mortality, observed_mortality, colour = factor(event_year))
) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = .55, size = 1.5) +
  geom_smooth(method = 'loess', se = FALSE, linewidth = .8) +
  facet_wrap(~candidate) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = 'Held-out calibration of prospective INLA and INLA–BRT',
    subtitle = 'Event colours expose whether tail gains distort low-mortality years',
    x = 'Predicted relative mortality', y = 'Observed relative mortality',
    colour = 'Event'
  ) + theme_bw(base_size = 10) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_calibration, 'Fig-ENS-03_heldout_calibration', calibration_data,
  'Observed versus predicted mortality for prospective explicit-hazard INLA and residual-BRT ensembles under leave-one-event-out validation.',
  'Provides a direct visual check for underprediction of severe mortality and inflation of low-mortality events.',
  'Smooths are descriptive diagnostics and do not represent an additional fitted calibration layer.',
  'inla_brt_residual_ensemble', 'INLA-BRT', 'calibration',
  'operational_candidate', root, TRUE,
  code_source = 'src/evaluation/test_inla_brt_residual_ensemble.R',
  width = 13, height = 7
)

focal_data <- severe_audit |>
  filter(str_detect(str_to_lower(ReefName), 'gannett') |
           abs(residual) >= .3) |>
  mutate(label = paste0(ReefName, ' (', event_year, ')'))
p_focal <- ggplot(
  focal_data,
  aes(predicted_mortality, label, colour = candidate)
) + geom_point(size = 2.7, position = position_dodge(width = .55)) +
  geom_point(aes(x = observed_mortality), shape = 4, colour = 'black',
             size = 3, stroke = 1) +
  labs(
    title = 'Severe-tail controls after explicit hazards',
    subtitle = 'Black crosses are observations; Gannett is the COTS positive control',
    x = 'Held-out predicted mortality', y = NULL, colour = 'Model'
  ) + theme_bw(base_size = 10) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_focal, 'Fig-ENS-04_severe_controls', focal_data,
  'Held-out predictions for Gannett and remaining absolute residuals of at least 0.30 under prospective INLA and residual-BRT ensembles.',
  'Shows whether the severe-tail learner actually addresses focal misses rather than only improving aggregate scores.',
  'Focal observations may combine thermal, COTS, cyclone or freshwater processes; the figure is diagnostic rather than causal attribution.',
  'inla_brt_residual_ensemble', 'INLA-BRT', 'severe controls',
  'operational_candidate', root, TRUE,
  code_source = 'src/evaluation/test_inla_brt_residual_ensemble.R',
  width = 12, height = 7.5
)

importance_plot_data <- influence |>
  mutate(variable = factor(variable, levels = rev(variable)))
p_importance <- ggplot(
  importance_plot_data,
  aes(relative_influence, variable)
) + geom_col(fill = '#7B3294') +
  labs(
    title = 'Residual BRT variable importance',
    subtitle = 'Cross-fitted prospective-INLA residual is the response',
    x = 'Relative influence (%)', y = NULL
  ) + theme_bw(base_size = 11)
save_figure_bundle(
  p_importance, 'Fig-ENS-05_residual_brt_importance', importance_plot_data,
  'Relative influence in the diagnostic residual BRT fitted to cross-fitted prospective-INLA residuals.',
  'Identifies which operational predictors contribute remaining nonlinear structure after explicit hazards.',
  'Full-data importance is interpretive; nested outer-fold predictions determine ensemble skill.',
  'inla_brt_residual_ensemble', 'BRT', 'variable importance',
  'operational_diagnostic', root, TRUE,
  code_source = 'src/evaluation/test_inla_brt_residual_ensemble.R',
  width = 10, height = 7.5
)

p_pdp <- ggplot(pdp, aes(value, residual_correction)) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey45') +
  geom_line(colour = '#7B3294', linewidth = .9) +
  facet_wrap(~variable, scales = 'free_x', ncol = 3) +
  labs(
    title = 'Residual BRT partial dependence',
    subtitle = 'Positive values add mortality to the prospective INLA prediction',
    x = 'Predictor value', y = 'Partial residual correction'
  ) + theme_bw(base_size = 10)
save_figure_bundle(
  p_pdp, 'Fig-ENS-06_residual_brt_partial_dependence', pdp,
  'Partial-dependence curves for the eight most influential residual-BRT predictors.',
  'Shows the direction and shape of candidate nonlinear corrections remaining after explicit INLA hazards.',
  'These full-data marginal curves are descriptive and can be affected by interactions; promotion depends on nested held-out validation.',
  'inla_brt_residual_ensemble', 'BRT', 'partial dependence',
  'operational_diagnostic', root, TRUE,
  code_source = 'src/evaluation/test_inla_brt_residual_ensemble.R',
  width = 11.5, height = 8
)
write_figure_readme(root)
message('INLA-BRT residual ensemble outputs written to: ', out_dir)
