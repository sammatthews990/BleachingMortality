# Standalone two-part BRT and environmental-envelope INLA ensemble.
# The BRT predicts occurrence and positive magnitude independently of INLA.
# Every outer assessment is an unseen event or reef block. Environmental
# applicability uses predictors only and is recalculated from each training set.

suppressPackageStartupMessages({
  library(dplyr)
  library(gbm)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source('scripts/model_registry.R')
source('scripts/model_diagnostics.R')

root <- project_root()
out_dir <- file.path(root, 'output', 'standalone_brt_envelope_ensemble')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260903L)

selected_id <- read_lines(file.path(
  root, 'output', 'prospective_cots_nowcast',
  'selected_prospective_candidate.txt'
))[[1]]
inla_predictions <- read_csv(file.path(
  root, 'output', 'prospective_cots_nowcast', 'cv_predictions.csv'
), show_col_types = FALSE) |>
  filter(candidate == selected_id) |>
  mutate(
    source_observation_id = as.character(source_observation_id),
    fold = as.character(fold)
  )

context <- read_csv(file.path(
  root, 'output', 'explanatory_event_dhw', 'event_dhw_brt_data.csv'
), show_col_types = FALSE) |>
  mutate(source_observation_id = as.character(source_observation_id)) |>
  distinct(source_observation_id, .keep_all = TRUE)

cots_state <- read_csv(file.path(
  root, 'data', 'processed', 'cots_event_start_features.csv'
), show_col_types = FALSE) |>
  mutate(ReefID = str_to_upper(str_trim(ReefID))) |>
  distinct(ReefID, event_year, .keep_all = TRUE)

rows <- context |>
  mutate(
    ReefID = str_to_upper(str_trim(ReefID)),
    observed_mortality = mortality_prop,
    observed_occurrence = as.numeric(mortality_prop > 0),
    programme_factor = droplevels(factor(programme_key)),
    region_factor = droplevels(factor(region_block)),
    rrn_event_excess_raw = pmax(coalesce(cot_interval_idw_max, 0) - .22, 0),
    log1p_cyclone_wave_hours = log1p(pmax(cyc_interval_maxHrs4mw, 0)),
    disease_applicable_numeric = as.numeric(coalesce(disease_risk_applicable, FALSE)),
    wqc_current_missing = as.numeric(!is.finite(wqc_freqcc12))
  ) |>
  left_join(cots_state, by = c('ReefID', 'event_year'),
            relationship = 'many-to-one') |>
  mutate(
    across(c(manta_latest_excess_raw, manta_peak_excess_raw,
             manta_years_since_peak, manta_supported,
             cull_removed_per_dive, log1p_cull_dives), ~ coalesce(.x, 0))
  )

stopifnot(
  nrow(rows) == n_distinct(rows$source_observation_id),
  all(inla_predictions$source_observation_id %in% rows$source_observation_id)
)

# Trees learn nonlinear thresholds directly; explicit hinges and hand-built
# interactions are not duplicated. Priority resolves non-structural |rho|>.8.
numeric_candidates <- c(
  'ann_maxdhw', 'prop_acropora_pre', 'observed_pre_cover',
  'dhw_novelty10', 'dhw_events_since2016_n6',
  'secc3m_p10', 'cloudp_90', 'mcur_90', 'sst_summer_skewness',
  'log_coastal_rain30', 'wqc_freqcc12', 'wqc_prior10_percentile',
  'era5_wind_calm_fraction', 'log1p_cyclone_wave_hours',
  'tc_interval_wind_distance_index', 'cots_outbreak_probability',
  'rrn_event_excess_raw', 'manta_latest_excess_raw',
  'manta_peak_excess_raw', 'manta_years_since_peak',
  'disease_risk_max', 'depth', 'manta_supported',
  'disease_applicable_numeric', 'wqc_current_missing'
)
factor_candidates <- c('programme_factor', 'region_factor')
stopifnot(all(c(numeric_candidates, factor_candidates) %in% names(rows)))

# Global predictor correlations are outcome-free and all operational mapping
# covariates are observable. A fixed screen keeps the contract identical in
# every fold while medians, factors, tuning and envelope PCs remain fold-local.
cor_matrix <- cor(
  rows |> select(all_of(numeric_candidates)) |>
    mutate(across(everything(), as.numeric)),
  use = 'pairwise.complete.obs', method = 'spearman'
)
retained_numeric <- character()
screen <- tibble()
for (candidate in numeric_candidates) {
  if (!length(retained_numeric)) {
    retained_numeric <- candidate
    screen <- bind_rows(screen, tibble(
      predictor = candidate, retained = TRUE,
      competing_predictor = NA_character_, absolute_correlation = NA_real_
    ))
    next
  }
  correlations <- abs(cor_matrix[candidate, retained_numeric])
  finite_correlations <- which(is.finite(correlations))
  if (!length(finite_correlations)) {
    competitor <- NA_character_
    maximum <- NA_real_
    keep <- TRUE
  } else {
    local_maximum <- finite_correlations[[which.max(
      correlations[finite_correlations]
    )]]
    competitor <- retained_numeric[[local_maximum]]
    maximum <- correlations[[local_maximum]]
    keep <- maximum <= .8
  }
  if (keep) retained_numeric <- c(retained_numeric, candidate)
  screen <- bind_rows(screen, tibble(
    predictor = candidate, retained = keep,
    competing_predictor = competitor, absolute_correlation = maximum
  ))
}
predictors <- c(retained_numeric, factor_candidates)
write_csv(screen, file.path(out_dir, 'collinearity_screen.csv'))
write_csv(as.data.frame(cor_matrix) |> as_tibble(rownames = 'predictor'),
          file.path(out_dir, 'predictor_correlations.csv'))

prepare_data <- function(training, assessment) {
  training <- as.data.frame(training)
  assessment <- as.data.frame(assessment)
  for (feature in retained_numeric) {
    replacement <- median(training[[feature]][is.finite(training[[feature]])],
                          na.rm = TRUE)
    if (!is.finite(replacement)) replacement <- 0
    training[[feature]][!is.finite(training[[feature]])] <- replacement
    assessment[[feature]][!is.finite(assessment[[feature]])] <- replacement
  }
  for (feature in factor_candidates) {
    training[[feature]] <- droplevels(factor(training[[feature]]))
    assessment[[feature]] <- factor(
      as.character(assessment[[feature]]), levels = levels(training[[feature]])
    )
    if (anyNA(assessment[[feature]])) {
      stop('Assessment contains an unseen factor level for ', feature)
    }
  }
  list(training = training, assessment = assessment)
}

reef_event_weights <- function(data) {
  reef_event_n <- table(data$reef_event_key)
  programme_n <- table(data$programme_key)
  weight <- 1 / as.numeric(reef_event_n[data$reef_event_key]) /
    as.numeric(programme_n[data$programme_key])
  weight / mean(weight)
}

fit_component <- function(data, response, distribution, weights, settings,
                          seed) {
  outcome <- data[[response]]
  if (length(unique(outcome)) < 2L || !is.finite(sd(outcome)) || sd(outcome) == 0) {
    return(list(type = 'constant', value = weighted.mean(outcome, weights)))
  }
  set.seed(seed)
  model <- gbm(
    reformulate(predictors, response = response),
    data = data, distribution = distribution, weights = weights,
    n.trees = settings$n_trees,
    interaction.depth = settings$depth,
    shrinkage = settings$shrinkage,
    n.minobsinnode = min(settings$minobs,
                         max(2L, floor(nrow(data) / 8L))),
    bag.fraction = .70, train.fraction = 1,
    keep.data = FALSE, verbose = FALSE
  )
  list(type = 'gbm', model = model)
}

predict_component <- function(component, newdata, n_trees) {
  if (component$type == 'constant') {
    return(rep(component$value, nrow(newdata)))
  }
  predict(component$model, newdata, n.trees = n_trees, type = 'response')
}

fit_two_part <- function(training, assessment, settings, seed) {
  prepared <- prepare_data(training, assessment)
  training <- prepared$training
  assessment <- prepared$assessment
  training$has_loss <- as.numeric(training$observed_mortality > 0)
  weights <- reef_event_weights(training)
  occurrence <- fit_component(
    training, 'has_loss', 'bernoulli', weights, settings, seed
  )
  positive <- training[training$has_loss == 1, , drop = FALSE]
  positive_weights <- weights[training$has_loss == 1]
  epsilon <- min(.01, .5 / nrow(positive))
  positive$positive_logit <- qlogis(pmin(pmax(
    positive$observed_mortality, epsilon
  ), 1 - epsilon))
  magnitude <- fit_component(
    positive, 'positive_logit', 'gaussian', positive_weights,
    settings, seed + 1L
  )
  direct <- fit_component(
    training, 'observed_mortality', 'gaussian', weights,
    settings, seed + 2L
  )
  predicted_occurrence <- pmin(pmax(predict_component(
    occurrence, assessment, settings$n_trees
  ), 0), 1)
  predicted_positive <- plogis(predict_component(
    magnitude, assessment, settings$n_trees
  ))
  predicted_direct <- pmin(pmax(predict_component(
    direct, assessment, settings$n_trees
  ), 0), 1)
  list(
    prediction = tibble(
      brt_occurrence = predicted_occurrence,
      brt_positive_mortality = predicted_positive,
      brt_prediction = predicted_occurrence * predicted_positive,
      direct_brt_prediction = predicted_direct
    ),
    model = list(occurrence = occurrence, magnitude = magnitude,
                 direct = direct),
    training = training, assessment = assessment
  )
}

metric_summary <- function(data, prediction = 'predicted_mortality',
                           occurrence = 'predicted_occurrence') {
  predicted <- data[[prediction]]
  predicted_occurrence <- data[[occurrence]]
  observed <- data$observed_mortality
  severe <- observed >= .5
  tibble(
    n = length(observed),
    rmse = sqrt(mean((observed - predicted)^2)),
    mae = mean(abs(observed - predicted)),
    predictive_r2 = 1 - sum((observed - predicted)^2) /
      sum((observed - mean(observed))^2),
    severe_rmse = if (any(severe)) sqrt(mean(
      (observed[severe] - predicted[severe])^2
    )) else NA_real_,
    occurrence_brier = mean((data$observed_occurrence -
                               predicted_occurrence)^2),
    false_extreme_rate = mean(predicted[observed < .3] >= .3),
    mean_bias = mean(observed - predicted)
  )
}

tuning_grid <- tribble(
  ~setting_id, ~n_trees, ~depth, ~shrinkage, ~minobs,
  'additive', 800L, 1L, .02, 8L,
  'pairwise', 800L, 2L, .02, 8L,
  'three_way', 800L, 3L, .02, 8L
)

tune_two_part <- function(training, inner_group, outer_key) {
  scores <- tibble()
  inner_levels <- sort(unique(inner_group))
  for (setting_index in seq_len(nrow(tuning_grid))) {
    setting <- tuning_grid[setting_index, ]
    inner_predictions <- tibble()
    for (inner_index in seq_along(inner_levels)) {
      held <- inner_levels[[inner_index]]
      analysis <- training[inner_group != held, , drop = FALSE]
      assessment <- training[inner_group == held, , drop = FALSE]
      fitted <- fit_two_part(
        analysis, assessment, setting,
        seed = 20260903L + outer_key * 100L +
          setting_index * 10L + inner_index
      )
      inner_predictions <- bind_rows(
        inner_predictions,
        assessment |> select(source_observation_id, observed_mortality,
                             observed_occurrence) |>
          bind_cols(fitted$prediction)
      )
    }
    metrics <- metric_summary(
      inner_predictions, 'brt_prediction', 'brt_occurrence'
    )
    direct_metrics <- metric_summary(
      inner_predictions, 'direct_brt_prediction', 'brt_occurrence'
    )
    severe_penalty <- ifelse(is.finite(metrics$severe_rmse),
                             metrics$severe_rmse, metrics$rmse)
    scores <- bind_rows(scores, setting |>
      mutate(
        inner_rmse = metrics$rmse,
        inner_severe_rmse = metrics$severe_rmse,
        inner_occurrence_brier = metrics$occurrence_brier,
        tuning_score = inner_rmse + .20 * severe_penalty +
          .05 * inner_occurrence_brier,
        direct_inner_rmse = direct_metrics$rmse,
        direct_inner_severe_rmse = direct_metrics$severe_rmse,
        direct_tuning_score = direct_inner_rmse + .20 * ifelse(
          is.finite(direct_inner_severe_rmse),
          direct_inner_severe_rmse, direct_inner_rmse
        ) + .05 * inner_occurrence_brier
      ))
  }
  list(
    best = scores |> arrange(tuning_score, inner_rmse) |> slice(1),
    best_direct = scores |>
      arrange(direct_tuning_score, direct_inner_rmse) |> slice(1),
    scores = scores
  )
}

environmental_envelope <- function(training, assessment, k = 10L) {
  training_matrix <- as.matrix(training[, retained_numeric, drop = FALSE])
  assessment_matrix <- as.matrix(assessment[, retained_numeric, drop = FALSE])
  centre <- colMeans(training_matrix)
  spread <- apply(training_matrix, 2, sd)
  spread[!is.finite(spread) | spread < 1e-8] <- 1
  training_scaled <- sweep(sweep(training_matrix, 2, centre), 2, spread, '/')
  assessment_scaled <- sweep(sweep(assessment_matrix, 2, centre), 2, spread, '/')
  variable <- apply(training_scaled, 2, sd) > 1e-8
  training_scaled <- training_scaled[, variable, drop = FALSE]
  assessment_scaled <- assessment_scaled[, variable, drop = FALSE]
  pca <- prcomp(training_scaled, center = FALSE, scale. = FALSE)
  cumulative <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  components <- min(which(cumulative >= .90), 8L, ncol(training_scaled))
  train_scores <- predict(pca, training_scaled)[, seq_len(components), drop = FALSE]
  assessment_scores <- predict(pca, assessment_scaled)[, seq_len(components), drop = FALSE]
  pc_scale <- pmax(pca$sdev[seq_len(components)], 1e-8)
  train_scores <- sweep(train_scores, 2, pc_scale, '/')
  assessment_scores <- sweep(assessment_scores, 2, pc_scale, '/')

  all_scores <- rbind(train_scores, assessment_scores)
  distance_matrix <- as.matrix(dist(all_scores))
  n_training <- nrow(train_scores)
  train_distance <- distance_matrix[seq_len(n_training), seq_len(n_training), drop = FALSE]
  diag(train_distance) <- Inf
  effective_k <- min(k, n_training - 1L)
  training_knn <- apply(train_distance, 1, function(x) {
    mean(sort(x, partial = seq_len(effective_k))[seq_len(effective_k)])
  })
  assessment_distance <- distance_matrix[
    n_training + seq_len(nrow(assessment_scores)), seq_len(n_training),
    drop = FALSE
  ]
  assessment_knn <- apply(assessment_distance, 1, function(x) {
    mean(sort(x, partial = seq_len(effective_k))[seq_len(effective_k)])
  })
  threshold <- quantile(training_knn, .95, na.rm = TRUE)
  transition <- max(IQR(training_knn, na.rm = TRUE), .15 * threshold, 1e-6)

  training_min <- apply(training_matrix, 2, min)
  training_max <- apply(training_matrix, 2, max)
  outside_fraction <- rowMeans(
    sweep(assessment_matrix, 2, training_min, '<') |
      sweep(assessment_matrix, 2, training_max, '>')
  )
  applicability <- plogis((threshold - assessment_knn) / transition) *
    exp(-4 * outside_fraction)
  tibble(
    envelope_distance = assessment_knn,
    envelope_threshold = as.numeric(threshold),
    outside_range_fraction = outside_fraction,
    applicability = pmin(pmax(applicability, 0), 1),
    inside_envelope = assessment_knn <= threshold & outside_fraction <= .05,
    envelope_components = components
  )
}

all_predictions <- tibble()
tuning_results <- tibble()
for (scheme_name in unique(inla_predictions$scheme)) {
  folds <- unique(inla_predictions$fold[inla_predictions$scheme == scheme_name])
  for (fold_index in seq_along(folds)) {
    fold_name <- folds[[fold_index]]
    base <- inla_predictions |>
      filter(scheme == scheme_name, fold == fold_name)
    assessment <- rows[match(base$source_observation_id,
                             rows$source_observation_id), , drop = FALSE]
    if (anyNA(assessment$source_observation_id)) {
      stop('Could not match all assessment observations')
    }
    if (scheme_name == 'leave_one_event_out') {
      held_event <- as.integer(fold_name)
      training <- rows |> filter(event_year != held_event)
      inner_group <- training$event_year
    } else {
      held_reefs <- unique(assessment$ReefID)
      training <- rows |> filter(!ReefID %in% held_reefs)
      inner_group <- training$joint_reef_fold
    }
    outer_key <- match(scheme_name, unique(inla_predictions$scheme)) * 10L +
      fold_index
    tuned <- tune_two_part(training, inner_group, outer_key)
    tuning_results <- bind_rows(
      tuning_results,
      tuned$scores |> mutate(scheme = scheme_name, fold = fold_name,
        selected_two_part = setting_id == tuned$best$setting_id,
        selected_direct = setting_id == tuned$best_direct$setting_id)
    )
    fitted <- fit_two_part(
      training, assessment, tuned$best,
      seed = 20261903L + outer_key
    )
    direct_fitted <- if (tuned$best_direct$setting_id ==
                         tuned$best$setting_id) fitted else fit_two_part(
      training, assessment, tuned$best_direct,
      seed = 20262903L + outer_key
    )
    envelope <- environmental_envelope(
      fitted$training, fitted$assessment
    )
    scored <- assessment |>
      transmute(
        source_observation_id, programme_key, ReefID, ReefName,
        event_year, observed_mortality, observed_occurrence,
        ann_maxdhw, scheme = scheme_name, fold = fold_name
      ) |>
      bind_cols(fitted$prediction, envelope) |>
      mutate(
        inla_prediction = base$predicted_mortality,
        inla_occurrence = base$predicted_occurrence,
        hard_weight = as.numeric(inside_envelope),
        smooth_weight = .5 * applicability,
        champion_weight = applicability,
        direct_brt_prediction =
          direct_fitted$prediction$direct_brt_prediction
      )
    candidate_rows <- bind_rows(
      scored |> transmute(
        across(everything()), candidate = 'INLA explicit hazards',
        predicted_mortality = inla_prediction,
        predicted_occurrence = inla_occurrence,
        brt_weight = 0
      ),
      scored |> transmute(
        across(everything()), candidate = 'Standalone two-part BRT',
        predicted_mortality = brt_prediction,
        predicted_occurrence = brt_occurrence,
        brt_weight = 1
      ),
      scored |> transmute(
        across(everything()), candidate = 'Standalone direct BRT',
        predicted_mortality = direct_brt_prediction,
        predicted_occurrence = brt_occurrence,
        brt_weight = 1
      ),
      scored |> transmute(
        across(everything()), candidate = 'Two-part hard envelope champion',
        predicted_mortality = hard_weight * brt_prediction +
          (1 - hard_weight) * inla_prediction,
        predicted_occurrence = hard_weight * brt_occurrence +
          (1 - hard_weight) * inla_occurrence,
        brt_weight = hard_weight
      ),
      scored |> transmute(
        across(everything()), candidate = 'Two-part smooth 50% ensemble',
        predicted_mortality = smooth_weight * brt_prediction +
          (1 - smooth_weight) * inla_prediction,
        predicted_occurrence = smooth_weight * brt_occurrence +
          (1 - smooth_weight) * inla_occurrence,
        brt_weight = smooth_weight
      ),
      scored |> transmute(
        across(everything()), candidate = 'Two-part applicability champion',
        predicted_mortality = champion_weight * brt_prediction +
          (1 - champion_weight) * inla_prediction,
        predicted_occurrence = champion_weight * brt_occurrence +
          (1 - champion_weight) * inla_occurrence,
        brt_weight = champion_weight
      ),
      scored |> transmute(
        across(everything()), candidate = 'Direct hard envelope champion',
        predicted_mortality = hard_weight * direct_brt_prediction +
          (1 - hard_weight) * inla_prediction,
        predicted_occurrence = hard_weight * brt_occurrence +
          (1 - hard_weight) * inla_occurrence,
        brt_weight = hard_weight
      ),
      scored |> transmute(
        across(everything()), candidate = 'Direct smooth 50% ensemble',
        predicted_mortality = smooth_weight * direct_brt_prediction +
          (1 - smooth_weight) * inla_prediction,
        predicted_occurrence = smooth_weight * brt_occurrence +
          (1 - smooth_weight) * inla_occurrence,
        brt_weight = smooth_weight
      ),
      scored |> transmute(
        across(everything()), candidate = 'Direct applicability champion',
        predicted_mortality = champion_weight * direct_brt_prediction +
          (1 - champion_weight) * inla_prediction,
        predicted_occurrence = champion_weight * brt_occurrence +
          (1 - champion_weight) * inla_occurrence,
        brt_weight = champion_weight
      )
    )
    all_predictions <- bind_rows(all_predictions, candidate_rows)
    message('Completed ', scheme_name, ' fold ', fold_name,
            ' with two-part ', tuned$best$setting_id,
            ' and direct ', tuned$best_direct$setting_id)
  }
}

model_comparison <- all_predictions |>
  group_by(candidate, scheme) |>
  group_modify(~ metric_summary(.x)) |>
  ungroup()
event_metrics <- all_predictions |>
  filter(scheme == 'leave_one_event_out') |>
  group_by(candidate, event_year) |>
  group_modify(~ metric_summary(.x)) |>
  ungroup()
envelope_metrics <- all_predictions |>
  filter(candidate %in% c(
    'INLA explicit hazards', 'Standalone two-part BRT',
    'Standalone direct BRT'
  )) |>
  mutate(domain = if_else(inside_envelope, 'Inside envelope', 'Outside envelope')) |>
  group_by(candidate, scheme, domain) |>
  group_modify(~ metric_summary(.x)) |>
  ungroup()

write_csv(all_predictions, file.path(out_dir, 'cv_predictions.csv'))
write_csv(model_comparison, file.path(out_dir, 'model_comparison.csv'))
write_csv(event_metrics, file.path(out_dir, 'event_metrics.csv'))
write_csv(envelope_metrics, file.path(out_dir, 'envelope_metrics.csv'))
write_csv(tuning_results, file.path(out_dir, 'nested_tuning.csv'))

base_primary <- model_comparison |>
  filter(candidate == 'INLA explicit hazards',
         scheme == 'leave_one_event_out')
base_reef <- model_comparison |>
  filter(candidate == 'INLA explicit hazards',
         scheme == 'reef_blocked_5fold')
base_events <- event_metrics |>
  filter(candidate == 'INLA explicit hazards', event_year %in% c(2020, 2022)) |>
  select(event_year, base_event_rmse = rmse)
promotion <- model_comparison |>
  filter(
    scheme == 'leave_one_event_out',
    candidate != 'INLA explicit hazards'
  ) |>
  select(candidate, rmse, severe_rmse, false_extreme_rate) |>
  left_join(
    model_comparison |>
      filter(scheme == 'reef_blocked_5fold') |>
      select(candidate, reef_rmse = rmse, reef_severe_rmse = severe_rmse,
             reef_false_extreme_rate = false_extreme_rate),
    by = 'candidate'
  ) |>
  left_join(
    event_metrics |>
      filter(event_year %in% c(2020, 2022)) |>
      select(candidate, event_year, event_rmse = rmse) |>
      left_join(base_events, by = 'event_year') |>
      mutate(event_delta = event_rmse - base_event_rmse) |>
      group_by(candidate) |>
      summarise(max_lanina_rmse_increase = max(event_delta), .groups = 'drop'),
    by = 'candidate'
  ) |>
  mutate(
    improves_event_rmse = rmse < base_primary$rmse,
    improves_event_severe = severe_rmse < base_primary$severe_rmse,
    false_extremes_safe = false_extreme_rate <=
      base_primary$false_extreme_rate + .005,
    reef_transfer_safe = reef_rmse <= base_reef$rmse + .002,
    lanina_safe = max_lanina_rmse_increase <= .005,
    passes = improves_event_rmse & improves_event_severe &
      false_extremes_safe & reef_transfer_safe & lanina_safe
  ) |>
  arrange(desc(passes), rmse)
selected_candidate <- if (any(promotion$passes)) {
  promotion$candidate[[which(promotion$passes)[[1]]]]
} else 'INLA explicit hazards'
write_csv(promotion, file.path(out_dir, 'promotion_assessment.csv'))
write_lines(selected_candidate, file.path(out_dir, 'selected_candidate.txt'))

severe_controls <- all_predictions |>
  filter(scheme == 'leave_one_event_out') |>
  group_by(candidate, ReefID, ReefName, event_year) |>
  summarise(
    observed_mortality = mean(observed_mortality),
    predicted_mortality = mean(predicted_mortality),
    residual = observed_mortality - predicted_mortality,
    applicability = mean(applicability),
    inside_envelope = all(inside_envelope), .groups = 'drop'
  ) |>
  group_by(ReefID, event_year) |>
  mutate(max_abs_residual = max(abs(residual))) |>
  ungroup() |>
  filter(observed_mortality >= .5 | max_abs_residual >= .3 |
           str_detect(str_to_lower(ReefName), 'gannett|mackay|penrith'))
write_csv(severe_controls, file.path(out_dir, 'severe_control_audit.csv'))

# Full-data diagnostic fit. This does not enter held-out model selection.
direct_setting_id <- tuning_results |>
  filter(selected_direct) |>
  count(setting_id, sort = TRUE) |>
  slice(1) |>
  pull(setting_id)
two_part_setting_id <- tuning_results |>
  filter(selected_two_part) |>
  count(setting_id, sort = TRUE) |>
  slice(1) |>
  pull(setting_id)
direct_setting <- tuning_grid |> filter(setting_id == direct_setting_id)
two_part_setting <- tuning_grid |> filter(setting_id == two_part_setting_id)
full_direct <- fit_two_part(rows, rows, direct_setting, 20263903L)
full_two_part <- if (direct_setting_id == two_part_setting_id) {
  full_direct
} else fit_two_part(rows, rows, two_part_setting, 20264903L)
saveRDS(list(
  direct = full_direct$model$direct,
  occurrence = full_two_part$model$occurrence,
  magnitude = full_two_part$model$magnitude,
  direct_setting = direct_setting,
  two_part_setting = two_part_setting,
  predictors = predictors
), file.path(out_dir, 'full_brt_models.rds'))

component_importance <- function(component, component_name, n_trees) {
  if (component$type == 'constant') return(tibble())
  summary(component$model, n.trees = n_trees, plotit = FALSE) |>
    as_tibble() |>
    transmute(component = component_name, variable = var,
              relative_influence = rel.inf)
}
importance <- bind_rows(
  component_importance(full_direct$model$direct, 'Direct mortality',
                       direct_setting$n_trees),
  component_importance(full_two_part$model$occurrence, 'Occurrence',
                       two_part_setting$n_trees),
  component_importance(full_two_part$model$magnitude, 'Positive magnitude',
                       two_part_setting$n_trees)
)
write_csv(importance, file.path(out_dir, 'variable_importance.csv'))

direct_top <- importance |>
  filter(component == 'Direct mortality', variable %in% retained_numeric) |>
  arrange(desc(relative_influence)) |>
  slice_head(n = 9) |>
  pull(variable)
grid_size <- 35L
bootstrap_n <- 25L
grids <- setNames(lapply(direct_top, function(variable) {
  bounds <- quantile(rows[[variable]], c(.02, .98), na.rm = TRUE)
  seq(bounds[[1]], bounds[[2]], length.out = grid_size)
}), direct_top)
bootstrap_values <- setNames(lapply(direct_top, function(variable) {
  matrix(NA_real_, nrow = bootstrap_n, ncol = grid_size)
}), direct_top)
cluster_ids <- unique(rows$reef_event_key)
for (bootstrap_index in seq_len(bootstrap_n)) {
  sampled <- sample(cluster_ids, length(cluster_ids), replace = TRUE)
  bootstrap_rows <- bind_rows(lapply(seq_along(sampled), function(copy) {
    rows |>
      filter(reef_event_key == sampled[[copy]]) |>
      mutate(reef_event_key = paste0(reef_event_key, '__', copy))
  }))
  fitted <- fit_two_part(
    bootstrap_rows, rows, direct_setting,
    seed = 20265903L + bootstrap_index
  )
  for (variable in direct_top) {
    for (grid_index in seq_along(grids[[variable]])) {
      newdata <- fitted$assessment
      newdata[[variable]] <- grids[[variable]][[grid_index]]
      bootstrap_values[[variable]][bootstrap_index, grid_index] <- mean(
        pmin(pmax(predict_component(
          fitted$model$direct, newdata, direct_setting$n_trees
        ), 0), 1)
      )
    }
  }
}
pdp <- map_dfr(direct_top, function(variable) {
  values <- bootstrap_values[[variable]]
  tibble(
    variable = variable,
    predictor_value = grids[[variable]],
    median = apply(values, 2, median, na.rm = TRUE),
    lower95 = apply(values, 2, quantile, .025, na.rm = TRUE),
    upper95 = apply(values, 2, quantile, .975, na.rm = TRUE),
    bootstrap_fits = bootstrap_n
  )
})
write_csv(pdp, file.path(out_dir, 'bootstrapped_partial_dependence.csv'))

candidate_levels <- c(
  'INLA explicit hazards', 'Standalone direct BRT',
  'Standalone two-part BRT', 'Direct smooth 50% ensemble',
  'Two-part smooth 50% ensemble', 'Direct applicability champion',
  'Two-part applicability champion', 'Direct hard envelope champion',
  'Two-part hard envelope champion'
)
validation_plot_data <- model_comparison |>
  filter(candidate %in% candidate_levels) |>
  select(candidate, scheme, rmse, severe_rmse, predictive_r2,
         false_extreme_rate) |>
  pivot_longer(c(rmse, severe_rmse, predictive_r2, false_extreme_rate),
               names_to = 'metric', values_to = 'value') |>
  mutate(
    candidate = factor(candidate, levels = rev(candidate_levels)),
    metric = recode(metric, rmse = 'RMSE', severe_rmse = 'Severe RMSE',
                    predictive_r2 = 'Predictive R2',
                    false_extreme_rate = 'False-extreme rate')
  )
p_validation <- ggplot(
  validation_plot_data, aes(value, candidate, colour = scheme)
) +
  geom_point(size = 2.7, position = position_dodge(width = .55)) +
  facet_wrap(~metric, scales = 'free_x') +
  labs(
    title = 'Standalone BRT and environmental-envelope ensemble validation',
    subtitle = 'Identical event-held-out and reef-blocked assessments',
    x = NULL, y = NULL, colour = 'Validation'
  ) + theme_bw(base_size = 10) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_validation, 'Fig-FULLBRT-01_validation', validation_plot_data,
  'Matched validation of the operational INLA composite, direct and two-part standalone BRTs, and pre-specified environmental-envelope ensembles.',
  'Separates new-event transfer from interpolation among reefs in represented events.',
  'Envelope weights use predictors only; none of the BRT candidates passes the event-held-out promotion safeguards unless explicitly stated.',
  selected_id, 'INLA-BRT comparison', 'model validation',
  'operational_candidate', root, TRUE,
  code_source = 'scripts/fit_standalone_brt_envelope_ensemble.R',
  width = 13, height = 8
)

event_plot_data <- event_metrics |>
  filter(candidate %in% c(
    'INLA explicit hazards', 'Standalone direct BRT',
    'Standalone two-part BRT', 'Direct smooth 50% ensemble',
    'Two-part smooth 50% ensemble'
  )) |>
  mutate(event_year = factor(event_year))
p_event <- ggplot(event_plot_data,
                  aes(event_year, rmse, colour = candidate, group = candidate)) +
  geom_line(linewidth = .8) + geom_point(size = 2.2) +
  labs(
    title = 'New-event transfer by bleaching year',
    subtitle = 'Leave-one-event-out RMSE; lower is better',
    x = 'Held-out event', y = 'RMSE', colour = 'Candidate'
  ) + theme_bw(base_size = 10) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_event, 'Fig-FULLBRT-02_event_transfer', event_plot_data,
  'Event-specific RMSE for operational INLA, standalone BRTs and smooth envelope ensembles.',
  'Shows whether interpolation gains survive the distinct environmental states of 2016, 2017, 2020, 2022 and 2024.',
  'Only five events are available; event identity is never used as an operational predictor.',
  selected_id, 'INLA-BRT comparison', 'event transfer',
  'operational_candidate', root, TRUE,
  code_source = 'scripts/fit_standalone_brt_envelope_ensemble.R',
  width = 11, height = 7
)

domain_plot_data <- envelope_metrics |>
  filter(candidate %in% c('INLA explicit hazards',
                          'Standalone direct BRT',
                          'Standalone two-part BRT'))
p_domain <- ggplot(
  domain_plot_data, aes(candidate, rmse, fill = candidate)
) + geom_col(show.legend = FALSE) +
  facet_grid(domain ~ scheme, scales = 'free_y') +
  labs(
    title = 'Performance inside and outside the environmental envelope',
    subtitle = 'Envelope is learned from training predictors only',
    x = NULL, y = 'RMSE'
  ) + theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_figure_bundle(
  p_domain, 'Fig-FULLBRT-03_envelope_performance', domain_plot_data,
  'INLA and standalone-BRT RMSE inside and outside each fold-specific multivariate environmental envelope.',
  'Directly tests the hypothesis that BRT is the stronger learner within supported environmental space.',
  'Applicability is relative to each training fold and does not imply ecological causality.',
  selected_id, 'INLA-BRT comparison', 'applicability domain',
  'operational_candidate', root, TRUE,
  code_source = 'scripts/fit_standalone_brt_envelope_ensemble.R',
  width = 12, height = 8
)

calibration_data <- all_predictions |>
  filter(
    scheme == 'leave_one_event_out',
    candidate %in% c('INLA explicit hazards', 'Standalone direct BRT',
                     'Standalone two-part BRT',
                     'Direct smooth 50% ensemble')
  )
p_calibration <- ggplot(
  calibration_data, aes(observed_mortality, predicted_mortality)
) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(colour = factor(event_year)), alpha = .55, size = 1.3) +
  geom_smooth(method = 'loess', se = FALSE, colour = 'black', linewidth = .7) +
  facet_wrap(~candidate, ncol = 2) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = 'Leave-one-event-out calibration',
    subtitle = 'Independent BRT predictions, not residual corrections',
    x = 'Observed mortality', y = 'Predicted mortality', colour = 'Event'
  ) + theme_bw(base_size = 10) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_calibration, 'Fig-FULLBRT-04_calibration', calibration_data,
  'Observed versus predicted mortality under leave-one-event-out validation for INLA, standalone BRTs and the best fixed smooth envelope blend.',
  'Displays mortality-floor and severe-tail calibration that aggregate metrics can obscure.',
  'Loess lines are descriptive and are not used to recalibrate predictions.',
  selected_id, 'INLA-BRT comparison', 'calibration',
  'operational_candidate', root, TRUE,
  code_source = 'scripts/fit_standalone_brt_envelope_ensemble.R',
  width = 11, height = 10
)

applicability_data <- all_predictions |>
  filter(candidate == 'INLA explicit hazards') |>
  mutate(
    direct_error_gain = abs(observed_mortality - inla_prediction) -
      abs(observed_mortality - direct_brt_prediction),
    two_part_error_gain = abs(observed_mortality - inla_prediction) -
      abs(observed_mortality - brt_prediction)
  ) |>
  select(scheme, source_observation_id, applicability,
         direct_error_gain, two_part_error_gain) |>
  pivot_longer(c(direct_error_gain, two_part_error_gain),
               names_to = 'brt', values_to = 'absolute_error_gain') |>
  mutate(brt = recode(brt,
    direct_error_gain = 'Direct BRT',
    two_part_error_gain = 'Two-part BRT'))
p_applicability <- ggplot(
  applicability_data, aes(applicability, absolute_error_gain, colour = brt)
) + geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = .25, size = 1) +
  geom_smooth(method = 'loess', se = TRUE, linewidth = .9) +
  facet_wrap(~scheme) +
  labs(
    title = 'Does BRT gain skill as environmental applicability increases?',
    subtitle = 'Positive values mean lower absolute error than INLA',
    x = 'Training-envelope applicability',
    y = 'INLA absolute error minus BRT absolute error', colour = 'BRT'
  ) + theme_bw(base_size = 10) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_applicability, 'Fig-FULLBRT-05_local_competence', applicability_data,
  'Observation-level BRT absolute-error gain over INLA against training-fold environmental applicability.',
  'Tests whether a denser supported environment is sufficient to justify greater BRT ensemble weight.',
  'The smoother is diagnostic only; the evaluated envelope weights are outcome-free and pre-specified.',
  selected_id, 'INLA-BRT comparison', 'local competence',
  'operational_candidate', root, TRUE,
  code_source = 'scripts/fit_standalone_brt_envelope_ensemble.R',
  width = 11, height = 7
)

labels <- c(
  ann_maxdhw = 'Local-first maximum DHW',
  prop_acropora_pre = 'Pre-event Acropora',
  observed_pre_cover = 'Starting coral cover',
  dhw_novelty10 = 'DHW novelty',
  dhw_events_since2016_n6 = 'Events above 6 DHW since 2016',
  secc3m_p10 = 'Low-tail Secchi depth', cloudp_90 = 'High cloud cover',
  mcur_90 = 'High current speed', sst_summer_skewness = 'SST skewness',
  log_coastal_rain30 = 'Coastal rainfall',
  wqc_freqcc12 = 'Coloured-water frequency',
  wqc_prior10_percentile = 'Relative coloured-water exposure',
  era5_wind_calm_fraction = 'Calm-wind fraction',
  log1p_cyclone_wave_hours = 'Damaging-wave hours',
  tc_interval_wind_distance_index = 'Cyclone wind-distance index',
  cots_outbreak_probability = 'COTS outbreak probability',
  rrn_event_excess_raw = 'RRN COTS excess',
  manta_latest_excess_raw = 'Latest Manta COTS excess',
  manta_peak_excess_raw = 'Prior Manta COTS peak',
  manta_years_since_peak = 'Years since Manta peak',
  disease_risk_max = 'NOAA disease risk', depth = 'Survey depth',
  programme_factor = 'Survey programme', region_factor = 'GBR region'
)
importance_plot_data <- importance |>
  group_by(component) |>
  slice_max(relative_influence, n = 12, with_ties = FALSE) |>
  ungroup() |>
  mutate(label = coalesce(unname(labels[variable]), variable),
         label = reorder(label, relative_influence))
p_importance <- ggplot(
  importance_plot_data, aes(relative_influence, label)
) + geom_col(fill = '#2A9D8F') + facet_wrap(~component, scales = 'free_y') +
  labs(
    title = 'Full standalone BRT variable importance',
    subtitle = 'Direct mortality and two-part occurrence/magnitude components',
    x = 'Relative influence (%)', y = NULL
  ) + theme_bw(base_size = 10)
save_figure_bundle(
  p_importance, 'Fig-FULLBRT-06_variable_importance', importance_plot_data,
  'Relative predictor influence for the full-data direct BRT and two-part occurrence and positive-magnitude components.',
  'Shows how the independently fitted BRT frameworks allocate split improvement across operational predictors.',
  'Full-data importance is descriptive; correlated mechanisms can share influence and held-out validation governs model promotion.',
  selected_id, 'Standalone BRT', 'variable importance',
  'operational_diagnostic', root, TRUE,
  code_source = 'scripts/fit_standalone_brt_envelope_ensemble.R',
  width = 13, height = 10
)

pdp_plot_data <- pdp |>
  mutate(label = coalesce(unname(labels[variable]), variable))
p_pdp <- ggplot(pdp_plot_data, aes(predictor_value, median)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95),
              fill = '#56B4E9', alpha = .28) +
  geom_line(colour = '#0072B2', linewidth = .9) +
  facet_wrap(~label, scales = 'free_x', ncol = 3) +
  labs(
    title = 'Direct BRT partial dependence with 95% bootstrap intervals',
    subtitle = paste0('Reef-event cluster bootstrap; B = ', bootstrap_n),
    x = 'Predictor value', y = 'Predicted relative mortality'
  ) + theme_bw(base_size = 10)
save_figure_bundle(
  p_pdp, 'Fig-FULLBRT-07_bootstrap_pdp', pdp_plot_data,
  'Mortality-scale partial-dependence curves with 95% reef-event cluster-bootstrap intervals for the nine most influential direct-BRT predictors.',
  'Provides visual checks of nonlinear thresholds, slopes and unstable tails in the independently fitted direct BRT.',
  'These are marginal full-data diagnostics, not prediction intervals or causal effects.',
  selected_id, 'Standalone direct BRT', 'partial dependence',
  'operational_diagnostic', root, TRUE,
  code_source = 'scripts/fit_standalone_brt_envelope_ensemble.R',
  width = 13, height = 10
)

severe_plot_data <- severe_controls |>
  filter(candidate %in% c(
    'INLA explicit hazards', 'Standalone direct BRT',
    'Standalone two-part BRT', 'Direct smooth 50% ensemble'
  )) |>
  mutate(label = paste0(ReefName, ' (', event_year, ')')) |>
  group_by(label) |>
  mutate(priority = max(abs(residual))) |>
  ungroup() |>
  arrange(desc(priority)) |>
  filter(dense_rank(desc(priority)) <= 16) |>
  mutate(label = factor(label, levels = rev(unique(label))))
p_severe <- ggplot(
  severe_plot_data, aes(predicted_mortality, label, colour = candidate)
) + geom_point(size = 2.3, position = position_dodge(width = .55)) +
  geom_point(aes(x = observed_mortality), colour = 'black', shape = 4,
             size = 2.8, stroke = 1) +
  labs(
    title = 'Severe and focal reef-event controls',
    subtitle = 'Black crosses are observed mortality',
    x = 'Held-out mortality', y = NULL, colour = 'Candidate'
  ) + theme_bw(base_size = 10) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_severe, 'Fig-FULLBRT-08_severe_controls', severe_plot_data,
  'Held-out predictions at severe and high-residual reef-events for INLA, standalone BRTs and the best fixed smooth envelope blend.',
  'Shows whether an aggregate validation change addresses Gannett, Mackay, Penrith and other ecological controls.',
  'Repeated programme/depth observations are summarised to reef-event means.',
  selected_id, 'INLA-BRT comparison', 'severe controls',
  'operational_candidate', root, TRUE,
  code_source = 'scripts/fit_standalone_brt_envelope_ensemble.R',
  width = 13, height = 9
)
write_figure_readme(root)
message('Standalone BRT and environmental-envelope outputs written to: ',
        out_dir)
