# Standalone occurrence diagnostic for raw versus local-first DHW.
#
# This deliberately does not modify or select the production model. The
# binomial response is whether relative mortality is greater than zero; survey
# mortality proportions are not treated as literal colony-level trial counts.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

set.seed(20260911L)

event_years <- c(2016L, 2017L, 2020L, 2022L, 2024L, 2025L)
programme_keys <- c('manta', 'ltmp', 'mmp')
output_dir <- file.path('output', 'bleaching_only_dhw_binomial')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

historical <- map_dfr(programme_keys, function(programme_key) {
  readRDS(file.path(
    'data', 'processed',
    paste0('validation_rows_', programme_key, '.rds')
  ))
})
forward_2025 <- map_dfr(programme_keys, function(programme_key) {
  path <- file.path(
    'data', 'processed',
    paste0('initial_forecast_rows_2025_', programme_key, '.rds')
  )
  if (!file.exists(path)) {
    stop('Missing locked 2025 assessment rows: ', path)
  }
  readRDS(path)
})

correction_path <- file.path(
  'data', 'processed',
  'noaa_dhw_correction_layer_local_first_validation.csv'
)
correction <- read_csv(correction_path, show_col_types = FALSE) |>
  select(
    ReefID, event_year, local_first_correction, correction_source,
    correction_sd, nearest_local_logger_km, effective_local_loggers,
    local_loggers_used
  )

all_rows <- bind_rows(historical, forward_2025) |>
  filter(event_year %in% .env$event_years) |>
  mutate(
    explicit_nonthermal =
      coalesce(disturbance_has_cots, FALSE) |
      coalesce(disturbance_has_cyclone, FALSE) |
      coalesce(disturbance_has_flood, FALSE),
    attribution_class = case_when(
      explicit_nonthermal ~ 'explicit_nonthermal_excluded',
      coalesce(disturbance_has_bleaching, FALSE) ~ 'explicit_bleaching',
      TRUE ~ 'unlabelled_no_explicit_nonthermal'
    )
  ) |>
  left_join(
    correction,
    by = c('ReefID', 'event_year'),
    relationship = 'many-to-one'
  ) |>
  mutate(
    raw_dhw = ann_maxdhw,
    local_first_correction = coalesce(local_first_correction, 0),
    correction_source = case_when(
      event_year == 2025L & is.na(correction_source) ~
        'unavailable_2025_no_adjustment',
      is.na(correction_source) ~ 'not_in_correction_grid_zero',
      TRUE ~ correction_source
    ),
    logger_adjusted_dhw = pmax(raw_dhw + local_first_correction, 0),
    has_mortality = mortality_prop > 0
  )

write_csv(
  all_rows |>
    count(event_year, programme_key, attribution_class, name = 'observations'),
  file.path(output_dir, 'eligibility_audit.csv')
)

rows <- all_rows |>
  filter(!explicit_nonthermal) |>
  mutate(.diagnostic_row_id = row_number())

if (any(!is.finite(rows$raw_dhw)) || any(!is.finite(rows$mortality_prop))) {
  stop('Eligible rows contain missing DHW or mortality outcomes')
}
if (n_distinct(rows$event_year) != length(event_years)) {
  stop('Expected all six bleaching events in the diagnostic data')
}

# The canonical reef folds were constructed separately by programme. This
# pooled diagnostic needs a single common assignment so the same reef cannot
# leak through another programme. Allocate reefs greedily by row count with
# ReefID as the deterministic tie-break, and save the resulting assignment.
reef_fold_lookup <- rows |>
  count(ReefID, name = 'observations') |>
  arrange(desc(observations), ReefID)
fold_load <- rep(0L, 5L)
reef_fold_lookup$reef_fold <- NA_integer_
for (reef_index in seq_len(nrow(reef_fold_lookup))) {
  selected_fold <- which(fold_load == min(fold_load))[1]
  reef_fold_lookup$reef_fold[[reef_index]] <- selected_fold
  fold_load[[selected_fold]] <-
    fold_load[[selected_fold]] + reef_fold_lookup$observations[[reef_index]]
}
reef_fold_lookup <- reef_fold_lookup |> select(ReefID, reef_fold)
write_csv(
  reef_fold_lookup,
  file.path(output_dir, 'pooled_reef_fold_assignments.csv')
)
rows <- rows |>
  select(-any_of('reef_fold')) |>
  left_join(reef_fold_lookup, by = 'ReefID', relationship = 'many-to-one')
stopifnot(all(rows$reef_fold %in% seq_len(5L)))

variant_table <- tribble(
  ~variant, ~dhw_column, ~label,
  'raw_noaa_dhw', 'raw_dhw', 'Raw NOAA DHW',
  'local_first_logger_adjusted_dhw', 'logger_adjusted_dhw',
  'Local-first logger-adjusted DHW'
)

fold_table <- bind_rows(
  tibble(scheme = 'apparent_full_data', fold = 'all'),
  tibble(scheme = 'leave_one_event_out', fold = as.character(event_years)),
  tibble(
    scheme = 'leave_one_sector_out',
    fold = sort(unique(rows$SECTOR))
  ),
  tibble(scheme = 'reef_blocked_5fold', fold = as.character(seq_len(5L))),
  tibble(
    scheme = 'leave_one_programme_out',
    fold = sort(unique(rows$programme_key))
  )
)

assessment_mask <- function(data, scheme, fold) {
  switch(
    scheme,
    apparent_full_data = rep(TRUE, nrow(data)),
    leave_one_event_out = data$event_year == as.integer(fold),
    leave_one_sector_out = data$SECTOR == fold,
    reef_blocked_5fold = data$reef_fold == as.integer(fold),
    leave_one_programme_out = data$programme_key == fold,
    stop('Unknown validation scheme: ', scheme)
  )
}

fit_binomial <- function(analysis, assessment, dhw_column) {
  analysis <- analysis |>
    mutate(dhw_input = .data[[dhw_column]])
  assessment <- assessment |>
    mutate(dhw_input = .data[[dhw_column]])
  if (n_distinct(analysis$has_mortality) < 2L) {
    stop('Training fold contains only one outcome class')
  }
  fit <- suppressWarnings(glm(
    has_mortality ~ dhw_input,
    data = analysis,
    family = binomial(link = 'logit')
  ))
  probability <- as.numeric(predict(
    fit, newdata = assessment, type = 'response'
  ))
  if (any(!is.finite(probability))) {
    stop('Non-finite binomial predictions')
  }
  list(fit = fit, probability = pmin(pmax(probability, 1e-8), 1 - 1e-8))
}

prediction_parts <- list()
fit_parts <- list()
part_index <- 0L
for (fold_index in seq_len(nrow(fold_table))) {
  scheme <- fold_table$scheme[[fold_index]]
  fold <- fold_table$fold[[fold_index]]
  held_out <- assessment_mask(rows, scheme, fold)
  assessment <- rows[held_out, , drop = FALSE]
  analysis <- if (identical(scheme, 'apparent_full_data')) {
    rows
  } else {
    rows[!held_out, , drop = FALSE]
  }
  for (variant_index in seq_len(nrow(variant_table))) {
    variant <- variant_table$variant[[variant_index]]
    dhw_column <- variant_table$dhw_column[[variant_index]]
    result <- fit_binomial(analysis, assessment, dhw_column)
    coefficients <- coef(summary(result$fit))
    part_index <- part_index + 1L
    prediction_parts[[part_index]] <- assessment |>
      transmute(
        .diagnostic_row_id, source_observation_id, ReefID, ReefName,
        programme_key, SECTOR, event_year, mortality_prop, has_mortality,
        raw_dhw, local_first_correction, logger_adjusted_dhw,
        correction_source, attribution_class, scheme, fold, variant,
        predicted_probability = result$probability
      )
    fit_parts[[part_index]] <- tibble(
      scheme, fold, variant,
      analysis_n = nrow(analysis), assessment_n = nrow(assessment),
      analysis_events = n_distinct(analysis$event_year),
      analysis_reefs = n_distinct(analysis$ReefID),
      converged = isTRUE(result$fit$converged),
      intercept = unname(coef(result$fit)[['(Intercept)']]),
      dhw_slope = unname(coef(result$fit)[['dhw_input']]),
      dhw_slope_se = coefficients['dhw_input', 'Std. Error'],
      aic = AIC(result$fit)
    )
  }
}
predictions <- bind_rows(prediction_parts)
fit_audit <- bind_rows(fit_parts)

binary_auc <- function(observed, probability) {
  positive <- observed == 1
  n_positive <- sum(positive)
  n_negative <- sum(!positive)
  if (n_positive == 0L || n_negative == 0L) return(NA_real_)
  ranks <- rank(probability, ties.method = 'average')
  (sum(ranks[positive]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

metric_row <- function(data) {
  observed <- as.integer(data$has_mortality)
  probability <- data$predicted_probability
  tibble(
    n = nrow(data),
    reefs = n_distinct(data$ReefID),
    events = n_distinct(data$event_year),
    observed_occurrence = mean(observed),
    predicted_occurrence = mean(probability),
    bias_predicted_minus_observed = mean(probability - observed),
    brier = mean((probability - observed)^2),
    rmse = sqrt(mean((probability - observed)^2)),
    mae = mean(abs(probability - observed)),
    predictive_r2 = 1 - sum((probability - observed)^2) /
      sum((observed - mean(observed))^2),
    log_loss = -mean(
      observed * log(probability) + (1 - observed) * log(1 - probability)
    ),
    auc = binary_auc(observed, probability),
    false_positive_rate_050 = if (any(observed == 0)) {
      mean(probability[observed == 0] >= 0.5)
    } else {
      NA_real_
    }
  )
}

fold_metrics <- predictions |>
  group_by(scheme, fold, variant) |>
  group_modify(~ metric_row(.x)) |>
  ungroup()
overall_metrics <- predictions |>
  group_by(scheme, variant) |>
  group_modify(~ metric_row(.x)) |>
  ungroup()
paired_comparison <- overall_metrics |>
  select(scheme, variant, brier, rmse, mae, log_loss, auc) |>
  pivot_wider(
    names_from = variant,
    values_from = c(brier, rmse, mae, log_loss, auc)
  ) |>
  mutate(
    delta_brier_adjusted_minus_raw =
      brier_local_first_logger_adjusted_dhw - brier_raw_noaa_dhw,
    delta_rmse_adjusted_minus_raw =
      rmse_local_first_logger_adjusted_dhw - rmse_raw_noaa_dhw,
    delta_log_loss_adjusted_minus_raw =
      log_loss_local_first_logger_adjusted_dhw - log_loss_raw_noaa_dhw,
    delta_auc_adjusted_minus_raw =
      auc_local_first_logger_adjusted_dhw - auc_raw_noaa_dhw
  )

# Retrospective event-specific curves. Each event and DHW representation gets
# its own intercept and slope. These describe the observed event; they are not
# event-held-out forecasts and event identity is never used as a predictor.
curve_grid <- seq(0, 20, by = 0.1)
event_curve_parts <- list()
event_coefficient_parts <- list()
curve_index <- 0L
for (event_year_value in event_years) {
  event_rows <- rows |> filter(event_year == .env$event_year_value)
  for (variant_index in seq_len(nrow(variant_table))) {
    variant <- variant_table$variant[[variant_index]]
    dhw_column <- variant_table$dhw_column[[variant_index]]
    event_rows_fit <- event_rows |>
      mutate(dhw_input = .data[[dhw_column]])
    fit <- suppressWarnings(glm(
      has_mortality ~ dhw_input,
      data = event_rows_fit,
      family = binomial(link = 'logit')
    ))
    link <- predict(
      fit,
      newdata = tibble(dhw_input = curve_grid),
      type = 'link', se.fit = TRUE
    )
    curve_index <- curve_index + 1L
    event_curve_parts[[curve_index]] <- tibble(
      event_year = event_year_value,
      variant,
      dhw = curve_grid,
      predicted_probability = plogis(as.numeric(link$fit)),
      lower95 = plogis(as.numeric(link$fit) - 1.96 * as.numeric(link$se.fit)),
      upper95 = plogis(as.numeric(link$fit) + 1.96 * as.numeric(link$se.fit)),
      observed_dhw_min = min(event_rows_fit$dhw_input),
      observed_dhw_max = max(event_rows_fit$dhw_input)
    )
    coefficient_values <- coef(fit)
    coefficient_summary <- coef(summary(fit))
    slope <- unname(coefficient_values[['dhw_input']])
    slope_se <- coefficient_summary['dhw_input', 'Std. Error']
    slope_p_value <- coefficient_summary['dhw_input', 'Pr(>|z|)']
    intercept <- unname(coefficient_values[['(Intercept)']])
    event_coefficient_parts[[curve_index]] <- tibble(
      event_year = event_year_value,
      variant,
      n = nrow(event_rows_fit),
      mortality_occurrences = sum(event_rows_fit$has_mortality),
      corrected_rows = sum(event_rows_fit$local_first_correction != 0),
      correction_coverage = mean(event_rows_fit$local_first_correction != 0),
      intercept,
      dhw_slope = slope,
      dhw_slope_se = slope_se,
      dhw_slope_lower95 = slope - 1.96 * slope_se,
      dhw_slope_upper95 = slope + 1.96 * slope_se,
      dhw_slope_p_value = slope_p_value,
      odds_ratio_per_dhw = exp(slope),
      probability_050_dhw = if (abs(slope) > 1e-10) -intercept / slope else NA_real_,
      transition_width_p10_to_p90 = if (abs(slope) > 1e-10) {
        2 * qlogis(0.9) / abs(slope)
      } else {
        NA_real_
      },
      converged = isTRUE(fit$converged)
    )
  }
}
event_curves <- bind_rows(event_curve_parts)
event_coefficients <- bind_rows(event_coefficient_parts)
event_shape_comparison <- event_coefficients |>
  select(
    event_year, variant, n, mortality_occurrences, corrected_rows,
    correction_coverage, intercept, dhw_slope, odds_ratio_per_dhw,
    dhw_slope_se, dhw_slope_lower95, dhw_slope_upper95, dhw_slope_p_value,
    probability_050_dhw, transition_width_p10_to_p90
  ) |>
  pivot_wider(
    names_from = variant,
    values_from = c(
      n, mortality_occurrences, corrected_rows, correction_coverage,
      intercept, dhw_slope, odds_ratio_per_dhw, probability_050_dhw,
      dhw_slope_se, dhw_slope_lower95, dhw_slope_upper95, dhw_slope_p_value,
      transition_width_p10_to_p90
    )
  ) |>
  mutate(
    delta_slope_adjusted_minus_raw =
      dhw_slope_local_first_logger_adjusted_dhw - dhw_slope_raw_noaa_dhw,
    slope_ratio_adjusted_over_raw =
      dhw_slope_local_first_logger_adjusted_dhw / dhw_slope_raw_noaa_dhw,
    delta_midpoint_dhw_adjusted_minus_raw =
      probability_050_dhw_local_first_logger_adjusted_dhw -
      probability_050_dhw_raw_noaa_dhw
  )

coverage <- rows |>
  group_by(event_year) |>
  summarise(
    observations = n(), reefs = n_distinct(ReefID),
    mortality_occurrences = sum(has_mortality),
    explicit_bleaching_rows = sum(attribution_class == 'explicit_bleaching'),
    unlabelled_compatible_rows =
      sum(attribution_class == 'unlabelled_no_explicit_nonthermal'),
    nonzero_corrections = sum(local_first_correction != 0),
    nonzero_correction_fraction = mean(local_first_correction != 0),
    negative_corrections = sum(local_first_correction < 0),
    positive_corrections = sum(local_first_correction > 0),
    correction_median = median(local_first_correction),
    correction_q10 = quantile(local_first_correction, 0.10),
    correction_q90 = quantile(local_first_correction, 0.90),
    raw_dhw_min = min(raw_dhw), raw_dhw_max = max(raw_dhw),
    adjusted_dhw_min = min(logger_adjusted_dhw),
    adjusted_dhw_max = max(logger_adjusted_dhw),
    .groups = 'drop'
  )

focal_probability_contrasts <- event_curves |>
  filter(
    event_year %in% c(2016L, 2024L),
    round(dhw, 1) %in% c(4, 8, 12)
  ) |>
  select(event_year, dhw, variant, predicted_probability, lower95, upper95) |>
  pivot_wider(
    names_from = variant,
    values_from = c(predicted_probability, lower95, upper95)
  ) |>
  mutate(
    delta_probability_adjusted_minus_raw =
      predicted_probability_local_first_logger_adjusted_dhw -
      predicted_probability_raw_noaa_dhw
  )

write_csv(predictions, file.path(output_dir, 'predictions.csv'))
write_csv(fit_audit, file.path(output_dir, 'fit_audit.csv'))
write_csv(fold_metrics, file.path(output_dir, 'fold_metrics.csv'))
write_csv(overall_metrics, file.path(output_dir, 'overall_metrics.csv'))
write_csv(paired_comparison, file.path(output_dir, 'paired_comparison.csv'))
write_csv(event_curves, file.path(output_dir, 'event_curves.csv'))
write_csv(event_coefficients, file.path(output_dir, 'event_coefficients.csv'))
write_csv(
  event_shape_comparison,
  file.path(output_dir, 'event_shape_comparison.csv')
)
write_csv(coverage, file.path(output_dir, 'coverage.csv'))
write_csv(
  focal_probability_contrasts,
  file.path(output_dir, 'focal_2016_2024_probability_contrasts.csv')
)

curve_plot <- event_curves |>
  left_join(
    variant_table |> select(variant, label),
    by = 'variant', relationship = 'many-to-one'
  ) |>
  ggplot(aes(dhw, predicted_probability, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ event_year, ncol = 3) +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 1)) +
  scale_colour_manual(values = c(
    'Raw NOAA DHW' = '#2b6cb0',
    'Local-first logger-adjusted DHW' = '#c05621'
  )) +
  scale_fill_manual(values = c(
    'Raw NOAA DHW' = '#2b6cb0',
    'Local-first logger-adjusted DHW' = '#c05621'
  )) +
  labs(
    title = 'Event-specific mortality-occurrence curves',
    subtitle = paste(
      'Bleaching-compatible rows; separate simple binomial fit per event.',
      'Shading is a model-based 95% interval.'
    ),
    x = 'DHW supplied to the model',
    y = 'Probability of observed relative mortality > 0',
    colour = NULL, fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom')
ggsave(
  file.path(output_dir, 'event_specific_curves.png'),
  curve_plot, width = 11, height = 7.5, dpi = 300
)
ggsave(
  file.path(output_dir, 'event_specific_curves.pdf'),
  curve_plot, width = 11, height = 7.5
)

focal_observations <- rows |>
  filter(event_year %in% c(2016L, 2024L)) |>
  select(
    event_year, has_mortality, raw_dhw, logger_adjusted_dhw
  ) |>
  pivot_longer(
    c(raw_dhw, logger_adjusted_dhw),
    names_to = 'dhw_column', values_to = 'dhw_input'
  ) |>
  left_join(
    variant_table,
    by = 'dhw_column', relationship = 'many-to-one'
  )

focal_curve_data <- event_curves |>
  filter(event_year %in% c(2016L, 2024L)) |>
  left_join(
    variant_table |> select(variant, label),
    by = 'variant', relationship = 'many-to-one'
  )

focal_plot <- ggplot() +
  geom_ribbon(
    data = focal_curve_data,
    aes(dhw, ymin = lower95, ymax = upper95, fill = label),
    alpha = 0.12, colour = NA
  ) +
  geom_line(
    data = focal_curve_data,
    aes(dhw, predicted_probability, colour = label),
    linewidth = 1
  ) +
  geom_jitter(
    data = focal_observations,
    aes(dhw_input, as.numeric(has_mortality)),
    width = 0, height = 0.035, alpha = 0.45, size = 1.5
  ) +
  facet_grid(event_year ~ label) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.05, 1.05)) +
  scale_colour_manual(values = c(
    'Raw NOAA DHW' = '#2b6cb0',
    'Local-first logger-adjusted DHW' = '#c05621'
  )) +
  scale_fill_manual(values = c(
    'Raw NOAA DHW' = '#2b6cb0',
    'Local-first logger-adjusted DHW' = '#c05621'
  )) +
  labs(
    title = 'Focal 2016 and 2024 event-specific DHW curves',
    subtitle = paste(
      'Points are observed mortality occurrence at the DHW supplied to each',
      'separate event-specific binomial model.'
    ),
    x = 'DHW supplied to the model',
    y = 'Observed or fitted probability of relative mortality > 0'
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'none')
ggsave(
  file.path(output_dir, 'focal_2016_2024_curves.png'),
  focal_plot, width = 10.5, height = 7.2, dpi = 300
)
ggsave(
  file.path(output_dir, 'focal_2016_2024_curves.pdf'),
  focal_plot, width = 10.5, height = 7.2
)

message(
  'Wrote standalone raw-versus-adjusted DHW binomial diagnostic for ',
  nrow(rows), ' bleaching-compatible observations across ',
  n_distinct(rows$event_year), ' events'
)
