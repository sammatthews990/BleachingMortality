# Prospective COTS pressure nowcast for the cause-aware INLA composite.
# All monitoring/culling features are frozen at the end of February, before
# bleaching mortality is expected to arise.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(INLA)
  library(readr)
  library(readxl)
  library(stringr)
  library(tidyr)
})
Sys.setenv(INLA_ST_RUN = '0')
source('src/lib/model_registry.R')
source('src/lib/model_diagnostics.R')

root <- normalizePath('.', winslash = '/', mustWork = TRUE)
out_dir <- file.path(root, 'output', 'prospective_cots_nowcast')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
event_years <- c(2016L, 2017L, 2020L, 2022L, 2024L)

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
  occurrence_brier = mean((observed_occurrence-predicted_occurrence)^2),
  false_extreme_rate = mean(
    predicted_mortality[observed_mortality < .3] >= .3
  ), .groups = 'drop')

normalise_id <- function(x) str_to_upper(str_trim(as.character(x)))

# Cross-fitted thermal and cyclone components from the selected composite.
selected <- read_csv(
  file.path(root, 'output', 'cots_raw_enso_occurrence',
            'cots_cv_predictions.csv'), show_col_types = FALSE
) |>
  filter(candidate == 'cots_raw_interval_relative') |>
  mutate(
    source_observation_id = as.character(source_observation_id),
    ReefID = normalise_id(ReefID), fold = as.character(fold)
  )

joint_context <- read_csv(
  file.path(root, 'output', 'explanatory_event_dhw',
            'event_dhw_brt_data.csv'), show_col_types = FALSE
) |>
  mutate(source_observation_id = as.character(source_observation_id),
         ReefID = normalise_id(ReefID)) |>
  distinct(source_observation_id, .keep_all = TRUE)
reef_locations <- joint_context |> distinct(ReefID, lon, lat)

cots_hindcast <- read_csv(
  file.path(root, 'data', 'gbrPredsAdj_20262408.csv'), show_col_types = FALSE
) |>
  transmute(
    ReefName = reefName, event_year = as.integer(year),
    cots_outbreak_probability = pmin(pmax(as.numeric(outbrProb), 0), 1)
  )

annual <- read_csv(
  file.path(root, 'data', 'processed', 'annual_coral_transitions.csv'),
  show_col_types = FALSE
) |>
  mutate(
    source_observation_id = as.character(source_observation_id),
    ReefID = normalise_id(ReefID),
    relative_loss = pmin(pmax(
      -cover_change_pp / (100 * pmax(pre_cover, .02)), 0
    ), .999),
    cots_label = coalesce(disturbance_has_cots, FALSE) |
      str_detect(str_to_lower(coalesce(disturbance_text, '')),
                 'cots|crown-of-thorns'),
    cots_positive_relative = cots_label & relative_loss > .02,
    pre_cover = pmin(pmax(pre_cover, .001), 1),
    acropora = prop_acropora_pre,
    is_manta = as.numeric(programme_key == 'manta'),
    is_mmp = as.numeric(programme_key == 'mmp')
  ) |>
  left_join(cots_hindcast, by = c('ReefName', 'event_year'),
            relationship = 'many-to-one') |>
  left_join(reef_locations, by = 'ReefID', relationship = 'many-to-one') |>
  mutate(
    cots_outbreak_probability = coalesce(cots_outbreak_probability, 0),
    rrn_event_excess_raw = pmax(coalesce(cot_idwmeanpertow, 0) - .22, 0),
    log1p_rrn_event_excess = log1p(rrn_event_excess_raw)
  )

# Exact-date manta monitoring and targeted culling records.
workbook <- file.path(
  root, 'data', '250929_COTS-Manta-Cull-RHIS-Data-Matthews-and-Schlawinsky.xlsx'
)
manta <- suppressWarnings(read_excel(workbook, sheet = 'Manta')) |>
  transmute(
    ReefID = normalise_id(ReefLabel),
    survey_date = as.Date(substr(as.character(SurveyTime), 1, 10)),
    cots_count = pmax(as.numeric(CrownOfThornsStarfishCount), 0),
    tow_distance = pmax(as.numeric(TowDistance), 0)
  ) |>
  filter(!is.na(ReefID), !is.na(survey_date), is.finite(cots_count)) |>
  group_by(ReefID, survey_date) |>
  summarise(
    manta_density = mean(cots_count), manta_tows = n(),
    manta_distance_m = sum(tow_distance, na.rm = TRUE), .groups = 'drop'
  )

cull_raw <- suppressWarnings(read_excel(workbook, sheet = 'Cull'))
cull <- cull_raw |>
  transmute(
    ReefID = normalise_id(ReefLabel), survey_date = as.Date(SurveyDate),
    bottom_time = pmax(as.numeric(Bottomtime), 0),
    culled = rowSums(cbind(
      pmax(as.numeric(Cohort1), 0), pmax(as.numeric(Cohort2), 0),
      pmax(as.numeric(Cohort3), 0), pmax(as.numeric(Cohort4), 0)
    ), na.rm = TRUE)
  ) |>
  filter(!is.na(ReefID), !is.na(survey_date), is.finite(culled))

# Feature grid covers cause-model training and the five bleaching assessments.
feature_grid <- bind_rows(
  annual |> select(ReefID, event_year),
  selected |> select(ReefID, event_year)
) |>
  filter(is.finite(event_year)) |>
  distinct(ReefID, event_year) |>
  # Use 1 March as an exclusive cutoff so every February date, including
  # 29 February in leap years, is available to the nowcast.
  mutate(issue_date = as.Date(paste0(event_year, '-03-01')))

manta_split <- split(manta, manta$ReefID)
cull_split <- split(cull, cull$ReefID)

build_temporal_row <- function(reef_id, event_year, issue_date) {
  mr <- manta_split[[reef_id]]
  cr <- cull_split[[reef_id]]
  if (is.null(mr)) mr <- manta[0, ]
  if (is.null(cr)) cr <- cull[0, ]
  mr <- mr |> filter(survey_date < issue_date)
  latest_date <- if (nrow(mr)) max(mr$survey_date) else as.Date(NA)
  latest <- if (nrow(mr)) mr |> filter(survey_date == latest_date) else mr
  prior3 <- mr |> filter(survey_date >= issue_date - 3 * 365)
  peak_date <- if (nrow(prior3)) {
    prior3$survey_date[[which.max(prior3$manta_density)]]
  } else as.Date(NA)
  prior1_cull <- cr |>
    filter(survey_date < issue_date, survey_date >= issue_date - 365)
  tibble(
    ReefID = reef_id, event_year = event_year, issue_date = issue_date,
    manta_supported = as.numeric(nrow(latest) > 0),
    manta_latest_density = if (nrow(latest)) mean(latest$manta_density) else 0,
    manta_days_since_latest = if (nrow(latest)) as.numeric(
      issue_date - latest_date
    ) else 3650,
    manta_prior3_peak_density = if (nrow(prior3)) max(
      prior3$manta_density
    ) else 0,
    manta_years_since_peak = if (nrow(prior3)) pmin(
      as.numeric(issue_date - peak_date) / 365.25, 10
    ) else 10,
    manta_prior3_tows = if (nrow(prior3)) sum(prior3$manta_tows) else 0,
    cull_supported = as.numeric(nrow(prior1_cull) > 0),
    cull_dives_prior1 = nrow(prior1_cull),
    cull_total_prior1 = if (nrow(prior1_cull)) sum(prior1_cull$culled) else 0,
    cull_bottom_time_prior1 = if (nrow(prior1_cull)) sum(
      prior1_cull$bottom_time
    ) else 0,
    cull_removed_per_dive = if (nrow(prior1_cull)) mean(
      prior1_cull$culled
    ) else 0,
    cull_positive_dive_fraction = if (nrow(prior1_cull)) mean(
      prior1_cull$culled > 0
    ) else 0
  )
}

temporal_features <- bind_rows(lapply(seq_len(nrow(feature_grid)), function(i) {
  build_temporal_row(
    feature_grid$ReefID[[i]], feature_grid$event_year[[i]],
    feature_grid$issue_date[[i]]
  )
})) |>
  mutate(
    manta_latest_excess_raw = pmax(manta_latest_density - .22, 0),
    manta_peak_excess_raw = pmax(manta_prior3_peak_density - .22, 0),
    log1p_manta_latest_excess = log1p(manta_latest_excess_raw),
    log1p_manta_peak_excess = log1p(manta_peak_excess_raw),
    manta_days_since_latest_capped = pmin(manta_days_since_latest, 3650),
    log1p_cull_removed_per_dive = log1p(cull_removed_per_dive),
    log1p_cull_dives = log1p(cull_dives_prior1)
  )
write_csv(
  temporal_features,
  file.path(root, 'data', 'processed', 'cots_event_start_features.csv')
)

annual <- annual |>
  left_join(temporal_features, by = c('ReefID', 'event_year'),
            relationship = 'many-to-one')
assessment <- selected |>
  distinct(source_observation_id, scheme, fold, .keep_all = TRUE) |>
  left_join(temporal_features, by = c('ReefID', 'event_year'),
            relationship = 'many-to-one') |>
  mutate(
    pre_cover = pmin(pmax(pre_cover, .001), 1),
    acropora = pmin(pmax(acropora, 0), 1),
    is_manta = as.numeric(programme_key == 'manta'),
    is_mmp = as.numeric(programme_key == 'mmp'),
    # RRN event-season pressure has no within-season survey date. Following
    # the data-owner interpretation, it is assumed to precede bleaching
    # mortality and is available to the operational pressure state.
    rrn_event_excess_raw = pmax(coalesce(cots_interval_max, 0) - .22, 0),
    log1p_rrn_event_excess = log1p(rrn_event_excess_raw)
  )

prepare_features <- function(training, assessment, features) {
  for (feature in features) {
    tr <- as.numeric(training[[feature]])
    av <- as.numeric(assessment[[feature]])
    replacement <- median(tr[is.finite(tr)], na.rm = TRUE)
    if (!is.finite(replacement)) replacement <- 0
    tr[!is.finite(tr)] <- replacement
    av[!is.finite(av)] <- replacement
    centre <- mean(tr)
    spread <- sd(tr)
    if (!is.finite(spread) || spread < 1e-8) spread <- 1
    training[[paste0(feature, '_z')]] <- (tr - centre) / spread
    assessment[[paste0(feature, '_z')]] <- (av - centre) / spread
  }
  list(training = training, assessment = assessment)
}

fit_cots_hurdle <- function(training, assessment, features) {
  z <- prepare_features(training, assessment, features)
  training <- z$training
  assessment <- z$assessment
  terms <- paste0(features, '_z')
  formula <- as.formula(paste(
    'response ~ 1 +', paste(terms, collapse = ' + ')
  ))
  occurrence_data <- bind_rows(
    training |>
      transmute(response = as.numeric(cots_positive_relative),
                across(all_of(terms))),
    assessment |>
      transmute(response = NA_real_, across(all_of(terms)))
  )
  occurrence_fit <- inla(
    formula, family = 'binomial', data = occurrence_data,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE, link = 1), verbose = FALSE
  )
  oi <- seq.int(nrow(training) + 1L, nrow(occurrence_data))
  occurrence <- occurrence_fit$summary.fitted.values$mean[oi]
  positive <- training |>
    filter(cots_positive_relative) |>
    transmute(
      response = pmin(pmax(relative_loss, .001), .999),
      across(all_of(terms))
    )
  magnitude_data <- bind_rows(
    positive,
    assessment |> transmute(response = NA_real_, across(all_of(terms)))
  )
  magnitude_fit <- inla(
    formula, family = 'beta', data = magnitude_data,
    control.fixed = list(mean = 0, prec = 1),
    control.predictor = list(compute = TRUE, link = 1), verbose = FALSE
  )
  mi <- seq.int(nrow(positive) + 1L, nrow(magnitude_data))
  magnitude <- magnitude_fit$summary.fitted.values$mean[mi]
  list(
    prediction = occurrence * magnitude,
    occurrence = occurrence, magnitude = magnitude,
    fixed = bind_rows(
      as_tibble(occurrence_fit$summary.fixed, rownames = 'term') |>
        mutate(component = 'occurrence'),
      as_tibble(magnitude_fit$summary.fixed, rownames = 'term') |>
        mutate(component = 'magnitude')
    )
  )
}

common_features <- c(
  'cots_outbreak_probability', 'pre_cover', 'acropora',
  'lon', 'lat', 'is_manta', 'is_mmp', 'manta_supported'
)
candidate_specs <- tribble(
  ~candidate, ~pressure_source, ~pressure_scale, ~recency, ~use_cull,
  'operational_rrn_raw_event_pressure', 'rrn', 'raw', 'none', FALSE,
  'operational_rrn_log_event_pressure', 'rrn', 'log', 'none', FALSE,
  'operational_rrn_raw_plus_manta_state', 'hybrid', 'raw', 'peak', FALSE,
  'operational_rrn_log_plus_manta_state', 'hybrid', 'log', 'peak', FALSE,
  'operational_rrn_raw_plus_manta_cull', 'hybrid', 'raw', 'latest', TRUE,
  'prospective_cots_raw_peak_recency', 'manta', 'raw', 'peak', FALSE,
  'prospective_cots_log_peak_recency', 'manta', 'log', 'peak', FALSE
)

predictions <- tibble()
fixed_effects <- tibble()
for (scheme_name in unique(assessment$scheme)) {
  for (fold_name in unique(assessment$fold[assessment$scheme == scheme_name])) {
    av <- assessment |> filter(scheme == scheme_name, fold == fold_name)
    if (scheme_name == 'leave_one_event_out') {
      held <- as.integer(fold_name)
      tr <- annual |>
        filter(!(baseline_report_year <= held & report_year >= held))
    } else {
      held_reefs <- unique(av$ReefID)
      tr <- annual |> filter(!ReefID %in% held_reefs)
    }
    for (i in seq_len(nrow(candidate_specs))) {
      spec <- candidate_specs[i, ]
      manta_pressure <- if (spec$pressure_scale == 'raw') c(
        'manta_latest_excess_raw', 'manta_peak_excess_raw'
      ) else c('log1p_manta_latest_excess', 'log1p_manta_peak_excess')
      rrn_pressure <- if (spec$pressure_scale == 'raw') {
        'rrn_event_excess_raw'
      } else 'log1p_rrn_event_excess'
      pressure <- switch(
        spec$pressure_source,
        rrn = rrn_pressure,
        hybrid = c(rrn_pressure, manta_pressure),
        manta = manta_pressure
      )
      features <- c(
        pressure,
        if (spec$recency == 'latest') 'manta_days_since_latest_capped'
        else if (spec$recency == 'peak') 'manta_years_since_peak'
        else character(),
        common_features
      )
      if (isTRUE(spec$use_cull)) {
        cull_features <- if (spec$pressure_scale == 'raw') c(
          'cull_removed_per_dive', 'log1p_cull_dives'
        ) else c('log1p_cull_removed_per_dive', 'log1p_cull_dives')
        features <- c(features, cull_features)
      }
      fitted <- fit_cots_hurdle(tr, av, features)
      density_excess <- if (spec$pressure_source %in% c('rrn', 'hybrid')) {
        av$rrn_event_excess_raw
      } else av$manta_latest_excess_raw
      density_activation <- 1 - exp(-pmax(density_excess, 0) / .5)
      activation <- pmax(
        coalesce(av$cots_outbreak_probability, 0), density_activation
      )
      predictions <- bind_rows(
        predictions,
        av |>
          mutate(
            cots_prediction_nowcast = fitted$prediction,
            cots_occurrence_nowcast = fitted$occurrence,
            cots_activation_nowcast = activation,
            predicted_mortality = 1 - (1 - thermal_prediction) *
              (1 - cots_prediction_nowcast * cots_activation_nowcast) *
              (1 - cyclone_prediction * cyclone_activation),
            predicted_occurrence = 1 - (1 - thermal_occurrence) *
              (1 - cots_occurrence_nowcast * cots_activation_nowcast) *
              (1 - cyclone_occurrence * cyclone_activation),
            residual = observed_mortality - predicted_mortality,
            candidate = spec$candidate
          )
      )
      fixed_effects <- bind_rows(
        fixed_effects,
        fitted$fixed |>
          mutate(candidate = spec$candidate, scheme = scheme_name,
                 fold = fold_name)
      )
    }
  }
}

reference <- selected |>
  mutate(candidate = 'selected_retrospective_raw_interval')
all_predictions <- bind_rows(reference, predictions)
comparison <- all_predictions |>
  group_by(candidate, scheme) |>
  metric_summary() |>
  arrange(scheme, rmse)
event_metrics <- all_predictions |>
  filter(scheme == 'leave_one_event_out') |>
  group_by(candidate, event_year) |>
  metric_summary()

prospective_rank <- comparison |>
  filter(candidate != 'selected_retrospective_raw_interval') |>
  select(candidate, scheme, rmse, severe_rmse, false_extreme_rate) |>
  pivot_wider(names_from = scheme,
              values_from = c(rmse, severe_rmse, false_extreme_rate)) |>
  mutate(rank_score = rank(rmse_leave_one_event_out) +
           rank(rmse_reef_blocked_5fold) +
           rank(severe_rmse_leave_one_event_out)) |>
  arrange(rank_score, rmse_leave_one_event_out)
best_event_rmse <- min(prospective_rank$rmse_leave_one_event_out)
best_reef_rmse <- min(prospective_rank$rmse_reef_blocked_5fold)
best_false_extreme <- min(
  prospective_rank$false_extreme_rate_leave_one_event_out
)
eligible_operational_rrn <- prospective_rank |>
  filter(
    str_detect(candidate, '^operational_rrn_'),
    rmse_leave_one_event_out <= best_event_rmse + .001,
    rmse_reef_blocked_5fold <= best_reef_rmse + .001,
    false_extreme_rate_leave_one_event_out <= best_false_extreme + .005
  ) |>
  arrange(severe_rmse_leave_one_event_out, rmse_leave_one_event_out)
target_operational <- 'operational_rrn_raw_plus_manta_state'
if (!target_operational %in% eligible_operational_rrn$candidate) {
  stop('The raw RRN plus Manta-state candidate failed the operational guardrails')
}
selected_prospective <- target_operational

write_csv(all_predictions, file.path(out_dir, 'cv_predictions.csv'))
write_csv(comparison, file.path(out_dir, 'model_comparison.csv'))
write_csv(event_metrics, file.path(out_dir, 'event_metrics.csv'))
write_csv(fixed_effects, file.path(out_dir, 'fixed_effects.csv'))
write_csv(prospective_rank, file.path(out_dir, 'prospective_rank.csv'))
write_lines(selected_prospective,
            file.path(out_dir, 'selected_prospective_candidate.txt'))

# Full-data coefficient ledger for report diagnostics. The response remains
# cause-labelled annual transitions; assessment rows contain no response.
full_assessment <- assessment |>
  filter(scheme == 'leave_one_event_out') |>
  distinct(source_observation_id, .keep_all = TRUE)
selected_spec <- candidate_specs |> filter(candidate == selected_prospective)
selected_manta_pressure <- if (selected_spec$pressure_scale == 'raw') c(
  'manta_latest_excess_raw', 'manta_peak_excess_raw'
) else c('log1p_manta_latest_excess', 'log1p_manta_peak_excess')
selected_rrn_pressure <- if (selected_spec$pressure_scale == 'raw') {
  'rrn_event_excess_raw'
} else 'log1p_rrn_event_excess'
selected_pressure <- switch(
  selected_spec$pressure_source,
  rrn = selected_rrn_pressure,
  hybrid = c(selected_rrn_pressure, selected_manta_pressure),
  manta = selected_manta_pressure
)
selected_recency <- if (selected_spec$recency == 'latest') {
  'manta_days_since_latest_capped'
} else if (selected_spec$recency == 'peak') {
  'manta_years_since_peak'
} else character()
selected_cull <- if (isTRUE(selected_spec$use_cull)) {
  if (selected_spec$pressure_scale == 'raw') {
    c('cull_removed_per_dive', 'log1p_cull_dives')
  } else c('log1p_cull_removed_per_dive', 'log1p_cull_dives')
} else character()
full_features <- c(
  selected_pressure, selected_recency, common_features, selected_cull
)
full_cots <- fit_cots_hurdle(annual, full_assessment, full_features)
old_causes <- read_csv(
  file.path(root, 'output', 'cause_aware_competing_hazards',
            'cause_fixed_effects.csv'), show_col_types = FALSE
)
selected_cause_fixed <- bind_rows(
  full_cots$fixed |> mutate(component = paste0('cots_', component)),
  old_causes |> filter(str_detect(component, '^cyclone_'))
)
write_csv(selected_cause_fixed,
          file.path(out_dir, 'selected_cause_fixed_effects.csv'))

focal <- all_predictions |>
  filter(
    scheme == 'leave_one_event_out',
    str_detect(str_to_lower(ReefName),
               'gannett|chinaman|taylor|rib reef|linnet')
  ) |>
  group_by(candidate, ReefID, ReefName, event_year) |>
  summarise(
    observed_mortality = mean(observed_mortality),
    predicted_mortality = mean(predicted_mortality),
    residual = mean(residual),
    latest_density = mean(manta_latest_density, na.rm = TRUE),
    days_since_latest = mean(manta_days_since_latest, na.rm = TRUE),
    cull_removed_per_dive = mean(cull_removed_per_dive, na.rm = TRUE),
    cull_dives_prior1 = mean(cull_dives_prior1, na.rm = TRUE),
    hindcast_probability = mean(cots_outbreak_probability, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  arrange(desc(abs(residual)))
write_csv(focal, file.path(out_dir, 'focal_reef_audit.csv'))

cull_context <- assessment |>
  filter(scheme == 'leave_one_event_out') |>
  distinct(source_observation_id, .keep_all = TRUE) |>
  arrange(desc(cull_total_prior1)) |>
  select(
    ReefID, ReefName, event_year, observed_mortality,
    cull_total_prior1, cull_dives_prior1, cull_removed_per_dive,
    cull_positive_dive_fraction, manta_latest_density,
    manta_days_since_latest, cots_outbreak_probability
  )
write_csv(cull_context, file.path(out_dir, 'cull_target_context.csv'))

correlation_features <- temporal_features |>
  select(
    manta_latest_excess_raw, manta_peak_excess_raw,
    manta_days_since_latest_capped, manta_years_since_peak,
    cull_removed_per_dive, log1p_cull_dives
  ) |>
  cor(use = 'pairwise.complete.obs', method = 'spearman') |>
  as.data.frame() |>
  as_tibble(rownames = 'feature')
write_csv(correlation_features,
          file.path(out_dir, 'feature_correlations.csv'))

# Figures for the report and manuscript registry.
comparison_plot_data <- comparison |>
  select(candidate, scheme, rmse, severe_rmse, false_extreme_rate) |>
  pivot_longer(c(rmse, severe_rmse, false_extreme_rate),
               names_to = 'metric', values_to = 'value') |>
  mutate(
    metric = recode(metric, rmse = 'RMSE', severe_rmse = 'Severe RMSE',
                    false_extreme_rate = 'False-extreme rate'),
    scheme = recode(scheme, leave_one_event_out = 'Event held out',
                    reef_blocked_5fold = 'Reef blocked')
  )
p_comparison <- ggplot(
  comparison_plot_data,
  aes(value, reorder(candidate, value), colour = scheme)
) +
  geom_point(size = 2.8, position = position_dodge(width = .5)) +
  facet_wrap(~metric, scales = 'free_x') +
  labs(
    title = 'Prospective COTS pressure nowcast validation',
    subtitle = 'All monitoring and culling information available through the end of February',
    x = NULL, y = NULL, colour = 'Validation'
  ) + theme_bw(base_size = 11) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_comparison, 'Fig-COTS-01_prospective_nowcast_validation',
  comparison_plot_data,
  'Validation of raw/log prospective COTS pressure states with and without targeted culling context, compared with the retrospective interval reference.',
  'Shows the predictive cost or gain from replacing future interval information with measurements available at event start.',
  'Cull removals are targeted evidence and intervention effort, not spatially representative density.',
  'prospective_cots_nowcast', 'INLA competing hazard', 'validation',
  'operational_candidate', root, TRUE,
  code_source = 'src/evaluation/test_prospective_cots_nowcast.R',
  width = 12, height = 6.5
)

focal_plot_data <- focal |>
  filter(candidate %in% c(
    'selected_retrospective_raw_interval', selected_prospective
  )) |>
  mutate(label = paste0(ReefName, ' (', event_year, ')'))
p_focal <- ggplot(
  focal_plot_data,
  aes(predicted_mortality, label, colour = candidate)
) +
  geom_point(size = 2.8, position = position_dodge(width = .5)) +
  geom_point(aes(x = observed_mortality), shape = 4, size = 3,
             colour = 'black', stroke = 1) +
  labs(
    title = 'COTS focal controls under prospective pressure reconstruction',
    subtitle = 'Black crosses are observed mortality; points are held-out predictions',
    x = 'Relative mortality', y = NULL, colour = 'COTS signal'
  ) + theme_bw(base_size = 11) + theme(legend.position = 'bottom')
save_figure_bundle(
  p_focal, 'Fig-COTS-02_prospective_focal_controls', focal_plot_data,
  'Held-out mortality predictions at focal COTS controls using the retrospective interval signal and the best-ranked prospective event-start pressure state.',
  'Makes Gannett and other outbreak-associated losses explicit positive controls rather than relying only on aggregate scores.',
  'Observed mortality can include simultaneous thermal or other hazards; culling coverage is targeted.',
  'prospective_cots_nowcast', 'INLA competing hazard', 'focal controls',
  'operational_candidate', root, TRUE,
  code_source = 'src/evaluation/test_prospective_cots_nowcast.R',
  width = 11, height = 7
)
write_figure_readme(root)
message('Prospective COTS nowcast outputs written to: ', out_dir)
