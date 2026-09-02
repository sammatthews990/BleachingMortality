# Spatially independent, leakage-safe within-event bleaching update.
# This script does not alter the initial operational forecast. It tests whether
# aerial or rapid in-water bleaching observations available by March/April can
# update the occurrence component at reefs in held-out sectors/spatial blocks.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(INLA)
  library(readr)
  library(readxl)
  library(sf)
  library(stringr)
  library(tidyr)
})
Sys.setenv(INLA_ST_RUN = '0')
source('scripts/model_registry.R')
source('scripts/model_diagnostics.R')

root <- normalizePath('.', winslash = '/', mustWork = TRUE)
out_dir <- file.path(root, 'output', 'spatial_early_bleaching_update')
fig_dir <- file.path(root, 'output', 'fig')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

event_years <- c(2016L, 2017L, 2020L, 2022L, 2024L)
max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

metric_summary <- function(x) {
  severe <- x$observed_mortality >= .5
  tibble(
    n = nrow(x),
    n_severe = sum(severe),
    rmse = sqrt(mean((x$observed_mortality - x$predicted_mortality)^2)),
    mae = mean(abs(x$observed_mortality - x$predicted_mortality)),
    predictive_r2 = 1 - sum((x$observed_mortality - x$predicted_mortality)^2) /
      sum((x$observed_mortality - mean(x$observed_mortality))^2),
    severe_rmse = if (any(severe)) sqrt(mean(
      (x$observed_mortality[severe] - x$predicted_mortality[severe])^2
    )) else NA_real_,
    occurrence_brier = mean(
      (x$observed_occurrence - x$predicted_occurrence)^2
    ),
    false_extreme_rate = mean(
      x$predicted_mortality[x$observed_mortality < .3] >= .3
    )
  )
}

normalise_reef_id <- function(x) str_to_upper(str_trim(as.character(x)))

# ---------------------------------------------------------------------------
# Reef geography and the selected initial forecast
# ---------------------------------------------------------------------------
reef_reference <- read_csv(
  file.path(root, 'data', 'AIMS-Reef_Reference.csv'), show_col_types = FALSE
) |>
  transmute(
    ReefID = normalise_reef_id(ReefID),
    reef_name_reference = ReefName,
    SECTOR = SECTOR,
    sector_name = SECT_NAME
  ) |>
  distinct(ReefID, .keep_all = TRUE)

reef_features <- st_read(
  file.path(root, 'data', 'Great_Barrier_Reef_Features',
            'Great_Barrier_Reef_Features.shp'), quiet = TRUE
) |>
  st_drop_geometry() |>
  transmute(
    ReefID = normalise_reef_id(LABEL_ID),
    lon = as.numeric(X_COORD), lat = as.numeric(Y_COORD)
  ) |>
  filter(between(lon, 142, 155), between(lat, -26, -9)) |>
  group_by(ReefID) |>
  summarise(lon = mean(lon), lat = mean(lat), .groups = 'drop') |>
  left_join(reef_reference, by = 'ReefID')

block_breaks <- quantile(
  reef_features$lat, probs = seq(0, 1, length.out = 6), na.rm = TRUE
)
block_breaks[1] <- block_breaks[1] - 1e-6
block_breaks[length(block_breaks)] <- block_breaks[length(block_breaks)] + 1e-6
assign_block <- function(lat) as.character(cut(
  lat, breaks = block_breaks, include.lowest = TRUE,
  labels = paste0('lat_block_', 1:5)
))
reef_features <- reef_features |>
  mutate(spatial_block = assign_block(lat))

joint_context <- read_csv(
  file.path(root, 'output', 'explanatory_event_dhw',
            'event_dhw_brt_data.csv'), show_col_types = FALSE
) |>
  mutate(source_observation_id = as.character(source_observation_id)) |>
  distinct(source_observation_id, .keep_all = TRUE) |>
  select(source_observation_id, SECTOR, survey_date)

base <- read_csv(
  file.path(root, 'output', 'cots_raw_enso_occurrence',
            'cots_cv_predictions.csv'), show_col_types = FALSE
) |>
  filter(
    candidate == 'cots_raw_interval_relative',
    scheme == 'leave_one_event_out', event_year %in% event_years
  ) |>
  mutate(
    source_observation_id = as.character(source_observation_id),
    ReefID = normalise_reef_id(ReefID),
    predicted_occurrence = pmin(pmax(predicted_occurrence, .001), .999),
    base_conditional_magnitude = pmin(
      predicted_mortality / predicted_occurrence, .999
    )
  ) |>
  left_join(joint_context, by = 'source_observation_id') |>
  mutate(
    SECTOR = coalesce(SECTOR, 'UNKNOWN'),
    spatial_block = assign_block(lat),
    target_key = paste(ReefID, event_year, sep = '__')
  )

targets <- base |>
  distinct(target_key, ReefID, event_year, lon, lat, SECTOR, spatial_block)

# ---------------------------------------------------------------------------
# RHIS conversion: requested 0--4 severity scale, community weighting and an
# explicitly separate recently-dead indicator.
# ---------------------------------------------------------------------------
rhis_raw <- read_excel(
  file.path(root, 'data',
            '250929_COTS-Manta-Cull-RHIS-Data-Matthews-and-Schlawinsky.xlsx'),
  sheet = 'RHIS'
)

morphologies <- c(
  'Soft', 'Branching', 'Bushy', 'Plate/Table',
  'Vase/Foliose', 'Encrusting', 'Mushroom', 'Massive'
)
cover_cols <- paste0(morphologies, ' Coral % (Live and Recently Dead)')
severity_cols <- paste0(morphologies, ' Bleaching Severity')
bleached_cols <- paste0(morphologies, ' Bleached %')

severity_map <- c(
  'None' = 0,
  'Bleached only on upper surface' = 1,
  'Pale/fluoro (very light or yellowish)' = 2,
  'Totally bleached white' = 3,
  'Recently dead coral lightly covered in algae' = 4
)
severity_matrix <- sapply(severity_cols, function(nm) {
  value <- as.character(rhis_raw[[nm]])
  unname(severity_map[value])
})
cover_matrix <- sapply(cover_cols, function(nm) as.numeric(rhis_raw[[nm]]))
bleached_matrix <- sapply(bleached_cols, function(nm) {
  pmin(pmax(as.numeric(rhis_raw[[nm]]) / 100, 0), 1)
})

valid_severity <- is.finite(severity_matrix) & is.finite(cover_matrix) &
  cover_matrix > 0
severity_numerator <- rowSums(ifelse(
  valid_severity, cover_matrix * severity_matrix / 4, 0
), na.rm = TRUE)
severity_denominator <- rowSums(ifelse(valid_severity, cover_matrix, 0),
                                na.rm = TRUE)
severity_fallback <- rowMeans(severity_matrix / 4, na.rm = TRUE)
severity_fallback[!is.finite(severity_fallback)] <- NA_real_
community_severity <- ifelse(
  severity_denominator > 0,
  severity_numerator / severity_denominator,
  severity_fallback
)

valid_burden <- valid_severity & is.finite(bleached_matrix)
burden_numerator <- rowSums(ifelse(
  valid_burden,
  cover_matrix * (severity_matrix / 4) * bleached_matrix, 0
), na.rm = TRUE)
burden_denominator <- rowSums(ifelse(valid_burden, cover_matrix, 0),
                              na.rm = TRUE)
community_burden <- ifelse(
  burden_denominator > 0, burden_numerator / burden_denominator,
  community_severity
)

bleaching_present <- rhis_raw[['Bleaching Present']] %in% TRUE
community_severity[!bleaching_present] <- 0
community_burden[!bleaching_present] <- 0
recent_dead <- pmin(pmax(
  as.numeric(rhis_raw[['%Benthos_Recently_Dead_Coral']]) / 100, 0
), 1)
recent_dead[!bleaching_present | !is.finite(recent_dead)] <- 0

rhis_observations <- rhis_raw |>
  transmute(
    SurveyId = as.character(SurveyId),
    survey_date = as.Date(SurveyTime),
    ReefID = normalise_reef_id(str_match(
      ReefName, '([0-9]{2}-[0-9]{3}[A-Za-z]?)'
    )[, 2]),
    ReefName = ReefName,
    lon = as.numeric(Longitude), lat = as.numeric(Latitude),
    bleaching_present = as.numeric(bleaching_present),
    rhis_severity = pmin(pmax(community_severity, 0), 1),
    rhis_burden = pmin(pmax(community_burden, 0), 1),
    rhis_recent_dead = recent_dead
  ) |>
  mutate(
    event_year = as.integer(format(survey_date, '%Y')),
    cutoff_april = survey_date <= as.Date(paste0(event_year, '-04-30')),
    cutoff_march = survey_date <= as.Date(paste0(event_year, '-03-31'))
  ) |>
  filter(event_year %in% event_years, !is.na(ReefID)) |>
  left_join(
    reef_features |>
      select(ReefID, ref_lon = lon, ref_lat = lat, SECTOR, spatial_block),
    by = 'ReefID'
  ) |>
  mutate(
    lon = coalesce(lon, ref_lon), lat = coalesce(lat, ref_lat),
    SECTOR = coalesce(SECTOR, 'UNKNOWN'),
    spatial_block = coalesce(spatial_block, assign_block(lat))
  ) |>
  select(-ref_lon, -ref_lat)

aggregate_rhis <- function(rows, cutoff_name) {
  rows |>
    filter(.data[[cutoff_name]], is.finite(lon), is.finite(lat)) |>
    group_by(ReefID, event_year) |>
    summarise(
      ReefName = first(ReefName), lon = mean(lon), lat = mean(lat),
      SECTOR = first(SECTOR), spatial_block = first(spatial_block),
      rhis_bleaching_present = max(bleaching_present),
      rhis_severity = max_or_na(rhis_severity),
      rhis_burden = max_or_na(rhis_burden),
      rhis_recent_dead = max_or_na(rhis_recent_dead),
      first_survey = min(survey_date), last_survey = max(survey_date),
      n_rhis_records = n(), .groups = 'drop'
    ) |>
    mutate(cutoff = recode(
      cutoff_name, cutoff_april = 'April 30', cutoff_march = 'March 31'
    ))
}
rhis_sources <- bind_rows(
  aggregate_rhis(rhis_observations, 'cutoff_april'),
  aggregate_rhis(rhis_observations, 'cutoff_march')
)

write_csv(
  rhis_sources,
  file.path(root, 'data', 'processed', 'rhis_early_bleaching_reef_event.csv')
)

# Hughes et al. aerial response is a binary severe-bleaching score. Survey
# dates are absent from the supplied CSV, so this source is treated as an
# event-time update with timing provenance explicitly marked as unavailable.
aerial_sources <- read_csv(
  file.path(root, 'data', 'Hughes2021', 'data', 'bleach_all.csv'),
  show_col_types = FALSE
) |>
  transmute(
    ReefID = normalise_reef_id(ReefID), event_year = as.integer(year),
    aerial_severe_bleaching = as.numeric(bin.score), aerial_dhw = DHW
  ) |>
  filter(event_year %in% c(2016L, 2017L, 2020L)) |>
  left_join(
    reef_features |>
      select(ReefID, lon, lat, SECTOR, spatial_block, reef_name_reference),
    by = 'ReefID'
  ) |>
  filter(is.finite(lon), is.finite(lat)) |>
  distinct(ReefID, event_year, .keep_all = TRUE) |>
  mutate(timing_provenance = 'Event-time aerial survey; date absent in source file')

write_csv(
  aerial_sources,
  file.path(root, 'data', 'processed', 'aerial_early_bleaching_reef_event.csv')
)

# ---------------------------------------------------------------------------
# Kernel summaries. The target sector/block is excluded before calculating
# each signal; no target mortality response enters these summaries.
# ---------------------------------------------------------------------------
haversine_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  dlon <- (lon2 - lon1) * rad
  dlat <- (lat2 - lat1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(1 - a, 0)))
}

kernel_signal <- function(target_rows, source_rows, value, exclusion,
                          bandwidth_km = 150, maximum_km = 450) {
  output <- vector('list', nrow(target_rows))
  for (i in seq_len(nrow(target_rows))) {
    target <- target_rows[i, ]
    sources <- source_rows |>
      filter(event_year == target$event_year, is.finite(.data[[value]]))
    if (exclusion == 'sector') {
      sources <- sources |>
        filter(SECTOR != target$SECTOR, SECTOR != 'UNKNOWN')
    } else {
      sources <- sources |>
        filter(spatial_block != target$spatial_block,
               !is.na(spatial_block))
    }
    if (!nrow(sources)) {
      output[[i]] <- tibble(signal = NA_real_, n_source = 0L,
                            effective_n = 0, nearest_km = NA_real_)
      next
    }
    distance <- haversine_km(
      target$lon, target$lat, sources$lon, sources$lat
    )
    keep <- is.finite(distance) & distance <= maximum_km
    if (!any(keep)) {
      output[[i]] <- tibble(signal = NA_real_, n_source = 0L,
                            effective_n = 0, nearest_km = NA_real_)
      next
    }
    distance <- distance[keep]
    values <- sources[[value]][keep]
    weight <- exp(-.5 * (distance / bandwidth_km)^2)
    effective_n <- sum(weight)^2 / sum(weight^2)
    output[[i]] <- tibble(
      signal = if (effective_n >= 3) weighted.mean(values, weight) else NA_real_,
      n_source = length(values), effective_n = effective_n,
      nearest_km = min(distance)
    )
  }
  bind_cols(target_rows |> select(target_key), bind_rows(output))
}

signal_specs <- tribble(
  ~source, ~cutoff, ~field, ~prefix,
  'RHIS', 'April 30', 'rhis_severity', 'rhis_severity_apr',
  'RHIS', 'April 30', 'rhis_burden', 'rhis_burden_apr',
  'RHIS', 'April 30', 'rhis_recent_dead', 'rhis_recent_dead_apr',
  'RHIS', 'March 31', 'rhis_severity', 'rhis_severity_mar',
  'RHIS', 'March 31', 'rhis_burden', 'rhis_burden_mar',
  'RHIS', 'March 31', 'rhis_recent_dead', 'rhis_recent_dead_mar',
  'Aerial', NA_character_, 'aerial_severe_bleaching', 'aerial_severe'
)

feature_table <- targets
coverage <- tibble()
for (design in c('sector', 'block')) {
  for (i in seq_len(nrow(signal_specs))) {
    spec <- signal_specs[i, ]
    sources <- if (spec$source == 'RHIS') {
      rhis_sources |> filter(cutoff == spec$cutoff)
    } else aerial_sources
    result <- kernel_signal(targets, sources, spec$field, design)
    stem <- paste(spec$prefix, design, sep = '_')
    names(result)[2:5] <- paste0(
      stem, c('', '_n_source', '_effective_n', '_nearest_km')
    )
    feature_table <- feature_table |>
      left_join(result, by = 'target_key', relationship = 'one-to-one')
    coverage <- bind_rows(
      coverage,
      result |>
        left_join(targets |> select(target_key, event_year), by = 'target_key') |>
        group_by(event_year) |>
        summarise(
          source = spec$source, cutoff = coalesce(spec$cutoff, 'not supplied'),
          field = spec$field, design = design,
          target_reefs = n(), supported_reefs = sum(is.finite(.data[[stem]])),
          coverage = mean(is.finite(.data[[stem]])),
          median_effective_n = median(
            .data[[paste0(stem, '_effective_n')]], na.rm = TRUE
          ),
          median_nearest_km = median(
            .data[[paste0(stem, '_nearest_km')]], na.rm = TRUE
          ), .groups = 'drop'
        )
    )
  }
}
write_csv(feature_table, file.path(out_dir, 'spatial_early_features.csv'))
write_csv(coverage, file.path(out_dir, 'spatial_source_coverage.csv'))

source_summary <- bind_rows(
  rhis_sources |>
    group_by(event_year, cutoff) |>
    summarise(
      source = 'RHIS', reefs = n(), records = sum(n_rhis_records),
      bleaching_prevalence = mean(rhis_bleaching_present),
      mean_severity = mean(rhis_severity, na.rm = TRUE),
      mean_burden = mean(rhis_burden, na.rm = TRUE),
      mean_recent_dead = mean(rhis_recent_dead, na.rm = TRUE),
      .groups = 'drop'
    ),
  aerial_sources |>
    group_by(event_year) |>
    summarise(
      source = 'Aerial', cutoff = 'date absent', reefs = n(), records = n(),
      bleaching_prevalence = mean(aerial_severe_bleaching),
      mean_severity = NA_real_, mean_burden = NA_real_,
      mean_recent_dead = NA_real_, .groups = 'drop'
    )
)
write_csv(source_summary, file.path(out_dir, 'early_source_summary.csv'))

# ---------------------------------------------------------------------------
# Leave-one-event-out occurrence calibration. Features are standardised only
# on training events. The base conditional magnitude is retained unchanged.
# ---------------------------------------------------------------------------
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

fit_occurrence_update <- function(training, assessment, features) {
  z <- prepare_features(training, assessment, features)
  training <- z$training
  assessment <- z$assessment
  terms <- paste0(features, '_z')
  formula <- as.formula(paste(
    'response ~ 1 +', paste(terms, collapse = ' + '), '+ offset(base_logit)'
  ))
  combined <- bind_rows(
    training |>
      transmute(
        response = observed_occurrence,
        base_logit = qlogis(predicted_occurrence), across(all_of(terms))
      ),
    assessment |>
      transmute(
        response = NA_real_,
        base_logit = qlogis(predicted_occurrence), across(all_of(terms))
      )
  )
  fit <- inla(
    formula, family = 'binomial', data = combined,
    control.fixed = list(
      mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
    ),
    control.predictor = list(compute = TRUE, link = 1), verbose = FALSE
  )
  index <- seq.int(nrow(training) + 1L, nrow(combined))
  list(
    occurrence = fit$summary.fitted.values$mean[index],
    fixed = as_tibble(fit$summary.fixed, rownames = 'term')
  )
}

fit_magnitude_update <- function(training, assessment, features) {
  positive_training <- training |>
    filter(observed_occurrence == 1) |>
    mutate(observed_magnitude = pmin(pmax(observed_mortality, .001), .999))
  z <- prepare_features(positive_training, assessment, features)
  positive_training <- z$training
  assessment <- z$assessment
  terms <- paste0(features, '_z')
  formula <- as.formula(paste(
    'response ~ 1 +', paste(terms, collapse = ' + '),
    '+ offset(base_magnitude_logit)'
  ))
  combined <- bind_rows(
    positive_training |>
      transmute(
        response = observed_magnitude,
        base_magnitude_logit = qlogis(pmin(pmax(
          base_conditional_magnitude, .001
        ), .999)),
        across(all_of(terms))
      ),
    assessment |>
      transmute(
        response = NA_real_,
        base_magnitude_logit = qlogis(pmin(pmax(
          base_conditional_magnitude, .001
        ), .999)),
        across(all_of(terms))
      )
  )
  fit <- inla(
    formula, family = 'beta', data = combined,
    control.fixed = list(
      mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
    ),
    control.predictor = list(compute = TRUE, link = 1), verbose = FALSE
  )
  index <- seq.int(nrow(positive_training) + 1L, nrow(combined))
  list(
    magnitude = fit$summary.fitted.values$mean[index],
    fixed = as_tibble(fit$summary.fixed, rownames = 'term')
  )
}

candidate_specs <- tribble(
  ~candidate, ~scope, ~feature_key, ~update_components,
  'RHIS severity only (April 30)', 'all_events', 'rhis_severity_apr',
  'occurrence',
  'RHIS burden + recent dead (April 30)', 'all_events',
  'rhis_burden_apr;rhis_recent_dead_apr',
  'occurrence',
  'RHIS burden + recent dead (March 31)', 'all_events',
  'rhis_burden_mar;rhis_recent_dead_mar',
  'occurrence',
  'Aerial severe bleaching', 'aerial_events', 'aerial_severe', 'occurrence',
  'Aerial + RHIS burden/recent dead', 'aerial_events',
  'aerial_severe;rhis_burden_apr;rhis_recent_dead_apr',
  'occurrence',
  'RHIS burden/recent dead: two-part update', 'all_events',
  'rhis_burden_apr;rhis_recent_dead_apr',
  'occurrence_and_magnitude',
  'Aerial + RHIS: two-part update', 'aerial_events',
  'aerial_severe;rhis_burden_apr;rhis_recent_dead_apr',
  'occurrence_and_magnitude'
)

model_rows <- base |>
  left_join(
    feature_table |> select(-ReefID, -event_year, -lon, -lat, -SECTOR,
                            -spatial_block),
    by = 'target_key', relationship = 'many-to-one'
  )
predictions <- tibble()
fixed_effects <- tibble()

for (design in c('sector', 'block')) {
  for (i in seq_len(nrow(candidate_specs))) {
    spec <- candidate_specs[i, ]
    feature_stems <- str_split_1(spec$feature_key, ';')
    features <- paste(feature_stems, design, sep = '_')
    rows <- model_rows
    if (spec$scope == 'aerial_events') {
      rows <- rows |> filter(event_year %in% c(2016L, 2017L, 2020L))
    }
    for (held_event in sort(unique(rows$event_year))) {
      training <- rows |> filter(event_year != held_event)
      assessment <- rows |> filter(event_year == held_event)
      fitted <- fit_occurrence_update(training, assessment, features)
      updated_occurrence <- pmin(pmax(fitted$occurrence, .001), .999)
      magnitude_fitted <- if (spec$update_components ==
                              'occurrence_and_magnitude') {
        fit_magnitude_update(training, assessment, features)
      } else NULL
      updated_magnitude <- if (is.null(magnitude_fitted)) {
        assessment$base_conditional_magnitude
      } else pmin(pmax(magnitude_fitted$magnitude, .001), .999)
      predictions <- bind_rows(
        predictions,
        assessment |>
          transmute(
            source_observation_id, programme_key, ReefID, ReefName,
            event_year, lon, lat, SECTOR, spatial_block,
            observed_mortality, observed_occurrence,
            initial_mortality = predicted_mortality,
            initial_occurrence = predicted_occurrence,
            predicted_occurrence = updated_occurrence,
            predicted_magnitude = updated_magnitude,
            predicted_mortality = updated_occurrence * updated_magnitude,
            candidate = spec$candidate, scope = spec$scope,
            update_components = spec$update_components,
            design = design, held_event = held_event,
            across(all_of(features))
          )
      )
      fixed_effects <- bind_rows(
        fixed_effects,
        fitted$fixed |>
          mutate(
            component = 'occurrence',
            candidate = spec$candidate, scope = spec$scope,
            design = design, held_event = held_event
          )
      )
      if (!is.null(magnitude_fitted)) {
        fixed_effects <- bind_rows(
          fixed_effects,
          magnitude_fitted$fixed |>
            mutate(
              component = 'conditional_magnitude',
              candidate = spec$candidate, scope = spec$scope,
              design = design, held_event = held_event
            )
        )
      }
    }
  }
}

write_csv(predictions, file.path(out_dir, 'cv_predictions.csv'))
write_csv(fixed_effects, file.path(out_dir, 'fixed_effects.csv'))

comparison <- predictions |>
  group_by(candidate, scope, design) |>
  group_modify(~ {
    updated <- metric_summary(.x)
    initial <- metric_summary(.x |>
      mutate(
        predicted_mortality = initial_mortality,
        predicted_occurrence = initial_occurrence
      ))
    bind_cols(
      updated |> rename_with(~ paste0('updated_', .x)),
      initial |> rename_with(~ paste0('initial_', .x))
    ) |>
      mutate(
        delta_rmse = updated_rmse - initial_rmse,
        delta_predictive_r2 = updated_predictive_r2 - initial_predictive_r2,
        delta_severe_rmse = updated_severe_rmse - initial_severe_rmse,
        delta_occurrence_brier = updated_occurrence_brier -
          initial_occurrence_brier
      )
  }) |>
  ungroup() |>
  arrange(delta_rmse)

event_metrics <- predictions |>
  group_by(candidate, scope, design, event_year) |>
  group_modify(~ {
    updated <- metric_summary(.x)
    initial <- metric_summary(.x |>
      mutate(
        predicted_mortality = initial_mortality,
        predicted_occurrence = initial_occurrence
      ))
    bind_cols(
      updated |> rename_with(~ paste0('updated_', .x)),
      initial |> rename_with(~ paste0('initial_', .x))
    ) |>
      mutate(
        delta_rmse = updated_rmse - initial_rmse,
        delta_severe_rmse = updated_severe_rmse - initial_severe_rmse,
        delta_occurrence_brier = updated_occurrence_brier -
          initial_occurrence_brier
      )
  }) |>
  ungroup()

write_csv(comparison, file.path(out_dir, 'model_comparison.csv'))
write_csv(event_metrics, file.path(out_dir, 'event_metrics.csv'))

signal_correlation <- feature_table |>
  select(target_key, event_year, matches('^(rhis|aerial).+_(sector|block)$')) |>
  pivot_longer(-c(target_key, event_year), names_to = 'feature',
               values_to = 'value') |>
  filter(!str_detect(feature, '_n_source|_effective_n|_nearest_km')) |>
  separate(feature, into = c('stem', 'design'), sep = '_(?=[^_]+$)') |>
  pivot_wider(names_from = stem, values_from = value) |>
  group_by(design) |>
  summarise(
    rhis_severity_burden = cor(
      rhis_severity_apr, rhis_burden_apr,
      use = 'pairwise.complete.obs', method = 'spearman'
    ),
    burden_recent_dead = cor(
      rhis_burden_apr, rhis_recent_dead_apr,
      use = 'pairwise.complete.obs', method = 'spearman'
    ), .groups = 'drop'
  )
write_csv(signal_correlation, file.path(out_dir, 'signal_correlation.csv'))

# ---------------------------------------------------------------------------
# Diagnostic and paper-figure bundles.
# ---------------------------------------------------------------------------
p_sources <- source_summary |>
  filter(source == 'RHIS', cutoff == 'April 30') |>
  ggplot(aes(factor(event_year), bleaching_prevalence)) +
  geom_col(fill = '#2C7FB8', width = .72) +
  geom_text(aes(label = paste0(reefs, ' reefs')), vjust = -0.35, size = 3.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, .12))) +
  labs(
    title = 'Early RHIS bleaching observations available by April 30',
    subtitle = 'Bleaching prevalence; mortality surveys after the event are not used',
    x = 'Bleaching event', y = 'RHIS reefs with bleaching present'
  ) +
  theme_bw(base_size = 12)

save_figure_bundle(
  p_sources, 'Fig-NOWCAST-01_early_source_coverage', source_summary,
  caption = paste(
    'Availability and prevalence of RHIS bleaching observations by April 30',
    'for each event; aerial data are summarised separately because dates are',
    'not present in the supplied file.'
  ),
  interpretation = paste(
    'Rapid in-water observations exist for all five events and can support a',
    'within-event occurrence update.'
  ),
  caveats = paste(
    'April 30 is a common operational cutoff rather than a reef-specific heat',
    'stress peak; March 31 is tested as a stricter sensitivity.'
  ),
  model_id = 'spatial_early_bleaching_update', framework = 'INLA occurrence update',
  figure_type = 'source coverage', analysis_role = 'operational update validation',
  root = root, manuscript_candidate = FALSE,
  code_source = 'scripts/test_spatial_early_bleaching_update.R',
  width = 9, height = 5.5
)

comparison_plot_data <- comparison |>
  select(candidate, design, delta_rmse, delta_severe_rmse) |>
  pivot_longer(starts_with('delta_'), names_to = 'metric', values_to = 'delta') |>
  mutate(
    metric = recode(metric, delta_rmse = 'Overall RMSE',
                    delta_severe_rmse = 'Severe RMSE'),
    design = recode(design, sector = 'Held-out sector',
                    block = 'Held-out latitude block')
  )
p_comparison <- ggplot(
  comparison_plot_data,
  aes(delta, reorder(candidate, delta), colour = design)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = 'grey45') +
  geom_point(size = 2.8, position = position_dodge(width = .55)) +
  facet_wrap(~metric, scales = 'free_x') +
  labs(
    title = 'Spatially independent early-event update validation',
    subtitle = 'Change from the same initial forecast; negative RMSE change is better',
    x = 'Updated minus initial RMSE', y = NULL, colour = 'Source exclusion'
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom')

save_figure_bundle(
  p_comparison, 'Fig-NOWCAST-02_spatial_update_validation', comparison_plot_data,
  caption = paste(
    'Change in overall and severe-event RMSE after applying early aerial/RHIS',
    'bleaching indicators, with every target sector or latitude block excluded',
    'from construction of its update signal.'
  ),
  interpretation = paste(
    'Negative values indicate that independent early-event observations improve',
    'the selected initial operational forecast.'
  ),
  caveats = paste(
    'Aerial comparisons cover 2016, 2017 and 2020 only; 2022 and 2024 aerial',
    'scores have not yet been supplied.'
  ),
  model_id = 'spatial_early_bleaching_update', framework = 'INLA occurrence update',
  figure_type = 'held-out validation', analysis_role = 'operational update validation',
  root = root, manuscript_candidate = TRUE,
  code_source = 'scripts/test_spatial_early_bleaching_update.R',
  width = 12, height = 6.5
)

event_plot_data <- event_metrics |>
  mutate(
    design = recode(design, sector = 'Held-out sector',
                    block = 'Held-out latitude block'),
    event_year = factor(event_year)
  )
p_events <- ggplot(
  event_plot_data,
  aes(event_year, delta_rmse, colour = candidate, group = candidate)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey45') +
  geom_line(linewidth = .65) + geom_point(size = 2) +
  facet_wrap(~design) +
  labs(
    title = 'Does the early-event update transfer among bleaching events?',
    subtitle = 'Event-held-out RMSE change; negative values improve the initial forecast',
    x = 'Held-out event', y = 'Updated minus initial RMSE', colour = 'Update signal'
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom')

save_figure_bundle(
  p_events, 'Fig-NOWCAST-03_event_transfer', event_plot_data,
  caption = paste(
    'Event-specific RMSE change for spatially independent early-event occurrence',
    'updates. Each event is predicted from update coefficients fitted to other events.'
  ),
  interpretation = paste(
    'The figure identifies whether gains are general or dominated by a single',
    'event such as 2020.'
  ),
  caveats = paste(
    'With five events, event-level heterogeneity remains a major source of',
    'uncertainty; aerial candidates have only three event folds.'
  ),
  model_id = 'spatial_early_bleaching_update', framework = 'INLA occurrence update',
  figure_type = 'event transfer', analysis_role = 'operational update validation',
  root = root, manuscript_candidate = TRUE,
  code_source = 'scripts/test_spatial_early_bleaching_update.R',
  width = 12, height = 7
)

write_figure_readme(root)
message('Spatial early-bleaching update outputs written to: ', out_dir)
