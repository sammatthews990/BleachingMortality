# Reconstruct and lock a 2025 initial forecast before testing aerial/RHIS
# updates. The fit sees only pre-2025 outcomes; 2025 responses are attached for
# scoring after prediction and never enter an INLA response vector.

suppressPackageStartupMessages({
    library(dplyr)
    library(INLA)
    library(readr)
    library(readxl)
    library(sf)
    library(stringr)
})

Sys.setenv(INLA_ST_RUN = '0')
source('src/models/fit_inla_spatiotemporal_screen.R')

root <- normalizePath('.', winslash = '/', mustWork = TRUE)
out_dir <- file.path(root, 'output', 'initial_forecast_2025')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
normalise_id <- function(x) toupper(trimws(as.character(x)))
processed <- file.path(root, 'data', 'processed')
forecast_files <- setNames(file.path(
    processed,
    paste0('initial_forecast_rows_2025_', c('manta', 'ltmp', 'mmp'), '.rds')
), c('manta', 'ltmp', 'mmp'))
forecast_weather <- file.path(processed, 'era5_weather_reef_year_2025.csv')
if (any(!file.exists(c(forecast_files, forecast_weather)))) {
    stop('Run Rscript src/data/build_2025_initial_forecast_rows.R first')
}
if (any(data$event_year >= 2025L)) {
    stop('The locked initial-forecast training data include 2025 outcomes')
}

# Build the assessment with the same joins and transformations as the selected
# INLA screen. Missing 2025 IMOS chlorophyll/disease features are deliberately
# left for pre-2025 median imputation inside prepare_fold_predictors().
assessment_base <- load_joint_compound_rows(
    mortality_files = forecast_files, weather_file = forecast_weather,
    require_complete_pressures = FALSE
)
cyclone_interval_2025 <- assessment_base |>
    select(.mortality_row_id, ReefID, baseline_survey_date, survey_date) |>
    inner_join(
        cyclone_all |>
            transmute(
                ReefID, pressure_year = event_year,
                tc_min_distance_km, tc_nearest_name,
                tc_nearest_max_wind_ms,
                annual_wind_distance_index = pmax(
                    tc_nearest_max_wind_ms - 17, 0
                ) * exp(-pmin(tc_min_distance_km, 1000) / 100)
            ),
        by = 'ReefID', relationship = 'many-to-many'
    ) |>
    mutate(
        pressure_start = as.Date(paste0(pressure_year - 1L, '-11-01')),
        pressure_end = as.Date(paste0(pressure_year, '-04-30'))
    ) |>
    filter(pressure_end >= baseline_survey_date,
           pressure_start <= pmin(survey_date, as.Date('2025-04-30'))) |>
    group_by(.mortality_row_id) |>
    summarise(
        tc_interval_min_distance_km = min(tc_min_distance_km, na.rm = TRUE),
        tc_interval_wind_distance_index = max(
            annual_wind_distance_index, na.rm = TRUE
        ),
        tc_interval_peak_name = tc_nearest_name[
            which.max(replace(annual_wind_distance_index,
                             !is.finite(annual_wind_distance_index), -Inf))
        ][1],
        .groups = 'drop'
    )

assessment <- assessment_base |>
    left_join(new_features, by = c('ReefID', 'event_year' = 'year'),
              relationship = 'many-to-one') |>
    left_join(cyclone, by = c('ReefID', 'event_year'),
              relationship = 'many-to-one') |>
    left_join(cyclone_interval_2025, by = '.mortality_row_id',
              relationship = 'many-to-one') |>
    left_join(cots_hindcast, by = c('ReefName', 'event_year'),
              relationship = 'many-to-one') |>
    left_join(disease, by = c('ReefID', 'event_year' = 'year'),
              relationship = 'many-to-one') |>
    left_join(dhw_correction, by = c('ReefID', 'event_year'),
              relationship = 'many-to-one') |>
    left_join(local_correction, by = c('ReefID', 'event_year'),
              relationship = 'many-to-one') |>
    mutate(
        tc_proximity100 = exp(-pmin(tc_min_distance_km, 1000) / 100),
        tc_wind_distance_index = pmax(tc_nearest_max_wind_ms - 17, 0) *
            tc_proximity100,
        tc_interval_proximity100 = exp(
            -pmin(tc_interval_min_distance_km, 1000) / 100
        ),
        correction_source = coalesce(correction_source, 'none'),
        correction_sd = coalesce(correction_sd, 0),
        correction_uncertainty_method = coalesce(
            correction_uncertainty_method, 'no_2025_local_logger_correction'
        ),
        nearest_local_logger_km = coalesce(nearest_local_logger_km, Inf),
        effective_local_loggers = coalesce(effective_local_loggers, 0),
        local_loggers_used = coalesce(local_loggers_used, 0),
        cots_outbreak_probability = coalesce(cots_outbreak_probability, 0),
        event_year = 2025L,
        reef_event_key = paste(ReefID, event_year, sep = '__')
    )

assessment_sf <- assessment |>
    distinct(ReefID, lon, lat) |>
    st_as_sf(coords = c('lon', 'lat'), crs = 4326, remove = FALSE) |>
    st_transform(3577)
assessment_xy <- st_coordinates(assessment_sf) / 1000
assessment_locations <- tibble(
    ReefID = assessment_sf$ReefID,
    x_km = assessment_xy[, 1], y_km = assessment_xy[, 2]
)
assessment <- assessment |>
    select(-any_of(c('x_km', 'y_km'))) |>
    left_join(assessment_locations, by = 'ReefID', relationship = 'many-to-one')

# Add an explicit prediction-only sixth event level. This changes neither the
# training rows nor model specification; the constrained event effect receives
# its zero-centred new-level prior.
event_years <- c(2016L, 2017L, 2020L, 2022L, 2024L, 2025L)
assessment$event_index <- match(assessment$event_year, event_years)

thermal_core_id <-
    'persistent_rw1_local_dhw_decomposed_hazards_freshwater_partial_pool'
thermal_row <- candidates |> filter(candidate == thermal_core_id)
if (nrow(thermal_row) != 1L) stop('Selected thermal candidate is unavailable')
cause_terms <- c(
    'log1p_cyc_interval_maxHrs4mw_z',
    'tc_interval_proximity100_z', 'tc_interval_wind_distance_index_z',
    'log1p_cot_interval_idw_max_z', 'cots_outbreak_probability_z'
)
thermal_terms <- setdiff(candidate_shared_terms(thermal_row), cause_terms)
thermal_candidate <- thermal_row |>
    mutate(candidate = 'locked_initial_2025_thermal',
           feature_set = 'thermal_freshwater')
base_candidate_shared_terms <- candidate_shared_terms
candidate_shared_terms <- function(candidate) {
    if (identical(as.character(candidate$feature_set), 'thermal_freshwater')) {
        thermal_terms
    } else base_candidate_shared_terms(candidate)
}
nonthermal_training <- coalesce(data$disturbance_has_cots, FALSE) |
    coalesce(data$disturbance_has_cyclone, FALSE) |
    coalesce(data$disturbance_has_flood, FALSE)
thermal_fit <- fit_candidate(
    data[!nonthermal_training, , drop = FALSE], assessment,
    thermal_candidate, compute_criteria = FALSE, seed = 20250910L
)
saveRDS(thermal_fit, file.path(out_dir, 'locked_thermal_fit.rds'))
thermal_predictions <- thermal_fit$predictions |>
    rename(
        thermal_prediction = predicted_mortality,
        thermal_occurrence = predicted_occurrence,
        thermal_magnitude = predicted_positive_mortality
    )

# Cause-labelled annual transitions are supervised only by pre-2025 rows.
cots_hindcast_annual <- read_csv(
    file.path(root, 'data', 'gbrPredsAdj_20262408.csv'), show_col_types = FALSE
) |>
    transmute(
        ReefName = reefName, event_year = as.integer(year),
        cots_outbreak_probability = pmin(pmax(as.numeric(outbrProb), 0), 1)
    )
cyclone_annual <- read_csv(
    file.path(processed, 'bom_cyclone_reef_year.csv'), show_col_types = FALSE
) |>
    transmute(
        ReefID, event_year,
        lon = reef_longitude, lat = reef_latitude,
        tc_proximity100 = exp(-pmin(tc_min_distance_km, 1000) / 100),
        tc_wind_distance_index = pmax(tc_nearest_max_wind_ms - 17, 0) *
            tc_proximity100
    )
annual <- read_csv(
    file.path(processed, 'annual_coral_transitions.csv'), show_col_types = FALSE
) |>
    filter(event_year < 2025L) |>
    left_join(cots_hindcast_annual, by = c('ReefName', 'event_year'),
              relationship = 'many-to-one') |>
    left_join(cyclone_annual, by = c('ReefID', 'event_year'),
              relationship = 'many-to-one') |>
    mutate(
        relative_loss = pmin(pmax(
            -cover_change_pp / (100 * pmax(pre_cover, .02)), 0
        ), .999),
        cots_label = coalesce(disturbance_has_cots, FALSE) |
            str_detect(str_to_lower(coalesce(disturbance_text, '')),
                       'cots|crown-of-thorns'),
        cyclone_label = coalesce(disturbance_has_cyclone, FALSE) |
            coalesce(disturbance_has_flood, FALSE) |
            str_detect(str_to_lower(coalesce(disturbance_text, '')),
                       'cyclone|storm|flood'),
        cots_positive = cots_label & relative_loss > .02,
        cyclone_positive = cyclone_label & relative_loss > .02,
        wave_log = log1p(pmax(cyc_maxHrs4mw, 0)),
        wind_distance = tc_wind_distance_index,
        rain_log = log_coastal_rain30,
        pre_cover = pmin(pmax(pre_cover, .001), 1),
        acropora = pmin(pmax(prop_acropora_pre, 0), 1),
        is_manta = as.numeric(programme_key == 'manta'),
        is_mmp = as.numeric(programme_key == 'mmp'),
        cots_outbreak_probability = coalesce(cots_outbreak_probability, 0),
        rrn_event_excess_raw = pmax(coalesce(cot_idwmeanpertow, 0) - .22, 0)
    )

assessment_cause <- assessment |>
    transmute(
        source_observation_id = as.character(source_observation_id),
        programme_key, ReefID, ReefName, event_year,
        observed_mortality = mortality_prop,
        observed_occurrence = as.numeric(mortality_prop > 0),
        wave_log = log1p(pmax(cyc_interval_maxHrs4mw, 0)),
        wind_distance = tc_interval_wind_distance_index,
        rain_log = log_coastal_rain30,
        pre_cover = pmin(pmax(observed_pre_cover, .001), 1),
        acropora = pmin(pmax(prop_acropora_pre, 0), 1),
        lon, lat,
        is_manta = as.numeric(programme_key == 'manta'),
        is_mmp = as.numeric(programme_key == 'mmp'),
        cots_interval_max = cot_interval_idw_max,
        cots_outbreak_probability = coalesce(cots_outbreak_probability, 0),
        cyclone_wave_hours = coalesce(cyc_interval_maxHrs4mw, 0),
        cyclone_name = coalesce(tc_interval_peak_name, 'none'),
        ann_maxdhw
    )

prepare_features <- function(training, prediction, features) {
    for (feature in features) {
        tr <- as.numeric(training[[feature]])
        av <- as.numeric(prediction[[feature]])
        replacement <- median(tr[is.finite(tr)], na.rm = TRUE)
        if (!is.finite(replacement)) replacement <- 0
        tr[!is.finite(tr)] <- replacement
        av[!is.finite(av)] <- replacement
        centre <- mean(tr)
        spread <- sd(tr)
        if (!is.finite(spread) || spread < 1e-8) spread <- 1
        training[[paste0(feature, '_z')]] <- (tr - centre) / spread
        prediction[[paste0(feature, '_z')]] <- (av - centre) / spread
    }
    list(training = training, prediction = prediction)
}

fit_hurdle <- function(training, prediction, positive_field, features) {
    prepared <- prepare_features(training, prediction, features)
    training <- prepared$training
    prediction <- prepared$prediction
    terms <- paste0(features, '_z')
    formula <- as.formula(paste('response ~ 1 +', paste(terms, collapse = ' + ')))
    occurrence_data <- bind_rows(
        training |> transmute(
            response = as.numeric(.data[[positive_field]]),
            across(all_of(terms))
        ),
        prediction |> transmute(response = NA_real_, across(all_of(terms)))
    )
    occurrence_fit <- inla(
        formula, family = 'binomial', data = occurrence_data,
        control.predictor = list(compute = TRUE, link = 1), verbose = FALSE
    )
    oi <- seq.int(nrow(training) + 1L, nrow(occurrence_data))
    positive <- training |>
        filter(.data[[positive_field]]) |>
        transmute(response = pmin(pmax(relative_loss, .001), .999),
                  across(all_of(terms)))
    magnitude_data <- bind_rows(
        positive,
        prediction |> transmute(response = NA_real_, across(all_of(terms)))
    )
    magnitude_fit <- inla(
        formula, family = 'beta', data = magnitude_data,
        control.predictor = list(compute = TRUE, link = 1), verbose = FALSE
    )
    mi <- seq.int(nrow(positive) + 1L, nrow(magnitude_data))
    list(
        prediction = occurrence_fit$summary.fitted.values$mean[oi] *
            magnitude_fit$summary.fitted.values$mean[mi],
        occurrence = occurrence_fit$summary.fitted.values$mean[oi],
        magnitude = magnitude_fit$summary.fitted.values$mean[mi]
    )
}

cyclone_features <- c(
    'wave_log', 'wind_distance', 'rain_log', 'pre_cover',
    'acropora', 'lon', 'lat', 'is_manta', 'is_mmp'
)
cyclone_fit <- fit_hurdle(
    annual, assessment_cause, 'cyclone_positive', cyclone_features
)

# Exact-date Manta COTS state is frozen at 1 March in every event.
workbook <- file.path(
    root, 'data', '250929_COTS-Manta-Cull-RHIS-Data-Matthews-and-Schlawinsky.xlsx'
)
manta <- suppressWarnings(read_excel(
    workbook, sheet = 'Manta', guess_max = 50000
)) |>
    transmute(
        ReefID = normalise_id(ReefLabel),
        survey_date = as.Date(substr(as.character(SurveyTime), 1, 10)),
        cots_count = pmax(as.numeric(CrownOfThornsStarfishCount), 0)
    ) |>
    filter(!is.na(ReefID), !is.na(survey_date), is.finite(cots_count)) |>
    group_by(ReefID, survey_date) |>
    summarise(manta_density = mean(cots_count), manta_tows = n(), .groups = 'drop')
manta_split <- split(manta, manta$ReefID)

feature_grid <- bind_rows(
    annual |> select(ReefID, event_year),
    assessment_cause |> select(ReefID, event_year)
) |>
    distinct(ReefID, event_year) |>
    mutate(issue_date = as.Date(paste0(event_year, '-03-01')))
temporal_features <- bind_rows(lapply(seq_len(nrow(feature_grid)), function(i) {
    reef_id <- feature_grid$ReefID[[i]]
    issue_date <- feature_grid$issue_date[[i]]
    observations <- manta_split[[reef_id]]
    if (is.null(observations)) observations <- manta[0, ]
    observations <- observations |> filter(survey_date < issue_date)
    latest_date <- if (nrow(observations)) max(observations$survey_date) else as.Date(NA)
    latest <- if (nrow(observations)) {
        observations |> filter(survey_date == latest_date)
    } else observations
    prior3 <- observations |> filter(survey_date >= issue_date - 3 * 365)
    peak_date <- if (nrow(prior3)) {
        prior3$survey_date[[which.max(prior3$manta_density)]]
    } else as.Date(NA)
    tibble(
        ReefID = reef_id, event_year = feature_grid$event_year[[i]],
        manta_supported = as.numeric(nrow(latest) > 0),
        manta_latest_excess_raw = if (nrow(latest)) {
            pmax(mean(latest$manta_density) - .22, 0)
        } else 0,
        manta_peak_excess_raw = if (nrow(prior3)) {
            pmax(max(prior3$manta_density) - .22, 0)
        } else 0,
        manta_years_since_peak = if (nrow(prior3)) {
            pmin(as.numeric(issue_date - peak_date) / 365.25, 10)
        } else 10
    )
}))
annual <- annual |>
    left_join(temporal_features, by = c('ReefID', 'event_year'),
              relationship = 'many-to-one')
assessment_cause <- assessment_cause |>
    left_join(temporal_features, by = c('ReefID', 'event_year'),
              relationship = 'many-to-one') |>
    mutate(rrn_event_excess_raw = pmax(coalesce(cots_interval_max, 0) - .22, 0))
cots_features <- c(
    'rrn_event_excess_raw', 'manta_latest_excess_raw',
    'manta_peak_excess_raw', 'manta_years_since_peak',
    'cots_outbreak_probability', 'pre_cover', 'acropora',
    'lon', 'lat', 'is_manta', 'is_mmp', 'manta_supported'
)
cots_fit <- fit_hurdle(annual, assessment_cause, 'cots_positive', cots_features)

locked <- thermal_predictions |>
    select(
        programme_key, source_observation_id, ReefID, ReefName, event_year,
        observed_mortality, observed_occurrence,
        thermal_prediction, thermal_occurrence, thermal_magnitude,
        ann_maxdhw, log_coastal_rain30, wqc_freqcc12,
        wqc_prior10_percentile, wqc_10yr_sum
    ) |>
    left_join(
        assessment_cause |>
            select(source_observation_id, lon, lat, pre_cover, acropora,
                   cots_interval_max, cots_outbreak_probability,
                   cyclone_wave_hours, cyclone_name, rrn_event_excess_raw,
                   manta_supported, manta_latest_excess_raw,
                   manta_peak_excess_raw, manta_years_since_peak),
        by = 'source_observation_id', relationship = 'one-to-one'
    ) |>
    mutate(
        cots_prediction = cots_fit$prediction,
        cots_occurrence = cots_fit$occurrence,
        cots_activation = pmax(
            cots_outbreak_probability,
            1 - exp(-pmax(rrn_event_excess_raw, 0) / .5)
        ),
        cyclone_prediction = cyclone_fit$prediction,
        cyclone_occurrence = cyclone_fit$occurrence,
        cyclone_activation = plogis((cyclone_wave_hours - 20) / 5),
        predicted_mortality = 1 - (1 - thermal_prediction) *
            (1 - cots_prediction * cots_activation) *
            (1 - cyclone_prediction * cyclone_activation),
        predicted_occurrence = 1 - (1 - thermal_occurrence) *
            (1 - cots_occurrence * cots_activation) *
            (1 - cyclone_occurrence * cyclone_activation),
        base_conditional_magnitude = pmin(
            predicted_mortality / pmax(predicted_occurrence, .001), .999
        ),
        residual = observed_mortality - predicted_mortality,
        candidate = 'operational_rrn_raw_plus_manta_state',
        forecast_stage = 'initial_environmental_forecast',
        fit_outcomes_through = 2024L,
        aerial_or_rhis_used = FALSE,
        reconstruction_status = 'retrospectively_locked_before_update_fit'
    )
if (any(!is.finite(locked$predicted_mortality)) ||
    any(!between(locked$predicted_mortality, 0, 1))) {
    stop('Locked 2025 predictions are invalid')
}
write_csv(locked, file.path(out_dir, 'locked_predictions.csv'))

metrics <- locked |>
    summarise(
        n = n(), reefs = n_distinct(ReefID),
        rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
        mae = mean(abs(observed_mortality - predicted_mortality)),
        bias = mean(predicted_mortality - observed_mortality),
        occurrence_brier = mean((observed_occurrence - predicted_occurrence)^2),
        false_extreme_rate = mean(
            predicted_mortality[observed_mortality < .3] >= .3
        ),
        observed_mean = mean(observed_mortality),
        predicted_mean = mean(predicted_mortality),
        observed_occurrence = mean(observed_occurrence),
        predicted_occurrence = mean(predicted_occurrence)
    )
write_csv(metrics, file.path(out_dir, 'locked_metrics.csv'))

manifest <- tibble(
    contract_version = 1L,
    event_year = 2025L,
    selected_model_id = 'operational_rrn_raw_plus_manta_state',
    training_rows = nrow(data),
    maximum_training_event = max(data$event_year),
    assessment_rows = nrow(locked),
    assessment_responses_enter_fit = FALSE,
    aerial_or_rhis_used = FALSE,
    cots_state_issue_date = '2025-03-01',
    prediction_file = 'output/initial_forecast_2025/locked_predictions.csv',
    prediction_md5 = as.character(tools::md5sum(
        file.path(out_dir, 'locked_predictions.csv')
    )),
    row_input_audit_md5 = as.character(tools::md5sum(file.path(
        processed, 'initial_forecast_2025_input_audit.csv'
    ))),
    model_code_md5 = as.character(tools::md5sum(
        file.path(root, 'src', 'models', 'fit_inla_spatiotemporal_screen.R')
    )),
    caveat = paste(
        'This is a retrospective reconstruction with inputs locked before',
        'the aerial/RHIS update fit, not an archived forecast issued in 2025.'
    )
)
write_csv(manifest, file.path(out_dir, 'locked_manifest.csv'))
message('Locked 2025 initial forecast written to: ', out_dir)
