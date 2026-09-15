# Audit where the selected model's leave-one-event-out prediction error remains.
# This script describes cross-fitted residuals; it does not fit or select a model.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

prediction_file <- 'output/prospective_cots_nowcast/cv_predictions.csv'
context_file <- 'output/explanatory_event_dhw/event_dhw_brt_data.csv'
output_dir <- 'output/selected_model_residual_gap_audit'
selected_candidate <- 'operational_rrn_raw_plus_manta_state'

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(prediction_file) || !file.exists(context_file)) {
    stop('Run the selected-model pipeline before the residual-gap audit.')
}

normalise_id <- function(x) toupper(trimws(as.character(x)))

predictions <- read_csv(prediction_file, show_col_types = FALSE) |>
    filter(
        candidate == selected_candidate,
        scheme == 'leave_one_event_out'
    ) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        ReefID = normalise_id(ReefID),
        residual = observed_mortality - predicted_mortality,
        absolute_error = abs(residual),
        squared_error = residual^2,
        severe_observed = observed_mortality >= 0.5,
        cots_effective_hazard = coalesce(
            cots_prediction_nowcast * cots_activation_nowcast, 0
        ),
        cyclone_effective_hazard = coalesce(
            cyclone_prediction * cyclone_activation, 0
        ),
        predicted_without_cots = 1 - (1 - thermal_prediction) *
            (1 - cyclone_effective_hazard),
        predicted_without_cyclone = 1 - (1 - thermal_prediction) *
            (1 - cots_effective_hazard),
        cots_marginal_uplift = predicted_mortality - predicted_without_cots,
        cyclone_marginal_uplift =
            predicted_mortality - predicted_without_cyclone
    )

if (nrow(predictions) == 0L) {
    stop('No selected leave-one-event-out predictions were found.')
}
if (anyDuplicated(predictions$source_observation_id)) {
    stop('Selected leave-one-event-out predictions are not observation-unique.')
}

context <- read_csv(context_file, show_col_types = FALSE) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        ReefID = normalise_id(ReefID)
    ) |>
    distinct(programme_key, source_observation_id, .keep_all = TRUE) |>
    select(
        programme_key, source_observation_id, region_block, SECTOR,
        DISTURBANCE_TYPE, disturbance_text,
        disturbance_has_bleaching, disturbance_has_cyclone,
        disturbance_has_flood, disturbance_has_cots,
        prop_acropora_pre, acropora_cover_pre,
        observed_pre_cover, log_coastal_rain30, wqc_freqcc12,
        wqc_prior10_percentile, wqc_10yr_sum, secc3m_p10, cloudp_90,
        mcur_90, sst_summer_skewness, sst_summer_excess_kurtosis,
        chla_wetseason_median, cyc_maxHrs4mw,
        cyc_interval_maxHrs4mw, tc_interval_min_distance_km,
        tc_interval_wind_distance_index, tc_interval_peak_name,
        cot_idwmeanpertow, cot_interval_idw_max
    )

rows <- predictions |>
    left_join(
        context,
        by = c('programme_key', 'source_observation_id'),
        relationship = 'many-to-one'
    ) |>
    mutate(
        explicit_cause_count = rowSums(cbind(
            coalesce(disturbance_has_flood, FALSE),
            coalesce(disturbance_has_cyclone, FALSE),
            coalesce(disturbance_has_cots, FALSE)
        )),
        cause_group = case_when(
            explicit_cause_count > 1L ~ 'Multiple non-thermal labels',
            coalesce(disturbance_has_flood, FALSE) ~ 'Freshwater/flood label',
            coalesce(disturbance_has_cyclone, FALSE) ~ 'Cyclone label',
            coalesce(disturbance_has_cots, FALSE) ~ 'COTS label',
            coalesce(disturbance_has_bleaching, FALSE) ~ 'Bleaching label only',
            TRUE ~ 'Unspecified/background'
        )
    )

if (anyNA(rows$region_block)) {
    stop('Region context is missing after joining selected predictions.')
}

# Event-external thresholds prevent an extreme event defining its own outliers.
threshold_keys <- rows |>
    distinct(programme_key, event_year) |>
    rowwise() |>
    mutate(
        calibration_n = sum(
            rows$programme_key == programme_key &
                rows$event_year != event_year
        ),
        residual_q05 = quantile(
            rows$residual[
                rows$programme_key == programme_key &
                    rows$event_year != event_year
            ],
            0.05, na.rm = TRUE, names = FALSE
        ),
        residual_q95 = quantile(
            rows$residual[
                rows$programme_key == programme_key &
                    rows$event_year != event_year
            ],
            0.95, na.rm = TRUE, names = FALSE
        )
    ) |>
    ungroup()

rows <- rows |>
    left_join(threshold_keys, by = c('programme_key', 'event_year')) |>
    mutate(
        underprediction_outlier = residual > pmax(0.10, residual_q95),
        overprediction_outlier = residual < pmin(-0.10, residual_q05)
    )

global_sse <- sum(rows$squared_error)

error_summary <- function(data, groups) {
    data |>
        group_by(across(all_of(groups))) |>
        summarise(
            observations = n(),
            reefs = n_distinct(ReefID),
            observed_mean = mean(observed_mortality),
            predicted_mean = mean(predicted_mortality),
            bias_observed_minus_predicted = mean(residual),
            rmse = sqrt(mean(squared_error)),
            mae = mean(absolute_error),
            severe_observations = sum(severe_observed),
            severe_rmse = if (any(severe_observed)) sqrt(mean(
                squared_error[severe_observed]
            )) else NA_real_,
            underprediction_outliers = sum(underprediction_outlier),
            overprediction_outliers = sum(overprediction_outlier),
            squared_error = sum(squared_error),
            squared_error_share = squared_error / global_sse,
            thermal_prediction_mean = mean(thermal_prediction),
            cots_marginal_uplift_mean = mean(cots_marginal_uplift),
            cyclone_marginal_uplift_mean = mean(cyclone_marginal_uplift),
            .groups = 'drop'
        ) |>
        arrange(desc(squared_error))
}

overview <- tibble(
    candidate = selected_candidate,
    scheme = 'leave_one_event_out',
    observations = nrow(rows),
    reefs = n_distinct(rows$ReefID),
    events = n_distinct(rows$event_year),
    rmse = sqrt(mean(rows$squared_error)),
    mae = mean(rows$absolute_error),
    predictive_r2 = 1 - sum(rows$squared_error) /
        sum((rows$observed_mortality - mean(rows$observed_mortality))^2),
    severe_rmse = sqrt(mean(rows$squared_error[rows$severe_observed])),
    underprediction_outliers = sum(rows$underprediction_outlier),
    overprediction_outliers = sum(rows$overprediction_outlier)
)

write_csv(overview, file.path(output_dir, 'overview.csv'))
write_csv(error_summary(rows, 'event_year'),
          file.path(output_dir, 'event_summary.csv'))
write_csv(error_summary(rows, 'programme_key'),
          file.path(output_dir, 'programme_summary.csv'))
write_csv(error_summary(rows, 'region_block'),
          file.path(output_dir, 'region_summary.csv'))
write_csv(error_summary(rows, c('event_year', 'region_block')),
          file.path(output_dir, 'event_region_summary.csv'))
write_csv(error_summary(rows, 'cause_group'),
          file.path(output_dir, 'cause_summary.csv'))

reef_summary <- rows |>
    group_by(ReefID, ReefName, region_block) |>
    summarise(
        observations = n(),
        events = n_distinct(event_year),
        programmes = paste(sort(unique(programme_key)), collapse = '; '),
        observed_mean = mean(observed_mortality),
        predicted_mean = mean(predicted_mortality),
        mean_residual = mean(residual),
        rmse = sqrt(mean(squared_error)),
        maximum_absolute_error = max(absolute_error),
        underprediction_outliers = sum(underprediction_outlier),
        overprediction_outliers = sum(overprediction_outlier),
        squared_error = sum(squared_error),
        squared_error_share = squared_error / global_sse,
        .groups = 'drop'
    ) |>
    arrange(desc(squared_error))
write_csv(reef_summary, file.path(output_dir, 'reef_summary.csv'))

# Equal-weighted reef-event summaries prevent repeatedly sampled reef-events
# from automatically becoming the highest priority.
reef_events <- rows |>
    group_by(ReefID, ReefName, event_year, region_block) |>
    summarise(
        observations = n(),
        programmes = paste(sort(unique(programme_key)), collapse = '; '),
        observed_mortality = mean(observed_mortality),
        predicted_mortality = mean(predicted_mortality),
        residual = observed_mortality - predicted_mortality,
        absolute_error = abs(residual),
        rmse_across_observations = sqrt(mean(squared_error)),
        severe_observed = observed_mortality >= 0.5,
        disturbance_text = paste(
            sort(unique(na.omit(disturbance_text))), collapse = '; '
        ),
        flood_label = any(coalesce(disturbance_has_flood, FALSE)),
        cyclone_label = any(coalesce(disturbance_has_cyclone, FALSE)),
        cots_label = any(coalesce(disturbance_has_cots, FALSE)),
        ann_maxdhw = mean(ann_maxdhw, na.rm = TRUE),
        log_coastal_rain30 = mean(log_coastal_rain30, na.rm = TRUE),
        wqc_freqcc12 = mean(wqc_freqcc12, na.rm = TRUE),
        wqc_prior10_percentile = mean(
            wqc_prior10_percentile, na.rm = TRUE
        ),
        cyclone_wave_hours = mean(cyclone_wave_hours, na.rm = TRUE),
        cyclone_activation = mean(cyclone_activation, na.rm = TRUE),
        cyclone_name = paste(sort(unique(na.omit(cyclone_name))),
                             collapse = '; '),
        cyclone_track_peak = paste(
            sort(unique(na.omit(tc_interval_peak_name))), collapse = '; '
        ),
        cots_activation = mean(cots_activation_nowcast, na.rm = TRUE),
        cots_pressure = mean(rrn_event_excess_raw, na.rm = TRUE),
        prop_acropora_pre = mean(prop_acropora_pre, na.rm = TRUE),
        thermal_prediction = mean(thermal_prediction),
        cots_marginal_uplift = mean(cots_marginal_uplift),
        cyclone_marginal_uplift = mean(cyclone_marginal_uplift),
        .groups = 'drop'
    ) |>
    group_by(event_year) |>
    mutate(
        rain_event_percentile = percent_rank(log_coastal_rain30),
        water_colour_event_percentile = percent_rank(wqc_freqcc12),
        cyclone_wave_event_percentile = percent_rank(cyclone_wave_hours),
        cots_pressure_event_percentile = percent_rank(cots_pressure)
    ) |>
    ungroup() |>
    mutate(
        high_rain = coalesce(rain_event_percentile >= 0.8, FALSE),
        high_water_colour = coalesce(
            water_colour_event_percentile >= 0.8, FALSE
        ),
        activated_cyclone = coalesce(cyclone_activation >= 0.5, FALSE),
        activated_cots = coalesce(cots_activation >= 0.5, FALSE),
        high_heat = coalesce(ann_maxdhw >= 8, FALSE)
    ) |>
    rowwise() |>
    mutate(
        candidate_mechanism_flags = paste(c(
            if (flood_label) 'explicit flood/freshwater label',
            if (cyclone_label) 'explicit cyclone label',
            if (cots_label) 'explicit COTS label',
            if (high_rain) 'upper-quintile event rainfall proxy',
            if (high_water_colour) 'upper-quintile current water colour',
            if (activated_cyclone) 'cyclone hazard activated',
            if (activated_cots) 'COTS hazard activated',
            if (high_heat) 'at least 8 DHW'
        ), collapse = '; '),
        candidate_mechanism_flags = if_else(
            candidate_mechanism_flags == '',
            'no strong measured mechanism flag',
            candidate_mechanism_flags
        )
    ) |>
    ungroup()

reef_event_sse <- sum(reef_events$residual^2)
reef_events <- reef_events |>
    mutate(
        squared_error = residual^2,
        squared_error_share = squared_error / reef_event_sse
    ) |>
    arrange(desc(squared_error))
write_csv(reef_events, file.path(output_dir, 'reef_event_summary.csv'))

error_decomposition_groups <- rows |>
    group_by(ReefID, event_year) |>
    summarise(
        observations = n(),
        mean_residual = mean(residual),
        between_reef_event_error = observations * mean_residual^2,
        within_reef_event_disagreement = sum(
            (residual - mean_residual)^2
        ),
        .groups = 'drop'
    )
error_decomposition <- error_decomposition_groups |>
    summarise(
        observations = sum(observations),
        reef_events = n(),
        total_squared_error = sum(
            between_reef_event_error + within_reef_event_disagreement
        ),
        reef_event_mean_error = sum(between_reef_event_error),
        within_reef_event_disagreement = sum(
            within_reef_event_disagreement
        ),
        reef_event_mean_error_share =
            reef_event_mean_error / total_squared_error,
        within_reef_event_disagreement_share =
            within_reef_event_disagreement / total_squared_error
    )
write_csv(error_decomposition,
          file.path(output_dir, 'error_decomposition.csv'))

programme_reef_events <- rows |>
    group_by(ReefID, ReefName, event_year, region_block, programme_key) |>
    summarise(
        observations = n(),
        observed_mortality = mean(observed_mortality),
        predicted_mortality = mean(predicted_mortality),
        residual = mean(residual),
        .groups = 'drop'
    )
programme_discordance <- programme_reef_events |>
    group_by(ReefID, ReefName, event_year, region_block) |>
    filter(n_distinct(programme_key) >= 2L) |>
    summarise(
        programmes = paste(sort(programme_key), collapse = '; '),
        observed_minimum = min(observed_mortality),
        observed_maximum = max(observed_mortality),
        observed_range = observed_maximum - observed_minimum,
        residual_minimum = min(residual),
        residual_maximum = max(residual),
        residual_range = residual_maximum - residual_minimum,
        residual_directions_disagree =
            residual_minimum < 0 & residual_maximum > 0,
        .groups = 'drop'
    ) |>
    arrange(desc(observed_range))
write_csv(programme_discordance,
          file.path(output_dir, 'cross_programme_discordance.csv'))

reef_persistence <- reef_events |>
    group_by(ReefID, ReefName, region_block) |>
    filter(n() >= 2L) |>
    summarise(
        event_years = paste(sort(event_year), collapse = '; '),
        events = n(),
        mean_event_residual = mean(residual),
        event_residual_rmse = sqrt(mean(residual^2)),
        maximum_event_absolute_error = max(absolute_error),
        underpredicted_events = sum(residual > 0.10),
        overpredicted_events = sum(residual < -0.10),
        same_direction_fraction = max(
            mean(residual > 0), mean(residual < 0)
        ),
        .groups = 'drop'
    ) |>
    arrange(desc(event_residual_rmse))
write_csv(reef_persistence,
          file.path(output_dir, 'reef_persistence_summary.csv'))

priority_reef_events <- bind_rows(
    reef_events |>
        slice_max(residual, n = 20, with_ties = FALSE) |>
        mutate(direction = 'underpredicted'),
    reef_events |>
        slice_min(residual, n = 12, with_ties = FALSE) |>
        mutate(direction = 'overpredicted')
) |>
    arrange(factor(direction, c('underpredicted', 'overpredicted')),
            desc(absolute_error))
write_csv(priority_reef_events,
          file.path(output_dir, 'priority_reef_events.csv'))

safe_cor <- function(x, y, minimum_n = 8L) {
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < minimum_n || length(unique(x[keep])) < 3L ||
        length(unique(y[keep])) < 3L) {
        return(NA_real_)
    }
    suppressWarnings(cor(x[keep], y[keep], method = 'spearman'))
}

feature_catalogue <- tribble(
    ~feature, ~label, ~mechanism, ~model_status,
    'ann_maxdhw', 'Local-first maximum DHW', 'thermal', 'already represented',
    'log_coastal_rain30', 'Coastal 30-day rainfall', 'freshwater', 'already represented',
    'wqc_freqcc12', 'Current coloured-water frequency', 'freshwater', 'already represented',
    'wqc_prior10_percentile', 'Prior water-colour percentile', 'freshwater', 'already represented',
    'wqc_10yr_sum', 'Cumulative water colour', 'freshwater', 'already represented',
    'secc3m_p10', 'Low water clarity', 'optical water quality', 'already represented',
    'chla_wetseason_median', 'Wet-season chlorophyll', 'optical water quality', 'already represented',
    'cyc_interval_maxHrs4mw', 'Maximum damaging-wave hours', 'cyclone', 'already represented',
    'tc_interval_min_distance_km', 'Minimum cyclone-track distance', 'cyclone', 'already represented',
    'tc_interval_wind_distance_index', 'Cyclone wind-distance index', 'cyclone', 'already represented',
    'cot_interval_idw_max', 'Interval COTS density', 'COTS', 'already represented',
    'cots_outbreak_probability', 'COTS outbreak probability', 'COTS', 'already represented',
    'prop_acropora_pre', 'Pre-event Acropora proportion', 'susceptibility', 'already represented',
    'observed_pre_cover', 'Pre-event coral cover', 'susceptibility', 'already represented',
    'cloudp_90', 'Upper-tail cloud', 'heat modulation', 'already represented',
    'mcur_90', 'Upper-tail current speed', 'heat modulation', 'already represented',
    'sst_summer_skewness', 'Summer SST skewness', 'thermal shape', 'already represented',
    'sst_summer_excess_kurtosis', 'Summer SST excess kurtosis', 'thermal shape', 'already represented',
    'cyclone_wave_hours', 'Cyclone component wave hours', 'cyclone', 'already represented',
    'cyclone_activation', 'Cyclone activation', 'cyclone', 'derived model diagnostic',
    'cyclone_marginal_uplift', 'Cyclone marginal prediction uplift', 'cyclone', 'derived model diagnostic',
    'cots_activation_nowcast', 'COTS activation', 'COTS', 'derived model diagnostic',
    'cots_marginal_uplift', 'COTS marginal prediction uplift', 'COTS', 'derived model diagnostic'
)

association_for_feature <- function(feature_name) {
    x <- rows[[feature_name]]
    complete <- is.finite(x) & is.finite(rows$residual)
    ranked <- rows |>
        transmute(
            programme_key, event_year,
            feature_value = x,
            residual, absolute_error
        ) |>
        filter(is.finite(feature_value), is.finite(residual)) |>
        group_by(programme_key, event_year) |>
        filter(n() >= 8L) |>
        mutate(
            feature_rank = rank(feature_value, ties.method = 'average'),
            residual_rank = rank(residual, ties.method = 'average'),
            feature_rank_centered = feature_rank - mean(feature_rank),
            residual_rank_centered = residual_rank - mean(residual_rank)
        ) |>
        ungroup()
    grouped <- rows |>
        transmute(
            programme_key, event_year,
            feature_value = x, residual
        ) |>
        group_by(programme_key, event_year) |>
        summarise(
            n = sum(is.finite(feature_value) & is.finite(residual)),
            spearman_rho = safe_cor(feature_value, residual),
            .groups = 'drop'
        ) |>
        filter(is.finite(spearman_rho))
    tibble(
        feature = feature_name,
        n = sum(complete),
        overall_residual_spearman = safe_cor(x, rows$residual),
        overall_absolute_error_spearman = safe_cor(x, rows$absolute_error),
        pooled_within_event_rank_r = if (nrow(ranked) >= 8L) cor(
            ranked$feature_rank_centered,
            ranked$residual_rank_centered
        ) else NA_real_,
        event_programme_groups = nrow(grouped),
        median_group_spearman = if (nrow(grouped)) {
            median(grouped$spearman_rho)
        } else NA_real_,
        positive_groups = sum(grouped$spearman_rho > 0),
        negative_groups = sum(grouped$spearman_rho < 0)
    )
}

associations <- bind_rows(lapply(
    feature_catalogue$feature, association_for_feature
)) |>
    left_join(feature_catalogue, by = 'feature') |>
    select(
        feature, label, mechanism, model_status, everything()
    ) |>
    arrange(desc(abs(pooled_within_event_rank_r)))
write_csv(associations, file.path(output_dir, 'residual_associations.csv'))

association_groups <- bind_rows(lapply(feature_catalogue$feature, function(
    feature_name
) {
    rows |>
        transmute(
            programme_key, event_year,
            feature_value = .data[[feature_name]], residual
        ) |>
        group_by(programme_key, event_year) |>
        summarise(
            feature = feature_name,
            n = sum(is.finite(feature_value) & is.finite(residual)),
            spearman_rho = safe_cor(feature_value, residual),
            .groups = 'drop'
        ) |>
        filter(is.finite(spearman_rho))
})) |>
    left_join(feature_catalogue, by = 'feature') |>
    select(
        feature, label, mechanism, programme_key, event_year,
        n, spearman_rho
    ) |>
    arrange(feature, programme_key, event_year)
write_csv(association_groups,
          file.path(output_dir, 'residual_associations_by_group.csv'))

coverage <- rows |>
    group_by(event_year) |>
    summarise(
        observations = n(),
        rainfall_available = sum(is.finite(log_coastal_rain30)),
        current_water_colour_available = sum(is.finite(wqc_freqcc12)),
        prior_water_colour_available = sum(is.finite(wqc_prior10_percentile)),
        cyclone_wave_available = sum(is.finite(cyc_interval_maxHrs4mw)),
        cyclone_track_available = sum(is.finite(tc_interval_min_distance_km)),
        cots_density_available = sum(is.finite(cot_interval_idw_max)),
        direct_salinity_in_selected_model = 0L,
        flood_labels = sum(coalesce(disturbance_has_flood, FALSE)),
        cyclone_labels = sum(coalesce(disturbance_has_cyclone, FALSE)),
        cots_labels = sum(coalesce(disturbance_has_cots, FALSE)),
        .groups = 'drop'
    )
write_csv(coverage, file.path(output_dir, 'mechanism_coverage.csv'))

write_csv(
    rows |>
        select(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, region_block, observed_mortality,
            predicted_mortality, residual, absolute_error, squared_error,
            severe_observed, underprediction_outlier,
            overprediction_outlier, cause_group,
            thermal_prediction, cots_marginal_uplift,
            cyclone_marginal_uplift
        ),
    file.path(output_dir, 'selected_cross_fitted_residuals.csv')
)

print(overview)
print(read_csv(file.path(output_dir, 'event_summary.csv'),
               show_col_types = FALSE))
print(read_csv(file.path(output_dir, 'region_summary.csv'),
               show_col_types = FALSE))
