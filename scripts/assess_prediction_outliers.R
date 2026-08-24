# Identify ecologically informative prediction failures from event-held-out
# predictions. These are diagnostics, not automatic exclusion rules.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

prediction_file <- 'output/formal_models/model_comparison_predictions.csv'
output_dir <- 'output/prediction_outliers'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(prediction_file)) {
    stop('Run scripts/build_formal_ensemble.R before outlier assessment.')
}

predictions <- read_csv(prediction_file, show_col_types = FALSE) |>
    mutate(fold = as.character(fold))
binomial_file <- 'output/formal_models/binomial_blocked_predictions.csv'
if (file.exists(binomial_file)) {
    predictions <- bind_rows(
        predictions,
        read_csv(binomial_file, show_col_types = FALSE) |>
            mutate(fold = as.character(fold))
    )
}
predictions <- predictions |>
    filter(scheme == 'leave_one_event_out') |>
    mutate(source_observation_id = as.character(source_observation_id))

context <- bind_rows(lapply(c('manta', 'ltmp', 'mmp'), function(programme) {
    readRDS(paste0('data/processed/validation_rows_', programme, '.rds')) |>
        transmute(
            programme_key,
            source_observation_id = as.character(source_observation_id),
            ann_maxdhw,
            prop_acropora_pre,
            acropora_cover_pre,
            observed_pre_cover,
            post_cover,
            DISTURBANCE_TYPE,
            description,
            lon,
            lat
        )
}))

predictions <- predictions |>
    left_join(
        context,
        by = c('programme_key', 'source_observation_id'),
        relationship = 'many-to-one'
    ) |>
    mutate(
        era5_grid_lat = floor(lat / 0.25 + 0.5) * 0.25,
        era5_grid_lon = floor(lon / 0.25 + 0.5) * 0.25,
        residual = observed_mortality - predicted_mortality,
        absolute_error = abs(residual),
        severe_observed = observed_mortality >= 0.2,
        low_thermal_stress = ann_maxdhw < 4
    )

weather_file <- 'data/processed/era5_weather_reef_year.csv'
if (file.exists(weather_file)) {
    weather_context <- read_csv(weather_file, show_col_types = FALSE) |>
        group_by(year, grid_lat, grid_lon) |>
        summarise(
            era5_coastal_rain_max_30day = median(
                era5_coastal_rain_dec_mar_max_30day
            ),
            era5_coastal_rain_q1 = median(
                era5_coastal_rain_q1_total
            ),
            era5_coastal_distance_km = median(
                coastal_grid_distance_km
            ),
            era5_wind_mean = median(era5_reef_wind_q1_mean),
            era5_wind_calm_fraction = median(
                era5_reef_wind_q1_fraction_below_3
            ),
            .groups = 'drop'
        )
    predictions <- predictions |>
        left_join(
            weather_context,
            by = c(
                'event_year' = 'year',
                'era5_grid_lat' = 'grid_lat',
                'era5_grid_lon' = 'grid_lon'
            ),
            relationship = 'many-to-one'
        )
}

freshwater_risk_file <-
    'data/processed/ereefs_salinity_2024_rainfall_sensitivity.csv'
if (file.exists(freshwater_risk_file)) {
    freshwater_context <- read_csv(
        freshwater_risk_file, show_col_types = FALSE
    ) |>
        transmute(
            ReefID = as.character(LABEL_ID),
            event_year = as.integer(year),
            freshwater_risk30_percentile = risk30_percentile,
            partial_ereefs_flood30,
            partial_ereefs_salinity_min = salinity_min
        ) |>
        group_by(ReefID, event_year) |>
        summarise(
            freshwater_risk30_percentile = median(
                freshwater_risk30_percentile, na.rm = TRUE
            ),
            partial_ereefs_flood30 = any(
                partial_ereefs_flood30, na.rm = TRUE
            ),
            partial_ereefs_salinity_min = min(
                partial_ereefs_salinity_min, na.rm = TRUE
            ),
            .groups = 'drop'
        )
    predictions <- predictions |>
        left_join(
            freshwater_context,
            by = c('ReefID', 'event_year'),
            relationship = 'many-to-one'
        )
}

enso_file <- 'data/processed/enso_event_context.csv'
if (!file.exists(enso_file)) {
    stop('Run scripts/extract_enso_indices.py before outlier assessment.')
}
enso_context <- read_csv(enso_file, show_col_types = FALSE) |>
    select(
        event_year, roni_djf, roni_jfm, roni_fma,
        roni_bleaching_summer_mean, roni_bleaching_summer_max_abs,
        enso_phase, enso_strength, enso_category,
        soi_dec_mar_mean, soi_dec_mar_min, soi_dec_mar_max, soi_phase,
        enso_ocean_atmosphere_state
    )
predictions <- predictions |>
    left_join(enso_context, by = 'event_year', relationship = 'many-to-one')

# Threshold every event using only cross-fitted residuals from other events.
groups <- split(
    predictions,
    interaction(predictions$programme_key, predictions$learner, drop = TRUE)
)
assessed <- bind_rows(lapply(groups, function(group) {
    bind_rows(lapply(sort(unique(group$event_year)), function(year) {
        assessment <- filter(group, event_year == year)
        calibration <- filter(group, event_year != year)
        if (nrow(calibration) < 20) {
            calibration <- group
        }
        lower <- quantile(calibration$residual, 0.05, na.rm = TRUE)
        upper <- quantile(calibration$residual, 0.95, na.rm = TRUE)
        absolute <- quantile(calibration$absolute_error, 0.90, na.rm = TRUE)
        assessment |>
            mutate(
                calibration_residual_q05 = lower,
                calibration_residual_q95 = upper,
                calibration_absolute_error_q90 = absolute,
                residual_outlier = residual < lower | residual > upper,
                outlier_direction = case_when(
                    residual > upper ~ 'underestimated',
                    residual < lower ~ 'overestimated',
                    TRUE ~ 'within_other_event_90pct_range'
                ),
                severe_2024_underprediction =
                    event_year == 2024 & severe_observed &
                    residual > pmax(0.10, upper),
                low_dhw_overprediction =
                    low_thermal_stress & -residual > pmax(0.10, -lower)
            )
    }))
}))

core_consensus <- assessed |>
    filter(learner %in% c('brms', 'brt', 'formal_ensemble')) |>
    group_by(programme_key, source_observation_id, event_year) |>
    summarise(
        models_underestimating = sum(outlier_direction == 'underestimated'),
        models_overestimating = sum(outlier_direction == 'overestimated'),
        brms_residual = residual[learner == 'brms'][1],
        brt_residual = residual[learner == 'brt'][1],
        formal_ensemble_residual = residual[learner == 'formal_ensemble'][1],
        .groups = 'drop'
    )

assessed <- assessed |>
    left_join(
        core_consensus,
        by = c('programme_key', 'source_observation_id', 'event_year')
    ) |>
    mutate(
        core_model_consensus = case_when(
            models_underestimating == 3 ~ 'all_three_underestimate',
            models_overestimating == 3 ~ 'all_three_overestimate',
            TRUE ~ 'mixed_or_not_outlying'
        )
    )

primary <- assessed |>
    filter(learner == 'formal_ensemble') |>
    arrange(desc(absolute_error))

summary <- primary |>
    group_by(programme_key, event_year) |>
    summarise(
        n = n(),
        unique_reefs = n_distinct(ReefID),
        enso_category = first(enso_category),
        roni_bleaching_summer_mean = first(roni_bleaching_summer_mean),
        soi_dec_mar_mean = first(soi_dec_mar_mean),
        observed_mean = mean(observed_mortality),
        predicted_mean = mean(predicted_mortality),
        bias_observed_minus_predicted = mean(residual),
        mae = mean(absolute_error),
        underestimated_outliers = sum(outlier_direction == 'underestimated'),
        overestimated_outliers = sum(outlier_direction == 'overestimated'),
        severe_2024_underpredictions = sum(severe_2024_underprediction),
        severe_2024_unique_reefs = n_distinct(
            ReefID[severe_2024_underprediction]
        ),
        low_dhw_overpredictions = sum(low_dhw_overprediction),
        low_dhw_unique_reefs = n_distinct(
            ReefID[low_dhw_overprediction]
        ),
        .groups = 'drop'
    )

family_comparison <- assessed |>
    filter(learner %in% c('brms', 'brms_binomial_olre')) |>
    group_by(programme_key, learner, event_year) |>
    summarise(
        n = n(),
        observed_mean = mean(observed_mortality),
        predicted_mean = mean(predicted_mortality),
        bias_observed_minus_predicted = mean(residual),
        mae = mean(absolute_error),
        rmse = sqrt(mean(residual^2)),
        severe_mae = mean(
            absolute_error[observed_mortality >= 0.2], na.rm = TRUE
        ),
        low_dhw_bias = mean(
            residual[ann_maxdhw < 4], na.rm = TRUE
        ),
        .groups = 'drop'
    )

priority_rows <- primary |>
    filter(
        residual_outlier | severe_2024_underprediction |
            low_dhw_overprediction
    ) |>
    arrange(
        desc(severe_2024_underprediction),
        desc(low_dhw_overprediction),
        desc(absolute_error)
    )

# Several MMP reefs have multiple survey observations in an event. Preserve
# those in a full audit table, but provide one representative record per
# programme/reef/event for ecological review so counts are not inflated.
priority <- priority_rows |>
    group_by(programme_key, ReefID, event_year) |>
    mutate(
        n_flagged_observations = n(),
        reef_year_severe_2024_underprediction = any(
            severe_2024_underprediction
        ),
        reef_year_low_dhw_overprediction = any(low_dhw_overprediction),
        reef_year_residual_outlier = any(residual_outlier),
        reef_year_all_three_underestimate = any(
            core_model_consensus == 'all_three_underestimate'
        ),
        reef_year_all_three_overestimate = any(
            core_model_consensus == 'all_three_overestimate'
        )
    ) |>
    slice_max(absolute_error, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
        severe_2024_underprediction =
            reef_year_severe_2024_underprediction,
        low_dhw_overprediction = reef_year_low_dhw_overprediction,
        residual_outlier = reef_year_residual_outlier,
        core_model_consensus = case_when(
            reef_year_all_three_underestimate ~ 'all_three_underestimate',
            reef_year_all_three_overestimate ~ 'all_three_overestimate',
            TRUE ~ core_model_consensus
        )
    ) |>
    select(-starts_with('reef_year_')) |>
    arrange(
        desc(severe_2024_underprediction),
        desc(low_dhw_overprediction),
        desc(absolute_error)
    )

priority_key <- paste(
    priority$programme_key, priority$ReefID, priority$event_year,
    sep = '::'
)
if (anyDuplicated(priority_key)) {
    stop('Priority ecological review table is not unique by reef-year')
}

write_csv(assessed, file.path(output_dir, 'all_cross_fitted_residuals.csv'))
write_csv(primary, file.path(output_dir, 'formal_ensemble_residuals.csv'))
write_csv(
    priority_rows,
    file.path(output_dir, 'priority_observation_rows.csv')
)
write_csv(summary, file.path(output_dir, 'event_outlier_summary.csv'))
write_csv(
    family_comparison,
    file.path(output_dir, 'beta_binomial_family_comparison.csv')
)

priority_file <- file.path(output_dir, 'priority_reef_years.csv')
tryCatch(
    write_csv(priority, priority_file),
    error = function(error) {
        warning(
            'Could not refresh ', priority_file,
            '; the file may be open in another application. ',
            conditionMessage(error),
            call. = FALSE
        )
    }
)

print(summary)
cat('Priority observation rows:', nrow(priority_rows), '\n')
cat('Priority programme-reef-years:', nrow(priority), '\n')
