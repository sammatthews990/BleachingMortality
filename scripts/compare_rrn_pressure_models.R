# Compare the RRN pressure ablations with the existing core and compound
# weather BRTs on exactly matched outer-fold predictions.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

output_dir <- 'output/rrn_pressure_assessment'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prediction_files <- c(
    core = 'output/joint_compound_models/joint_brt_core_predictions.csv',
    weather = 'output/joint_compound_models/joint_brt_predictions.csv',
    rrn_wq = 'output/joint_compound_models/joint_brt_rrn_wq_predictions.csv',
    rrn_cyclone = 'output/joint_compound_models/joint_brt_rrn_cyclone_predictions.csv',
    rrn_both = 'output/joint_compound_models/joint_brt_rrn_predictions.csv',
    weather_rrn = 'output/joint_compound_models/joint_brt_compound_rrn_predictions.csv',
    weather_rrn_relative = paste0(
        'output/joint_compound_models/',
        'joint_brt_compound_rrn_relative_predictions.csv'
    ),
    weather_rrn_full = paste0(
        'output/joint_compound_models/',
        'joint_brt_compound_rrn_full_predictions.csv'
    )
)

missing_files <- prediction_files[!file.exists(prediction_files)]
if (length(missing_files) > 0L) {
    stop('Missing RRN comparison predictions: ', paste(missing_files, collapse = ', '))
}

rrn <- read_csv(
    'data/processed/rrn_pressure_reef_year.csv', show_col_types = FALSE
) |>
    select(
        LABEL_ID, event_year, wqc_freqcc12,
        wqc_prior10_percentile, wqc_prior10_delta,
        cyc_maxHrs4mw, log1p_cyc_maxHrs4mw,
        cot_meanpertow, cot_idwmeanpertow,
        log1p_cot_idwmeanpertow
    )

predictions <- bind_rows(lapply(names(prediction_files), function(model_name) {
    model_rows <- read_csv(
        prediction_files[[model_name]], show_col_types = FALSE
    )
    required_rrn <- setdiff(names(rrn), c('LABEL_ID', 'event_year'))
    if (!all(required_rrn %in% names(model_rows))) {
        model_rows <- model_rows |> select(-any_of(required_rrn))
        model_rows <- model_rows |>
            left_join(
                rrn,
                by = c('ReefID' = 'LABEL_ID', 'event_year'),
                relationship = 'many-to-one'
            )
    }
    model_rows |>
        mutate(
            model = model_name,
            source_observation_id = as.character(source_observation_id),
            fold = as.character(fold)
        )
}))

overall_metrics <- predictions |>
    group_by(model, scheme) |>
    summarise(
        observations = n(),
        rmse = sqrt(mean(residual^2)),
        mae = mean(absolute_error),
        severe_mae = mean(
            absolute_error[observed_mortality >= 0.2], na.rm = TRUE
        ),
        predictive_r2 = 1 - sum(residual^2) /
            sum((observed_mortality - mean(observed_mortality))^2),
        bias_observed_minus_predicted = mean(residual),
        .groups = 'drop'
    )
write_csv(overall_metrics, file.path(output_dir, 'mortality_brt_metrics.csv'))

programme_metrics <- predictions |>
    group_by(model, scheme, programme_key) |>
    summarise(
        observations = n(),
        rmse = sqrt(mean(residual^2)),
        mae = mean(absolute_error),
        severe_mae = mean(
            absolute_error[observed_mortality >= 0.2], na.rm = TRUE
        ),
        predictive_r2 = 1 - sum(residual^2) /
            sum((observed_mortality - mean(observed_mortality))^2),
        bias_observed_minus_predicted = mean(residual),
        .groups = 'drop'
    )
write_csv(
    programme_metrics,
    file.path(output_dir, 'mortality_brt_programme_metrics.csv')
)

key_columns <- c(
    'scheme', 'fold', 'programme_key', 'source_observation_id',
    'ReefID', 'ReefName', 'event_year', 'region_block',
    'observed_mortality', 'ann_maxdhw', 'prop_acropora_pre',
    'wqc_freqcc12', 'wqc_prior10_percentile', 'wqc_prior10_delta',
    'cyc_maxHrs4mw', 'cot_meanpertow', 'cot_idwmeanpertow',
    'DISTURBANCE_TYPE', 'storm_name', 'description',
    'tooltip', 'disturbance_text'
)

paired <- predictions |>
    select(all_of(key_columns), model, predicted_mortality, absolute_error) |>
    pivot_wider(
        names_from = model,
        values_from = c(predicted_mortality, absolute_error)
    ) |>
    mutate(
        wq_gain_vs_core = absolute_error_core - absolute_error_rrn_wq,
        cyclone_gain_vs_core = absolute_error_core - absolute_error_rrn_cyclone,
        rrn_gain_vs_core = absolute_error_core - absolute_error_rrn_both,
        weather_gain_vs_core = absolute_error_core - absolute_error_weather,
        weather_rrn_gain_vs_weather =
            absolute_error_weather - absolute_error_weather_rrn,
        relative_wq_gain_vs_weather_rrn =
            absolute_error_weather_rrn - absolute_error_weather_rrn_relative,
        cots_gain_vs_relative =
            absolute_error_weather_rrn_relative - absolute_error_weather_rrn_full,
        full_gain_vs_weather =
            absolute_error_weather - absolute_error_weather_rrn_full
    )
write_csv(paired, file.path(output_dir, 'mortality_brt_paired_rows.csv'))

strata <- paired |>
    mutate(
        severe = observed_mortality >= 0.2,
        cyclone_exposed = cyc_maxHrs4mw > 0,
        high_wq = wqc_freqcc12 >= 0.5,
        high_relative_wq = wqc_prior10_percentile >= 0.9,
        cots_pressure = cot_idwmeanpertow > 0
    ) |>
    group_by(
        scheme, programme_key, event_year, severe, cyclone_exposed,
        high_wq, high_relative_wq, cots_pressure
    ) |>
    summarise(
        observations = n(),
        mean_wq_gain_vs_core = mean(wq_gain_vs_core),
        mean_cyclone_gain_vs_core = mean(cyclone_gain_vs_core),
        mean_rrn_gain_vs_core = mean(rrn_gain_vs_core),
        mean_weather_rrn_gain_vs_weather =
            mean(weather_rrn_gain_vs_weather),
        mean_relative_wq_gain = mean(relative_wq_gain_vs_weather_rrn),
        mean_cots_gain = mean(cots_gain_vs_relative),
        mean_full_gain_vs_weather = mean(full_gain_vs_weather),
        .groups = 'drop'
    )
write_csv(strata, file.path(output_dir, 'mortality_brt_pressure_strata.csv'))

severe_2024 <- paired |>
    filter(
        scheme == 'reef_blocked_2024',
        event_year == 2024L,
        observed_mortality >= 0.2
    ) |>
    mutate(
        core_underprediction = observed_mortality - predicted_mortality_core,
        wq_underprediction = observed_mortality - predicted_mortality_rrn_wq,
        cyclone_underprediction =
            observed_mortality - predicted_mortality_rrn_cyclone,
        weather_underprediction =
            observed_mortality - predicted_mortality_weather,
        weather_rrn_underprediction =
            observed_mortality - predicted_mortality_weather_rrn,
        relative_wq_underprediction =
            observed_mortality - predicted_mortality_weather_rrn_relative,
        full_underprediction =
            observed_mortality - predicted_mortality_weather_rrn_full
    ) |>
    arrange(desc(core_underprediction))
write_csv(severe_2024, file.path(output_dir, 'severe_2024_reef_diagnostics.csv'))

event_exposure <- paired |>
    filter(scheme == 'reef_blocked_2024') |>
    group_by(event_year, region_block) |>
    summarise(
        observations = n(),
        reefs = n_distinct(ReefID),
        median_wq = median(wqc_freqcc12),
        median_relative_wq = median(wqc_prior10_percentile),
        high_wq_observations = sum(wqc_freqcc12 >= 0.5),
        high_relative_wq_observations = sum(
            wqc_prior10_percentile >= 0.9
        ),
        cyclone_exposed_observations = sum(cyc_maxHrs4mw > 0),
        maximum_cyclone_hours = max(cyc_maxHrs4mw),
        median_modelled_cots = median(cot_idwmeanpertow),
        observed_mortality = mean(observed_mortality),
        core_prediction = mean(predicted_mortality_core),
        wq_prediction = mean(predicted_mortality_rrn_wq),
        cyclone_prediction = mean(predicted_mortality_rrn_cyclone),
        weather_prediction = mean(predicted_mortality_weather),
        weather_rrn_prediction = mean(predicted_mortality_weather_rrn),
        weather_rrn_relative_prediction = mean(
            predicted_mortality_weather_rrn_relative
        ),
        weather_rrn_full_prediction = mean(
            predicted_mortality_weather_rrn_full
        ),
        .groups = 'drop'
    )
write_csv(event_exposure, file.path(output_dir, 'event_pressure_predictions.csv'))

print(overall_metrics |> arrange(scheme, rmse))
print(
    programme_metrics |>
        filter(scheme == 'reef_blocked_2024') |>
        arrange(programme_key, rmse)
)
