# Test whether ERA5 rainfall, wind, and IMOS optics improve reconstruction of
# eReefs low-salinity exposure. Thirty PSU is the primary screening threshold;
# 28 PSU is retained as a severe-tail check.

suppressPackageStartupMessages({
    library(dplyr)
    library(mgcv)
    library(readr)
})

salinity_file <- 'data/processed/ereefs_salinity_imos_kd490_reef_year.csv'
weather_file <- 'data/processed/era5_weather_reef_year.csv'
output_dir <- 'output/salinity_rainfall'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(file.exists(c(salinity_file, weather_file)))) {
    stop('Run the eReefs salinity and ERA5 weather extractors first.')
}

keys <- c('LABEL_ID', 'LOC_NAME_S', 'lon', 'lat', 'year')
salinity <- read_csv(salinity_file, show_col_types = FALSE)
weather <- read_csv(weather_file, show_col_types = FALSE)
combined <- salinity |>
    left_join(weather, by = keys, relationship = 'one-to-one')
if (anyNA(combined$era5_reef_rain_q1_total)) {
    stop('ERA5 weather failed to join one-to-one with salinity rows')
}

# Multiple registry features can share one 4-km eReefs cell. Use one row per
# cell-year so duplicated salinity values do not create pseudo-replication.
cell_year <- combined |>
    group_by(ereefs_lat_index, ereefs_lon_index, year) |>
    summarise(
        lon = first(ereefs_longitude),
        lat = first(ereefs_latitude),
        n_reef_features = n(),
        salinity_min = first(salinity_min),
        freshwater_exposure_30 = first(freshwater_exposure_30),
        freshwater_exposure_28 = first(freshwater_exposure_28),
        days_below_30 = first(days_below_30),
        days_below_28 = first(days_below_28),
        salinity_complete_q1 = first(salinity_complete_q1),
        salinity_forcing_warning = first(salinity_forcing_warning),
        k490_q90 = median(k490_q90, na.rm = TRUE),
        coastal_distance_km = median(coastal_grid_distance_km),
        coastal_rain_december = median(
            era5_coastal_rain_december_total
        ),
        coastal_rain_q1 = median(era5_coastal_rain_q1_total),
        coastal_rain_max_7day = median(era5_coastal_rain_q1_max_7day),
        coastal_rain_max_30day = median(
            era5_coastal_rain_dec_mar_max_30day
        ),
        reef_rain_q1 = median(era5_reef_rain_q1_total),
        reef_wind_mean = median(era5_reef_wind_q1_mean),
        reef_wind_calm3 = median(era5_reef_wind_q1_fraction_below_3),
        .groups = 'drop'
    ) |>
    mutate(
        flood30 = freshwater_exposure_30 > 0,
        flood28 = freshwater_exposure_28 > 0,
        log_kd = log(k490_q90),
        log_coastal_rain30 = log1p(coastal_rain_max_30day),
        log_coastal_rain7 = log1p(coastal_rain_max_7day),
        log_reef_rain = log1p(reef_rain_q1)
    )

historical <- cell_year |>
    filter(
        year %in% c(2016, 2017, 2020),
        salinity_complete_q1,
        !salinity_forcing_warning,
        if_all(
            c(
                log_kd, log_coastal_rain30, log_coastal_rain7,
                log_reef_rain, reef_wind_mean
            ),
            is.finite
        )
    )

scale_fold <- function(training, testing, variables) {
    for (variable in variables) {
        centre <- mean(training[[variable]])
        spread <- sd(training[[variable]])
        if (!is.finite(spread) || spread == 0) spread <- 1
        training[[paste0(variable, '_z')]] <-
            (training[[variable]] - centre) / spread
        testing[[paste0(variable, '_z')]] <-
            (testing[[variable]] - centre) / spread
    }
    list(training = training, testing = testing)
}

average_precision <- function(truth, estimate) {
    ordered <- order(estimate, decreasing = TRUE)
    truth <- as.integer(truth[ordered])
    estimate <- estimate[ordered]
    if (sum(truth) == 0) return(NA_real_)
    ends <- c(which(diff(estimate) != 0), length(estimate))
    true_positive <- cumsum(truth)[ends]
    recall <- true_positive / sum(truth)
    precision <- true_positive / ends
    sum(diff(c(0, recall)) * precision)
}

top_fraction_recall <- function(truth, estimate, fraction) {
    if (sum(truth) == 0) return(NA_real_)
    selected <- order(estimate, decreasing = TRUE)[
        seq_len(max(1, ceiling(length(truth) * fraction)))
    ]
    sum(truth[selected]) / sum(truth)
}

variables <- c(
    'log_coastal_rain30', 'log_coastal_rain7', 'log_reef_rain',
    'log_kd', 'reef_wind_mean'
)
occurrence_predictions <- list()
magnitude_predictions <- list()
index <- 1L

for (held_out_year in c(2016, 2017, 2020)) {
    training <- filter(historical, year != held_out_year)
    testing <- filter(historical, year == held_out_year)
    prepared <- scale_fold(training, testing, variables)
    training <- prepared$training
    testing <- prepared$testing

    occurrence_models <- list(
        'Intercept only' = glm(
            flood30 ~ 1, data = training, family = binomial()
        ),
        'Spatial prior' = bam(
            flood30 ~ s(lon, lat, k = 40),
            data = training, family = binomial(), method = 'fREML',
            select = TRUE, discrete = TRUE
        ),
        'Rainfall + spatial' = bam(
            flood30 ~ s(log_coastal_rain30_z, k = 4, bs = 'cr') +
                s(log_coastal_rain7_z, k = 4, bs = 'cr') +
                s(lon, lat, k = 40),
            data = training, family = binomial(), method = 'fREML',
            select = TRUE, discrete = TRUE
        ),
        'Rainfall + optics + wind + spatial' = bam(
            flood30 ~ s(log_coastal_rain30_z, k = 4, bs = 'cr') +
                s(log_coastal_rain7_z, k = 4, bs = 'cr') +
                s(log_kd_z, k = 4, bs = 'cr') +
                s(reef_wind_mean_z, k = 4, bs = 'cr') +
                s(lon, lat, k = 40),
            data = training, family = binomial(), method = 'fREML',
            select = TRUE, discrete = TRUE
        )
    )
    for (model_name in names(occurrence_models)) {
        occurrence_predictions[[index]] <- testing |>
            transmute(
                held_out_year,
                model = model_name,
                flood30,
                flood28,
                observed_exposure_30 = freshwater_exposure_30,
                predicted = as.numeric(predict(
                    occurrence_models[[model_name]],
                    newdata = testing,
                    type = 'response'
                ))
            )
        index <- index + 1L
    }

    magnitude_models <- list(
        'Spatial prior' = bam(
            salinity_min ~ s(lon, lat, k = 40),
            data = training, method = 'fREML', select = TRUE,
            discrete = TRUE
        ),
        'Rainfall + optics + wind + spatial' = bam(
            salinity_min ~ s(log_coastal_rain30_z, k = 4, bs = 'cr') +
                s(log_coastal_rain7_z, k = 4, bs = 'cr') +
                s(log_kd_z, k = 4, bs = 'cr') +
                s(reef_wind_mean_z, k = 4, bs = 'cr') +
                s(lon, lat, k = 40),
            data = training, method = 'fREML', select = TRUE,
            discrete = TRUE
        )
    )
    for (model_name in names(magnitude_models)) {
        magnitude_predictions[[index]] <- testing |>
            transmute(
                held_out_year,
                model = model_name,
                flood30,
                flood28,
                observed = salinity_min,
                predicted = as.numeric(predict(
                    magnitude_models[[model_name]], newdata = testing
                ))
            )
        index <- index + 1L
    }
}

occurrence_predictions <- bind_rows(occurrence_predictions)
magnitude_predictions <- bind_rows(magnitude_predictions)

occurrence_metrics <- occurrence_predictions |>
    group_by(model, held_out_year) |>
    summarise(
        cells = n(),
        flood30_cells = sum(flood30),
        flood28_cells = sum(flood28),
        prevalence30 = mean(flood30),
        average_precision30 = average_precision(flood30, predicted),
        balanced_brier30 = 0.5 * (
            mean((predicted[flood30] - 1)^2) +
                mean(predicted[!flood30]^2)
        ),
        top_1pct_recall30 = top_fraction_recall(flood30, predicted, 0.01),
        top_5pct_recall30 = top_fraction_recall(flood30, predicted, 0.05),
        severe28_top_1pct_recall = top_fraction_recall(
            flood28, predicted, 0.01
        ),
        severe28_top_5pct_recall = top_fraction_recall(
            flood28, predicted, 0.05
        ),
        .groups = 'drop'
    ) |>
    mutate(held_out_year = as.character(held_out_year))

occurrence_event_mean <- occurrence_metrics |>
    group_by(model) |>
    summarise(
        held_out_year = 'Event mean',
        across(
            c(
                prevalence30, average_precision30, balanced_brier30,
                top_1pct_recall30, top_5pct_recall30,
                severe28_top_1pct_recall, severe28_top_5pct_recall
            ),
            ~ mean(.x, na.rm = TRUE)
        ),
        cells = sum(cells),
        flood30_cells = sum(flood30_cells),
        flood28_cells = sum(flood28_cells),
        .groups = 'drop'
    )

magnitude_metrics <- magnitude_predictions |>
    group_by(model, held_out_year) |>
    summarise(
        all_cell_mae = mean(abs(predicted - observed)),
        flood30_cell_mae = mean(
            abs(predicted[flood30] - observed[flood30]), na.rm = TRUE
        ),
        flood30_cell_bias = mean(
            predicted[flood30] - observed[flood30], na.rm = TRUE
        ),
        flood30_detection = mean(
            predicted[flood30] < 30, na.rm = TRUE
        ),
        flood28_detection = mean(
            predicted[flood28] < 28, na.rm = TRUE
        ),
        .groups = 'drop'
    )

# Train the full historical occurrence model and create a 2024 sensitivity.
prepared <- scale_fold(historical, filter(cell_year, year == 2024), variables)
training <- prepared$training
future <- prepared$testing
full_model <- bam(
    flood30 ~ s(log_coastal_rain30_z, k = 4, bs = 'cr') +
        s(log_coastal_rain7_z, k = 4, bs = 'cr') +
        s(log_kd_z, k = 4, bs = 'cr') +
        s(reef_wind_mean_z, k = 4, bs = 'cr') +
        s(lon, lat, k = 40),
    data = training, family = binomial(), method = 'fREML',
    select = TRUE, discrete = TRUE
)
future$risk30_probability <- as.numeric(predict(
    full_model, newdata = future, type = 'response'
))
future$risk30_percentile <- rank(
    future$risk30_probability, ties.method = 'average'
) / nrow(future)

risk_2024 <- combined |>
    filter(year == 2024) |>
    left_join(
        future |>
            select(
                ereefs_lat_index, ereefs_lon_index,
                risk30_probability, risk30_percentile
            ),
        by = c('ereefs_lat_index', 'ereefs_lon_index'),
        relationship = 'many-to-one'
    ) |>
    mutate(
        partial_ereefs_flood30 = freshwater_exposure_30 > 0,
        use_in_primary_mortality_model = FALSE,
        reconstruction_status = paste(
            'event-held-out validation required; rainfall uses nearest',
            'ERA5 land-cell proxy, not catchment runoff'
        )
    )

write_csv(cell_year, file.path(output_dir, 'cell_year_inputs.csv'))
write_csv(
    occurrence_predictions,
    file.path(output_dir, 'occurrence_predictions.csv')
)
write_csv(
    bind_rows(occurrence_metrics, occurrence_event_mean),
    file.path(output_dir, 'occurrence_metrics.csv')
)
write_csv(
    magnitude_predictions,
    file.path(output_dir, 'magnitude_predictions.csv')
)
write_csv(
    magnitude_metrics,
    file.path(output_dir, 'magnitude_metrics.csv')
)
write_csv(
    risk_2024,
    'data/processed/ereefs_salinity_2024_rainfall_sensitivity.csv'
)

print(occurrence_event_mean |>
    arrange(desc(average_precision30)))
print(magnitude_metrics)
