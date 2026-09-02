# Transparent first benchmark for generalising measured AIMS salinity to GBR
# reefs. The main validation keeps 2024 in the training pool but holds out whole
# reefs; a separate train-pre-2024 sensitivity shows why 2024 cannot be the only
# operational holdout for this compound event.

suppressPackageStartupMessages({
    library(dplyr)
    library(mgcv)
    library(readr)
})

input_file <- 'data/processed/aims_freshwater_training_points.csv'
output_dir <- 'output/aims_freshwater_proxy'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read_csv(input_file, show_col_types = FALSE) |>
    filter(
        event_year %in% c(2016, 2017, 2020, 2022, 2024),
        reef_proximal_20km
    ) |>
    mutate(
        flood30 = salinity_min_psu < 30,
        freshwater_deficit30 = pmax(30 - salinity_min_psu, 0),
        log_deficit30 = log1p(freshwater_deficit30),
        log_rain30 = log1p(era5_coastal_rain_dec_mar_max30),
        log_kd490 = log(imos_kd490_q90),
        cyclone_proximity = exp(-pmin(tc_min_distance_km, 1000) / 150),
        observation_depth_class = factor(observation_depth_class)
    ) |>
    filter(if_all(c(
        salinity_min_psu, ereefs_salinity_min, log_rain30, log_kd490,
        rrn_coloured_water, cyclone_proximity, reef_longitude, reef_latitude
    ), is.finite))

reef_folds <- data |>
    distinct(ReefID) |>
    arrange(ReefID) |>
    mutate(fold = (row_number() - 1L) %% 5L + 1L)
data <- data |>
    left_join(reef_folds, by = 'ReefID', relationship = 'many-to-one')

occurrence_formulas <- list(
    'Spatial prior' = flood30 ~ observation_depth_class +
        s(reef_longitude, reef_latitude, k = 30),
    'eReefs calibrated' = flood30 ~ observation_depth_class +
        s(ereefs_salinity_min, k = 5, bs = 'cr') +
        s(reef_longitude, reef_latitude, k = 30),
    'Operational proxies' = flood30 ~ observation_depth_class +
        s(log_rain30, k = 4, bs = 'cr') +
        s(log_kd490, k = 4, bs = 'cr') +
        s(rrn_coloured_water, k = 4, bs = 'cr') +
        s(cyclone_proximity, k = 4, bs = 'cr') +
        s(reef_longitude, reef_latitude, k = 30),
    'Hybrid' = flood30 ~ observation_depth_class +
        s(ereefs_salinity_min, k = 5, bs = 'cr') +
        s(log_rain30, k = 4, bs = 'cr') +
        s(log_kd490, k = 4, bs = 'cr') +
        s(rrn_coloured_water, k = 4, bs = 'cr') +
        s(cyclone_proximity, k = 4, bs = 'cr') +
        s(reef_longitude, reef_latitude, k = 30)
)
magnitude_formulas <- list(
    'Spatial prior' = log_deficit30 ~ observation_depth_class +
        s(reef_longitude, reef_latitude, k = 15),
    'eReefs calibrated' = log_deficit30 ~ observation_depth_class +
        s(ereefs_salinity_min, k = 4, bs = 'cr') +
        s(reef_longitude, reef_latitude, k = 15),
    'Operational proxies' = log_deficit30 ~ observation_depth_class +
        log_rain30 + log_kd490 + rrn_coloured_water +
        cyclone_proximity + s(reef_longitude, reef_latitude, k = 15),
    'Hybrid' = log_deficit30 ~ observation_depth_class +
        s(ereefs_salinity_min, k = 4, bs = 'cr') +
        log_rain30 + log_kd490 + rrn_coloured_water +
        cyclone_proximity + s(reef_longitude, reef_latitude, k = 15)
)

occurrence_predictions <- list()
magnitude_predictions <- list()
prediction_index <- 1L

for (held_out_fold in 1:5) {
    training <- filter(data, fold != held_out_fold)
    testing <- filter(data, fold == held_out_fold)
    magnitude_training <- filter(training, flood30)
    magnitude_testing <- filter(testing, flood30)

    for (model_name in names(occurrence_formulas)) {
        occurrence_model <- gam(
            occurrence_formulas[[model_name]], data = training,
            family = binomial(), method = 'REML', select = TRUE
        )
        magnitude_model <- gam(
            magnitude_formulas[[model_name]], data = magnitude_training,
            method = 'REML', select = TRUE
        )
        occurrence_predictions[[prediction_index]] <- testing |>
            transmute(
                ReefID, event_year, fold = held_out_fold,
                model = model_name, flood30, salinity_min_psu,
                observation_depth_class,
                predicted = as.numeric(predict(
                    occurrence_model, newdata = testing, type = 'response'
                ))
            )
        predicted_log_deficit <- as.numeric(predict(
            magnitude_model, newdata = magnitude_testing, type = 'response'
        ))
        magnitude_predictions[[prediction_index]] <- magnitude_testing |>
            transmute(
                ReefID, event_year, fold = held_out_fold,
                model = model_name, flood30, salinity_min_psu,
                observation_depth_class,
                predicted = pmax(
                    0, pmin(30, 30 - expm1(predicted_log_deficit))
                )
            )
        prediction_index <- prediction_index + 1L
    }
}

occurrence_predictions <- bind_rows(occurrence_predictions)
magnitude_predictions <- bind_rows(magnitude_predictions)

average_precision <- function(truth, estimate) {
    ordered <- order(estimate, decreasing = TRUE)
    truth <- as.integer(truth[ordered])
    if (sum(truth) == 0L) return(NA_real_)
    cumulative_true <- cumsum(truth)
    recall <- cumulative_true / sum(truth)
    precision <- cumulative_true / seq_along(truth)
    sum(diff(c(0, recall)) * precision)
}

occurrence_metrics <- occurrence_predictions |>
    group_by(model, evaluation = if_else(
        event_year == 2024, 'Reef-blocked 2024', 'All reef-blocked rows'
    )) |>
    filter(
        evaluation == 'All reef-blocked rows' |
            event_year == 2024
    ) |>
    summarise(
        rows = n(),
        flood30_rows = sum(flood30),
        prevalence = mean(flood30),
        average_precision = average_precision(flood30, predicted),
        brier = mean((predicted - flood30)^2),
        balanced_brier = 0.5 * (
            mean((predicted[flood30] - 1)^2) +
                mean(predicted[!flood30]^2)
        ),
        top_10pct_recall = sum(flood30[
            order(predicted, decreasing = TRUE)[
                seq_len(ceiling(0.1 * n()))
            ]
        ]) / sum(flood30),
        .groups = 'drop'
    )

# Recalculate the overall rows explicitly because the grouped expression above
# creates the 2024 subset and the non-2024 complement.
overall_occurrence <- occurrence_predictions |>
    group_by(model) |>
    summarise(
        evaluation = 'All reef-blocked rows', rows = n(),
        flood30_rows = sum(flood30), prevalence = mean(flood30),
        average_precision = average_precision(flood30, predicted),
        brier = mean((predicted - flood30)^2),
        balanced_brier = 0.5 * (
            mean((predicted[flood30] - 1)^2) +
                mean(predicted[!flood30]^2)
        ),
        top_10pct_recall = sum(flood30[
            order(predicted, decreasing = TRUE)[seq_len(ceiling(0.1 * n()))]
        ]) / sum(flood30),
        .groups = 'drop'
    )
occurrence_metrics <- occurrence_metrics |>
    filter(evaluation == 'Reef-blocked 2024') |>
    bind_rows(overall_occurrence)

magnitude_metrics <- magnitude_predictions |>
    mutate(evaluation = if_else(
        event_year == 2024, 'Reef-blocked 2024', 'Other reef-blocked rows'
    )) |>
    group_by(model, evaluation) |>
    summarise(
        rows = n(),
        conditional_rmse_psu = sqrt(mean(
            (predicted - salinity_min_psu)^2
        )),
        conditional_mae_psu = mean(abs(predicted - salinity_min_psu)),
        conditional_bias_psu = mean(predicted - salinity_min_psu),
        .groups = 'drop'
    )
overall_magnitude <- magnitude_predictions |>
    group_by(model) |>
    summarise(
        evaluation = 'All reef-blocked rows', rows = n(),
        conditional_rmse_psu = sqrt(mean(
            (predicted - salinity_min_psu)^2
        )),
        conditional_mae_psu = mean(abs(predicted - salinity_min_psu)),
        conditional_bias_psu = mean(predicted - salinity_min_psu),
        .groups = 'drop'
    )
magnitude_metrics <- bind_rows(magnitude_metrics, overall_magnitude)

write_csv(
    occurrence_predictions,
    file.path(output_dir, 'reef_blocked_occurrence_predictions.csv')
)
write_csv(
    magnitude_predictions,
    file.path(output_dir, 'reef_blocked_magnitude_predictions.csv')
)
write_csv(
    occurrence_metrics,
    file.path(output_dir, 'occurrence_metrics.csv')
)
write_csv(
    magnitude_metrics,
    file.path(output_dir, 'magnitude_metrics.csv')
)

cat('Calibration rows:', nrow(data), '\n')
cat('Unique reefs:', n_distinct(data$ReefID), '\n')
cat('Below-30 rows:', sum(data$flood30), '\n')
print(filter(occurrence_metrics, evaluation == 'All reef-blocked rows'))
print(filter(magnitude_metrics, evaluation == 'All reef-blocked rows'))
