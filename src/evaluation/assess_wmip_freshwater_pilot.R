# Map and screen the isolated WMIP freshwater-source pilot. This script does
# not fit or alter a mortality model. Its response is direct AIMS salinity.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
})

input_file <- 'data/processed/aims_freshwater_training_points.csv'
reef_file <- 'data/processed/freshwater_reef_event.csv'
output_dir <- 'output/wmip_freshwater_pilot'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(file.exists(c(input_file, reef_file)))) {
    stop('Run fetch_wmip_discharge.py and build_aims_freshwater_calibration.R first.')
}

safe_scale <- function(training, testing, variables) {
    for (variable in variables) {
        centre <- mean(training[[variable]], na.rm = TRUE)
        spread <- sd(training[[variable]], na.rm = TRUE)
        if (!is.finite(spread) || spread == 0) spread <- 1
        name <- paste0(variable, '_z')
        training[[name]] <- (training[[variable]] - centre) / spread
        testing[[name]] <- (testing[[variable]] - centre) / spread
    }
    list(training = training, testing = testing)
}

average_precision <- function(truth, estimate) {
    keep <- is.finite(estimate) & !is.na(truth)
    truth <- as.integer(truth[keep])
    estimate <- estimate[keep]
    if (length(truth) == 0 || sum(truth) == 0) return(NA_real_)
    ordered <- order(estimate, decreasing = TRUE)
    truth <- truth[ordered]
    estimate <- estimate[ordered]
    ends <- c(which(diff(estimate) != 0), length(estimate))
    recall <- cumsum(truth)[ends] / sum(truth)
    precision <- cumsum(truth)[ends] / ends
    sum(diff(c(0, recall)) * precision)
}

balanced_brier <- function(truth, estimate) {
    positive <- truth
    if (!any(positive) || !any(!positive)) return(NA_real_)
    0.5 * (
        mean((estimate[positive] - 1)^2) +
            mean(estimate[!positive]^2)
    )
}

aims <- read_csv(input_file, show_col_types = FALSE) |>
    mutate(
        log_rain = log1p(era5_coastal_rain_dec_mar_max30),
        log_kd490 = log(pmax(imos_kd490_q90, 1e-6)),
        log_discharge = log1p(wmip_routed_discharge_total_ml),
        low_salinity = salinity_min_psu < 30
    ) |>
    filter(
        reef_proximal_20km,
        wmip_connection_weight_sum > 0,
        if_all(c(log_rain, log_kd490, log_discharge), is.finite),
        !is.na(low_salinity)
    )

coverage <- aims |>
    group_by(event_year, wmip_dominant_source_id) |>
    summarise(
        aims_rows = n(),
        aims_sites = n_distinct(LOCATION_NAME),
        low_salinity_rows = sum(low_salinity),
        minimum_salinity_psu = min(salinity_min_psu),
        median_routed_discharge_ml = median(wmip_routed_discharge_total_ml),
        .groups = 'drop'
    )
write_csv(coverage, file.path(output_dir, 'validation_coverage.csv'))

variables <- c('log_rain', 'log_kd490', 'log_discharge')
predictions <- list()
prediction_index <- 1L
for (held_out_event in sort(unique(aims$event_year))) {
    training <- filter(aims, event_year != held_out_event)
    testing <- filter(aims, event_year == held_out_event)
    if (
        nrow(testing) == 0 || nrow(training) == 0 ||
        n_distinct(training$low_salinity) < 2
    ) next
    prepared <- safe_scale(training, testing, variables)
    training <- prepared$training
    testing <- prepared$testing
    models <- list(
        `Rainfall + Kd490` = glm(
            low_salinity ~ log_rain_z + log_kd490_z,
            data = training, family = binomial()
        ),
        `Rainfall + Kd490 + routed discharge` = glm(
            low_salinity ~ log_rain_z + log_kd490_z + log_discharge_z,
            data = training, family = binomial()
        )
    )
    for (model_name in names(models)) {
        predictions[[prediction_index]] <- testing |>
            transmute(
                held_out_event,
                model = model_name,
                event_year,
                ReefID,
                LOCATION_NAME,
                source = wmip_dominant_source_id,
                observed_salinity_psu = salinity_min_psu,
                low_salinity,
                predicted = as.numeric(predict(
                    models[[model_name]], newdata = testing,
                    type = 'response'
                ))
            )
        prediction_index <- prediction_index + 1L
    }
}

predictions <- bind_rows(predictions)
if (nrow(predictions) > 0) {
    metrics <- predictions |>
        group_by(model, held_out_event) |>
        summarise(
            rows = n(),
            low_salinity_rows = sum(low_salinity),
            average_precision = average_precision(low_salinity, predicted),
            balanced_brier = balanced_brier(low_salinity, predicted),
            .groups = 'drop'
        )
    event_mean <- metrics |>
        group_by(model) |>
        summarise(
            held_out_event = 'event_mean',
            rows = sum(rows),
            low_salinity_rows = sum(low_salinity_rows),
            average_precision = mean(average_precision, na.rm = TRUE),
            balanced_brier = mean(balanced_brier, na.rm = TRUE),
            .groups = 'drop'
        )
    metrics <- bind_rows(
        mutate(metrics, held_out_event = as.character(held_out_event)),
        event_mean
    )
} else {
    metrics <- tibble(
        model = character(), held_out_event = character(), rows = integer(),
        low_salinity_rows = integer(), average_precision = numeric(),
        balanced_brier = numeric()
    )
}
write_csv(predictions, file.path(output_dir, 'event_held_out_predictions.csv'))
write_csv(metrics, file.path(output_dir, 'event_held_out_metrics.csv'))

reef_2024 <- read_csv(reef_file, show_col_types = FALSE) |>
    filter(event_year == 2024, cutoff_name == 'march')
component_map <- reef_2024 |>
    select(
        reef_longitude, reef_latitude,
        `Routed discharge (log ML)` = wmip_routed_discharge_total_ml,
        `ERA5 coastal rain (mm)` = era5_coastal_rain_dec_mar_max30,
        `IMOS Kd490 q90` = imos_kd490_q90,
        `IMOS chlorophyll median` = chla_wetseason_median
    ) |>
    mutate(`Routed discharge (log ML)` = log1p(`Routed discharge (log ML)`)) |>
    pivot_longer(
        -c(reef_longitude, reef_latitude),
        names_to = 'component', values_to = 'value'
    ) |>
    group_by(component) |>
    mutate(component_percentile = percent_rank(value)) |>
    ungroup()

map_plot <- ggplot(
    component_map,
    aes(reef_longitude, reef_latitude, colour = component_percentile)
) +
    geom_point(size = 0.55) +
    coord_equal() +
    facet_wrap(~ component) +
    scale_colour_viridis_c(option = 'C', na.value = 'grey90') +
    labs(
        x = 'Longitude', y = 'Latitude', colour = 'Within-layer\npercentile',
        title = 'March 2024 freshwater pilot components',
        subtitle = paste(
            'Discharge uses provisional distance routing;',
            'components remain separate and are not salinity'
        )
    ) +
    theme_bw(base_size = 9)
ggsave(
    file.path(output_dir, 'component_map_2024.png'), map_plot,
    width = 10, height = 8, dpi = 180
)

decision <- tibble(
    decision = 'not_promoted',
    claims_salinity = FALSE,
    mortality_refit_allowed = FALSE,
    validation_rows = nrow(aims),
    validation_events = n_distinct(aims$event_year),
    low_salinity_rows = sum(aims$low_salinity),
    reason = paste(
        'Pilot ingestion and component mapping only; outlet positions, ratings,',
        'regulation and routing require review before proxy promotion.'
    )
)
write_csv(decision, file.path(output_dir, 'decision.csv'))

print(coverage)
print(metrics)
print(decision)
