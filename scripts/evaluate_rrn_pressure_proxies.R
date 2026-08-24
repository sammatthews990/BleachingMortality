# Evaluate GBRMPA water-colour exposure as a freshwater/plume proxy against
# complete eReefs salinity years. Cyclone wave exposure is summarised here but
# is not treated as a salinity proxy.

suppressPackageStartupMessages({
    library(dplyr)
    library(mgcv)
    library(readr)
    library(tidyr)
})

output_dir <- 'output/rrn_pressure_assessment'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rrn <- read_csv(
    'data/processed/rrn_pressure_reef_year.csv', show_col_types = FALSE
)
salinity <- read_csv(
    'data/processed/ereefs_salinity_imos_kd490_reef_year.csv',
    show_col_types = FALSE
)
weather <- read_csv(
    'data/processed/era5_weather_reef_year.csv', show_col_types = FALSE
) |>
    group_by(year, grid_lat, grid_lon) |>
    summarise(
        coastal_rain30 = median(
            era5_coastal_rain_dec_mar_max_30day, na.rm = TRUE
        ),
        .groups = 'drop'
    )

proxy_data <- salinity |>
    filter(salinity_complete_q1 %in% TRUE) |>
    mutate(
        grid_lat = floor(lat / 0.25 + 0.5) * 0.25,
        grid_lon = floor(lon / 0.25 + 0.5) * 0.25
    ) |>
    left_join(
        rrn,
        by = c('LABEL_ID', 'year' = 'event_year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        weather, by = c('year', 'grid_lat', 'grid_lon'),
        relationship = 'many-to-one'
    ) |>
    mutate(
        flood30 = freshwater_exposure_30 > 0,
        log_coastal_rain30 = log1p(coastal_rain30),
        log_kd490 = log1p(1 / pmax(secc3m_p10, 0.05))
    ) |>
    filter(complete.cases(
        salinity_min, flood30, lon, lat, wqc_freqcc12,
        log_coastal_rain30, log_kd490
    ))

if (nrow(proxy_data) == 0L) stop('No complete RRN/eReefs proxy rows')

coverage <- salinity |>
    filter(salinity_complete_q1 %in% TRUE) |>
    left_join(
        rrn, by = c('LABEL_ID', 'year' = 'event_year'),
        relationship = 'many-to-one'
    ) |>
    group_by(year) |>
    summarise(
        ereefs_reef_cells = n(),
        rrn_wq_cells = sum(!is.na(wqc_freqcc12)),
        rrn_cyclone_cells = sum(!is.na(cyc_maxHrs4mw)),
        overlap_percent = 100 * rrn_wq_cells / ereefs_reef_cells,
        .groups = 'drop'
    )
write_csv(coverage, file.path(output_dir, 'proxy_coverage.csv'))

rank_comparison <- proxy_data |>
    group_by(year) |>
    summarise(
        n = n(),
        flood30_cells = sum(flood30),
        wq_vs_salinity_min = cor(
            wqc_freqcc12, salinity_min, method = 'spearman'
        ),
        wq_vs_days_below30 = cor(
            wqc_freqcc12, days_below_30, method = 'spearman'
        ),
        rainfall_vs_salinity_min = cor(
            log_coastal_rain30, salinity_min, method = 'spearman'
        ),
        kd490_vs_salinity_min = cor(
            log_kd490, salinity_min, method = 'spearman'
        ),
        .groups = 'drop'
    ) |>
    bind_rows(
        proxy_data |>
            summarise(
                year = 0L,
                n = n(),
                flood30_cells = sum(flood30),
                wq_vs_salinity_min = cor(
                    wqc_freqcc12, salinity_min, method = 'spearman'
                ),
                wq_vs_days_below30 = cor(
                    wqc_freqcc12, days_below_30, method = 'spearman'
                ),
                rainfall_vs_salinity_min = cor(
                    log_coastal_rain30, salinity_min, method = 'spearman'
                ),
                kd490_vs_salinity_min = cor(
                    log_kd490, salinity_min, method = 'spearman'
                )
            )
    ) |>
    mutate(year = if_else(year == 0L, 'all_complete_years', as.character(year)))
write_csv(rank_comparison, file.path(output_dir, 'proxy_rank_correlations.csv'))

candidate_terms <- list(
    spatial = 's(lon, lat, k = 40)',
    rainfall = paste(
        's(lon, lat, k = 40)',
        's(log_coastal_rain30, k = 4)', sep = ' + '
    ),
    imos_kd490 = paste(
        's(lon, lat, k = 40)', 's(log_kd490, k = 4)', sep = ' + '
    ),
    rrn_wq = paste(
        's(lon, lat, k = 40)', 's(wqc_freqcc12, k = 4)', sep = ' + '
    ),
    rrn_wq_rainfall = paste(
        's(lon, lat, k = 40)', 's(wqc_freqcc12, k = 4)',
        's(log_coastal_rain30, k = 4)', sep = ' + '
    ),
    all_proxies = paste(
        's(lon, lat, k = 40)', 's(wqc_freqcc12, k = 4)',
        's(log_coastal_rain30, k = 4)', 's(log_kd490, k = 4)', sep = ' + '
    )
)

average_precision <- function(observed, predicted) {
    ordering <- order(predicted, decreasing = TRUE)
    observed <- as.logical(observed[ordering])
    if (!any(observed)) return(NA_real_)
    mean(cumsum(observed)[observed] / which(observed))
}

rank_auc <- function(observed, predicted) {
    observed <- as.logical(observed)
    positives <- sum(observed)
    negatives <- sum(!observed)
    if (positives == 0L || negatives == 0L) return(NA_real_)
    ranks <- rank(predicted, ties.method = 'average')
    (sum(ranks[observed]) - positives * (positives + 1) / 2) /
        (positives * negatives)
}

predictions <- tibble()
for (held_out_year in sort(unique(proxy_data$year))) {
    training <- filter(proxy_data, year != held_out_year)
    assessment <- filter(proxy_data, year == held_out_year)
    for (model_name in names(candidate_terms)) {
        rhs <- candidate_terms[[model_name]]
        occurrence <- bam(
            as.formula(paste('flood30 ~', rhs)),
            data = training, family = binomial(), method = 'fREML',
            discrete = TRUE
        )
        minimum_salinity <- bam(
            as.formula(paste('salinity_min ~', rhs)),
            data = training, method = 'fREML', discrete = TRUE
        )
        predictions <- bind_rows(
            predictions,
            assessment |>
                transmute(
                    LABEL_ID, year, salinity_min, flood30,
                    wqc_freqcc12, cyc_maxHrs4mw,
                    model = model_name,
                    predicted_flood30 = as.numeric(predict(
                        occurrence, newdata = assessment, type = 'response'
                    )),
                    predicted_salinity_min = as.numeric(predict(
                        minimum_salinity, newdata = assessment
                    ))
                )
        )
    }
}

metrics <- predictions |>
    group_by(model) |>
    summarise(
        n = n(),
        events = n_distinct(year),
        flood30_cells = sum(flood30),
        average_precision = average_precision(flood30, predicted_flood30),
        roc_auc = rank_auc(flood30, predicted_flood30),
        balanced_brier = 0.5 * (
            mean((predicted_flood30[flood30] - 1)^2) +
                mean(predicted_flood30[!flood30]^2)
        ),
        salinity_rmse = sqrt(mean(
            (salinity_min - predicted_salinity_min)^2
        )),
        salinity_mae = mean(abs(salinity_min - predicted_salinity_min)),
        salinity_r2 = 1 - sum(
            (salinity_min - predicted_salinity_min)^2
        ) / sum((salinity_min - mean(salinity_min))^2),
        .groups = 'drop'
    ) |>
    arrange(desc(average_precision), balanced_brier)
write_csv(predictions, file.path(output_dir, 'proxy_loyo_predictions.csv'))
write_csv(metrics, file.path(output_dir, 'proxy_loyo_metrics.csv'))

mortality <- read_csv(
    'output/joint_compound_models/joint_compound_rows.csv',
    show_col_types = FALSE
) |>
    left_join(
        rrn,
        by = c('ReefID' = 'LABEL_ID', 'event_year'),
        relationship = 'many-to-one'
    )

mortality_coverage <- mortality |>
    group_by(programme_key, event_year) |>
    summarise(
        observations = n(),
        reefs = n_distinct(ReefID),
        wq_available = sum(!is.na(wqc_freqcc12)),
        cyclone_available = sum(!is.na(cyc_maxHrs4mw)),
        cyclone_exposed = sum(cyc_maxHrs4mw > 0, na.rm = TRUE),
        maximum_cyclone_hours = max(cyc_maxHrs4mw, na.rm = TRUE),
        .groups = 'drop'
    )
write_csv(
    mortality_coverage,
    file.path(output_dir, 'mortality_pressure_coverage.csv')
)

cyclone_events <- mortality |>
    filter(cyc_maxHrs4mw > 0) |>
    group_by(event_year, region_block) |>
    summarise(
        observations = n(), reefs = n_distinct(ReefID),
        median_cyclone_hours = median(cyc_maxHrs4mw),
        maximum_cyclone_hours = max(cyc_maxHrs4mw),
        mean_observed_mortality = mean(mortality_prop),
        .groups = 'drop'
    )
write_csv(cyclone_events, file.path(output_dir, 'cyclone_event_summary.csv'))

print(coverage)
print(metrics)
print(mortality_coverage)

