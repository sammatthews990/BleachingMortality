# Construct validation of the satellite day-to-following-night contrast
# against subdaily in-situ temperature variability at AIMS WQ logger sites.
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

logger_file <- 'data/processed/aims_wq_logger_hourly_event_summary.csv'
satellite_file <- 'data/processed/imos_thermal_metrics_aims_loggers.csv'
output_dir <- 'output/imos_thermal_screen'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

logger <- read_csv(logger_file, show_col_types = FALSE)
satellite <- read_csv(satellite_file, show_col_types = FALSE)
paired <- logger |>
    inner_join(
        satellite |>
            select(
                ReefID, year, imos_day_to_following_night_drop_c,
                imos_night_relief_fraction,
                imos_hot_day_night_pair_count,
                imos_hot_day_night_pair_coverage
            ),
        by = c('site' = 'ReefID', 'event_year' = 'year'),
        relationship = 'one-to-one'
    ) |>
    group_by(event_year) |>
    mutate(
        satellite_drop_event_centred =
            imos_day_to_following_night_drop_c -
            mean(imos_day_to_following_night_drop_c, na.rm = TRUE),
        logger_drop_event_centred =
            logger_peak30_median_nocturnal_drop_c -
            mean(logger_peak30_median_nocturnal_drop_c, na.rm = TRUE),
        logger_dtr_event_centred = logger_peak30_mean_dtr_c -
            mean(logger_peak30_mean_dtr_c, na.rm = TRUE)
    ) |>
    ungroup()

safe_spearman <- function(x, y) {
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 4L) return(NA_real_)
    cor(x[keep], y[keep], method = 'spearman')
}
comparison <- tibble(
    comparison = c(
        'satellite drop vs logger nocturnal drop',
        'satellite drop vs logger DTR',
        'event-centred satellite vs logger nocturnal drop',
        'event-centred satellite vs logger DTR'
    ),
    n = c(
        sum(is.finite(paired$imos_day_to_following_night_drop_c) &
            is.finite(paired$logger_peak30_median_nocturnal_drop_c)),
        sum(is.finite(paired$imos_day_to_following_night_drop_c) &
            is.finite(paired$logger_peak30_mean_dtr_c)),
        sum(is.finite(paired$satellite_drop_event_centred) &
            is.finite(paired$logger_drop_event_centred)),
        sum(is.finite(paired$satellite_drop_event_centred) &
            is.finite(paired$logger_dtr_event_centred))
    ),
    spearman_rho = c(
        safe_spearman(
            paired$imos_day_to_following_night_drop_c,
            paired$logger_peak30_median_nocturnal_drop_c
        ),
        safe_spearman(
            paired$imos_day_to_following_night_drop_c,
            paired$logger_peak30_mean_dtr_c
        ),
        safe_spearman(
            paired$satellite_drop_event_centred,
            paired$logger_drop_event_centred
        ),
        safe_spearman(
            paired$satellite_drop_event_centred,
            paired$logger_dtr_event_centred
        )
    )
)

coverage <- paired |>
    summarise(
        logger_site_events = n(),
        logger_sites = n_distinct(site),
        event_years = n_distinct(event_year),
        finite_satellite_cooling = sum(is.finite(
            imos_day_to_following_night_drop_c
        )),
        finite_logger_nocturnal_drop = sum(is.finite(
            logger_peak30_median_nocturnal_drop_c
        )),
        finite_logger_dtr = sum(is.finite(logger_peak30_mean_dtr_c))
    )

write_csv(paired, file.path(output_dir, 'logger_cooling_pairs.csv'))
write_csv(comparison, file.path(output_dir, 'logger_cooling_validation.csv'))
write_csv(coverage, file.path(output_dir, 'logger_cooling_coverage.csv'))
print(coverage)
print(comparison)
