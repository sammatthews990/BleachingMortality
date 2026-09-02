# Spatial and temporal audit of freshwater evidence around the Lizard Island
# mortality misses. This distinguishes continuous loggers from discrete AIMS
# water-quality samples; the latter can show a pulse but cannot establish its
# duration or minimum salinity.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

output_dir <- 'output/low_dhw_sensitivity'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

haversine_km <- function(lon1, lat1, lon2, lat2) {
    rad <- pi / 180
    6371 * 2 * asin(sqrt(
        sin((lat2 - lat1) * rad / 2)^2 +
            cos(lat1 * rad) * cos(lat2 * rad) *
            sin((lon2 - lon1) * rad / 2)^2
    ))
}

# Use the North West Lizard reef observation as the spatial reference.
outliers <- read_csv(
    'output/prediction_outliers/priority_observation_rows.csv',
    show_col_types = FALSE
)
lizard_reference <- outliers |>
    filter(grepl('Lizard Island Reef \\(North West\\)', ReefName)) |>
    slice(1)
if (nrow(lizard_reference) != 1L) {
    stop('Could not identify the Lizard Island North West reference row.')
}
lizard_lon <- lizard_reference$lon[[1]]
lizard_lat <- lizard_reference$lat[[1]]

lizard_cluster <- outliers |>
    filter(event_year == 2024) |>
    mutate(
        distance_to_lizard_nw_km = haversine_km(
            lizard_lon, lizard_lat, lon, lat
        )
    ) |>
    filter(distance_to_lizard_nw_km <= 100) |>
    select(
        programme_key, ReefID, ReefName, lon, lat,
        distance_to_lizard_nw_km, ann_maxdhw,
        observed_mortality, predicted_mortality, residual,
        era5_coastal_rain_max_30day, freshwater_risk30_percentile,
        partial_ereefs_salinity_min, partial_ereefs_flood30,
        DISTURBANCE_TYPE, description
    ) |>
    arrange(desc(residual))

logger_sites <- read_csv(
    'data/processed/aims_wq_logger_sites.csv', show_col_types = FALSE
) |>
    mutate(
        distance_to_lizard_nw_km = haversine_km(
            lizard_lon, lizard_lat, site_longitude, site_latitude
        )
    ) |>
    arrange(distance_to_lizard_nw_km)

logger_coverage <- read_csv(
    'data/processed/aims_wq_logger_hourly_coverage.csv',
    show_col_types = FALSE
)
logger_audit <- logger_sites |>
    left_join(logger_coverage, by = c(
        'site', 'site_longitude', 'site_latitude'
    )) |>
    select(
        site, site_longitude, site_latitude, distance_to_lizard_nw_km,
        first_time, last_time, hourly_rows, salinity_hours
    )

samples <- read_csv(
    'data/processed/aims_wq_discrete_samples_2015_2025.csv',
    show_col_types = FALSE
) |>
    mutate(
        sample_datetime = as.POSIXct(sample_datetime, tz = 'UTC'),
        distance_to_lizard_nw_km = haversine_km(
            lizard_lon, lizard_lat, LONGITUDE, LATITUDE
        )
    )

# Retain the wet-season samples close enough to resolve Lizard, Linnet and Eyrie.
lizard_samples <- samples |>
    filter(
        distance_to_lizard_nw_km <= 30,
        sample_datetime >= as.POSIXct('2023-11-01', tz = 'UTC'),
        sample_datetime <= as.POSIXct('2024-04-30 23:59:59', tz = 'UTC')
    ) |>
    mutate(sample_date = as.Date(sample_datetime)) |>
    group_by(
        PROJECT, LOCATION_NAME, SHORT_NAME, sample_date,
        LATITUDE, LONGITUDE, distance_to_lizard_nw_km
    ) |>
    summarise(
        depths_sampled = n_distinct(SAMPLE_DEPTH[is.finite(SAMPLE_DEPTH)]),
        salinity_min_psu = if (all(!is.finite(SAL))) NA_real_ else min(SAL, na.rm = TRUE),
        salinity_mean_psu = if (all(!is.finite(SAL))) NA_real_ else mean(SAL, na.rm = TRUE),
        chlorophyll_mean = if (all(!is.finite(CHL))) NA_real_ else mean(CHL, na.rm = TRUE),
        chlorophyll_max = if (all(!is.finite(CHL))) NA_real_ else max(CHL, na.rm = TRUE),
        suspended_solids_max = if (all(!is.finite(SS))) NA_real_ else max(SS, na.rm = TRUE),
        secchi_min_m = if (all(!is.finite(SECCHI_DEPTH))) NA_real_ else min(SECCHI_DEPTH, na.rm = TRUE),
        .groups = 'drop'
    ) |>
    arrange(sample_date, distance_to_lizard_nw_km)

# Before/after contrasts are possible at the two repeatedly sampled island sites.
repeated_site_contrast <- lizard_samples |>
    filter(grepl('Lizard Island north|MacGillivray', LOCATION_NAME, ignore.case = TRUE)) |>
    group_by(LOCATION_NAME) |>
    arrange(sample_date, .by_group = TRUE) |>
    summarise(
        first_date = first(sample_date),
        last_date = last(sample_date),
        first_salinity = first(salinity_min_psu),
        last_salinity = last(salinity_min_psu),
        salinity_change_psu = last_salinity - first_salinity,
        first_chlorophyll = first(chlorophyll_mean),
        last_chlorophyll = last(chlorophyll_mean),
        chlorophyll_ratio = last_chlorophyll / first_chlorophyll,
        .groups = 'drop'
    )

write_csv(lizard_cluster, file.path(output_dir, 'lizard_cluster_outliers.csv'), na = '')
write_csv(logger_audit, file.path(output_dir, 'lizard_logger_distance_audit.csv'), na = '')
write_csv(lizard_samples, file.path(output_dir, 'lizard_wet_season_discrete_samples.csv'), na = '')
write_csv(repeated_site_contrast, file.path(output_dir, 'lizard_repeated_site_contrast.csv'), na = '')

print(lizard_cluster)
print(head(logger_audit, 3))
print(lizard_samples)
print(repeated_site_contrast)
